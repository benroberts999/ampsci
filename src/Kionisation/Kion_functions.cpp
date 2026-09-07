#include "Kion_functions.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracODE/include.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHFcomplex.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Maths.hpp"
#include "qip/String.hpp"
#include "qip/Widgets.hpp"
#include "qip/omp.hpp"
#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <utility>

namespace Kion {

//==============================================================================
AtomicMethod parseStatesMethod(const std::string &in_method) {
  if (qip::ci_compare(in_method, "HF") ||
      qip::ci_compare(in_method, "HartreeFock"))
    return AtomicMethod::HF;
  if (qip::ci_compare(in_method, "Zeff"))
    return AtomicMethod::Zeff;
  if (qip::ci_wc_compare(in_method, "ZeffAn*"))
    return AtomicMethod::ZeffAnalytic;
  if (qip::ci_compare(in_method, "RPA") ||
      qip::ci_wc_compare(in_method, "TDHF*"))
    return AtomicMethod::RPA;
  std::cout << "Warning: AtomicMethod: " << in_method
            << " ?? Defaulting to HF\n";
  return AtomicMethod::HF;
}

std::string parseStatesMethod(const AtomicMethod &in_method) {
  switch (in_method) {
  case AtomicMethod::HF:
    return "HF";
  case AtomicMethod::Zeff:
    return "Zeff";
  case AtomicMethod::ZeffAnalytic:
    return "ZeffAnalytic";
  case AtomicMethod::RPA:
    return "RPA";
  }
  assert(false);
}

//==============================================================================
LinAlg::Matrix<double>
calculateK_nk(const HF::HartreeFock *vHF, const DiracSpinor &Fnk, int max_L,
              const Grid &Egrid, const DiracOperator::jL *jl,
              bool force_rescale, bool hole_particle, bool force_orthog,
              bool zeff_cont, bool zeff_bound, double ec_cut) {
  assert(vHF != nullptr && "Hartree-Fock potential must not be null");
  assert(jl != nullptr && "jl operator must not be null");

  const auto &qgrid = jl->q_grid();
  const auto qsteps = qgrid.num_points();

  LinAlg::Matrix Knk_Eq(Egrid.num_points(), qgrid.num_points());

  if (std::abs(Fnk.en()) > Egrid.r().back()) {
    return Knk_Eq;
  }

  // Same Zeff as used by DarkARC (eqn B35 of arxiv:1912.08204):
  // Zeff = sqrt{I_{njl} eV / 13.6 eV} * n
  // au: Zeff = sqrt{2 * I_{njl}} * n
  // const double Zeff = std::sqrt(-2.0 * Fnk.en()) * Fnk.n();

  const double Zeff = std::sqrt(-2.0 * Fnk.en()) * Fnk.n();
  if (zeff_cont || zeff_bound) {
    std::cout << Fnk << " E = " << Fnk.en() << ", Zeff = " << Zeff << "\n";
  }
  // nb: Fnk_zeff will have not exact right energy; use HF energy
  // This is either real state, or Zeff version.
  // Use this when calculating matrix elements (but not energies)
  const auto &Fnk_t =
    zeff_bound ? DiracSpinor::exactHlike(Fnk.n(), Fnk.kappa(), Fnk.grid_sptr(),
                                         Zeff, vHF->alpha()) :
                 Fnk;

  // Definition of matrix element:
  // matrix element defined such that:
  // K(E,q) = (2L+1) * |me|^2
  // me = <a||jL||e>

  // Find first energy grid point for which Fnk is accessible:
  const auto idE_first_accessible = std::size_t(std::distance(
    Egrid.begin(), std::find_if(Egrid.begin(), Egrid.end(),
                                [&](auto e) { return e > -Fnk.en(); })));
  const auto num_accessible_E_steps = Egrid.num_points() - idE_first_accessible;

  // decide what to parallelise over:
  const bool parallelise_E =
    num_accessible_E_steps >
    std::min(qsteps, (std::size_t)omp_get_max_threads());

  (void)parallelise_E; //suppress unused variable warning clang, when no OMP
#pragma omp parallel for if (parallelise_E)
  for (std::size_t idE = idE_first_accessible; idE < Egrid.num_points();
       ++idE) {
    const auto dE = Egrid(idE);

    // Convert energy deposition to contimuum state energy:
    const double ec = dE + Fnk.en();
    if (ec <= 0.0 || ec > ec_cut)
      continue;

    const int l = Fnk.l();
    const int lc_max = l + max_L;
    const int lc_min = std::max(l - max_L, 0);
    // occupancy fraction. Usually 1. = N(j)/(2j+1)
    const double x_ocf = Fnk.occ_frac();

    // create cntm object [survives locally only]
    ContinuumOrbitals cntm(vHF);
    if (zeff_cont) {
      cntm.solveContinuumZeff(ec, lc_min, lc_max, Zeff, &Fnk_t, force_orthog);
    } else {
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fnk_t, force_rescale,
                            hole_particle, force_orthog);
    }

// Generate AK for each L, lc, and q
// L and lc are summed, not stored individually
#pragma omp parallel for if (!parallelise_E)
    for (std::size_t iq = 0; iq < qsteps; iq++) {
      for (std::size_t L = 0; L <= std::size_t(max_L); L++) {
        for (const auto &Fe : cntm.orbitals) {
          if (jl->is_zero(Fe, Fnk_t, L))
            continue;
          const auto q = jl->q_grid().r(iq);
          const auto me = jl->rme(Fe, Fnk_t, L, q);
          Knk_Eq(idE, iq) += double(2 * L + 1) * me * me * x_ocf;
        }
      }
    }
  }
  return Knk_Eq;
}

//==============================================================================
FormFactorSet allocate_formFactors(std::size_t E_steps, std::size_t q_steps,
                                   bool vectorQ, bool axialQ, bool scalarQ,
                                   bool pseudoscalarQ, bool spatialQ) {
  // Order: V_T, V_E, V_M, V_L, X, A_T, A_E, A_M, A_L, Y, Z, S, P
  const std::array<bool, 13> requested{vectorQ,
                                       vectorQ && spatialQ,
                                       vectorQ && spatialQ,
                                       vectorQ && spatialQ,
                                       vectorQ && spatialQ,
                                       axialQ,
                                       axialQ && spatialQ,
                                       axialQ && spatialQ,
                                       axialQ && spatialQ,
                                       axialQ && spatialQ,
                                       vectorQ && axialQ && spatialQ,
                                       scalarQ,
                                       pseudoscalarQ};
  FormFactorSet K_factors;
  for (std::size_t i = 0; i < K_factors.size(); ++i) {
    if (requested[i]) {
      K_factors[i].resize(E_steps, q_steps, 0.0);
    }
  }
  return K_factors;
}

//==============================================================================
std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>
multipole_operators(const Grid &grid, bool low_q,
                    const SphericalBessel::JL_table *jK_tab, bool vectorQ,
                    bool axialQ, bool scalarQ, bool pseudoscalarQ,
                    bool spatialQ) {
  const auto make =
    [&](char type, char comp,
        bool include) -> std::unique_ptr<DiracOperator::TensorOperator> {
    return include ? DiracOperator::MultipoleOperator(grid, 0, 0.0, type, comp,
                                                      low_q, jK_tab) :
                     nullptr;
  };
  return {make('V', 'T', vectorQ),
          make('V', 'E', vectorQ && spatialQ),
          make('V', 'M', vectorQ && spatialQ),
          make('V', 'L', vectorQ && spatialQ),
          make('A', 'T', axialQ),
          make('A', 'E', axialQ && spatialQ),
          make('A', 'M', axialQ && spatialQ),
          make('A', 'L', axialQ && spatialQ),
          make('S', 'T', scalarQ),
          make('P', 'T', pseudoscalarQ)};
}

//==============================================================================
void accumulate_formFactors(FormFactorSet &K_factors, std::size_t iE,
                            std::size_t iq, double tkp1_x,
                            const ChannelAmplitudes &A) {
  const auto &[t, E, M, L, t5, E5, M5, L5, S, S5] = A;

  const auto add = [&](std::size_t i, double value) {
    if (!K_factors[i].empty()) {
      K_factors[i](iE, iq) += tkp1_x * value;
    }
  };
  // Interference of two amplitudes in the same channel: Re(A A'^*)
  const auto cross = [](std::complex<double> a, std::complex<double> b) {
    return std::real(a * std::conj(b));
  };

  // Vector: temporal, electric, magnetic, longitudinal, and t-L cross term
  add(0, std::norm(t));
  add(1, std::norm(E));
  add(2, std::norm(M));
  add(3, std::norm(L));
  add(4, cross(t, L));
  // Axial (gamma^5) partners
  add(5, std::norm(t5));
  add(6, std::norm(E5));
  add(7, std::norm(M5));
  add(8, std::norm(L5));
  add(9, cross(t5, L5));
  // Vector-axial spatial interference
  add(10, cross(E5, M) - cross(E, M5));
  // Scalar, pseudoscalar
  add(11, std::norm(S));
  add(12, std::norm(S5));
}

//==============================================================================
std::vector<FormFactorSet> calculate_formFactors(
  const HF::HartreeFock *vHF,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  double ec_max, bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, AtomicMethod method,
  double zeff_constant) {

  assert(vHF != nullptr);
  if (diagonal_Eq) {
    assert(qgrid.size() == 1);
  }

  const auto &core = vHF->core();
  const auto n_core = core.size();
  const auto E_steps = Egrid.size();
  const auto q_steps = qgrid.size();

  // H-like (Zeff) approximation: constant Zeff (if given), else "real"
  // Zeff = n*sqrt(-2*en) from the binding energy (same as DarkARC).
  // Bound state used in the matrix elements: the real state, or its H-like
  // (Zeff, pointlike) version, solved analytically or numerically (DiracODE).
  // nb: Zeff state will not have exactly right energy; continuum energies
  // (ec = E + en) and occupation always use the real (HF) orbital.
  std::vector<double> Zeff(n_core, 0.0);
  std::vector<DiracSpinor> bound;
  bound.reserve(n_core);
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    const auto &Fa = core[ia];
    Zeff[ia] = zeff_constant > 0.0 ? zeff_constant : Zeff_real(Fa.en(), Fa.n());
    if (method == AtomicMethod::ZeffAnalytic) {
      bound.push_back(DiracSpinor::exactHlike(
        Fa.n(), Fa.kappa(), Fa.grid_sptr(), Zeff[ia], vHF->alpha()));
    } else if (method == AtomicMethod::Zeff) {
      const auto v_z =
        Nuclear::sphericalNuclearPotential(Zeff[ia], 0.0, vHF->grid().r());
      const auto e0 =
        AtomData::diracen(Zeff[ia], Fa.n(), Fa.kappa(), vHF->alpha());
      bound.push_back(DiracODE::boundState(
        Fa.n(), Fa.kappa(), e0, Fa.grid_sptr(), v_z, {}, vHF->alpha()));
    } else {
      bound.push_back(Fa);
    }
  }

  // Output: the factors of each core orbital
  std::vector<FormFactorSet> K_nk(
    n_core, allocate_formFactors(E_steps, q_steps, vectorQ, axialQ, scalarQ,
                                 pseudoscalarQ, spatialQ));

  // Work list: every (orbital, E) with the ejected electron energy within the
  // limits. Flattened for load balance: a deep shell is ionised at only the
  // few highest energies.
  std::vector<std::pair<std::size_t, std::size_t>> tasks;
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    for (std::size_t iE = 0; iE < E_steps; ++iE) {
      const auto ec = Egrid[iE] + core[ia].en();
      if (ec > ec_min && ec <= ec_max) {
        tasks.emplace_back(ia, iE);
      }
    }
  }

  // One operator set per thread: the rank and frequency are set in place
  // inside the loop. Operators not requested are null.
  const auto prototypes =
    multipole_operators(vHF->grid(), low_q, &jK_tab, vectorQ, axialQ, scalarQ,
                        pseudoscalarQ, spatialQ);
  const auto max_threads = std::size_t(omp_get_max_threads());
  std::vector<std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>>
    thread_operators(max_threads);
  for (auto &operators : thread_operators) {
    for (std::size_t i = 0; i < operators.size(); ++i) {
      operators[i] = prototypes[i] ? prototypes[i]->clone() : nullptr;
    }
  }

  // Each task owns one row (iE) of one orbital's factors: no reduction
  qip::ProgressBar bar(tasks.size());
#pragma omp parallel for schedule(dynamic)
  for (std::size_t it = 0; it < tasks.size(); ++it) {
    const auto ia = tasks[it].first;
    const auto iE = tasks[it].second;
    const auto &Fa = core[ia];
    const auto &Fa_t = bound[ia];
    const auto ec = Egrid[iE] + Fa.en();
    // Continuum l reached by multipoles up to Kmax (either parity):
    // j_e = j_a -/+ Kmax, l_e = j_e -/+ 1/2; within the lc_minmax limits
    auto lc_min = std::max((Fa.twoj() - 2 * Kmax - 1) / 2, 0);
    auto lc_max = (Fa.twoj() + 2 * Kmax + 1) / 2;
    if (lc_minmax) {
      lc_min = std::max(lc_min, lc_minmax->at(0));
      lc_max = std::min(lc_max, lc_minmax->at(1));
    }

    ContinuumOrbitals cntm(vHF);
    if (method == AtomicMethod::Zeff) {
      cntm.solveContinuumZeff(ec, lc_min, lc_max, Zeff[ia], &Fa_t,
                              force_orthog);
    } else if (method == AtomicMethod::ZeffAnalytic) {
      cntm.solveContinuumZeffAnalytic(ec, lc_min, lc_max, Zeff[ia], &Fa_t,
                                      force_orthog);
    } else {
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa_t, force_rescale,
                            hole_particle, force_orthog);
    }

    auto &operators = thread_operators[std::size_t(omp_get_thread_num())];
    for (int k = Kmin; k <= Kmax; ++k) {
      const auto tkp1_x = (2.0 * k + 1.0) * Fa.occ_frac();
      for (std::size_t iq = 0; iq < q_steps; ++iq) {
        // The operators take qc as their "frequency" (or E itself in the
        // diagonal, massless-absorption, case). Rank first (sets the
        // parity), then frequency (fills the Bessel vectors)
        const auto qc = diagonal_Eq ? Egrid[iE] : qgrid[iq] * PhysConst::c;
        for (auto &h : operators) {
          if (h) {
            h->updateRank(k);
            h->updateFrequency(qc);
          }
        }
        for (const auto &Fe : cntm.orbitals) {
          ChannelAmplitudes A{};
          for (std::size_t i = 0; i < operators.size(); ++i) {
            if (operators[i]) {
              A[i] = operators[i]->reducedME(Fe, Fa_t);
            }
          }
          accumulate_formFactors(K_nk[ia], iE, iq, tkp1_x, A);
        }
      }
    }
    bar.update();
  }

  return K_nk;
}

//==============================================================================
FormFactorsRPA calculate_formFactors_RPA(
  const HF::HartreeFock *vHF,
  const std::optional<std::array<int, 2>> &lc_minmax, double ec_min,
  double ec_max, bool force_rescale, bool hole_particle, bool force_orthog,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, const RPAOptions &rpa_options) {

  assert(vHF != nullptr);
  if (diagonal_Eq) {
    assert(qgrid.size() == 1);
  }

  const auto &core = vHF->core();
  const auto n_core = core.size();
  const auto E_steps = Egrid.size();
  const auto q_steps = qgrid.size();

  FormFactorsRPA out;
  out.bare.assign(n_core,
                  allocate_formFactors(E_steps, q_steps, vectorQ, axialQ,
                                       scalarQ, pseudoscalarQ, spatialQ));
  out.rpa = out.bare;
  out.eps.resize(E_steps, q_steps, 0.0);

  // The requested operators (the others are null); per-thread clones are
  // made per (E, K, operator) block below
  const auto prototypes =
    multipole_operators(vHF->grid(), low_q, &jK_tab, vectorQ, axialQ, scalarQ,
                        pseudoscalarQ, spatialQ);
  std::vector<std::size_t> active_operators;
  for (std::size_t i = 0; i < prototypes.size(); ++i) {
    if (prototypes[i]) {
      active_operators.push_back(i);
    }
  }
  if (active_operators.empty()) {
    return out;
  }

  // The operators take qc as their "frequency" (or E itself in the diagonal,
  // massless-absorption, case)
  const auto qc_at = [&](std::size_t iE, std::size_t iq) {
    return diagonal_Eq ? Egrid[iE] : qgrid[iq] * PhysConst::c;
  };

  // Enough q points to keep every thread busy: parallel over q, each thread
  // owning its solver (the solver's own parallel regions are then nested,
  // hence inert). Otherwise q runs serially, with the parallelism inside
  // the solver.
  const bool parallel_q = q_steps >= std::size_t(omp_get_max_threads());

  // A first-order solve (max_its <= 1) is what was asked for: its eps is
  // the size of the correction, not a convergence measure, and is not tested
  const bool test_convergence = rpa_options.max_its > 1;
  const auto converged = [=](double eps) {
    return !test_convergence ||
           (!std::isnan(eps) && eps < rpa_options.eps_fail);
  };
  // The worse of two eps values (nan is the worst)
  const auto worst_eps = [](double a, double b) {
    return (std::isnan(a) || std::isnan(b)) ? std::nan("") : std::max(a, b);
  };

  const auto n_K = std::size_t(Kmax - Kmin + 1);
  const auto solves_per_E = n_K * active_operators.size() * q_steps;
  fmt::print("RPA: {} solves per energy (K x operators x q = {} x {} x {}); "
             " parallel over {}\n",
             solves_per_E, n_K, active_operators.size(), q_steps,
             rpa_options.eps, rpa_options.eps_fail,
             parallel_q ? "q" : "RPA channels");

  for (std::size_t iE = 0; iE < E_steps; ++iE) {
    const auto omega = Egrid[iE];

    // List of ionised (energetically accessible) orbitals
    std::vector<std::size_t> ionised;
    for (std::size_t ia = 0; ia < n_core; ++ia) {
      const auto ec = omega + core[ia].en();
      if (ec > ec_min && ec <= ec_max) {
        ionised.push_back(ia);
      }
    }
    if (ionised.empty()) {
      continue;
    }

    // Continuum states of the ejected electron of each
    std::vector<ContinuumOrbitals> cntm(ionised.size(), ContinuumOrbitals(vHF));
#pragma omp parallel for schedule(dynamic)
    for (std::size_t j = 0; j < ionised.size(); ++j) {
      const auto &Fa = core[ionised[j]];
      // Continuum l reached by multipoles up to Kmax (either parity):
      // j_e = j_a -/+ Kmax, l_e = j_e -/+ 1/2; within the lc_minmax limits
      auto lc_min = std::max((Fa.twoj() - 2 * Kmax - 1) / 2, 0);
      auto lc_max = (Fa.twoj() + 2 * Kmax + 1) / 2;
      if (lc_minmax) {
        lc_min = std::max(lc_min, lc_minmax->at(0));
        lc_max = std::min(lc_max, lc_minmax->at(1));
      }
      cntm[j].solveContinuumHF(omega + Fa.en(), lc_min, lc_max, &Fa,
                               force_rescale, hole_particle, force_orthog);
    }

    // Every (hole orbital, ejected state) channel at this energy
    std::vector<std::pair<std::size_t, const DiracSpinor *>> channels;
    for (std::size_t j = 0; j < ionised.size(); ++j) {
      for (const auto &Fe : cntm[j].orbitals) {
        channels.emplace_back(ionised[j], &Fe);
      }
    }

    fmt::print("E = {:.6g} eV: {} orbital(s) ionised, {} channels\n",
               omega * PhysConst::Hartree_eV, ionised.size(), channels.size());
    std::cout << std::flush;

    qip::ProgressBar bar(solves_per_E);
    for (int k = Kmin; k <= Kmax; ++k) {
      // Amplitudes of every channel at every q, [channel][iq], filled one
      // operator at a time, then accumulated
      std::vector<std::vector<ChannelAmplitudes>> A_bare(
        channels.size(), std::vector<ChannelAmplitudes>(q_steps));
      auto A_rpa = A_bare;

      for (const auto i_op : active_operators) {
#pragma omp parallel if (parallel_q)
        {
          // Per thread: its own operator (rank set here, frequency per q)
          // and solver. Static schedule: consecutive q on one thread, so
          // each solve warm starts from the neighbouring q
          auto h = prototypes[i_op]->clone();
          h->updateRank(k);
          h->updateFrequency(qc_at(iE, 0));
          ExternalField::TDHFcntm rpa(h.get(), vHF);
          rpa.eps_target() = rpa_options.eps;
#pragma omp for schedule(static)
          for (std::size_t iq = 0; iq < q_steps; ++iq) {
            h->updateFrequency(qc_at(iE, iq));

            rpa.solve_core(omega, rpa_options.max_its, false);
            const auto eps_solve = rpa.last_eps();
            const bool use_rpa = converged(eps_solve);
            if (!use_rpa) {
              // prevents contaminating next q iteration
              rpa.clear();
            }

            out.eps(iE, iq) = worst_eps(out.eps(iE, iq), eps_solve);

            for (std::size_t ic = 0; ic < channels.size(); ++ic) {
              const auto &Fa = core[channels[ic].first];
              const auto &Fe = *channels[ic].second;
              if (h->isZero(Fe, Fa))
                continue;
              const auto t0 = h->reducedME(Fe, Fa);
              A_bare[ic][iq][i_op] = t0;
              A_rpa[ic][iq][i_op] = use_rpa ? t0 + rpa.dV_complex(Fe, Fa) : t0;
            }
            bar.update();
          }
        }
      }

      for (std::size_t ic = 0; ic < channels.size(); ++ic) {
        const auto ia = channels[ic].first;
        const auto tkp1_x = (2.0 * k + 1.0) * core[ia].occ_frac();
        for (std::size_t iq = 0; iq < q_steps; ++iq) {
          accumulate_formFactors(out.bare[ia], iE, iq, tkp1_x, A_bare[ic][iq]);
          accumulate_formFactors(out.rpa[ia], iE, iq, tkp1_x, A_rpa[ic][iq]);
        }
      }
    }
  }

  return out;
}

//==============================================================================
std::pair<std::size_t, std::size_t>
interpolate_failed_rpa(FormFactorsRPA &K_rpa, double eps_fail) {

  const auto E_steps = K_rpa.eps.rows();
  const auto q_steps = K_rpa.eps.cols();

  const auto failed = [&](std::size_t iE, std::size_t iq) {
    const auto eps = K_rpa.eps(iE, iq);
    return std::isnan(eps) || eps > eps_fail;
  };

  // Count number that failed
  std::size_t n_failed = 0;
  for (std::size_t iE = 0; iE < E_steps; ++iE) {
    for (std::size_t iq = 0; iq < q_steps; ++iq) {
      if (failed(iE, iq)) {
        ++n_failed;
      }
    }
  }
  // No neighbour to interpolate from (a diagonal E-q calculation, say)
  if (q_steps < 2) {
    return {n_failed, 0};
  }

  std::size_t n_corrected = 0;
  for (std::size_t iE = 0; iE < E_steps; ++iE) {
    for (std::size_t iq = 0; iq < q_steps; ++iq) {
      if (!failed(iE, iq))
        continue;
      const bool have_below = iq > 0 && !failed(iE, iq - 1);
      const bool have_above = iq + 1 < q_steps && !failed(iE, iq + 1);
      if (!have_below && !have_above)
        continue;
      for (std::size_t ia = 0; ia < K_rpa.rpa.size(); ++ia) {
        for (std::size_t i = 0; i < K_rpa.rpa[ia].size(); ++i) {
          auto &K_factor = K_rpa.rpa[ia][i];
          const auto &K_bare = K_rpa.bare[ia][i];
          if (K_factor.empty())
            continue;
          // Relative shift of each converged neighbour in q. A zero bare
          // factor has no shift to speak of
          double sum_shift = 0.0;
          int n_sides = 0;
          if (have_below && K_bare(iE, iq - 1) != 0.0) {
            sum_shift += K_factor(iE, iq - 1) / K_bare(iE, iq - 1);
            ++n_sides;
          }
          if (have_above && K_bare(iE, iq + 1) != 0.0) {
            sum_shift += K_factor(iE, iq + 1) / K_bare(iE, iq + 1);
            ++n_sides;
          }
          if (n_sides == 0)
            continue;
          // Both sides: their mean. One side only: half its correction,
          // since the shift is unconstrained on the other side
          const auto shift =
            n_sides == 2 ? 0.5 * sum_shift : 1.0 + 0.5 * (sum_shift - 1.0);
          K_factor(iE, iq) = shift * K_bare(iE, iq);
        }
      }
      ++n_corrected;
    }
  }
  return {n_failed, n_corrected};
}

//==============================================================================
bool check_radial_grid(double Emax_au, double qmax_au, const Grid &rgrid,
                       double alpha) {
  bool ok = true;

  // Check grid type: only loglinear is reasonable for this module
  if (rgrid.type() != GridType::loglinear) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning:\n");
    std::cout << "This module unlikely to work with grid type: "
              << GridParameters::parseType(rgrid.type())
              << "; consider changing to loglinear\n";
  }

  // *Very* rough estimate of good aximum q range.
  const auto r_q = 1.0;
  std::cout << "\n"
            << "Check grid for q_max:\n";
  const auto i = rgrid.getIndex(r_q);
  const auto dr_q = rgrid.drdu(i) * rgrid.du();
  const auto qmax_targ = 2.0 * M_PI / (3.0 * dr_q);
  fmt::print("Very rough guess for maximum safe q: {:.0f} au = {:.2f} MeV\n",
             qmax_targ, qmax_targ * UnitConv::Momentum_au_to_MeV);
  if (qmax_au > qmax_targ) {
    // dr required at r_q for qmax (same formula, inverted):
    const auto dr_targ = 2.0 * M_PI / (8.0 * qmax_au);
    // num_points required (same r0, rmax, b):
    const auto n_q = Grid::calc_num_points_from_du(
      rgrid.r0(), rgrid.rmax(), dr_targ / rgrid.drdu(i), rgrid.type(),
      rgrid.loglin_b());
    // smallest sufficient b (same num_points): dr(r_q) shrinks as b grows
    // (larger b = denser at low r) - opposite direction to the E case.
    // Negative means no b is sufficient (need more points):
    const double L = std::log(rgrid.rmax() / rgrid.r0());
    const double tn = dr_targ * double(rgrid.num_points() - 1);
    const double b_q =
      r_q * (tn - (rgrid.rmax() - rgrid.r0())) / (L * r_q - tn);
    fmt2::styled_print(fg(fmt::color::orange), "Warning: ");
    fmt::print("Grid may not be dense enough for q = {:.0f} au = {:.2f} MeV\n",
               qmax_au, qmax_au * UnitConv::Momentum_au_to_MeV);
    if (b_q > 0.0) {
      fmt::print("Try:\n - increasing num_points to {}; or\n - increasing b "
                 "to {:.1f} (denser at low r)\n",
                 n_q, b_q);
    } else {
      fmt::print("Try increasing num_points to {}\n", n_q);
    }
    ok = false;
    std::cout << "(nb: Rough: high q might not contribute, in which case "
                 "numerical error there might not matter. Always check)\n";
  }

  // E check: grid must resolve the continuum oscillations at Emax
  // (RequiredContinuumGrid: 15 points-per-wavelength at the coarsest
  // point, as required to store the state pointwise to rmax):
  const auto req = DiracODE::RequiredContinuumGrid(Emax_au, rgrid, 10.0, alpha);

  fmt::print("\nGrid check for continuum states: E_max = {:.3g} au:\n"
             "Require num_points ~ {} (same r0, rmax, b); have {}\n",
             Emax_au, req.num_points, rgrid.num_points());

  if (rgrid.num_points() < req.num_points) {
    fmt2::styled_print(fg(fmt::color::red), "Warning: ");
    fmt::print("Grid may not be dense enough for continuum states with "
               "E = {:.2f} au\n",
               Emax_au);
    fmt::print(
      "Try:\n - increasing num_points to {}; or\n - reduce b to {:.2f} and "
      "increase num_points to {}\n",
      req.num_points, req.b, req.num_points_b);
    std::cout << "Program will continue; results may be inaccurate\n";
    ok = false;
  }

  return ok;
}

//==============================================================================
void write_to_file_xyz(const std::string &filename,
                       const std::vector<double> &E_grid,
                       const std::vector<double> &q_grid,
                       const std::vector<std::string> &titles,
                       const std::vector<std::string> &descriptions,
                       std::vector<LinAlg::Matrix_view<const double>> factors,
                       Units units, int num_digits, bool diagonal) {

  assert(titles.size() == factors.size() && "Each factor must have a title");

  std::ofstream out_file(filename);

  // Just write to screen:
  fmt::print("Writing factors to file: {}\n", filename);
  fmt::print("Using {} units for E/q (K dimensionless)\n",
             units == Units::Atomic ? "atomic" : "eV");

  // check array sizes:
  for (const auto &K : factors) {
    if (K.size() == 0)
      continue;
    assert(K.rows() == E_grid.size() && "Factors: rows must match energy grid");
    assert(K.cols() == q_grid.size() &&
           "Factors: cols must match momentum grid");
  }

  // Require num_digits be between 3 and 16:
  num_digits = std::clamp(num_digits, 3, 16);

  const auto unit_E = units == Units::Atomic ? 1.0 : UnitConv::Energy_au_to_eV;
  const auto unit_q =
    units == Units::Atomic ? 1.0 : UnitConv::Momentum_au_to_eV;
  const auto unit_str = units == Units::Particle ? "eV" : "au";

  out_file << "# ampsci Kion form factors output data file: " << filename
           << " - xyz format\n";
  fmt::print(out_file, "# Units: ");
  if (units == Units::Atomic) {
    fmt::print(out_file, "Atomic units. [q] = [1/a0], [E] = [E_H], [K] = 1\n");
  } else if (units == Units::Particle) {
    fmt::print(out_file, "Particle units. [q] = eV, [E] = eV, [K] = 1\n");
  } else {
    std::cout << "units error\n";
  }
  fmt::print(out_file, "# nb: E_H = m_e (c*α)^2 = ~27.21 eV\n");
  fmt::print(out_file,
             "# nb: 1/a0 = m_e*c*α/hbar = E_H / (c*α*hbar) = ~3729 eV\n");

  // Map short column headers to useful descriptions
  out_file << "# Columns:\n";
  fmt::print(out_file, "# {:<4} : {} in {}\n", "E", "Energy exchange",
             unit_str);
  fmt::print(out_file, "# {:<4} : {} in {}\n", "q", "Momentum transfer",
             unit_str);
  for (std::size_t i = 0; i < factors.size(); ++i) {
    if (factors.at(i).size() == 0)
      continue;
    fmt::print(out_file, "# {:<4} : {}\n", titles[i], descriptions[i]);
  }
  out_file << "################################################################"
              "################\n";

  // Add titles (column headers)
  const auto width = 7 + num_digits;
  fmt::print(out_file, "{:<{}}  ", "E", width);
  fmt::print(out_file, "{:<{}} ", "q", width);
  for (std::size_t i = 0; i < factors.size(); ++i) {
    if (factors.at(i).size() == 0)
      continue;
    fmt::print(out_file, "{:<{}} ", titles[i], width);
  }
  out_file << "\n";

  // fmt::print(out_file, "{:<{}} {:<{}.{}f}", str, n, value, n, m);
  for (std::size_t iE = 0; iE < E_grid.size(); ++iE) {
    const auto E = E_grid.at(iE) * unit_E;
    for (std::size_t iq = 0; iq < q_grid.size(); ++iq) {
      const auto q = diagonal ? E_grid.at(iE) * PhysConst::alpha * unit_q :
                                q_grid.at(iq) * unit_q;
      fmt::print(out_file, "{:+{}.{}e} ", E, width, num_digits);
      fmt::print(out_file, "{:+{}.{}e} ", q, width, num_digits);
      for (const auto &K : factors) {
        if (K.size() == 0)
          continue;
        fmt::print(out_file, "{:+{}.{}e} ", K(iE, iq), width, num_digits);
      }
      out_file << "\n";
    }
    // Print new line between each new energy (so long as >1 q).
    // This makes gnuplot happy, and doesn't impact pyplot
    if (q_grid.size() > 1 && iE + 1 < E_grid.size()) {
      out_file << "\n";
    }
  }
}

//==============================================================================
void write_to_file_xyz_13(const std::string &filename,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::vector<std::string> &titles,
                          const std::vector<std::string> &descriptions,
                          const FormFactorSet &K_factors, Units units,
                          int num_digits, bool diagonal) {

  const auto &K_VT = K_factors[0];  // Vector: temporal
  const auto &K_VE = K_factors[1];  // Vector: electric
  const auto &K_VM = K_factors[2];  // Vector: magnetic
  const auto &K_VL = K_factors[3];  // Vector: longitudinal
  const auto &K_X = K_factors[4];   // Vector: v-v cross
  const auto &K_T5 = K_factors[5];  // Axial: temporal
  const auto &K_E5 = K_factors[6];  // Axial: electric
  const auto &K_M5 = K_factors[7];  // Axial: magnetic
  const auto &K_L5 = K_factors[8];  // Axial: longitudinal
  const auto &K_X5 = K_factors[9];  // Axial: a-a cross
  const auto &K_Z = K_factors[10];  // Vector-Axial interference
  const auto &K_S = K_factors[11];  // Scalar
  const auto &K_S5 = K_factors[12]; // Pseudo-scalar

  std::vector<LinAlg::Matrix_view<const double>> factors{
    K_VT, K_VE, K_VM, K_VL, K_X, K_T5, K_E5, K_M5, K_L5, K_X5, K_Z, K_S, K_S5};

  return write_to_file_xyz(filename, E_grid, q_grid, titles, descriptions,
                           factors, units, num_digits, diagonal);
}

//==============================================================================
void write_to_file_matrix(const LinAlg::Matrix<double> &K,
                          const std::vector<double> &E_grid,
                          const std::vector<double> &q_grid,
                          const std::string &filename, int num_digits,
                          Units units) {
  // optional format argument?
  assert(K.rows() == E_grid.size());
  assert(K.cols() == q_grid.size());
  std::ofstream out_file(filename);

  const auto unit_E = units == Units::Atomic ? 1.0 : UnitConv::Energy_au_to_eV;
  const auto unit_q =
    units == Units::Atomic ? 1.0 : UnitConv::Momentum_au_to_eV;

  out_file << "# Kion output data file: " << filename << " - matrix format\n";

  if (units == Units::Atomic) {
    fmt::print(out_file,
               "# Atomic units. [q] = [1/a0], [E] = [E_H], [K] = 1\n");
  } else if (units == Units::Particle) {
    fmt::print(out_file, "# Particle units. [q] = eV, [E] = eV, [K] = 1\n");
  } else {
    std::cout << "units error\n";
  }
  fmt::print(out_file, "# E_H = m_e (c*α)^2 = ~27.21 eV\n");
  fmt::print(out_file, "# 1/a0 = m_e*c*α/hbar = E_H / (c*α*hbar) = ~3729 eV\n");

  out_file << "\n# E values:\n";
  for (auto E : E_grid) {
    fmt::print(out_file, "{:+.{}e} ", E * unit_E, num_digits);
  }
  out_file << "\n\n# q values:\n";
  for (auto q : q_grid) {
    fmt::print(out_file, "{:+.{}e} ", q * unit_q, num_digits);
  }
  out_file << "\n\n";

  out_file << "# K values K(E,q). Each new row is new E, each col is new q\n";
  for (std::size_t iE = 0; iE < E_grid.size(); ++iE) {
    for (std::size_t iq = 0; iq < q_grid.size(); ++iq) {
      fmt::print(out_file, "{:+.{}e} ", K(iE, iq), num_digits);
    }
    out_file << '\n';
  }
}

} // namespace Kion