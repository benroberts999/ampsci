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
#include "qip/Vector.hpp"
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
void accumulate_formFactors(FormFactorSet *K_factors, std::size_t iE,
                            std::size_t iq, double tkp1_x,
                            const ChannelAmplitudes &A) {
  assert(K_factors != nullptr);
  const auto &[t, E, M, L, t5, E5, M5, L5, S, S5] = A;

  const auto add = [&](std::size_t i, double value) {
    auto &K_factor = (*K_factors)[i];
    if (!K_factor.empty()) {
      K_factor(iE, iq) += tkp1_x * value;
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
  std::vector<DiracSpinor> bound_states;
  bound_states.reserve(n_core);
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    const auto &Fa_hf = core[ia];
    Zeff[ia] =
      zeff_constant > 0.0 ? zeff_constant : Zeff_real(Fa_hf.en(), Fa_hf.n());
    if (method == AtomicMethod::ZeffAnalytic) {
      bound_states.push_back(DiracSpinor::exactHlike(
        Fa_hf.n(), Fa_hf.kappa(), Fa_hf.grid_sptr(), Zeff[ia], vHF->alpha()));
    } else if (method == AtomicMethod::Zeff) {
      const auto v_z =
        Nuclear::sphericalNuclearPotential(Zeff[ia], 0.0, vHF->grid().r());
      const auto e0 =
        AtomData::diracen(Zeff[ia], Fa_hf.n(), Fa_hf.kappa(), vHF->alpha());
      bound_states.push_back(DiracODE::boundState(Fa_hf.n(), Fa_hf.kappa(), e0,
                                                  Fa_hf.grid_sptr(), v_z, {},
                                                  vHF->alpha()));
    } else {
      bound_states.push_back(Fa_hf);
    }
  }

  // Output: the factors of each core orbital
  std::vector<FormFactorSet> K_nk(
    n_core, allocate_formFactors(E_steps, q_steps, vectorQ, axialQ, scalarQ,
                                 pseudoscalarQ, spatialQ));

  // Work list: every (orbital, E) with the ejected electron energy within the
  // limits. Flattened for load balance: a deep shell is ionised at only the
  // few highest energies.
  std::vector<std::pair<std::size_t, std::size_t>> orbital_E_pairs;
  for (std::size_t ia = 0; ia < n_core; ++ia) {
    for (std::size_t iE = 0; iE < E_steps; ++iE) {
      const auto ec = Egrid[iE] + core[ia].en();
      if (ec > ec_min && ec <= ec_max) {
        orbital_E_pairs.emplace_back(ia, iE);
      }
    }
  }

  // One operator set per thread: the rank and frequency are set in place
  // inside the loop. Operators not requested are null.
  const auto multipoles =
    multipole_operators(vHF->grid(), low_q, &jK_tab, vectorQ, axialQ, scalarQ,
                        pseudoscalarQ, spatialQ);
  const auto max_threads = std::size_t(omp_get_max_threads());
  std::vector<std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>>
    thread_multipoles(max_threads);
  for (auto &clones : thread_multipoles) {
    for (std::size_t i = 0; i < clones.size(); ++i) {
      clones[i] = multipoles[i] ? multipoles[i]->clone() : nullptr;
    }
  }

  // Each pair owns one row (iE) of one orbital's factors: no reduction
  qip::ProgressBar bar(orbital_E_pairs.size());
#pragma omp parallel for schedule(dynamic)
  for (std::size_t it = 0; it < orbital_E_pairs.size(); ++it) {
    const auto [ia, iE] = orbital_E_pairs[it];
    // Fa is the state in the matrix elements; Fa_hf sets the energies
    // and occupation
    const auto &Fa_hf = core[ia];
    const auto &Fa = bound_states[ia];
    const auto ec = Egrid[iE] + Fa_hf.en();
    // Continuum l reached by multipoles up to Kmax (either parity):
    // j_e = j_a -/+ Kmax, l_e = j_e -/+ 1/2; within the lc_minmax limits
    auto lc_min = std::max((Fa_hf.twoj() - 2 * Kmax - 1) / 2, 0);
    auto lc_max = (Fa_hf.twoj() + 2 * Kmax + 1) / 2;
    if (lc_minmax) {
      lc_min = std::max(lc_min, lc_minmax->at(0));
      lc_max = std::min(lc_max, lc_minmax->at(1));
    }

    ContinuumOrbitals cntm(vHF);
    if (method == AtomicMethod::Zeff) {
      cntm.solveContinuumZeff(ec, lc_min, lc_max, Zeff[ia], &Fa, force_orthog);
    } else if (method == AtomicMethod::ZeffAnalytic) {
      cntm.solveContinuumZeffAnalytic(ec, lc_min, lc_max, Zeff[ia], &Fa,
                                      force_orthog);
    } else {
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa, force_rescale,
                            hole_particle, force_orthog);
    }

    auto &own_multipoles = thread_multipoles[std::size_t(omp_get_thread_num())];
    for (int k = Kmin; k <= Kmax; ++k) {
      const auto tkp1_x = (2.0 * k + 1.0) * Fa_hf.occ_frac();
      for (std::size_t iq = 0; iq < q_steps; ++iq) {
        // The operators take qc as their "frequency" (or E itself in the
        // diagonal, massless-absorption, case). Rank first (sets the
        // parity), then frequency (fills the Bessel vectors)
        const auto qc = diagonal_Eq ? Egrid[iE] : qgrid[iq] * PhysConst::c;
        for (auto &h : own_multipoles) {
          if (h) {
            h->updateRank(k);
            h->updateFrequency(qc);
          }
        }
        for (const auto &Fe : cntm.orbitals) {
          ChannelAmplitudes A{};
          for (std::size_t i = 0; i < own_multipoles.size(); ++i) {
            if (own_multipoles[i]) {
              A[i] = own_multipoles[i]->reducedME(Fe, Fa);
            }
          }
          accumulate_formFactors(&K_nk[ia], iE, iq, tkp1_x, A);
        }
      }
    }
    bar.update();
  }

  return K_nk;
}

//==============================================================================
std::pair<int, int>
continuum_l_range(const DiracSpinor &Fa, int Kmax,
                  const std::optional<std::array<int, 2>> &lc_minmax) {
  // j_e = j_a -/+ Kmax, l_e = j_e -/+ 1/2 (either parity)
  auto lc_min = std::max((Fa.twoj() - 2 * Kmax - 1) / 2, 0);
  auto lc_max = (Fa.twoj() + 2 * Kmax + 1) / 2;
  if (lc_minmax) {
    lc_min = std::max(lc_min, lc_minmax->at(0));
    lc_max = std::min(lc_max, lc_minmax->at(1));
  }
  return {lc_min, lc_max};
}

//==============================================================================
std::vector<IonisedOrbital> solve_ionised_orbitals_at_omega(
  const HF::HartreeFock *vHF, double omega, double ec_min, double ec_max,
  int Kmax, const std::optional<std::array<int, 2>> &lc_minmax,
  bool force_rescale, bool hole_particle, bool force_orthog) {
  assert(vHF != nullptr);
  const auto &core = vHF->core();

  // The orbitals whose ejected electron energy lies within the limits
  std::vector<IonisedOrbital> ionised;
  for (std::size_t ia = 0; ia < core.size(); ++ia) {
    const auto ec = omega + core[ia].en();
    if (ec > ec_min && ec <= ec_max) {
      ionised.push_back({ia, ContinuumOrbitals(vHF)});
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t j = 0; j < ionised.size(); ++j) {
    const auto &Fa = core[ionised[j].core_index];
    const auto [lc_min, lc_max] = continuum_l_range(Fa, Kmax, lc_minmax);
    ionised[j].ejected.solveContinuumHF(omega + Fa.en(), lc_min, lc_max, &Fa,
                                        force_rescale, hole_particle,
                                        force_orthog);
  }
  return ionised;
}

//==============================================================================
std::vector<IonisationChannel>
construct_channels(const std::vector<DiracSpinor> &core,
                   const std::vector<IonisedOrbital> &ionised) {
  std::vector<IonisationChannel> channels;
  for (const auto &orbital : ionised) {
    for (const auto &Fe : orbital.ejected.orbitals) {
      channels.push_back({orbital.core_index, &core[orbital.core_index], &Fe});
    }
  }
  return channels;
}

//==============================================================================
std::vector<std::size_t> active_operators(
  const std::array<std::unique_ptr<DiracOperator::TensorOperator>, 10>
    &operators) {
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < operators.size(); ++i) {
    if (operators[i]) {
      active.push_back(i);
    }
  }
  return active;
}

//==============================================================================
double solve_channel_amplitudes(const DiracOperator::TensorOperator &h,
                                std::size_t i_op, ExternalField::TDHFcntm *rpa,
                                double omega, const RPAOptions &rpa_options,
                                const std::vector<IonisationChannel> &channels,
                                std::vector<ChannelAmplitudes> *A_bare,
                                std::vector<ChannelAmplitudes> *A_rpa) {
  assert(rpa != nullptr && A_bare != nullptr && A_rpa != nullptr);
  assert(A_bare->size() == channels.size() && A_rpa->size() == channels.size());

  rpa->solve_core(omega, rpa_options.max_its, false);
  const auto eps = rpa->last_eps();
  // A first-order solve (max_its <= 1) is what was asked for: its eps is the
  // size of the correction, not a convergence measure
  const bool use_rpa = rpa_options.max_its <= 1 ||
                       (!std::isnan(eps) && eps < rpa_options.eps_fail);
  // A failed solve is not used, and is not warm started from
  if (!use_rpa) {
    rpa->clear();
  }

  for (std::size_t ic = 0; ic < channels.size(); ++ic) {
    const auto &Fa = *channels[ic].hole;
    const auto &Fe = *channels[ic].ejected;
    if (h.isZero(Fe, Fa))
      continue;
    const auto bare = h.reducedME(Fe, Fa);
    (*A_bare)[ic][i_op] = bare;
    (*A_rpa)[ic][i_op] = use_rpa ? bare + rpa->dV_complex(Fe, Fa) : bare;
  }
  return eps;
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
  using namespace qip::overloads;

  assert(vHF != nullptr);
  if (diagonal_Eq) {
    assert(qgrid.size() == 1);
  }

  const auto &core = vHF->core();
  const auto n_core = core.size();
  const auto E_steps = Egrid.size();
  const auto q_steps = qgrid.size();

  FormFactorsRPA factors;
  factors.bare.assign(n_core,
                      allocate_formFactors(E_steps, q_steps, vectorQ, axialQ,
                                           scalarQ, pseudoscalarQ, spatialQ));
  factors.rpa = factors.bare;
  factors.eps.resize(E_steps, q_steps, 0.0);

  // The requested operators (the others are null); per-thread clones are
  // made per (E, K, operator) block below
  const auto multipoles =
    multipole_operators(vHF->grid(), low_q, &jK_tab, vectorQ, axialQ, scalarQ,
                        pseudoscalarQ, spatialQ);
  const auto active_multipoles = active_operators(multipoles);
  if (active_multipoles.empty()) {
    return factors;
  }

  // The energies at which some orbital is ionised
  // highest first (better OMP load balance)
  std::vector<std::size_t> active_E;
  for (std::size_t i = 0; i < E_steps; ++i) {
    const auto iE = E_steps - 1 - i;
    for (const auto &Fa : core) {
      const auto ec = Egrid[iE] + Fa.en();
      if (ec > ec_min && ec <= ec_max) {
        active_E.push_back(iE);
        break;
      }
    }
  }

  // Decide which (E or q) to parallelise over
  // Neither: is parallelised over RPA instead
  const auto n_threads = std::size_t(omp_get_max_threads());
  const bool parallel_E = n_threads > 1 && active_E.size() >= n_threads;
  const bool parallel_q = !parallel_E && q_steps >= n_threads;

  const auto num_Ks = std::size_t(Kmax - Kmin + 1);
  const auto solves_per_E = num_Ks * active_multipoles.size() * q_steps;
  fmt::print("RPA: {} energies with an ionised orbital, {} solves per energy "
             "(K x operators x q = {} x {} x {}); parallel over {}\n",
             active_E.size(), solves_per_E, num_Ks, active_multipoles.size(),
             q_steps,
             parallel_E ? "E" :
             parallel_q ? "q" :
                          "RPA channels");

  // Progress: one bar over the run when E is parallel (energies then finish
  // out of order), else one per energy
  qip::ProgressBar run_bar(active_E.size() * solves_per_E, parallel_E);

#pragma omp parallel for schedule(dynamic) if (parallel_E)
  for (std::size_t i = 0; i < active_E.size(); ++i) {
    const auto iE = active_E[i];
    const auto omega = Egrid[iE];

    const auto ionised = solve_ionised_orbitals_at_omega(
      vHF, omega, ec_min, ec_max, Kmax, lc_minmax, force_rescale, hole_particle,
      force_orthog);
    assert(!ionised.empty());
    const auto channels = construct_channels(core, ionised);
    // The operators take qc as their "frequency" (omega itself in the
    // diagonal, massless-absorption, case)
    const auto frequencies =
      diagonal_Eq ? std::vector<double>(q_steps, omega) : qgrid * PhysConst::c;

    if (!parallel_E) {
      fmt::print("E = {:.6g} eV: {} orbital(s) ionised, {} channels\n",
                 omega * PhysConst::Hartree_eV, ionised.size(),
                 channels.size());
      std::cout << std::flush;
    }
    qip::ProgressBar E_bar(solves_per_E, !parallel_E);
    auto &bar = parallel_E ? run_bar : E_bar;
    for (int k = Kmin; k <= Kmax; ++k) {
      // Amplitudes of every channel at every q, [iq][channel], filled one
      // operator at a time, then accumulated
      std::vector<std::vector<ChannelAmplitudes>> A_bare(
        q_steps, std::vector<ChannelAmplitudes>(channels.size()));
      auto A_rpa = A_bare;

      for (const auto i_op : active_multipoles) {
#pragma omp parallel if (parallel_q)
        {
          // Per thread: its own operator (rank set here, frequency per q)
          // and solver. Static schedule: consecutive q on one thread, so
          // each solve warm starts from the neighbouring q
          auto h = multipoles[i_op]->clone();
          h->updateRank(k);
          ExternalField::TDHFcntm rpa(h.get(), vHF);
          rpa.eps_target() = rpa_options.eps;
#pragma omp for schedule(static)
          for (std::size_t iq = 0; iq < q_steps; ++iq) {
            h->updateFrequency(frequencies[iq]);
            const auto eps_solve =
              solve_channel_amplitudes(*h, i_op, &rpa, omega, rpa_options,
                                       channels, &A_bare[iq], &A_rpa[iq]);
            // Worst over K and operators (nan is the worst); (iE, iq) is
            // visited by one thread within this block
            auto &eps_Eq = factors.eps(iE, iq);
            if (std::isnan(eps_solve) || eps_solve > eps_Eq) {
              eps_Eq = eps_solve;
            }
            bar.update();
          }
        }
      }

      for (std::size_t iq = 0; iq < q_steps; ++iq) {
        for (std::size_t ic = 0; ic < channels.size(); ++ic) {
          const auto &channel = channels[ic];
          const auto tkp1_x = (2.0 * k + 1.0) * channel.hole->occ_frac();
          accumulate_formFactors(&factors.bare[channel.hole_index], iE, iq,
                                 tkp1_x, A_bare[iq][ic]);
          accumulate_formFactors(&factors.rpa[channel.hole_index], iE, iq,
                                 tkp1_x, A_rpa[iq][ic]);
        }
      }
    }
  }

  return factors;
}

//==============================================================================
std::pair<std::size_t, std::size_t>
interpolate_failed_rpa(FormFactorsRPA *K_rpa, double eps_fail) {
  assert(K_rpa != nullptr);

  const auto E_steps = K_rpa->eps.rows();
  const auto q_steps = K_rpa->eps.cols();

  const auto failed = [&](std::size_t iE, std::size_t iq) {
    const auto eps = K_rpa->eps(iE, iq);
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
      for (std::size_t ia = 0; ia < K_rpa->rpa.size(); ++ia) {
        for (std::size_t i = 0; i < K_rpa->rpa[ia].size(); ++i) {
          auto &K_factor = K_rpa->rpa[ia][i];
          const auto &K_bare = K_rpa->bare[ia][i];
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