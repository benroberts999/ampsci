#include "Kion_functions.hpp"
#include "DiracODE/include.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHFcntm.hpp"
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
#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <optional>
#include <type_traits>

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
std::array<LinAlg::Matrix<double>, 13> calculate_formFactors_nk(
  const HF::HartreeFock *vHF, const DiracSpinor &Fa, int lc_min, int lc_max,
  double ec_min, double ec_max, bool force_rescale, bool hole_particle,
  bool force_orthog, const std::vector<double> &Egrid,
  const std::vector<double> &qgrid, bool diagonal_Eq, bool low_q,
  const SphericalBessel::JL_table &jK_tab, int Kmin, int Kmax, bool vectorQ,
  bool axialQ, bool scalarQ, bool pseudoscalarQ, bool spatialQ,
  AtomicMethod method, double zeff_constant) {

  // can never get this far, but debugging check
  assert(vHF != nullptr);

  if (diagonal_Eq) {
    assert(qgrid.size() == 1);
  }

  // H-like (Zeff) approximation: constant Zeff (if given), else "real"
  // Zeff = n*sqrt(-2*en) from the binding energy (same as DarkARC):
  const bool use_zeff = method != AtomicMethod::HF;
  const double Zeff =
    zeff_constant > 0.0 ? zeff_constant : Zeff_real(Fa.en(), Fa.n());

  // Bound state used in the matrix elements: the real state, or its H-like
  // (Zeff, pointlike) version, solved analytically or numerically (DiracODE).
  // nb: Zeff state will not have exactly right energy; continuum energies
  // (ec = E + en) and occupation always use the real (input) Fa.
  std::optional<DiracSpinor> Fa_zeff{};
  if (method == AtomicMethod::ZeffAnalytic) {
    Fa_zeff = DiracSpinor::exactHlike(Fa.n(), Fa.kappa(), Fa.grid_sptr(), Zeff,
                                      vHF->alpha());
  } else if (method == AtomicMethod::Zeff) {
    const auto v_z =
      Nuclear::sphericalNuclearPotential(Zeff, 0.0, vHF->grid().r());
    const auto e0 = AtomData::diracen(Zeff, Fa.n(), Fa.kappa(), vHF->alpha());
    Fa_zeff = DiracODE::boundState(Fa.n(), Fa.kappa(), e0, Fa.grid_sptr(), v_z,
                                   {}, vHF->alpha());
  }
  const auto &Fa_t = use_zeff ? *Fa_zeff : Fa;

  const auto E_steps = Egrid.size();
  const auto q_steps = qgrid.size();

  // Factors (output) for specific bound state, Fa
  std::array<LinAlg::Matrix<double>, 13> K_factors;
  auto &K_VT = K_factors[0];  // Vector: temporal
  auto &K_VE = K_factors[1];  // Vector: electric
  auto &K_VM = K_factors[2];  // Vector: magnetic
  auto &K_VL = K_factors[3];  // Vector: longitudinal
  auto &K_X = K_factors[4];   // Vector: v-v cross
  auto &K_T5 = K_factors[5];  // Axial: temporal
  auto &K_E5 = K_factors[6];  // Axial: electric
  auto &K_M5 = K_factors[7];  // Axial: magnetic
  auto &K_L5 = K_factors[8];  // Axial: longitudinal
  auto &K_X5 = K_factors[9];  // Axial: a-a cross
  auto &K_Z = K_factors[10];  // Vector-Axial interference
  auto &K_S = K_factors[11];  // Scalar
  auto &K_S5 = K_factors[12]; // Pseudo-scalar

  // Re-size output arrays (if required)
  if (vectorQ) {
    K_VT.resize(E_steps, q_steps);
    if (spatialQ) {
      K_VE.resize(E_steps, q_steps);
      K_VM.resize(E_steps, q_steps);
      K_VL.resize(E_steps, q_steps);
      K_X.resize(E_steps, q_steps);
    }
  }
  if (axialQ) {
    K_T5.resize(E_steps, q_steps);
    if (spatialQ) {
      K_E5.resize(E_steps, q_steps);
      K_M5.resize(E_steps, q_steps);
      K_L5.resize(E_steps, q_steps);
      K_X5.resize(E_steps, q_steps);
    }
  }
  if (vectorQ && axialQ && spatialQ) {
    K_Z.resize(E_steps, q_steps);
  }
  if (scalarQ) {
    K_S.resize(E_steps, q_steps);
  }
  if (pseudoscalarQ) {
    K_S5.resize(E_steps, q_steps);
  }

  // Find which parts of E grid can contribute (more efficient parallelisation)
  const auto idE_0_tmp = std::distance(
    Egrid.begin(), std::find_if(Egrid.begin(), Egrid.end(),
                                [&](auto e) { return e + Fa.en() > ec_min; }));
  assert(idE_0_tmp >= 0);
  const auto idE_0 = std::size_t(idE_0_tmp);

  const auto idE_max = std::size_t(std::distance(
    Egrid.begin(), std::find_if(Egrid.begin(), Egrid.end(),
                                [&](auto e) { return e + Fa.en() > ec_max; })));
  assert(idE_max <= Egrid.size());

  const auto max_threads = std::size_t(omp_get_max_threads());

  // Build one operator set per thread. updateRank() and updateFrequency() are
  // called inside the loop before use; operators will be null for invalid
  // type/comp combinations, so all uses are null-guarded below.
  auto build_thread_ops = [&](char type, char comp, bool include)
    -> std::vector<std::unique_ptr<DiracOperator::TensorOperator>> {
    const auto op0 =
      include ? DiracOperator::MultipoleOperator(vHF->grid(), 0, 0.0, type,
                                                 comp, low_q, &jK_tab) :
                nullptr;
    std::vector<std::unique_ptr<DiracOperator::TensorOperator>> ops;
    ops.reserve(max_threads);
    for (auto t = 0ul; t < max_threads; ++t)
      ops.push_back(op0 ? op0->clone() : nullptr);
    return ops;
  };

  auto Phik_ops = build_thread_ops('V', 'T', vectorQ);
  auto Ek_ops = build_thread_ops('V', 'E', vectorQ && spatialQ);
  auto Mk_ops = build_thread_ops('V', 'M', vectorQ && spatialQ);
  auto Lk_ops = build_thread_ops('V', 'L', vectorQ && spatialQ);
  auto Phi5k_ops = build_thread_ops('A', 'T', axialQ);
  auto E5k_ops = build_thread_ops('A', 'E', axialQ && spatialQ);
  auto M5k_ops = build_thread_ops('A', 'M', axialQ && spatialQ);
  auto L5k_ops = build_thread_ops('A', 'L', axialQ && spatialQ);
  auto Sk_ops = build_thread_ops('S', 'T', scalarQ);
  auto S5k_ops = build_thread_ops('P', 'T', pseudoscalarQ);

  // nb: sum into K(iE,iq).
  // Therefore, each thread has unique iE, and no reduction required if //-isation over iE/iq
  // BUT if we change this (e.g., over k), then this will change
  qip::ProgressBar bar(idE_max - idE_0);
#pragma omp parallel for schedule(dynamic)
  for (std::size_t iE = idE_0; iE < idE_max; ++iE) {
    const auto omega = Egrid.at(iE);

    const auto ec = omega + Fa.en();
    assert(ec >= ec_min);
    assert(ec <= ec_max);

    ContinuumOrbitals cntm(vHF);
    if (method == AtomicMethod::Zeff) {
      cntm.solveContinuumZeff(ec, lc_min, lc_max, Zeff, &Fa_t, force_orthog);
    } else if (method == AtomicMethod::ZeffAnalytic) {
      cntm.solveContinuumZeffAnalytic(ec, lc_min, lc_max, Zeff, &Fa_t,
                                      force_orthog);
    } else {
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa_t, force_rescale,
                            hole_particle, force_orthog);
    }

    const auto tid = std::size_t(omp_get_thread_num());
    auto Phik = Phik_ops[tid].get();
    auto Ek = Ek_ops[tid].get();
    auto Mk = Mk_ops[tid].get();
    auto Lk = Lk_ops[tid].get();
    auto Phi5k = Phi5k_ops[tid].get();
    auto E5k = E5k_ops[tid].get();
    auto M5k = M5k_ops[tid].get();
    auto L5k = L5k_ops[tid].get();
    auto Sk = Sk_ops[tid].get();
    auto S5k = S5k_ops[tid].get();

    for (int k = Kmin; k <= Kmax; ++k) {
      const auto tkp1_x = (2.0 * k + 1.0) * Fa.occ_frac();

      for (std::size_t iq = 0; iq < qgrid.size(); ++iq) {

        // Use qc as expected for "omega" in operators (just units)
        const auto qc =
          diagonal_Eq ? Egrid.at(iE) : qgrid.at(iq) * PhysConst::c;

        // Update rank (adjusts parity), then frequency (resets Bessel vectors)
        if (Phik) {
          Phik->updateRank(k);
          Phik->updateFrequency(qc);
        }
        if (Ek) {
          Ek->updateRank(k);
          Ek->updateFrequency(qc);
        }
        if (Mk) {
          Mk->updateRank(k);
          Mk->updateFrequency(qc);
        }
        if (Lk) {
          Lk->updateRank(k);
          Lk->updateFrequency(qc);
        }
        if (Phi5k) {
          Phi5k->updateRank(k);
          Phi5k->updateFrequency(qc);
        }
        if (E5k) {
          E5k->updateRank(k);
          E5k->updateFrequency(qc);
        }
        if (M5k) {
          M5k->updateRank(k);
          M5k->updateFrequency(qc);
        }
        if (L5k) {
          L5k->updateRank(k);
          L5k->updateFrequency(qc);
        }
        if (Sk) {
          Sk->updateRank(k);
          Sk->updateFrequency(qc);
        }
        if (S5k) {
          S5k->updateRank(k);
          S5k->updateFrequency(qc);
        }

        for (const auto &Fe : cntm.orbitals) {

          // vector:
          const auto t = Phik ? Phik->reducedME(Fe, Fa_t) : 0.0;
          const auto E = Ek ? Ek->reducedME(Fe, Fa_t) : 0.0;
          const auto M = Mk ? Mk->reducedME(Fe, Fa_t) : 0.0;
          const auto L = Lk ? Lk->reducedME(Fe, Fa_t) : 0.0;
          // axial:
          const auto t5 = Phi5k ? Phi5k->reducedME(Fe, Fa_t) : 0.0;
          const auto E5 = E5k ? E5k->reducedME(Fe, Fa_t) : 0.0;
          const auto M5 = M5k ? M5k->reducedME(Fe, Fa_t) : 0.0;
          const auto L5 = L5k ? L5k->reducedME(Fe, Fa_t) : 0.0;
          // Scalar, Pseudoscalar
          const auto S = Sk ? Sk->reducedME(Fe, Fa_t) : 0.0;
          const auto S5 = S5k ? S5k->reducedME(Fe, Fa_t) : 0.0;

          // Vector operators
          if (vectorQ) {
            K_VT(iE, iq) += tkp1_x * qip::pow(t, 2);
            if (spatialQ) {
              K_VE(iE, iq) += tkp1_x * qip::pow(E, 2);
              K_VM(iE, iq) += tkp1_x * qip::pow(M, 2);
              K_VL(iE, iq) += tkp1_x * qip::pow(L, 2);
              K_X(iE, iq) += tkp1_x * t * L;
            }
          }

          // Axial (γ^5) operators
          if (axialQ) {
            K_T5(iE, iq) += tkp1_x * qip::pow(t5, 2);
            if (spatialQ) {
              K_E5(iE, iq) += tkp1_x * qip::pow(E5, 2);
              K_M5(iE, iq) += tkp1_x * qip::pow(M5, 2);
              K_L5(iE, iq) += tkp1_x * qip::pow(L5, 2);
              K_X5(iE, iq) += tkp1_x * t5 * L5;
            }
          }

          // Vector-Axial Spatial Interference:
          if (vectorQ && axialQ && spatialQ) {
            K_Z(iE, iq) += tkp1_x * (E5 * M - E * M5);
          }

          // Scalar and Pseudoscalar
          if (scalarQ) {
            K_S(iE, iq) += tkp1_x * qip::pow(S, 2);
          }
          if (pseudoscalarQ) {
            K_S5(iE, iq) += tkp1_x * qip::pow(S5, 2);
          }
        }
      }
    }
    bar.update();
  }

  return K_factors;
}

//==============================================================================
RPAFormFactors calculate_formFactors_rpa(
  const HF::HartreeFock *vHF, double ec_min, double ec_max,
  const std::vector<double> &Egrid, const std::vector<double> &qgrid,
  bool diagonal_Eq, bool low_q, const SphericalBessel::JL_table &jK_tab,
  int Kmin, int Kmax, bool vectorQ, bool axialQ, bool scalarQ,
  bool pseudoscalarQ, bool spatialQ, const RPAOptions &options, bool print) {

  assert(vHF != nullptr);
  if (diagonal_Eq) {
    assert(qgrid.size() == 1);
  }

  const auto &core = vHF->core();
  const auto n_core = core.size();
  const auto E_steps = Egrid.size();
  const auto q_steps = qgrid.size();

  // Operators (null if not included), as calculate_formFactors_nk. The rank
  // is set per k; the frequency is set per q, on per-thread clones.
  const auto &grid = vHF->grid();
  const auto Phik = vectorQ ? DiracOperator::MultipoleOperator(
                                grid, Kmin, 0.0, 'V', 'T', low_q, &jK_tab) :
                              nullptr;
  const auto Ek = vectorQ && spatialQ ?
                    DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'V', 'E',
                                                     low_q, &jK_tab) :
                    nullptr;
  const auto Mk = vectorQ && spatialQ ?
                    DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'V', 'M',
                                                     low_q, &jK_tab) :
                    nullptr;
  const auto Lk = vectorQ && spatialQ ?
                    DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'V', 'L',
                                                     low_q, &jK_tab) :
                    nullptr;
  const auto Phi5k = axialQ ? DiracOperator::MultipoleOperator(
                                grid, Kmin, 0.0, 'A', 'T', low_q, &jK_tab) :
                              nullptr;
  const auto E5k = axialQ && spatialQ ?
                     DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'A', 'E',
                                                      low_q, &jK_tab) :
                     nullptr;
  const auto M5k = axialQ && spatialQ ?
                     DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'A', 'M',
                                                      low_q, &jK_tab) :
                     nullptr;
  const auto L5k = axialQ && spatialQ ?
                     DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'A', 'L',
                                                      low_q, &jK_tab) :
                     nullptr;
  const auto Sk = scalarQ ? DiracOperator::MultipoleOperator(
                              grid, Kmin, 0.0, 'S', 'T', low_q, &jK_tab) :
                            nullptr;
  const auto S5k = pseudoscalarQ ?
                     DiracOperator::MultipoleOperator(grid, Kmin, 0.0, 'P', 'T',
                                                      low_q, &jK_tab) :
                     nullptr;
  // In a fixed order; the amplitudes below are stored by the same index
  const std::array<DiracOperator::TensorOperator *, 10> operators{
    Phik.get(), Ek.get(),  Mk.get(),  Lk.get(), Phi5k.get(),
    E5k.get(),  M5k.get(), L5k.get(), Sk.get(), S5k.get()};

  // Operator "frequency": the multipoles take q = alpha*omega_op, so pass
  // qc (or E itself in the diagonal, massless-absorption, case)
  const auto qc_at = [&](std::size_t iE, std::size_t iq) {
    return diagonal_Eq ? Egrid.at(iE) : qgrid.at(iq) * PhysConst::c;
  };

  // Output: the 13 factors of each core orbital, allocated (empty if not
  // requested) as calculate_formFactors_nk
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
  RPAFormFactors out;
  out.K_nk.resize(n_core);
  for (auto &K_factors : out.K_nk) {
    for (std::size_t i = 0; i < K_factors.size(); ++i) {
      if (requested[i]) {
        K_factors[i].resize(E_steps, q_steps);
      }
    }
  }
  out.rpa_eps.assign(E_steps, 0.0);
  out.rpa_its.assign(E_steps, 0.0);
  out.KpiD_dev.assign(E_steps, 0.0);
  out.Kbar_asym.assign(E_steps, 0.0);

  // The RPA solver: one instance, its core-only tables built once here. Its
  // operator is switched per (rank, parity) class below; per-thread copies
  // carry its omega-only caches (channel pairs; K matrix, or the
  // outgoing-wave incident waves). Either the standing-wave + K-matrix
  // solver (TDHFcntm) or the outgoing-wave one (TDHFcomplex), which has no
  // standing-wave amplitudes: unitarised always.
  DiracOperator::TensorOperator *any_operator = nullptr;
  for (auto *h : operators) {
    if (h != nullptr) {
      any_operator = h;
    }
  }
  if (any_operator == nullptr) {
    return out;
  }
  const bool unitarise = options.unitarise || options.outgoing_wave;

  // The scan, generic in the solver type (TDHFcntm: standing-wave solve +
  // K matrix; TDHFcomplex: outgoing-wave solve). The K-matrix-specific
  // calls exist only for TDHFcntm and are compiled only for it.
  const auto run_scan = [&](auto &rpa) {
    using Solver = std::decay_t<decltype(rpa)>;
    constexpr bool has_kmatrix =
      std::is_same_v<Solver, ExternalField::TDHFcntm>;
    rpa.eps_target() = options.eps;
    if constexpr (has_kmatrix) {
      rpa.set_unitarise(unitarise);
    }

    if (print) {
      fmt::print("\nRPA form factors: {} energies, {} q points, K = {}-{}, {} "
                 "amplitudes ({})\n",
                 E_steps, q_steps, Kmin, Kmax,
                 unitarise ? "unitarised" : "standing-wave",
                 options.outgoing_wave ? "outgoing-wave solve" :
                                         "standing-wave solve + K matrix");
      fmt::print("{:>9s} {:>5s} {:>8s} {:>4s} {:>8s} {:>8s}\n", "E/eV", "chans",
                 "eps", "its", "KpiD", "Kasym");
    }

    // Energies run serially: threads are used inside (over the K-matrix
    // columns, and over q or inside each driven solve)
    for (std::size_t iE = 0; iE < E_steps; ++iE) {
      const auto omega = Egrid.at(iE);
      if (options.E_max > 0.0 && omega > options.E_max) {
        continue;
      }

      // Nothing to output unless some orbital is ionised within the ec limits
      bool any_open = false;
      for (const auto &Fa : core) {
        const auto ec = omega + Fa.en();
        if (ec > 0.0 && ec >= ec_min && ec <= ec_max) {
          any_open = true;
        }
      }
      if (!any_open) {
        continue;
      }

      // Standing-wave option: the orthogonalised V^{N-1} HF continuum bra of
      // each ionised orbital, at every l reachable from it by a rank up to
      // Kmax (the unitarised amplitudes need no bra)
      std::vector<ContinuumOrbitals> bra(n_core, ContinuumOrbitals(vHF));
      if (!unitarise) {
#pragma omp parallel for schedule(dynamic)
        for (std::size_t ib = 0; ib < n_core; ++ib) {
          const auto &Fa = core[ib];
          const auto ec = omega + Fa.en();
          if (ec <= 0.0 || ec < ec_min || ec > ec_max) {
            continue;
          }
          const auto lc_min = std::max((Fa.twoj() - 2 * Kmax - 1) / 2, 0);
          const auto lc_max = (Fa.twoj() + 2 * Kmax + 1) / 2;
          bra[ib].solveContinuumHF(ec, lc_min, lc_max, &Fa, false, true, true);
        }
      }

      std::size_t max_channels = 0;

      for (int k = Kmin; k <= Kmax; ++k) {
        for (auto *h : operators) {
          if (h != nullptr) {
            h->updateRank(k);
          }
        }

        // The operators of each (k, parity) class share the channel
        // structure, hence the channel caches and the K matrix: built once
        // per class, on the master instance (K-matrix columns in parallel),
        // with any operator of the class
        for (const int parity : {1, -1}) {
          DiracOperator::TensorOperator *class_operator = nullptr;
          for (auto *h : operators) {
            if (h != nullptr && h->parity() == parity) {
              class_operator = h;
            }
          }
          if (class_operator == nullptr) {
            continue;
          }
          rpa.set_operator(class_operator);
          if constexpr (has_kmatrix) {
            rpa.prepare(omega, options.max_its);
          } else {
            rpa.prepare(omega);
          }
          const auto channels = rpa.channel_list();
          const auto n_channels = channels.size();
          if (n_channels == 0) {
            continue;
          }
          max_channels = std::max(max_channels, n_channels);
          if constexpr (has_kmatrix) {
            if (rpa.kmatrix()) {
              out.Kbar_asym[iE] =
                std::max(out.Kbar_asym[iE], rpa.kmatrix()->asymmetry);
            }
          }

          // Amplitudes A/pi of each operator, [iq][channel], indexed as
          // `operators`; zero for the operators not in this class
          std::array<std::vector<std::vector<std::complex<double>>>, 10>
            amplitudes;
          for (auto &A : amplitudes) {
            A.assign(q_steps,
                     std::vector<std::complex<double>>(n_channels, 0.0));
          }

          for (std::size_t i_op = 0; i_op < operators.size(); ++i_op) {
            const auto *h = operators[i_op];
            if (h == nullptr || h->parity() != parity) {
              continue;
            }
            auto &A = amplitudes[i_op];
#pragma omp parallel if (options.parallel_q)
            {
              // Per thread: its own operator (frequency set per q), and a
              // copy of the prepared solver (carries the channel caches and
              // K matrix; its own corrections and Anderson history). With
              // parallel_q the solve's own parallel regions are nested,
              // hence inert; otherwise this region has one thread and they
              // are active.
              auto h_thread = h->clone();
              auto rpa_thread = rpa;
              rpa_thread.set_operator(h_thread.get());

            // Static schedule: consecutive q on one thread, so each solve
            // warm-starts from the neighbouring q
#pragma omp for schedule(static)
              for (std::size_t iq = 0; iq < q_steps; ++iq) {
                h_thread->updateFrequency(qc_at(iE, iq));
                rpa_thread.solve_core(omega, options.max_its, false);
                if (unitarise) {
                  A[iq] = rpa_thread.A_phys();
                  assert(A[iq].size() == n_channels);
                } else {
                  // Standing-wave amplitude D = <e|h + dV|a> with the HF bra
                  for (std::size_t i = 0; i < n_channels; ++i) {
                    const auto &Fa = core[channels[i].i_core];
                    for (const auto &Fe : bra[channels[i].i_core].orbitals) {
                      if (Fe.kappa() == channels[i].kappa &&
                          Fe.norm2() != 0.0) {
                        A[iq][i] =
                          h_thread->reducedME(Fe, Fa) + rpa_thread.dV(Fe, Fa);
                      }
                    }
                  }
                }
#pragma omp critical(kion_rpa_diagnostics)
                {
                  out.rpa_eps[iE] =
                    std::max(out.rpa_eps[iE], rpa_thread.last_eps());
                  out.rpa_its[iE] =
                    std::max(out.rpa_its[iE], rpa_thread.last_its());
                  out.KpiD_dev[iE] =
                    std::max(out.KpiD_dev[iE], rpa_thread.KpiD_dev());
                }
              }
            }
          }

          // Accumulate as calculate_formFactors_nk, with |A|^2 for the
          // squares and Re(A conj(A')) for the interference terms (every
          // interference pair lies within one parity class). Channels outside
          // the ec limits are left out of the output only.
          for (std::size_t iq = 0; iq < q_steps; ++iq) {
            for (std::size_t i = 0; i < n_channels; ++i) {
              const auto &channel = channels[i];
              if (channel.en < ec_min || channel.en > ec_max) {
                continue;
              }
              const auto tkp1_x =
                (2.0 * k + 1.0) * core[channel.i_core].occ_frac();
              auto &K_factors = out.K_nk[channel.i_core];
              const auto t = amplitudes[0][iq][i];
              const auto E = amplitudes[1][iq][i];
              const auto M = amplitudes[2][iq][i];
              const auto L = amplitudes[3][iq][i];
              const auto t5 = amplitudes[4][iq][i];
              const auto E5 = amplitudes[5][iq][i];
              const auto M5 = amplitudes[6][iq][i];
              const auto L5 = amplitudes[7][iq][i];
              const auto S = amplitudes[8][iq][i];
              const auto S5 = amplitudes[9][iq][i];

              // Vector operators
              if (vectorQ) {
                K_factors[0](iE, iq) += tkp1_x * std::norm(t); // VT
                if (spatialQ) {
                  K_factors[1](iE, iq) += tkp1_x * std::norm(E); // VE
                  K_factors[2](iE, iq) += tkp1_x * std::norm(M); // VM
                  K_factors[3](iE, iq) += tkp1_x * std::norm(L); // VL
                  K_factors[4](iE, iq) +=
                    tkp1_x * std::real(t * std::conj(L)); // X
                }
              }

              // Axial operators
              if (axialQ) {
                K_factors[5](iE, iq) += tkp1_x * std::norm(t5); // T5
                if (spatialQ) {
                  K_factors[6](iE, iq) += tkp1_x * std::norm(E5); // E5
                  K_factors[7](iE, iq) += tkp1_x * std::norm(M5); // M5
                  K_factors[8](iE, iq) += tkp1_x * std::norm(L5); // L5
                  K_factors[9](iE, iq) +=
                    tkp1_x * std::real(t5 * std::conj(L5)); // X5
                }
              }

              // Vector-Axial spatial interference
              if (vectorQ && axialQ && spatialQ) {
                K_factors[10](iE, iq) +=
                  tkp1_x *
                  std::real(E5 * std::conj(M) - E * std::conj(M5)); // Z
              }

              // Scalar and Pseudoscalar
              if (scalarQ) {
                K_factors[11](iE, iq) += tkp1_x * std::norm(S); // S
              }
              if (pseudoscalarQ) {
                K_factors[12](iE, iq) += tkp1_x * std::norm(S5); // S5
              }
            }
          }
        }
      }

      if (print) {
        fmt::print("{:9.2f} {:5d} {:8.1e} {:4.0f} {:8.1e} {:8.1e}\n",
                   omega * PhysConst::Hartree_eV, max_channels, out.rpa_eps[iE],
                   out.rpa_its[iE], out.KpiD_dev[iE], out.Kbar_asym[iE]);
        std::cout << std::flush;
      }
    }
  };

  if (options.outgoing_wave) {
    ExternalField::TDHFcomplex rpa(any_operator, vHF);
    run_scan(rpa);
  } else {
    ExternalField::TDHFcntm rpa(any_operator, vHF);
    run_scan(rpa);
  }

  return out;
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
void write_to_file_xyz_13(
  const std::string &filename, const std::vector<double> &E_grid,
  const std::vector<double> &q_grid, const std::vector<std::string> &titles,
  const std::vector<std::string> &descriptions,
  const std::array<LinAlg::Matrix<double>, 13> K_factors, Units units,
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