#include "Kion_functions.hpp"
#include "DiracOperator/include.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Maths.hpp"
#include "qip/Methods.hpp"
#include "qip/omp.hpp"
#include <iostream>
#include <memory>

namespace Kion {

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
  bool axialQ, bool scalarQ, bool pseudoscalarQ, bool spatialQ) {

  // can never get this far, but debugging check
  assert(vHF != nullptr);

  if (diagonal_Eq)
    assert(qgrid.size() == 1);

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
#pragma omp parallel for
  for (std::size_t iE = idE_0; iE < idE_max; ++iE) {
    const auto omega = Egrid.at(iE);

    const auto ec = omega + Fa.en();
    assert(ec >= ec_min);
    assert(ec <= ec_max);

    ContinuumOrbitals cntm(vHF);
    cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa, force_rescale, hole_particle,
                          force_orthog);

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
          const auto t = Phik ? Phik->reducedME(Fe, Fa) : 0.0;
          const auto E = Ek ? Ek->reducedME(Fe, Fa) : 0.0;
          const auto M = Mk ? Mk->reducedME(Fe, Fa) : 0.0;
          const auto L = Lk ? Lk->reducedME(Fe, Fa) : 0.0;
          // axial:
          const auto t5 = Phi5k ? Phi5k->reducedME(Fe, Fa) : 0.0;
          const auto E5 = E5k ? E5k->reducedME(Fe, Fa) : 0.0;
          const auto M5 = M5k ? M5k->reducedME(Fe, Fa) : 0.0;
          const auto L5 = L5k ? L5k->reducedME(Fe, Fa) : 0.0;
          // Scalar, Pseudoscalar
          const auto S = Sk ? Sk->reducedME(Fe, Fa) : 0.0;
          const auto S5 = S5k ? S5k->reducedME(Fe, Fa) : 0.0;

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
  }

  return K_factors;
}

//==============================================================================
bool check_radial_grid(double Emax_au, double qmax_au, const Grid &rgrid) {
  // Check grid type: only loglinear is reasonable for this module
  if (rgrid.type() != GridType::loglinear) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning:\n");
    std::cout << "This module unlikely to work with grid type: "
              << GridParameters::parseType(rgrid.type())
              << "; consider changing to loglinear\n";
  }

  // *Very* rough estimate of good aximum q range
  const auto i = rgrid.getIndex(1.0);
  const auto dr_q = rgrid.drdu(i) * rgrid.du();
  const auto qmax_targ = 1.5 * 2.0 * M_PI / (3.0 * dr_q);
  fmt::print("\nVery rough guess at maximum safe q: {:.0f} au = {:.2f} MeV\n",
             qmax_targ, qmax_targ * UnitConv::Momentum_au_to_MeV);
  if (qmax_au > qmax_targ) {
    fmt::print("Warning: Grid may not be dense enough for highest q\n\n");
  }

  // Check if grid dense enough. If not, offers advice
  const double approx_wavelength = 2.0 * M_PI / std::sqrt(2.0 * Emax_au);
  const auto dr = rgrid.du(); //rough
  const int N_ppw = 15;
  const double dr_target = approx_wavelength / N_ppw;
  fmt::print("Approx continuum wavelength for max E: {:.3f}.\n"
             "Approx dr at large r: {:.4f}\n",
             approx_wavelength, dr);
  if (dr > dr_target) {
    fmt2::styled_print(fg(fmt::color::orange), "Warning: ");
    fmt::print(
      "Grid may not be dense enough for continuum state with e={:.2f}\n",
      Emax_au);
    fmt::print("Have dr~{:.3f} at large r, but need dr~{:.3f}\n", dr,
               dr_target);

    // For given b, find required num_points:
    const auto n_target =
      Grid::calc_num_points_from_du(rgrid.r0(), rgrid.rmax(), 0.9 * dr_target,
                                    GridType::loglinear, rgrid.loglin_b());
    // For given num_points, find required b [by solving: du(b)-du_targ = 0]
    const auto f = [&](double b) {
      return Grid::calc_du_from_num_points(rgrid.r0(), rgrid.rmax(),
                                           rgrid.num_points(),
                                           GridType::loglinear, b) -
             0.9 * dr_target;
    };
    const auto [b_target, db] =
      qip::Newtons(f, 1.0, {0.05, 100.0}, 1.0e-1, 0.01);
    fmt::print("Try increasing num_points to {}; or decreasing b to {:.2f}\n",
               n_target, b_target);
    std::cout << "(Program will continue, but results may be inaccurate)\n\n";
    return false;
  }
  return true;
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