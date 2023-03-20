#include "Kion_functions.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/DiagramRPA0_jL.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Methods.hpp"
#include "qip/omp.hpp"
#include <iostream>
#include <memory>

namespace Kion {

//==============================================================================
LinAlg::Matrix<double>
calculateK_nk(const HF::HartreeFock *vHF, const DiracSpinor &Fnk, int max_L,
              const Grid &Egrid, const DiracOperator::jL *jl, bool subtract_1,
              bool force_rescale, bool hole_particle, bool force_orthog,
              bool zeff_cont, bool use_rpa0,
              const std::vector<DiracSpinor> &basis) {
  assert(vHF != nullptr && "Hartree-Fock potential must not be null");

  const auto &qgrid = jl->q_grid();
  const auto qsteps = qgrid.num_points();

  LinAlg::Matrix Knk_Eq(Egrid.num_points(), qgrid.num_points());

  std::unique_ptr<ExternalField::DiagramRPA0_jL> rpa{nullptr};
  if (use_rpa0)
    rpa =
        std::make_unique<ExternalField::DiagramRPA0_jL>(jl, basis, vHF, max_L);

  if (std::abs(Fnk.en()) > Egrid.r().back()) {
    return Knk_Eq;
  }

  // Definition of matrix element:
  // matrix element defined such that:
  // K(E,q) = (2L+1) * |me|^2
  // me = <a||jL||e>
  // Now that we have 'me' directly (rather than me^2), we can -1 easily:
  // <a| jL - 1 |e> = <a| jL |e> - <a|e>
  // note: only works for vector/scalar, since for pseudo-cases,
  // we factored out factor of i from me?
  if (subtract_1 && (jl->name() != "jL" && jl->name() != "g0jL")) {
    std::cout
        << "\nWARNING: subtract 1 option currently only checked for vector "
           "and scalar operator (due to factoring out i)\n";
  }

  // Find first energy grid point for which Fnk is accessible:
  const auto idE_first_accessible = std::size_t(std::distance(
      Egrid.begin(), std::find_if(Egrid.begin(), Egrid.end(),
                                  [&](auto e) { return e > -Fnk.en(); })));
  const auto num_accessible_E_steps = Egrid.num_points() - idE_first_accessible;

  // decide what to parallelise over:
  const bool parallelise_E =
      num_accessible_E_steps >
      std::min(qsteps, (std::size_t)omp_get_max_threads());

#pragma omp parallel for if (parallelise_E)
  for (std::size_t idE = idE_first_accessible; idE < Egrid.num_points();
       ++idE) {
    const auto dE = Egrid(idE);

    // Convert energy deposition to contimuum state energy:
    double ec = dE + Fnk.en();
    if (ec <= 0.0)
      continue;

    const int l = Fnk.l();
    const int lc_max = l + max_L;
    const int lc_min = std::max(l - max_L, 0);
    // occupancy fraction. Usually 1. = N(j)/(2j+1)
    const double x_ocf = Fnk.occ_frac();

    // create cntm object [survives locally only]
    ContinuumOrbitals cntm(vHF);
    if (zeff_cont) {
      // Same Zeff as used by DarkARC (eqn B35 of arxiv:1912.08204):
      // Zeff = sqrt{I_{njl} eV / 13.6 eV} * n
      // au: Zeff = sqrt{2 * I_{njl}} * n
      const double Zeff = std::sqrt(-2.0 * Fnk.en()) * Fnk.n();
      cntm.solveContinuumZeff(ec, lc_min, lc_max, Zeff, &Fnk, force_orthog);
    } else {
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fnk, force_rescale,
                            hole_particle, force_orthog);
    }

// Generate AK for each L, lc, and q
// L and lc are summed, not stored individually
#pragma omp parallel for if (!parallelise_E)
    for (std::size_t iq = 0; iq < qsteps; iq++) {
      for (std::size_t L = 0; L <= std::size_t(max_L); L++) {
        for (const auto &Fe : cntm.orbitals) {
          if (jl->is_zero(Fe, Fnk, L))
            continue;
          const auto q = jl->q_grid().r(iq);
          auto me = jl->rme(Fe, Fnk, L, q);
          if (rpa) {
            // nb: only first-order core-pol
            me += rpa->dV_diagram_jL(Fe, Fnk, jl, L, q);
          }
          if (subtract_1 && (L == 0 && Fe.kappa() == Fnk.kappa())) {
            me -= Fe * Fnk;
          }
          Knk_Eq(idE, iq) += double(2 * L + 1) * me * me * x_ocf;
        }
      }
    }
  }
  return Knk_Eq;
}

//==============================================================================
std::vector<LinAlg::Matrix<double>> calculateK_nk_rpa(
    const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &core, int max_L,
    const Grid &Egrid, DiracOperator::jL *jl, bool subtract_1,
    bool force_rescale, bool hole_particle, bool force_orthog,
    const std::vector<DiracSpinor> &basis, const std::string &atom) {
  assert(vHF != nullptr && "Hartree-Fock potential must not be null");

  const auto &qgrid = jl->q_grid();
  const auto qsteps = qgrid.num_points();

  const int max_rpa_its = 128;

  std::vector<LinAlg::Matrix<double>> K_nk_Eq(
      core.size(), {Egrid.num_points(), qgrid.num_points()});

  // Find index of first accessible core state:
  const auto Emax = Egrid.back();
  const auto i_first_accessible_nk = std::size_t(std::distance(
      core.begin(), std::find_if(core.begin(), core.end(), [Emax](auto &Fc) {
        return std::abs(Fc.en()) < Emax;
      })));

  // This could be done more efficiently, but this seems fine.
  // Don't need to run RPA for large number of q anyway
  std::cout << "\nInlcuding RPA. RPA equations solved for each L and q\n";
  for (std::size_t L = 0; L <= std::size_t(max_L); L++) {
    jl->set_L_q(L, 0.0);
    auto rpa = ExternalField::DiagramRPA(jl, basis, vHF, atom);
    for (std::size_t iq = 0; iq < qsteps; iq++) {
      const auto q = jl->q_grid().at(iq);
      jl->set_L_q(L, q);
      rpa.update_t0s(); // Udate t0, since JL(q) changed
      rpa.solve_core(0.0, max_rpa_its, false);
      fmt::print("RPA: L={}, q={:7.2f} au [{}/{}], eps={:.1e}, its={} \r", L, q,
                 iq + 1, qsteps, rpa.get_eps(), rpa.get_its());
      if (rpa.get_eps() > 1.0e-5)
        std::cout << " **\n";
      std::cout << std::flush;
#pragma omp parallel for collapse(2)
      for (std::size_t ink = i_first_accessible_nk; ink < core.size(); ink++) {
        for (std::size_t idE = 0; idE < Egrid.num_points(); ++idE) {
          const auto &Fnk = core.at(ink);
          const auto dE = Egrid(idE);
          const auto ec = dE + Fnk.en();
          if (ec <= 0.0)
            continue;
          const auto [lc_max, lc_min] =
              std::pair{Fnk.l() + max_L, std::max(Fnk.l() - max_L, 0)};
          const double x_ocf = Fnk.occ_frac();
          ContinuumOrbitals cntm(vHF);
          cntm.solveContinuumHF(ec, lc_min, lc_max, &Fnk, force_rescale,
                                hole_particle, force_orthog);
          double K_tmp_Fe = 0.0;
          for (const auto &Fe : cntm.orbitals) {
            if (jl->isZero(Fe, Fnk))
              continue;
            auto me = jl->reducedME(Fe, Fnk) + rpa.dV(Fe, Fnk);
            if (subtract_1 && (L == 0 && Fe.kappa() == Fnk.kappa())) {
              me -= Fe * Fnk;
            }
            K_tmp_Fe += double(2 * L + 1) * me * me * x_ocf;
          }
          K_nk_Eq.at(ink).at(idE, iq) += K_tmp_Fe;
        }
      }
    }
    std::cout << "\n";
  }

  return K_nk_Eq;
}

//==============================================================================
LinAlg::Matrix<double> calculateK_nk_approx(
    const HF::HartreeFock *vHF, const std::vector<DiracSpinor> &core, int max_L,
    const DiracOperator::jL *jl, bool subtract_1, bool force_rescale,
    bool hole_particle, bool force_orthog, bool zeff_cont, bool use_rpa,
    const std::vector<DiracSpinor> &basis) {

  assert(vHF != nullptr && "Hartree-Fock potential must not be null");

  const auto &qgrid = jl->q_grid();
  const auto qsteps = qgrid.num_points();

  LinAlg::Matrix K_nk_q(core.size(), qgrid.num_points());

  std::unique_ptr<ExternalField::DiagramRPA0_jL> rpa;
  if (use_rpa)
    rpa =
        std::make_unique<ExternalField::DiagramRPA0_jL>(jl, basis, vHF, max_L);

  // Definition of matrix element:
  // matrix element defined such that:
  // K(E,q) = (2L+1) * |me|^2
  // me = <a||jL||e>
  // Now that we have 'me' directly (rather than me^2), we can -1 easily:
  // <a| jL - 1 |e> = <a| jL |e> - <a|e>
  // note: only works for vector/scalar, since for pseudo-cases,
  // we factored out factor of i from me?
  if (subtract_1 && (jl->name() != "jL" && jl->name() != "g0jL")) {
    std::cout
        << "\nWARNING: subtract 1 option currently only checked for vector "
           "and scalar operator (due to factoring out i)\n";
  }

  for (std::size_t ink = 0; ink < core.size(); ink++) {
    const auto &Fnk = core.at(ink);

    double e_cntm = 0.5; // solve all at 0.5 au - roughly OK?

    const int l = Fnk.l();
    const int lc_max = l + max_L;
    const int lc_min = std::max(l - max_L, 0);
    // occupancy fraction. Usually 1. = N(j)/(2j+1)
    const double x_ocf = Fnk.occ_frac();

    // create cntm object [survives locally only]
    ContinuumOrbitals cntm(vHF);
    if (zeff_cont) {
      // Same Zeff as used by DarkARC (eqn B35 of arxiv:1912.08204):
      // Zeff = sqrt{I_{njl} eV / 13.6 eV} * n
      // au: Zeff = sqrt{2 * I_{njl}} * n
      const double Zeff = std::sqrt(-2.0 * Fnk.en()) * Fnk.n();
      cntm.solveContinuumZeff(e_cntm, lc_min, lc_max, Zeff, &Fnk, force_orthog);
    } else {
      cntm.solveContinuumHF(e_cntm, lc_min, lc_max, &Fnk, force_rescale,
                            hole_particle, force_orthog);
    }

    // Generate AK for each L, lc, and q
    // L and lc are summed, not stored individually
#pragma omp parallel for
    for (std::size_t iq = 0; iq < qsteps; iq++) {
      for (std::size_t L = 0; L <= std::size_t(max_L); L++) {
        for (const auto &Fe : cntm.orbitals) {
          if (jl->is_zero(Fe, Fnk, L))
            continue;
          const auto q = jl->q_grid().r(iq);
          auto me = jl->rme(Fe, Fnk, L, q);
          if (rpa) {
            // nb: only first-order core-pol
            me += rpa->dV_diagram_jL(Fe, Fnk, jl, L, q);
          }
          if (subtract_1 && (L == 0 && Fe.kappa() == Fnk.kappa())) {
            me -= Fe * Fnk;
          }
          K_nk_q(ink, iq) += double(2 * L + 1) * me * me * x_ocf;
        }
      }
    }
  }
  return K_nk_q;
}

//==============================================================================
LinAlg::Matrix<double>
convert_K_nk_approx_to_std(const LinAlg::Matrix<double> &Kaprx,
                           const Grid &Egrid,
                           const std::vector<DiracSpinor> &core) {

  // K(E,q) = \sum_nk Theta[|e_nk|<E]* K_nk(q)

  assert(Kaprx.rows() == core.size());
  const auto qsteps = Kaprx.cols();

  LinAlg::Matrix Knk_E_q(Egrid.size(), qsteps);

  for (std::size_t iE = 0; iE < Egrid.size(); ++iE) {
    const auto dE = Egrid(iE);
    for (std::size_t ink = 0; ink < core.size(); ++ink) {
      const auto &Fnk = core.at(ink);
      if (std::abs(Fnk.en()) > dE)
        continue;
      for (std::size_t iq = 0; iq < qsteps; ++iq) {
        Knk_E_q(iE, iq) += Kaprx(ink, iq);
      }
    }
  }

  return Knk_E_q;
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
void write_to_file_xyz(const LinAlg::Matrix<double> &K, const Grid &E_grid,
                       const Grid &q_grid, const std::string &filename,
                       int num_digits, Units units) {
  // optional format argument?
  assert(K.rows() == E_grid.num_points());
  assert(K.cols() == q_grid.num_points());
  std::ofstream out_file(filename + "_xyz.txt");

  const auto unit_E = units == Units::Atomic ? 1.0 : UnitConv::Energy_au_to_keV;
  const auto unit_q =
      units == Units::Atomic ? 1.0 : UnitConv::Momentum_au_to_MeV;

  out_file << "# Kion output data file: " << filename << " - xyz format\n";
  fmt::print(out_file, "# Units: ");

  if (units == Units::Atomic) {
    fmt::print(out_file, "Atomic units. [q] = [1/a0], [E] = [E_H], [K] = 1\n");
    fmt::print(out_file, "#   E_H = m_e (c*alpha)^2 = ~0.027 keV\n");
    fmt::print(out_file, "#   1/a0 = m_e*c*alpha/hbar = E_H / "
                         "(c*alpha*hbar) = ~0.0037 MeV\n");
  } else if (units == Units::Particle) {
    fmt::print(out_file, "Particle units. [q] = MeV, [E] = keV, [K] = 1\n");
  } else {
    std::cout << "units error\n";
  }

  fmt::print(out_file, "# E           q            K\n");
  for (std::size_t iE = 0; iE < E_grid.num_points(); ++iE) {
    for (std::size_t iq = 0; iq < q_grid.num_points(); ++iq) {
      fmt::print(out_file, "{:+.{}e} {:+.{}e} {:+.{}e}\n", E_grid(iE) * unit_E,
                 num_digits, q_grid(iq) * unit_q, num_digits, K(iE, iq),
                 num_digits);
    }
  }
}

//==============================================================================
void write_to_file_matrix(const LinAlg::Matrix<double> &K, const Grid &E_grid,
                          const Grid &q_grid, const std::string &filename,
                          int num_digits, Units units) {
  // optional format argument?
  assert(K.rows() == E_grid.num_points());
  assert(K.cols() == q_grid.num_points());
  std::ofstream out_file(filename + "_mat.txt");

  const auto unit_E = units == Units::Atomic ? 1.0 : UnitConv::Energy_au_to_keV;
  const auto unit_q =
      units == Units::Atomic ? 1.0 : UnitConv::Momentum_au_to_MeV;

  out_file << "# Kion output data file: " << filename << " - matrix format\n";
  fmt::print(out_file, "# Units: ");

  if (units == Units::Atomic) {
    fmt::print(out_file, "Atomic units. [q] = [1/a0], [E] = [E_H], [K] = 1\n");
    fmt::print(out_file, "#   E_H = m_e (c*alpha)^2 = ~0.027 keV\n");
    fmt::print(out_file, "#   1/a0 = m_e*c*alpha/hbar = E_H / "
                         "(c*alpha*hbar) = ~0.0037 MeV\n");
  } else if (units == Units::Particle) {
    fmt::print(out_file, "Particle units. [q] = MeV, [E] = keV, [K] = 1\n");
  } else {
    std::cout << "units error\n";
  }

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
  for (std::size_t iE = 0; iE < E_grid.num_points(); ++iE) {
    for (std::size_t iq = 0; iq < q_grid.num_points(); ++iq) {
      fmt::print(out_file, "{:+.{}e} ", K(iE, iq), num_digits);
    }
    out_file << '\n';
  }
}

//==============================================================================
void write_to_file_gnuplot(const LinAlg::Matrix<double> &K, const Grid &E_grid,
                           const Grid &q_grid, const std::string &filename,
                           int num_digits, Units units) {
  // optional format argument?
  assert(K.rows() == E_grid.num_points());
  assert(K.cols() == q_grid.num_points());
  std::ofstream out_file(filename + "_gnu.txt");

  const auto unit_E = units == Units::Atomic ? 1.0 : UnitConv::Energy_au_to_keV;
  const auto unit_q =
      units == Units::Atomic ? 1.0 : UnitConv::Momentum_au_to_MeV;

  out_file << "# Kion output data file: " << filename << " - gnuplot format\n";
  fmt::print(out_file, "# Units: ");

  if (units == Units::Atomic) {
    fmt::print(out_file, "Atomic units. [q] = [1/a0], [E] = [E_H], [K] = 1\n");
    fmt::print(out_file, "#   E_H = m_e (c*alpha)^2 = ~0.027 keV\n");
    fmt::print(out_file, "#   1/a0 = m_e*c*alpha/hbar = E_H / "
                         "(c*alpha*hbar) = ~0.0037 MeV\n");
  } else if (units == Units::Particle) {
    fmt::print(out_file, "Particle units. [q] = MeV, [E] = keV, [K] = 1\n");
  } else {
    std::cout << "units error\n";
  }

  fmt::print(out_file, "# K_ion(E,q).\n"
                       "# Each column is new E; each row is new q.\n"
                       "# First column is q values; first row is E values.\n");

  fmt::print(out_file, "{:{}s} ", "q\\E", num_digits + 7);
  for (const auto &E : E_grid) {
    fmt::print(out_file, "{:+.{}e} ", E * unit_E, num_digits);
  }
  out_file << "\n";

  // nb: write out is 'transposed' order from array: this is to make plotting easier
  for (std::size_t iq = 0; iq < q_grid.num_points(); ++iq) {
    fmt::print(out_file, "{:+.{}e} ", q_grid(iq) * unit_q, num_digits);
    for (std::size_t iE = 0; iE < E_grid.num_points(); ++iE) {
      fmt::print(out_file, "{:+.{}e} ", K(iE, iq), num_digits);
    }
    fmt::print(out_file, "\n");
  }
}

//==============================================================================
void write_approxTable_to_file(const LinAlg::Matrix<double> &K,
                               const std::vector<DiracSpinor> &core,
                               const Grid &q_grid, const std::string &filename,
                               int num_digits, Units units) {

  assert(K.rows() == core.size());
  assert(K.cols() == q_grid.num_points());
  std::ofstream out_file(filename + "_approxTable.txt");

  const auto unit_q =
      units == Units::Atomic ? 1.0 : UnitConv::Momentum_au_to_MeV;

  out_file << "# Kion Approximate Tables output data file: " << filename
           << "\n";
  fmt::print(out_file, "# Units: ");
  if (units == Units::Atomic) {
    fmt::print(out_file, "Atomic units. [q] = [1/a0], [E] = [E_H], [K] = 1\n");
    fmt::print(out_file, "#   E_H = m_e (c*alpha)^2 = ~0.027 keV\n");
    fmt::print(out_file, "#   1/a0 = m_e*c*alpha/hbar = E_H / "
                         "(c*alpha*hbar) = ~0.0037 MeV\n");
  } else if (units == Units::Particle) {
    fmt::print(out_file, "Particle units. [q] = MeV, [E] = keV, [K] = 1\n");
  } else {
    std::cout << "units error\n";
  }

  out_file << "# Each column is a new state (first col is q values)\n";
  out_file << "# Each row is a new q\n";
  out_file << "# First row is state labels\n";

  fmt::print(out_file, "q\\state ");
  for (const auto &Fnk : core) {
    fmt::print(out_file, "{} ", Fnk.shortSymbol());
  }
  out_file << "\n";

  for (std::size_t iq = 0; iq < q_grid.num_points(); ++iq) {
    fmt::print(out_file, "{:+.{}e} ", q_grid(iq) * unit_q, num_digits);
    for (std::size_t ink = 0; ink < core.size(); ++ink) {
      fmt::print(out_file, "{:+.{}e} ", K(ink, iq), num_digits);
    }
    fmt::print(out_file, "\n");
  }
}

//==============================================================================
void write_to_file(const std::vector<OutputFormat> &formats,
                   const LinAlg::Matrix<double> &K, const Grid &E_grid,
                   const Grid &q_grid, const std::string &filename,
                   int num_digits, Units units) {
  // Update: Always use atomic units for mat and xyz
  for (const auto &format : formats) {
    if (format == OutputFormat::gnuplot)
      write_to_file_gnuplot(K, E_grid, q_grid, filename, num_digits, units);
    if (format == OutputFormat::matrix)
      write_to_file_matrix(K, E_grid, q_grid, filename, num_digits,
                           Units::Atomic);
    if (format == OutputFormat::xyz)
      write_to_file_xyz(K, E_grid, q_grid, filename, num_digits, Units::Atomic);
  }
}

} // namespace Kion