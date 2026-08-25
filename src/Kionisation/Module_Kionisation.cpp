#include "DiracODE/ContinuumState.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHFcntm.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Kionisation/Kion_functions.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Modules/Modules.hpp"
#include "Physics/DiracContinuum.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Maths.hpp"
#include "qip/Methods.hpp"
#include "qip/Widgets.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <memory>

namespace Module {

// Declare, register, then define below.
void Kionisation(const IO::InputBlock &input, const Wavefunction &wf);
void photo(const IO::InputBlock &input, const Wavefunction &wf);
void photoRPA(const IO::InputBlock &input, const Wavefunction &wf);
void formFactors(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_Kionisation{
  "Kionisation", "Calculate atomic ionisation form-factors", &Kionisation};
const Register r_photo{
  "photo", "Calculate atomic photo-ionisation form-factors", &photo};
const Register r_photoRPA{
  "photoRPA", "Photo-ionisation cross-section with RPA (continuum TDHF)",
  &photoRPA};
const Register r_formFactors{"formFactors",
                             "Calculate general atomic ionisation form-factors",
                             &formFactors};
} // namespace

//==============================================================================

void Kionisation(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("Kionisation");

  input.check(
    {{"", "Calculates atomic ionisation factors, K.\n"
          "Specifically, it calculates the temporal part of the form factor."
          "This module is superceded by the newer 'formFactors' module, "
          "which calculates more gereral form factors."
          "The temporal component of those should match that from here.\n"
          "Output formats:\n"
          "  - xyz:      Standard output. Each row is in form: 'E q K(E,q)'\n"
          "  - matrix:   Outputs entire K matrix in table form, with E and q "
          "grids printed prior. Legacy. This is form expected by 'dmex'"
          "program. 'dmex' program also expects output in Atomic units.\n"},
     {"E_range",
      "List (2). Minimum, maximum energy transfer (dE), in keV [0.1,0.1]"},
     {"E_steps", "Numer of steps along dE grid (logarithmic grid) [1]"},
     {"q_range",
      "List (2). Minimum, maximum momentum transfer (q), in MeV [0.01,0.01]"},
     {"q_steps", "Number of steps along q grid (logarithmic grid) [1]"},
     {"max_L", "Maximum multipolarity used in exp(iqr) expansion [6]"},
     {"ec_max", "Cut-off (in au) for continuum energy. [inf]"},
     {"Zeff_bound", "Use Zeff for bound state [false]"},
     {"Zeff_cont", "Use Zeff for continuum state. False by default, unless "
                   "Zeff_bound=true, in which case true by default. [false]"},
     {"label", "optional extra label for output files"},
     {"subtract_1", "Replace e^(iqr) -> e^(iqr)-1 [false]"},
     {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
     {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                       "hole-particle interaction) [true]"},
     {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
     {"coupling", "Vector (1), Scalar (g0), Pseudovector (g5), Pseudoscalar "
                  "(g0g5) [Vector]"},
     {"each_state", "bool. If true, will output K(E,q) seperately for each "
                    "(accessible) bound state [false]"},
     {"units", "Units for output: Particle (eV) or Atomic (E_H,1/a0). Use "
               "atomic units for old dmex program. [Particle]"}});

  // Don't run module if requesting help:
  if (input.has_option("help")) {
    return;
  }

  //----------------------------------------------------------------------------

  // Read in energy-deposit/momentum-exchange input options:
  auto [Emin_keV, Emax_keV] = input.get("E_range", std::array{0.1, 0.1});
  auto E_steps = input.get<std::size_t>("E_steps", 1);
  if (E_steps <= 1) {
    E_steps = 1;
    Emax_keV = Emin_keV;
  }
  const auto Emin_au = Emin_keV * UnitConv::Energy_keV_to_au;
  const auto Emax_au =
    Emax_keV < Emin_keV ? Emin_au : Emax_keV * UnitConv::Energy_keV_to_au;

  std::cout << "\nSummary of inputs:\n";
  fmt::print(
    "Energy  : [{:.2f}, {:.2f}] keV  = [{:.1f}, {:.1f}] au, in {} steps\n",
    Emin_keV, Emax_keV, Emin_au, Emax_au, E_steps);

  auto [qmin_MeV, qmax_MeV] = input.get("q_range", std::array{0.01, 0.01});
  auto q_steps = input.get<std::size_t>("q_steps", 1);
  if (q_steps <= 1) {
    q_steps = 1;
    qmax_MeV = qmin_MeV;
  }
  const auto qmin_au = qmin_MeV * UnitConv::Momentum_MeV_to_au;
  const auto qmax_au = qmax_MeV * UnitConv::Momentum_MeV_to_au;
  const auto max_L = input.get("max_L", 6);
  const auto label = input.get("label", std::string{""});

  const auto ec_max = input.get("ec_max", 1.0 / 0.0);

  fmt::print(
    "Momentum: [{:.3f}, {:.3f}] MeV = [{:.1f}, {:.1f}] au, in {} steps\n",
    qmin_MeV, qmax_MeV, qmin_au, qmax_au, q_steps);

  // Set up the E and q grids
  const Grid Egrid({E_steps, Emin_au, Emax_au, 0, GridType::logarithmic});
  const Grid qgrid({q_steps, qmin_au, qmax_au, 0, GridType::logarithmic});

  // Check to see if grid is reasonable for maximum energy:
  Kion::check_radial_grid(std::min(ec_max, Emax_au), qmax_au, wf.grid(),
                          wf.alpha());

  //----------------------------------------------------------------------------
  // Read in and parse options:

  fmt::print("\nMax L   : {}  (multipolarity in e^iqr expansion)\n", max_L);

  if (!label.empty()) {
    fmt::print("Label   : {}\n", label);
  }

  //----------------------------------------------------------------------------

  // Other methods options
  const auto subtract_1 = input.get("subtract_1", false);
  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);

  // Summarise input options
  std::cout << "\nOptions:\n";
  if (subtract_1) {
    std::cout << "Subtract 1: Replacing: e^iqr -> e^iqr - 1\n";
  }
  if (force_rescale) {
    std::cout << "Force rescale: Enforcing V(r) ~ -Z_ion/r at large r\n";
  }
  if (hole_particle) {
    std::cout << "Subtracting HF self-interaction (account for hole-particle "
                 "interaction)\n";
  }
  if (force_orthog) {
    std::cout << "Explicitely enforcing orthogonality between bound and "
                 "continuum states\n";
  }

  // Perform checks, print possible warnings
  if (force_rescale && !(force_orthog || subtract_1)) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Force rescale will ruin orthogonality; suggest use "
               "force_orthog or subtract_1\n");
  }
  if (!force_rescale && !hole_particle && wf.Zion() == 0) {
    fmt::print("\nWarning: Long-range behaviour of V(r) may be incorrect. "
               "Suggest to either "
               "force rescaling, or include hole-particle interaction.\n");
  }
  if (force_rescale && hole_particle) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Should not force rescaling of V(r) if also subtracting hole "
               "particle (self-interaction) potential.\n");
  }

  //----------------------------------------------------------------------------

  // Print core information for convenience:
  // This assumes core is in energy order; always true but not guarenteed.
  // Has no impact on result, just what is printed.
  std::cout << "\n\nCore orbitals:\n";
  std::cout << "        E(au)   E(keV)  N_el\n";
  bool reached_accessible = false;
  if (std::abs(wf.core().front().en()) > Emax_au) {
    std::cout << "------- Inaccessible -------\n";
  }
  for (const auto &Fnk : wf.core()) {
    if (std::abs(Fnk.en()) <= Emax_au && !reached_accessible) {
      reached_accessible = true;
      std::cout << "-------- Accessible --------\n";
    }
    fmt::print("{:3s}  {:8.2f}  {:7.3f}   {}\n", Fnk.shortSymbol(), Fnk.en(),
               Fnk.en() * UnitConv::Energy_au_to_keV, Fnk.num_electrons());
  }
  std::cout << '\n';

  //----------------------------------------------------------------------------

  // DM-electron couplings (allow common variations)
  using qip::ci_compare;
  const auto tcoupling = input.get<std::string>("coupling", "vector");
  const auto coupling =
    ci_compare(tcoupling, "vector")       ? Kion::Coupling::Vector :
    ci_compare(tcoupling, "v")            ? Kion::Coupling::Vector :
    ci_compare(tcoupling, "1")            ? Kion::Coupling::Vector :
    ci_compare(tcoupling, "scalar")       ? Kion::Coupling::Scalar :
    ci_compare(tcoupling, "s")            ? Kion::Coupling::Scalar :
    ci_compare(tcoupling, "g0")           ? Kion::Coupling::Scalar :
    ci_compare(tcoupling, "pseudovector") ? Kion::Coupling::AxialVector :
    ci_compare(tcoupling, "axial")        ? Kion::Coupling::AxialVector :
    ci_compare(tcoupling, "axialvector")  ? Kion::Coupling::AxialVector :
    ci_compare(tcoupling, "a")            ? Kion::Coupling::AxialVector :
    ci_compare(tcoupling, "g5")           ? Kion::Coupling::AxialVector :
    ci_compare(tcoupling, "pseudoscalar") ? Kion::Coupling::PseudoScalar :
    ci_compare(tcoupling, "p")            ? Kion::Coupling::PseudoScalar :
    ci_compare(tcoupling, "g0g5")         ? Kion::Coupling::PseudoScalar :
                                            Kion::Coupling::Error;
  if (coupling == Kion::Coupling::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 212: ");
    fmt::print("Coupling option: `{}' unknown. Options are: vector, scalar, "
               "axialvector, pseudoscalar\n",
               tcoupling);
    return;
  }

  // Construct the effective electron-coupling operator.
  // Uses pointers, since we use polymorphism to swap between coupling types
  std::unique_ptr<DiracOperator::jL> jl = nullptr;
  if (coupling == Kion::Coupling::Vector) {
    jl = std::make_unique<DiracOperator::jL>(wf.grid(), qgrid,
                                             std::size_t(max_L), subtract_1);
  } else if (coupling == Kion::Coupling::Scalar) {
    jl = std::make_unique<DiracOperator::g0jL>(wf.grid(), qgrid,
                                               std::size_t(max_L), subtract_1);
  } else if (coupling == Kion::Coupling::AxialVector) {
    jl = std::make_unique<DiracOperator::ig5jL>(wf.grid(), qgrid,
                                                std::size_t(max_L));
  } else if (coupling == Kion::Coupling::PseudoScalar) {
    jl = std::make_unique<DiracOperator::ig0g5jL>(wf.grid(), qgrid,
                                                  std::size_t(max_L));
  }
  assert(jl != nullptr && "Error in coupling type");
  std::cout << "Operator: " << jl->name() << "\n";

  // For testing: Use H-like with Z-eff for bound and/or continuum:
  const bool use_Zeff_bound = input.get("Zeff_bound", false);
  if (use_Zeff_bound) {
    std::cout << "Using Zeff for bound wavefunction\n";
  }
  const bool use_Zeff_cont = input.get("Zeff_cont", use_Zeff_bound);
  if (use_Zeff_cont) {
    std::cout << "Using Zeff for continuum wavefunction\n";
  }

  //----------------------------------------------------------------------------

  // Create output file-name template

  const auto hf_text =
    use_Zeff_bound ? "Zeff" : HF::parseMethod_short(wf.vHF()->method());

  const std::string coupling_text =
    coupling == Kion::Coupling::Vector       ? "v" :
    coupling == Kion::Coupling::Scalar       ? "s" :
    coupling == Kion::Coupling::AxialVector  ? "a" :
    coupling == Kion::Coupling::PseudoScalar ? "p" :
                                               "???";

  // doesn't include suffix
  std::string oname =
    "K_" + wf.atomicSymbol() + "_" + hf_text + "_" + coupling_text + "_";
  oname += fmt::format("{}_", max_L);
  if (force_rescale)
    oname += "rescale_";
  if (hole_particle)
    oname += "hp_";
  if (force_orthog)
    oname += "orth_";
  if (subtract_1)
    oname += "sub1_";
  if (use_Zeff_cont)
    oname += "Zeffcont_";
  if (label != "")
    oname += label + "_";
  if (oname.back() == '_')
    oname.pop_back();

  //----------------------------------------------------------------------------

  // Output format:
  const bool write_each_state = input.get("each_state", false);

  // Units:
  const auto tunits = input.get<std::string>("units", "Particle");
  auto units = ci_compare(tunits, "eV")       ? Kion::Units::Particle :
               ci_compare(tunits, "particle") ? Kion::Units::Particle :
               ci_compare(tunits, "au")       ? Kion::Units::Atomic :
               ci_compare(tunits, "atomic")   ? Kion::Units::Atomic :
                                                Kion::Units::Error;
  if (units == Kion::Units::Error) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Output units option: `{}' unknown. Options are: Particle, "
               "Atomic. Defaulting to atomic\n",
               tunits);
    units = Kion::Units::Atomic;
  }
  std::cout << "Will use: " << tunits << " units\n";

  //----------------------------------------------------------------------------

  std::cout << "\nCalculating K(E,q) - ionisation factor\n" << std::flush;
  const int num_output_digits = 5;

  // Kion stored K(dE, q)
  LinAlg::Matrix<double> Kion(Egrid.num_points(), qgrid.num_points());

  for (const auto &Fnk : wf.core()) {
    const auto accessible = std::abs(Fnk.en()) < Emax_au;
    if (!accessible)
      continue;
    std::cout << Fnk << ", " << std::flush;

    const auto K_nk = Kion::calculateK_nk(
      wf.vHF(), Fnk, max_L, Egrid, jl.get(), force_rescale, hole_particle,
      force_orthog, use_Zeff_cont, use_Zeff_bound, ec_max);
    if (write_each_state) {
      const auto oname_nk = oname + "_" + Fnk.shortSymbol();
      std::cout << "Written to file: " << oname_nk << "\n";
      write_to_file_xyz(oname_nk + "_xyz.txt", Egrid.r(), qgrid.r(), {"K"},
                        {tcoupling + " (temporal)"}, {K_nk}, units,
                        num_output_digits);
      // Old format: kept for legacy:
      write_to_file_matrix(K_nk, Egrid.r(), qgrid.r(), oname_nk + "_mat.txt",
                           num_output_digits, units);
    }
    Kion += K_nk;
  }

  std::cout << "\nWritten to file: " << oname << "\n";
  write_to_file_xyz(oname + "_xyz.txt", Egrid.r(), qgrid.r(), {"K"},
                    {tcoupling + " (temporal)"}, {Kion}, units,
                    num_output_digits);
  // Old format: kept for legacy:
  write_to_file_matrix(Kion, Egrid.r(), qgrid.r(), oname + "_mat.txt",
                       num_output_digits, units);

  std::cout << '\n';
}

//==============================================================================
void photo(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("photo");

  input.check({
    {"", "For calculating photoionisation cross-section, including comparison "
         "of beyond-dipole approximations. These can also be reconstructed "
         "from the vector form factors from formFactors{} module, which can "
         "serve as a check."},
    {"E_range",
     "List (2). Minimum, maximum energy transfer (dE), in eV [10, 1000]"},
    {"E_steps", "Numer of steps along dE grid (logarithmic grid) [50]"},
    {"E_threshold", "Numer of extra E steps to add in -15% range on either side"
                    " of each threshold. If <2, will add no new points [0]"},
    {"E_extra", "List (comma separated) extra energies (in eV) to add 10 "
                "points around. Useful for specific regions we want more "
                "resolution in."},
    {"oname", "oname"},
    {"ec_max", "Cut-off (in au) for continuum energy. [1e99]"},
    {"K_minmax", "List (2). Minimum, maximum K [1, 1]"},
    {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
    {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                      "hole-particle interaction) [true]"},
    {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
  });
  if (input.has_option("help")) {
    return;
  }

  // Set up energy grid:
  auto [Emin_eV, Emax_eV] = input.get("E_range", std::array{10.0, 1000.0});
  auto E_steps = input.get<std::size_t>("E_steps", 50);
  auto E_threshold = input.get<std::size_t>("E_threshold", 0);
  if (E_steps <= 1) {
    E_steps = 1;
    Emax_eV = Emin_eV;
  }
  // Convert to atomic units for calculations:
  const auto Emin_au = Emin_eV / PhysConst::Hartree_eV;
  const auto Emax_au =
    Emax_eV < Emin_eV ? Emin_au : Emax_eV / PhysConst::Hartree_eV;

  // const Grid Egrid({E_steps, Emin_au, Emax_au, 0, GridType::logarithmic});
  // auto energies = Egrid.r();

  // Instead of using "grid" - specificly add extra points around
  auto energies = qip::logarithmic_range(Emin_au, Emax_au, E_steps);

  std::cout << "\nCore ionisation energies, in MeV\n";
  for (const auto &Fc : wf.core()) {
    fmt::print("{:3} : {:.4e}\n", Fc.shortSymbol(),
               -1 * Fc.en() * PhysConst::Hartree_eV / 1.0e6);
  }
  std::cout << "\n";

  // Add extra energy points near thresholds:
  if (E_threshold > 1) {
    for (const auto &Fc : wf.core()) {

      // just below thresholds:
      const auto extra1 =
        qip::uniform_range(-0.85 * Fc.en(), -0.999 * Fc.en(), E_threshold);

      // Just above thresholds (note: careful, since hard to solve
      // Dirac equation for cntm states with very small energy)
      const auto e0 = 0.01; // smallest energy can calculate well for cntm
      const auto extra2 =
        qip::uniform_range(-Fc.en() + e0, 1.15 * (-Fc.en() + e0), E_threshold);
      energies = qip::merge(energies, extra1, extra2);
    }
  }

  // Add extra points around specific energies
  const auto E_extra = input.get("E_extra", std::vector<double>{});
  for (const auto &Em_eV : E_extra) {
    const auto Em = Em_eV / PhysConst::Hartree_eV;
    const auto extra3 = qip::uniform_range(0.8 * Em, 1.2 * Em, 10);
    energies = qip::merge(energies, extra3);
  }

  // If added extra points, sort list:
  if (E_threshold > 1 || E_extra.size() > 0) {
    std::sort(energies.begin(), energies.end());
  }

  // "cut-off"/ceiling energy for continuum. Bad idea?
  const auto ec_max = input.get("ec_max", 1 / 0.0);

  std::cout << "\nSummary of inputs:\n";
  fmt::print(
    "Energy  : [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n\n",
    energies.front() * PhysConst::Hartree_eV,
    energies.back() * PhysConst::Hartree_eV, energies.front(), energies.back(),
    energies.size());

  const auto [Kmin, Kmax] = input.get("K_minmax", std::array{1, 1});
  // Error with structured bindings + clang++ with OpenMP:
  // error: capturing a structured binding is not yet supported in OpenMP
  // So, workaround: make local copy
  const auto Kmin_ = Kmin;
  const auto Kmax_ = Kmax;

  const auto label = input.get("label", std::string{""});

  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);

  auto oname = input.get("oname", std::string{"photoel-out.txt"});

  // "full" dipole operator
  const auto E1 = DiracOperator::E1(wf.grid());
  const auto E2 = DiracOperator::Ek(wf.grid(), 2);

  auto M1nr = DiracOperator::M1nr();

  // Store per-omega results in an array, then write to file after the parallel
  // loop.
  // Columns (per row, matching the existing file format):
  // Q_E1, Q_M1, Q_M1_nr, Q_E, Q_E_len, Q_M, Q_EM(=Q_E+Q_M), Q_E2, Q_Ek2, Q_Mk1
  // omega itself is recovered from energies[i_omega] at write time.
  constexpr std::size_t N_cols = 10;
  std::vector<std::array<double, N_cols>> results(energies.size());

  qip::ProgressBar prog(int(energies.size()));
#pragma omp parallel for schedule(dynamic)
  for (std::size_t i_omega = 0; i_omega < energies.size(); ++i_omega) {
    const auto omega = energies[i_omega];

    // Conversion factor from dimensionless Q absorption form factor to sigma
    const auto Ksigma = 4.0 * M_PI * M_PI * PhysConst::alpha *
                        PhysConst::aB_cm * PhysConst::aB_cm * omega;

    // Regular (length) E1 operator
    double Q_E1 = 0.0;
    // "Regular" M1
    double Q_M1 = 0.0;
    // Non-relativistic M1 operator
    double Q_M1_nr = 0.0;
    // Mk at k=1 (compare to regular M1)
    double Q_Mk1 = 0.0;
    // Ek at K=2 (compare to regular E2)
    double Q_Ek2 = 0.0;
    // Regular E2 (length)
    double Q_E2 = 0.0;
    // Full electric multipole
    double Q_E = 0.0;
    // Full magnetic multipole
    double Q_M = 0.0;
    // Full electric multipole (length form)
    double Q_E_len = 0.0;

    for (int k = Kmin_; k <= Kmax_; ++k) {

      // Electric, magnetic parts
      const auto Ek = DiracOperator::VEk(wf.grid(), k, omega);
      const auto Mk = DiracOperator::VMk(wf.grid(), k, omega);
      // Magnetic dipole
      const auto M1 = DiracOperator::M1(wf.grid(), PhysConst::alpha, omega);
      // "Length" form - for tests only
      const auto Ek_len = DiracOperator::VEk_Len(wf.grid(), k, omega);

      for (const auto &Fa : wf.core()) {
        const auto ec = omega + Fa.en();
        if (ec < 0.0 || ec > ec_max)
          continue;

        const int l = Fa.l();
        const int lc_max = l + k + 1;
        const int lc_min = std::max(l - k - 1, 0);

        ContinuumOrbitals cntm(wf.vHF());
        cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa, force_rescale,
                              hole_particle, force_orthog);

        for (const auto &Fe : cntm.orbitals) {

          const auto q = PhysConst::alpha * omega;

          const auto tkp1 = 2.0 * k + 1.0;
          const auto pol_av = 1.0 / 2.0;
          const auto f_Q =
            tkp1 * pol_av / qip::pow(PhysConst::alpha * omega, 2);

          // check!
          const auto f_Q_E1 = 1.0 / 3.0;
          const auto f_Q_M1 = 1.0 / 3.0 * qip::pow(PhysConst::muB_CGS, 2);

          if (k == 1) {
            Q_E1 += f_Q_E1 * qip::pow(E1.reducedME(Fe, Fa), 2);
            Q_M1 += f_Q_M1 * qip::pow(M1.reducedME(Fe, Fa), 2);
            Q_Mk1 += f_Q * qip::pow(Mk.reducedME(Fe, Fa), 2);
            Q_M1_nr += f_Q_M1 * qip::pow(M1nr.reducedME(Fe, Fa), 2);
          }
          if (k == 2) {
            // test with "actual" E2 as well!
            Q_Ek2 += f_Q * qip::pow(Ek.reducedME(Fe, Fa), 2);

            Q_E2 += f_Q_E1 * qip::pow(E2.reducedME(Fe, Fa), 2) / 20 * q * q;
          }

          Q_E += f_Q * qip::pow(Ek.reducedME(Fe, Fa), 2);
          Q_M += f_Q * qip::pow(Mk.reducedME(Fe, Fa), 2);

          Q_E_len += f_Q * qip::pow(Ek_len.reducedME(Fe, Fa), 2);
        }
      }
    }

    results[i_omega] = {Ksigma * Q_E1,        //
                        Ksigma * Q_M1,        //
                        Ksigma * Q_M1_nr,     //
                        Ksigma * Q_E,         //
                        Ksigma * Q_E_len,     //
                        Ksigma * Q_M,         //
                        Ksigma * (Q_E + Q_M), //
                        Ksigma * Q_E2,        //
                        Ksigma * Q_Ek2,       //
                        Ksigma * Q_Mk1};
    prog.update();
  }

  // Sequential write after the parallel loop completes.
  std::ofstream out_file(oname);
  out_file << "# Photoelectric effect::\n"
           << "# Cross section (cm^2):\n"
           << "# Columns:\n"
           << "# omega_MeV      : photon energy (MeV)\n"
           << "# sigma_E1       : E1 (length) dipole cross section\n"
           << "# sigma_M1       : M1 dipole\n"
           << "# sigma_M1_nr    : M1 non-relativistic\n"
           << "# sigma_E        : electric multipole (velocity)\n"
           << "# sigma_E_len    : electric multipole (length) (Ek) all K\n"
           << "# sigma_M        : magnetic multipole (Mk) all K\n"
           << "# sigma_EM       : sigma_E + sigma_M (total multipole)\n"
           << "# sigma_E2       : E2 (length)\n"
           << "# sigma_Ek2      : Ek at K=2\n"
           << "# sigma_Mk1      : Mk at K=1\n";
  out_file << "#\n"
           << "omega_MeV  sigma_E1  sigma_M1  sigma_M1_nr  sigma_E  "
              "sigma_E_len  sigma_M  "
              "sigma_EM  sigma_E2  sigma_Ek2  sigma_Mk1\n";

  for (std::size_t i_omega = 0; i_omega < energies.size(); ++i_omega) {
    const auto &r = results[i_omega];
    out_file << energies[i_omega] * PhysConst::Hartree_eV / 1e6 // omega (MeV)
             << " " << r[0] // s_E1      : E1 (length) dipole
             << " " << r[1] // s_M1      : M1 dipole
             << " " << r[2] // s_M1_nr   : M1 non-relativistic
             << " " << r[3] // s_E       : electric multipole (velocity)
             << " " << r[4] // s_E_len   : electric multipole (length)
             << " " << r[5] // s_M       : magnetic multipole
             << " " << r[6] // s_EM      : Q_E + Q_M (total multipole)
             << " " << r[7] // s_E2      : E2 (length)
             << " " << r[8] // s_Ek2     : Ek at K=2
             << " " << r[9] // s_Mk1     : Mk at K=1
             << "\n";
  }
}

//==============================================================================
void photoRPA(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("photoRPA");

  input.check({
    {"", "E1 photoionisation cross-section, with core polarisation (RPA) "
         "included via the continuum TDHF method (Johnson RRPA), in BOTH "
         "length and velocity gauges (RRPA-L = RRPA-V is the gauge test). "
         "Columns per gauge: bare (V^{N-1} tree), fully-iterated RPA "
         "(standing-wave), and unitarised RPA, plus convergence "
         "diagnostics. The photoelectron is always treated in the V^{N-1} "
         "(residual ion) potential, orthogonalised to the core."},
    {"E_range",
     "List (2). Minimum, maximum photon energy (omega), in eV [10, 1000]"},
    {"E_steps", "Number of steps along omega grid (logarithmic) [50]"},
    {"ec_max", "Cut-off (in au) for continuum energy. [1e99]"},
    {"max_its", "Maximum RPA (TDHF) iterations per omega [60]"},
    {"eps", "Convergence target for the RPA iterations [1e-5]"},
    {"eta", "Damping factor for the RPA iterations (larger = heavier "
            "damping; helps near autoionising resonances) [0.4]"},
    {"unitarise", "Unitarise the RPA amplitudes (on-shell rescattering, "
                  "Johnson 1979 appendix): the standing-wave amplitudes "
                  "have real poles (e.g. just above deep thresholds) where "
                  "the physical cross-section is finite. Costs one extra "
                  "RPA-type solve per open channel per energy. [true]"},
    {"each_shell", "Add per-shell RPA cross-section columns [false]"},
    {"rpa_max", "Maximum photon energy (eV) at which the RPA is solved. Core "
                "polarisation is a low-energy effect: above a few hundred eV "
                "the correction is well under a percent, while the solve is "
                "at its most expensive (every shell open). Above this, the "
                "RPA columns repeat the bare values. Set 0 (or negative) for "
                "no limit [1000]"},
    {"oname", "Output file name [photoRPA-out.txt]"},
  });
  if (input.has_option("help")) {
    return;
  }

  // Set up energy grid:
  const auto [Emin_eV, Emax_eV] =
    input.get("E_range", std::array{10.0, 1000.0});
  const auto E_steps = input.get<std::size_t>("E_steps", 50);
  const auto Emin_au = Emin_eV / PhysConst::Hartree_eV;
  const auto Emax_au =
    Emax_eV < Emin_eV ? Emin_au : Emax_eV / PhysConst::Hartree_eV;
  const auto energies = qip::logarithmic_range(Emin_au, Emax_au, E_steps);

  const auto ec_max = input.get("ec_max", 1 / 0.0);
  const auto max_its = input.get("max_its", 60);
  const auto eps_target = input.get("eps", 1.0e-5);
  const auto eta = input.get("eta", 0.4);
  const auto unitarise = input.get("unitarise", true);
  const auto each_shell = input.get("each_shell", false);
  const auto rpa_max_eV = input.get("rpa_max", 1000.0);
  const auto rpa_max_au =
    rpa_max_eV > 0.0 ? rpa_max_eV / PhysConst::Hartree_eV : 1.0 / 0.0;
  const auto oname = input.get("oname", std::string{"photoRPA-out.txt"});

  fmt::print("\nomega : [{:.1f}, {:.1f}] eV, {} steps (logarithmic)\n", Emin_eV,
             Emax_eV, energies.size());
  if (rpa_max_au < energies.back()) {
    fmt::print("RPA solved up to {:.0f} eV; bare above\n", rpa_max_eV);
  }

  // Grid-resolution guidance for the requested range: the channel solves
  // and the (c, K) extraction need ~40 points per wavelength on the main
  // grid (the pair itself is fine-grid solved past that; see
  // solveContinuumForward). Report the omega up to which the current grid
  // meets this everywhere, and what the full range would need. Purely
  // informative: higher omegas run on the reduced extraction window, and
  // failures are flagged by the K = pi*D warnings / KpiD_dev column.
  if (!wf.core().empty()) {
    const auto e_least = std::max_element(wf.core().cbegin(), wf.core().cend(),
                                          [](const auto &a, const auto &b) {
                                            return a.en() < b.en();
                                          })
                           ->en();
    const auto en_top = energies.back() + e_least;
    const auto N_ppw = 40.0;
    if (en_top > 0.0) {
      const auto req =
        DiracODE::RequiredContinuumGrid(en_top, wf.grid(), N_ppw, wf.alpha());
      if (req.num_points > wf.grid().num_points()) {
        // largest en with dr_max <= lambda(en)/N_ppw, from
        // k^2 = en*(2 + alpha^2 en):
        const auto dr_max = wf.grid().drdu().back() * wf.grid().du();
        const auto k_ok = 2.0 * M_PI / (N_ppw * dr_max);
        const auto al2 = wf.alpha() * wf.alpha();
        const auto en_ok = (std::sqrt(1.0 + al2 * k_ok * k_ok) - 1.0) / al2;
        fmt::print("\nNote: this grid resolves the continuum solves ({:.0f} "
                   "pts/wavelength) up to omega ~ {:.0f} eV; the full range "
                   "would need num_points ~ {} (or b ~ {:.3f}). Above that, "
                   "watch the K = pi*D warnings (KpiD_dev column).\n",
                   N_ppw, (en_ok - e_least) * PhysConst::Hartree_eV,
                   req.num_points, req.b);
      }
    }
  }

  const auto E1 = DiracOperator::E1(wf.grid());

  // Per-omega results: {sigma_bare, sigma_rpa, sigma_rpaU} for E1 (length)
  // then E1v (velocity), eps, its, KpiD_dev, Kbar_asym (diagnostics = worst
  // over the two gauges), then optional per-shell (unitarised, length) sigmas
  const auto n_shells = wf.core().size();
  const auto n_cols = 10ul + (each_shell ? n_shells : 0ul);
  std::vector<std::vector<double>> results(energies.size(),
                                           std::vector<double>(n_cols, 0.0));

  // Scan-wide diagnostics: collected here and reported once at the end,
  // rather than warned on per omega.
  double KpiD_worst = 0.0;
  double KpiD_worst_at = 0.0;
  std::string KpiD_worst_ch{};
  std::size_t n_excluded_omega = 0;

  fmt::print("\n{:>9s}  {:>10s}  {:>10s}{}  {:>8s} {:>4s}  {:>8s}\n",
             "omega/eV", "bare", "RPA", unitarise ? "       RPA_U" : "", "eps",
             "its", "KpiD");

  // Energies run serially; threads are used inside each solve (the channel
  // tasks, and the rescattering seed columns). With unitarisation the cost
  // concentrates in the few high-omega points (many open channels = many
  // seeded solves), which parallelise well over the seed columns --
  // omega-parallel would leave the most expensive omega serial on one
  // thread (and hold an Anderson history per thread).
  for (std::size_t i_omega = 0; i_omega < energies.size(); ++i_omega) {
    const auto omega = energies[i_omega];

    // Below every ionisation threshold (or every open shell beyond ec_max):
    // no contributing channels, sigma = 0 -- nothing to solve. (Same shell
    // condition as the matrix-element loop below.)
    const auto any_open =
      std::any_of(wf.core().cbegin(), wf.core().cend(), [&](const auto &Fa) {
        const auto ec = omega + Fa.en();
        return ec >= 0.0 && ec <= ec_max;
      });
    if (!any_open) {
      continue;
    }

    // Above rpa_max the RPA is not solved: the columns repeat the bare
    // values (see the rpa_max option).
    const bool do_rpa = omega <= rpa_max_au;
    const bool do_unitarise = unitarise && do_rpa;

    // Conversion: (1/3) sum |<e||E1||a>|^2 -> cross-section (cm^2)
    const auto Ksigma = 4.0 * M_PI * M_PI * PhysConst::alpha *
                        PhysConst::aB_cm * PhysConst::aB_cm * omega / 3.0;

    // Velocity-form operators: frequency-dependent, so per-omega instances,
    // updated to (+omega, -omega) for the (t_+, t_-) pair.
    auto E1v = DiracOperator::E1v(wf.alpha(), omega);
    auto E1v_minus = DiracOperator::E1v(wf.alpha(), -omega);

    // Solve the RPA (continuum TDHF) at this omega, both gauges:
    auto rpa_L = ExternalField::TDHFcntm(&E1, wf.vHF());
    auto rpa_V = ExternalField::TDHFcntm(&E1v, wf.vHF(), &E1v_minus);
    if (do_rpa) {
      rpa_L.eps_target() = eps_target;
      rpa_L.set_eta(eta);
      rpa_L.solve_core(omega, max_its, false);

      rpa_V.eps_target() = eps_target;
      rpa_V.set_eta(eta);
      rpa_V.solve_core(omega, max_its, false);
    }

    // Matrix elements: D = D0 + <e|dV|a> (standing-wave), and the
    // unitarised (physical) amplitudes D_phys = |A|/pi [A = (1-i*Kbar)^{-1}
    // pi*D: finite through the standing-wave poles above thresholds]:
    double sigma_E1 = 0.0;
    double sigma_E1_rpa = 0.0;
    double sigma_E1_rpaU = 0.0;
    double sigma_E1v = 0.0;
    double sigma_E1v_rpa = 0.0;
    double sigma_E1v_rpaU = 0.0;
    for (std::size_t i_shell = 0; i_shell < n_shells; ++i_shell) {
      const auto &Fa = wf.core()[i_shell];
      const auto ec = omega + Fa.en();
      if (ec < 0.0 || ec > ec_max)
        continue;
      ContinuumOrbitals cntm(wf.vHF());
      cntm.solveContinuumHF(ec, std::max(Fa.l() - 1, 0), Fa.l() + 1, &Fa, false,
                            true, true);
      for (const auto &Fe : cntm.orbitals) {
        if (E1.isZero(Fe, Fa) || Fe.norm2() == 0.0)
          continue;
        const auto D0 = E1.reducedME(Fe, Fa);
        const auto D = do_rpa ? D0 + rpa_L.dV_cntm(Fe, Fa) : D0;
        const auto D0v = E1v.reducedME(Fe, Fa);
        const auto Dv = do_rpa ? D0v + rpa_V.dV_cntm(Fe, Fa) : D0v;
        sigma_E1 += Ksigma * D0 * D0;
        sigma_E1_rpa += Ksigma * D * D;
        sigma_E1v += Ksigma * D0v * D0v;
        sigma_E1v_rpa += Ksigma * Dv * Dv;
        // Unitarised: the physical amplitudes. Without the RPA solve there
        // is nothing to rescatter, so the bare value stands.
        const auto DU = do_unitarise ? rpa_L.D_phys(Fa, Fe.kappa()) : D0;
        const auto DUv = do_unitarise ? rpa_V.D_phys(Fa, Fe.kappa()) : D0v;
        if (unitarise) {
          sigma_E1_rpaU += Ksigma * DU * DU;
          sigma_E1v_rpaU += Ksigma * DUv * DUv;
          if (each_shell) {
            results[i_omega][10 + i_shell] += Ksigma * DU * DU;
          }
        }
      }
    }

    // Convergence / internal-consistency diagnostics, worst over gauges
    // [KpiD: the K = pi*D identity, computed by solve_core]
    const auto KpiD_dev = std::max(rpa_L.KpiD_dev(), rpa_V.KpiD_dev());
    if (KpiD_dev > KpiD_worst) {
      KpiD_worst = KpiD_dev;
      KpiD_worst_at = omega;
      KpiD_worst_ch = rpa_L.KpiD_dev() >= rpa_V.KpiD_dev() ?
                        rpa_L.KpiD_worst_channel() :
                        rpa_V.KpiD_worst_channel();
    }
    for (const auto *rpa : {&rpa_L, &rpa_V}) {
      if (!rpa->excluded_channels().empty()) {
        n_excluded_omega++;
        break;
      }
    }

    results[i_omega][0] = sigma_E1;
    results[i_omega][1] = sigma_E1_rpa;
    results[i_omega][2] = sigma_E1_rpaU;
    results[i_omega][3] = sigma_E1v;
    results[i_omega][4] = sigma_E1v_rpa;
    results[i_omega][5] = sigma_E1v_rpaU;
    // Convergence of the RPA solve; zero where it was not solved (above
    // rpa_max), rather than the unsolved default.
    results[i_omega][6] =
      do_rpa ? std::max(rpa_L.last_eps(), rpa_V.last_eps()) : 0.0;
    results[i_omega][7] =
      do_rpa ? std::max(rpa_L.last_its(), rpa_V.last_its()) : 0.0;
    results[i_omega][8] = KpiD_dev;
    results[i_omega][9] =
      std::max(rpa_L.rescattering_asymmetry(), rpa_V.rescattering_asymmetry());

    fmt::print("{:9.2f}  {:10.3e}  {:10.3e}{}  {:8.1e} {:4.0f}  {:8.1e}{}\n",
               omega * PhysConst::Hartree_eV, sigma_E1, sigma_E1_rpa,
               unitarise ? fmt::format("  {:10.3e}", sigma_E1_rpaU) : "",
               results[i_omega][6], results[i_omega][7], KpiD_dev,
               do_rpa ? "" : "  (bare)");
    std::cout << std::flush;
  }

  // Scan-wide quality report: printed once, in place of per-omega warnings.
  if (KpiD_worst > 0.05) {
    fmt::print("\nK = pi*D worst: {:.1e} [{}] at w = {:.1f} eV. Above ~5% the "
               "channel solve is unreliable on this grid (increase "
               "num_points); check the KpiD_dev column.\n",
               KpiD_worst, KpiD_worst_ch,
               KpiD_worst_at * PhysConst::Hartree_eV);
  }
  if (n_excluded_omega != 0) {
    fmt::print("\n{} of {} energies had open channels dropped (continuum solve "
               "returned zero: grid resolves nothing at that energy). Increase "
               "num_points to include them.\n",
               n_excluded_omega, energies.size());
  }

  std::ofstream out_file(oname);
  out_file << "# E1 photoionisation cross-section (cm^2), with RPA; "
              "length (E1) and velocity (E1v) gauges\n"
           << "# Columns:\n"
           << "# omega_eV       : photon energy (eV)\n"
           << "# sigma_E1       : bare (V^{N-1} tree, no RPA), length\n"
           << "# sigma_E1_rpa   : full (iterated) RPA, standing-wave "
              "amplitudes (poles above thresholds!), length\n"
           << "# sigma_E1_rpaU  : full RPA, unitarised (physical) "
              "amplitudes, length\n"
           << "# sigma_E1v      : bare, velocity\n"
           << "# sigma_E1v_rpa  : full RPA, standing-wave, velocity\n"
           << "# sigma_E1v_rpaU : full RPA, unitarised, velocity\n"
           << "# rpa_eps        : RPA convergence achieved (worst gauge)\n"
           << "# rpa_its        : RPA iterations used (worst gauge)\n"
           << "# KpiD_dev       : worst |K - pi*D| / max|K| internal "
              "consistency\n"
           << "# Kbar_asym      : rescattering-matrix asymmetry "
              "(consistency)\n";
  if (rpa_max_au < energies.back()) {
    out_file << "# nb: RPA solved only up to " << rpa_max_eV
             << " eV; above that the rpa/rpaU columns repeat the bare\n";
  }
  if (each_shell) {
    out_file << "# sigma_<shell>  : per-shell unitarised-RPA cross-section "
                "(length)\n";
  }
  out_file << "#\n"
           << "omega_eV  sigma_E1  sigma_E1_rpa  sigma_E1_rpaU  sigma_E1v  "
              "sigma_E1v_rpa  sigma_E1v_rpaU  rpa_eps  rpa_its  KpiD_dev  "
              "Kbar_asym";
  if (each_shell) {
    for (const auto &Fc : wf.core()) {
      out_file << "  sigma_" << Fc.shortSymbol();
    }
  }
  out_file << "\n";

  for (std::size_t i_omega = 0; i_omega < energies.size(); ++i_omega) {
    const auto &r = results[i_omega];
    out_file << energies[i_omega] * PhysConst::Hartree_eV;
    for (std::size_t ic = 0; ic < r.size(); ++ic) {
      out_file << " " << r[ic];
    }
    out_file << "\n";
  }
  std::cout << "Written to: " << oname << "\n";
}

//==============================================================================
void formFactors(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("formFactors");

  input.check(
    {{"",
      "Calculates generalised atomic ionisation (scattering) form factors.\n"
      "See paper <arXiv:xxxx.xxxxx> for details.\n"
      "These depend on energy exchange, E, and momentum exchange q. Be careful "
      "to distinguish between energy exchange, E, and final-state energy e_f.\n"
      "Output file will be in the 'XYZ' format, with each requested "
      "formfactor in new column (with header)\n"
      "    E1 q1 K1(E1,q1) K2(E1,q1) ... KN(E1,q1)\n"
      "    E1 q2 K1(E1,q2) ...\n"
      "    ....               \n"
      "    E1 qM K1(E1,qM) ...\n"
      "    E2 q1 K1(E2,q1) ...\n"
      "    ....               \n"
      "    EN qM K1(E1,q1) ... KN(E1,q1)\n"
      "q and E are given in eV; factors K are dimensionless.\n"
      "Output filename will be in form, e.g., "
      "Xe0_method_kmin-kmax_VASP.txt\n\n"},
     {"E_range", "List (2), comma separated. Minimum, maximum energy transfer "
                 "(E), in eV [10.0, 1.0e4]"},
     {"E_steps", "Numer of steps along dE grid (logarithmic) [1]"},
     {"E_set", "List: set of specific E values to calculate for (in eV) - "
               "Will override E_range/E_set if set"},
     {"q_range", "List (2). Minimum, maximum momentum transfer (q), in eV "
                 "(hbar=c=1). For reference, 1/a0 ~ 3730 eV. [1.0e4, 1.0e7]"},
     {"q_steps", "Numer of steps along q grid (logarithmic) [1]"},
     {"q_set", "List: set of specific q values to calculate for (in eV) - "
               "Will override q_range/q_set if set"},
     {"diagonal",
      "true/false. If true, forces momentum exchange to be equal (in hbar=c=1 "
      "units) to energy exchange (e.g., for absorption of massless particle). "
      "Overrides q_steps/q_range/q_set if true. [false]"},
     {"label", "Extra label appended to output file name (not usually nedded)"},
     {"operators",
      "String. Which factors to calculate. Any combination of:\n"
      "   'V' (vector), \n"
      "   'A' (axial-vector), \n"
      "   'S' (scalar), \n"
      "   'P' (pseudoscalar).\n"
      "If 'V' (e.g.), output will include: temporal V0, electric VE, magnetic "
      "VM, longitudinal VL, and cross X terms."
      "Cross-terms X,Y,Z will be calculated automatically depending on input."},
     {"temporal_only", "If only non-relativistic scattering is required, save "
                       "time by not calculating the spatial terms. [false]"},
     {"K_minmax", "List (2). Minimum, maximum multipolarity K [0, 6]"},
     {"lc_minmax", "List (2). Minimum and maximum orbital quantum number l to "
                   "include in continuum states. By "
                   "default, will include all allowed values (limitted by "
                   "K_max). This should usually be left blank"},
     {"ec_minmax",
      "List of floats (2). Minimum and maximum continuum state energy to "
      "include (eV). "
      "Minimum energy can be used, e.g., to model recombination of "
      "low-E ionised electron. Maximum energy may be used to simplify "
      "calculations, in cases where lower shells dominate (and high-energy "
      "electron ionised from high shell doesn't contribute to cross section). "
      "Both must be used with care, and should not normally be different from "
      "default. By default: minimum=0.0, maximum=infinity. For infinity, input "
      "very large number. [0.0, 1e99]"},
     {"each_state",
      "true/false. If true, will output a separate K file for each "
      "(accessible) bound state. Output file will have same name as total, but "
      "with short-form orbital label appended (i.e., 3p-=3p_1/2, "
      "3p+=3p_3/2). [false]"},
     {"low_q", "Explicitly use low-q form of operators. These are only valid "
               "at low q only, and are used for numerical tests. (nb: All are "
               "zero for K>2.) [false]"},
     {"force_rescale", "Rescale atomic potential V(r) at large r when solving "
                       "continuum orbitals. Should be false for local "
                       "potentials or if hole-particle is included. [false]"},
     {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                       "hole-particle interaction) [true]"},
     {"force_orthog", "Enforce orthogonality of the continuum orbitals [true]"},
     {"method",
      "Method for bound and continuum states: HF (standard), Zeff (H-like, "
      "solved numerically with DiracODE), ZeffAnalytic (H-like, exact "
      "analytic functions; requires FLINT). [HF; Zeff if Zeff option is set]"},
     {"Zeff",
      "Effective charge for the Zeff/ZeffAnalytic methods. If set (to "
      "anything), the default method becomes Zeff. Set to 'true' or <=0 to "
      "use 'real' Zeff = n*sqrt(-2*en) from each binding energy; "
      "set to a positive value to use that constant Zeff for all states."}});
  if (input.has_option("help")) {
    return;
  }

  // Read in energy grid. Input is in eV:
  auto [Emin_eV, Emax_eV] = input.get("E_range", std::array{10.0, 1.0e3});
  auto E_steps = input.get<std::size_t>("E_steps", 1);
  if (E_steps <= 1 || Emax_eV <= Emin_eV) {
    E_steps = 1;
    Emax_eV = Emin_eV;
  }

  // Option for specific set of energies
  const auto E_set_eV = input.get("E_set", std::vector<double>{});
  if (!E_set_eV.empty()) {
    // override above
    Emin_eV = E_set_eV.front();
    Emax_eV = E_set_eV.back();
    E_steps = E_set_eV.size();
  }

  // Convert to atomic units (from keV):
  const auto Emin_au = Emin_eV / PhysConst::Hartree_eV;
  const auto Emax_au = Emax_eV / PhysConst::Hartree_eV;

  // Form the actual energy grid (in atomic units)
  using namespace qip::overloads;
  const auto Egrid = E_set_eV.empty() ?
                       qip::logarithmic_range(Emin_au, Emax_au, E_steps) :
                       E_set_eV / PhysConst::Hartree_eV;

  //----------------------------------------------------------------------------

  // Momentum-transfer range (in eV)
  const auto diagonal_Eq = input.get("diagonal", false);
  auto [qmin_eV, qmax_eV] = input.get("q_range", std::array{1.0, 1.0e4});
  auto q_steps = diagonal_Eq ? 1 : input.get<std::size_t>("q_steps", 1);
  if (q_steps <= 1 || qmax_eV <= qmin_eV) {
    q_steps = 1;
    qmax_eV = qmin_eV;
  }

  // Option for specific set of momentum exchange
  const auto q_set_eV = input.get("q_set", std::vector<double>{});
  if (!q_set_eV.empty() && !diagonal_Eq) {
    // override above
    qmin_eV = q_set_eV.front();
    qmax_eV = q_set_eV.back();
    q_steps = q_set_eV.size();
  }

  // Convert momentum from keV to atomic units.
  const auto q_min = qmin_eV * UnitConv::Momentum_eV_to_au;
  const auto q_max = qmax_eV * UnitConv::Momentum_eV_to_au;

  // Set up the q grid
  const auto qgrid = qip::logarithmic_range(q_min, q_max, q_steps);

  const auto alpha_ratio = wf.alpha() / PhysConst::alpha;
  if (std::abs(alpha_ratio - 1) > 0.001) {
    fmt::print("\nEffective speed of light: c_eff = {} c\n", 1.0 / alpha_ratio);
    fmt::print("(dα^2 = {})", wf.dalpha2());
  }

  std::cout << "\nEnergy/Momentum exchange grids:\n";
  fmt::print(
    "Energy  : [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n",
    Emin_eV, Emax_eV, Emin_au, Emax_au, E_steps);
  if (diagonal_Eq) {
    std::cout << "Fixing momentum exchange equal to energy exchange: qc = E\n";
  } else {
    fmt::print(
      "Momentum: [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n",
      qmin_eV, qmax_eV, q_min, q_max, q_steps);
  }

  //----------------------------------------------------------------------------

  // Check to see if grid is reasonable for maximum energy/momentum.
  // Just prints warning to screen if not the case
  Kion::check_radial_grid(Emax_au, q_max, wf.grid(), wf.alpha());

  //----------------------------------------------------------------------------

  // Print core information for convenience:
  // This assumes core is in energy order; always true but not guarenteed.
  // Has no impact on result, just what is printed.
  std::cout << "\nBound orbitals:\n";
  std::cout << "state  κ l   j       E(au)      E(eV)  N_el\n";
  bool reached_accessible = false;
  for (const auto &Fnk : wf.core()) {
    if (std::abs(Fnk.en()) <= Emax_au && !reached_accessible) {
      reached_accessible = true;
      std::cout << "-------------------------------------------\n";
    }
    fmt::print("{:3s}   {:+2} {:1} {:1}/2   {:9.3f}  {:9.2f}     {}\n",
               Fnk.shortSymbol(), Fnk.kappa(), Fnk.l(), Fnk.twoj(), Fnk.en(),
               Fnk.en() * UnitConv::Energy_au_to_eV, Fnk.num_electrons());
  }

  const auto each_state = input.get("each_state", false);
  // if (each_state)

  //----------------------------------------------------------------------------

  // Which operators to calculate
  // Case insensitive, and only check first letter
  const auto operators = input.get("operators", std::string{"V"});
  const auto temporal_only = input.get("temporal_only", false);
  const auto spatialQ = !temporal_only;

  bool vectorQ{false}, axialQ{false}, scalarQ{false}, pseudoscalarQ{false};

  // Interference terms:
  // vector-vector spatial-temporal interference: auto-include with vector
  bool XvvQ{false};
  // axial-axial spatial-temporal interference: auto-include with axial
  bool YaaQ{false};
  // axial-vector spatial interference: auto-include if axial and vector
  bool ZvaQ{false};

  // Check all characters in input string for v,a,s,p [case insensitive]
  for (auto &w : operators) {
    switch (std::tolower(w)) {
    case 'v':
      vectorQ = true;
      if (spatialQ)
        XvvQ = true;
      break;
    case 'a':
      axialQ = true;
      if (spatialQ)
        YaaQ = true;
      break;
    case 's':
      scalarQ = true;
      break;
    case 'p':
      pseudoscalarQ = true;
      break;
    }
  }
  // automatically include if doing both
  ZvaQ = (vectorQ && axialQ && spatialQ);

  fmt::print("\nComputing operators: "
             "{}{}{}{}\n",
             vectorQ ? "Vector; " : "", axialQ ? "Axial; " : "",
             scalarQ ? "Scalar; " : "", pseudoscalarQ ? "Pseudoscalar; " : "");
  if (XvvQ) {
    std::cout << "- and spatial-temporal vector interference term X\n";
  }
  if (YaaQ) {
    std::cout << "- and spatial-temporal axial interference term Y\n";
  }
  if (ZvaQ) {
    std::cout << "- and spatial vector-axial interference term Z\n";
  }
  if (!spatialQ) {
    std::cout << "Only calculating temporal parts (non-relativistic)\n";
  }

  const auto low_q = input.get("low_q", false);
  if (low_q)
    std::cout << "\nExplicitely using low-q form of operators.\n"
                 "Valid only for q << 1/a0 ~ 1e-3 MeV\n\n";

  // Multipolarity:
  const auto [Kmin, Kmax] = input.get("K_minmax", std::array{0, 6});
  fmt::print("\nIncluding K = {} - {}\n", Kmin, Kmax);

  // Optional:
  const auto lc_minmax = input.get<std::array<int, 2>>("lc_minmax");
  if (lc_minmax) {
    fmt::print("Limitting continuum state orbital L to: {} <= L <= {}\n",
               lc_minmax->at(0), lc_minmax->at(1));
  }

  const auto ec_minmax_eV = input.get<std::array<double, 2>>("ec_minmax");
  const double ec_min =
    ec_minmax_eV ? ec_minmax_eV->at(0) / PhysConst::Hartree_eV : 0.0;
  const double ec_max =
    ec_minmax_eV ? ec_minmax_eV->at(1) / PhysConst::Hartree_eV : 1.0 / 0.0;
  if (ec_minmax_eV) {
    fmt::print("Limitting continuum state energy to: {} <= E/eV <= "
               "{} ({:.1e} <= E/au <= {:.1e})\n",
               ec_minmax_eV->at(0), ec_minmax_eV->at(1), ec_min, ec_max);
  }

  // Method for continuum states:
  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);

  // H-like (Zeff) approximation (for testing/comparison):
  // If the Zeff option is set (to anything), the default method is Zeff.
  // Zeff = 'true' or <= 0 (or not set): "real" Zeff = n*sqrt(-2*en);
  // Zeff = positive value: use that constant Zeff for all states.
  const auto zeff_input = input.get<std::string>("Zeff");
  const auto t_method =
    input.get<std::string>("method", zeff_input ? "Zeff" : "HF");

  const auto states_method = Kion::parseStatesMethod(t_method);
  const bool use_Zeff = states_method != Kion::AtomicMethod::HF;
  const bool Zeff_analytic = states_method == Kion::AtomicMethod::ZeffAnalytic;

  // Zeff value: positive number -> constant; 'true'/non-numeric/<=0 -> real:
  double Zeff_constant = 0.0;
  if (zeff_input && !qip::ci_compare(*zeff_input, "true")) {
    // strtod returns 0.0 for non-numeric input (-> real Zeff)
    const double z = std::strtod(zeff_input->c_str(), nullptr);
    Zeff_constant = z > 0.0 ? z : 0.0;
  }

  if (zeff_input && !use_Zeff) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Zeff option has no effect for method=HF; ignoring\n");
  }
  if (Zeff_analytic && !DiracContinuum::available) {
    fmt2::styled_print(fg(fmt::color::red), "\nFail: ");
    fmt::print("method=ZeffAnalytic requires ampsci be compiled with FLINT; "
               "use method=Zeff (solves with DiracODE)\n");
    return;
  }

  std::cout << "\n";
  if (use_Zeff) {
    std::cout << "Using H-like (Zeff) bound and continuum states";
    if (Zeff_constant != 0.0) {
      fmt::print(": constant Zeff = {:g}\n", Zeff_constant);
    } else {
      std::cout << ": Zeff = n*sqrt(-2*en) for each state\n";
    }
    std::cout << (Zeff_analytic ?
                    "Using exact analytic H-like functions\n" :
                    "Solving H-like states numerically (DiracODE)\n");
    std::cout << "(force_rescale/hole_particle have no effect for Zeff "
                 "states)\n";
  }
  if (force_rescale && !use_Zeff) {
    std::cout << "Force rescale: Enforcing V(r) ~ -Z_ion/r at large r\n";
  }
  if (hole_particle && !use_Zeff) {
    std::cout << "Subtracting HF self-interaction (account for hole-particle "
                 "interaction)\n";
  }
  if (force_orthog) {
    std::cout << "Explicitely enforcing orthogonality between bound and "
                 "continuum states\n";
  }
  if (!use_Zeff && !force_rescale && !hole_particle &&
      wf.vHF()->method() == HF::Method::HartreeFock) {
    fmt::print("\n"
               "Warning: Long-range behaviour of V(r) may be incorrect.\n"
               "Suggest to either force rescaling, or include hole-particle "
               "interaction.\n");
  }
  if (!use_Zeff && force_rescale && hole_particle) {
    // meaningless to include
    fmt2::styled_print(fg(fmt::color::red), "\nFail: ");
    fmt::print("Do not force rescaling of V(r) if also subtracting hole "
               "particle (self-interaction) potential. Meaningless.\n");
    return;
  }
  std::cout << "\n";

  // Spherical Bessel lookup table.
  // Each operator stores a pointer to this table.
  std::cout << "Filling jL spherical Bessel table.." << std::flush;
  const SphericalBessel::JL_table jK_tab(Kmax + 1, qgrid, wf.grid().r());
  std::cout << "..done\n" << std::flush;

  // Matrices for each K(E,q) form factor
  // These will be empty (size==0), unless we calculate for relevant operator
  // nb: Order matters, so that output file columns match expected.
  // Always:
  // Vector (V_T, V_E, V_M, V_L, X),
  // Axial (A_T, A_E, A_M, A_L, Y),
  // Scalar (S),
  // Pseudoscalar (P)
  std::array<LinAlg::Matrix<double>, 13> K_factors;

  // Note: important that these are in the same order at the form-factor arrays
  // Always: Vector, Axial, Scalar, Pseudoscalar
  const auto titles =
    std::vector<std::string>{"V_T", "V_E", "V_M", "V_L", "X", "A_T", "A_E",
                             "A_M", "A_L", "Y",   "Z",   "S", "P"};
  const auto descriptions =
    std::vector<std::string>{"Vector (temporal)",
                             "Vector (electric)",
                             "Vector (magnetic)",
                             "Vector (longitudinal)",
                             "Vector (spatial-temporal cross term)",
                             "Axial (temporal)",
                             "Axial (electric)",
                             "Axial (magnetic)",
                             "Axial (longitudinal)",
                             "Axial (spatial-temporal cross term)",
                             "Vector-Axial interference",
                             "Scalar",
                             "Pseudoscalar"};

  //----------------------------------------------------------------------------
  // Output format: Filename and descriptions:

  // Output file name: depends on approximation, options
  // (force_rescale/hole_particle have no effect for Zeff states)
  const auto base_method =
    use_Zeff ? (Zeff_constant != 0.0 ? fmt::format("Zeff{:g}", Zeff_constant) :
                                       std::string{"Zeff"}) +
                 (Zeff_analytic ? "an" : "") :
               HF::parseMethod_short(wf.vHF()->method());
  const auto method =
    base_method + (force_rescale && !use_Zeff ? "_rescale" : "") +
    (hole_particle && !use_Zeff ? "_hp" : "") + (force_orthog ? "_orth" : "");

  const auto units = Kion::Units::Particle;
  const int num_digits = 6;

  // optional extra label
  const auto label = input.get("label", std::string{""});

  std::string ofname_prefix =
    wf.identity() + "_"                                       //
    + method + "_"                                            //
    + std::to_string(Kmin) + "-" + std::to_string(Kmax) + "_" //
    + (lc_minmax ? std::to_string(lc_minmax->at(0)) + "-" +
                     std::to_string(lc_minmax->at(1)) + "_" :
                   "")               //
    + (ec_minmax_eV ? "eclim_" : "") //
    + (low_q ? "lowq_" : "");        //

  if (vectorQ)
    ofname_prefix += "V";
  if (axialQ)
    ofname_prefix += "A";
  if (scalarQ)
    ofname_prefix += "S";
  if (pseudoscalarQ)
    ofname_prefix += "P";

  if (!label.empty()) {
    ofname_prefix += "_" + label;
  }

  //-------------------------------------------------------------------------
  std::cout << "\nCalculating ionisation factors for each bound electron:\n";

  int count = 1;
  for (const auto &Fa : wf.core()) {

    // Find min/max allowed l for continuum states
    // l_min = j_min - 1/2
    //       = j - K_max - 1/2
    //       = (2j - 2*Kmax - 1) / 2
    // nb: 2j is always odd! And allow for either parity.
    const auto lc_min_tmp = std::max((Fa.twoj() - 2 * Kmax - 1) / 2, 0);
    const auto lc_max_tmp = (Fa.twoj() + 2 * Kmax + 1) / 2;
    const auto lc_min =
      lc_minmax ? std::max(lc_min_tmp, lc_minmax->at(0)) : lc_min_tmp;
    const auto lc_max =
      lc_minmax ? std::min(lc_max_tmp, lc_minmax->at(1)) : lc_max_tmp;

    fmt::print("[{:2}/{:2}]: {:4s} -> {} - {}  (L = {} -> {} - {})\n", count++,
               wf.core().size(), Fa.shortSymbol(), AtomData::l_symbol(lc_min),
               AtomData::l_symbol(lc_max), Fa.l(), lc_min, lc_max);
    if (use_Zeff) {
      fmt::print("         Zeff = {:.4f}\n",
                 Zeff_constant != 0.0 ? Zeff_constant :
                                        Kion::Zeff_real(Fa.en(), Fa.n()));
    }
    std::cout << std::flush;

    // Form factors for specific bound state, Fa:
    const auto K_factors_nk = Kion::calculate_formFactors_nk(
      wf.vHF(), Fa, lc_min, lc_max, ec_min, ec_max, force_rescale,
      hole_particle, force_orthog, Egrid, qgrid, diagonal_Eq, low_q, jK_tab,
      Kmin, Kmax, vectorQ, axialQ, scalarQ, pseudoscalarQ, spatialQ,
      states_method, Zeff_constant);

    assert(K_factors_nk.size() == K_factors.size());

    // Optionally: write all K factors for each bound state to disk
    if (each_state) {
      Kion::write_to_file_xyz_13(
        Fa.shortSymbol() + "." + ofname_prefix + ".txt", Egrid, qgrid, titles,
        descriptions, K_factors_nk, units, num_digits, diagonal_Eq);
    }

    for (std::size_t i = 0; i < K_factors.size(); ++i) {
      if (K_factors_nk[i].empty())
        continue;
      if (K_factors[i].empty()) {
        // resize if first non-zero nk
        K_factors[i].resize(E_steps, q_steps);
      }
      K_factors[i] += K_factors_nk[i];
    }
  }
  std::cout << "done\n\n";

  std::cout << "Calculated: ";
  for (std::size_t i = 0; i < titles.size(); ++i) {
    if (!K_factors[i].empty()) {
      std::cout << titles[i] << ", ";
    }
  }
  std::cout << "\n\n";

  //----------------------------------------

  // Write total (summed) form factors to disk
  Kion::write_to_file_xyz_13(ofname_prefix + ".txt", Egrid, qgrid, titles,
                             descriptions, K_factors, units, num_digits,
                             diagonal_Eq);
}

} // namespace Module
