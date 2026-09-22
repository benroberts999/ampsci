#include "DiracODE/ContinuumState.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHFcomplex.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Kionisation/Kion_functions.hpp"
#include "Kionisation/Kion_ridge.hpp"
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
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include "qip/Widgets.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
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
  "photoRPA", "Photo-ionisation cross-section with RPA (outgoing-wave TDHF)",
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
    {"ec_max", "Cut-off (in au) for continuum energy. [1e99]"},
    {"K_minmax", "List (2). Minimum, maximum K [1, 1]"},
    {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
    {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                      "hole-particle interaction) [true]"},
    {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
    {"label", "Optional extra label appended to output file name"},
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

  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);
  const auto label = input.get("label", std::string{""});

  // Output file name: identity, method, continuum options, optional label
  const auto oname =
    wf.identity() + "_photo_" +
    (wf.vHF()->method() == HF::Method::HartreeFock ? "hf" : "local") +
    (force_rescale ? "_rescale" : "") + (hole_particle ? "_hp" : "") +
    (force_orthog ? "_orth" : "") + (label.empty() ? "" : "_" + label) + ".txt";

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
           << "# sigma_Mk1      : Mk at K=1\n"
           << "#\n"
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
    {"", "Photoionisation cross-section for one or more operators, without "
         "or with RPA (core polarisation) from the outgoing-wave TDHF "
         "(TDHFcntm)."},
    {"operator", "List. Operators: any of E1, E1v, M1, E2, VEk, VEk_Len, VMk. "
                 "Each is written to own column in output file. [E1]"},
    {"method", "HF (Hartree-Fock cross-section only) or RPA (with "
               "RPA/core-polarisation corrections; requires the HF method). "
               "As formFactors, without the Zeff options [RPA]"},
    {"K_minmax", "List (2). Minimum, maximum multipolarity K [1, 1]"},
    {"E_range", "List (2). Minimum, maximum photon energy, in eV [10, 1000]"},
    {"E_steps", "Number of photon energies (logarithmic grid) [64]"},
    {"E_threshold", "Number of extra energies to add in the 15% range on "
                    "either side of each ionisation threshold. If <2, none "
                    "are added [0]"},
    {"E_extra", "List (comma separated) of extra energies (in eV) to add 10 "
                "points around (20% either side); for specific regions that "
                "need more resolution"},
    {"max_its", "Maximum RPA iterations per omega [60]"},
    {"eps", "RPA convergence target [1e-10]"},
    {"eps_fail", "RPA solutions whose final eps is above this (or nan) are "
                 "discarded: the RPA shift is interpolated from the "
                 "neighbouring energies, else the no-RPA value is used "
                 "[1e-3]"},
    {"rpa_E_max", "Solve the RPA only for photon energies up to this (eV); "
                  "the no-RPA value is used above it [no limit]"},
    {"rpa_K_max", "Solve the RPA only for multipoles K up to this; the no-RPA "
                  "value is used for the higher K [K_max]"},
    {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                      "hole-particle interaction) [true]"},
    {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
    {"label", "Optional extra label appended to output file name"},
  });
  if (input.has_option("help")) {
    return;
  }

  // Operator names, matched case-insensitively against the DiracOperator
  // names. A list keeps its brackets in the parsed entries: stripped here
  const std::vector<std::string> supported_operators{
    "E1", "E1v", "M1", "E2", "VEk", "VEk_Len", "VMk"};
  std::vector<std::string> operators;
  for (auto name : input.get("operator", std::vector<std::string>{"E1"})) {
    name.erase(std::remove_if(name.begin(), name.end(),
                              [](char c) { return c == '[' || c == ']'; }),
               name.end());
    const auto match =
      std::find_if(supported_operators.begin(), supported_operators.end(),
                   [&name](const std::string &supported) {
                     return qip::ci_compare(supported, name);
                   });
    if (match == supported_operators.end()) {
      fmt2::styled_print(fg(fmt::color::red), "\nFail: ");
      fmt::print("photoRPA: unknown or unsupported operator {}\n", name);
      return;
    }
    operators.push_back(*match);
  }
  const auto n_ops = operators.size();

  const auto [Kmin, Kmax] = input.get("K_minmax", std::array{1, 1});
  auto [Emin_eV, Emax_eV] = input.get("E_range", std::array{10.0, 1000.0});
  auto E_steps = input.get<std::size_t>("E_steps", 64);
  const auto max_its = input.get("max_its", 60);
  const auto eps = input.get("eps", 1.0e-10);
  const auto eps_fail = input.get("eps_fail", 1.0e-3);
  const auto hole_particle = input.get("hole_particle", true);
  const auto force_orthog = input.get("force_orthog", true);
  const auto label = input.get("label", std::string{""});

  // Method: HF (no RPA) or RPA; the Zeff methods of formFactors do not apply
  const auto method =
    Kion::parseStatesMethod(input.get("method", std::string{"RPA"}));
  if (method != Kion::AtomicMethod::HF && method != Kion::AtomicMethod::RPA) {
    fmt2::styled_print(fg(fmt::color::red), "\nFail: ");
    fmt::print("photoRPA: method must be HF or RPA (have {})\n",
               Kion::parseStatesMethod(method));
    return;
  }
  const bool use_rpa = method == Kion::AtomicMethod::RPA;
  if (use_rpa && wf.vHF()->method() != HF::Method::HartreeFock) {
    fmt2::styled_print(fg(fmt::color::red), "\nFail: ");
    fmt::print("method=RPA requires a Hartree-Fock core (have {}); RPA is "
               "meaningless for a local potential\n",
               HF::parseMethod_short(wf.vHF()->method()));
    return;
  }

  // Limits on where the RPA is solved (the no-RPA value is used elsewhere)
  const auto rpa_E_max_eV = input.get<double>("rpa_E_max");
  const auto rpa_E_max =
    rpa_E_max_eV ? *rpa_E_max_eV / PhysConst::Hartree_eV : 1.0e99;
  const auto rpa_K_max = input.get("rpa_K_max", Kmax);

  // Output file name: identity, method, operators, K, continuum options,
  // label
  const auto oname =
    wf.identity() + "_photo_" + Kion::parseStatesMethod(method) + "_" +
    qip::concat(operators, "-") + fmt::format("_{}-{}", Kmin, Kmax) +
    (hole_particle ? "_hp" : "") + (force_orthog ? "_orth" : "") +
    (label.empty() ? "" : "_" + label) + ".txt";

  if (E_steps <= 1) {
    E_steps = 1;
    Emax_eV = Emin_eV;
  }
  const auto Emin_au = Emin_eV / PhysConst::Hartree_eV;
  const auto Emax_au =
    Emax_eV < Emin_eV ? Emin_au : Emax_eV / PhysConst::Hartree_eV;
  auto energies = qip::logarithmic_range(Emin_au, Emax_au, E_steps);

  // Extra energies near each ionisation threshold: just below, and just
  // above (the continuum states cannot be solved well at very small energy,
  // so start e0 above the threshold)
  const auto E_threshold = input.get<std::size_t>("E_threshold", 0);
  if (E_threshold > 1) {
    for (const auto &Fc : wf.core()) {
      const auto below =
        qip::uniform_range(-0.85 * Fc.en(), -0.999 * Fc.en(), E_threshold);
      const auto e0 = 0.01;
      const auto above =
        qip::uniform_range(-Fc.en() + e0, 1.15 * (-Fc.en() + e0), E_threshold);
      energies = qip::merge(energies, below, above);
    }
  }
  // Extra energies around specific values (20% either side)
  const auto E_extra_eV = input.get("E_extra", std::vector<double>{});
  for (const auto Em_eV : E_extra_eV) {
    const auto Em = Em_eV / PhysConst::Hartree_eV;
    energies = qip::merge(energies, qip::uniform_range(0.8 * Em, 1.2 * Em, 10));
  }
  if (E_threshold > 1 || !E_extra_eV.empty()) {
    std::sort(energies.begin(), energies.end());
  }
  const auto n_E = energies.size();

  // Few energies: solve them in order (each warm starts from the previous)
  // with the parallelism inside the RPA solver. Many energies (or no RPA):
  // parallelise over the energies instead, each thread owning its own
  // solver.
  constexpr std::size_t parallel_omega_threshold = 32;
  const bool parallel_omega = !use_rpa || n_E > parallel_omega_threshold;

  // Operator by name; the frequency-dependent ones are updated per omega.
  // k is the multipolarity of the VEk, VEk_Len, VMk family; the others have
  // a fixed rank
  const auto make_operator =
    [&](const std::string &name, int k,
        double omega) -> std::unique_ptr<DiracOperator::TensorOperator> {
    if (name == "E1")
      return std::make_unique<DiracOperator::E1>(wf.grid());
    if (name == "E1v")
      return std::make_unique<DiracOperator::E1v>(wf.alpha(), omega);
    if (name == "M1")
      return std::make_unique<DiracOperator::M1>(wf.grid(), wf.alpha(), omega);
    if (name == "E2")
      return std::make_unique<DiracOperator::E2>(wf.grid());
    if (name == "VEk")
      return std::make_unique<DiracOperator::VEk>(wf.grid(), k, omega);
    if (name == "VEk_Len")
      return std::make_unique<DiracOperator::VEk_Len>(wf.grid(), k, omega);
    // VMk
    return std::make_unique<DiracOperator::VMk>(wf.grid(), k, omega);
  };
  const auto variable_rank = [](const std::string &name) {
    return name == "VEk" || name == "VEk_Len" || name == "VMk";
  };

  // The (operator, K) blocks: one RPA solve per block per energy, with the
  // convergence bookkeeping per block. A fixed-rank operator contributes at
  // its own K only (E1 at K = 1, E2 at K = 2, ...)
  struct Block {
    std::size_t i_op;
    int k;
  };
  std::vector<Block> blocks;
  for (int k = Kmin; k <= Kmax; ++k) {
    for (std::size_t i_op = 0; i_op < n_ops; ++i_op) {
      const auto &name = operators[i_op];
      if (variable_rank(name) ||
          make_operator(name, k, energies.front())->rank() == k) {
        blocks.push_back({i_op, k});
      }
    }
  }

  fmt::print("\nCore ionisation energies:\n");
  fmt::print("{:>4} {:>12} {:>12}\n", "", "au", "eV");
  for (const auto &Fc : wf.core()) {
    fmt::print("{:>4} {:12.6f} {:12.3f}\n", Fc.shortSymbol(), -Fc.en(),
               -Fc.en() * PhysConst::Hartree_eV);
  }
  fmt::print("\nOperators: {}; K = {} to {}\n", qip::concat(operators, ", "),
             Kmin, Kmax);
  fmt::print("Energy   : [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} "
             "steps\n",
             energies.front() * PhysConst::Hartree_eV,
             energies.back() * PhysConst::Hartree_eV, energies.front(),
             energies.back(), n_E);
  if (!use_rpa) {
    std::cout << "Method   : HF (no RPA); parallel over energies\n";
  } else {
    fmt::print("Method   : RPA, eps target {:.1e}, discarded if eps > {:.1e}; "
               "parallel over {}\n",
               eps, eps_fail, parallel_omega ? "energies" : "RPA channels");
  }
  if (use_rpa && (rpa_E_max_eV || rpa_K_max < Kmax)) {
    std::cout << "RPA solved only for:";
    if (rpa_E_max_eV) {
      fmt::print(" E <= {:.4g} eV;", *rpa_E_max_eV);
    }
    if (rpa_K_max < Kmax) {
      fmt::print(" K <= {};", rpa_K_max);
    }
    std::cout << " no-RPA value elsewhere\n";
  }
  if (blocks.empty()) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    std::cout << "no operator contributes in this K range; nothing to do\n";
    return;
  }

  // The operator pair for the RPA solver: t_+ and t_-. E1v depends on the
  // sign of omega, so t_- is E1v at -omega. Every other operator depends on
  // |omega| only, so t_- = t_+^dagger (automatic: nullptr)
  using OperatorPtr = std::unique_ptr<DiracOperator::TensorOperator>;
  const auto make_operator_pair =
    [&](const Block &block) -> std::pair<OperatorPtr, OperatorPtr> {
    const auto &name = operators[block.i_op];
    return {make_operator(name, block.k, energies.front()),
            name == "E1v" ? make_operator(name, block.k, -energies.front()) :
                            nullptr};
  };

  // Angular/polarisation factor: dimensionless absorption form factor Q is
  // f * |<e||t||a>|^2 summed over channels (as in photo{})
  const auto Q_factor = [&](const Block &block, double omega) {
    const auto &name = operators[block.i_op];
    const auto q = PhysConst::alpha * omega;
    if (name == "E1" || name == "E1v")
      return 1.0 / 3.0;
    if (name == "M1")
      return 1.0 / 3.0 * qip::pow(PhysConst::muB_CGS, 2);
    if (name == "E2")
      return 1.0 / 3.0 / 20.0 * q * q;
    return (2.0 * block.k + 1.0) / 2.0 / (q * q);
  };

  // Below the lowest threshold no channel is open: nothing to do
  const auto any_open = [&wf](double omega) {
    return std::any_of(
      wf.core().begin(), wf.core().end(),
      [omega](const auto &Fa) { return omega + Fa.en() > 0.0; });
  };

  // A first-order solve (max_its <= 1) is what was asked for: its eps is
  // the size of the correction, not a convergence measure, and is not tested
  const bool test_convergence = max_its > 1;
  const auto converged = [=](double eps_rpa) {
    return !test_convergence || (!std::isnan(eps_rpa) && eps_rpa < eps_fail);
  };

  // Cross-sections without and with RPA, and the final RPA eps, per block
  // and photon energy [block][energy]. Zero below every threshold; eps zero
  // where the RPA was not solved (outside the limits)
  LinAlg::Matrix<double> sigma_0(blocks.size(), n_E, 0.0);
  LinAlg::Matrix<double> sigma_rpa(blocks.size(), n_E, 0.0);
  LinAlg::Matrix<double> eps_rpa(blocks.size(), n_E, 0.0);

  // Solves the RPA for block ib at energies[i_omega] (unless there is no
  // solver, or outside the RPA limits), and accumulates the cross-sections
  // over the open channels of every core orbital. Returns the final RPA eps
  // (zero if not solved). An unconverged solve is retried once from a
  // cleared state (the warm start may be poor near a resonance); if that
  // also fails, the no-RPA value is used for sigma_rpa, and the solver is
  // cleared so the next energy does not warm start from a broken state.
  const auto solve_energy = [&](std::size_t ib, std::size_t i_omega,
                                DiracOperator::TensorOperator *h,
                                DiracOperator::TensorOperator *h_minus,
                                ExternalField::TDHFcntm *rpa, bool print) {
    const auto &block = blocks[ib];
    const auto omega = energies[i_omega];

    if (h->freqDependantQ()) {
      h->updateFrequency(omega);
      if (h_minus) {
        h_minus->updateFrequency(-omega);
      }
    }

    const bool solve_rpa =
      rpa != nullptr && block.k <= rpa_K_max && omega <= rpa_E_max;
    double eps_final = 0.0;
    bool dressed = false;
    if (solve_rpa) {
      rpa->solve_core(omega, max_its, print);
      if (!converged(rpa->last_eps())) {
        rpa->clear();
        rpa->solve_core(omega, max_its, print);
      }
      eps_final = rpa->last_eps();
      dressed = converged(eps_final);
      if (!dressed) {
        rpa->clear();
      }
    }

    // Conversion factor from dimensionless Q absorption form factor to sigma
    const auto Ksigma = 4.0 * M_PI * M_PI * PhysConst::alpha *
                        PhysConst::aB_cm * PhysConst::aB_cm * omega;
    const auto f_Q = Q_factor(block, omega);

    for (const auto &Fa : wf.core()) {
      const auto ec = omega + Fa.en();
      if (ec < 0.0)
        continue;

      // Continuum states of the ejected electron, in V^(N-1) of hole Fa
      const int lc_max = Fa.l() + h->rank() + 1;
      const int lc_min = std::max(Fa.l() - h->rank() - 1, 0);
      ContinuumOrbitals cntm(wf.vHF());
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa, false, hole_particle,
                            force_orthog);

      for (const auto &Fe : cntm.orbitals) {
        if (h->isZero(Fe, Fa))
          continue;
        const auto t0 = h->reducedME(Fe, Fa);
        const auto t0_sq = t0 * t0;
        const auto t_rpa_sq =
          dressed ? std::norm(t0 + rpa->dV_complex(Fe, Fa)) : t0_sq;
        sigma_0(ib, i_omega) += Ksigma * f_Q * t0_sq;
        sigma_rpa(ib, i_omega) += Ksigma * f_Q * t_rpa_sq;
      }
    }
    return eps_final;
  };

  for (std::size_t ib = 0; ib < blocks.size(); ++ib) {
    const auto &block = blocks[ib];
    fmt::print("{} K={}\n", operators[block.i_op], block.k);

    qip::ProgressBar bar(n_E, parallel_omega);
    // Each thread owns its operators and solver (none without RPA): the
    // frequency-dependent operators are updated in place, and the solver
    // warm starts from the last energy the thread solved. In serial (team
    // of one) the energies are visited in order, so each warm starts from
    // the previous one.
#pragma omp parallel if (parallel_omega)
    {
      // Bug(?): OMP cannot work with structured bindings?
      const auto operator_pair = make_operator_pair(block);
      const auto &h = operator_pair.first;
      const auto &h_minus = operator_pair.second;
      std::unique_ptr<ExternalField::TDHFcntm> rpa;
      if (use_rpa) {
        rpa = std::make_unique<ExternalField::TDHFcntm>(h.get(), wf.vHF(),
                                                        h_minus.get());
        rpa->eps_target() = eps;
      }
#pragma omp for schedule(dynamic)
      for (std::size_t i_omega = 0; i_omega < n_E; ++i_omega) {
        if (any_open(energies[i_omega])) {
          eps_rpa(ib, i_omega) = solve_energy(
            ib, i_omega, h.get(), h_minus.get(), rpa.get(), !parallel_omega);
        }
        bar.update();
      }
    }
  }

  // Failed solves: the relative RPA shift of the block is interpolated from
  // the neighbouring energies (as for the form factors); with no converged
  // neighbour the no-RPA value stays. Not after a first-order solve, whose
  // eps is not a convergence measure
  if (use_rpa && test_convergence) {
    const auto [n_failed, n_interpolated] =
      Kion::count_failed_rpa(eps_rpa, eps_fail);
    Kion::interpolate_failed_rpa(eps_rpa, eps_fail, sigma_0, &sigma_rpa);
    if (n_failed > 0) {
      fmt::print("\nNote: RPA not converged (eps > {:.0e}) at {} of {} "
                 "(energy, operator, K) points: dRPA interpolated in energy "
                 "for {}; no-RPA used for {}\n",
                 eps_fail, n_failed, blocks.size() * n_E, n_interpolated,
                 n_failed - n_interpolated);
    }
  }

  // Per operator: summed over its K blocks
  LinAlg::Matrix<double> sigma_op(n_ops, n_E, 0.0);
  LinAlg::Matrix<double> sigma_rpa_op(n_ops, n_E, 0.0);
  for (std::size_t ib = 0; ib < blocks.size(); ++ib) {
    const auto i_op = blocks[ib].i_op;
    for (std::size_t i_omega = 0; i_omega < n_E; ++i_omega) {
      sigma_op(i_op, i_omega) += sigma_0(ib, i_omega);
      sigma_rpa_op(i_op, i_omega) += sigma_rpa(ib, i_omega);
    }
  }

  std::ofstream out_file(oname);
  out_file << "# Photoionisation cross section (cm^2), K = " << Kmin << " to "
           << Kmax << "\n# Per operator, summed over K: sigma (no RPA)";
  if (use_rpa) {
    out_file << ", sigma_rpa (with RPA; where a solve failed, eps > "
             << eps_fail
             << " or nan, the RPA shift was interpolated from the "
                "neighbouring energies, else sigma_rpa = sigma)";
  }
  out_file << "\n# omega_eV";
  for (const auto &name : operators) {
    out_file << "  sigma_" << name;
    if (use_rpa) {
      out_file << "  sigma_rpa_" << name;
    }
  }
  out_file << "\n";
  for (std::size_t i_omega = 0; i_omega < n_E; ++i_omega) {
    fmt::print(out_file, "{:.6e}", energies[i_omega] * PhysConst::Hartree_eV);
    for (std::size_t i_op = 0; i_op < n_ops; ++i_op) {
      fmt::print(out_file, " {:.6e}", sigma_op(i_op, i_omega));
      if (use_rpa) {
        fmt::print(out_file, " {:.6e}", sigma_rpa_op(i_op, i_omega));
      }
    }
    out_file << "\n";
  }
  fmt::print("\nWritten to {}\n", oname);
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
      "    EN qM K1(EN,qM) ... KN(EN,qM)\n"
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
     {"ridge_correction",
      "true/false. Bethe-ridge correction: completes the multipole sum above "
      "K_max with free (plane-wave) states for the ejected electron. A sum "
      "truncated at K_max misses the response near the quasi-free ridge "
      "E ~ q^2/2m, where multipoles up to K ~ q*r contribute; off the ridge "
      "the correction vanishes. Adds '_ridge' to the output file name. "
      "Skipped (with a warning) if K_max < 2*l_max of the core, or with "
      "diagonal, lc_minmax, or low_q. [false]"},
     {"force_rescale", "Rescale atomic potential V(r) at large r when solving "
                       "continuum orbitals. Should be false for local "
                       "potentials or if hole-particle is included. [false]"},
     {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                       "hole-particle interaction) [true]"},
     {"force_orthog", "Enforce orthogonality of the continuum orbitals [true]"},
     {"method",
      "Method for bound and continuum states: HF (standard), RPA (HF states, "
      "with RPA/core-polarisation corrections to every amplitude from the "
      "outgoing-wave TDHF; requires the HF method; writes the bare HF file "
      "as well as the RPA file), Zeff (H-like, solved numerically with "
      "DiracODE), ZeffAnalytic (H-like, exact analytic functions; requires "
      "FLINT). [HF; Zeff if Zeff option is set]"},
     {"rpa_max_its", "RPA: maximum iterations per solve; 1 gives the "
                     "first-order correction [60]"},
     {"rpa_eps", "RPA: convergence target [1e-10]"},
     {"rpa_eps_fail", "RPA: a solve whose final eps is above this (or nan) is "
                      "discarded, and the no-RPA value used there [1e-3]"},
     {"rpa_E_max", "RPA: solve the RPA only for energy transfers E up to this "
                   "(eV); above it, the RPA factors are the bare (HF) ones. "
                   "[no limit]"},
     {"rpa_q_max", "RPA: solve the RPA only for momentum transfers q up to "
                   "this (eV); above it, the RPA factors are the bare ones. In "
                   "the diagonal case, this limits E (q = E/c). [no limit]"},
     {"rpa_K_max", "RPA: solve the RPA only for multipoles K up to this; the "
                   "higher multipoles contribute their bare factors. [K_max]"},
     {"Zeff", "true/false. Use H-like (Zeff) bound and continuum states, with "
              "Zeff = n*sqrt(-2*en) from each binding energy; the same as "
              "method = Zeff. [false]"}});
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
  // Bug(?) capturing structured bindings in lambdas with clang
  const auto K_minmax = input.get("K_minmax", std::array{0, 6});
  const auto Kmin = K_minmax.at(0);
  const auto Kmax = K_minmax.at(1);
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

  // Bethe-ridge (plane-wave) completion of the multipole sum above Kmax.
  // Skipped at the q where the multipoles above Kmax carry less than
  // ridge_eps of the norm (see Kion::calculate_ridge_correction)
  auto ridge_correction = input.get("ridge_correction", false);
  const double ridge_eps = 1.0e-4;
  if (ridge_correction) {
    int l_max_core = 0;
    for (const auto &Fa : wf.core()) {
      l_max_core = std::max(l_max_core, Fa.l());
    }
    if (diagonal_Eq || lc_minmax || low_q) {
      fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
      fmt::print("ridge_correction skipped: not available with diagonal, "
                 "lc_minmax, or low_q\n");
      ridge_correction = false;
    } else if (Kmax < 2 * l_max_core) {
      fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
      fmt::print("ridge_correction skipped: requires K_max >= 2*l_max(core) "
                 "= {} (have {})\n",
                 2 * l_max_core, Kmax);
      ridge_correction = false;
    } else {
      fmt::print("Ridge correction: completing K > {} with free (plane-wave) "
                 "states, at the q where K <= {} carries less than 1 - {:.0e} "
                 "of the norm of exp(iq.r)psi\n",
                 Kmax, Kmax, ridge_eps);
    }
  }

  // Method for continuum states:
  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);

  // H-like (Zeff) approximation (for testing/comparison), with
  // Zeff = n*sqrt(-2*en) for each state. Zeff = true is method = Zeff
  const auto zeff_input = input.get("Zeff", false);
  const auto t_method =
    input.get<std::string>("method", zeff_input ? "Zeff" : "HF");

  const auto method = Kion::parseStatesMethod(t_method);
  // RPA is a method for the amplitudes: the states themselves are those of HF
  const bool use_rpa = method == Kion::AtomicMethod::RPA;
  const auto states_method = use_rpa ? Kion::AtomicMethod::HF : method;
  const bool use_Zeff = states_method != Kion::AtomicMethod::HF;
  const bool Zeff_analytic = states_method == Kion::AtomicMethod::ZeffAnalytic;

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
  if (use_rpa && wf.vHF()->method() != HF::Method::HartreeFock) {
    fmt2::styled_print(fg(fmt::color::red), "\nFail: ");
    fmt::print("method=RPA requires a Hartree-Fock core (have {}); RPA is "
               "meaningless for a local potential\n",
               HF::parseMethod_short(wf.vHF()->method()));
    return;
  }
  Kion::RPAOptions rpa_options;
  rpa_options.max_its = input.get("rpa_max_its", 60);
  rpa_options.eps = input.get("rpa_eps", 1.0e-10);
  rpa_options.eps_fail = input.get("rpa_eps_fail", 1.0e-3);
  // Limits on where the RPA is solved (bare factors elsewhere); input in eV
  const auto rpa_E_max_eV = input.get<double>("rpa_E_max");
  const auto rpa_q_max_eV = input.get<double>("rpa_q_max");
  const auto rpa_K_max = input.get<int>("rpa_K_max");
  if (rpa_E_max_eV) {
    rpa_options.E_max = *rpa_E_max_eV * UnitConv::Energy_eV_to_au;
  }
  if (rpa_q_max_eV) {
    rpa_options.q_max = *rpa_q_max_eV * UnitConv::Momentum_eV_to_au;
  }
  if (rpa_K_max) {
    rpa_options.K_max = *rpa_K_max;
  }

  std::cout << "\n";
  if (use_Zeff) {
    std::cout << "Using H-like (Zeff) bound and continuum states: "
                 "Zeff = n*sqrt(-2*en) for each state\n";
    std::cout << (Zeff_analytic ?
                    "Using exact analytic H-like functions\n" :
                    "Solving H-like states numerically (DiracODE)\n");
    std::cout << "(force_rescale/hole_particle have no effect for Zeff "
                 "states)\n";
  }
  if (use_rpa) {
    fmt::print("Including RPA (core polarisation) via TDHF:\n"
               "  max_its = {}, eps = {:.1e}, discarded if eps > {:.1e}\n",
               rpa_options.max_its, rpa_options.eps, rpa_options.eps_fail);
    if (rpa_E_max_eV || rpa_q_max_eV || rpa_K_max) {
      std::cout << "  RPA solved only for:";
      if (rpa_E_max_eV) {
        fmt::print(" E <= {:.4g} eV;", *rpa_E_max_eV);
      }
      if (rpa_q_max_eV) {
        fmt::print(" q <= {:.4g} eV;", *rpa_q_max_eV);
      }
      if (rpa_K_max) {
        fmt::print(" K <= {};", *rpa_K_max);
      }
      std::cout << " bare factors elsewhere\n";
    }
    if (!hole_particle) {
      fmt2::styled_print(fg(fmt::color::orange), "Warning: ");
      fmt::print("The RPA amplitudes assume continuum states of the residual "
                 "ion (hole_particle=true).\n");
    }
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

  // Spherical Bessel lookup table, on the q grid (in the diagonal case, on
  // q = E/c at each energy). Each operator stores a pointer to this table.
  std::cout << "Filling jL spherical Bessel table.." << std::flush;
  const auto q_table = diagonal_Eq ? Egrid * PhysConst::alpha : qgrid;
  const SphericalBessel::JL_table jK_tab(Kmax + 1, q_table, wf.grid().r());
  std::cout << "..done\n" << std::flush;

  // Titles/descriptions of each K(E,q) form factor, in the order of
  // Kion::FormFactorSet (which is the order of the output file columns):
  // Vector (V_T, V_E, V_M, V_L, X),
  // Axial (A_T, A_E, A_M, A_L, Y),
  // Vector-Axial interference (Z),
  // Scalar (S),
  // Pseudoscalar (P)
  // Always (order): Vector, Axial, Scalar, Pseudoscalar
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
  const auto base_method = use_Zeff ?
                             std::string{"Zeff"} + (Zeff_analytic ? "an" : "") :
                             HF::parseMethod_short(wf.vHF()->method());
  const auto method_suffix =
    std::string{} + (force_rescale && !use_Zeff ? "_rescale" : "") +
    (hole_particle && !use_Zeff ? "_hp" : "") + (force_orthog ? "_orth" : "") +
    (ridge_correction ? "_ridge" : "");

  const auto units = Kion::Units::Particle;
  const int num_digits = 6;

  // optional extra label
  const auto label = input.get("label", std::string{""});

  // e.g., Xe0_HF_hp_orth_0-6_VA; the RPA factors go to Xe0_RPA_hp_orth_0-6_VA
  const auto output_prefix = [&](const std::string &method_text) {
    std::string prefix =
      wf.identity() + "_"                                       //
      + method_text + method_suffix + "_"                       //
      + std::to_string(Kmin) + "-" + std::to_string(Kmax) + "_" //
      + (lc_minmax ? std::to_string(lc_minmax->at(0)) + "-" +
                       std::to_string(lc_minmax->at(1)) + "_" :
                     "")               //
      + (ec_minmax_eV ? "eclim_" : "") //
      + (low_q ? "lowq_" : "")         //
      + (temporal_only ? "T_" : "");   //
    if (vectorQ)
      prefix += "V";
    if (axialQ)
      prefix += "A";
    if (scalarQ)
      prefix += "S";
    if (pseudoscalarQ)
      prefix += "P";
    if (!label.empty()) {
      prefix += "_" + label;
    }
    return prefix;
  };
  const auto ofname_prefix = output_prefix(base_method);
  const auto ofname_prefix_rpa = output_prefix("RPA");

  //-------------------------------------------------------------------------
  if (use_Zeff) {
    std::cout << "\nZeff for each bound electron:\n";
    for (const auto &Fa : wf.core()) {
      fmt::print("{:4s}: Zeff = {:.4f}\n", Fa.shortSymbol(),
                 Kion::Zeff_nonrel(Fa.en(), Fa.n()));
    }
    std::cout << std::flush;
  }

  // Bound states used in the matrix elements: the HF orbitals, or (Zeff
  // methods) their H-like versions
  const auto bound_states = Kion::model_bound_states(*wf.vHF(), states_method);

  // Writes the per-orbital files (if requested), then the total
  const auto write_factors = [&](const std::string &prefix,
                                 const std::vector<Kion::FormFactorSet> &K_nk) {
    if (each_state) {
      for (std::size_t ia = 0; ia < wf.core().size(); ++ia) {
        Kion::write_to_file_xyz_13(
          wf.core()[ia].shortSymbol() + "." + prefix + ".txt", Egrid, qgrid,
          titles, descriptions, K_nk[ia], units, num_digits, diagonal_Eq);
      }
    }
    // Total over the orbitals (factors not calculated stay empty)
    Kion::FormFactorSet total;
    for (const auto &K_factors : K_nk) {
      for (std::size_t i = 0; i < total.size(); ++i) {
        if (K_factors[i].empty())
          continue;
        if (total[i].empty()) {
          total[i] = K_factors[i];
        } else {
          total[i] += K_factors[i];
        }
      }
    }
    std::cout << "Calculated: ";
    for (std::size_t i = 0; i < titles.size(); ++i) {
      if (!total[i].empty()) {
        std::cout << titles[i] << ", ";
      }
    }
    std::cout << "\n";
    Kion::write_to_file_xyz_13(prefix + ".txt", Egrid, qgrid, titles,
                               descriptions, total, units, num_digits,
                               diagonal_Eq);
    std::cout << "\n";
  };

  // Bethe-ridge correction: plane waves and the bound states only, so the
  // same for every amplitude method (added to both the bare and RPA factors)
  const auto ridge = [&]() {
    return Kion::calculate_ridge_correction(
      wf.vHF(), bound_states, ec_min, ec_max, Egrid, qgrid, jK_tab, Kmax,
      vectorQ, axialQ, scalarQ, pseudoscalarQ, spatialQ, ridge_eps);
  };

  if (!use_rpa) {
    auto K_nk = Kion::calculate_formFactors(
      wf.vHF(), bound_states, lc_minmax, ec_min, ec_max, force_rescale,
      hole_particle, force_orthog, Egrid, qgrid, diagonal_Eq, low_q, jK_tab,
      Kmin, Kmax, vectorQ, axialQ, scalarQ, pseudoscalarQ, spatialQ,
      states_method);
    std::cout << "done\n\n";
    if (ridge_correction) {
      Kion::add_formFactors(&K_nk, ridge());
      std::cout << "done\n\n";
    }
    write_factors(ofname_prefix, K_nk);
    return;
  }

  // When we include RPA: Write both (with + without) to disk seperately,
  // Since we have to calculate both anyway

  auto result = Kion::calculate_formFactors_RPA(
    wf.vHF(), lc_minmax, ec_min, ec_max, force_rescale, hole_particle,
    force_orthog, Egrid, qgrid, diagonal_Eq, low_q, jK_tab, Kmin, Kmax, vectorQ,
    axialQ, scalarQ, pseudoscalarQ, spatialQ, rpa_options);
  std::cout << "done\n\n";

  // Try to smoothly interpolate
  if (rpa_options.max_its > 1) {
    // Isolated failed solves leave a step in an otherwise smooth factor
    const auto [n_failed, n_interpolated] =
      Kion::interpolate_failed_rpa(&result, rpa_options.eps_fail);
    if (n_failed > 0) {
      fmt::print(
        "\nNote: RPA not converged (eps > {:.0e}) at {} of {} "
        "(E,q) points: dRPA interpolated in q for {}; no-RPA used for {}\n\n",
        rpa_options.eps_fail, n_failed, E_steps * q_steps, n_interpolated,
        n_failed - n_interpolated);
    }
  }

  if (ridge_correction) {
    const auto dK = ridge();
    std::cout << "done\n\n";
    Kion::add_formFactors(&result.bare, dK);
    Kion::add_formFactors(&result.rpa, dK);
  }
  std::cout << "Without RPA:\n";
  write_factors(ofname_prefix, result.bare);
  std::cout << "With RPA:\n";
  write_factors(ofname_prefix_rpa, result.rpa);
}

} // namespace Module
