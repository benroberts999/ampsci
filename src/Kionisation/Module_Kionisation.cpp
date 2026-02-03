#include "Kionisation/Module_Kionisation.hpp"
#include "DiracOperator/include.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Kionisation/Kion_functions.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/Maths.hpp"
#include "qip/Methods.hpp"
#include <cassert>
#include <iostream>
#include <memory>
#include "Maths/SphericalBessel.hpp"
//
#include "ExternalField/DiagramRPA.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "DiracOperator/Operators/EM_multipole_nr.hpp"

static const std::string Kionisation_description_text{R"(
This module calculates atomic ionisation factors.
Some further detail is given here.

Method:
  - 'hf' is the standard choice, calculated continuum wavefunctions in 
    Hartree-Fock potential of bound states. "hole_particle" option, which is 
    true by default, corrects the Hartree-Fock potential to account for missing 
    ionised electron.
  - RPA means include core-polarisation (many-body) corrections.
    'RPA0' will include only lowest-order RPA corrections, and is quite fast
    'RPA' will include all-orders rpa. This is slow, since the RPA equations
    need to be iterated for each L and q seperately. Therefore, it is advised 
    only to use full RPA for a small subset of E/q grids.
  - Other methods (Zeff etc.) are mainly used for tests, and to compare with 
    other less accurate codes. These are not accurate methods to use.

Output format:
  - xyz:      For easy 2D interpolation. Each row is in form: 'E q K(E,q)'
  - gnuplot:  For easy plotting. Each column is new E
  - matrix:   Outputs entire matrix in table form, with E and q grids 
              printed prior. This is form expected by 'dmex' program

    )"};

namespace Module {

//==============================================================================

void Kionisation(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("Kionisation");

  input.check(
      {{"E_range",
        "List (2). Minimum, maximum energy transfer (dE), in keV [0.1,0.1]"},
       {"E_steps", "Numer of steps along dE grid (logarithmic grid) [1]"},
       {"q_range",
        "List (2). Minimum, maximum momentum transfer (q), in MeV [0.01,0.01]"},
       {"q_steps", "Number of steps along q grid (logarithmic grid) [1]"},
       {"max_L", "Maximum multipolarity used in exp(iqr) expansion [6]"},
       {"ec_cut", "Cut-off (in au) for continuum energy. [1000.0]"},
       {"label", "optional extra label for output files"},
       {"method", "'hf' (relativistic Hartree-Fock), "
                  "'rpa0' (lowst-order RPA), "
                  "'rpa' (all-orders RPA), "
                  "'zeff' (Z_eff for continuum), "
                  "'approx' (step function). [hf]"},
       {"subtract_1", "Replace e^(iqr) -> e^(iqr)-1 [false]"},
       {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
       {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                         "hole-particle interaction) [true]"},
       {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
       {"coupling", "Vector, Scalar (g0), Pseudovector (g5), Pseudoscalar "
                    "(g0g5) [Vector]"},
       {"output_format", "List: Format for output. List any of: gnuplot, xyz, "
                         "matrix (comma-separated). [gnuplot]"},
       {"each_state", "bool. If true, will output K(E,q) for each "
                      "(accessible) core-state [false]"},
       {"units",
        "Units for 'gnuplot' output: Particle (keV/MeV) or Atomic (E_H,1/a0). "
        "Only affects _gnu output format, all _mat and _xyz are "
        "always in atomic units. [Particle]"}});
  if (input.has_option("help")) {
    std::cout << Kionisation_description_text;
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

  const auto ec_cut = input.get("ec_cut", 1000.0);

  fmt::print(
      "Momentum: [{:.3f}, {:.3f}] MeV = [{:.1f}, {:.1f}] au, in {} steps\n",
      qmin_MeV, qmax_MeV, qmin_au, qmax_au, q_steps);

  // Set up the E and q grids
  const Grid Egrid({E_steps, Emin_au, Emax_au, 0, GridType::logarithmic});
  const Grid qgrid({q_steps, qmin_au, qmax_au, 0, GridType::logarithmic});

  // Check to see if grid is reasonable for maximum energy:
  Kion::check_radial_grid(std::min(ec_cut, Emax_au), qmax_au, wf.grid());

  //----------------------------------------------------------------------------
  // Read in and parse options:

  fmt::print("\nMax L   : {}  (multipolarity in e^iqr expansion)\n", max_L);

  if (!label.empty()) {
    fmt::print("Label   : {}\n", label);
  }

  // Method for cntm states:
  using qip::ci_compare;    // case-insensitive string comparison
  using qip::ci_wc_compare; // case-insensitive string comparison with *
  const auto tmethod = input.get<std::string>("method", "hf");
  const auto method = ci_wc_compare(tmethod, "h*")      ? Kion::Method::HF :
                      ci_compare(tmethod, "rpa0")       ? Kion::Method::RPA0 :
                      ci_compare(tmethod, "rpa")        ? Kion::Method::RPA :
                      ci_compare(tmethod, "zeff")       ? Kion::Method::Zeff :
                      ci_wc_compare(tmethod, "approx*") ? Kion::Method::Approx :
                                                          Kion::Method::Error;
  const auto use_Zeff_cont = method == Kion::Method::Zeff;
  if (method == Kion::Method::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 104: ");
    fmt::print("Method option: `{}' unknown. Options are: standard, approx, "
               "zeff\n",
               tmethod);
    return;
  }
  std::cout << "Using   : " << tmethod << " method\n";

  const bool use_rpa0 = method == Kion::Method::RPA0;
  const bool use_rpa_ao = method == Kion::Method::RPA;
  if (use_rpa0 || use_rpa_ao) {
    std::cout << "Using RPA " << (use_rpa0 ? "(lowest-order)" : "(all-orders)")
              << ", with basis: " << DiracSpinor::state_config(wf.basis())
              << "\n";
  }
  if ((use_rpa0 || use_rpa_ao) && wf.basis().empty()) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("RPA method required a basis to work. See `ampsci -a "
               "Basis`\nRPA will evaluate to zero without a basis\n\n");
  }
  // sanity check for RPA - meaningless unless we used HF method
  if ((use_rpa0 || use_rpa_ao) && wf.vHF() &&
      wf.vHF()->method() != HF::Method::HartreeFock) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 113: ");
    fmt::print("It's only meaningful to include RPA if we start with "
               "HartreeFock method (for core states).\n");
    return;
  }

  //----------------------------------------------------------------------------

  // Other methods options
  const auto subtract_1 = input.get("subtract_1", false);
  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale =
      use_Zeff_cont ? false : input.get("force_rescale", false);
  const auto hole_particle =
      use_Zeff_cont ? false : input.get("hole_particle", true);

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
  if (use_Zeff_cont) {
    std::cout
        << "Using Z_eff model for continuum; hole_particle and force_rescale "
           "have no effect. Caution: this should only be used for tests\n";
  }

  // Perform checks, print possible warnings
  if (force_rescale && !(force_orthog || subtract_1)) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Force rescale will ruin orthogonality; suggest use "
               "force_orthog or subtract_1\n");
  }
  if (use_Zeff_cont && !(force_orthog || subtract_1)) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("using Zeff for continuum will ruin orthogonality; suggest use "
               "force_orthog or subtract_1\n");
  }
  if (!force_rescale && !hole_particle && !use_Zeff_cont && wf.Zion() == 0) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print(
        "Long-range behaviour of V(r) may be incorrect. Suggest to either "
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

  // DM-electron couplings
  const auto tcoupling = input.get<std::string>("coupling", "vector");
  const auto coupling =
      ci_compare(tcoupling, "vector")       ? Kion::Coupling::Vector :
      ci_compare(tcoupling, "scalar")       ? Kion::Coupling::Scalar :
      ci_compare(tcoupling, "pseudovector") ? Kion::Coupling::PseudoVector :
      ci_compare(tcoupling, "pseudoscalar") ? Kion::Coupling::PseudoScalar :
                                              Kion::Coupling::Error;
  if (coupling == Kion::Coupling::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 212: ");
    fmt::print("Coupling option: `{}' unknown. Options are: vector, scalar, "
               "pseudovector, pseudoscalar\n",
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
  } else if (coupling == Kion::Coupling::PseudoVector) {
    jl = std::make_unique<DiracOperator::ig5jL>(wf.grid(), qgrid,
                                                std::size_t(max_L));
  } else if (coupling == Kion::Coupling::PseudoScalar) {
    jl = std::make_unique<DiracOperator::ig0g5jL>(wf.grid(), qgrid,
                                                  std::size_t(max_L));
  }
  assert(jl != nullptr && "Error in coupling type");
  std::cout << "Operator: " << jl->name() << "\n";

  //----------------------------------------------------------------------------

  // Create output file-name template

  const std::string rpa_text =
      (use_rpa0 ? "rpa0_" : "rpa_") + DiracSpinor::state_config(wf.basis());

  const auto hf_text = wf.vHF() == nullptr                           ? "??" :
                       wf.vHF()->method() == HF::Method::HartreeFock ? "hf" :
                       wf.vHF()->method() == HF::Method::Hartree     ? "ha" :
                       wf.vHF()->method() == HF::Method::KohnSham    ? "ks" :
                       wf.vHF()->method() == HF::Method::ApproxHF    ? "ahf" :
                       wf.vHF()->method() == HF::Method::Local       ? "loc" :
                                                                       "??";

  const std::string method_text = method == Kion::Method::HF     ? hf_text :
                                  method == Kion::Method::RPA0   ? rpa_text :
                                  method == Kion::Method::RPA    ? rpa_text :
                                  method == Kion::Method::Approx ? "aprx" :
                                  method == Kion::Method::Zeff   ? "zeff" :
                                                                   "???";

  const std::string coupling_text =
      coupling == Kion::Coupling::Vector       ? "v" :
      coupling == Kion::Coupling::Scalar       ? "s" :
      coupling == Kion::Coupling::PseudoVector ? "pv" :
      coupling == Kion::Coupling::PseudoScalar ? "ps" :
                                                 "???";

  // doesn't include suffix
  std::string oname =
      "K_" + wf.identity() + wf.ion_symbol(0) + "_" + method_text + "_";
  if (wf.Zion() != 0)
    oname += fmt::format("{}+_", wf.Zion());
  oname += coupling_text + "_";
  oname += fmt::format("{}_", max_L);
  if (force_rescale)
    oname += "rescale_";
  if (hole_particle)
    oname += "hp_";
  if (force_orthog)
    oname += "orth_";
  if (subtract_1)
    oname += "sub1_";
  if (label != "")
    oname += label + "_";
  if (oname.back() == '_')
    oname.pop_back();

  //----------------------------------------------------------------------------

  // Output format:

  const bool write_each_state = input.get("each_state", false);

  const auto toutput =
      input.get<std::vector<std::string>>("output_format", {"gnuplot"});
  // Use vector, since allow outputting in multiple formats
  std::vector<Kion::OutputFormat> output_formats;
  for (auto &each : toutput) {
    const auto output =
        ci_wc_compare(each, "gnu*") ? Kion::OutputFormat::gnuplot :
        ci_wc_compare(each, "xyz")  ? Kion::OutputFormat::xyz :
        ci_wc_compare(each, "mat*") ? Kion::OutputFormat::matrix :
                                      Kion::OutputFormat::Error;
    if (output == Kion::OutputFormat::Error) {
      fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
      fmt::print(
          "Output Format option: `{}' unknown. Options are: gnuplot, xyz, "
          "matrix\n",
          each);
    } else {
      // only add if not already in list
      if (std::find(output_formats.begin(), output_formats.end(), output) ==
          output_formats.end())
        output_formats.push_back(output);
    }
  }
  if (output_formats.empty()) {
    fmt2::styled_print(fg(fmt::color::red), "\nFail 300: ");
    fmt::print("No output formats? No output will be generated!\n");
    return;
  }

  // Units (only for gnuplot-style):
  const auto tunits = input.get<std::string>("units", "Particle");
  auto units = ci_compare(tunits, "particle") ? Kion::Units::Particle :
               ci_compare(tunits, "atomic")   ? Kion::Units::Atomic :
                                                Kion::Units::Error;
  if (units == Kion::Units::Error) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Output units option: `{}' unknown. Options are: Particle, "
               "Atomic. Defaulting to atomic\n",
               tunits);
    units = Kion::Units::Atomic;
  }
  std::cout << "_gnu output file (if requested) will use: " << tunits
            << " units\n";
  std::cout << "(_mat and _xyz output files always use atomic units)\n";

  //----------------------------------------------------------------------------

  std::cout << "\nCalculating K(E,q) - ionisation factor\n" << std::flush;
  const int num_output_digits = 5;

  // Kion stored K(dE, q)
  LinAlg::Matrix<double> Kion(Egrid.num_points(), qgrid.num_points());

  if (method == Kion::Method::HF || method == Kion::Method::RPA0) {
  }

  if (method == Kion::Method::RPA) {

    /// XXX print each nk
    const auto Knks = Kion::calculateK_nk_rpa(
        wf.vHF(), wf.core(), max_L, Egrid, jl.get(), force_rescale,
        hole_particle, force_orthog, wf.basis(), wf.identity());

    assert(Knks.size() == wf.core().size());

    for (std::size_t ic = 0; ic < wf.core().size(); ++ic) {
      const auto &Fnk = wf.core().at(ic);
      const auto Fc_accessible = std::abs(Fnk.en()) <= Emax_au;
      if (write_each_state && Fc_accessible) {
        const auto oname_nk = oname + "_" + Fnk.shortSymbol();
        std::cout << "Written to file: " << oname_nk << "\n";
        Kion::write_to_file(output_formats, Knks.at(ic), Egrid, qgrid, oname_nk,
                            num_output_digits, units);
      }
      Kion += Knks.at(ic);
    }

  } else if (method == Kion::Method::Approx) {

    std::cout << "Using   : approx (step-function) method\n";
    const auto K_approx = Kion::calculateK_nk_approx(
        wf.vHF(), wf.core(), max_L, jl.get(), force_rescale, hole_particle,
        force_orthog, use_Zeff_cont, use_rpa0, wf.basis());
    Kion::write_approxTable_to_file(K_approx, wf.core(), qgrid, oname,
                                    num_output_digits, units);
    // convert to "standard" form, for easy comparison
    Kion = Kion::convert_K_nk_approx_to_std(K_approx, Egrid, wf.core());

  } else {
    // all other methods (including standard)

    for (const auto &Fnk : wf.core()) {
      const auto accessible = std::abs(Fnk.en()) < Emax_au;
      if (!accessible)
        continue;
      std::cout << Fnk << ", " << std::flush;
      const auto K_nk = Kion::calculateK_nk(
          wf.vHF(), Fnk, max_L, Egrid, jl.get(), force_rescale, hole_particle,
          force_orthog, use_Zeff_cont, use_rpa0, wf.basis(), ec_cut);
      if (write_each_state) {
        const auto oname_nk = oname + "_" + Fnk.shortSymbol();
        std::cout << "Written to file: " << oname_nk << "\n";
        Kion::write_to_file(output_formats, K_nk, Egrid, qgrid, oname_nk,
                            num_output_digits, units);
      }
      Kion += K_nk;
    }
  }

  std::cout << "\nWritten to file: " << oname << "\n";
  Kion::write_to_file(output_formats, Kion, Egrid, qgrid, oname,
                      num_output_digits, units);

  std::cout << '\n';
}

//==============================================================================
void photo(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("photo");

  input.check({
      {"E_range",
       "List (2). Minimum, maximum energy transfer (dE), in eV [10, 1000]"},
      {"E_steps", "Numer of steps along dE grid (logarithmic grid) [50]"},
      {"E_threshold",
       "Numer of extra E steps to add in -15% range on either side"
       " of each threshold. If <2, will add no new points [0]"},
      {"E_extra", "List (comma separated) extra energies (in eV) to add 10 "
                  "points around. Useful for specific regions we want more "
                  "resolution in."},
      {"oname", "oname"},
      {"ec_cut", "Cut-off (in au) for continuum energy. [inf]"},
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
      const auto extra2 = qip::uniform_range(
          -Fc.en() + e0, 1.15 * (-Fc.en() + e0), E_threshold);
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
  const auto ec_cut = input.get("ec_cut", 1 / 0.0);

  std::cout << "\nSummary of inputs:\n";
  fmt::print(
      "Energy  : [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n\n",
      energies.front() * PhysConst::Hartree_eV,
      energies.back() * PhysConst::Hartree_eV, energies.front(),
      energies.back(), energies.size());

  const auto [Kmin, Kmax] = input.get("K_minmax", std::array{1, 1});
  const auto label = input.get("label", std::string{""});

  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);

  auto oname = input.get("oname", std::string{"out.txt"});
  std::ofstream out_file(oname);

  // "full" dipole operator
  const auto E1 = DiracOperator::E1(wf.grid());
  const auto E2 = DiracOperator::Ek(wf.grid(), 2);

  auto M1nr = DiracOperator::M1nr();

  // Note: This is parallelised very ineficiently.
  // Could be improved quite a bit probably..
  // Also: perhaps faster to use jL lookup table? Maybe not the bottleneck

  int count = 0;
  for (const auto omega : energies) {

    // First, loop through and just find list of what we shall do.
    // THEN parellelise over that!
    std::vector<std::size_t> iclist;
    for (std::size_t i = 0; i < wf.core().size(); ++i) {
      const auto ec = omega + wf.core()[i].en();
      if (ec < 0.0)
        continue;
      iclist.push_back(i);
    }

    fmt::print("{:3} {:9.2f} eV ; {:3} shells accessible\n", count++,
               omega * PhysConst::Hartree_eV, iclist.size());

    const auto Ksigma = 4.0 * M_PI * M_PI * PhysConst::alpha *
                        PhysConst::aB_cm * PhysConst::aB_cm * omega;

    double Q_E1 = 0.0;    // Regular (length) E1 operator
    double Q_M1 = 0.0;    // "Regular" M1
    double Q_M1_nr = 0.0; // Non-relativistic M1 operator
    double Q_Mk1 = 0.0;   // Mk at k=1

    double Q_Ek2 = 0.0; // Ek at K=2
    double Q_E2 = 0.0;  // Regular E2 (length)

    double Q_E = 0.0; // Full electric multipole
    double Q_M = 0.0; // Full magnetic multipole

    double Q_E_len = 0.0; // Full electric multipole (length form)

    for (int k = Kmin; k <= Kmax; ++k) {

      // Electric, magnetic parts
      const auto Ek = DiracOperator::Ek_w(wf.grid(), k, omega);
      const auto Mk = DiracOperator::Mk_w(wf.grid(), k, omega);
      // Magnetic dipole
      const auto M1 = DiracOperator::M1(wf.grid(), PhysConst::alpha, omega);
      // "Length" form - for tests only
      const auto Ek_len = DiracOperator::Ek_w_L(wf.grid(), k, omega);

#pragma omp parallel for reduction(+ : Q_E1, Q_M1, Q_Mk1, Q_M1_nr, Q_Ek2,      \
                                       Q_E2, Q_E, Q_M, Q_E_len)
      for (const auto ic : iclist) {
        const auto &Fa = wf.core()[ic];
        const auto ec = omega + Fa.en();
        if (ec < 0.0)
          continue;
        const auto ec_t = std::min(ec, ec_cut);

        const int l = Fa.l();
        const int lc_max = l + k + 1;
        const int lc_min = std::max(l - k - 1, 0);

        ContinuumOrbitals cntm(wf.vHF());
        cntm.solveContinuumHF(ec_t, lc_min, lc_max, &Fa, force_rescale,
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

    out_file << omega * PhysConst::Hartree_eV / 1e6 << " " << Q_E1 * Ksigma
             << " " << Q_M1 * Ksigma << " " << Q_M1_nr * Ksigma << " "
             << Q_E * Ksigma << " " << Q_E_len * Ksigma << " " << Q_M * Ksigma
             << " " << (Q_E + Q_M) * Ksigma << " " << Q_E2 * Ksigma << " "
             << Q_Ek2 * Ksigma << " " << Q_Mk1 * Ksigma << "\n";
  }
}

//==============================================================================
void formFactors(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("formFactors");

  input.check(
      {{"E_range",
        "List (2). Minimum, maximum energy transfer (dE), in eV [10.0, 1.0e4]"},
       {"E_steps", "Numer of steps along dE grid (logarithmic) [1]"},
       // {"E_threshold",
       //  "Numer of extra E steps to add in -15% range on either side"
       //  " of each threshold. If <2, will add no new points [0]"},
       // {"E_extra", "List (comma separated) extra energies (in eV) to add 10 "
       //             "points around. Useful for specific regions we want more "
       //             "resolution in."},
       {"q_range", "List (2). Minimum, maximum momentum transfer (q), in eV "
                   "(hbar=c=1). For reference, 1/a0 ~ 3730 eV. [1.0e4, 1.0e7]"},
       {"q_steps", "Numer of steps along q grid (logarithmic) [1]"},
       {"oname", "oname"},
       {"operators", "List, comma separated. Any of 'V' (vector), 'A' "
                     "(axial/pseudo-vector), 'S' "
                     "(scalar), 'P' (pseudoscalar), 'C' (cross-terms)"},
       {"K_minmax", "List (2). Minimum, maximum K [0, 5]"},
       {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
       {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                         "hole-particle interaction) [true]"},
       {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
       {"output_format",
        "List: Format for output. List any of: gnuplot, gnuplot_E, xyz, "
        "matrix (comma-separated). gnuplot is for easy plotting as function of "
        "q, gnuplot_E is for easy plotting as function of E [gnuplot]"},
       {"print_all", "Print all spatial parts (E, M, L) seperately? [false]"},
       {"units",
        "Units for 'gnuplot' output: Particle (keV/MeV) or Atomic (E_H,1/a0). "
        "Only affects _gnu output format, all _mat and _xyz are "
        "always in atomic units. [Particle]"}});
  if (input.has_option("help")) {
    return;
  }

  // Read in energy grid. Input is in keV:
  auto [Emin_eV, Emax_eV] = input.get("E_range", std::array{10.0, 1.0e3});
  auto E_steps = input.get<std::size_t>("E_steps", 1);
  if (E_steps <= 1 || Emax_eV <= Emin_eV) {
    E_steps = 1;
    Emax_eV = Emin_eV;
  }
  // Convert to atomic units (from keV):
  const auto Emin_au = Emin_eV / PhysConst::Hartree_eV;
  const auto Emax_au = Emax_eV / PhysConst::Hartree_eV;

  //----------------------------------------------------------------------------

  // Print core information for convenience:
  // This assumes core is in energy order; always true but not guarenteed.
  // Has no impact on result, just what is printed.
  std::cout << "\nCore orbitals:\n";
  std::cout << "        E(au)      E(eV)  N_el\n";
  bool reached_accessible = false;
  if (std::abs(wf.core().front().en()) > Emax_au) {
    std::cout << "------- Inaccessible ---------\n";
  }
  for (const auto &Fnk : wf.core()) {
    if (std::abs(Fnk.en()) <= Emax_au && !reached_accessible) {
      reached_accessible = true;
      std::cout << "-------- Accessible ----------\n";
    }
    fmt::print("{:3s}  {:8.2f}  {:9.2f}     {}\n", Fnk.shortSymbol(), Fnk.en(),
               Fnk.en() * UnitConv::Energy_au_to_eV, Fnk.num_electrons());
  }

  // Momentum-transfer range (in eV)
  auto [qmin_eV, qmax_eV] = input.get("q_range", std::array{1.0, 1.0e4});
  auto q_steps = input.get<std::size_t>("q_steps", 1);
  if (q_steps <= 1 || qmax_eV <= qmin_eV) {
    q_steps = 1;
    qmax_eV = qmin_eV;
  }

  // Convert momentum from keV to atomic units.
  const auto q_min = qmin_eV * UnitConv::Momentum_eV_to_au;
  const auto q_max = qmax_eV * UnitConv::Momentum_eV_to_au;

  // Set up the E and q grids
  const Grid Egrid({E_steps, Emin_au, Emax_au, 0, GridType::logarithmic});
  const Grid qgrid({q_steps, q_min, q_max, 0, GridType::logarithmic});

  // Check to see if grid is reasonable for maximum energy:
  Kion::check_radial_grid(Emax_au, q_max, wf.grid());

  std::cout << "\nEnergy/Momentum exchange grids:\n";
  fmt::print(
      "Energy  : [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n",
      Emin_eV, Emax_eV, Emin_au, Emax_au, E_steps);
  fmt::print(
      "Momentum: [{:.1e}, {:.1e}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n",
      qmin_eV, qmax_eV, q_min, q_max, q_steps);

  // Which operators to calculate
  // Case insensitive, and only check first letter
  const auto operators = input.get("operators", std::vector<std::string>{"V"});
  bool vectorQ{false}, axialQ{false}, scalarQ{false}, pseudoscalarQ{false}, crossQ{false};
  for (auto &w : operators) {
    if (!w.empty()) {
      switch (std::tolower(w[0])) {
      case 'v':
        vectorQ = true;
        break;
      case 'a':
        axialQ = true;
        break;
      case 's':
        scalarQ = true;
        break;
      case 'p':
        pseudoscalarQ = true;
        break;
      case 'c':
        crossQ = true;
      break;
      }
    }
  }
  fmt::print("\nComputing operators: "
             "{}{}{}{}\n",
             vectorQ ? "Vector; " : "", axialQ ? "Axial; " : "",
             scalarQ ? "Scalar; " : "", pseudoscalarQ ? "Pseudoscalar; " : "", crossQ ? "Cross-term; " : "");

  const auto [Kmin, Kmax] = input.get("K_minmax", std::array{0, 5});
  fmt::print("\nIncluding K = {} - {}\n", Kmin, Kmax);

  // Method for continuum states:
  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);
  std::cout << "\n";
  if (force_rescale) {
    std::cout << "Force rescale: Enforcing V(r) ~ -Z_ion/r at large r\n";
  }
  if (hole_particle) {
    std::cout << "Subtracting HF self-interaction (account for hole-particle "
                 "interaction)\n";
  }
  if (force_orthog) {
    std::cout << "Explicitly enforcing orthogonality between bound and "
                 "continuum states\n";
  }
  if (!force_rescale && !hole_particle &&
      wf.vHF()->method() == HF::Method::HartreeFock) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print(
        "Long-range behaviour of V(r) may be incorrect. Suggest to either "
        "force rescaling, or include hole-particle interaction.\n");
  }
  if (force_rescale && hole_particle) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Should not force rescaling of V(r) if also subtracting hole "
               "particle (self-interaction) potential.\n");
  }
  std::cout << "\n";

  // Spherical Bessel lookup table
  std::cout << "Filling jL spherical Bessel table.." << std::flush;
  const SphericalBessel::JL_table jK_tab(Kmax + 1, qgrid.r(), wf.grid().r());
  std::cout << "..done\n" << std::flush;

  LinAlg::Matrix<double>        // Matrices for each V(E,q) form factor
      Q_Phi(E_steps, q_steps),  // Vector: temporal
      Q_E(E_steps, q_steps),    // Vector: electric
      Q_M(E_steps, q_steps),    // Vector: magnetic
      Q_L(E_steps, q_steps),    // Vector: longitudinal
      Q_Phi5(E_steps, q_steps), // Axial: temporal
      Q_E5(E_steps, q_steps),   // Axial: electric
      Q_M5(E_steps, q_steps),   // Axial: magnetic
      Q_L5(E_steps, q_steps),   // Axial: longitudinal
      Q_S(E_steps, q_steps),    // Scalar
      Q_S5(E_steps, q_steps),   // Pseudo-scalar

      // Cross terms
      Q_Z(E_steps, q_steps),
      Q_X(E_steps, q_steps),
      Q_Y(E_steps, q_steps),

      // Checking NR limits
      Q_Phi_nr(E_steps, q_steps),   // Vector: temporal (NR limit)
      Q_E_nr(E_steps, q_steps),    // Vector: electric (NR limit)
      Q_M_nr(E_steps, q_steps),    // Vector: magnetic (NR limit)
      Q_L_nr(E_steps, q_steps),    // Vector: longitudinal (NR Limit)
      Q_Phi5_nr(E_steps, q_steps),   // Axial: temporal (NR limit)
      Q_E5_nr(E_steps, q_steps),    // Axial: electric (NR limit)
      Q_M5_nr(E_steps, q_steps),    // Axial: magnetic (NR limit)
      Q_L5_nr(E_steps, q_steps),   // Axial: longitudinal (NR Limit)
      Q_S_nr(E_steps, q_steps),     // Scalar (NR limit)
      Q_S5_nr(E_steps, q_steps);     // Pseudo-Scalar (NR limit)


  //-------------------------------------------------------------------------
  std::cout << "\nCalculating ionisation factors:\n";
  for (const auto &Fa : wf.core()) {
    std::cout << Fa << ", " << std::flush;

    // First lowest energy able to ionise Fa:
    const auto idE_0 = std::size_t(std::distance(
        Egrid.begin(), std::find_if(Egrid.begin(), Egrid.end(),
                                    [&Fa](auto e) { return e > -Fa.en(); })));

#pragma omp parallel for
    for (std::size_t iE = idE_0; iE < Egrid.size(); ++iE) {
      const auto omega = Egrid.at(iE);

      const auto ec = omega + Fa.en();
      assert(ec > 0.0);

      const int l = Fa.l();
      const int lc_max = l + Kmax;
      const int lc_min = std::max(l - Kmax, 0);

      ContinuumOrbitals cntm(wf.vHF());
      cntm.solveContinuumHF(ec, lc_min, lc_max, &Fa, force_rescale,
                            hole_particle, force_orthog);

      for (std::size_t iq = 0; iq < qgrid.size(); ++iq) {
        for (int k = Kmin; k <= Kmax; ++k) {

          const auto tkp1_x = (2.0 * k + 1.0) * Fa.occ_frac();
          // Use qc for as expected for "omega" in operators
          const auto qc = qgrid.at(iq) * PhysConst::c;

          const auto Phik = DiracOperator::Phik_w(wf.grid(), k, qc, &jK_tab);
          const auto Ek = DiracOperator::Ek_w(wf.grid(), k, qc, &jK_tab);
          const auto Mk = DiracOperator::Mk_w(wf.grid(), k, qc, &jK_tab);
          const auto Lk = DiracOperator::Lk_w(wf.grid(), k, qc, &jK_tab);

          const auto Phi5k = DiracOperator::Phi5k_w(wf.grid(), k, qc, &jK_tab);
          const auto E5k = DiracOperator::E5k_w(wf.grid(), k, qc, &jK_tab);
          const auto M5k = DiracOperator::M5k_w(wf.grid(), k, qc, &jK_tab);
          const auto L5k = DiracOperator::L5k_w(wf.grid(), k, qc, &jK_tab);

          const auto Sk = DiracOperator::Sk_w(wf.grid(), k, qc, &jK_tab);
          const auto S5k = DiracOperator::S5k_w(wf.grid(), k, qc, &jK_tab);

          const auto Phik_nr = DiracOperator::Phi_nr(wf.grid(), k, qc);
          const auto Phi5k_nr = DiracOperator::Phi5_nr(wf.grid(), k, qc);

          const auto Ek_nr = DiracOperator::E_nr(wf.grid(), qc);
          const auto Mk_nr = DiracOperator::M_nr(wf.grid(), qc);
          const auto Lk_nr = DiracOperator::L_nr(wf.grid(), qc);

          const auto E5k_nr = DiracOperator::E5_nr(wf.grid(), qc);
          const auto M5k_nr = DiracOperator::M5_nr(wf.grid(), qc);
          const auto L5k_nr = DiracOperator::L5_nr(wf.grid(), k, qc);

          const auto Sk_nr = DiracOperator::S_nr(wf.grid(), k, qc);
          const auto S5k_nr = DiracOperator::S5_nr(wf.grid(), k, qc);

          for (const auto &Fe : cntm.orbitals) {
            // Vector operators
            if (vectorQ) {
              Q_Phi(iE, iq) += tkp1_x * qip::pow(Phik.reducedME(Fe, Fa), 2);
              Q_E(iE, iq) += tkp1_x * qip::pow(Ek.reducedME(Fe, Fa), 2);
              Q_M(iE, iq) += tkp1_x * qip::pow(Mk.reducedME(Fe, Fa), 2);
              Q_L(iE, iq) += tkp1_x * qip::pow(Lk.reducedME(Fe, Fa), 2);

              if (k == 0) {
                // Small qr / 'dipole' limit
                Q_Phi_nr(iE, iq) += tkp1_x * qip::pow(Phik_nr.reducedME(Fe, Fa), 2);
              }

              if (k == 1) {
                // Small qr / 'dipole' limit
                Q_Phi_nr(iE, iq) += tkp1_x * qip::pow(Phik_nr.reducedME(Fe, Fa), 2);
                Q_E_nr(iE, iq) += tkp1_x * qip::pow(Ek_nr.reducedME(Fe, Fa), 2);
                Q_M_nr(iE, iq) += tkp1_x * qip::pow(Mk_nr.reducedME(Fe, Fa), 2);
                Q_L_nr(iE, iq) += tkp1_x * qip::pow(Lk_nr.reducedME(Fe, Fa), 2);
              }
            }

            // Axial (γ^5) operators
            if (axialQ) {
              Q_Phi5(iE, iq) += tkp1_x * qip::pow(Phi5k.reducedME(Fe, Fa), 2);
              Q_E5(iE, iq) += tkp1_x * qip::pow(E5k.reducedME(Fe, Fa), 2);
              Q_M5(iE, iq) += tkp1_x * qip::pow(M5k.reducedME(Fe, Fa), 2);
              Q_L5(iE, iq) += tkp1_x * qip::pow(L5k.reducedME(Fe, Fa), 2);

              if (k == 0){
                Q_Phi5_nr(iE, iq) += tkp1_x * qip::pow(Phi5k_nr.reducedME(Fe, Fa), 2);
                Q_L5_nr(iE, iq) += tkp1_x * qip::pow(L5k_nr.reducedME(Fe, Fa), 2);
              }

              if (k == 1){
                // Small qr / 'dipole' limit
                Q_Phi5_nr(iE, iq) += tkp1_x * qip::pow(Phi5k_nr.reducedME(Fe, Fa), 2);
                Q_E5_nr(iE, iq) += tkp1_x * qip::pow(E5k_nr.reducedME(Fe, Fa), 2);
                Q_M5_nr(iE, iq) += tkp1_x * qip::pow(M5k_nr.reducedME(Fe, Fa), 2);
                Q_L5_nr(iE, iq) += tkp1_x * qip::pow(L5k_nr.reducedME(Fe, Fa), 2);
              }
            }

            // Scalar and Pseudoscalar
            if (scalarQ){
              Q_S(iE, iq) += tkp1_x * qip::pow(Sk.reducedME(Fe, Fa), 2);

              if (k == 1){
                Q_S_nr(iE, iq) += tkp1_x * qip::pow(Sk_nr.reducedME(Fe, Fa), 2);
              }

              if (k == 0){
                Q_S_nr(iE, iq) += tkp1_x * qip::pow(Sk_nr.reducedME(Fe, Fa), 2);
              }
            }

            if (pseudoscalarQ){
              Q_S5(iE, iq) += tkp1_x * qip::pow(S5k.reducedME(Fe, Fa), 2);

              if (k == 1){
                Q_S5_nr(iE, iq) += tkp1_x * qip::pow(S5k_nr.reducedME(Fe, Fa), 2);
              }

              if (k == 0){
                Q_S5_nr(iE, iq) += tkp1_x * qip::pow(S5k_nr.reducedME(Fe, Fa), 2);
              }

            }

            if (crossQ){

              Q_Z(iE, iq) += tkp1_x * (E5k.reducedME(Fe, Fa)*Mk.reducedME(Fe, Fa) - Ek.reducedME(Fe, Fa)*M5k.reducedME(Fe, Fa));
              Q_X(iE, iq) += tkp1_x * Phik.reducedME(Fe, Fa)*Lk.reducedME(Fe, Fa);
              Q_Y(iE, iq) += tkp1_x * Phi5k.reducedME(Fe, Fa)*L5k.reducedME(Fe, Fa);

            }

          }
        }
      }
    }
  }
  std::cout << "\ndone\n\n";

  //----------------------------------------------------------------------------

  // Output format:

  const auto print_all = input.get("print_all", false);

  const auto toutput =
      input.get<std::vector<std::string>>("output_format", {"gnuplot"});
  // Use vector, since allow outputting in multiple formats
  std::vector<Kion::OutputFormat> output_formats;
  using namespace qip;
  for (auto &each : toutput) {
    const auto output =
        ci_wc_compare(each, "gnuplot") ? Kion::OutputFormat::gnuplot :
        ci_wc_compare(each, "gnu*_E")  ? Kion::OutputFormat::gnuplot_E :
        ci_wc_compare(each, "xyz")     ? Kion::OutputFormat::xyz :
        ci_wc_compare(each, "mat*")    ? Kion::OutputFormat::matrix :
                                         Kion::OutputFormat::Error;
    if (output == Kion::OutputFormat::Error) {
      fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
      fmt::print(
          "Output Format option: `{}' unknown. Options are: gnuplot, xyz, "
          "matrix\n",
          each);
    } else {
      // only add if not already in list
      if (std::find(output_formats.begin(), output_formats.end(), output) ==
          output_formats.end())
        output_formats.push_back(output);
    }
  }
  if (output_formats.empty()) {
    fmt2::styled_print(fg(fmt::color::red), "\nFail 300: ");
    fmt::print("No output formats? No output will be generated!\n");
    // return;
    output_formats.push_back(Kion::OutputFormat::gnuplot);
  }

  // Units (only for gnuplot-style):
  const auto tunits = input.get<std::string>("units", "Particle");
  auto units = ci_compare(tunits, "particle") ? Kion::Units::Particle :
               ci_compare(tunits, "atomic")   ? Kion::Units::Atomic :
                                                Kion::Units::Error;
  if (units == Kion::Units::Error) {
    fmt2::styled_print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Output units option: `{}' unknown. Options are: Particle, "
               "Atomic. Defaulting to atomic\n",
               tunits);
    units = Kion::Units::Atomic;
  }
  std::cout << "_gnu output file (if requested) will use: " << tunits
            << " units\n";
  std::cout << "(_mat and _xyz output files always use atomic units)\n";

  //----------------------------------------------------------------------------

  const auto method =
      HF::parseMethod_short(wf.vHF()->method()) + (hole_particle ? "_hp" : "") +
      (force_orthog ? "_orth" : "") + (force_rescale ? "_rescale" : "");

  std::string prefix = wf.identity() + "_" + method + "_" +
                       std::to_string(Kmin) + "-" + std::to_string(Kmax) + "_";

  std::cout << "\nWriting to files: " << prefix << "...\n";
  if (vectorQ) {
    Kion::write_to_file(output_formats, Q_Phi, Egrid, qgrid, prefix + "Phi", 8,
                        units);
    Kion::write_to_file(output_formats, Q_E + Q_M + Q_L, Egrid, qgrid,
                        prefix + "V", 8, units);
    if (print_all) {
      Kion::write_to_file(output_formats, Q_E, Egrid, qgrid, prefix + "V^e", 8,
                          units);
      Kion::write_to_file(output_formats, Q_M, Egrid, qgrid, prefix + "V^m", 8,
                          units);
      Kion::write_to_file(output_formats, Q_L, Egrid, qgrid, prefix + "V^l", 8,
                          units);
      Kion::write_to_file(output_formats, Q_E_nr, Egrid, qgrid, prefix + "V^e_nr", 8,
                          units);
      Kion::write_to_file(output_formats, Q_M_nr, Egrid, qgrid, prefix + "V^m_nr", 8,
                          units);
      Kion::write_to_file(output_formats, Q_L_nr, Egrid, qgrid, prefix + "V^l_nr", 8,
                          units);
      Kion::write_to_file(output_formats, Q_Phi_nr, Egrid, qgrid, prefix + "Phi_nr", 8,
                        units);
    }
  }
  if (axialQ) {
    Kion::write_to_file(output_formats, Q_Phi5, Egrid, qgrid, prefix + "Phi5",
                        8, units);
    Kion::write_to_file(output_formats, Q_E5 + Q_M5 + Q_L5, Egrid, qgrid,
                        prefix + "V5", 8, units);
    if (print_all) {
      Kion::write_to_file(output_formats, Q_E5, Egrid, qgrid, prefix + "V5^e",
                          8, units);
      Kion::write_to_file(output_formats, Q_M5, Egrid, qgrid, prefix + "V5^m",
                          8, units);
      Kion::write_to_file(output_formats, Q_L5, Egrid, qgrid, prefix + "V5^l",
                          8, units);
      Kion::write_to_file(output_formats, Q_E5_nr, Egrid, qgrid, prefix + "V5^e_nr",
                          8, units);
      Kion::write_to_file(output_formats, Q_M5_nr, Egrid, qgrid, prefix + "V5^m_nr",
                          8, units);
      Kion::write_to_file(output_formats, Q_L5_nr, Egrid, qgrid, prefix + "V5^l_nr",
                          8, units);
      Kion::write_to_file(output_formats, Q_Phi5_nr, Egrid, qgrid, prefix + "Phi5_nr", 8,
                        units);
    }
  }
  if (scalarQ) {
    Kion::write_to_file(output_formats, Q_S, Egrid, qgrid, prefix + "S", 8,
                        units);
    Kion::write_to_file(output_formats, Q_S_nr, Egrid, qgrid, prefix + "S_nr", 8,
                        units);
  }
  if (pseudoscalarQ) {
    Kion::write_to_file(output_formats, Q_S5, Egrid, qgrid, prefix + "S5", 8,
                        units);
    Kion::write_to_file(output_formats, Q_S5_nr, Egrid, qgrid, prefix + "S5_nr", 8,
                        units);
  }

  if (crossQ) {
    Kion::write_to_file(output_formats, Q_Z, Egrid, qgrid, prefix + "Z", 8,
                        units);
    Kion::write_to_file(output_formats, Q_X, Egrid, qgrid, prefix + "X", 8,
                        units);
    Kion::write_to_file(output_formats, Q_Y, Egrid, qgrid, prefix + "Y", 8,
                        units);
  }
}

} // namespace Module