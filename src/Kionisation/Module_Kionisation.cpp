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
#include "qip/Methods.hpp"
#include "qip/String.hpp"
#include <cassert>
#include <iostream>
#include <memory>
//
#include "ExternalField/DiagramRPA.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"

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
       "List (2). Minimum, maximum energy transfer (dE), in eV [10, 100]"},
      {"E_steps", "Numer of steps along dE grid (logarithmic grid) [1]"},
      {"oname", "oname"},
      {"ec_cut", "Cut-off (in au) for continuum energy. [100.0]"},
      {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
      {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                        "hole-particle interaction) [true]"},
      {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
  });
  if (input.has_option("help")) {
    return;
  }

  auto [Emin_eV, Emax_eV] = input.get("E_range", std::array{10.0, 100.0});
  auto E_steps = input.get<std::size_t>("E_steps", 1);
  if (E_steps <= 1) {
    E_steps = 1;
    Emax_eV = Emin_eV;
  }
  const auto Emin_au = Emin_eV / PhysConst::Hartree_eV;
  const auto Emax_au =
      Emax_eV < Emin_eV ? Emin_au : Emax_eV / PhysConst::Hartree_eV;

  std::cout << "\nSummary of inputs:\n";
  fmt::print(
      "Energy  : [{:.1f}, {:.1f}] eV  = [{:.1e}, {:.1e}] au, in {} steps\n",
      Emin_eV, Emax_eV, Emin_au, Emax_au, E_steps);
  const auto ec_cut = input.get("ec_cut", 100.0);

  const auto K = input.get("K", 1);
  const auto label = input.get("label", std::string{""});

  const auto force_orthog = input.get("force_orthog", true);
  const auto force_rescale = input.get("force_rescale", false);
  const auto hole_particle = input.get("hole_particle", true);

  // Set up the E and q grids
  const Grid Egrid({E_steps, Emin_au, Emax_au, 0, GridType::logarithmic});

  int k = K;

  const auto E1 = DiracOperator::E1(wf.grid());

  auto oname = input.get("oname", std::string{"out.txt"});
  std::ofstream out_file(oname);

  auto Ek = DiracOperator::Ekv_omega(wf.grid(), k, wf.alpha(), 0.0, true);
  ExternalField::DiagramRPA dV0(&Ek, wf.basis(), wf.vHF(), wf.identity());
  ExternalField::DiagramRPA dV(&Ek, wf.basis(), wf.vHF(), wf.identity());

  int count = 0;
  for (const auto omega : Egrid.r()) {
    std::cout << count++ << " " << omega << " " << omega * PhysConst::Hartree_eV
              << "\n";

    const auto Ksigma = 4.0 * M_PI * M_PI * wf.alpha() * PhysConst::aB_cm *
                        PhysConst::aB_cm * omega;

    Ek.updateFrequency(omega);
    // dV0.update_t0s(&Ek); // required??
    // dV.solve_core(0.0);
    // dV.solve_core(omega);

    const auto EkL =
        DiracOperator::Ek_omega(wf.grid(), k, wf.alpha(), omega, true);

    const auto Mk =
        DiracOperator::Mk_omega(wf.grid(), k, wf.alpha(), omega, true);

    const auto E1v = DiracOperator::E1v(wf.alpha(), omega);
    const auto M1 = DiracOperator::M1(wf.grid(), wf.alpha(), omega);

    const auto ialpha = DiracOperator::ialpha();

    double Q_E = 0.0;
    double Q_El = 0.0;
    double Q_E1 = 0.0;
    double Q_E1v = 0.0;
    double Q_ialpha = 0.0;

    // First, loop through and just find list of what we shall do.
    // THEN parellelise over that!

    std::vector<std::size_t> iclist;
    for (std::size_t i = 0; i < wf.core().size(); ++i) {
      const auto ec = omega + wf.core()[i].en();
      if (ec < 0.0)
        continue;
      iclist.push_back(i);
    }

#pragma omp parallel for reduction(+ : Q_E, Q_El, Q_E1, Q_E1v, Q_ialpha)
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

        const auto f_mp = 1.0 / (wf.alpha() * omega) / (wf.alpha() * omega);
        const auto f_1o3 = 1.0 / 3.0;

        const auto rme_E1 = E1.reducedME(Fe, Fa);
        Q_E1 += f_1o3 * rme_E1 * rme_E1;

        const auto rme_E1v = E1v.reducedME(Fe, Fa);
        Q_E1v += f_1o3 * rme_E1v * rme_E1v;

        const auto rme_ialpha = ialpha.reducedME(Fe, Fa);
        Q_ialpha += f_mp * f_1o3 * rme_ialpha * rme_ialpha;

        const auto tkp1 = 2.0 * k + 1.0;
        const auto pol_av = 1.0 / 2.0;
        const auto rme_E = Ek.reducedME(Fe, Fa);
        Q_E += f_mp * rme_E * rme_E * tkp1 * pol_av;

        const auto rme_El = EkL.reducedME(Fe, Fa);
        Q_El += f_mp * rme_El * rme_El * tkp1 * pol_av;
      }
    }
    // std::cout << "\n";

    out_file << omega * PhysConst::Hartree_eV / 1e6 << " " << Q_E * Ksigma
             << " " << Q_El * Ksigma << " " << Q_E1 * Ksigma << " "
             << Q_E1v * Ksigma << " " << Q_ialpha * Ksigma << "\n";
  }
}

} // namespace Module