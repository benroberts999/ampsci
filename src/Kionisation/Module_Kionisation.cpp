#include "Kionisation/Module_Kionisation.hpp"
#include "DiracOperator/DiracOperator.hpp"
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
#include <iostream>
#include <memory>

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
       {"label", "label for output files"},
       {"method", "Method for continuum: standard (rel. HF), zeff (for "
                  "continuum), approx (step function). [standard]"},
       {"rpa", "bool. Use RPA (rel HF, with first-order RPA)? [false]"},
       {"subtract_1", "Replace e^(iqr) -> e^(iqr)-1 [false]"},
       {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
       {"hole_particle", "Subtract Hartree-Fock self-interaction (account for "
                         "hole-particle interaction) [true]"},
       {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
       {"coupling", "Vector, Scalar (g0), Pseudovector (g5), Pseudoscalar "
                    "(g0g5) [Vector]"},
       {"output_format", "List: Format for output. List any of: gnuplot, xyz, "
                         "matrix (comma-separated). [gnuplot]\n\n"
                         "xyz: For easy 2D interpolation. list formmated with "
                         "each row 'E q K(E,q)'\n"
                         "gnuplot: For easy plotting. Each column is new E\n"
                         "matrix: Outputs entire matrix in table form. E and q "
                         "grids printed prior\n"},
       {"each_state", "bool. If true, will output K(E,q) for each "
                      "(accessible) core-state [false]"},
       {"units", "Units for output: Particle (keV/MeV) or Atomic (E_H,1/a0) "
                 "[Particle]"}});
  if (input.has_option("help"))
    return;

  //----------------------------------------------------------------------------

  // Read in energy-deposit/momentum-exchange input options:
  const auto [Emin_keV, Emax_keV] = input.get("E_range", std::array{0.1, 0.1});
  const auto E_steps = input.get<std::size_t>("E_steps", 1);
  const auto Emin_au = Emin_keV * UnitConv::Energy_keV_to_au;
  const auto Emax_au =
      Emax_keV < Emin_keV ? Emin_au : Emax_keV * UnitConv::Energy_keV_to_au;

  std::cout << "\nSummary of inputs:\n";
  fmt::print(
      "Energy  : [{:.2f}, {:.2f}] keV  = [{:.1f}, {:.1f}] au, in {} steps\n",
      Emin_keV, Emax_keV, Emin_au, Emax_au, E_steps);

  const auto [qmin_keV, qmax_keV] =
      input.get("q_range", std::array{0.01, 0.01});
  const auto q_steps = input.get<std::size_t>("q_steps", 1);
  const auto qmin_au = qmin_keV * UnitConv::Momentum_MeV_to_au;
  const auto qmax_au = qmax_keV * UnitConv::Momentum_MeV_to_au;
  const auto max_L = input.get("max_L", 6);
  const auto label = input.get("label", std::string{""});

  fmt::print(
      "Momentum: [{:.3f}, {:.3f}] MeV = [{:.1f}, {:.1f}] au, in {} steps\n",
      qmin_keV, qmax_keV, qmin_au, qmax_au, q_steps);

  // Set up the E and q grids
  const Grid Egrid({E_steps, Emin_au, Emax_au, 0, GridType::logarithmic});
  const Grid qgrid({q_steps, qmin_au, qmax_au, 0, GridType::logarithmic});

  // Check to see if grid is reasonable for maximum energy:
  Kion::check_radial_grid(Emax_au, qmax_au, wf.grid());

  //----------------------------------------------------------------------------
  // Read in and parse options:

  fmt::print("\nMax L   : {}  (multipolarity in e^iqr expansion)\n", max_L);

  if (!label.empty()) {
    fmt::print("Label   : {}\n", label);
  }

  // Method for cntm states:
  using qip::ci_compare;          // case-insensitive string comparison
  using qip::ci_wildcard_compare; // case-insensitive string comparison with *
  const auto tmethod = input.get<std::string>("method", "standard");
  const auto method = ci_compare(tmethod, "standard") ? Kion::Method::Standard :
                      ci_compare(tmethod, "approx")   ? Kion::Method::Approx :
                      ci_compare(tmethod, "zeff")     ? Kion::Method::Zeff :
                                                        Kion::Method::Error;
  const auto use_Zeff_cont = method == Kion::Method::Zeff;
  if (method == Kion::Method::Error) {
    fmt::print(fg(fmt::color::red), "\nError 104: ");
    fmt::print("Method option: `{}' unknown. Options are: standard, approx, "
               "zeff\n",
               tmethod);
    return;
  }

  // Use RPA?
  const bool use_rpa = input.get("rpa", false);
  if (use_rpa) {
    std::cout << "Using RPA (first-order), with basis: "
              << DiracSpinor::state_config(wf.basis()) << "\n";
  }

  // sanity check for RPA - meaningless unless we used HF method
  if (use_rpa && wf.vHF() && wf.vHF()->method() != HF::Method::HartreeFock) {
    fmt::print(fg(fmt::color::red), "\nError 113: ");
    fmt::print("It's only meaningful to include RPA if we start with "
               "HartreeFock method (for core states).\n");
    return;
  }
  if (use_rpa && wf.basis().empty()) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("RPA method required a basis to work. See `ampsci -a "
               "Basis`\nRPA will evaluate to zero without a basis\n\n");
  }
  std::cout << "Using   : " << tmethod << " method\n";

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
  if (use_rpa && use_Zeff_cont) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Using RPA with Z_eff model doesn't make sense.\n");
  }

  // Perform checks, print possible warnings
  if (force_rescale && !(force_orthog || subtract_1)) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Force rescale will ruin orthogonality; suggest use "
               "force_orthog or subtract_1\n");
  }
  if (use_Zeff_cont && !(force_orthog || subtract_1)) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("using Zeff for continuum will ruin orthogonality; suggest use "
               "force_orthog or subtract_1\n");
  }
  if (!force_rescale && !hole_particle && !use_Zeff_cont && wf.Zion() == 0) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print(
        "Long-range behaviour of V(r) may be incorrect. Suggest to either "
        "force rescaling, or include hole-particle interaction.\n");
  }
  if (force_rescale && hole_particle) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
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
    fmt::print(fg(fmt::color::red), "\nError 212: ");
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
                                             std::size_t(max_L));
  } else if (coupling == Kion::Coupling::Scalar) {
    jl = std::make_unique<DiracOperator::g0jL>(wf.grid(), qgrid,
                                               std::size_t(max_L));
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
  const std::string method_text = method == Kion::Method::Standard ? "" :
                                  method == Kion::Method::Approx   ? "aprx_" :
                                  method == Kion::Method::Zeff     ? "zeff_" :
                                                                     "???_";

  const std::string rpa_text =
      use_rpa ? "rpa" + DiracSpinor::state_config(wf.basis()) + "_" : "";

  const std::string coupling_text =
      coupling == Kion::Coupling::Vector       ? "v" :
      coupling == Kion::Coupling::Scalar       ? "s" :
      coupling == Kion::Coupling::PseudoVector ? "pv" :
      coupling == Kion::Coupling::PseudoScalar ? "ps" :
                                                 "???";
  // doesn't include suffix
  std::string oname = "K_" + method_text + rpa_text + wf.atomicSymbol() + "_";
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
        ci_wildcard_compare(each, "gnu*") ? Kion::OutputFormat::gnuplot :
        ci_wildcard_compare(each, "xyz")  ? Kion::OutputFormat::xyz :
        ci_wildcard_compare(each, "mat*") ? Kion::OutputFormat::matrix :
                                            Kion::OutputFormat::Error;
    if (output == Kion::OutputFormat::Error) {
      fmt::print(fg(fmt::color::orange), "\nWarning: ");
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
    fmt::print(fg(fmt::color::red), "\nFail 300: ");
    fmt::print("No output formats? No output will be generated!\n");
    return;
  }

  // Units:
  const auto tunits = input.get<std::string>("units", "particle");
  auto units = ci_compare(tunits, "particle") ? Kion::Units::Particle :
               ci_compare(tunits, "atomic")   ? Kion::Units::Atomic :
                                                Kion::Units::Error;
  if (units == Kion::Units::Error) {
    fmt::print(fg(fmt::color::orange), "\nWarning: ");
    fmt::print("Output units option: `{}' unknown. Options are: Particle, "
               "Atomic. Defaulting to particle\n",
               tunits);
    units = Kion::Units::Particle;
  }
  std::cout << "Output files will use: " << tunits << " units\n";
  const int num_output_digits = 5;

  //----------------------------------------------------------------------------

  std::cout << "Calculating K(E,q) - ionisation factor\n" << std::flush;

  // Kion stored K(dE, q)
  LinAlg::Matrix<double> Kion(Egrid.num_points(), qgrid.num_points());

  if (method != Kion::Method::Approx) {
    std::cout << "Using   : " << tmethod << " method\n";
    for (const auto &Fnk : wf.core()) {
      const auto accessible = std::abs(Fnk.en()) < Emax_au;
      if (!accessible)
        continue;
      std::cout << Fnk << ", " << std::flush;
      const auto K_nk = Kion::calculateK_nk(
          wf.vHF(), Fnk, max_L, Egrid, jl.get(), subtract_1, force_rescale,
          hole_particle, force_orthog, use_Zeff_cont, use_rpa, wf.basis());
      if (write_each_state) {
        const auto oname_nk = oname + "_" + Fnk.shortSymbol();
        std::cout << "Written to file: " << oname_nk << "\n";
        Kion::write_to_file(output_formats, K_nk, Egrid, qgrid, oname_nk,
                            num_output_digits, units);
      }
      Kion += K_nk;
    }
    std::cout << "\n";
  } else {
    std::cout << "Using   : approx (step-function) method\n";
    const auto K_approx = Kion::calculateK_nk_approx(
        wf.vHF(), wf.core(), max_L, jl.get(), subtract_1, force_rescale,
        hole_particle, force_orthog, use_Zeff_cont, use_rpa, wf.basis());

    Kion::write_approxTable_to_file(K_approx, wf.core(), qgrid, oname,
                                    num_output_digits, units);

    // convert to "standard" form, for easy comparison
    Kion = Kion::convert_K_nk_approx_to_std(K_approx, Egrid, wf.core());
  }

  std::cout << "\nWritten to file: " << oname << "\n";
  Kion::write_to_file(output_formats, Kion, Egrid, qgrid, oname,
                      num_output_digits, units);

  std::cout << '\n';
}

} // namespace Module