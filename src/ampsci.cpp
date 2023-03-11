#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/InputBlock.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp" //for 'ExtraPotential'
#include "Modules/modules_list.hpp"
#include "Modules/runModules.hpp"
#include "Physics/include.hpp"
#include "Physics/periodicTable.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include "qip/omp.hpp"
#include "version/EasterEgg.hpp"
#include "version/version.hpp"
#include <iostream>
#include <memory>
#include <string>

//! General usage instructions
const std::string ampsci_help{R"(
ampsci
Atomic Many-body Perturbation theory in the Screened Coulomb Interaction. 
Benjamin M. Roberts (https://broberts.io/), University of Queensland, Australia.

OPTIONS
    <filename>
        Runs ampsci taking options specified in file "filename" (eg, ./ampsci filename). See documentation (or option -a) for input file format options.
        Example:
        ./ampsci input.in
            -Runs ampsci taking input options from file 'input.in'
    
    <At> <Core> <Valence>
        For quick and simple HF calculation. If core is not given, guesses core configuration and runs using V^N approximation.
        Examples:
        ./ampsci Cs
            - Runs ampsci for Cs using Hartree Fock (V^N) approximation
        ./ampsci Cs [Xe] 6sd5d
            - Runs ampsci for Cs using Hartree Fock with Xe-like core and valence states up to n=6 for s,p-states and n=5 for d-states

    -v (--version)
        Prints ampsci version (and git commit) details

    -l (--libs, --libraries)
        Prints ampsci version details for libaries

    -h (--help, -?)
        Print help info, including some detail on input options

    -a <BlockName> (--ampsci)
        Prints list of available top-level ampsci options. BlockName is optional; if given it will print options for given ampsci Block. You may list any number of blocks (space separated)
        Example:
        ./ampsci -a Atom HartreeFock

    -m <ModuleName> (--modules)
        Prints list of available Modules. ModuleName is optional; if given, will list avaiable options for that Module
        Example:
        ./ampsci -m MatrixElements

    -o <OperatorName> (--operators)
        Prints list of available operators. OperatorName is optional; if given, will list avaiable options for that operator (most operators take no options).
        Example:
        ./ampsci -o E1

    -p <Atom> <Isotope> (--periodicTable)
        Prints textual periodic table with electronic + nuclear information.
        Atom and Isotope are optional; if given, will print info for that isotope. Atom should be atomic symbol (eg Cs), or Z (55).
        If Isotope is blank, will print for 'default' isotope. Can also list 'all' known isotope info
        Examples:
        ./ampsci -p Cs
        ./ampsci -p Cs 131
        ./ampsci -p Cs all
    
    -c (--constants)
        Prints some handy physical constants
)"};

//! Calculates wavefunction and runs optional modules
void ampsci(const IO::InputBlock &input);

//==============================================================================
int main(int argc, char *argv[]) {
  using namespace std::string_literals;

  // Parse input text into strings:
  const std::string input_text = (argc > 1) ? argv[1] : "";
  const std::string core_text = (argc > 2) ? argv[2] : "";
  const std::string valence_text = (argc > 3) ? argv[3] : "";

  // check for special commands
  if (input_text == "") {
    std::cout << ampsci_help << '\n';
    return 0;
  } else if (input_text == "-v" || input_text == "--version") {
    std::cout << "AMPSCI v: " << version::version() << '\n';
    std::cout << "Libraries:\n" << version::libraries() << '\n';
    std::cout << "Compiled: " << version::compiled() << '\n';
    std::cout << "Benjamin M. Roberts (https://broberts.io/), University of "
                 "Queensland, Australia\n";
    return 0;
  } else if (input_text == "-l" || input_text.substr(0, 5) == "--lib") {
    std::cout << "Libraries:\n" << version::libraries() << '\n';
    return 0;
  } else if (input_text == "-h" || input_text == "--help" ||
             input_text == "-?") {
    std::cout << ampsci_help << '\n';
    return 0;
  } else if (input_text == "-m" || input_text == "--modules") {
    std::cout << "Available modules: \n";
    Module::list_modules();
    const std::string module_name = (argc > 2) ? argv[2] : "";
    if (!module_name.empty()) {
      // make an arbitrary required as input to runModule
      Wavefunction wf{{1, 1.0, 1.0}, {1, 1}};
      // run the module, with option 'help' set. This will trigger the helper
      // to print the details for the available options in that module
      Module::runModule(IO::InputBlock{"Module::"s + module_name, {"help;"}},
                        wf);
    }
    return 0;
  } else if (input_text == "-o" || input_text == "--operators") {
    std::cout << "Available operators: \n";
    DiracOperator::list_operators();
    const std::string op_name = (argc > 2) ? argv[2] : "";
    if (!op_name.empty()) {
      Wavefunction wf{{1, 1.0, 1.0}, {1, 1}};
      DiracOperator::generate(op_name, IO::InputBlock{op_name, {"help;"}}, wf);
    }
    return 0;
  } else if (input_text == "-a" || input_text == "--ampsci") {
    auto temp_input = IO::InputBlock{"ampsci", {"help;"}};
    for (int i_in = 2; i_in < argc; ++i_in) {
      const std::string block_name = (argc > i_in) ? argv[i_in] : "";
      if (!block_name.empty()) {
        temp_input.add(IO::InputBlock{block_name, {"help;"}});
      }
    }
    ampsci(temp_input);
    return 0;
  } else if (input_text == "-p" || input_text == "--periodicTable") {
    std::string z_str = (argc > 2) ? argv[2] : "";
    std::string a_str = (argc > 3) ? argv[3] : "";
    AtomData::periodicTable(z_str, a_str);
    return 0;
  } else if (input_text == "-c" || input_text == "--constants") {
    AtomData::printConstants();
    return 0;
  } else if (!input_text.empty() && input_text.front() == '-') {
    std::cout << "Unrecognised option: " << input_text << '\n';
    std::cout << ampsci_help << '\n';
    return 0;
  }

  // Print git/version info to screen:
  std::cout << '\n';
  IO::print_line();
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Parallel: " << qip::omp_details() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';
  std::cout << "Run time: " << IO::time_date() << '\n';

  // If we are not given a valid input text file, assume input is in form:
  // <At> <core> <valence> (e.g., "Cs [Xe] 6sp")
  // All optional, but must appear in order. ("Core" may be skipped for H)
  const auto atom = AtomData::atomicSymbol(AtomData::atomic_Z(input_text));
  const auto core =
      atom == "H" ? "" : (core_text == "" ? "[" + atom + "]" : core_text);
  // allow core to be skipped for Hydrogen
  const auto valence = (atom == "H" && argc == 3) ? core_text : valence_text;
  const std::string default_input = "Atom{Z=" + atom + ";}" +
                                    "HartreeFock { core = " + core +
                                    "; valence = " + valence + ";}";

  // nb: std::filesystem not available in g++-7 (getafix version)
  const auto fstream = std::fstream(input_text);
  auto input = fstream.good() ? IO::InputBlock("ampsci", fstream) :
                                IO::InputBlock("ampsci", default_input);

  // Run program. Add option to run multiple times
  ampsci(input);

  return 0;
}

//==============================================================================
//==============================================================================
//==============================================================================

void ampsci(const IO::InputBlock &input) {
  IO::ChronoTimer timer("\nampsci");

  using namespace std::string_literals;
  // vector overloads:
  using namespace qip::overloads;

  //check if running in 'help' mode: if so, suppress some output
  // (i.e., check if 'help' has been included, even if not set)
  const auto help_mode = input.has_option("help");
  if (!help_mode) {
    std::cout << '\n';
    IO::print_line();
    input.print();
  }

  // Top-level input blocks
  input.check(
      {{"",
        "Format for descriptions are:\n Description [default_value]\n Blocks "
        "end with '{}', options end with ';'"},
       {"Atom{}", "Which atom to run for"},
       {"Nucleus{}", "Set nuclear parameters"},
       {"Grid{}", "Set radial grid parameters"},
       {"HartreeFock{}", "Options for Solving atomic system"},
       {"RadPot{}", "Inlcude QED radiative potential"},
       {"ExtraPotential{}", "Include an extra potential. Rarely used."},
       {"Basis{}", "Basis of HF eigenstates used for MBPT"},
       {"Correlations{}", "Options for correlations"},
       {"Spectrum{}",
        "Like basis, but includes correlations. Used for sum-over-states"},
       {"Module::*{}", "Run any number of modules (* -> module name). `ampsci "
                       "-m` to see available modules"}});
  input.check({"EasterEgg"}, {{"", EasterEgg::get_egg()}});

  // Atom: Get + setup atom parameters
  input.check({"Atom"},
              {{"Z", "Atomic number or symbol (e.g., 55 or Cs). [H]"},
               {"A", "Atomic mass number, for nuclear parameters including "
                     "finite nuclear size. Default based on Z."},
               {"varAlpha2",
                "Fractional variation of the fine-structure constant, alpha^2: "
                "d(a^2)/a_0^2. Use to enforce the non-relativistic limit "
                "(c->infinity => alpha->0), or calculate sensitivity to "
                "variation of alpha. [1.0]"}});

  const auto atom_Z = AtomData::atomic_Z(input.get({"Atom"}, "Z", "H"s));
  const auto atom_A = input.get({"Atom"}, "A", AtomData::defaultA(atom_Z));
  const auto var_alpha = [&]() {
    const auto varAlpha2 = input.get({"Atom"}, "varAlpha2", 1.0);
    // cannot explicitely set alpha to zero - so make it very small
    return (varAlpha2 > 0.0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Nucleus: Get + setup nuclear parameters
  input.check({"Nucleus"},
              {{"", "Options for nuclear potential (finite nuclear size). All "
                    "are optional. Default is a Fermi-like nucleus, with "
                    "parameters chosen according to isotope (see Atom{A;})"},
               {"rrms", "Root-mean-square charge radius, in fm "
                        "[default depends on Z and A]"},
               {"c", "Half-density radius, in fm (will over-ride rms) [default "
                     "depends on Z and A]"},
               {"t", "Nuclear skin thickness, in fm [2.3]"},
               {"type", "Fermi, spherical, pointlike, Gaussian [Fermi]"}});

  // Set nuclear type. If given
  const auto nuc_type = input.get({"Nucleus"}, "type", "Fermi"s);
  // Get default nucleus:
  auto nucleus = Nuclear::Nucleus{atom_Z, atom_A, nuc_type};
  // over-ride default options
  const auto rrms = input.get<double>({"Nucleus"}, "rrms");
  const auto t = input.get<double>({"Nucleus"}, "t");
  const auto c_hdr = input.get<double>({"Nucleus"}, "c");
  if (t) {
    nucleus.t() = *t;
  }
  if (rrms) {
    nucleus.r_rms() = *rrms;
  }
  if (c_hdr) {
    // this will over-ride given rms
    nucleus.r_rms() = Nuclear::rrms_formula_c_t(*c_hdr, nucleus.t());
  }
  // If A or given rrms are zero, explicitely set to pointlike nucleus
  // This isn't required, but makes output more explicit
  if (nucleus.a() == 0.0 || nucleus.r_rms() == 0.0) {
    nucleus.t() = 0.0;
    nucleus.r_rms() = 0.0;
    nucleus.type() = Nuclear::ChargeDistro::point;
  }

  // Grid: Get + setup grid parameters
  input.check(
      {"Grid"},
      {{"", "Options for radial grid (lattice) used for integrations, solving "
            "equations and storing oritals. All relevant quantities are in "
            "units of Bohr radius (aB)."},
       {"r0", "Initial grid point, in aB [1.0e-6]"},
       {"rmax", "Finial grid point [120.0]"},
       {"num_points", "Number of grid points [2000]"},
       {"type", "Type of grid: loglinear, logarithmic, linear [loglinear]"},
       {"b", "Only used for loglinear: grid is ~ logarithmic for r<b, "
             "linear for r>b [rmax/3]"},
       {"du", "du is uniform grid step size; set this instead of "
              "num_points - will override num_points [default set by "
              "num_points]. Rarely used."}});

  // Radial grid. Shared resource used by all wavefunctions/orbitals etc
  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);
  const auto du = input.get<double>({"Grid"}, "du");
  // if num_points = 0, uses du.
  const auto num_points = du ? 0ul : input.get({"Grid"}, "num_points", 2000ul);
  const auto b = input.get({"Grid"}, "b", rmax / 3.0);
  const auto grid_type = (b <= r0 || b >= rmax) ?
                             "logarithmic" :
                             input.get({"Grid"}, "type", "loglinear"s);
  const auto radial_grid = std::make_shared<const Grid>(
      GridParameters{num_points, r0, rmax, b, grid_type, du ? *du : 0});

  // Create wavefunction object
  Wavefunction wf(radial_grid, std::move(nucleus), var_alpha);

  if (!help_mode) {
    std::cout << "\nRunning for " << wf.atom() << '\n';
    if (std::abs(var_alpha - 1.0) > 1.0e-8) {
      if (var_alpha > 1.0e-4) {
        std::cout << "With variation of alpha: ";
      } else {
        std::cout << "Non-relativistic limit: ";
      }
      std::cout << "a/a0 = " << var_alpha
                << " (a/a0)^2 = " << var_alpha * var_alpha << "\n";
    }
    std::cout << wf.nucleus() << '\n'
              << wf.grid().gridParameters() << '\n'
              << "========================================================\n";
  }

  // Parse input for Hartree-Fock
  input.check(
      {"HartreeFock"},
      {{"", "Options for solving lowest-order atomic wavefunction"},
       {"core", "Core configuration. Either list entire core, or use [At] "
                "short-hand. e.g., [He] equivilant to 1s2; [Xe],6s1 equivilant "
                "to [Cs] and to 1s2,2s2,...,5p6,6s1. [blank by default]"},
       {"valence", "e.g., 7sp5d will include valence states up to "
                   "n=7 for s and p, but n=5 for d states. Automatically "
                   "excludes states "
                   "in the core. [blank by default]"},
       {"eps", "HF convergance goal [1.0e-13]"},
       {"method", "HartreeFock, Hartree, KohnSham, Local [HartreeFock]"},
       {"Breit", "Scale for factor for Breit Hamiltonian. Usially 0.0 (no "
                 "Breit) or 1.0 (full Breit), but can take any value. [0.0]"}});

  const auto core = input.get({"HartreeFock"}, "core", "[]"s);
  const auto HF_method = input.get({"HartreeFock"}, "method", "HartreeFock"s);
  const auto eps_HF = input.get({"HartreeFock"}, "eps", 1.0e-13);
  const auto x_Breit = input.get({"HartreeFock"}, "Breit", 0.0);
  const auto valence = input.get({"HartreeFock"}, "valence", ""s);

  // Set up the Hartree Fock potential/method (does not solve)
  // (Must set HF before adding RadPot - but must add RadPot before solving HF)
  wf.set_HF(HF_method, x_Breit, core, eps_HF, true);

  // Forms QED radiative potential, if RadPot{} block is present.
  // Note: input options are parsed inside radiativePotential()
  const auto qed_input = input.getBlock("RadPot");
  if (qed_input != std::nullopt) {
    std::cout << "Including QED into Hartree-Fock core (and valence)\n\n";
    wf.radiativePotential(*qed_input, true, true);
  }

  // Inlcude extra potential. Either read in from text file.
  // or "polarisation operator": V(r) = -0.5/(r^4 + r_cut^4)
  // Note: If read in, it is interpolated onto grid, but NOT extrapolated
  input.check(
      {"ExtraPotential"},
      {{"",
        "Adds an extra potential (to Vnuc), before HF solved. Either effective "
        "polarisation potential, V(r) = -0.5/(r^4 + r_cut^4), or read in from "
        "a file."},
       {"filename", "Read potential from file (r v(r)) - will be interpolated. "
                    "If not given, will use pol. potential"},
       {"r_cut", "Radial cut-off parameter for effective pol. potential [=1]"},
       {"scale",
        "Overall scaling factor for potential is scaled by this value [1]"}});
  const auto include_extra_potential =
      input.getBlock("ExtraPotential") != std::nullopt;
  if (include_extra_potential) {
    const auto ep_fname = input.get({"ExtraPotential"}, "filename", ""s);
    const auto ep_scale = input.get({"ExtraPotential"}, "scale", 1.0);

    if (ep_fname != "") {
      std::cout << "Including extra potential, read in from file: " << ep_fname
                << ", with scale factor: " << ep_scale << "\n";
      const auto &[x, y] = IO::FRW::readFile_xy_PoV(ep_fname);
      if (x.size() == y.size() && x.size() > 1) {
        const auto Vextra = Interpolator::interpolate(x, y, wf.grid().r());
        wf.update_Vnuc(wf.vnuc() + ep_scale * Vextra);
      } else {
        std::cout << "Error 391: Bad file. Not included.\n";
      }
    } else {
      const auto rc = input.get({"ExtraPotential"}, "r_cut", 1.0);
      std::cout << "Adding effective polarisation potential [-0.5 * a * (r^4 + "
                   "rc^4)]: a="
                << ep_scale << ", rc=" << rc << '\n';
      const auto a4 = rc * rc * rc * rc;
      auto dV = [=](auto r) { return -0.5 / (r * r * r * r + a4); };
      const auto dv = qip::apply_to(dV, wf.grid().r());
      wf.update_Vnuc(wf.vnuc() + ep_scale * dv);
    }
  }

  // Solve Hartree equations for the core:
  wf.solve_core(true);

  // Solve for the valence states:
  wf.solve_valence(valence);

  // Output Hartree Fock energies:
  if (!help_mode) {
    std::cout << '\n' << wf.identity() << "-" << wf.Anuc() << '\n';
    wf.printCore();
    printf("E_c = %.6f\n", wf.coreEnergyHF());
    wf.printValence();
  }

  // Construct B-spline basis:
  input.check(
      {"Basis"},
      {{"number", "Number of splines used in expansion [0]"},
       {"order", "order of splines ~7-9 [7]"},
       {"r0", "minimum cavity radius (first internal knot) [1.0e-4]"},
       {"r0_eps", "Select cavity radius r0 for each l by position "
                  "where |psi(r0)/psi_max| falls below r0_eps [1.0e-3]"},
       {"rmax", "maximum cavity radius [Grid{rmax}]"},
       {"states", "states to keep (e.g., 30spdf20ghi)"},
       {"orthogonalise", "Force orthogonal to core [false]"},
       {"print", "Print all spline energies (for testing) [false]"},
       {"positron", "Include -ve energy states [false]]"},
       {"type", "Derevianko (DKB) or Johnson [Derevianko]"}});

  const auto basis_input = input.getBlock("Basis");
  if (basis_input) {
    wf.formBasis(*basis_input);
    if (input.get({"Basis"}, "print", false) && !wf.basis().empty()) {
      std::cout << "Basis:\n";
      wf.printBasis(wf.basis());
    }
  }

  // Correlations: read in options
  // This is a mess - will re-do correlations part
  const auto Sigma_ok = input.check(
      {"Correlations"},
      {{"", "Options for inclusion of correlations (correlation potential "
            "method). It's become a bit of a mess, and will be refactored "
            "~soon~"},
       {"Brueckner", "Form Brueckner orbitals [false]"},
       {"energyShifts", "Calculate MBPT2 shift [false]"},
       {"n_min_core", "Minimum core n to polarise [1]"},
       {"fitTo_cm", "List of binding energies (in cm^-1) to scale Sigma for. "
                    "Must be in same order as valence states"},
       {"lambda_kappa",
        "Scaling factors for Sigma. Must be in same order as valence states"},
       {"read", "Filename to read in Sigma [false=don't read]"},
       {"write", "Filename to write Sigma to [false=don't write]"},
       {"rmin", "minimum radius to calculate sigma for [1.0e-4]"},
       {"rmax", "maximum radius to calculate sigma for [30.0]"},
       {"stride", "Only calculate Sigma every <stride> points"},
       {"each_valence", "Different Sigma for each valence states? [false]"},
       {"ek",
        "Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127, "
        "7s+=-0.552;}. Blank => HF energies"},
       {"Feynman", "Use Feynman method [false]"},
       {"fk",
        "List of doubles. Screening factors for effective all-order "
        "exchange. In Feynman method, used in exchange+ladder only; "
        "Goldstone, "
        "used direct also. If blank, will calculate them from scratch. []"},
       {"eta", "List of doubles. Hole-Particle factors. In Feynman method, "
               "used in ladder only; Goldstone, used direct also. []"},
       {"screening", "bool. Include Screening [false]"},
       {"holeParticle", "Include hole-particle interaction [false]"},
       {"ladder", "Experimental feature. Filename for ladder diagram file "
                  "(generated in the "
                  "ladder Module). If blank, ladder "
                  "not included. Only in Feynman. []"},
       {"lmax", "Maximum l used for Feynman method [6]"},
       {"basis_for_Green", "Use basis for Feynman Greens function [false]"},
       {"basis_for_pol", "Use basis for Feynman polarisation op [false]"},
       {"real_omega", "[worked out by default]"},
       {"imag_omega", "w0, wratio for Im(w) grid [0.01, 1.5]"},
       {"include_G", "Inlcude lower g-part into Sigma [false]"}});

  const bool do_energyShifts =
      input.get({"Correlations"}, "energyShifts", false);
  const bool do_brueckner = input.get({"Correlations"}, "Brueckner", false);
  const auto n_min_core = input.get({"Correlations"}, "n_min_core", 1);
  const auto sigma_rmin = input.get({"Correlations"}, "rmin", 1.0e-4);
  const auto sigma_rmax = input.get({"Correlations"}, "rmax", 30.0);
  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto stride =
        int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (stride <= 2) ? 2 : stride;
  }();
  const auto sigma_stride =
      input.get({"Correlations"}, "stride", default_stride);

  // Feynman method:
  const auto sigma_Feynman = input.get({"Correlations"}, "Feynman", false);
  const auto sigma_Screening = input.get({"Correlations"}, "screening", false);
  const auto hole_particle = input.get({"Correlations"}, "holeParticle", false);
  const auto sigma_lmax = input.get({"Correlations"}, "lmax", 6);
  const auto GreenBasis = input.get({"Correlations"}, "basis_for_Green", false);
  const auto PolBasis = input.get({"Correlations"}, "basis_for_pol", false);
  const auto each_valence = input.get({"Correlations"}, "each_valence", false);
  const auto include_G = input.get({"Correlations"}, "include_G", false);
  const auto ladder_file = input.get({"Correlations"}, "ladder", ""s);
  // force sigma_omre to be always -ve
  const auto sigma_omre = -std::abs(
      input.get({"Correlations"}, "real_omega", -0.33 * wf.energy_gap()));

  // Imaginary omegagrid params (only used for Feynman)
  double w0 = 0.01;
  double wratio = 1.5;
  {
    const auto imag_om =
        input.get({"Correlations"}, "imag_omega", std::vector{w0, wratio});
    if (imag_om.size() != 2) {
      std::cout << "ERROR: imag_omega must be a list of 2: omega_0 (first "
                   "step), and omega_ratio (ratio for log w grid)\n";
    } else {
      w0 = imag_om[0];
      wratio = imag_om[1];
    }
  }

  // Solve Sigma at specific energies
  // e.g., ek{6s+=-0.127, 7s+=-0.552;}
  const auto ek_Sig = input.getBlock({"Correlations"}, "ek");

  // Read/write Sigma to file:
  auto sigma_write = input.get({"Correlations"}, "write", ""s);
  // By default,  try to  read  from  write  file  (if it exists)
  const auto sigma_read = input.get({"Correlations"}, "read", sigma_write);
  // don't  write to default filename when reading from another file
  if (sigma_read != "" && sigma_write == "")
    sigma_write = "false";

  // To fit Sigma to energies:
  // (nb: energies given in cm^-1, convert to au on input)
  auto fit_energies =
      (1.0 / PhysConst::Hartree_invcm) *
      input.get({"Correlations"}, "fitTo_cm", std::vector<double>{});
  const auto lambda_k =
      input.get({"Correlations"}, "lambda_kappa", std::vector<double>{});

  const auto fk = input.get({"Correlations"}, "fk", std::vector<double>{});
  const auto etak = input.get({"Correlations"}, "eta", std::vector<double>{});

  // Form correlation potential:
  if ((do_energyShifts || do_brueckner) && Sigma_ok) {
    IO::ChronoTimer time("Sigma");
    wf.formSigma(n_min_core, do_brueckner, sigma_rmin, sigma_rmax, sigma_stride,
                 each_valence, include_G, lambda_k, fk, etak, sigma_read,
                 sigma_write, ladder_file, sigma_Feynman, sigma_Screening,
                 hole_particle, sigma_lmax, GreenBasis, PolBasis, sigma_omre,
                 w0, wratio, ek_Sig);
  }

  // Calculate + print second-order energy shifts
  if (!wf.valence().empty() && do_energyShifts && Sigma_ok) {
    wf.SOEnergyShift();
  }

  // Solve Brueckner orbitals (optionally, fit Sigma to exp energies)
  if (!wf.valence().empty() && do_brueckner && Sigma_ok) {
    std::cout << '\n';
    IO::ChronoTimer time("Brueckner");
    if (!fit_energies.empty())
      wf.fitSigma_hfBrueckner(valence, fit_energies);
    else
      wf.hartreeFockBrueckner();
  }
  // Print out info for new "Brueckner" valence orbitals:
  if (!wf.valence().empty() && do_brueckner && Sigma_ok) {
    std::cout << "\nBrueckner orbitals:\n";
    wf.printValence();
  }

  // Construct B-spline Spectrum:
  input.check({"Spectrum"},
              {{"", "Options for 'spectrum', Spectrum is the same as 'Basis', "
                    "but includes correlations. Spectrum is used for "
                    "sum-over-states (while basis is used for MBPT)."},
               {"number", "Number of splines used in expansion"},
               {"order", "order of splines ~7-9"},
               {"r0", "minimum cavity radius"},
               {"r0_eps", "Select cavity radius r0 for each l by position "
                          "where |psi(r0)/psi_max| falls below r0_eps"},
               {"rmax", "maximum cavity radius"},
               {"states", "states to keep (e.g., 30spdf20ghi)"},
               {"orthogonalise", "Force orthogonal to valence [false]"},
               {"print", "Print all spline energies (for testing)"},
               {"positron", "Include -ve energy states (true/false)"},
               {"type", "Derevianko (DKB) or Johnson [Derevianko]"}});

  const auto spectrum_in = input.getBlock("Spectrum");
  if (spectrum_in) {
    wf.formSpectrum(*spectrum_in);
    if (input.get({"Spectrum"}, "print", false) && !wf.spectrum().empty()) {
      std::cout << "Spectrum:\n";
      wf.printBasis(wf.spectrum());
    }
  }

  // run each of the modules with the calculated wavefunctions
  Module::runModules(input, wf);
}
