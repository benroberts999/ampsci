#include "ampsci.hpp"
#include "DiracODE/BoundState.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/InputBlock.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp" //for 'ExtraPotential'
#include "Modules/modules_list.hpp"
#include "Modules/runModules.hpp"
#include "Physics/include.hpp"
#include "Wavefunction/Wavefunction.hpp"

Wavefunction ampsci(const IO::InputBlock &input) {
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

  // Check Top-level input blocks
  input.check(
    {{"", "These are the top-level ampsci input blocks and options. Default "
          "values are given in square brackets following the description: "
          "[default_value]. Blocks end with '{}', options end with ';'. run "
          "`ampsci -a BlockName` (for any of the following blocks) to see "
          "all the options available for that block."},
     {"Atom{}", "Which atom to run for"},
     {"Nucleus{}", "Set nuclear parameters"},
     {"Grid{}", "Set radial grid/lattice parameters"},
     {"HartreeFock{}", "Options for solving atomic system"},
     {"RadPot{}",
      "Options for the QED radiative potential (usually defaults suffice)"},
     {"ExtraPotential{}", "Include an extra effective potential. Rarely used."},
     {"Basis{}", "Basis of HF eigenstates used for MBPT"},
     {"Correlations{}", "Options for MBPT and correlation corrections"},
     {"Spectrum{}",
      "Like basis, but includes correlations. Used for sum-over-states"},
     {"Exotic{}", "Option for including `exotic` (e.g., muonic) atom states."},
     {"CI{}", "Configuration Interaction"},
     {"Module::*{}", "Run any number of modules (* -> module name). `ampsci "
                     "-m` to see available modules"}});

  //----------------------------------------------------------------------------
  // Atom: Get + setup atom parameters
  input.check(
    {"Atom"},
    {{"Z", "Atomic number or symbol (e.g., 55 or Cs). [H]"},
     {"A", "Atomic mass number, for nuclear parameters including "
           "finite nuclear size. Default based on Z."},
     {"varAlpha2",
      "Fractional variation of the fine-structure constant, alpha^2: "
      "(a/a0)^2. Use to enforce the non-relativistic limit "
      "(c->infinity => alpha->0), or calculate sensitivity to "
      "variation of alpha. [1.0]"},
     {"run_label", "Optional label for output identity - for distinguishing "
                   "outputs with different parameters"},
     {"json_out", "Write (partial) wavefunction details to json file? "
                  "true/false [false]"}});

  const auto atom_block = input.get_block("Atom");

  // Z is taken as a string, so can write "Cs" or "55"
  const auto atom_Z = AtomData::atomic_Z(atom_block.get("Z", "H"s));
  // A: is a "std::optional<int>" (if blank, nucleus will use default value)
  const auto atom_A = atom_block.get<int>("A");

  // Variation of fine structure constant
  const auto var_alpha = [&]() {
    const auto varAlpha2 = atom_block.get("varAlpha2", 1.0);
    // cannot explicitely set alpha to zero - so make it very small
    return (varAlpha2 > 0.0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();
  const auto run_label = atom_block.get("run_label", ""s);

  //----------------------------------------------------------------------------
  // Nucleus: parse nuclear options
  const auto nucleus_options = input.get_block("Nucleus");

  const auto nucleus = Nuclear::form_nucleus(atom_Z, atom_A, nucleus_options);

  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  // Create wavefunction object
  Wavefunction wf(radial_grid, std::move(nucleus), var_alpha, run_label);

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

  //----------------------------------------------------------------------------
  // Parse input for Hartree-Fock
  input.check(
    {"HartreeFock"},
    {{"", "Options for solving lowest-order atomic wavefunction"},
     {"core", "Core configuration. Either list entire core, or use [At] "
              "short-hand. e.g., [He] equivilant to 1s2; [Xe],6s1 equivilant "
              "to [Cs] and to 1s2,2s2,...,5p6,6s1. Instead of one of the "
              "commas, you may use a ':' - states above this are included "
              "into the core, but not excluded from the valence list. Use "
              "this method for KohnSham, for example. [blank by default]"},
     {"valence",
      "Valence configuration in `basis string' format. e.g., 7sp5df will "
      "include valence states up to  n=7 for s and p, and up to n=5 for d "
      "and f states. Automatically excludes states in the core (except "
      "those above the optional ':'). [blank by default]"},
     {"eps", "HF convergance goal [1.0e-13]"},
     {"method", "Method for mean-field approximation: HartreeFock, Hartree, "
                "KohnSham, Local [HartreeFock]"},
     {"Breit", "Include Breit into HF? true/false, or scale factor. Scale "
               "factor for Breit Hamiltonian is usially 0.0 (no "
               "Breit) or 1.0 (full Breit), but can take any value. [false]"},
     {"QED",
      "Include QED? Three options: true, false, valence. If 'valence, will "
      "include QED only into valence states, but not the core. Detailed QED "
      "options are set within the RadPot{} block - if that block is not set, "
      "defaults will be used. By default, this option is false, unless the "
      "RadPot{} block exists, in which case it is true"}});

  const auto core = input.get({"HartreeFock"}, "core", "[]"s);
  const auto HF_method = input.get({"HartreeFock"}, "method", "HartreeFock"s);
  const auto eps_HF = input.get({"HartreeFock"}, "eps", 1.0e-13);
  const auto tf_Breit = input.get({"HartreeFock"}, "Breit", false);
  const auto x_Breit =
    tf_Breit ? 1.0 : input.get({"HartreeFock"}, "Breit", 0.0);
  const auto valence = input.get({"HartreeFock"}, "valence", ""s);

  // Decide if to include QED into core+valence, just core, or not at all
  const auto qed_input = input.getBlock("RadPot");
  const auto [include_qed, qed_core] = [&]() {
    const auto default_qed = qed_input ? "true"s : "false"s;
    const auto qed_option = input.get({"HartreeFock"}, "QED", default_qed);
    const auto valence_qed = qip::ci_compare(qed_option, "valence");
    const auto t_include_qed =
      qip::ci_compare(qed_option, "true") || valence_qed;
    return std::pair{t_include_qed, !valence_qed};
  }();

  // Set up the Hartree Fock potential/method (does not solve)
  // (Must set HF before adding RadPot - but must add RadPot before solving HF)
  wf.set_HF(HF_method, x_Breit, core, eps_HF, true);

  // Forms QED radiative potential, if RadPot{} block is present.
  // Note: input options are parsed inside radiativePotential()
  if (include_qed && qed_core) {
    std::cout << "\nIncluding QED into core (and valence)\n";
    wf.radiativePotential(qed_input ? *qed_input : IO::InputBlock{}, true,
                          true);
    std::cout << "\n";
  }

  //----------------------------------------------------------------------------
  // Inlcude extra potential. Either read in from text file.
  // or "polarisation operator": V(r) = -0.5/(r^4 + r_cut^4)
  // Note: If read in, it is interpolated onto grid, but NOT extrapolated
  input.check(
    {"ExtraPotential"},
    {{"", "Adds an extra potential (to Vnuc), before HF solved. Either "
          "effective polarisation potential, V(r) = -0.5/(r^4 + r_cut^4), or "
          "read in from a file."},
     {"filename", "Read potential from file (r v(r)) - will be interpolated. "
                  "If not given, will use pol. potential"},
     {"r_cut", "Radial cut-off parameter for effective pol. potential [=1]"},
     {"scale", "Overall scaling factor for potential is scaled by this "
               "value [1]"}});

  const auto extra_in = input.getBlock("ExtraPotential");
  if (extra_in && !extra_in->has_option("help")) {
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
      std::cout << "Adding effective polarisation potential "
                   "[-0.5 * a * (r^4 + rc^4)]: a="
                << ep_scale << ", rc=" << rc << '\n';
      const auto a4 = rc * rc * rc * rc;
      auto dV = [=](auto r) { return -0.5 / (r * r * r * r + a4); };
      const auto dv = qip::apply_to(dV, wf.grid().r());
      wf.update_Vnuc(wf.vnuc() + ep_scale * dv);
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // Solve Hartree equations for the core:
  wf.solve_core(true);

  if (include_qed && !qed_core) {
    std::cout << "\nIncluding QED into valence (not included in core)\n";
    wf.radiativePotential(qed_input ? *qed_input : IO::InputBlock{}, true,
                          true);
    std::cout << "\n";
  }

  // Solve for the valence states:
  wf.solve_valence(valence);

  // Output Hartree Fock energies:
  if (!help_mode) {
    std::cout << '\n' << wf.atomicSymbol() << "-" << wf.Anuc() << '\n';
    wf.printCore();
    printf("E_c = %.6f\n\n", wf.coreEnergyHF());

    wf.printValence();
  }

  //----------------------------------------------------------------------------

  // Construct B-spline basis:
  input.check(
    {"Basis"},
    {{"number", "Number of splines used in expansion [30]"},
     {"order", "order of splines ~7-9 [7]"},
     {"r0", "minimum cavity radius (first internal knot) [1.0e-4]"},
     {"r0_eps", "Select cavity radius r0 for each l by position "
                "where |psi(r0)/psi_max| falls below r0_eps [0.0]"},
     {"rmax", "maximum cavity radius [40.0]"},
     {"states", "states to keep (e.g., 30spdf20ghi)"},
     {"orthogonalise", "Force orthogonal to core [false]"},
     {"print", "Print all spline energies (for testing) [false]"},
     {"positron",
      "Basis string for negative energy states (same format as states). []"},
     {"type", "Derevianko (DKB) or Johnson [Derevianko]"}});

  const auto basis_input = input.getBlock("Basis");
  if (basis_input) {
    wf.formBasis(*basis_input);
    if (input.get({"Basis"}, "print", false) && !wf.basis().empty()) {
      std::cout << "Basis:\n";
      wf.printBasis(wf.basis());
    }
  }

  //----------------------------------------------------------------------------

  // Correlations: read in options
  // This is a mess - will re-do correlations part
  const auto Sigma_ok = input.check(
    {"Correlations"},
    {{"", "Options for inclusion of correlations (correlation potential "
          "method)."},
     {"n_min_core", "Minimum core n to polarise [1]"},
     //  {"n_min_core_F",
     //   "Minimum core n to polarise in Feynman method. By default, same as "
     //   "n_min_core. If n_min_core_F>n_min_core, code will use Goldstone "
     //   "method for lowest ns"},
     {"each_valence",
      "Construct seperate Sigma for each valence state? [false]"},
     {"fitTo_cm", "List of binding energies (in cm^-1) to scale Sigma for. "
                  "Must be in same order as valence states"},
     {"lambda_kappa",
      "Scaling factors for Sigma. Must be in same order as valence states"},
     {"read", "Read/write correlation potential to disk [true]"},
     {"filename",
      "Filename prefix (i.e., not including extension) to "
      "read/write correlation matrix to/from. "
      "Normally, this option should remain unset. By default it will be "
      "<identity>.sigX.abf, where X details parameters of Sigma."},
     {"rmin", "minimum radius to calculate sigma for [1.0e-4]"},
     {"rmax", "maximum radius to calculate sigma for [30.0]"},
     {"stride", "Only calculate Sigma every <stride> points. Default such "
                "that there are 150 points between (1e-4, 30)"},
     {"ek{}", "Block: Explicit list of energies to solve for. e.g., "
              "ek{6s+=-0.127; 7s+=-0.552;}. Blank => HF energies. Takes "
              "precidence over each_valence. [blank]"},
     {"all_order", "Use all-orders method (implies Feynman=true; "
                   "screening=true; hole_particle=true;) [false]"},
     {"Feynman", "Use Feynman method [false]"},
     {"fk", "List of doubles. Screening factors for effective all-order "
            "exchange. In Feynman method, used in exchange only (and G-part); "
            "Goldstone, used direct also. If blank, will calculate them from "
            "scratch. []"},
     {"eta", "List of doubles. Hole-Particle factors. In Feynman method, "
             "used only for G part; Goldstone, used in direct also. []"},
     {"screening", "Include all-orders screening. Only applicable for "
                   "Feynman method [false]"},
     {"hole_particle",
      "Include all-orders hole-particle interaction. Only applicable for "
      "Feynman method [false]"},
     {"lmax", "Maximum l used for internal lines in Feynman method [6]"},
     {"real_omega", "Real part of frequency used in contour integral. By "
                    "Default, ~1/3 of the core/valence energy gap"},
     {"imag_omega",
      "Pair of comma-separated doubles: w0, wratio. Initial point, and "
      "ratio, for logarithimg Im(w) grid [0.01, 1.5]"},
     {"include_G", "Inlcude lower g-part into Sigma [false]"},
     {"include_Breit", "Inlcude two-body Breit corrections into Sigma [false]"},
     {"n_max_Breit",
      "Maximum n for excited states to include in two-body Breit "
      "correction to Correlation potential [<=0, means entire basis]"}});

  const bool do_brueckner = input.getBlock({"Correlations"}) != std::nullopt;
  const auto n_min_core = input.get({"Correlations"}, "n_min_core", 1);
  // const auto n_min_core_F =
  //     input.get({"Correlations"}, "n_min_core_F", n_min_core);
  const auto sigma_rmin = input.get({"Correlations"}, "rmin", 1.0e-4);
  const auto sigma_rmax = input.get({"Correlations"}, "rmax", 30.0);
  const auto each_valence = input.get({"Correlations"}, "each_valence", false);
  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto stride =
      int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;
    return (stride <= 2) ? 2 : stride;
  }();
  const auto sigma_stride =
    input.get({"Correlations"}, "stride", default_stride);

  // Feynman method:
  const auto all_order = input.get({"Correlations"}, "all_order", false);
  const auto sigma_Feynman = input.get({"Correlations"}, "Feynman", all_order);
  const auto sigma_Screening =
    input.get({"Correlations"}, "screening", all_order);
  const auto hole_particle =
    input.get({"Correlations"}, "hole_particle", all_order);
  const auto sigma_lmax = input.get({"Correlations"}, "lmax", 6);
  const auto include_G = input.get({"Correlations"}, "include_G", false);
  const auto include_Breit =
    input.get({"Correlations"}, "include_Breit", false);
  const auto n_max_Breit = input.get({"Correlations"}, "n_max_Breit", -1);
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

  const auto sigma_readwrite = input.get({"Correlations"}, "read", true);
  const auto sigma_filename = input.get({"Correlations"}, "filename", ""s);

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
  if (Sigma_ok && do_brueckner) {
    IO::ChronoTimer time("Sigma");
    wf.formSigma(n_min_core, sigma_rmin, sigma_rmax, sigma_stride, each_valence,
                 include_G, include_Breit, n_max_Breit, lambda_k, fk, etak,
                 sigma_readwrite, sigma_filename, sigma_Feynman,
                 sigma_Screening, hole_particle, sigma_lmax, sigma_omre, w0,
                 wratio, ek_Sig);
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

  //----------------------------------------------------------------------------

  // Construct B-spline Spectrum:
  input.check(
    {"Spectrum"},
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
     {"positron", "Basis string (as states) to include -ve energy states"},
     {"type", "Derevianko (DKB) or Johnson [Derevianko]"}});

  const auto spectrum_in = input.getBlock("Spectrum");
  if (spectrum_in && !spectrum_in->has_option("help")) {
    wf.formSpectrum(*spectrum_in);
    if (input.get({"Spectrum"}, "print", false) && !wf.spectrum().empty()) {
      std::cout << "Spectrum:\n";
      wf.printBasis(wf.spectrum());
    }
  }

  //----------------------------------------------------------------------------
  const auto exotic = input.getBlock("Exotic");
  input.check(
    {"Exotic"},
    {{"", "Option for including `exotic` (e.g., muonic) atom states.\n"
          "Adds them to end of valence list; usually valence should be empty.\n"
          "Includes screening. Use muon module as test"},
     {"mass", "Mass (in au=m_e) of exotic lepton [M_muon = 206.7682827]"},
     {"mass_MeV", "Mass (in MeV) of exotic lepton [M_muon = 105.6583755 "
                  "MeV]; will be overridden by the above au version"},
     {"states", "Which states to calculate [1s]"}});
  if (exotic && !exotic->has_option("help")) {

    const auto states_str = exotic->get("states", std::string{"1s"});

    const auto mass_MeV = exotic->get<double>("mass_MeV");
    const auto mass = exotic->get(
      "mass", mass_MeV ? *mass_MeV / PhysConst::m_e_MeV : PhysConst::m_muon);

    wf.solve_exotic(states_str, mass, true);
  }

  //----------------------------------------------------------------------------
  const auto CI_in = input.getBlock("CI");
  if (CI_in) {
    wf.ConfigurationInteraction(*CI_in);
  }

  const auto json_out = input.get({"Atom"}, "json_out", false);
  if (json_out) {
    const std::string json_out_name = wf.identity() + ".json";
    wf.output_to_json(json_out_name);
  }

  // run each of the modules with the calculated wavefunctions
  Module::runModules(input, wf);

  return wf;
}
