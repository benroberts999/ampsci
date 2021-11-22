#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/InputBlock.hpp"
#include "Maths/Interpolator.hpp" //for 'ExtraPotential'
#include "Modules/runModules.hpp"
#include "Physics/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "git.info"
#include "qip/Vector.hpp"
#include <iostream>
#include <string>

void ampsci(const IO::InputBlock &input);
int main(int argc, char *argv[]) {
  std::cout << "\n";
  IO::print_line();
  std::cout << "ampsci git:" << GitInfo::gitversion << " ("
            << GitInfo::gitbranch << ")\n";
  std::cout << IO::time_date() << '\n';

  const std::string input_text = (argc > 1) ? argv[1] : "ampsci.in";

  // std::filesystem not available in g++-7 (getafix version)
  // Reading from a file? Or from command-line?
  const auto fstream = std::fstream(input_text);
  const auto symb = AtomData::atomicSymbol(AtomData::atomic_Z(input_text));
  const auto core = symb == "H" ? "" : symb;
  const std::string default_input = (input_text.size() <= 3) ?
                                        "Atom{Z=" + symb + ";}" +
                                            "HartreeFock { core = [" + core +
                                            "]; valence = 2sp;}" :
                                        input_text;

  const auto input = fstream.good() ? IO::InputBlock("ampsci", fstream) :
                                      IO::InputBlock("ampsci", default_input);

  // Run program. Add option to run multiple times
  ampsci(input);

  return 0;
}

//******************************************************************************
void ampsci(const IO::InputBlock &input) {
  using namespace std::string_literals;
  IO::ChronoTimer timer("\nampsci");
  std::cout << "\n";
  IO::print_line();
  input.print();

  input.check({{"Atom", "Which atom to run for"},
               {"Grid", "Set radial grid parameters"},
               {"HartreeFock", "Expl"},
               {"Nucleus", "Set nuclear parameters"},
               {"RadPot", "Inlcude QED radiative potential"},
               {"Basis", "Basis used for MBPT"},
               {"Spectrum", "Like basis; used for sum-over-states"},
               {"Correlations", "Options for correlations"},
               {"ExtraPotential", "Include an extra potential"},
               {"dVpol", "Approximate correlation (polarisation) potential"},
               {"Module::*", "Run any number of modules (* -> module name)"}});

  // Atom: Get + setup atom parameters
  auto input_ok =
      input.check({"Atom"}, {{"Z", "Atomic symbol/number (int or string)"},
                             {"A", "Atomic mass number (blank for default)"},
                             {"varAlpha2", "d(a^2)/a_0^2 (1 by default)"}});

  const auto atom_Z = AtomData::atomic_Z(input.get({"Atom"}, "Z", "H"s));
  const auto atom_A = input.get({"Atom"}, "A", AtomData::defaultA(atom_Z));
  const auto var_alpha = [&]() {
    const auto varAlpha2 = input.get({"Atom"}, "varAlpha2", 1.0);
    return (varAlpha2 > 0.0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Grid: Get + setup grid parameters
  input_ok &= input.check(
      {"Grid"},
      {{"r0", "Initial grid point, in au (~1e-6)"},
       {"rmax", "Finial grid point ~100.0"},
       {"num_points", "Number of grid points ~1000"},
       {"type", "loglinear or logarithmic"},
       {"b", "only loglinear: roughly logarithmic for r<b, linear for r>b"},
       {"du", "du is uniform grid step size; set this instead of "
              "num_points (~0.1)"}});

  const auto r0 = input.get({"Grid"}, "r0", 1.0e-6);
  const auto rmax = input.get({"Grid"}, "rmax", 120.0);
  // du>0 means calc num_points
  const auto du_option = input.get<double>({"Grid"}, "du");
  const auto du = du_option ? *du_option : -1.0;
  const auto num_points =
      (du > 0.0) ? 0ul : input.get({"Grid"}, "num_points", 1600ul);
  const auto b = input.get({"Grid"}, "b", 0.33 * rmax);
  const auto grid_type =
      (b <= r0 || b >= rmax) ?
          "logarithmic" :
          input.get<std::string>({"Grid"}, "type", "loglinear");

  // Nucleus: Get + setup nuclear parameters
  input_ok &= input.check({"Nucleus"},
                          {{"rrms", "root-mean-square charge radius, in fm "
                                    "(blank means will look up default)"},
                           {"c", "Half-density radius (use instead of rrms)"},
                           {"t", "Nuclear skin thickness; default = 2.3"},
                           {"type", "Fermi, spherical, pointlike, Gaussian"}});

  // usually, give rrms. Giving c will over-ride rrms
  const auto c_hdr = input.get<double>({"Nucleus"}, "c");
  auto dflt_rrms = Nuclear::find_rrms(atom_Z, atom_A);
  if (dflt_rrms <= 0.0 && atom_A != 0 && !c_hdr) {
    dflt_rrms = Nuclear::approximate_r_rms(atom_A);
    std::cout << "\nWARNING: isotope Z=" << atom_Z << ", A=" << atom_A
              << " - cannot find rrms. Using approx formula: rrms=" << dflt_rrms
              << "\n";
  }
  const auto t_skin = input.get({"Nucleus"}, "t", Nuclear::default_t);
  // c (half density radius) takes precidence if c and r_rms are given.
  const auto rrms = c_hdr ? Nuclear::rrms_formula_c_t(*c_hdr, t_skin) :
                            input.get({"Nucleus"}, "rrms", dflt_rrms);

  // Set nuc. type explicitely to 'pointlike' if A=0, or r_rms = 0.0
  const auto nuc_type =
      (atom_A == 0 || rrms == 0.0) ?
          "pointlike" :
          input.get<std::string>({"Nucleus"}, "type", "Fermi");

  // Create wavefunction object
  Wavefunction wf({num_points, r0, rmax, b, grid_type, du},
                  {atom_Z, atom_A, nuc_type, rrms, t_skin}, var_alpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid->gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  input_ok &= input.check(
      {"HartreeFock"},
      {{"core", "Core configuration. e.g., [Xe] for Cs"},
       {"valence", "Which valence states? e.g., 7sp5d"},
       {"convergence", "HF convergance goal, 1e-12"},
       {"method", "HartreeFock(default), Hartree, KohnSham"},
       {"Breit", "Scale for Breit. 0.0 default (no Breit), 1.0 include Breit"},
       {"sortOutput", "Sort energy tables by energy? (default=false)"}});

  if (!input_ok) {
    std::cout
        << "\nProgram halted due to input errors; review above warnings\n";
    return;
  }

  const auto str_core = input.get<std::string>({"HartreeFock"}, "core", "[]");
  const auto eps_HF = input.get({"HartreeFock"}, "convergence", 1.0e-12);
  const auto HF_method =
      input.get<std::string>({"HartreeFock"}, "method", "HartreeFock");
  if (HF_method == "Hartree")
    std::cout << "Using Hartree Method (no Exchange)\n";
  else if (HF_method == "ApproxHF")
    std::cout << "Using approximate HF Method (approx Exchange)\n";
  else if (HF_method == "KohnSham") {
    std::cout << "Using Kohn-Sham Method.\n"
              << "Note: You should include first valence state into the core:\n"
                 "Kohn-Sham is NOT a V^N-1 method!\n";
  } else if (HF_method == "Local") {
    std::cout << "Using local potential\n";
  } else if (HF_method != "HartreeFock") {
    std::cout << "\n⚠️  WARNING unkown method: " << HF_method
              << "\nDefaulting to HartreeFock method.\n";
  }

  // Breit:
  const auto x_Breit = input.get({"HartreeFock"}, "Breit", 0.0);
  // Can only include Breit within HF
  if (HF_method == "HartreeFock" && x_Breit != 0.0) {
    std::cout << "Including Breit (scale = " << x_Breit << ")\n";
  } else if (HF_method != "HartreeFock" && x_Breit != 0.0) {
    std::cout << "\n⚠️  WARNING can only include Breit in Hartree-Fock "
                 "method. Breit will not be included.\n";
  }

  // Inlcude QED radiatve potential
  const auto qed_ok = input.check(
      {"RadPot"},
      {{"",
        "(QED Radiative potential will be included if this block is present)"},
       {"Ueh", " for Uehling typical [0.0, 1.0], Default = 1"},
       {"SE_h", " for self-energy high-freq electric. Default = 1"},
       {"SE_l", " for self-energy low-freq electric. Default = 1"},
       {"SE_m", " self-energy magnetic. Default = 1"},
       {"WK", " Wickman-Kroll. Default = 0"},
       {"rcut", "Maximum r to calculate Rad Pot (~5)"},
       {"scale_rN", "Nuclear size. 0 for pointlike, 1 for typical"},
       {"scale_l", "Extra scaling factor for each l e.g., (1,1,1)"},
       {"core_qed", "Include rad pot into core Hartree-Fock (default=true)"}});

  const auto include_qed = input.getBlock("RadPot") != std::nullopt;
  const auto x_Ueh = input.get({"RadPot"}, "Ueh", 1.0);
  const auto x_SEe_h = input.get({"RadPot"}, "SE_h", 1.0);
  const auto x_SEe_l = input.get({"RadPot"}, "SE_l", 1.0);
  const auto x_SEm = input.get({"RadPot"}, "SE_m", 1.0);
  const auto x_wk = input.get({"RadPot"}, "WK", 0.0);
  const auto rcut = input.get({"RadPot"}, "rcut", 5.0);
  const auto scale_rN = input.get({"RadPot"}, "scale_rN", 1.0);
  const auto x_spd = input.get({"RadPot"}, "scale_l", std::vector{1.0});
  const bool core_qed = input.get({"RadPot"}, "core_qed", true);

  if (include_qed && qed_ok && core_qed) {
    wf.radiativePotential({x_Ueh, x_SEe_h, x_SEe_l, x_SEm, x_wk}, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Hartree-Fock core (and valence)\n\n";
  }

  // Inlcude extra potential (read in from text file):
  // Note: interpolated onto grid, but NOT extrapolated
  // (zero outside region!)
  const auto extra_ok = input.check(
      {"ExtraPotential"},
      {{"filename", ""},
       {"factor", "potential is scaled by this value [default=1]"},
       {"beforeHF", "include before HF (into core states). default=false"}});
  const auto ep_fname =
      input.get<std::string>({"ExtraPotential"}, "filename", "");
  const auto ep_factor = input.get({"ExtraPotential"}, "factor", 1.0);
  const auto ep_beforeHF = input.get({"ExtraPotential"}, "beforeHF", false);
  const auto extra_pot =
      ep_fname != "" && std::abs(ep_factor) > 0.0 && extra_ok;
  std::vector<double> Vextra;
  if (extra_pot) {
    const auto &[x, y] = IO::FRW::readFile_xy_PoV("testIn.txt");
    Vextra = Interpolator::interpolate(x, y, wf.rgrid->r());
    qip::scale(&Vextra, ep_factor);
  }

  // Add "extra potential", before HF (core + valence)
  if (extra_pot && ep_beforeHF) {
    qip::add(&wf.vnuc, Vextra);
  }

  { // Solve Hartree equations for the core:
    IO::ChronoTimer t(" core");
    wf.solve_core(HF_method, x_Breit, str_core, eps_HF);
  }

  if (include_qed && qed_ok && !core_qed) {
    wf.radiativePotential({x_Ueh, x_SEe_h, x_SEe_l, x_SEm, x_wk}, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Valence only\n\n";
  }

  // Add "extra potential", after HF (only valence)
  if (extra_pot && !ep_beforeHF) {
    wf.add_to_Vdir(Vextra);
  }

  // Adds effective polarision potential to direct
  // potential (After HF core, before HF valence)
  const auto Vpol_ok = input.check(
      {"dVpol"}, {{"a_eff", "scale factor for effective pol. potential [1]"},
                  {"r_cut", "cut-off parameter [=1]"}});
  const auto a_eff = input.get({"dVpol"}, "a_eff", 0.0);
  if (std::abs(a_eff) > 0.0 && Vpol_ok) {
    const auto r_cut = input.get({"dVpol"}, "r_cut", 1.0);
    const auto a4 = r_cut * r_cut * r_cut * r_cut;
    auto dV = [=](auto x) { return -0.5 * a_eff / (x * x * x * x + a4); };
    std::vector<double> dv(wf.rgrid->num_points());
    for (auto i = 0u; i < wf.rgrid->num_points(); ++i) {
      dv[i] = dV(wf.rgrid->r()[i]);
    }
    wf.add_to_Vdir(dv);
  }

  // Solve for the valence states:
  const auto valence_list =
      (wf.Ncore() < wf.Znuc() || HF_method == "KohnSham") ?
          input.get<std::string>({"HartreeFock"}, "valence", "") :
          "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad
    // happens if below are called
    IO::ChronoTimer t("  val");
    wf.solve_valence(valence_list);
  }

  // Output Hartree Fock energies:
  std::cout << "\nHartree Fock: " << wf.identity() << "-" << wf.Anuc() << "\n";
  const auto sorted = input.get({"HartreeFock"}, "sortOutput", false);
  wf.printCore(sorted);
  wf.printValence(sorted);

  // Construct B-spline basis:
  const auto basis_ok = input.check(
      {"Basis"}, {{"number", "Number of splines used in expansion"},
                  {"order", "order of splines ~7-9"},
                  {"r0", "minimum cavity radius (first internal knot)"},
                  {"r0_eps", "Select cavity radius r0 for each l by position "
                             "where |psi(r0)/psi_max| falls below r0_eps"},
                  {"rmax", "maximum cavity radius"},
                  {"states", "states to keep (e.g., 30spdf20ghi)"},
                  {"print", "Print all spline energies (for testing)"},
                  {"positron", "Include -ve energy states (true/false)"},
                  {"type", "Derevianko (DKB) or Johnson"}});
  if (basis_ok) {
    const auto basis_in = input.getBlock("Basis");
    if (basis_in)
      wf.formBasis(*basis_in);
    if (input.get({"Basis"}, "print", false) && !wf.basis.empty()) {
      std::cout << "Basis:\n";
      wf.printBasis(wf.basis);
    }
  }

  // Correlations: read in options
  const auto Sigma_ok = input.check(
      {"Correlations"},
      {{"Brueckner", "Form Brueckner orbitals [false]"},
       {"energyShifts", "Calculate MBPT2 shift [false]"},
       {"n_min_core", "Minimum core n to polarise [1]"},
       {"fitTo_cm", "List of binding energies (in cm^-1) to scale Sigma for. "
                    "Must be in same order as valence states"},
       {"lambda_kappa",
        "Scaling factors for Sigma. Must be in same order as valence states"},
       {"fk", "Screening factors for effective all-order exchange"},
       {"read", "Filename to read in Sigma [false=don't read]"},
       {"write", "Filename to write Sigma to [false=don't write]"},
       {"rmin", "minimum radius to calculate sigma for [1.0e-4]"},
       {"rmax", "maximum radius to calculate sigma for [1.0e-4]"},
       {"stride", "Only calculate Sigma every <stride> points"},
       {"each_valence", "Different Sigma for each valence states? [false]"},
       {"ek",
        "Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127, "
        "7s+=-0.552;}. Blank => HF energies"},
       {"Feynman", "Use Feynman method [false]"},
       {"screening", "Include Screening [false]"},
       {"holeParticle", "Include hole-particle interaction [false]"},
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
        int(wf.rgrid->getIndex(30.0) - wf.rgrid->getIndex(1.0e-4)) / 150;
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
  auto sigma_write = input.get<std::string>({"Correlations"}, "write", "");
  // By default,  try to  read  from  write  file  (if it exists)
  const auto sigma_read = input.get({"Correlations"}, "read", sigma_write);
  // don't  write to default filename when reading from another file
  if (sigma_read != "" && sigma_write == "")
    sigma_write = "false";

  // To fit Sigma to energies:
  auto fit_energies =
      input.get({"Correlations"}, "fitTo_cm", std::vector<double>{});
  // energies given in cm^-1, convert to au:
  qip::scale(&fit_energies, 1.0 / PhysConst::Hartree_invcm);
  const auto lambda_k =
      input.get({"Correlations"}, "lambda_kappa", std::vector<double>{});
  const auto fk = input.get({"Correlations"}, "fk", std::vector<double>{});

  // Form correlation potential:
  if ((do_energyShifts || do_brueckner) && Sigma_ok) {
    IO::ChronoTimer t("Sigma");
    wf.formSigma(n_min_core, do_brueckner, sigma_rmin, sigma_rmax, sigma_stride,
                 each_valence, include_G, lambda_k, fk, sigma_read, sigma_write,
                 sigma_Feynman, sigma_Screening, hole_particle, sigma_lmax,
                 GreenBasis, PolBasis, sigma_omre, w0, wratio, ek_Sig);
  }

  // Calculate + print second-order energy shifts
  if (!wf.valence.empty() && do_energyShifts && Sigma_ok) {
    IO::ChronoTimer t("de");
    wf.SOEnergyShift();
  }

  // Solve Brueckner orbitals (optionally, fit Sigma to exp energies)
  if (!wf.valence.empty() && do_brueckner && Sigma_ok) {
    std::cout << "\n";
    IO::ChronoTimer t("Br");
    if (!fit_energies.empty())
      wf.fitSigma_hfBrueckner(valence_list, fit_energies);
    else
      wf.hartreeFockBrueckner();
  }
  // Print out info for new "Brueckner" valence orbitals:
  if (!wf.valence.empty() && do_brueckner && Sigma_ok) {
    std::cout << "\nBrueckner orbitals:\n";
    wf.printValence(sorted);
  }

  // Construct B-spline Spectrum:
  const auto spectrum_ok =
      input.check({"Spectrum"},
                  {{"number", "Number of splines used in expansion"},
                   {"order", "order of splines ~7-9"},
                   {"r0", "minimum cavity radius"},
                   {"r0_eps", "Select cavity radius r0 for each l by position "
                              "where |psi(r0)/psi_max| falls below r0_eps"},
                   {"rmax", "maximum cavity radius"},
                   {"states", "states to keep (e.g., 30spdf20ghi)"},
                   {"print", "Print all spline energies (for testing)"},
                   {"positron", "Include -ve energy states (true/false)"},
                   {"type", "Derevianko (DKB) or Johnson"}});

  const auto spectrum_in = input.getBlock("Spectrum");
  if (spectrum_ok) {
    if (spectrum_in)
      wf.formSpectrum(*spectrum_in);
    if (input.get({"Spectrum"}, "print", false) && !wf.spectrum.empty()) {
      std::cout << "Spectrum:\n";
      wf.printBasis(wf.spectrum);
    }
  }

  // run each of the modules with the calculated wavefunctions
  Module::runModules(input, wf);
}

//******************************************************************************
