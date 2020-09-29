#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/UserInput.hpp"
#include "Maths/Interpolator.hpp" //for 'ExtraPotential'
#include "Modules/runModules.hpp"
#include "Physics/PhysConst_constants.hpp" //for fit_energies
#include "Wavefunction/Wavefunction.hpp"
#include "git.info"
#include "qip/Vector.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  IO::ChronoTimer timer("\ndiracSCAS");
  const std::string input_file = (argc > 1) ? argv[1] : "diracSCAS.in";
  IO::print_line();

  // Read in input options file
  std::cout << "Reading input from: " << input_file << "\n";
  const IO::UserInput input(input_file);
  std::cout << "diracSCAS git:" << GitInfo::gitversion << " ("
            << GitInfo::gitbranch << ")\n";
  input.print();

  // Atom: Get + setup atom parameters
  auto input_ok = input.check("Atom", {"Z", "A", "varAlpha2"});
  const auto atom_Z = input.get<std::string>("Atom", "Z");
  const auto atom_A = input.get("Atom", "A", -1);
  const auto var_alpha = [&]() {
    const auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
    return (varAlpha2 > 0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Grid: Get + setup grid parameters
  input_ok = input_ok && input.check("Grid", {"r0", "rmax", "num_points",
                                              "type", "b", "fixed_du"});
  const auto r0 = input.get("Grid", "r0", 1.0e-6);
  const auto rmax = input.get("Grid", "rmax", 120.0);
  // du_tmp>0 means calc num_points
  const auto du_tmp = input.get("Grid", "fixed_du", -1.0);
  const auto num_points =
      (du_tmp > 0) ? 0ul : input.get("Grid", "num_points", 1600ul);
  const auto b = input.get("Grid", "b", 0.33 * rmax);
  const auto grid_type =
      (b <= r0 || b >= rmax)
          ? "logarithmic"
          : input.get<std::string>("Grid", "type", "loglinear");

  // Nucleus: Get + setup nuclear parameters
  input_ok = input_ok && input.check("Nucleus", {"rrms", "c", "t", "type"});
  // if {rrms, t} < 0 means get default (depends on A)
  const auto c_hdr = input.get("Nucleus", "c", -1.0);
  const auto skint = input.get("Nucleus", "skin_t", Nuclear::default_t);
  // c (half density radius) takes precidence if c and r_rms are given.
  const auto rrms = c_hdr <= 0.0 ? input.get("Nucleus", "rrms", -1.0)
                                 : Nuclear::rrms_formula_c_t(c_hdr, skint);
  // Set nuc. type explicitely to 'pointlike' if A=0, or r_rms = 0.0
  const auto nuc_type =
      (atom_A == 0 || rrms == 0.0)
          ? "pointlike"
          : input.get<std::string>("Nucleus", "type", "Fermi");

  // Create wavefunction object
  Wavefunction wf({num_points, r0, rmax, b, grid_type, du_tmp},
                  {atom_Z, atom_A, nuc_type, rrms, skint}, var_alpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid->gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  input_ok =
      input_ok && input.check("HartreeFock", {"core", "valence", "convergence",
                                              "method", "Breit", "sortOutput"});
  if (!input_ok) {
    std::cout
        << "\nProgram halted due to input errors; review above warnings\n";
    return 1;
  }

  const auto str_core = input.get<std::string>("HartreeFock", "core", "[]");
  const auto eps_HF = input.get("HartreeFock", "convergence", 1.0e-12);
  const auto HF_method =
      input.get<std::string>("HartreeFock", "method", "HartreeFock");
  if (HF_method == "Hartree")
    std::cout << "Using Hartree Method (no Exchange)\n";
  else if (HF_method == "ApproxHF")
    std::cout << "Using approximate HF Method (approx Exchange)\n";
  else if (HF_method == "KohnSham") {
    std::cout << "Using Kohn-Sham Method.\n"
              << "Note: You should include first valence state into the core:\n"
                 "Kohn-Sham is NOT a V^N-1 method!\n";
  } else if (HF_method != "HartreeFock") {
    std::cout << "\n⚠️  WARNING unkown method: " << HF_method
              << "\nDefaulting to HartreeFock method.\n";
  }

  // Breit:
  const auto x_Breit = input.get("HartreeFock", "Breit", 0.0);
  // Can only include Breit within HF
  if (HF_method == "HartreeFock" && x_Breit != 0.0) {
    std::cout << "Including Breit (scale = " << x_Breit << ")\n";
  } else if (HF_method != "HartreeFock" && x_Breit != 0.0) {
    std::cout << "\n⚠️  WARNING can only include Breit in Hartree-Fock "
                 "method. Breit will not be included.\n";
  }

  // Inlcude QED radiatve potential
  const auto qed_ok =
      input.check("RadPot", {"RadPot", "Simple", "Ueh", "SE_h", "SE_l", "SE_m",
                             "rcut", "scale_rN", "scale_l", "core_qed"});
  const auto include_qed = input.get("RadPot", "RadPot", false);
  const auto x_Simple = input.get("RadPot", "Simple", 0.0);
  const auto xrp_dflt = (include_qed && x_Simple == 0.0) ? 1.0 : 0.0;
  const auto x_Ueh = input.get("RadPot", "Ueh", xrp_dflt);
  const auto x_SEe_h = input.get("RadPot", "SE_h", xrp_dflt);
  const auto x_SEe_l = input.get("RadPot", "SE_l", xrp_dflt);
  const auto x_SEm = input.get("RadPot", "SE_m", xrp_dflt);
  const auto rcut = input.get("RadPot", "rcut", 5.0);
  const auto scale_rN = input.get("RadPot", "scale_rN", 1.0);
  const auto x_spd = input.get_list("RadPot", "scale_l", std::vector{1.0});
  const bool core_qed = input.get("RadPot", "core_qed", true);

  if (include_qed && qed_ok && core_qed) {
    wf.radiativePotential(x_Simple, x_Ueh, x_SEe_h, x_SEe_l, x_SEm, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Hartree-Fock core (and valence)\n\n";
  }

  // Inlcude extra potential (read in from text file):
  // Note: interpolated onto grid, but NOT extrapolated (zero outside region!)
  const auto extra_ok =
      input.check("ExtraPotential", {"filename", "factor", "beforeHF"});
  const auto ep_fname =
      input.get<std::string>("ExtraPotential", "filename", "");
  const auto ep_factor = input.get("ExtraPotential", "factor", 0.0);
  const auto ep_beforeHF = input.get("ExtraPotential", "beforeHF", false);
  const auto extra_pot =
      ep_fname != "" && std::abs(ep_factor) > 0.0 && extra_ok;
  std::vector<double> Vextra;
  if (extra_pot) {
    const auto &[x, y] = IO::FRW::readFile_xy_PoV("testIn.txt");
    Vextra = Interpolator::interpolate(x, y, wf.rgrid->r);
    qip::scale(&Vextra, ep_factor);
  }

  // Add "extra potential", before HF (core + valence)
  if (extra_pot && ep_beforeHF) {
    qip::add(&wf.vnuc, Vextra);
  }

  { // Solve Hartree equations for the core:
    IO::ChronoTimer t(" core");
    wf.hartreeFockCore(HF_method, x_Breit, str_core, eps_HF);
  }

  if (include_qed && qed_ok && !core_qed) {
    wf.radiativePotential(x_Simple, x_Ueh, x_SEe_h, x_SEe_l, x_SEm, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Valence only\n\n";
  }

  // Add "extra potential", after HF (only valence)
  if (extra_pot && !ep_beforeHF) {
    qip::add(&wf.vdir, Vextra);
  }

  // Adds effective polarision potential to direct potential
  // (After HF core, before HF valence)
  const auto Vpol_ok = input.check("dVpol", {"a_eff", "r_cut"});
  const auto a_eff = input.get("dVpol", "a_eff", 0.0);
  if (std::abs(a_eff) > 0.0 && Vpol_ok) {
    const auto r_cut = input.get("dVpol", "r_cut", 1.0);
    const auto a4 = r_cut * r_cut * r_cut * r_cut;
    auto dV = [=](auto x) { return -0.5 * a_eff / (x * x * x * x + a4); };
    for (auto i = 0u; i < wf.rgrid->num_points; ++i) {
      wf.vdir[i] += dV(wf.rgrid->r[i]);
    }
  }

  // Solve for the valence states:
  const auto valence_list =
      (wf.Ncore() < wf.Znuc() || HF_method == "KohnSham")
          ? input.get<std::string>("HartreeFock", "valence", "")
          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    IO::ChronoTimer t("  val");
    if (HF_method == "KohnSham") {
      // Need different energy Guess for Kohn-sham!
      // Also: different way of reading valence list (since core cross-over!)
      wf.localValence(valence_list, true);
    } else if (wf.core.empty()) {
      wf.localValence(valence_list);
    } else {
      wf.hartreeFockValence(valence_list);
    }
  }

  // Output Hartree Fock energies:
  std::cout << "\nHartree Fock: " << wf.identity() << "-" << wf.Anuc() << "\n";
  const auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  // Construct B-spline basis:
  const auto basis_ok =
      input.check("Basis", {"number", "order", "r0", "r0_eps", "rmax", "states",
                            "print", "positron"});
  if (basis_ok)
    wf.formBasis({input.get("Basis")});
  if (input.get("Basis", "print", false) && !wf.basis.empty()) {
    std::cout << "Basis:\n";
    wf.printBasis(wf.basis);
  }

  // Correlations: read in options
  const auto Sigma_ok = input.check(
      "Correlations",
      {"Brueckner", "energyShifts", "n_min_core", "fitTo_cm", "lambda_k", "fk",
       "io_file", "rmin", "rmax", "stride", "Feynman", "screening",
       "holeParticle", "lmax", "basis_for_Green", "basis_for_pol", "real_omega",
       "imag_omega", "include_G"});
  const bool do_energyShifts = input.get("Correlations", "energyShifts", false);
  const bool do_brueckner = input.get("Correlations", "Brueckner", false);
  const auto n_min_core = input.get("Correlations", "n_min_core", 1);
  const auto sigma_rmin = input.get("Correlations", "rmin", 1.0e-4);
  const auto sigma_rmax = input.get("Correlations", "rmax", 30.0);
  const auto default_stride = [&]() {
    // By default, choose stride such that there is 150 points over [1e-4,30]
    const auto stride =
        int(wf.rgrid->getIndex(30.0) - wf.rgrid->getIndex(1.0e-4)) / 150;
    return (stride <= 2) ? 2 : stride;
  }();
  const auto sigma_stride = input.get("Correlations", "stride", default_stride);
  // Feynman method:
  const auto sigma_Feynman = input.get("Correlations", "Feynman", false);
  const auto sigma_Screening = input.get("Correlations", "screening", false);
  const auto hole_particle = input.get("Correlations", "holeParticle", false);
  const auto sigma_lmax = input.get("Correlations", "lmax", 6);
  const auto GreenBasis = input.get("Correlations", "basis_for_Green", false);
  const auto PolBasis = input.get("Correlations", "basis_for_pol", false);
  const auto include_G = input.get("Correlations", "include_G", false);
  // force sigma_omre to be always -ve
  const auto sigma_omre = -std::abs(
      input.get("Correlations", "real_omega", -0.33 * wf.energy_gap()));

  // Imaginary omegagrid params (onlu used for Feynman)
  double w0 = 0.01;
  double wratio = 1.5;
  {
    const auto imag_om =
        input.get_list("Correlations", "imag_omega", std::vector{w0, wratio});
    if (imag_om.size() != 2) {
      std::cout
          << "ERROR: imag_omega must be a list of 2: omega_0 (first step), "
             "and omegra_ratio (ratio for log w grid)\n";
    } else {
      w0 = imag_om[0];
      wratio = imag_om[1];
    }
  }

  // Read/write sigma matrix to file:
  const auto sigma_io = input.get("Correlations", "io_file", true);
  auto sigma_file =
      input.get<std::string>("Correlations", "io_file", wf.identity());
  if (sigma_file == "true")
    sigma_file = wf.identity();
  else if (!sigma_io)
    sigma_file = "";

  // To fit Sigma to energies:
  auto fit_energies =
      input.get_list("Correlations", "fitTo_cm", std::vector<double>{});
  // energies given in cm^-1, convert to au:
  qip::scale(&fit_energies, 1.0 / PhysConst::Hartree_invcm);
  const auto lambda_k =
      input.get_list("Correlations", "lambda_k", std::vector<double>{});
  const auto fk = input.get_list("Correlations", "fk", std::vector<double>{});

  // Form correlation potential:
  if ((do_energyShifts || do_brueckner) && Sigma_ok) {
    IO::ChronoTimer t("Sigma");
    wf.formSigma(n_min_core, true, sigma_rmin, sigma_rmax, sigma_stride,
                 include_G, lambda_k, fk, sigma_file, sigma_Feynman,
                 sigma_Screening, hole_particle, sigma_lmax, GreenBasis,
                 PolBasis, sigma_omre, w0, wratio);
  }

  // Calculate + print second-order energy shifts
  if (!wf.valence.empty() && do_energyShifts && Sigma_ok) {
    IO::ChronoTimer t("de");
    wf.SOEnergyShift();
  }

  // Solve Brueckner orbitals (optionally, fit Sigma to exp energies)
  if (!wf.valence.empty() && do_brueckner && Sigma_ok) {
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
      input.check("Spectrum", {"number", "order", "r0", "r0_eps", "rmax",
                               "states", "print", "positron"});
  if (spectrum_ok)
    wf.formSpectrum({input.get("Spectrum")});
  if (input.get("Spectrum", "print", false) && !wf.spectrum.empty()) {
    std::cout << "Spectrum:\n";
    wf.printBasis(wf.spectrum);
  }

  // run each of the modules with the calculated wavefunctions
  Module::runModules(input, wf);

  return 0;
}

//******************************************************************************
