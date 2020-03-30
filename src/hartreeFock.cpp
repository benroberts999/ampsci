#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/UserInput.hpp"
#include "Maths/Interpolator.hpp"          //for 'ExtraPotential'
#include "Maths/NumCalc_quadIntegrate.hpp" //for 'ExtraPotential'
#include "Modules/Module_runModules.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  IO::ChronoTimer timer("\nhartreeFock");
  const std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  IO::UserInput input(input_file);

  // Get + setup atom parameters
  auto input_ok = input.check("Atom", {"Z", "A", "varAlpha2"});
  const auto atom_Z = input.get<std::string>("Atom", "Z");
  auto atom_A = input.get("Atom", "A", -1);
  const auto var_alpha = [&]() {
    const auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
    return (varAlpha2 > 0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Get + setup Grid parameters
  input_ok = input_ok && input.check("Grid", {"r0", "rmax", "num_points",
                                              "type", "b", "fixed_du"});
  const auto r0 = input.get("Grid", "r0", 1.0e-5);
  const auto rmax = input.get("Grid", "rmax", 150.0);
  const auto du_tmp =
      input.get("Grid", "fixed_du", -1.0); // >0 means calc num_points
  const auto num_points =
      (du_tmp > 0) ? 0 : input.get("Grid", "num_points", 1600ul);
  const auto b = input.get("Grid", "b", 0.33 * rmax);
  const auto grid_type =
      (b <= r0 || b >= rmax)
          ? "logarithmic"
          : input.get<std::string>("Grid", "type", "loglinear");

  // Get + setup nuclear parameters
  input_ok =
      input_ok && input.check("Nucleus", {"A", "rrms", "skin_t", "type"});
  atom_A = input.get("Nucleus", "A", atom_A); // over-writes "atom" A
  const auto nuc_type = input.get<std::string>("Nucleus", "type", "Fermi");
  const auto rrms = input.get("Nucleus", "rrms", -1.0); // <0 means get default
  const auto skint = input.get("Nucleus", "skin_t", -1.0);

  // Create wavefunction object
  Wavefunction wf(atom_Z, {num_points, r0, rmax, b, grid_type, du_tmp},
                  {atom_Z, atom_A, nuc_type, rrms, skint}, var_alpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  input_ok =
      input_ok &&
      input.check("HartreeFock", {"core", "valence", "convergence", "method",
                                  "Green_H", "Green_d", "Tietz_g", "Tietz_t",
                                  "orthonormaliseValence", "sortOutput"});
  if (!input_ok)
    return 1;
  const auto str_core = input.get<std::string>("HartreeFock", "core", "[]");
  const auto eps_HF = input.get("HartreeFock", "convergence", 1.0e-12);
  const auto HF_method = HF::parseMethod(
      input.get<std::string>("HartreeFock", "method", "HartreeFock"));

  // For when using Hartree, or a parametric potential:
  // XXX Move into new block!?
  double H_d = 0.0, g_t = 0.0;
  if (HF_method == HF::Method::GreenPRM) {
    H_d = input.get("HartreeFock", "Green_H", 0.0);
    g_t = input.get("HartreeFock", "Green_d", 0.0);
    std::cout << "Using Greens Parametric Potential\n";
  } else if (HF_method == HF::Method::TietzPRM) {
    H_d = input.get("HartreeFock", "Tietz_g", 0.0);
    g_t = input.get("HartreeFock", "Tietz_t", 0.0);
    std::cout << "Using Tietz Parametric Potential\n";
  } else if (HF_method == HF::Method::Hartree) {
    std::cout << "Using Hartree Method (no Exchange)\n";
  }

  // Inlcude QED radiatve potential
  const auto qed_ok =
      input.check("RadPot", {"Simple", "Ueh", "SE_h", "SE_l", "SE_m", "rcut",
                             "scale_rN", "scale_l", "core_qed"});
  const auto x_Simple = input.get("RadPot", "Simple", 0.0);
  const auto x_Ueh = input.get("RadPot", "Ueh", 0.0);
  const auto x_SEe_h = input.get("RadPot", "SE_h", 0.0);
  const auto x_SEe_l = input.get("RadPot", "SE_l", 0.0);
  const auto x_SEm = input.get("RadPot", "SE_m", 0.0);
  const auto rcut = input.get("RadPot", "rcut", 1.0);
  const auto scale_rN = input.get("RadPot", "scale_rN", 1.0);
  const auto x_spd = input.get_list("RadPot", "scale_l", std::vector{1.0});
  const bool core_qed = input.get("RadPot", "core_qed", true);
  if (qed_ok && core_qed) {
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
    Vextra = Interpolator::interpolate(x, y, wf.rgrid.r);
    NumCalc::scaleVec(Vextra, ep_factor);
  }

  // Add "extra potential", before HF (core + valence)
  if (extra_pot && ep_beforeHF) {
    wf.vnuc = NumCalc::add_vectors(wf.vnuc, Vextra);
  }

  { // Solve Hartree equations for the core:
    IO::ChronoTimer t(" core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  if (qed_ok && !core_qed) {
    wf.radiativePotential(x_Simple, x_Ueh, x_SEe_h, x_SEe_l, x_SEm, rcut,
                          scale_rN, x_spd);
    std::cout << "Including QED into Valence only\n\n";
  }

  // Add "extra potential", after HF (only valence)
  if (extra_pot && !ep_beforeHF) {
    wf.vdir = NumCalc::add_vectors(wf.vdir, Vextra);
  }

  // Adds effective polarision potential to direct potential
  // (After HF core, before HF valence)
  const auto Vpol_ok = input.check("dVpol", {"a_eff", "r_cut"});
  const auto a_eff = input.get("dVpol", "a_eff", 0.0);
  if (std::abs(a_eff) > 0.0 && Vpol_ok) {
    const auto r_cut = input.get("dVpol", "r_cut", 1.0);
    const auto a4 = r_cut * r_cut * r_cut * r_cut;
    auto dV = [=](auto x) { return -0.5 * a_eff / (x * x * x * x + a4); };
    for (auto i = 0u; i < wf.rgrid.num_points; ++i) {
      wf.vdir[i] += dV(wf.rgrid.r[i]);
    }
  }

  // Solve for the valence states:
  const auto valence_list =
      (wf.Ncore() < wf.Znuc())
          ? input.get<std::string>("HartreeFock", "valence", "")
          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    IO::ChronoTimer t("  val");
    wf.hartreeFockValence(valence_list);
    if (input.get("HartreeFock", "orthonormaliseValence", false))
      wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);
  }

  // Output Hartree Fock energies:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  const auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  // Construct B-spline basis:
  const auto basis_ok =
      input.check("Basis", {"number", "order", "r0", "r0_eps", "rmax", "states",
                            "print", "positron"});
  const auto n_spl = input.get("Basis", "number", 0ul);
  const auto k_spl = input.get("Basis", "order", 0ul);
  const auto r0_spl = input.get("Basis", "r0", 0.0);
  const auto r0_eps = input.get("Basis", "r0_eps", 0.0);
  const auto rmax_spl = input.get("Basis", "rmax", 0.0);
  const auto basis_states = input.get<std::string>("Basis", "states", "");
  const auto print_basisQ = input.get("Basis", "print", false);
  const auto positronQ = input.get("Basis", "positron", false);
  if (n_spl > 0 && basis_ok) {
    wf.formBasis(basis_states, n_spl, k_spl, r0_spl, r0_eps, rmax_spl,
                 positronQ);
    if (print_basisQ)
      wf.printBasis();
  }

  // Correlations:
  const auto Sigma_ok =
      input.check("Correlations", {"Brueckner", "energyShifts", "n_min_core"});
  const bool do_energyShifts =
      Sigma_ok && input.get("Correlations", "energyShifts", false);
  const bool do_brueckner =
      Sigma_ok && input.get("Correlations", "Brueckner", false);
  const auto n_min_core = input.get("Correlations", "n_min_core", 1);
  // Just energy shifts
  if (!wf.valence_orbitals.empty() && do_energyShifts) {
    IO::ChronoTimer t(" de");
    wf.SOEnergyShift(n_min_core);
  }
  // Brueckner orbitals
  if (!wf.valence_orbitals.empty() && do_brueckner) {
    std::cout
        << "\nConstructing correlation potential for Brueckner orbitals:\n"
        << std::flush;
    IO::ChronoTimer t(" Br");
    wf.hartreeFockBrueckner(n_min_core);
  }
  // Print out info for new "Brueckner" valence orbitals:
  if (!wf.valence_orbitals.empty() && do_brueckner) {
    std::cout << "\nBrueckner orbitals:\n";
    wf.printValence(sorted);
  }

  // run each of the modules
  Module::runModules(input, wf);

  return 0;
}

//******************************************************************************
