#include "Wavefunction/Wavefunction.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FileIO_fileReadWrite.hpp" //for 'ExtraPotential'
#include "IO/UserInput.hpp"
#include "Maths/Interpolator.hpp"          //for 'ExtraPotential'
#include "Maths/NumCalc_quadIntegrate.hpp" //for 'ExtraPotential'
#include "Modules/Module_runModules.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  ChronoTimer timer("\nhartreeFock");
  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  UserInput input(input_file);

  // Get + setup atom parameters
  auto input_ok = input.check("Atom", {"Z", "A", "varAlpha2"});
  auto atom_Z = input.get<std::string>("Atom", "Z");
  auto atom_A = input.get("Atom", "A", -1);
  auto var_alpha = [&]() {
    auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
    return (varAlpha2 > 0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Get + setup Grid parameters
  input_ok = input_ok && input.check("Grid", {"r0", "rmax", "num_points",
                                              "type", "b", "fixed_du"});
  auto r0 = input.get("Grid", "r0", 1.0e-5);
  auto rmax = input.get("Grid", "rmax", 150.0);
  auto num_points = input.get("Grid", "num_points", 1600ul);
  auto du_tmp = input.get("Grid", "fixed_du", -1.0); // >0 means calc num_points
  if (du_tmp > 0)
    num_points = 0;
  auto b = input.get("Grid", "b", 0.33 * rmax);
  auto grid_type = input.get<std::string>("Grid", "type", "loglinear");
  if (b <= r0 || b >= rmax)
    grid_type = "logarithmic";

  // Get + setup nuclear parameters
  input_ok =
      input_ok && input.check("Nucleus", {"A", "rrms", "skin_t", "type"});
  atom_A = input.get("Nucleus", "A", atom_A); // over-writes "atom" A
  auto nuc_type = input.get<std::string>("Nucleus", "type", "Fermi");
  auto rrms = input.get("Nucleus", "rrms", -1.0); /*<0 means lookup default*/
  auto skint = input.get("Nucleus", "skin_t", -1.0);

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
  auto str_core = input.get<std::string>("HartreeFock", "core", "[]");
  auto eps_HF = input.get("HartreeFock", "convergence", 1.0e-12);
  auto HF_method = HartreeFock::parseMethod(
      input.get<std::string>("HartreeFock", "method", "HartreeFock"));

  // For when using Hartree, or a parametric potential:
  double H_d = 0.0, g_t = 0.0;
  if (HF_method == HFMethod::GreenPRM) {
    H_d = input.get("HartreeFock", "Green_H", 0.0);
    g_t = input.get("HartreeFock", "Green_d", 0.0);
    std::cout << "Using Greens Parametric Potential\n";
  } else if (HF_method == HFMethod::TietzPRM) {
    H_d = input.get("HartreeFock", "Tietz_g", 0.0);
    g_t = input.get("HartreeFock", "Tietz_t", 0.0);
    std::cout << "Using Tietz Parametric Potential\n";
  } else if (HF_method == HFMethod::Hartree) {
    std::cout << "Using Hartree Method (no Exchange)\n";
  }

  // Inlcude QED radiatve potential?
  input_ok = input_ok && input.check("RadPot", {"Ueh", "SE_h", "SE_l", "SE_m",
                                                "rcut", "scale_rN"});
  auto x_Ueh = input.get("RadPot", "Ueh", 0.0);
  auto x_SEe_h = input.get("RadPot", "SE_h", 0.0);
  auto x_SEe_l = input.get("RadPot", "SE_l", 0.0);
  auto x_SEm = input.get("RadPot", "SE_m", 0.0);
  auto rcut = input.get("RadPot", "rcut", 1.0);
  auto scale_rN = input.get("RadPot", "scale_rN", 1.0);
  if (input_ok)
    wf.radiativePotential(x_Ueh, x_SEe_h, x_SEe_l, x_SEm, rcut, scale_rN);

  // Inlcude extra potential (read in from text file):
  // Note: interpolated onto grid, but NOT extrapolated (zero outside region!)
  input_ok = input_ok &&
             input.check("ExtraPotential", {"filename", "factor", "beforeHF"});
  auto ep_fname = input.get<std::string>("ExtraPotential", "filename", "");
  auto ep_factor = input.get("ExtraPotential", "factor", 0.0);
  auto ep_beforeHF = input.get("ExtraPotential", "beforeHF", false);
  auto extra_pot = ep_fname != "" && std::abs(ep_factor) > 0.0;
  std::vector<double> Vextra;
  if (extra_pot) {
    const auto &[x, y] = FileIO::readFile_xy_PoV("testIn.txt");
    Vextra = Interpolator::interpolate(x, y, wf.rgrid.r);
    NumCalc::scaleVec(Vextra, ep_factor);
  }

  // Add "extra potential", before HF (core + valence)
  if (extra_pot && ep_beforeHF) {
    wf.vnuc = NumCalc::add_vectors(wf.vnuc, Vextra);
  }

  { // Solve Hartree equations for the core:
    ChronoTimer t(" core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  // Add "extra potential", after HF (only valence)
  if (extra_pot && !ep_beforeHF) {
    wf.vdir = NumCalc::add_vectors(wf.vdir, Vextra);
  }

  // Adds effective polarision potential to direct potential
  // (After HF core, before HF valence)
  auto a_eff = input.get("dVpol", "a_eff", 0.0);
  if (std::abs(a_eff) > 0.0) {
    auto r_cut = input.get("dVpol", "r_cut", 1.0);
    auto a4 = r_cut * r_cut * r_cut * r_cut;
    auto dV = [=](auto x) { return -0.5 * a_eff / (x * x * x * x + a4); };
    // auto dV = [=](auto x) { return -0.5 * a_eff / std::pow(x + r_cut, 4); };
    for (auto i = 0u; i < wf.rgrid.num_points; ++i) {
      wf.vdir[i] += dV(wf.rgrid.r[i]);
    }
  }

  // Solve for the valence states:
  auto valence_list = (wf.Ncore() < wf.Znuc())
                          ? input.get<std::string>("HartreeFock", "valence", "")
                          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    ChronoTimer t("  val");
    wf.hartreeFockValence(valence_list);
    if (input.get("HartreeFock", "orthonormaliseValence", false))
      wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);
  }

  // Output results:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  input.check("Basis",
              {"number", "order", "r0", "rmax", "states", "print", "positron"});
  auto n_spl = input.get("Basis", "number", 0ul);
  auto k_spl = input.get("Basis", "order", 0ul);
  auto r0_spl = input.get("Basis", "r0", 0.0);
  auto rmax_spl = input.get("Basis", "rmax", 0.0);
  auto basis_states = input.get<std::string>("Basis", "states", "");
  auto print = input.get("Basis", "print", false);
  auto positronQ = input.get("Basis", "positron", false);
  if (n_spl > 0) {
    wf.formBasis(basis_states, n_spl, k_spl, r0_spl, rmax_spl, positronQ);
    if (print)
      wf.printBasis();
  }

  // run each of the modules
  Module::runModules(input, wf);

  return 0;
}

//******************************************************************************
