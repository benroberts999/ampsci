#include "Dirac/Wavefunction.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "Maths/Grid.hpp"
#include "Modules/Module_runModules.hpp"
#include "Physics/AtomInfo.hpp"
#include "Physics/Nuclear.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  ChronoTimer timer("\nhartreeFock");
  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  UserInput input(input_file);

  // Get + setup atom parameters
  auto atom = input.get<std::string>("Atom", "Z");
  auto Z = AtomInfo::get_z(atom);
  auto A = input.get("Atom", "A", -1);
  auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
  if (varAlpha2 <= 0)
    varAlpha2 = 1.0e-25;
  auto varalpha = std::sqrt(varAlpha2);

  // Get + setup Grid parameters
  auto r0 = input.get("Grid", "r0", 1.0e-5);
  auto rmax = input.get("Grid", "rmax", 150.0);
  auto ngp = input.get("Grid", "ngp", 1600ul);
  auto b = input.get("Grid", "b", 4.0);
  auto grid_type = GridParameters::parseType(
      input.get<std::string>("Grid", "type", "loglinear"));
  GridParameters grid_params(ngp, r0, rmax, b, grid_type);

  // Get + setup nuclear parameters
  A = input.get("Nucleus", "A", A); // over-writes "atom" A
  auto nuc_type =
      Nuclear::parseType(input.get<std::string>("Nucleus", "type", "Fermi"));
  auto rrms = input.get("Nucleus", "rrms", -1.0); /*<0 means lookup default*/
  auto skint = input.get("Nucleus", "skin_t", -1.0);
  Nuclear::Parameters nuc_params(Z, A, nuc_type, rrms, skint);

  // create wavefunction object
  Wavefunction wf(Z, grid_params, nuc_params, varalpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  auto str_core = input.get<std::string>("HartreeFock", "core");
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
    std::cout << "Using Greens Tietz Potential\n";
  } else if (HF_method == HFMethod::Hartree) {
    std::cout << "Using Hartree Method (no Exchange)\n";
  }

  { // Solve Hartree equations for the core:
    ChronoTimer t("Core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  // Solve for the valence states:
  auto valence_list = (wf.Ncore() < wf.Znuc())
                          ? input.get<std::string>("HartreeFock", "valence", "")
                          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    ChronoTimer t("Valence");
    wf.hartreeFockValence(valence_list);
    if (input.get("HartreeFock", "orthonormaliseValence", false))
      wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);
  }

  // Output results:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  // run each of the modules
  Module::runModules(input, wf);

  //*********************************************************
  //               TESTS
  //*********************************************************

  bool test_hf_basis = false;
  if (test_hf_basis) {
    auto basis_lst = AtomInfo::listOfStates_nk("9spd8f");
    std::vector<DiracSpinor> basis = wf.core_orbitals;
    HartreeFock hfbasis(wf, basis, 1.0e-6);
    hfbasis.verbose = false;
    for (const auto &nk : basis_lst) {
      if (wf.isInCore(nk.n, nk.k))
        continue;
      basis.emplace_back(DiracSpinor(nk.n, nk.k, wf.rgrid));
      auto tmp_vex = std::vector<double>{};
      hfbasis.solveValence(basis.back(), tmp_vex);
    }
    wf.orthonormaliseOrbitals(basis, 2);
    wf.printValence(false, basis);
    std::cout << "\n Total time: " << timer.reading_str() << "\n";
  }

  return 0;
}

//******************************************************************************
