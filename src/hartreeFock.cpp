#include "ChronoTimer.hpp"
#include "Grid.hpp"
#include "Modules/Module_runModules.hpp"
#include "Operators.hpp"
#include "Physics/AtomInfo.hpp"
#include "Physics/Nuclear.hpp"
#include "UserInput.hpp"
#include "Wavefunction.hpp"
#include <iostream>
#include <string>
//
#include "NumCalc_quadIntegrate.hpp"

int main(int argc, char *argv[]) {
  ChronoTimer timer("hartreeFock");
  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  UserInput input(input_file);

  auto atom = input.get<std::string>("Atom", "Z");
  auto Z = AtomInfo::get_z(atom);
  auto A = input.get("Atom", "A", -1);
  auto varalpha = std::sqrt(input.get("Atom", "varAlpha2", 1.0));
  if (varalpha == 0)
    varalpha = 1.0e-16;

  auto r0 = input.get("Grid", "r0", 1.0e-5);
  auto rmax = input.get("Grid", "rmax", 150.0);
  auto ngp = input.get("Grid", "ngp", 1600ul);
  auto b = input.get("Grid", "b", 3.5);
  auto grid_type = GridParameters::parseType(
      input.get<std::string>("Grid", "type", "loglinear"));
  GridParameters grid_params(ngp, r0, rmax, b, grid_type);

  A = input.get("Nucleus", "A", A); // over-writes "atom" A
  auto nuc_type =
      Nuclear::parseType(input.get<std::string>("Nucleus", "type", "Fermi"));
  auto rrms = input.get("Nucleus", "rrms", -1.0); /*<0 means lookup default*/
  auto skint = input.get("Nucleus", "skin_t", -1.0);
  Nuclear::Parameters nuc_params(Z, A, nuc_type, rrms, skint);

  Wavefunction wf(Z, grid_params, nuc_params, varalpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  auto str_core = input.get<std::string>("HartreeFock", "core");
  auto eps_HF = input.get("HartreeFock", "convergence", 0.0);
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

  // Solve Hartree equations for the core:
  {
    ChronoTimer t("Core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  // Solve for the valence states:
  auto valence_list = (wf.Ncore() < wf.Znuc())
                          ? input.get<std::string>("HartreeFock", "valence", "")
                          : "";
  timer.start();
  wf.hartreeFockValence(valence_list);
  if (!wf.valence_orbitals.empty()) {
    if (input.get("HartreeFock", "orthonormaliseValence", false))
      wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);
    std::cout << "Valence: " << timer.lap_reading_str() << "\n";
  }

  // Output results:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  auto sorted = input.get("HartreeFock", "sortOutput", true);
  wf.printCore(sorted);
  wf.printValence(sorted);

  Module::runModules(input, wf);

  // DirectHamiltonian Hd(wf.vnuc, wf.vdir, wf.get_alpha());
  // for (const auto phi : wf.valence_orbitals) {
  //   // std::cout << (Hd.matrixEl(phi, phi) + phi * wf.get_VexPsi(phi) -
  //   phi.en)
  //   // /
  //   //                  phi.en
  //   //           << "\n";
  //   printf("%.7f\n", Hd.matrixEl(phi, phi) + phi * wf.get_VexPsi(phi));
  // }

  // for (const auto phi : wf.core_orbitals) {
  //   auto df = NumCalc::derivative(phi.f, wf.rgrid.drdu, wf.rgrid.du, 1);
  //   std::cout << (2.0 / wf.get_alpha()) *
  //                    NumCalc::integrate(phi.g, df, wf.rgrid.drdu,
  //                    wf.rgrid.du)
  //             << " ";
  //   std::cout << (2.0 / wf.get_alpha() / wf.get_alpha()) *
  //                    NumCalc::integrate(phi.g, phi.g, wf.rgrid.drdu,
  //                                       wf.rgrid.du)
  //             << " ";
  //   auto invr = wf.rgrid.inverse_r();
  //   std::cout << (2.0 / wf.get_alpha()) * phi.k *
  //                    NumCalc::integrate(phi.f, phi.g, invr, wf.rgrid.drdu,
  //                                       wf.rgrid.du)
  //             << "\n";
  // }

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

  bool testpnc = false;
  if (testpnc) {
    double t = Nuclear::default_t; // approximate_t_skin(wf.Anuc());
    auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
    double c = Nuclear::c_hdr_formula_rrms_t(r_rms);
    PNCnsiOperator hpnc(c, t, wf.rgrid, -wf.Nnuc());
    E1Operator he1(wf.rgrid);

    double Ac = 2. / 6.; // angular coef
    auto a6s_i = wf.getStateIndex(6, -1);
    auto a7s_i = wf.getStateIndex(7, -1);
    auto &a6s = wf.valence_orbitals[a6s_i];
    auto &a7s = wf.valence_orbitals[a7s_i];
    std::cout << "E_pnc: " << wf.Anuc() << "-"
              << AtomInfo::atomicSymbol(wf.Znuc()) << " " << a6s.symbol()
              << " -> " << a7s.symbol() << "\n";

    double pnc = 0;
    for (int i = 0; i < 2; i++) {
      auto &tmp_orbs = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
      for (auto &np : tmp_orbs) {
        if (np.k != 1)
          continue; // p_1/2 only
        // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
        double pnc1 =
            Ac * (a7s * (he1 * np)) * (np * (hpnc * a6s)) / (a6s.en - np.en);
        double pnc2 =
            Ac * (a7s * (hpnc * np)) * (np * (he1 * a6s)) / (a7s.en - np.en);
        std::cout << "n=" << np.n << " pnc= " << pnc1 << " + " << pnc2 << " = "
                  << pnc1 + pnc2 << "\n";
        pnc += pnc1 + pnc2;
      }
    }
    std::cout << "Total= " << pnc << "\n";
    std::cout << "\n Total time: " << timer.reading_str() << "\n";
  }

  return 0;
}

//******************************************************************************
