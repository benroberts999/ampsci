#include "ChronoTimer.hpp"
#include "DiracOperator.hpp"
#include "HartreeFockClass.hpp"
#include "Module_runModules.hpp"
#include "Operators.hpp"
#include "Physics/AtomInfo.hpp"
#include "Physics/Nuclear.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "UserInput.hpp"
#include "Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  ChronoTimer timer;

  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  UserInput input(input_file);
  auto atom = input.get<std::string>("Atom", "Z");
  auto varalpha = sqrt(input.get("Atom", "varAlpha2", 1.0));
  auto r0 = input.get("Grid", "r0", 1.0e-5);
  auto rmax = input.get("Grid", "rmax", 150.0);
  auto ngp = input.get("Grid", "ngp", 1600);
  auto A = input.get("Nucleus", "A", -1);

  auto Z = AtomInfo::get_z(atom);
  Wavefunction wf(Z, A, ngp, r0, rmax, varalpha);

  std::cout << "\nRunning for ";
  std::cout << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";

  // Parse input for HF method
  auto str_core = input.get<std::string>("HartreeFock", "Core");
  auto eps_HF = input.get("HartreeFock", "Convergence", 0.0);
  auto HF_method = HartreeFock::parseMethod(
      input.get<std::string>("HartreeFock", "Method", "HartreeFock"));

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
  timer.start();
  HartreeFock hf(HF_method, wf, str_core, eps_HF, H_d, g_t);

  double core_energy = hf.calculateCoreEnergy();
  std::cout << "core: " << timer.lap_reading_str() << "\n";

  // Create list of valence states to solve for
  auto in_valence_str = input.get<std::string>("HartreeFock", "valence", "");
  if (wf.Ncore() >= wf.Znuc())
    in_valence_str = "";
  auto val_lst = AtomInfo::listOfStates_nk(in_valence_str);

  // Solve for the valence states:
  timer.start();
  for (const auto &nk : val_lst) {
    if (!wf.isInCore(nk.n, nk.k))
      hf.solveNewValence(nk.n, nk.k);
  }
  if (val_lst.size() > 0)
    std::cout << "Valence: " << timer.lap_reading_str() << "\n";

  if (input.get("HartreeFock", "OrthonormaliseValence", false))
    wf.orthonormaliseOrbitals(wf.valence_orbitals, 2);

  // Output results:
  std::cout << "\nHartree Fock: " << wf.atom() << "\n";
  bool sorted = true;
  wf.printCore(sorted);
  std::cout << "E_core = " << core_energy
            << " au;  = " << core_energy * PhysConst::Hartree_invcm << "/cm\n";
  wf.printValence(sorted);

  {
    auto modules = input.module_list();
    for (const auto &module : modules) {
      timer.start();
      Module::runModule(module, input, wf, hf);
      // std::cout << module << " time: " << timer.lap_reading_str() << "\n";
    }
  }

  std::cout << "\nTotal time: " << timer.reading_str() << "\n";

  //*********************************************************
  //               TESTS
  //*********************************************************

  bool test_hf_basis = false;
  if (test_hf_basis) {
    auto basis_lst = AtomInfo::listOfStates_nk("9spdf");
    std::vector<DiracSpinor> basis = wf.core_orbitals;
    HartreeFock hfbasis(wf, basis, eps_HF);
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

  // bool test_hfs = false;
  // if (test_hfs) {
  //   std::cout << "\nHyperfine:\n";
  //
  //   double muN = Nuclear::find_mu(wf.Znuc(), wf.Anuc());
  //   double IN = Nuclear::find_spin(wf.Znuc(), wf.Anuc());
  //   auto r_rms_fm = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  //   // 4.1989
  //   auto r_rms = r_rms_fm / PhysConst::aB_fm;
  //   std::cout << "mu=" << muN << ", I=" << IN << " ,r=" << r_rms_fm << "\n";
  //
  //   std::cout << "Gridpoints below Rrms: " << wf.rgrid.getIndex(r_rms) <<
  //   "\n";
  //
  //   // example for using lambda
  //   auto l1 = [](double r, double) { return 1. / (r * r); };
  //   // auto l1 = [](double r, double rN) { return r > rN ? 1. / (r * r) : 0.;
  //   // };
  //   HyperfineOperator vhfs(muN, IN, r_rms, wf.rgrid, l1);
  //   for (auto phi : wf.valence_orbitals) {
  //     auto A_tmp = phi * (vhfs * phi);
  //     double j = phi.j();
  //     auto factor = PhysConst::Hartree_MHz * phi.k / (j * (j + 1.));
  //     std::cout << phi.symbol() << ": ";
  //     std::cout << A_tmp * factor << "\n";
  //   }
  // }

  return 0;
}

//******************************************************************************
