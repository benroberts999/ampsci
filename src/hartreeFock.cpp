#include "Dirac/Wavefunction.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "Modules/Module_runModules.hpp"
#include <iostream>
#include <string>
// #include "Physics/AtomData.hpp" //need for testing basis only
#include "Dirac/Operators.hpp"
#include "HF/ExternalField.hpp"

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
  auto b = input.get("Grid", "b", 4.0);
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

  // create wavefunction object
  Wavefunction wf(atom_Z, {num_points, r0, rmax, b, grid_type, du_tmp},
                  {atom_Z, atom_A, nuc_type, rrms, skint}, var_alpha);

  // Parse input for HF method
  input_ok =
      input_ok &&
      input.check("HartreeFock", {"core", "valence", "convergence", "method",
                                  "Green_H", "Green_d", "Tietz_g", "Tietz_t",
                                  "orthonormaliseValence", "sortOutput"});
  auto str_core = input.get<std::string>("HartreeFock", "core", "[]");
  auto eps_HF = input.get("HartreeFock", "convergence", 1.0e-12);
  auto HF_method = HartreeFock::parseMethod(
      input.get<std::string>("HartreeFock", "method", "HartreeFock"));

  if (!input_ok)
    return 1;

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n";
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

  { // Solve Hartree equations for the core:
    ChronoTimer t(" core");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  // Adds effective polarision potential to direct potential
  // (After HF core, before HF valence)
  auto a_eff = input.get("dV", "a_eff", 0.0);
  if (a_eff > 0) { // a=0.61 works well for Cs ns, n=6-18
    auto r_cut = input.get("dV", "r_cut", 1.0);
    auto dV = [=](double x) { return -0.5 * a_eff / (x * x * x * x + r_cut); };
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

  // run each of the modules
  Module::runModules(input, wf);

  auto h = E1Operator(wf.rgrid);
  // auto h = HyperfineOperator(1.0, 1.0, 4.0 / PhysConst::aB_fm, wf.rgrid);
  // auto h = PNCnsiOperator(5.0, 2.3, wf.rgrid);

  auto tdhf = ExternalField(&h, wf.core_orbitals,
                            NumCalc::add_vectors(wf.vnuc, wf.vdir),
                            wf.get_alpha(), 0.0);

  tdhf.solve_TDHFcore();

  auto psis = wf.getState(6, -1);
  auto psip = wf.getState(6, 1);

  auto me = h.reducedME(psis, psip);
  auto dv = tdhf.dV_ab(psis, psip);
  // / std::sqrt(h.rme3js(psis.twoj(), psip.twoj(), 1));
  // std::cout << tdhf.dV_ab(psis, psip) << "\n";
  // std::cout << tdhf.dV_ab(psis, psip) /
  //                  std::sqrt(h.rme3js(psis.twoj(), psip.twoj(), 1))
  //           << "\n";
  std::cout << me << " " << me + dv << "\n";

  // auto psis = wf.getState(3, -1);
  // auto psip = wf.getState(3, 1);
  // auto dpsis = tdhf.get_dPsi_x(psis, dPsiType::X, 1);
  // auto dpsisY = tdhf.get_dPsi_x(psis, dPsiType::Y, 1);
  //
  // auto metha = psip * dpsis;
  // auto methb = h.reducedME(psip, psis) / (psis.en - psip.en + 0.1);
  // std::cout << metha << " - " << methb << " => "
  //           << 2.0 * (metha - methb) / (metha + methb) << "\n";
  //
  // metha = dpsisY * psip;
  // methb = h.reducedME(psis, psip) / (psis.en - psip.en - 0.1);
  // std::cout << metha << " - " << methb << " => "
  //           << 2.0 * (metha - methb) / (metha + methb) << "\n";

  //*********************************************************
  //               TESTS
  //*********************************************************

  // needs: #include "Physics/AtomData.hpp" (for AtomData::listOfStates_nk)
  //
  // bool test_hf_basis = false;
  // if (test_hf_basis) {
  //   auto basis_lst = AtomData::listOfStates_nk("9spd8f");
  //   std::vector<DiracSpinor> basis = wf.core_orbitals;
  //   HartreeFock hfbasis(wf, basis, 1.0e-6);
  //   hfbasis.verbose = false;
  //   for (const auto &nk : basis_lst) {
  //     if (wf.isInCore(nk.n, nk.k))
  //       continue;
  //     basis.emplace_back(DiracSpinor(nk.n, nk.k, wf.rgrid));
  //     auto tmp_vex = std::vector<double>{};
  //     hfbasis.hf_valence_approx(basis.back(), tmp_vex);
  //   }
  //   wf.orthonormaliseOrbitals(basis, 2);
  //   wf.printValence(false, basis);
  //   std::cout << "\n Total time: " << timer.reading_str() << "\n";
  // }

  return 0;
}

//******************************************************************************
