#include "Adams/Adams_Greens.hpp"
#include "Dirac/Operators.hpp"
#include "Dirac/Wavefunction.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "Modules/Module_runModules.hpp"
#include <iostream>
#include <string>
// #include "Physics/AtomInfo.hpp" //need for testing basis only

int main(int argc, char *argv[]) {
  ChronoTimer timer("\nhartreeFock");
  std::string input_file = (argc > 1) ? argv[1] : "hartreeFock.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  UserInput input(input_file);

  // Get + setup atom parameters
  auto atom_Z = input.get<std::string>("Atom", "Z");
  auto atom_A = input.get("Atom", "A", -1);
  auto var_alpha = [&]() {
    auto varAlpha2 = input.get("Atom", "varAlpha2", 1.0);
    return (varAlpha2 > 0) ? std::sqrt(varAlpha2) : 1.0e-25;
  }();

  // Get + setup Grid parameters
  auto r0 = input.get("Grid", "r0", 1.0e-5);
  auto rmax = input.get("Grid", "rmax", 150.0);
  auto num_points = input.get("Grid", "num_points", 1600ul);
  auto b = input.get("Grid", "b", 4.0);
  auto grid_type = input.get<std::string>("Grid", "type", "loglinear");
  if (b <= r0 || b >= rmax)
    grid_type = "logarithmic";

  // Get + setup nuclear parameters
  atom_A = input.get("Nucleus", "A", atom_A); // over-writes "atom" A
  auto nuc_type = input.get<std::string>("Nucleus", "type", "Fermi");
  auto rrms = input.get("Nucleus", "rrms", -1.0); /*<0 means lookup default*/
  auto skint = input.get("Nucleus", "skin_t", -1.0);

  // create wavefunction object
  Wavefunction wf(atom_Z, {num_points, r0, rmax, b, grid_type},
                  {atom_Z, atom_A, nuc_type, rrms, skint}, var_alpha);

  std::cout << "\nRunning for " << wf.atom() << "\n"
            << wf.nuclearParams() << "\n"
            << wf.rgrid.gridParameters() << "\n"
            << "********************************************************\n\n";

  // Parse input for HF method
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

  { // Solve Hartree equations for the core:
    ChronoTimer t(" ");
    wf.hartreeFockCore(HF_method, str_core, eps_HF, H_d, g_t);
  }

  // Solve for the valence states:
  auto valence_list = (wf.Ncore() < wf.Znuc())
                          ? input.get<std::string>("HartreeFock", "valence", "")
                          : "";
  if (valence_list != "") {
    // 'if' is only for output format, nothing bad happens if below are called
    ChronoTimer t(" ");
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

  // // needs: #include "Physics/AtomInfo.hpp" (for AtomInfo::listOfStates_nk)
  // bool test_hf_basis = false;
  // if (test_hf_basis) {
  //   auto basis_lst = AtomInfo::listOfStates_nk("9spd8f");
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
  int nb = 7;
  int kb = -1;
  const auto &aA = wf.getState(6, -1);
  const auto &aB = wf.getState(nb, kb);
  const auto &a6p1 = wf.getState(6, 1);   // dummy
  const auto &a6p3 = wf.getState(6, -kb); // dummy

  // auto k1 = 1;
  // auto k2 = -2;
  auto alpha = wf.get_alpha();

  E1Operator he1(wf.rgrid);
  auto hpnc = PNCnsiOperator(5.67073, 2.3, wf.rgrid, -wf.Nnuc());

  auto v = wf.vnuc;
  for (unsigned i = 0; i < v.size(); i++) {
    v[i] += wf.vdir[i];
  }

  auto v1 = 0.0, v2 = 0.0;
  {
    double sa1 = -1.0 * std::pow(-1, (aA.twoj() - a6p1.twoj()) / 2);
    std::cout << "sa1=" << sa1 << "\n";
    auto hA_dag = sa1 * hpnc.reduced_rhs(a6p1, aA);
    /*minus 1 comes from cc of REDUCED matrix element! (have 'lhs' instead)*/
    auto hB = +1.0 * hpnc.reduced_rhs(a6p3, aB);

    double sa2 = std::pow(-1, (aA.twoj() - a6p3.twoj()) / 2);
    std::cout << "sa2=" << sa2 << "\n";
    auto e1A = sa2 * he1.reduced_rhs(a6p3, aA);
    auto e1B = he1.reduced_rhs(a6p1, aB);

    auto vxa = wf.m_pHF->get_vex(aA);
    auto vxb = wf.m_pHF->get_vex(aB);

    auto del_A_dag = HartreeFock::solveMixedState(aA, hA_dag.k, 0, v, alpha,
                                                  wf.core_orbitals, hA_dag, {});
    auto del_B = HartreeFock::solveMixedState(aB, hB.k, 0, v, alpha,
                                              wf.core_orbitals, hB, {});

    auto om = std::fabs(aA.en - aB.en);
    auto xA = HartreeFock::solveMixedState(aA, e1A.k, -om, v, alpha,
                                           wf.core_orbitals, e1A, {});
    auto yB = HartreeFock::solveMixedState(aB, e1B.k, om, v, alpha,
                                           wf.core_orbitals, e1B, {});

    auto tja = aA.twoj();
    auto tjb = aB.twoj();
    auto twom = std::min(tja, tjb);
    auto ss1 = std::pow(-1, (tja + tja - 2 * twom) / 2);
    auto ss2 = std::pow(-1, (tja + tjb - 2 * twom) / 2);
    std::cout << ss1 << " " << ss2 << " (ss)\n";
    auto c10 = Wigner::threej_2(tjb, 2, tja, -twom, 0, twom) *
               Wigner::threej_2(tja, 0, tja, -twom, 0, twom) * ss2;
    auto c01 = Wigner::threej_2(tjb, 0, tjb, -twom, 0, twom) *
               Wigner::threej_2(tjb, 2, tja, -twom, 0, twom) * ss1;

    // XXX "reducedME" not quite right! Already 3js in del_A_dag!
    auto pnc1_w = c10 * he1.reducedME(del_A_dag, aB);
    auto pnc2_w = c01 * he1.reducedME(aA, del_B);

    std::cout << "\n<dA |d| B> + <A |d| dB> = ";
    std::cout << pnc1_w << " + " << pnc2_w << " = " << pnc1_w + pnc2_w << "\n";

    auto pnc1_d = c01 * hpnc.reducedME(aA, yB);
    auto pnc2_d = c10 * hpnc.reducedME(xA, aB);
    std::cout << "<A |h|Y_B> + <X_A|h| B> = ";
    std::cout << pnc1_d << " + " << pnc2_d << " = " << pnc1_d + pnc2_d << "\n";

    v1 = pnc1_w + pnc2_w;
    v2 = pnc1_d + pnc2_d;
    std::cout << "eps = " << (v1 - v2) / (0.5 * (v1 + v2)) << "\n";

    // orthog wrt core:
    wf.orthogonaliseWrtCore(del_A_dag);
    wf.orthogonaliseWrtCore(del_B);
    // Core contribution:
    auto pnc1_c = pnc1_w - c10 * he1.reducedME(del_A_dag, aB);
    auto pnc2_c = pnc2_w - c01 * he1.reducedME(aA, del_B);
    // orthog wrt 'main':
    auto ncore = wf.maxCore_n();
    auto main_n = ncore + 4;
    for (const auto &phiv : wf.valence_orbitals) {
      if (phiv.n > main_n)
        continue;
      if (phiv.k == del_A_dag.k) {
        del_A_dag -= (del_A_dag * phiv) * phiv;
      }
      if (phiv.k == del_B.k) {
        del_B -= (del_B * phiv) * phiv;
      }
    }

    // Min contribution:
    auto pnc1_m = pnc1_w - pnc1_c - c10 * he1.reducedME(del_A_dag, aB);
    auto pnc2_m = pnc2_w - pnc2_c - c01 * he1.reducedME(aA, del_B);

    std::cout << "\n<dA'|d| B>  +  <A |d|dB'>  (by force orthog):\n";
    std::cout << "core : " << pnc1_c + pnc2_c << "\n";
    std::cout << "main : " << pnc1_m + pnc2_m << "\n";
    std::cout << "tail : "
              << pnc1_w - pnc1_m - pnc1_c + pnc2_w - pnc2_m - pnc2_c << "\n";
    std::cout << "Total: " << v1 << "\n";
  }
  return 0;
}

//******************************************************************************
