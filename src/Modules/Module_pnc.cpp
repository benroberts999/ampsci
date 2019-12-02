#include "Modules/Module_runModules.hpp"
//
#include "Dirac/DiracOperator.hpp"
#include "Dirac/Operators.hpp"
#include "Dirac/Wavefunction.hpp"
#include "IO/UserInput.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace Module {

//******************************************************************************
void Module_testPNC(const UserInputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::PNC";

  input.checkBlock({"t", "c", "transition", "nmain"});

  auto t_dflt = Nuclear::default_t;
  auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  auto c_dflt = Nuclear::c_hdr_formula_rrms_t(r_rms);
  auto t = input.get("t", t_dflt);
  auto c = input.get("c", c_dflt);

  auto transition_str = input.get<std::string>("transition");
  std::replace(transition_str.begin(), transition_str.end(), ',', ' ');
  auto ss = std::stringstream(transition_str);
  int na, ka, nb, kb;
  ss >> na >> ka >> nb >> kb;

  auto ncore = wf.maxCore_n();
  auto main_n = input.get("nmain", ncore + 4);

  PNCnsiOperator hpnc(c, t, wf.rgrid, -wf.Nnuc());
  E1Operator he1(wf.rgrid);
  auto alpha = wf.get_alpha();

  const auto &aA = wf.getState(na, ka);
  const auto &aB = wf.getState(nb, kb);
  std::cout << "\n********************************************** \n";
  std::cout << "E_pnc: " << wf.atom() << ":   A = " << aA.symbol()
            << "   -->   B = " << aB.symbol() << "\n\n";

  auto tja = aA.twoj();
  auto tjb = aB.twoj();
  auto twom = std::min(tja, tjb);
  auto c10 = he1.rme3js(tja, tjb, twom) * hpnc.rme3js(tjb, tjb, twom);
  auto c01 = hpnc.rme3js(tja, tja, twom) * he1.rme3js(tja, tjb, twom);

  bool main_ok = false;
  int max_n_main = 0;
  {
    std::cout << "Sum-over-states method (HF valence):\n";
    std::cout << " <A|d|n><n|hw|B>/dEB + <A|hw|n><n|d|B>/dEA\n";
    double pnc = 0, core = 0, main = 0;
    for (int i = 0; i < 2; i++) {
      auto &tmp_orbs = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
      for (auto &np : tmp_orbs) {
        // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
        if (np == aB || np == aA)
          continue;
        if (hpnc.isZero(np.k, aA.k) && hpnc.isZero(np.k, aB.k))
          continue;
        double pnc1 = c10 * he1.reducedME(aA, np) * hpnc.reducedME(np, aB) /
                      (aB.en - np.en);
        double pnc2 = c01 * hpnc.reducedME(aA, np) * he1.reducedME(np, aB) /
                      (aA.en - np.en);
        printf("%7s, pnc= %12.5e + %12.5e = %12.5e\n", np.symbol().c_str(),
               pnc1, pnc2, pnc1 + pnc2);
        pnc += pnc1 + pnc2;
        if (np.n <= main_n && np.n > ncore)
          main = pnc - core;
        if (np.n == main_n)
          main_ok = true;
        if (np.n > max_n_main)
          max_n_main = np.n;
      }
      if (i == 0)
        core = pnc;
    }
    std::cout << "Core = " << core << "\n";
    if (main_ok) {
      std::cout << "Main = " << main << "\n";
      std::cout << "Tail = " << pnc - main - core << "\n";
    } else {
      std::cout << "Rest*= " << pnc - core << " ";
      std::cout << "    *(didnt have full main, nmax=" << max_n_main << ")\n";
    }
    std::cout << "Total= " << pnc << "\n";
  }

  {
    auto v = NumCalc::add_vectors(wf.vnuc, wf.vdir);
    auto v1 = 0.0, v2 = 0.0;

    auto hA_dag = hpnc.reduced_lhs(-aA.k, aA);
    auto hB = hpnc.reduced_rhs(-aB.k, aB);

    auto e1A = he1.reduced_lhs(-aB.k, aA);
    auto e1B = he1.reduced_rhs(-aA.k, aB);

    auto del_A_dag = HartreeFock::solveMixedState(hA_dag.k, aA, 0, v, alpha,
                                                  wf.core_orbitals, hA_dag);
    auto del_B = HartreeFock::solveMixedState(hB.k, aB, 0, v, alpha,
                                              wf.core_orbitals, hB);

    auto omega = aB.en - aA.en;
    auto xA = HartreeFock::solveMixedState(e1A.k, aA, omega, v, alpha,
                                           wf.core_orbitals, e1A);
    auto yB = HartreeFock::solveMixedState(e1B.k, aB, -omega, v, alpha,
                                           wf.core_orbitals, e1B);

    auto pnc1_w = c01 * he1.reducedME(del_A_dag, aB);
    auto pnc2_w = c10 * he1.reducedME(aA, del_B);

    std::cout << "\nMixed states method: \n";
    std::cout << "<dA |d| B> + <A |d| dB> = ";
    std::cout << pnc1_w << " + " << pnc2_w << " = " << pnc1_w + pnc2_w << "\n";

    auto pnc1_d = c01 * hpnc.reducedME(aA, yB);
    auto pnc2_d = c10 * hpnc.reducedME(xA, aB);
    std::cout << "<A |h|Y_B> + <X_A|h| B> = ";
    std::cout << pnc1_d << " + " << pnc2_d << " = " << pnc1_d + pnc2_d << "\n";

    v1 = pnc1_w + pnc2_w;
    v2 = pnc1_d + pnc2_d;
    std::cout << "eps = " << (v1 - v2) / (0.5 * (v1 + v2)) << "\n";

    // Find the core+main contributions by forcing mixed-states to be orthoganal
    // to the core/main states:
    // orthog wrt core:
    wf.orthogonaliseWrt(del_A_dag, wf.core_orbitals);
    wf.orthogonaliseWrt(del_B, wf.core_orbitals);
    // Core contribution:
    auto pnc1_c = pnc1_w - c01 * he1.reducedME(del_A_dag, aB);
    auto pnc2_c = pnc2_w - c10 * he1.reducedME(aA, del_B);
    // Further orthog wrt 'main' part of valence (now orthog to core+main)
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
    // Main contribution:
    auto pnc1_m = pnc1_w - pnc1_c - c10 * he1.reducedME(del_A_dag, aB);
    auto pnc2_m = pnc2_w - pnc2_c - c01 * he1.reducedME(aA, del_B);

    std::cout << "\n<dA'|d| B>  +  <A |d|dB'>  (by force orthog):\n";
    std::cout << "core : " << pnc1_c + pnc2_c << "\n";
    if (main_ok) {
      std::cout << "main : " << pnc1_m + pnc2_m << "\n";
    } else {
      std::cout << "main*: " << pnc1_m + pnc2_m << " ";
      std::cout << "    *(didnt have full main, nmax=" << max_n_main << ")\n";
    }
    std::cout << "tail : "
              << pnc1_w - pnc1_m - pnc1_c + pnc2_w - pnc2_m - pnc2_c << "\n";
    std::cout << "Total: " << v1 << "\n";
  }
}

} // namespace Module
