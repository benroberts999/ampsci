#include "Modules/Module_runModules.hpp"
//
#include "Operators/DiracOperator.hpp"
#include "Operators/Operators.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "IO/UserInput.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
//
#include "HF/ExternalField.hpp"

namespace Module {

//******************************************************************************
void Module_testPNC(const UserInputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::PNC";

  input.checkBlock({"t", "c", "transition", "nmain", "rpa", "omega"});

  bool print_all = false;

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

  auto rpaQ = input.get("rpa", false);

  PNCnsiOperator hpnc(c, t, wf.rgrid, -wf.Nnuc());
  E1Operator he1(wf.rgrid);
  auto alpha = wf.get_alpha();

  auto pA = wf.getState(na, ka);
  auto pB = wf.getState(nb, kb);

  if (!pA || !pB) {
    std::cerr << "\nFAIL in Module:PNC\n"
              << "Couldn't find requested state: (" << transition_str
              << "). Is it in valence list?\n";
    return;
  }
  const auto &aA = *pA;
  const auto &aB = *pB;

  std::cout << "\n********************************************** \n";
  std::cout << "E_pnc (B->A): " << wf.atom() << ":   A = " << aA.symbol()
            << " ,  B = " << aB.symbol() << "\n"
            << "z-component, z=min(ja,jb). units: i(-Qw/N)10^-11."
            << "\n\n";

  auto dVE1 = ExternalField(&he1, wf.core_orbitals,
                            NumCalc::add_vectors(wf.vnuc, wf.vdir), alpha);
  auto dVpnc = ExternalField(&hpnc, wf.core_orbitals,
                             NumCalc::add_vectors(wf.vnuc, wf.vdir), alpha);
  if (rpaQ) {
    auto omega_dflt = std::abs(aA.en - aB.en);
    auto omega = input.get("omega", omega_dflt);
    dVE1.solve_TDHFcore(omega);
    dVpnc.solve_TDHFcore(0.0);
  }
  const bool dVconj = aA.en > aB.en ? true : false;

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
        auto coreQ = wf.isInCore(np);
        auto dAp = he1.reducedME(aA, np) + dVE1.dV_ab(aA, np, dVconj);
        auto hpB = hpnc.reducedME(np, aB) + dVpnc.dV_ab(np, aB, dVconj);
        auto hAp = hpnc.reducedME(aA, np) + dVpnc.dV_ab(aA, np, dVconj);
        auto dpB = he1.reducedME(np, aB) + dVE1.dV_ab(np, aB, dVconj);
        double pnc1 = c10 * dAp * hpB / (aB.en - np.en);
        double pnc2 = c01 * hAp * dpB / (aA.en - np.en);
        if (print_all)
          printf("%7s, pnc= %12.5e + %12.5e = %12.5e\n", np.symbol().c_str(),
                 pnc1, pnc2, pnc1 + pnc2);
        pnc += pnc1 + pnc2;
        if (np.n <= main_n && !coreQ)
          main += pnc1 + pnc2;
        if (np.n == main_n)
          main_ok = true;
        if (np.n > max_n_main)
          max_n_main = np.n; // largest n used (only used for printing)
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

  if (!wf.basis.empty()) {
    std::cout << "\nSum-over-states method (basis):\n";
    std::cout << " <A|d|n><n|hw|B>/dEB + <A|hw|n><n|d|B>/dEA\n";
    double pnc = 0, core = 0, main = 0;
    for (auto &np : wf.basis) {
      if (np == aB || np == aA)
        continue;
      if (hpnc.isZero(np.k, aA.k) && hpnc.isZero(np.k, aB.k))
        continue;
      auto coreQ = wf.isInCore(np);
      // Exclude np from dV ?? Brings closer to Derevianko result...???
      // auto excl = &np;
      auto excl = nullptr;
      auto dAp = he1.reducedME(aA, np) + dVE1.dV_ab(aA, np, dVconj, excl);
      auto hpB = hpnc.reducedME(np, aB) + dVpnc.dV_ab(np, aB, dVconj, excl);
      auto hAp = hpnc.reducedME(aA, np) + dVpnc.dV_ab(aA, np, dVconj, excl);
      auto dpB = he1.reducedME(np, aB) + dVE1.dV_ab(np, aB, dVconj, excl);
      double pnc1 = c10 * dAp * hpB / (aB.en - np.en);
      double pnc2 = c01 * hAp * dpB / (aA.en - np.en);
      if (np.n <= main_n && print_all)
        printf("%7s, pnc= %12.5e + %12.5e = %12.5e\n", np.symbol().c_str(),
               pnc1, pnc2, pnc1 + pnc2);
      pnc += pnc1 + pnc2;
      if (coreQ)
        core += pnc1 + pnc2;
      if (np.n <= main_n && !coreQ)
        main += pnc1 + pnc2;
      if (np.n == main_n)
        main_ok = true;
      if (np.n > max_n_main)
        max_n_main = np.n;
    }
    if (print_all)
      std::cout << "...(only printing core+main)\n";

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

    auto hA = hpnc.reduced_lhs(-aA.k, aA);
    auto hB = hpnc.reduced_rhs(-aB.k, aB);

    auto e1A = he1.reduced_lhs(-aB.k, aA);
    auto e1B = he1.reduced_rhs(-aA.k, aB);

    // XXX LHS! ?
    hA -= dVpnc.dV_ab_rhs(hA, aA, false);
    hB += dVpnc.dV_ab_rhs(hB, aB, false);
    e1A += dVE1.dV_ab_rhs(e1A, aA, !dVconj);
    e1B += dVE1.dV_ab_rhs(e1B, aB, dVconj);

    auto yA_w = HartreeFock::solveMixedState(hA.k, aA, 0, v, alpha,
                                             wf.core_orbitals, hA);
    auto xB_w = HartreeFock::solveMixedState(hB.k, aB, 0, v, alpha,
                                             wf.core_orbitals, hB);

    auto omegaSE = aA.en - aB.en;
    auto xB_d = HartreeFock::solveMixedState(e1B.k, aB, omegaSE, v, alpha,
                                             wf.core_orbitals, e1B);
    auto yA_d = HartreeFock::solveMixedState(e1A.k, aA, -omegaSE, v, alpha,
                                             wf.core_orbitals, e1A);

    auto pnc1_w = c01 * he1.reducedME(yA_w, aB);
    auto pnc2_w = c10 * he1.reducedME(aA, xB_w);
    pnc1_w += c01 * dVE1.dV_ab(yA_w, aB, dVconj);
    pnc2_w += c10 * dVE1.dV_ab(aA, xB_w, dVconj);

    std::cout << "\nMixed states method: \n";
    std::cout << "<yA_w|d| B> + <A |d|xB_w> = ";
    std::cout << pnc1_w << " + " << pnc2_w << " = " << pnc1_w + pnc2_w << "\n";

    auto pnc1_d = c01 * hpnc.reducedME(aA, xB_d);
    auto pnc2_d = c10 * hpnc.reducedME(yA_d, aB);
    pnc1_d += c01 * dVpnc.dV_ab(aA, xB_d, dVconj);
    pnc2_d += c10 * dVpnc.dV_ab(yA_d, aB, dVconj);
    std::cout << "<A |h|xB_d> + <yA_d|h| B> = ";
    std::cout << pnc1_d << " + " << pnc2_d << " = " << pnc1_d + pnc2_d << "\n";

    v1 = pnc1_w + pnc2_w;
    v2 = pnc1_d + pnc2_d;
    std::cout << "eps = " << (v1 - v2) / (0.5 * (v1 + v2)) << "\n";

    // Find the core+main contributions by forcing mixed-states to be orthoganal
    // to the core/main states:
    // orthog wrt core:
    wf.orthogonaliseWrt(yA_w, wf.core_orbitals);
    wf.orthogonaliseWrt(xB_w, wf.core_orbitals);
    // Core contribution:
    auto pnc1_c =
        pnc1_w - c01 * (he1.reducedME(yA_w, aB) + dVE1.dV_ab(yA_w, aB, dVconj));
    auto pnc2_c =
        pnc2_w - c10 * (he1.reducedME(aA, xB_w) + dVE1.dV_ab(aA, xB_w, dVconj));
    // Further orthog wrt 'main' part of valence (now orthog to core+main)
    for (const auto &phiv : wf.valence_orbitals) {
      if (phiv.n > main_n)
        continue;
      if (phiv.k == yA_w.k) {
        yA_w -= (yA_w * phiv) * phiv;
      }
      if (phiv.k == xB_w.k) {
        xB_w -= (xB_w * phiv) * phiv;
      }
    }
    // Main contribution:
    auto pnc1_m =
        pnc1_w - pnc1_c -
        c10 * (he1.reducedME(yA_w, aB) + dVE1.dV_ab(yA_w, aB, dVconj));
    auto pnc2_m =
        pnc2_w - pnc2_c -
        c01 * (he1.reducedME(aA, xB_w) + dVE1.dV_ab(aA, xB_w, dVconj));

    std::cout << "\n<dA'|d| B>  +  <A |d|dB'>  (by force orthog):\n";
    printf("core :%10.7f  (=%10.7f %+10.7f)\n", pnc1_c + pnc2_c, pnc1_c,
           pnc2_c);
    if (main_ok) {
      printf("main :%10.7f  (=%10.7f %+10.7f)\n", pnc1_m + pnc2_m, pnc1_m,
             pnc2_m);
    } else {
      printf("main*:%10.7f  (=%10.7f %+10.7f)", pnc1_m + pnc2_m, pnc1_m,
             pnc2_m);
      std::cout << "    *(didnt have full main, nmax=" << max_n_main << ")\n";
    }
    printf("tail :%10.7f  (=%10.7f %+10.7f)\n",
           pnc1_w - pnc1_m - pnc1_c + pnc2_w - pnc2_m - pnc2_c,
           pnc1_w - pnc1_m - pnc1_c, pnc2_w - pnc2_m - pnc2_c);
    printf("Total:%10.7f  (=%10.7f %+10.7f)\n", v1, pnc1_w, pnc2_w);
  }
}

} // namespace Module
