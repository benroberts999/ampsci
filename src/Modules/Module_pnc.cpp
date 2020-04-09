#include "Modules/Module_pnc.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "HF/MixedStates.hpp"
#include "IO/UserInput.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <tuple>

namespace Module {

//******************************************************************************
void polarisability(const IO::UserInputBlock &input, const Wavefunction &wf) {
  std::cout << "\n NOTE: In developement! Not finished!\n";

  const auto alpha = wf.alpha;

  input.checkBlock({"a", "b", "omega"});
  const auto a_str = input.get<std::string>("a", "");
  const auto w = input.get("omega", 0.00);
  // const auto b_str = input.get<std::string>("a", a_str);

  // [[maybe_unused]] // don't use e
  const auto [na, ka, e] = AtomData::listOfStates_singlen(a_str).front();
  const auto pFa = wf.getState(na, ka);

  const auto he1 = DiracOperator::E1(wf.rgrid);
  auto dVE1 = HF::ExternalField(&he1, wf.core, wf.get_Vlocal(), alpha);

  dVE1.solve_TDHFcore(w, 1); // 1 it means no RPA (for TDHF version)

  // auto twom = std::min(pFa->twoj(), pFa->twoj());
  // auto c10 = he1.rme3js(pFa->twoj(), pFa->twoj(), twom);

  auto fudge_factor = (4.0 / 3.0); // ???

  // core contribution:
  auto alpha_core = 0.0;
  for (const auto &Fb : wf.core) {
    const auto &Xb = dVE1.get_dPsis(Fb, HF::dPsiType::X);
    const auto &Yb = dVE1.get_dPsis(Fb, HF::dPsiType::Y);
    // const auto Xb = dVE1.solve_dPsis(Fb, w, HF::dPsiType::X);
    // const auto Yb = dVE1.solve_dPsis(Fb, w, HF::dPsiType::X);

    for (const auto &Xbeta : Xb) {
      alpha_core += he1.reducedME(Fb, Xbeta);
    }
    for (const auto &Ybeta : Yb) {
      alpha_core += he1.reducedME(Fb, Ybeta);
      // alpha_core += he1.reducedME(Ybeta, Fb);
    }
  }
  alpha_core *= (-1.0 / 3.0) * fudge_factor;
  std::cout << alpha_core << "\n";

  // Valence part:
  double alpha_valence = 0.0;
  if (pFa) {
    // omega?
    const auto Xb = dVE1.solve_dPsis(*pFa, w, HF::dPsiType::X, wf.getSigma());
    const auto Yb = dVE1.solve_dPsis(*pFa, w, HF::dPsiType::Y, wf.getSigma());
    for (const auto &Xbeta : Xb) {
      alpha_valence += he1.reducedME(*pFa, Xbeta);
    }
    for (const auto &Ybeta : Yb) {
      alpha_valence += he1.reducedME(Ybeta, *pFa);
    }
    alpha_valence *= (-1.0 / 3.0) * fudge_factor / (pFa->twoj() + 1);
  }
  std::cout << alpha_valence << "\n";

  std::cout << "\n";

  double alpha_2 = 0.0;
  for (const auto &Fb : wf.core) {
    for (const auto &Fn : wf.spectrum) {
      if (he1.isZero(Fb.k, Fn.k))
        continue;
      const auto d1 = he1.reducedME(Fb, Fn);
      const auto d2 = he1.reducedME(Fn, Fb) + dVE1.dV_ab(Fn, Fb);
      auto de = Fb.en - Fn.en;
      alpha_2 += fudge_factor * (-2.0 / 3.0) * d1 * d2 * de / (de * de - w * w);
    }
  }
  std::cout << alpha_2 << "\n";

  double alpha_3 = 0.0;
  if (pFa) {
    for (const auto &Fn : wf.spectrum) {
      if (he1.isZero(pFa->k, Fn.k))
        continue;
      const auto d1 = he1.reducedME(*pFa, Fn);
      const auto d2 = he1.reducedME(Fn, *pFa) + dVE1.dV_ab(Fn, *pFa);
      auto de = pFa->en - Fn.en;
      alpha_3 += fudge_factor * (-2.0 / 3.0) * d1 * d2 * de /
                 (de * de - w * w) / (pFa->twoj() + 1);
    }
  }
  std::cout << alpha_3 << "\n";
}

//******************************************************************************
void calculatePNC(const IO::UserInputBlock &input, const Wavefunction &wf) {
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

  auto rpaQ = input.get("rpa", true);

  DiracOperator::PNCnsi hpnc(c, t, wf.rgrid, -wf.Nnuc());
  DiracOperator::E1 he1(wf.rgrid);
  auto alpha = wf.alpha;

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

  auto dVE1 = HF::ExternalField(&he1, wf.core, wf.get_Vlocal(), alpha);
  auto dVpnc =
      HF::ExternalField(&hpnc, wf.core, wf.get_Vlocal(), alpha);
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
    std::cout << "Sum-over-states method (HF core+valence):\n";
    std::cout << " <A|d|n><n|hw|B>/dEB + <A|hw|n><n|d|B>/dEA\n";
    double pnc = 0, core = 0, main = 0;
    for (int i = 0; i < 2; i++) {
      auto &tmp_orbs = (i == 0) ? wf.core : wf.valence;
      for (auto &np : tmp_orbs) {
        // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
        if (np == aB || np == aA)
          continue;
        if (hpnc.isZero(np.k, aA.k) && hpnc.isZero(np.k, aB.k))
          continue;
        auto coreQ = wf.isInCore(np.n, np.k);
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

  if (!wf.spectrum.empty()) {
    std::cout << "\nSum-over-states method (Spectrum):\n";
    std::cout << " <A|d|n><n|hw|B>/dEB + <A|hw|n><n|d|B>/dEA\n";
    double pnc = 0, core = 0, main = 0;
    for (auto &np : wf.spectrum) {
      if (np == aB || np == aA)
        continue;
      if (hpnc.isZero(np.k, aA.k) && hpnc.isZero(np.k, aB.k))
        continue;
      auto coreQ = wf.isInCore(np.n, np.k);
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
    auto v1 = 0.0, v2 = 0.0;

    //<yA_w|d| B> + <A |d|xB_w>
    auto yA_w = dVpnc.solve_dPsi(aA, 0, HF::dPsiType::Y, -aA.k, wf.getSigma());
    auto xB_w = dVpnc.solve_dPsi(aB, 0, HF::dPsiType::X, -aB.k, wf.getSigma());
    const auto pnc1_w =
        c01 * (he1.reducedME(yA_w, aB) + dVE1.dV_ab(yA_w, aB, dVconj));
    const auto pnc2_w =
        c10 * (he1.reducedME(aA, xB_w) + dVE1.dV_ab(aA, xB_w, dVconj));

    //<A |h|xB_d> + <yA_d|h| B>
    const auto omegaSE = (aA.en - aB.en);
    auto xB_d =
        dVE1.solve_dPsi(aB, omegaSE, HF::dPsiType::X, -aA.k, wf.getSigma());
    auto yA_d =
        dVE1.solve_dPsi(aA, omegaSE, HF::dPsiType::Y, -aB.k, wf.getSigma());
    const auto pnc1_d =
        c01 * (hpnc.reducedME(aA, xB_d) + dVpnc.dV_ab(aA, xB_d, dVconj));
    const auto pnc2_d =
        c10 * (hpnc.reducedME(yA_d, aB) + dVpnc.dV_ab(yA_d, aB, dVconj));

    std::cout << "\nMixed states method: \n";
    std::cout << "<yA_w|d| B> + <A |d|xB_w> = ";
    std::cout << pnc1_w << " + " << pnc2_w << " = " << pnc1_w + pnc2_w << "\n";

    std::cout << "<A |h|xB_d> + <yA_d|h| B> = ";
    std::cout << pnc1_d << " + " << pnc2_d << " = " << pnc1_d + pnc2_d << "\n";

    v1 = pnc1_w + pnc2_w;
    v2 = pnc1_d + pnc2_d;
    std::cout << "eps = " << (v1 - v2) / (0.5 * (v1 + v2)) << "\n";

    // Find the core+main contributions by forcing mixed-states to be
    // orthoganal to the core/main states: orthog wrt core:
    wf.orthogonaliseWrt(yA_w, wf.core);
    wf.orthogonaliseWrt(xB_w, wf.core);
    // Core contribution:
    auto pnc1_c =
        pnc1_w - c01 * (he1.reducedME(yA_w, aB) + dVE1.dV_ab(yA_w, aB, dVconj));
    auto pnc2_c =
        pnc2_w - c10 * (he1.reducedME(aA, xB_w) + dVE1.dV_ab(aA, xB_w, dVconj));
    // Further orthog wrt 'main' part of valence (now orthog to core+main)
    for (const auto &phiv : wf.valence) {
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
