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

  DiracOperator::PNCnsi hpnc(c, t, *(wf.rgrid), -wf.Nnuc());
  DiracOperator::E1 he1(*(wf.rgrid));

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

  // Solve TDHF:
  auto dVE1 = HF::ExternalField(&he1, wf.getHF());
  auto dVpnc = HF::ExternalField(&hpnc, wf.getHF());
  if (rpaQ) {
    auto omega_dflt = std::abs(aA.en - aB.en);
    auto omega = input.get("omega", omega_dflt);
    dVE1.solve_TDHFcore(omega);
    dVpnc.solve_TDHFcore(0.0);
  }
  const bool dVconj = aA.en < aB.en ? true : false;

  // Angular factors (RME -> z-comp, z=min(ja,jb))
  auto tja = aA.twoj();
  auto tjb = aB.twoj();
  auto twom = std::min(tja, tjb);
  auto c10 = he1.rme3js(tja, tjb, twom) * hpnc.rme3js(tjb, tjb, twom);
  auto c01 = hpnc.rme3js(tja, tja, twom) * he1.rme3js(tja, tjb, twom);

  // Sum over states, use HF orbitals
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
        auto dAp = he1.reducedME(aA, np) + dVE1.dV(aA, np, dVconj);
        auto hpB = hpnc.reducedME(np, aB) + dVpnc.dV(np, aB, dVconj);
        auto hAp = hpnc.reducedME(aA, np) + dVpnc.dV(aA, np, dVconj);
        auto dpB = he1.reducedME(np, aB) + dVE1.dV(np, aB, dVconj);
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
      auto dAp = he1.reducedME(aA, np) + dVE1.dV(aA, np, dVconj, excl);
      auto hpB = hpnc.reducedME(np, aB) + dVpnc.dV(np, aB, dVconj, excl);
      auto hAp = hpnc.reducedME(aA, np) + dVpnc.dV(aA, np, dVconj, excl);
      auto dpB = he1.reducedME(np, aB) + dVE1.dV(np, aB, dVconj, excl);
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
        c01 * (he1.reducedME(yA_w, aB) + dVE1.dV(yA_w, aB, dVconj));
    const auto pnc2_w =
        c10 * (he1.reducedME(aA, xB_w) + dVE1.dV(aA, xB_w, dVconj));

    //<A |h|xB_d> + <yA_d|h| B>
    const auto omegaSE = (aA.en - aB.en);
    auto xB_d =
        dVE1.solve_dPsi(aB, omegaSE, HF::dPsiType::X, -aA.k, wf.getSigma());
    auto yA_d =
        dVE1.solve_dPsi(aA, omegaSE, HF::dPsiType::Y, -aB.k, wf.getSigma());
    const auto pnc1_d =
        c01 * (hpnc.reducedME(aA, xB_d) + dVpnc.dV(aA, xB_d, dVconj));
    const auto pnc2_d =
        c10 * (hpnc.reducedME(yA_d, aB) + dVpnc.dV(yA_d, aB, dVconj));

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
        pnc1_w - c01 * (he1.reducedME(yA_w, aB) + dVE1.dV(yA_w, aB, dVconj));
    auto pnc2_c =
        pnc2_w - c10 * (he1.reducedME(aA, xB_w) + dVE1.dV(aA, xB_w, dVconj));
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
    auto pnc1_m = pnc1_w - pnc1_c -
                  c10 * (he1.reducedME(yA_w, aB) + dVE1.dV(yA_w, aB, dVconj));
    auto pnc2_m = pnc2_w - pnc2_c -
                  c01 * (he1.reducedME(aA, xB_w) + dVE1.dV(aA, xB_w, dVconj));

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

//******************************************************************************
//******************************************************************************
//    Polarisability:
//******************************************************************************

inline static std::pair<double, double>
alpha_w_tdhf(const Wavefunction &wf, const DiracSpinor *const pFa,
             const DiracOperator::E1 &he1, HF::ExternalField &dVE1,
             double omega) {

  if (dVE1.get_eps() > 1.0e-5)
    dVE1.clear_dPsi();
  dVE1.solve_TDHFcore(omega, 99, false);

  double alpha = 0.0;
  const auto f1 = (-1.0 / 3.0);
  for (const auto &Fb : wf.core) {
    const auto Xb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::X);
    const auto Yb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::Y);
    for (const auto &Xbeta : Xb) {
      alpha += f1 * he1.reducedME(Xbeta, Fb);
    }
    for (const auto &Ybeta : Yb) {
      alpha += f1 * he1.reducedME(Fb, Ybeta);
    }
  }

  if (pFa) {
    const auto f = (-1.0 / 3.0) / (pFa->twoj() + 1);
    const auto Xb =
        dVE1.solve_dPsis(*pFa, omega, HF::dPsiType::X, wf.getSigma());
    const auto Yb =
        dVE1.solve_dPsis(*pFa, omega, HF::dPsiType::Y, wf.getSigma());
    for (const auto &Xbeta : Xb) {
      alpha += f * he1.reducedME(Xbeta, *pFa);
    }
    for (const auto &Ybeta : Yb) {
      alpha += f * he1.reducedME(*pFa, Ybeta);
    }
  }
  // alpha*eps is a very rough indication of num. error
  return {alpha, alpha * dVE1.get_eps()};
}

//------------------------------------------------------------------------------
static inline std::pair<double, double>
alpha_w_sos(const Wavefunction &wf, const DiracSpinor *const pFa,
            const DiracOperator::E1 &he1, HF::ExternalField &dVE1,
            double omega) {

  const auto &basis = wf.spectrum.empty() ? wf.basis : wf.spectrum;
  if (basis.empty())
    return {0, 0};

  if (dVE1.get_eps() > 1.0e-5)
    dVE1.clear_dPsi();
  dVE1.solve_TDHFcore(omega, 99, false);

  double alpha = 0.0;
  // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
  for (const auto &Fb : wf.core) {
    for (const auto &Fn : basis) {
      if (wf.isInCore(Fn.n, Fn.k))
        continue; // if core(HF) = core(basis), these cancel excactly
      if (he1.isZero(Fb.k, Fn.k))
        continue;
      const auto d1 = he1.reducedME(Fn, Fb);
      const auto d2 = he1.reducedME(Fn, Fb) + dVE1.dV(Fn, Fb);
      const auto de = Fb.en - Fn.en;
      const auto term =
          (-2.0 / 3.0) * std::abs(d1 * d2) * de / (de * de - omega * omega);
      alpha += term;
    }
  }

  if (pFa) {
    // valence part (and core-valence)
    for (const auto &Fn : basis) {

      if (he1.isZero(pFa->k, Fn.k))
        continue;
      const auto d1 = he1.reducedME(Fn, *pFa);
      const auto d2 = he1.reducedME(Fn, *pFa) + dVE1.dV(Fn, *pFa);
      const auto de = pFa->en - Fn.en;
      const auto term = (-2.0 / 3.0) * std::abs(d1 * d2) * de /
                        (de * de - omega * omega) / (pFa->twoj() + 1);
      alpha += term;
    }
  }

  return {alpha, alpha * dVE1.get_eps()};
}

//******************************************************************************
void polarisability(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"a", "rpa", "omega_max", "omega_steps"});

  const auto ank = input.get_list("a", std::vector<int>{0, 0});
  // const auto bnk = input.get_list("b", ank);
  const auto pFa = (ank.size() == 2) ? wf.getState(ank[0], ank[1]) : nullptr;
  const auto omega = 0.0;

  {
    std::cout << "\nDipole polarisability, alpha: " << wf.atom() << " ";
    if (pFa != nullptr)
      std::cout << pFa->symbol() << "\n";
    else
      std::cout << wf.coreConfiguration_nice() << "\n";
  }
  std::cout << "Omega = " << omega << "\n";

  const auto he1 = DiracOperator::E1(*(wf.rgrid));
  auto dVE1 = HF::ExternalField(&he1, wf.getHF());

  const auto rpaQ = input.get("rpa", true);
  if (rpaQ)
    dVE1.solve_TDHFcore(omega);

  std::size_t num_pws = 4; // s, p-, p+, rest:
  auto get_index = [&](const auto &Fb) {
    return std::size_t(Fb.k_index()) < num_pws ? std::size_t(Fb.k_index())
                                               : num_pws - 1;
  };

  std::cout << "\nTDHF method:\n";
  std::cout << "Inter:    s1/2    p1/2    p3/2   other |  sum\n";
  // core contribution:
  std::vector<double> core_pws(num_pws);
  auto alpha_core = 0.0;
  {
    const auto f = (-1.0 / 3.0);
    for (const auto &Fb : wf.core) {
      const auto Xb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::X);
      const auto Yb = dVE1.solve_dPsis(Fb, omega, HF::dPsiType::Y);
      for (const auto &Xbeta : Xb) {
        core_pws[get_index(Xbeta)] += f * he1.reducedME(Xbeta, Fb);
        alpha_core += f * he1.reducedME(Xbeta, Fb);
      }
      for (const auto &Ybeta : Yb) {
        core_pws[get_index(Ybeta)] += f * he1.reducedME(Fb, Ybeta);
        alpha_core += f * he1.reducedME(Fb, Ybeta);
      }
    }
    std::cout << "core : ";
    for (auto &pw : core_pws) {
      printf("%+7.1f ", pw);
    }
    printf("= %9.4f\n", alpha_core);
  }

  // // Valence part:
  double alpha_valence = 0.0;
  std::vector<double> val_pws(num_pws);
  if (pFa) {
    const auto f = (-1.0 / 3.0) / (pFa->twoj() + 1);
    const auto Xb =
        dVE1.solve_dPsis(*pFa, omega, HF::dPsiType::X, wf.getSigma());
    const auto Yb =
        dVE1.solve_dPsis(*pFa, omega, HF::dPsiType::Y, wf.getSigma());
    for (const auto &Xbeta : Xb) {
      val_pws[get_index(Xbeta)] += f * he1.reducedME(Xbeta, *pFa);
      alpha_valence += f * he1.reducedME(Xbeta, *pFa);
    }
    for (const auto &Ybeta : Yb) {
      val_pws[get_index(Ybeta)] += f * he1.reducedME(*pFa, Ybeta);
      alpha_valence += f * he1.reducedME(*pFa, Ybeta);
    }
    std::cout << "val  : ";
    for (auto &pw : val_pws) {
      printf("%+7.1f ", pw);
    }
    printf("= %8.3f\n", alpha_valence);
  }
  std::cout << "total: ";
  for (auto i = 0ul; i < core_pws.size(); ++i) {
    printf("%+7.1f ", core_pws[i] + val_pws[i]);
  }
  printf("= %8.3f\n", alpha_core + alpha_valence);

  const auto &basis = wf.spectrum.empty() ? wf.basis : wf.spectrum;
  if (!basis.empty()) {
    std::cout << "\nSum-over-states method: ";
    if (&basis == &wf.basis)
      std::cout << "(w/ basis)\n";
    else
      std::cout << "(w/ spectrum)\n";
    std::cout << "Inter:    s1/2    p1/2    p3/2   other |  sum\n";

    double alpha_2 = 0.0;
    std::vector<double> core_pws2(num_pws);
    {
      // core part: sum_{n,c} |<n|d|c>|^2, n excited states, c core states
      for (const auto &Fb : wf.core) {
        for (const auto &Fn : basis) {
          if (wf.isInCore(Fn.n, Fn.k))
            continue; // if core(HF) = core(basis), these cancel excactly
          if (he1.isZero(Fb.k, Fn.k))
            continue;
          const auto d1 = he1.reducedME(Fn, Fb);
          const auto d2 = he1.reducedME(Fn, Fb) + dVE1.dV(Fn, Fb);
          const auto de = Fb.en - Fn.en;
          const auto term =
              (-2.0 / 3.0) * std::abs(d1 * d2) * de / (de * de - omega * omega);
          alpha_2 += term;
          core_pws2[get_index(Fn)] += term;
        }
      }
      std::cout << "core : ";
      for (auto &pw : core_pws2) {
        printf("%+7.1f ", pw);
      }
      printf("= %9.4f\n", alpha_2);
    }

    double alpha_3 = 0.0;
    double alpha_cv = 0.0;
    std::vector<double> val_pws2(num_pws);
    std::vector<double> val_pws_cv(num_pws);
    if (pFa) {
      // valence part: sum_n |<n|d|A>|^2, n excited states
      for (const auto &Fn : basis) {
        if (wf.isInCore(Fn.n, Fn.k))
          continue; // just excited part
        if (he1.isZero(pFa->k, Fn.k))
          continue;
        const auto d1 = he1.reducedME(Fn, *pFa);
        const auto d2 = he1.reducedME(Fn, *pFa) + dVE1.dV(Fn, *pFa);
        const auto de = pFa->en - Fn.en;
        const auto term = (-2.0 / 3.0) * std::abs(d1 * d2) * de /
                          (de * de - omega * omega) / (pFa->twoj() + 1);
        alpha_3 += term;
        val_pws2[get_index(Fn)] += term;
      }
      std::cout << "val  : ";
      for (auto &pw : val_pws2) {
        printf("%+7.1f ", pw);
      }
      printf("= %8.3f\n", alpha_3);

      // core-valence part: sum_c |<c|d|A>|^2, c core states
      for (const auto &Fn : basis) {
        if (!wf.isInCore(Fn.n, Fn.k))
          continue; // just core part (compensates for val in basis from core)
        if (he1.isZero(pFa->k, Fn.k))
          continue;
        const auto d1 = he1.reducedME(Fn, *pFa);
        const auto d2 = he1.reducedME(Fn, *pFa) + dVE1.dV(Fn, *pFa);
        const auto de = pFa->en - Fn.en;
        const auto term = (-2.0 / 3.0) * std::abs(d1 * d2) * de /
                          (de * de - omega * omega) / (pFa->twoj() + 1);
        alpha_cv += term;
        val_pws_cv[get_index(Fn)] += term;
      }
      std::cout << "c-v  : ";
      for (auto &pw : val_pws_cv) {
        printf("%+7.1f ", pw);
      }
      printf("= %8.3f\n", alpha_cv);
    }
    std::cout << "total: ";
    for (auto i = 0ul; i < core_pws.size(); ++i) {
      printf("%+7.1f ", core_pws2[i] + val_pws2[i] + val_pws_cv[i]);
    }
    printf("= %8.3f\n", alpha_2 + alpha_3 + alpha_cv);
    const auto eps =
        1.0 - ((alpha_2 + alpha_3 + alpha_cv) / (alpha_core + alpha_valence));
    printf("eps=%.1e\n", eps);
  }

  // For fn of omega
  const auto omega_max = input.get("omega_max", 0.0);
  const auto omega_steps = input.get("omega_steps", 30);

  if (omega_max > 0) {
    std::cout << "\nDynamic polarisability:\n";
    // std::cout << "ww    a(w)[TDHF] error  a(w)[SOS]  error\n";
    std::cout << "ww      a(w)          error\n";
    const auto dw = omega_max / omega_steps;
    for (auto ww = 0.0; ww < omega_max; ww += dw) {
      // auto [a2, da2] = alpha_w_sos(wf, pFa, he1, dVE1, ww);
      auto [a1, da1] = alpha_w_tdhf(wf, pFa, he1, dVE1, ww);
      // printf("%4.2f  %8.4e %.0e  %8.4e %.0e\n", ww, a1, da1, a2, da2);
      printf("%6.4f  %11.6e  %.0e\n", ww, a1, da1);
    }
  }
}

} // namespace Module
