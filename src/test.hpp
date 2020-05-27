#include "Angular/Angular_369j.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracOperator/Operators.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <string>

void runRPA(const Wavefunction &wf) {

  std::cout << "\n";

  const auto &Fv = wf.valence[0];
  const auto &Fw = wf.valence[1];
  auto omega = 0.00; // fabs(Fv.en - Fw.en);

  DiracOperator::E1 he1(wf.rgrid);
  // DiracOperator::Hyperfine he1(1.0, 1.0, 0.0, wf.rgrid);
  const auto rank = he1.rank();

  std::vector<DiracSpinor> holes;
  std::vector<DiracSpinor> excited; // = wf.core;
  for (const auto &Fi : wf.basis) {
    if (wf.isInCore(Fi.n, Fi.k)) {
      holes.push_back(Fi);
    } else {
      excited.push_back(Fi);
    }
  }

  std::vector<std::vector<double>> t0_am;
  std::vector<std::vector<double>> t0_ma;
  for (const auto &Fa : holes) {
    std::vector<double> t0a_m;
    for (const auto &Fm : excited) {
      t0a_m.push_back(he1.reducedME(Fa, Fm));
    }
    t0_am.push_back(t0a_m);
  }
  for (const auto &Fm : excited) {
    std::vector<double> t0m_a;
    for (const auto &Fa : holes) {
      t0m_a.push_back(he1.reducedME(Fm, Fa));
    }
    t0_ma.push_back(t0m_a);
  }
  auto t_am = t0_am;
  auto t_ma = t0_ma;

  // RPA: store Z Coulomb integrals (used only for Core RPA its)
  std::cout << "Filling RPA MBPT Matrix:\n" << std::flush;
  std::vector<std::vector<std::vector<std::vector<double>>>> Zanmb;
  std::vector<std::vector<std::vector<std::vector<double>>>> Zabmn;
  Zanmb.resize(holes.size());
  Zabmn.resize(holes.size());
#pragma omp parallel for
  for (std::size_t i = 0; i < holes.size(); i++) {
    const auto &Fa = holes[i];
    std::vector<std::vector<std::vector<double>>> Za_nmb;
    std::vector<std::vector<std::vector<double>>> Za_bmn;
    Za_nmb.reserve(excited.size());
    Za_bmn.reserve(excited.size());
    for (const auto &Fn : excited) {
      std::vector<std::vector<double>> Zan_mb;
      std::vector<std::vector<double>> Zab_mn;
      Zan_mb.reserve(excited.size());
      Zab_mn.reserve(excited.size());
      for (const auto &Fm : excited) {
        std::vector<double> Zanm_b;
        std::vector<double> Zabm_n;
        Zanm_b.reserve(holes.size());
        Zabm_n.reserve(holes.size());
        for (const auto &Fb : holes) {
          auto x = Coulomb::Zk_abcd(Fa, Fn, Fm, Fb, rank);
          auto y = Coulomb::Zk_abcd(Fa, Fb, Fm, Fn, rank);
          Zanm_b.push_back(x);
          Zabm_n.push_back(y);
        }
        Zan_mb.push_back(Zanm_b);
        Zab_mn.push_back(Zabm_n);
      }
      Za_nmb.push_back(Zan_mb);
      Za_bmn.push_back(Zab_mn);
    }
    Zanmb[i] = Za_nmb;
    Zabmn[i] = Za_bmn;
  }
  std::vector<std::vector<std::vector<std::vector<double>>>> Zmnab;
  std::vector<std::vector<std::vector<std::vector<double>>>> Zmban;
  Zmnab.resize(excited.size());
  Zmban.resize(excited.size());
#pragma omp parallel for
  for (std::size_t i = 0; i < excited.size(); i++) {
    const auto &Fm = excited[i];
    std::vector<std::vector<std::vector<double>>> Za_nmb;
    std::vector<std::vector<std::vector<double>>> Za_bmn;
    for (const auto &Fn : excited) {
      std::vector<std::vector<double>> Zan_mb;
      std::vector<std::vector<double>> Zab_mn;
      for (const auto &Fa : holes) {
        std::vector<double> Zanm_b;
        std::vector<double> Zabm_n;
        for (const auto &Fb : holes) {
          auto x = Coulomb::Zk_abcd(Fm, Fn, Fa, Fb, rank);
          auto y = Coulomb::Zk_abcd(Fm, Fb, Fa, Fn, rank);
          Zanm_b.push_back(x);
          Zabm_n.push_back(y);
        }
        Zan_mb.push_back(Zanm_b);
        Zab_mn.push_back(Zabm_n);
      }
      Za_nmb.push_back(Zan_mb);
      Za_bmn.push_back(Zab_mn);
    }
    Zmnab[i] = Za_nmb;
    Zmban[i] = Za_bmn;
  }

  double tvw_0 = he1.reducedME(Fv, Fw);

  std::cout << "RPA itterations for core:\n" << std::flush;
  const int max_its = 99;
  const double eps_targ = 1.0e-6;
  for (int i = 1; i <= max_its; i++) {
    double max = 0.0;
    for (std::size_t ia = 0; ia < holes.size(); ia++) {
      const auto &Fa = holes[ia];
      for (std::size_t im = 0; im < excited.size(); im++) {
        const auto &Fm = excited[im];

        double sum_am = 0;
        double sum_ma = 0;
        auto f = (1.0 / (2 * rank + 1));
        std::size_t ib = 0;
        for (const auto &Fb : holes) {
          std::size_t in = 0;
          for (const auto &Fn : excited) {

            auto s1 = ((abs(Fb.twoj() - Fn.twoj()) + 2) % 4 == 0) ? 1 : -1;

            auto zanmb = Zanmb[ia][in][im][ib];
            auto zabmn = Zabmn[ia][in][im][ib];
            auto zmnab = Zmnab[im][in][ia][ib];
            auto zmban = Zmban[im][in][ia][ib];

            auto t_bn = t_am[ib][in];
            auto t_nb = t_ma[in][ib];
            auto A = t_bn * zanmb / (Fb.en - Fn.en - omega);
            auto B = t_nb * zabmn / (Fb.en - Fn.en + omega);
            auto C = t_bn * zmnab / (Fb.en - Fn.en - omega);
            auto D = t_nb * zmban / (Fb.en - Fn.en + omega);
            sum_am += s1 * (A + B);
            sum_ma += s1 * (C + D);

            ++in;
          }
          ++ib;
        }
        auto prev = t_am[ia][im];
        t_am[ia][im] = 0.5 * (t_am[ia][im] + t0_am[ia][im] + f * sum_am);
        t_ma[im][ia] = 0.5 * (t_ma[im][ia] + t0_ma[im][ia] + f * sum_ma);

        auto delta = std::abs((t_am[ia][im] - prev) / t_am[ia][im]);
        if (delta > max)
          max = fabs(delta);
      }
    }
    if ((i % 10 == 0 && i > 1) || max < eps_targ || i == max_its) {
      std::cout << "it=" << i << ", eps=" << max << "\n";
    }
    if (max < eps_targ)
      break;
  }

  // For transition:
  double sum = 0.0;
  double sum0 = 0.0;
  auto f = (1.0 / (2 * rank + 1));
  std::size_t ia = 0;
  for (const auto &Fa : holes) {
    std::size_t im = 0;
    for (const auto &Fm : excited) {

      auto s1 = ((abs(Fa.twoj() - Fm.twoj()) + 2) % 4 == 0) ? 1 : -1;

      auto Zwmva = Coulomb::Zk_abcd(Fw, Fm, Fv, Fa, rank);
      auto Zwavm = Coulomb::Zk_abcd(Fw, Fa, Fv, Fm, rank);

      auto tt_am = t_am[ia][im];
      auto tt_ma = t_ma[im][ia];
      auto tt0_am = t0_am[ia][im];
      auto tt0_ma = t0_ma[im][ia];

      auto A = tt_am * Zwmva / (Fa.en - Fm.en - omega);
      auto B = Zwavm * tt_ma / (Fa.en - Fm.en + omega);
      auto A0 = tt0_am * Zwmva / (Fa.en - Fm.en - omega);
      auto B0 = Zwavm * tt0_ma / (Fa.en - Fm.en + omega);

      sum += s1 * (A + B);
      sum0 += s1 * (A0 + B0);
      ++im;
    }
    ++ia;
  }
  std::cout << "\n";
  auto tvw_rpa = tvw_0 + f * sum;
  auto tvw_rpa1 = tvw_0 + f * sum0;
  std::cout << tvw_0 << " " << tvw_rpa1 << " " << tvw_rpa << "\n";
  // auto a = DiracOperator::Hyperfine::convertRMEtoA(Fw, Fv);
  // std::cout << tvw_0 * a << " " << tvw_rpa1 * a << " " << tvw_rpa * a <<
  // "\n";
}
