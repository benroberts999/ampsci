#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

namespace UnitTest {

//******************************************************************************
//! Unit tests for correlations (second-order MBPT correlation energy/potential)
bool DiagramRPA(std::ostream &obuff) {
  bool pass = true;

  Wavefunction wf({3000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.solve_valence("6sp5d4f");
  wf.formBasis({"40spdfgh", 50, 9, 1.0e-5, 0.0, 30.0, false});

  auto sorter = [](auto x, auto y) { return x.first < y.first; };
  auto compr = [](const auto &x, const auto &y) {
    return (x.second - y.second) / y.second;
  };

  {
    auto dE1 = DiracOperator::E1(*wf.rgrid);
    auto rpa = ExternalField::DiagramRPA(&dE1, wf.basis, wf.core, "");

    { // E1, ww=0
      using sp = std::pair<std::string, double>;
      auto e1VD = std::vector<sp>{
          {"6p-6s+", -0.30326}, {"6p+6s+", 0.41334},  {"6s+6p-", -0.30326},
          {"5d-6p-", -0.33843}, {"6s+6p+", -0.41334}, {"5d-6p+", 0.14597},
          {"5d+6p+", -0.43300}, {"6p-5d-", 0.33843},  {"6p+5d-", 0.14597},
          {"4f-5d-", -0.13890}, {"6p+5d+", 0.43300},  {"4f-5d+", 0.03692},
          {"4f+5d+", -0.16510}, {"5d-4f-", 0.13890},  {"5d+4f-", 0.03692},
          {"5d+4f+", 0.16510}};
      std::sort(begin(e1VD), end(e1VD), sorter);

      rpa.clear(); // dont' need here, but symmetric
      rpa.solve_core(0.0);
      std::vector<sp> e1me;
      for (const auto &Fv : wf.valence) {
        for (const auto &Fw : wf.valence) {
          if (dE1.isZero(Fv.k, Fw.k))
            continue;
          // nb: -ve sign account for |e| vs -|e| diff from VD
          e1me.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(),
                            -rpa.dV(Fv, Fw));
        }
      }
      std::sort(begin(e1me), end(e1me), sorter);
      for (auto i = 0ul; i < e1me.size(); ++i) {
        printf("%6s %8.5f [%8.5f]\n", e1me[i].first.c_str(), e1me[i].second,
               e1VD[i].second);
      }

      auto [eps, at] = qip::compare(e1me, e1VD, compr);
      pass &= qip::check_value(&obuff, "RPA(D) E1 w=0.00 " + at->first, eps,
                               0.0, 0.0004);
    }

    { // E1, w=0.05
      using sp = std::pair<std::string, double>;
      auto e1VD = std::vector<sp>{
          {"6p-6s+", -0.30321}, {"6p+6s+", 0.41294},  {"6s+6p-", -0.30321},
          {"5d-6p-", -0.34250}, {"6s+6p+", -0.41294}, {"5d-6p+", 0.14779},
          {"5d+6p+", -0.43950}, {"6p-5d-", 0.34250},  {"6p+5d-", 0.14779},
          {"4f-5d-", -0.13940}, {"6p+5d+", 0.43950},  {"4f-5d+", 0.03703},
          {"4f+5d+", -0.16560}, {"5d-4f-", 0.13940},  {"5d+4f-", 0.03703},
          {"5d+4f+", 0.16560}};
      std::sort(begin(e1VD), end(e1VD), sorter);

      rpa.clear(); // start from scratch
      rpa.solve_core(0.05);
      std::vector<sp> e1me;
      for (const auto &Fv : wf.valence) {
        for (const auto &Fw : wf.valence) {
          if (dE1.isZero(Fv.k, Fw.k))
            continue;
          // nb: -ve sign account for |e| vs -|e| diff from VD
          e1me.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(),
                            -rpa.dV(Fv, Fw));
        }
      }
      std::sort(begin(e1me), end(e1me), sorter);
      for (auto i = 0ul; i < e1me.size(); ++i) {
        printf("%6s %8.5f [%8.5f]\n", e1me[i].first.c_str(), e1me[i].second,
               e1VD[i].second);
      }

      auto [eps, at] = qip::compare(e1me, e1VD, compr);
      pass &= qip::check_value(&obuff, "RPA(D) E1 w=0.05 " + at->first, eps,
                               0.0, 0.0005);
    }
  }

  { // HFS (compare A, not dV)
    auto h = DiracOperator::HyperfineA(1.0, 1.0, 0.0, *wf.rgrid,
                                       DiracOperator::Hyperfine::pointlike_F());
    auto rpa = ExternalField::DiagramRPA(&h, wf.basis, wf.core, "");
    using sp = std::pair<std::string, double>;
    auto e1VD = std::vector<sp>{{"6s+", 2.342288e3},  {"6p-", 2.732209e2},
                                {"6p+", 5.808505e1},  {"5d-", 0.219042e2},
                                {"5d+", -0.330845e2}, {"4f-", 0.422569e-1},
                                {"4f+", -0.194970e-1}};
    std::sort(begin(e1VD), end(e1VD), sorter);

    rpa.solve_core(0.0);
    std::vector<sp> e1me;
    for (const auto &Fv : wf.valence) {
      auto a = DiracOperator::HyperfineA::convertRMEtoA(Fv, Fv);
      e1me.emplace_back(Fv.shortSymbol(),
                        a * (h.reducedME(Fv, Fv) + rpa.dV(Fv, Fv)));
    }
    std::sort(begin(e1me), end(e1me), sorter);
    for (auto i = 0ul; i < e1me.size(); ++i) {
      printf("%3s %12.5e [%12.5e]\n", e1me[i].first.c_str(), e1me[i].second,
             e1VD[i].second);
    }

    auto [eps, at] = qip::compare(e1me, e1VD, compr);
    pass &= qip::check_value(&obuff, "RPA(D) HFS(A) " + at->first, eps, 0.0,
                             0.0007);
  }

  return pass;
}
} // namespace UnitTest
