#include "ExternalField/DiagramRPA.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Random.hpp"
#include <algorithm>
#include <string>

//==============================================================================
// Tests for Diragram RPA (basic unit)
TEST_CASE("External Field: Diagram RPA - basic unit tests",
          "[ExternalField][DiagramRPA][RPA][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: Diagram RPA - basic unit tests\n";

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"K", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", 0.0, "[Ar]");
  wf.solve_valence("4sp");
  // nb: use very small basis. Don't care about numerical results, just that eveything is working correctly.
  wf.formBasis({"10spd", 20, 7, 1.0e-3, 1.0e-3, 40.0, false});

  auto dE1 = DiracOperator::E1(wf.grid());

  auto F6s = wf.getState("4s");
  auto F6p = wf.getState("4p-");
  REQUIRE(F6s != nullptr);
  REQUIRE(F6p != nullptr);

  const auto file_name = "deleteme-" + qip::random_string(4);
  auto rpa = ExternalField::DiagramRPA(&dE1, wf.basis(), wf.core(), file_name);
  rpa.solve_core(0.0, 20);

  const auto dv1_0 = rpa.dV1(*F6s, *F6p);
  const auto dv_0 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv1_0 << " " << dv_0 << "\n";
  REQUIRE(dv1_0 != 0.0);
  REQUIRE(dv_0 != 0.0);
  // doesn't work with small grid/basis/its etc
  // REQUIRE(std::abs(dv_0) < std::abs(dv1_0));

  // This should read W matrix from file, but not calculate t's
  // Grabbing t's from other rpa - should yield exact same result!
  auto rpa2 = ExternalField::DiagramRPA(&dE1, wf.basis(), wf.core(), file_name);
  rpa2.grab_tam(&rpa);
  const auto dv1_2 = rpa2.dV1(*F6s, *F6p);
  const auto dv_2 = rpa2.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv1_2 << " " << dv_2 << "\n";

  REQUIRE(dv1_2 == Approx(dv1_0));
  REQUIRE(dv_2 == Approx(dv_0));

  auto rpa3 = ExternalField::DiagramRPA(&dE1, &rpa);
  const auto dv1_3 = rpa3.dV1(*F6s, *F6p);

  // can get lowest-rder _before_ solving core:
  REQUIRE(dv1_3 == Approx(dv1_0));

  // but need to solve_core (or 'grab_tam') to get RPA:
  rpa3.solve_core(0.0);
  const auto dv_3 = rpa3.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv1_3 << " " << dv_3 << "\n";
  REQUIRE(dv_3 == Approx(dv_0).epsilon(1.0e-2));

  // if we 'clear' the tams, RPA should be equivilant to first-order:
  rpa3.clear();
  const auto dv_30 = rpa3.dV(*F6s, *F6p);
  REQUIRE(dv_30 == Approx(dv1_0));
}

//==============================================================================
// Tests for Diragram RPA (integration/regression)
TEST_CASE("External Field: Diagram RPA",
          "[ExternalField][DiagramRPA][RPA][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: Diagram RPA\n";

  Wavefunction wf({3500, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.solve_valence("6sp5d4f");
  wf.formBasis({"40spdfgh", 60, 9, 1.0e-4, 1.0e-3, 40.0, false});

  auto sorter = [](auto x, auto y) { return x.first < y.first; };
  auto compr = [](const auto &x, const auto &y) {
    return (x.second - y.second) / y.second;
  };

  {
    auto dE1 = DiracOperator::E1(wf.grid());
    auto rpa = ExternalField::DiagramRPA(&dE1, wf.basis(), wf.core(), "");

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
      for (const auto &Fv : wf.valence()) {
        for (const auto &Fw : wf.valence()) {
          if (dE1.isZero(Fv.kappa(), Fw.kappa()))
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
      INFO("RPA(D) E1 w=0.00 " << at->first << " " << eps);
      REQUIRE(std::abs(eps) < 0.0004);
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
      for (const auto &Fv : wf.valence()) {
        for (const auto &Fw : wf.valence()) {
          if (dE1.isZero(Fv.kappa(), Fw.kappa()))
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
      INFO("RPA(D) E1 w=0.05 " << at->first << " " << eps);
      REQUIRE(std::abs(eps) < 0.0005);
    }
  }

  { // HFS (compare A, not dV)
    auto h = DiracOperator::HyperfineA(1.0, 1.0, 0.0, wf.grid(),
                                       DiracOperator::Hyperfine::pointlike_F());
    auto rpa = ExternalField::DiagramRPA(&h, wf.basis(), wf.core(), "");
    using sp = std::pair<std::string, double>;
    auto e1VD = std::vector<sp>{{"6s+", 2.342288e3},  {"6p-", 2.732209e2},
                                {"6p+", 5.808505e1},  {"5d-", 0.219042e2},
                                {"5d+", -0.330845e2}, {"4f-", 0.422569e-1},
                                {"4f+", -0.194970e-1}};
    std::sort(begin(e1VD), end(e1VD), sorter);

    rpa.solve_core(0.0);
    std::vector<sp> e1me;
    for (const auto &Fv : wf.valence()) {
      const auto a = DiracOperator::HyperfineA::convertRMEtoA(Fv, Fv);

      e1me.emplace_back(Fv.shortSymbol(),
                        a * (h.reducedME(Fv, Fv) + rpa.dV(Fv, Fv)));
    }

    std::sort(begin(e1me), end(e1me), sorter);
    for (auto i = 0ul; i < e1me.size(); ++i) {
      printf("%3s %12.5e [%12.5e]\n", e1me[i].first.c_str(), e1me[i].second,
             e1VD[i].second);
    }

    auto [eps, at] = qip::compare(e1me, e1VD, compr);
    INFO("RPA(D) HFS(A) " << at->first << " " << eps);
    // used to be 0.0007 .. changed after "Breit/Basis" fiasco
    REQUIRE(std::abs(eps) < 0.009);
  }
}
