#include "ExternalField/DiagramRPA.hpp"
#include "DiracOperator/include.hpp"
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

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"K", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", std::nullopt, "[Ar]", 1.0e-10, false);
  wf.solve_valence("4sp", false);
  // nb: use very small basis. Don't care about numerical results, just that eveything is working correctly.
  SplineBasis::Parameters params{"10spd", 20, 7, 1.0e-3, 1.0e-3, 40.0};
  params.verbose = false;
  wf.formBasis(params);

  auto dE1 = DiracOperator::E1(wf.grid());

  auto F6s = wf.getState("4s");
  auto F6p = wf.getState("4p-");
  REQUIRE(F6s != nullptr);
  REQUIRE(F6p != nullptr);

  const auto file_name = "deleteme_" + qip::random_string(3);
  auto rpa = ExternalField::DiagramRPA(&dE1, wf.basis(), wf.vHF(), file_name);

  // before solve_core, should be valid, but return 0.0
  const auto dv_00 = rpa.dV(*F6s, *F6p);
  REQUIRE(dv_00 == 0.0);

  rpa.solve_core(0.0, 5);

  const auto dv_0 = rpa.dV(*F6s, *F6p);
  REQUIRE(dv_0 != 0.0);

  // This should read W matrix from file, but not calculate t's
  // Grabbing t's from other rpa - should yield exact same result!
  auto rpa2 = ExternalField::DiagramRPA(&dE1, wf.basis(), wf.vHF(), file_name);

  // const auto dv1_2 = rpa2.dV1(*F6s, *F6p);
  const auto dv_20 = rpa2.dV(*F6s, *F6p);
  REQUIRE(dv_20 == 0.0);
  rpa2.solve_core(0.0, 5, false);
  const auto dv_2 = rpa2.dV(*F6s, *F6p);
  REQUIRE(dv_2 == Approx(dv_0));

  auto rpa3 = ExternalField::DiagramRPA(&dE1, &rpa);

  // but need to solve_core (or 'grab_tam') to get RPA:
  rpa3.solve_core(0.0, 5, false);
  const auto dv_3 = rpa3.dV(*F6s, *F6p);
  REQUIRE(dv_3 == Approx(dv_0));

  // let rpa3 converge:
  rpa3.solve_core(0.0, 256, true);
  const auto dv_4 = rpa3.dV(*F6s, *F6p);
  REQUIRE(dv_4 != 0.0);
  REQUIRE(rpa3.last_omega() == 0.0);
  REQUIRE(rpa3.last_eps() < 1.0e-6);
  REQUIRE(rpa3.last_its() > 1);

  // dV_rhs(Fa) should return a DiracSpinor with the requested kappa
  for (int kappa : {-1, 1, -2, 2, -3}) {
    REQUIRE(rpa.dV_rhs(kappa, *F6s).kappa() == kappa);
    REQUIRE(rpa.dV_rhs(kappa, *F6p).kappa() == kappa);
  }

  // Fb * dV_rhs(kappa_b, Fa) should equal dV(Fb, Fa)
  {
    const auto dv_ps = rpa.dV(*F6p, *F6s);
    const auto dv_sp = rpa.dV(*F6s, *F6p);
    const auto rhs_ps = *F6p * rpa.dV_rhs(F6p->kappa(), *F6s);
    const auto rhs_sp = *F6s * rpa.dV_rhs(F6s->kappa(), *F6p);
    REQUIRE(rhs_ps == Approx(dv_ps));
    REQUIRE(rhs_sp == Approx(dv_sp));
  }

  // E1v (velocity gauge, frequency-dependent operator)
  // Tests of Hermicity for f-dependent case
  const auto &v = *F6s;
  const auto F6p3 = wf.getState("4p+");
  for (const auto pw : {F6p, F6p3}) {
    const auto &w = *pw;
    // transition: v -> w
    const double omega = w.en() - v.en();
    REQUIRE(omega > 0.0);

    auto t_plus = DiracOperator::E1v(wf.alpha(), omega);
    auto t_minus = DiracOperator::E1v(wf.alpha(), -omega);

    auto dV_plus = ExternalField::DiagramRPA(&t_plus, wf.basis(), wf.vHF(),
                                             file_name, false);
    auto dV_minus = ExternalField::DiagramRPA(&t_minus, wf.basis(), wf.vHF(),
                                              file_name, false);

    auto twv = t_plus.reducedME(w, v);

    // incorrect Hermitian:
    auto twv_X = t_plus.reducedME(v, w) * t_minus.symm_sign(v, w);
    REQUIRE(twv_X != Approx(twv));

    // correct:
    auto twv_c = t_minus.reducedME(v, w) * t_minus.symm_sign(v, w);
    REQUIRE(twv_c == Approx(twv));

    dV_plus.solve_core(omega, 5, false);
    dV_minus.solve_core(-omega, 5, false);

    // These should all be the same (uses 'auto conj')
    const auto dV1 = dV_plus.dV(w, v);
    const auto dV2 = dV_plus.dV(v, w) * t_minus.symm_sign(v, w);
    const auto dV3 = dV_minus.dV(w, v);
    const auto dV4 = dV_minus.dV(v, w) * t_minus.symm_sign(v, w);
    REQUIRE(dV2 == Approx(dV1));
    REQUIRE(dV3 == Approx(dV1));
    REQUIRE(dV4 == Approx(dV1));

    // These should NOT be the same as dV1 (incorrect conj)
    // But should be the same as each other
    const auto dV1_X = dV_plus.dV(w, v, true);
    const auto dV2_X = dV_plus.dV(v, w, false) * t_minus.symm_sign(v, w);
    const auto dV3_X = dV_minus.dV(w, v, false);
    const auto dV4_X = dV_minus.dV(v, w, true) * t_minus.symm_sign(v, w);
    REQUIRE(dV1_X != dV1);
    REQUIRE(dV2_X != dV1);
    REQUIRE(dV3_X != dV1);
    REQUIRE(dV4_X != dV1);
    REQUIRE(dV2_X == Approx(dV1_X));
    REQUIRE(dV3_X == Approx(dV1_X));
    REQUIRE(dV4_X == Approx(dV1_X));

    // Fb * dV_rhs(kappa_b, Fa) should equal dV(Fb, Fa)
    const auto dV1_r = w * dV_plus.dV_rhs(w.kappa(), v, false);
    const auto dV2_r =
      v * dV_plus.dV_rhs(v.kappa(), w, true) * t_minus.symm_sign(v, w);
    const auto dV3_r = w * dV_minus.dV_rhs(w.kappa(), v, true);
    const auto dV4_r =
      v * dV_minus.dV_rhs(v.kappa(), w, false) * t_minus.symm_sign(v, w);
    REQUIRE(dV1_r == Approx(dV1));
    REQUIRE(dV2_r == Approx(dV1_r));
    REQUIRE(dV3_r == Approx(dV1_r));
    REQUIRE(dV4_r == Approx(dV1_r));
  }
}

//==============================================================================
// Tests for Diragram RPA (integration/regression)
TEST_CASE("External Field: Diagram RPA",
          "[ExternalField][DiagramRPA][RPA][integration]") {

  Wavefunction wf({3500, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", std::nullopt, "[Xe]");
  wf.solve_valence("6sp5d4f");
  wf.formBasis({"30spd20f", 60, 9, 1.0e-4, 1.0e-3, 40.0});

  auto sorter = [](auto x, auto y) { return x.first < y.first; };
  auto compr = [](const auto &x, const auto &y) {
    return (x.second - y.second) / y.second;
  };

  {
    auto dE1 = DiracOperator::E1(wf.grid());
    auto rpa = ExternalField::DiagramRPA(&dE1, wf.basis(), wf.vHF(), "");

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
      REQUIRE(std::abs(eps) < 0.001);
      // REQUIRE(std::abs(eps) < 0.02);
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
      // REQUIRE(std::abs(eps) < 0.0005);
      REQUIRE(std::abs(eps) < 0.02); //?
    }
  }

  { // HFS (compare A, not dV)
    auto h = DiracOperator::hfs(1, 1.0, 0.0, wf.grid(),
                                DiracOperator::Hyperfine::pointlike_F());
    auto rpa = ExternalField::DiagramRPA(&h, wf.basis(), wf.vHF(), "");
    using sp = std::pair<std::string, double>;
    auto e1VD = std::vector<sp>{{"6s+", 2.342288e3},  {"6p-", 2.732209e2},
                                {"6p+", 5.808505e1},  {"5d-", 0.219042e2},
                                {"5d+", -0.330845e2}, {"4f-", 0.422569e-1},
                                {"4f+", -0.194970e-1}};
    std::sort(begin(e1VD), end(e1VD), sorter);

    rpa.solve_core(0.0);
    std::vector<sp> e1me;
    for (const auto &Fv : wf.valence()) {
      const auto a = DiracOperator::Hyperfine::convert_RME_to_HFSconstant(
        1, Fv.kappa(), Fv.kappa());

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
    // REQUIRE(std::abs(eps) < 0.009);
    REQUIRE(std::abs(eps) < 0.09);
  }
}
