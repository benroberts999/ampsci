#include "TDHF_DCP.hpp"
#include "Amplitudes/DCP.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

//==============================================================================
TEST_CASE("TDHF_DCP: unit tests", "[ExternalField][TDHF_DCP][DCP][unit]") {

  // nb: small grid, light atom. Don't care about numerical results, just
  // that everything is consistent.
  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);
  wf.solve_valence("4sp", false);

  const DiracOperator::E1 d(wf.grid());
  const DiracOperator::Ek q(wf.grid(), 2);
  const DiracOperator::PNCnsi hpnc(5.67073, 2.3, wf.grid());

  // Pair of operators (t, s), coupled to rank K
  struct OperatorPair {
    const DiracOperator::TensorOperator *t;
    const DiracOperator::TensorOperator *s;
    int K;
  };

  //----------------------------------------------------------------------------
  SECTION("parts add up and -- is the conjugate of ++") {

    // Between them: both signs of (-1)^{K+kt+ks}, and an imaginary operator
    const std::vector<OperatorPair> pairs{
      {&d, &d, 0}, {&d, &d, 1}, {&d, &q, 2}, {&d, &hpnc, 1}};

    // Both frequencies non-zero, and different: X and Y all differ
    const auto omega_t = 0.05;
    const auto omega_s = 0.02;

    for (const auto &[t, s, K] : pairs) {

      ExternalField::TDHF T(t, wf.vHF());
      ExternalField::TDHF S(s, wf.vHF());
      T.solve_core(omega_t, 30, false);
      S.solve_core(omega_s, 30, false);

      ExternalField::TDHF_DCP d2V(&T, &S, K, wf.vHF());
      d2V.solve_core(omega_t + omega_s, 30, false);
      REQUIRE(d2V.last_eps() < 1.0e-10);

      const auto iota_ts = d2V.imagQ() ? -1.0 : 1.0;

      int count = 0;
      for (const auto &w : wf.valence()) {
        for (const auto &v : wf.valence()) {
          if (w.parity() * v.parity() != d2V.parity() ||
              Angular::triangle(w.twoj(), v.twoj(), 2 * K) == 0)
            continue;
          ++count;

          // dV is the sum of the three parts (both components)
          for (const bool conj : {false, true}) {
            const auto d2V_11 = w * d2V.d2V_11_rhs(w.kappa(), v, conj);
            const auto d2V_Norm = w * d2V.d2V_Norm_rhs(w.kappa(), v);
            const auto d2V_2 = w * d2V.d2V_2_rhs(w.kappa(), v, conj);
            const auto delta_AK = d2V.dV(w, v, conj);
            const auto sign = conj ? iota_ts : 1.0;
            REQUIRE(delta_AK != 0.0);
            REQUIRE(delta_AK == Approx(sign * (d2V_11 + d2V_Norm + d2V_2)));
          }

          // <w||d2V_--||v> = iota_t*iota_s*(-1)^{jw-jv} <v||d2V_++||w>;
          // holds for each part. (The right-hand sides of -- carry the
          // factor iota_t*iota_s.)
          const auto sign_wv = Angular::neg1pow_2(w.twoj() - v.twoj());
          REQUIRE(w * d2V.d2V_11_rhs(w.kappa(), v, true) ==
                  Approx(sign_wv * (v * d2V.d2V_11_rhs(v.kappa(), w, false))));
          REQUIRE(w * d2V.d2V_Norm_rhs(w.kappa(), v) ==
                  Approx(sign_wv * (v * d2V.d2V_Norm_rhs(v.kappa(), w))));
          REQUIRE(w * d2V.d2V_2_rhs(w.kappa(), v, true) ==
                  Approx(sign_wv * (v * d2V.d2V_2_rhs(v.kappa(), w, false))));
          REQUIRE(d2V.dV(w, v, true) ==
                  Approx(iota_ts * sign_wv * d2V.dV(v, w, false)));
        }
      }
      REQUIRE(count != 0);
    }
  }

  //----------------------------------------------------------------------------
  SECTION("swapping t and s changes only the order of the coupling") {

    // A^K for (s, t) is (-1)^{kt+ks-K} times A^K for (t, s)
    const std::vector<OperatorPair> pairs{
      {&d, &d, 0}, {&d, &d, 1}, {&d, &q, 2}, {&d, &hpnc, 1}};

    // Both frequencies non-zero, and different
    const auto omega_t = 0.05;
    const auto omega_s = 0.02;

    for (const auto &[t, s, K] : pairs) {

      ExternalField::TDHF T(t, wf.vHF());
      ExternalField::TDHF S(s, wf.vHF());
      T.solve_core(omega_t, 30, false);
      S.solve_core(omega_s, 30, false);

      ExternalField::TDHF_DCP d2V_ts(&T, &S, K, wf.vHF());
      ExternalField::TDHF_DCP d2V_st(&S, &T, K, wf.vHF());
      d2V_ts.solve_core(omega_t + omega_s, 30, false);
      d2V_st.solve_core(omega_t + omega_s, 30, false);
      REQUIRE(d2V_ts.last_eps() < 1.0e-10);
      REQUIRE(d2V_st.last_eps() < 1.0e-10);

      const auto sign_swap = Angular::neg1pow(t->rank() + s->rank() - K);

      int count = 0;
      for (const auto &w : wf.valence()) {
        for (const auto &v : wf.valence()) {
          if (w.parity() * v.parity() != d2V_ts.parity() ||
              Angular::triangle(w.twoj(), v.twoj(), 2 * K) == 0)
            continue;
          ++count;
          REQUIRE(d2V_ts.dV(w, v) != 0.0);
          REQUIRE(d2V_st.dV(w, v) == Approx(sign_swap * d2V_ts.dV(w, v)));
        }
      }
      REQUIRE(count != 0);
    }
  }

  //----------------------------------------------------------------------------
  SECTION("static limit: solving for dY reproduces the static relation") {

    // In the static case dY is not solved for:
    // dY = (-1)^{K+kt+ks} * iota_t * iota_s * dX.
    // At a very small frequency both dX and dY are solved for, with all of
    // the -- terms; the result must be the same, up to terms linear in the
    // frequency.
    // (Not t = s with odd K: that amplitude is zero in the static limit.)
    const std::vector<OperatorPair> pairs{
      {&d, &d, 0}, {&d, &q, 2}, {&d, &hpnc, 1}};

    // The linear terms can be large: energy denominators between the
    // fine-structure partners of the core (2p-, 2p+) are small. So the
    // frequency is tiny, and differences are measured relative to the
    // largest correction (or amplitude), not channel by channel.
    const auto omega_small = 1.0e-8;
    const auto tolerance = 1.0e-4;

    for (const auto &[t, s, K] : pairs) {

      // Static: only dX solved for
      ExternalField::TDHF T_static(t, wf.vHF());
      ExternalField::TDHF S_static(s, wf.vHF());
      T_static.solve_core(0.0, 30, false);
      S_static.solve_core(0.0, 30, false);
      ExternalField::TDHF_DCP d2V_static(&T_static, &S_static, K, wf.vHF());
      d2V_static.solve_core(0.0, 30, false);
      REQUIRE(d2V_static.last_eps() < 1.0e-10);

      // Nearly static: dX and dY both solved for
      ExternalField::TDHF T(t, wf.vHF());
      ExternalField::TDHF S(s, wf.vHF());
      T.solve_core(omega_small, 30, false);
      S.solve_core(0.0, 30, false);
      ExternalField::TDHF_DCP d2V(&T, &S, K, wf.vHF());
      d2V.solve_core(omega_small, 30, false);
      REQUIRE(d2V.last_eps() < 1.0e-10);

      const auto sign_Y = Angular::neg1pow(K + t->rank() + s->rank()) *
                          (d2V.imagQ() ? -1.0 : 1.0);

      // dY = sign_Y * dX, for every core orbital and channel
      double largest_dX2 = 0.0;
      for (const auto &Fa : wf.core()) {
        for (const auto &dX : d2V.get_dPsis(Fa, ExternalField::dPsiType::X)) {
          largest_dX2 = std::max(largest_dX2, dX.norm2());
        }
      }
      REQUIRE(largest_dX2 != 0.0);
      for (const auto &Fa : wf.core()) {
        const auto &dX = d2V.get_dPsis(Fa, ExternalField::dPsiType::X);
        const auto &dY = d2V.get_dPsis(Fa, ExternalField::dPsiType::Y);
        REQUIRE(dX.size() == dY.size());
        for (std::size_t i = 0; i < dX.size(); ++i) {
          const auto difference2 = (dY[i] - sign_Y * dX[i]).norm2();
          REQUIRE(difference2 < tolerance * tolerance * largest_dX2);
        }
      }

      // and the same amplitudes
      double largest_dV = 0.0;
      for (const auto &w : wf.valence()) {
        for (const auto &v : wf.valence()) {
          largest_dV = std::max(largest_dV, std::abs(d2V_static.dV(w, v)));
        }
      }
      REQUIRE(largest_dV != 0.0);
      for (const auto &w : wf.valence()) {
        for (const auto &v : wf.valence()) {
          if (w.parity() * v.parity() != d2V.parity() ||
              Angular::triangle(w.twoj(), v.twoj(), 2 * K) == 0)
            continue;
          const auto difference = d2V.dV(w, v) - d2V_static.dV(w, v);
          REQUIRE(std::abs(difference) < tolerance * largest_dV);
        }
      }
    }
  }
}

//==============================================================================
TEST_CASE("TDHF_DCP: lowest order agrees with diagrams",
          "[ExternalField][TDHF_DCP][DCP]") {

  // The lowest-order DCP correction (a single iteration of the second-order
  // TDHF equations) is the sum of the Goldstone diagrams, with the same
  // core polarisation of t and s in both.

  Wavefunction wf({1600, 1.0e-6, 100.0, 40.0, "loglinear", -1.0},
                  {"Cs", 133, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Xe]");
  wf.solve_valence("7sp");
  SplineBasis::Parameters params{"40spd20fg", 45, 7, {1.0e-4}, 40.0};
  params.verbose = false;
  wf.formBasis(params);

  const auto [core, excited] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());

  const DiracOperator::E1 d(wf.grid());
  const DiracOperator::PNCnsi hpnc(5.67073, 2.3, wf.grid());

  // t carries the whole frequency of the transition; s is static
  struct Transition {
    std::string label;
    const DiracOperator::TensorOperator *t;
    const DiracOperator::TensorOperator *s;
    int K;
    std::string initial;
    std::string final;
  };
  const std::vector<Transition> transitions{
    {"PNC transition", &d, &hpnc, 1, "6s+", "7s+"},
    {"alpha transition", &d, &d, 0, "6s+", "7s+"},
    {"alpha static", &d, &d, 0, "6s+", "6s+"}};

  const auto tolerance = 1e-3;

  for (const auto &[label, t, s, K, initial, final] : transitions) {

    const auto *v_ptr = wf.getState(initial);
    const auto *w_ptr = wf.getState(final);
    REQUIRE(v_ptr != nullptr);
    REQUIRE(w_ptr != nullptr);
    const auto &v = *v_ptr;
    const auto &w = *w_ptr;

    const auto omega_t = w.en() - v.en();
    const auto omega_s = 0.0;

    std::cout << "\n" << label << "\n" << std::flush;
    fmt::print("operators: t = {} (w = {:.5f}), s = {} (w = {:.5f}), K = {}\n",
               t->name(), omega_t, s->name(), omega_s, K);
    fmt::print("transition: {} -> {}\n", v.shortSymbol(), w.shortSymbol());

    // First order: core polarisation of t and s
    ExternalField::TDHF T(t, wf.vHF());
    ExternalField::TDHF S(s, wf.vHF());
    T.solve_core(omega_t);
    S.solve_core(omega_s);

    // TDHF: lowest order, then all orders
    ExternalField::TDHF_DCP d2V(&T, &S, K, wf.vHF());
    d2V.solve_core(omega_t + omega_s, 1);
    const auto tdhf = d2V.dV(w, v);
    d2V.solve_core(omega_t + omega_s);
    const auto tdhf_all_orders = d2V.dV(w, v);

    // Diagrams
    std::cout << "Calculating DCP diagrams.." << std::flush;
    const auto T_me = Amplitudes::me_table(wf.basis(), t, &T, nullptr, omega_t);
    const auto S_me = Amplitudes::me_table(wf.basis(), s, &S, nullptr, omega_s);
    const auto diagram = Amplitudes::dcp(w, v, omega_t, omega_s, K, t->rank(),
                                         s->rank(), T_me, S_me, core, excited);
    std::cout << " done\n" << std::flush;

    const auto eps = std::abs((tdhf - diagram) / diagram);

    std::cout << "\n";
    fmt::print("{:18s} {:>13s} {:>13s} {:>8s}\n", "Term", "Diagram", "TDHF",
               "eps");
    fmt::print("{:18s} {:13.5e} {:13.5e} {:8.1e}\n", "Total (1)", diagram, tdhf,
               eps);
    fmt::print("{:18s} {:>13s} {:13.5e}\n", "Total (all orders)", "",
               tdhf_all_orders);

    REQUIRE(eps < tolerance);

    // Iterating converges, and changes the result
    REQUIRE(d2V.last_eps() < 1.0e-10);
    REQUIRE(tdhf_all_orders != Approx(tdhf));
  }
}
