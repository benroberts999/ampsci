#include "TDHF_DCP.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
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
  SECTION("parts add up, and -- is the conjugate of ++") {

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
