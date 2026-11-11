#include "Amplitudes/DCP.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <optional>
#include <vector>

//==============================================================================
TEST_CASE("DCP: diagrams basic unit tests", "[Amplitudes][DCP][unit]") {

  // nb: small grid, very small basis. Don't care about numerical results,
  // just that everything is consistent.
  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);
  wf.solve_valence("3sp", false);
  SplineBasis::Parameters params{"8sp6d", 15, 7, {1e-3}, 30.0};
  params.verbose = false;
  wf.formBasis(params);

  const auto [core, excited] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
  REQUIRE(!core.empty());
  REQUIRE(!excited.empty());

  // Operators, and their matrix elements (no core polarisation)
  const DiracOperator::E1 d(wf.grid());
  const DiracOperator::Ek q(wf.grid(), 2);
  const DiracOperator::PNCnsi hpnc(5.67073, 2.3, wf.grid());
  const auto d_me = Amplitudes::me_table(wf.basis(), &d);
  const auto q_me = Amplitudes::me_table(wf.basis(), &q);
  const auto hpnc_me = Amplitudes::me_table(wf.basis(), &hpnc);

  // Pair of operators (t, s), with their tables, coupled to rank K.
  // Between them: both signs of (-1)^{K+kt+ks}, and an imaginary operator
  struct OperatorPair {
    const DiracOperator::TensorOperator *t;
    const Coulomb::meTable<double> *T;
    const DiracOperator::TensorOperator *s;
    const Coulomb::meTable<double> *S;
    int K;
  };
  const std::vector<OperatorPair> pairs{{&d, &d_me, &d, &d_me, 0},
                                        {&d, &d_me, &q, &q_me, 2},
                                        {&d, &d_me, &hpnc, &hpnc_me, 1}};

  // Both frequencies non-zero, and different
  const auto omega_t = 0.03;
  const auto omega_s = 0.01;

  for (const auto &[t, T, s, S, K] : pairs) {
    const auto kt = t->rank();
    const auto ks = s->rank();
    const auto parity = t->parity() * s->parity();
    const auto iota_ts = t->imaginaryQ() != s->imaginaryQ() ? -1.0 : 1.0;

    int count = 0;
    for (const auto &w : wf.valence()) {
      for (const auto &v : wf.valence()) {
        if (w.parity() * v.parity() != parity ||
            Angular::triangle(w.twoj(), v.twoj(), 2 * K) == 0)
          continue;
        ++count;

        const auto C1 = Amplitudes::dcp_C1(w, v, omega_t, omega_s, K, kt, ks,
                                           *T, *S, core, excited);
        const auto C2 = Amplitudes::dcp_C2(w, v, omega_t, omega_s, K, kt, ks,
                                           *T, *S, core, excited);
        const auto R = Amplitudes::dcp_R12(w, v, omega_t, omega_s, K, kt, ks,
                                           *T, *S, core, excited);
        const auto L = Amplitudes::dcp_L12(w, v, omega_t, omega_s, K, kt, ks,
                                           *T, *S, core, excited);

        // Total is the sum of the terms
        const auto total = Amplitudes::dcp(w, v, omega_t, omega_s, K, kt, ks,
                                           *T, *S, core, excited);
        REQUIRE(total != 0.0);
        REQUIRE(total == Approx(C1 + C2 + R + L));

        // Hermitian conjugate: swap w and v, reverse the frequencies.
        // The right and left diagrams swap; the centre ones map to
        // themselves. Phase from the conjugate of the pair coupled to
        // rank K, and of an imaginary operator
        const auto phase = iota_ts * Angular::neg1pow(K + kt + ks) *
                           Angular::neg1pow_2(w.twoj() - v.twoj());
        REQUIRE(Amplitudes::dcp_L12(v, w, -omega_t, -omega_s, K, kt, ks, *T, *S,
                                    core, excited) == Approx(phase * R));
        REQUIRE(Amplitudes::dcp_R12(v, w, -omega_t, -omega_s, K, kt, ks, *T, *S,
                                    core, excited) == Approx(phase * L));
        REQUIRE(Amplitudes::dcp_C1(v, w, -omega_t, -omega_s, K, kt, ks, *T, *S,
                                   core, excited) == Approx(phase * C1));
        REQUIRE(Amplitudes::dcp_C2(v, w, -omega_t, -omega_s, K, kt, ks, *T, *S,
                                   core, excited) == Approx(phase * C2));

        // Swapping t and s (with their frequencies) changes the order of the
        // coupling only: each term changes by (-1)^{kt+ks-K}
        const auto sign_swap = Angular::neg1pow(kt + ks - K);
        REQUIRE(Amplitudes::dcp_C1(w, v, omega_s, omega_t, K, ks, kt, *S, *T,
                                   core, excited) == Approx(sign_swap * C1));
        REQUIRE(Amplitudes::dcp_C2(w, v, omega_s, omega_t, K, ks, kt, *S, *T,
                                   core, excited) == Approx(sign_swap * C2));
        REQUIRE(Amplitudes::dcp_R12(w, v, omega_s, omega_t, K, ks, kt, *S, *T,
                                    core, excited) == Approx(sign_swap * R));
        REQUIRE(Amplitudes::dcp_L12(w, v, omega_s, omega_t, K, ks, kt, *S, *T,
                                    core, excited) == Approx(sign_swap * L));
      }
    }
    REQUIRE(count != 0);
  }
}
