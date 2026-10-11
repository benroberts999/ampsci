#include "MixedStates.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <string>

//==============================================================================
// Tests the CI mixed states against the direct sum-over-states:
//   <A|dPsi_0> = <A||h||Psi_0> / (E_0 - E_A)
// which holds for each CI solution, A, in the block of the mixed state
TEST_CASE("CI: mixed states", "[CI][MixedStates][unit]") {

  std::cout << "CI mixed states (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, {1.0e-2}, 20.0));

  const auto ci_basis = CI::basis_subset(wf.basis(), "5spd");

  Coulomb::QkTable qk;
  {
    const Coulomb::YkTable yk(ci_basis);
    qk.fill(ci_basis, yk, 8, false);
  }
  const auto h1 = CI::calculate_h1_table(ci_basis, {}, {}, qk, false);

  // Reference state: the ground state, J=0, even parity
  CI::PsiJPi Psi0(0, +1, ci_basis);
  Psi0.solve(CI::construct_Hci(Psi0, h1, qk));
  const auto E0 = Psi0.energy(0);

  // Solve mixed state for operator h, and compare to the sum-over-states
  // Returns largest relative deviation, over all solutions of target block
  const auto test_mixed_state = [&](int twoJ, int parity,
                                    const DiracOperator::TensorOperator *h) {
    const auto h_tab = Amplitudes::me_table(ci_basis, h);

    // Solve CI for the block the mixed state lives in
    CI::PsiJPi target(twoJ, parity, ci_basis);
    const auto Hci = CI::construct_Hci(target, h1, qk);
    target.solve(Hci);

    const auto dPsi =
      CI::solve_mixed_state(Psi0, 0, target, Hci, h_tab, h->rank());
    const auto dc = dPsi.coefs(0);

    fmt::print("\n{}: 2J={}, pi={:+}, K={}\n", h->name(), twoJ, parity,
               h->rank());
    fmt::print("{:>3} {:>12} {:>16} {:>16} {:>9}\n", "A", "E0-EA",
               "<A||h||Psi>", "(E0-EA)<A|dPsi>", "eps");

    double worst = 0.0;
    double max_me = 0.0;
    for (std::size_t iA = 0; iA < target.num_solutions(); ++iA) {

      const auto dE = E0 - target.energy(iA);

      // Direct: sum over CSFs of both states
      const auto me =
        CI::ReducedME(target, iA, Psi0, 0, h_tab, h->rank(), h->parity());

      // From mixed state: <A|dPsi> * (E0 - EA)
      const auto cA = target.coefs(iA);
      double overlap = 0.0;
      for (std::size_t i = 0; i < cA.size(); ++i) {
        overlap += cA[i] * dc[i];
      }
      const auto me_ms = dE * overlap;

      // Psi0 itself is projected out of the mixed state (same J and parity)
      if (std::abs(dE) < 1.0e-10) {
        REQUIRE(std::abs(overlap) < 1.0e-10);
        continue;
      }

      max_me = std::max(max_me, std::abs(me));
      worst = std::max(worst, std::abs(me_ms - me));

      // Print first few, and the last, solution
      const auto last = target.num_solutions() - 1;
      if (iA == last && last > 5) {
        fmt::print("...\n");
      }
      if (iA < 5 || iA == last) {
        fmt::print("{:>3} {:12.6f} {:16.8e} {:16.8e} {:9.1e}\n", iA, dE, me,
                   me_ms, std::abs(me_ms - me));
      }
    }

    const auto eps = worst / max_me;
    fmt::print("worst: {:.1e}\n", eps);
    return eps;
  };

  //----------------------------------------------------------------------------
  // E1: rank 1, odd parity => 2J=2, odd
  const DiracOperator::E1 d{wf.grid()};
  REQUIRE(test_mixed_state(2, -1, &d) < 1.0e-10);

  // PNC: rank 0, odd parity => 2J=0, odd
  const DiracOperator::PNCnsi hpnc{1.0, Nuclear::default_t, wf.grid()};
  REQUIRE(test_mixed_state(0, -1, &hpnc) < 1.0e-10);

  // Scalar r^2: rank 0, even parity => same block as Psi0 (singular case)
  const DiracOperator::RadialF r2{wf.grid(), 2.0};
  REQUIRE(test_mixed_state(0, +1, &r2) < 1.0e-10);

  //----------------------------------------------------------------------------
  // As above, but at non-zero frequency, and with states projected out:
  //   <A|dPsi> = <A||h||Psi_0> / (E_0 + omega - E_A),
  // and zero for any A that was projected out
  {
    const auto d_tab = Amplitudes::me_table(ci_basis, &d);

    CI::PsiJPi target(2, -1, ci_basis);
    const auto Hci = CI::construct_Hci(target, h1, qk);
    target.solve(Hci);

    // Remove these solutions from the mixed state
    const std::vector<std::size_t> project{0, 2};

    const auto omega = 0.05;
    const auto dPsi = CI::project_out(
      CI::solve_mixed_state(Psi0, 0, target, Hci, d_tab, d.rank(), omega),
      target, project);
    const auto dc = dPsi.coefs(0);

    fmt::print("\n{}: 2J={}, pi={:+}, omega={}, projecting out 0 and 2\n",
               d.name(), 2, -1, omega);
    fmt::print("{:>3} {:>12} {:>16} {:>16} {:>9}\n", "A", "E0+w-EA",
               "<A||h||Psi>", "(E0+w-EA)<A|dPsi>", "eps");

    double worst = 0.0;
    double max_me = 0.0;
    for (std::size_t iA = 0; iA < target.num_solutions(); ++iA) {

      const auto dE = E0 + omega - target.energy(iA);
      const auto me =
        CI::ReducedME(target, iA, Psi0, 0, d_tab, d.rank(), d.parity());

      const auto cA = target.coefs(iA);
      double overlap = 0.0;
      for (std::size_t i = 0; i < cA.size(); ++i) {
        overlap += cA[i] * dc[i];
      }

      // Projected-out solutions have no part in the mixed state
      const auto projectedQ =
        std::find(project.begin(), project.end(), iA) != project.end();
      if (projectedQ) {
        REQUIRE(std::abs(overlap * dE) < 1.0e-10 * std::abs(me));
        continue;
      }

      const auto me_ms = dE * overlap;
      max_me = std::max(max_me, std::abs(me));
      worst = std::max(worst, std::abs(me_ms - me));

      const auto last = target.num_solutions() - 1;
      if (iA == last && last > 5) {
        fmt::print("...\n");
      }
      if (iA < 5 || iA == last) {
        fmt::print("{:>3} {:12.6f} {:16.8e} {:16.8e} {:9.1e}\n", iA, dE, me,
                   me_ms, std::abs(me_ms - me));
      }
    }
    fmt::print("worst: {:.1e}\n", worst / max_me);
    REQUIRE(worst / max_me < 1.0e-10);
  }

  //----------------------------------------------------------------------------
  // The overload that constructs the CI matrix itself must agree
  {
    const auto d_tab = Amplitudes::me_table(ci_basis, &d);

    CI::PsiJPi target(2, -1, ci_basis);
    const auto dPsi1 = CI::solve_mixed_state(
      Psi0, 0, target, CI::construct_Hci(target, h1, qk), d_tab, d.rank());

    const auto dPsi2 =
      CI::solve_mixed_state(Psi0, 0, 2, -1, ci_basis, d_tab, d.rank(), h1, qk);

    REQUIRE(dPsi2.twoJ() == 2);
    REQUIRE(dPsi2.parity() == -1);
    REQUIRE(dPsi2.num_solutions() == 1);
    REQUIRE(dPsi2.energy(0) == Approx(E0));
    REQUIRE(dPsi1.coefs(0).size() == dPsi2.coefs(0).size());
    for (std::size_t i = 0; i < dPsi1.coefs(0).size(); ++i) {
      REQUIRE(dPsi1.coefs(0)[i] == Approx(dPsi2.coefs(0)[i]));
    }
  }
}
