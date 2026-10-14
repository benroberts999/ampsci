#include "SecondOrder.hpp"
#include "Amplitudes/MatrixElements.hpp"
#include "Angular/include.hpp"
#include "CI_Integrals.hpp"
#include "Coulomb/include.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <sstream>
#include <vector>

//==============================================================================
// Wigner-Eckart factor for the q=0 component:
// <a,m|h^k_0|b,m> = we_factor * <a||h^k||b>  (cf TensorOperator::rme3js)
static double we_factor(int k, int twoja, int twojb, int two_m) {
  return Angular::neg1pow_2(twoja - two_m) *
         Angular::threej_2(twoja, 2 * k, twojb, -two_m, 0, two_m);
}

//==============================================================================
// Tests the angular coefficients of the second-order amplitude A^K.
// A^K is the reduced ME of the composite operator [t x s]^K, so for each
// intermediate state the sum over K of the z-component factors must reproduce
// the product of the two Wigner-Eckart factors of <b,m|t_0|n,m><n,m|s_0|a,m>
TEST_CASE("CI: second-order angular coefficients", "[CI][SecondOrder][unit]") {

  std::cout << "CI second-order amplitude: angular coefficients\n";

  double worst = 0.0;
  int num_tested = 0;

  for (int kt = 0; kt <= 2; ++kt) {
    for (int ks = 0; ks <= 2; ++ks) {

      for (int twoJa = 0; twoJa <= 6; twoJa += 2) {
        for (int twoJb = 0; twoJb <= 6; twoJb += 2) {
          for (int twoJn = 0; twoJn <= 8; twoJn += 2) {

            if (Angular::triangle(twoJb, 2 * kt, twoJn) == 0)
              continue;
            if (Angular::triangle(twoJn, 2 * ks, twoJa) == 0)
              continue;

            const auto two_m = std::min(twoJa, twoJb);

            // The two factors, for each of the two terms of A^K
            const auto zz_1 = we_factor(kt, twoJb, twoJn, two_m) *
                              we_factor(ks, twoJn, twoJa, two_m);
            const auto zz_2 = we_factor(ks, twoJb, twoJn, two_m) *
                              we_factor(kt, twoJn, twoJa, two_m);

            // The same, summed over the ranks K
            double sum_1 = 0.0;
            double sum_2 = 0.0;
            for (int K = std::abs(kt - ks); K <= kt + ks; ++K) {
              const auto [c1, c2] =
                CI::A_K_coefs(K, kt, ks, twoJb, twoJn, twoJa);
              const auto z = CI::z_component(K, kt, ks, twoJb, twoJa, two_m);
              sum_1 += z * c1;
              sum_2 += z * c2;
            }

            const auto scale = std::max({std::abs(zz_1), std::abs(zz_2), 1.0});
            worst = std::max({worst, std::abs(sum_1 - zz_1) / scale,
                              std::abs(sum_2 - zz_2) / scale});
            ++num_tested;
          }
        }
      }
    }
  }

  fmt::print("{} cases, worst: {:.1e}\n", num_tested, worst);
  REQUIRE(num_tested > 100);
  REQUIRE(worst < 1.0e-12);
}

//==============================================================================
// Reduced ME of the total spin between LS-coupled states of the same L and S:
//   <LSJb||S||LSJa> = (-1)^(L+S+Jb+1) sqrt([Ja][Jb]) {S Jb L; Ja S 1}
//                     * sqrt(S(S+1)(2S+1))
// A reference for CI::sigma_rme, which does not assume L or S are good
// quantum numbers (they are not, for relativistic CI states)
static double spin_rme_LS(int L, int twoS, int twoJb, int twoJa) {
  const auto rme_S = std::sqrt(double(twoS) * (twoS + 2) * (twoS + 1)) / 2;
  return Angular::neg1pow_2(2 * L + twoS + twoJb + 2) *
         std::sqrt(double(twoJa + 1) * (twoJb + 1)) *
         Angular::sixj_2(twoS, twoJb, 2 * L, twoJa, twoS, 2) * rme_S;
}

//==============================================================================
// The LS reference above must reproduce Angular::S_kk for a single electron
// (L = l, S = 1/2, J = j)
TEST_CASE("CI: spin reduced ME", "[CI][SecondOrder][unit]") {

  std::cout << "CI: <LSJ'||S||LSJ> against S_kk\n";

  double worst = 0.0;
  int num_tested = 0;

  for (int l = 0; l <= 4; ++l) {
    // kappa = -l-1 (j=l+1/2) and kappa = l (j=l-1/2)
    const std::vector<int> kappas =
      l == 0 ? std::vector{-1} : std::vector{-l - 1, l};
    for (const auto ka : kappas) {
      for (const auto kb : kappas) {
        const auto twoja = Angular::twoj_k(ka);
        const auto twojb = Angular::twoj_k(kb);
        // nb: S_kk(ka,kb) is <ka||S||kb>, i.e., 'a' is the bra
        const auto expected = Angular::S_kk(ka, kb);
        const auto found = spin_rme_LS(l, 1, twoja, twojb);
        worst = std::max(worst, std::abs(found - expected));
        ++num_tested;
        fmt::print("l={}: <{}/2||S||{}/2> = {:12.8f} [{:12.8f}]\n", l, twoja,
                   twojb, found, expected);
      }
    }
  }

  fmt::print("{} cases, worst: {:.1e}\n", num_tested, worst);
  REQUIRE(num_tested == 17);
  REQUIRE(worst < 1.0e-12);
}

//==============================================================================
// Tests the second-order amplitude for a two-electron CI:
//  - the two routes (mixed states of s, and of t) must agree
//  - the static scalar polarisability of the He ground state, which fixes the
//    overall sign and normalisation
TEST_CASE("CI: second-order amplitudes", "[CI][SecondOrder][unit]") {

  std::cout << "CI second-order amplitudes (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, {1.0e-2}, 20.0));

  CI::Integrals ints;
  ints.ci_basis = CI::basis_subset(wf.basis(), "5spd");
  {
    const Coulomb::YkTable yk(ints.ci_basis);
    ints.qk.fill(ints.ci_basis, yk, 8, false);
  }
  ints.h1 = CI::calculate_h1_table(ints.ci_basis, {}, {}, ints.qk, false);
  REQUIRE(ints.availableQ());

  // J=0 even (ground state and first excited), and J=1 even
  CI::PsiJPi Psi_0e(0, +1, ints.ci_basis);
  Psi_0e.solve(CI::construct_Hci(Psi_0e, ints));
  CI::PsiJPi Psi_1e(2, +1, ints.ci_basis);
  Psi_1e.solve(CI::construct_Hci(Psi_1e, ints));

  const DiracOperator::E1 d{wf.grid()};
  const DiracOperator::PNCnsi hpnc{1.0, Nuclear::default_t, wf.grid()};
  const auto d_me = Amplitudes::me_table(ints.ci_basis, &d);
  const auto pnc_me = Amplitudes::me_table(ints.ci_basis, &hpnc);

  //----------------------------------------------------------------------------
  // The spin matrix element that defines beta. The lowest J=1 even solution is
  // 1s2s 3S_1: two s-electrons coupled to J=1 can only be a triplet S, so it
  // is (almost) pure L=0, S=1, and <b||sigma||a> = 2<LSJ||S||LSJ>
  {
    const auto sigma = CI::sigma_rme(Psi_1e, 0, Psi_1e, 0, ints.ci_basis);
    const auto expected = 2.0 * spin_rme_LS(0, 2, 2, 2);
    fmt::print("\n<b||sigma||a> for 1s2s 3S_1: {:.6f} [{:.6f}]\n", sigma,
               expected);
    REQUIRE(std::abs(sigma - expected) / expected < 1.0e-3);
  }

  //----------------------------------------------------------------------------
  // Static polarisability of the ground state (all K), and the dynamic case,
  // for which s carries -omega. The only intermediate J/Pi is 1-, so
  // the direct sum over its complete CI spectrum is the reference
  {
    CI::PsiJPi Psi_1o(2, -1, ints.ci_basis);
    Psi_1o.solve(CI::construct_Hci(Psi_1o, ints));
    const auto [c1, c2] = CI::A_K_coefs(0, 1, 1, 0, 2, 0);

    for (const auto omega : {0.0, 0.05}) {
      const auto omega_s = -omega;
      fmt::print("\nE1-E1, J=0 -> J=0, omega = {:.2f}\n", omega);
      const auto [A_ket, A_bra] = CI::A_K(0, Psi_0e, 0, Psi_0e, 0, &d, d_me, &d,
                                          d_me, omega, omega_s, ints);
      const auto eps = std::abs(A_ket - A_bra) / std::abs(A_ket);
      fmt::print("A^0: {:.8e} (ket), {:.8e} (bra) [{:.1e}]\n", A_ket, A_bra,
                 eps);
      REQUIRE(eps < 1.0e-10);

      // Direct sum over states: c1/(E_a + omega_s - E_p) + c2/(E_a + omega -
      // E_p)
      double A_sos = 0.0;
      for (std::size_t ip = 0; ip < Psi_1o.num_solutions(); ++ip) {
        const auto dE = Psi_0e.energy(0) - Psi_1o.energy(ip);
        const auto d_bp = CI::ReducedME(Psi_0e, 0, Psi_1o, ip, d_me, 1, -1);
        const auto d_pa = CI::ReducedME(Psi_1o, ip, Psi_0e, 0, d_me, 1, -1);
        A_sos += d_bp * d_pa * (c1 / (dE + omega_s) + c2 / (dE + omega));
      }
      const auto eps_sos = std::abs(A_ket - A_sos) / std::abs(A_sos);
      fmt::print("     {:.8e} (sum over states)          [{:.1e}]\n", A_sos,
                 eps_sos);
      REQUIRE(eps_sos < 1.0e-10);

      if (omega == 0.0) {
        // alpha_0 = +A^0/sqrt(3[J]); exact for He is 1.383 au, but this basis
        // is small: only the sign and rough magnitude are meaningful
        const auto alpha = A_ket / std::sqrt(3.0);
        fmt::print("alpha_0 = {:.4f} au (expect ~1.4)\n", alpha);
        REQUIRE(alpha > 1.0);
        REQUIRE(alpha < 2.0);
      }
    }
  }

  //----------------------------------------------------------------------------
  // PNC-like amplitude: both operators odd, so a and b have the same parity.
  // Tests the case k_t != k_s (the two terms then have different J_n)
  {
    fmt::print("\nE1-pnc, J=0 -> J=1\n");
    const auto omega = Psi_1e.energy(0) - Psi_0e.energy(0);
    const auto [A_ket, A_bra] = CI::A_K(1, Psi_1e, 0, Psi_0e, 0, &d, d_me,
                                        &hpnc, pnc_me, omega, 0.0, ints);
    const auto eps = std::abs(A_ket - A_bra) / std::abs(A_ket);
    fmt::print("A^1: {:.8e} (ket), {:.8e} (bra) [{:.1e}]\n", A_ket, A_bra, eps);
    REQUIRE(eps < 1.0e-10);
    REQUIRE(std::abs(A_ket) > 1.0e-12);
  }

  //----------------------------------------------------------------------------
  // Transition amplitude between the two lowest J=0 even states: K = 0 only
  {
    fmt::print("\nE1-E1, transition between two J=0 even states\n");
    const auto omega = Psi_0e.energy(1) - Psi_0e.energy(0);
    const auto [A_ket, A_bra] =
      CI::A_K(0, Psi_0e, 1, Psi_0e, 0, &d, d_me, &d, d_me, omega, 0.0, ints);
    const auto eps = std::abs(A_ket - A_bra) / std::abs(A_ket);
    fmt::print("A^0: {:.8e} (ket), {:.8e} (bra) [{:.1e}]\n", A_ket, A_bra, eps);
    REQUIRE(eps < 1.0e-10);
  }

  //----------------------------------------------------------------------------
  // Projecting out intermediate levels: the amplitude must drop by exactly
  // their sum-over-states contributions
  {
    fmt::print("\nProjecting out intermediate levels\n");
    const auto A_all =
      CI::A_K(0, Psi_0e, 0, Psi_0e, 0, &d, d_me, &d, d_me, 0.0, 0.0, ints)
        .first;

    const std::vector<CI::Level> project{{2, -1, 0}, {2, -1, 2}};
    const auto A_cut = CI::A_K(0, Psi_0e, 0, Psi_0e, 0, &d, d_me, &d, d_me, 0.0,
                               0.0, ints, project)
                         .first;

    // The removed levels, by direct sum-over-states
    CI::PsiJPi Psi_1o(2, -1, ints.ci_basis);
    Psi_1o.solve(CI::construct_Hci(Psi_1o, ints));
    double removed = 0.0;
    for (const auto &level : project) {
      const auto ip = level.index;
      const auto dE = Psi_0e.energy(0) - Psi_1o.energy(ip);
      const auto d_bp = CI::ReducedME(Psi_0e, 0, Psi_1o, ip, d_me, 1, -1);
      const auto d_pa = CI::ReducedME(Psi_1o, ip, Psi_0e, 0, d_me, 1, -1);
      const auto [c1, c2] = CI::A_K_coefs(0, 1, 1, 0, 2, 0);
      removed += (c1 + c2) * d_bp * d_pa / dE;
    }

    const auto eps = std::abs(A_cut + removed - A_all) / std::abs(A_all);
    fmt::print("A^0: {:.8e} = {:.8e} + {:.8e} [{:.1e}]\n", A_all, A_cut,
               removed, eps);
    REQUIRE(eps < 1.0e-10);
  }
}

//==============================================================================
// With an energy-dependent CI Hamiltonian (the dSigma/dE correction), each J/Pi
// has its own reference energy E0. The mixed states of an intermediate J/Pi
// that was solved in the CI run must be formed at that E0, so that the
// amplitude is the sum over the levels of that run. Without the solved J/Pi,
// they are formed at the ground-state E0: a different Hamiltonian
TEST_CASE("CI: second-order amplitude energy-dependent Hamiltonian",
          "[CI][SecondOrder][unit]") {

  std::cout << "\nCI second-order: energy-dependent Hamiltonian\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, {1.0e-2}, 20.0));

  CI::Integrals ints;
  ints.ci_basis = CI::basis_subset(wf.basis(), "5spd");
  {
    const Coulomb::YkTable yk(ints.ci_basis);
    ints.qk.fill(ints.ci_basis, yk, 8, false);
  }
  ints.h1 = CI::calculate_h1_table(ints.ci_basis, {}, {}, ints.qk, false);

  // A synthetic Sigma_1 correction: constant Sigma_1 and dSigma/dE for each
  // orbital, with the reference energies e_sigma at zero, so that the
  // correction is active (dE < 0) and depends on E0 through every orbital
  for (const auto &a : ints.ci_basis) {
    ints.s1_corr.S1.add(a, a, -0.02);
    ints.s1_corr.dS1.add(a, a, -0.1);
    ints.s1_corr.e_sigma[a.kappa()] = 0.0;
    ints.s1_corr.en[a.nk_index()] = a.en();
  }

  // The CI run: each J/Pi solved at its own E0. The stored fallback E0 is the
  // ground-state energy
  std::ostringstream quiet;
  CI::PsiJPi Psi_0e(0, +1, ints.ci_basis);
  Psi_0e.solve(CI::iterate_E0(&Psi_0e, ints.s1_corr, ints.h1, ints.qk, nullptr,
                              nullptr, {}, nullptr, 0.0, quiet));
  ints.s1_corr.E0 = Psi_0e.energy(0);

  CI::PsiJPi Psi_1o(2, -1, ints.ci_basis);
  Psi_1o.solve(CI::iterate_E0(&Psi_1o, ints.s1_corr, ints.h1, ints.qk, nullptr,
                              nullptr, {}, nullptr, 0.0, quiet));

  // The same J/Pi, at the ground-state E0
  CI::PsiJPi Psi_1o_gs(2, -1, ints.ci_basis);
  Psi_1o_gs.solve(CI::construct_Hci(Psi_1o_gs, ints));

  fmt::print("E0: ground {:.6f}, 1- {:.6f}, 1- at ground E0 {:.6f}\n",
             Psi_0e.energy(0), Psi_1o.energy(0), Psi_1o_gs.energy(0));
  REQUIRE(std::abs(Psi_1o.energy(0) - Psi_1o_gs.energy(0)) > 1.0e-4);

  const DiracOperator::E1 d{wf.grid()};
  const auto d_me = Amplitudes::me_table(ints.ci_basis, &d);

  // Static scalar polarisability of the ground state: the direct sum over the
  // complete spectrum of the given 1- levels
  const auto sos = [&](const CI::PsiJPi &levels) {
    const auto [c1, c2] = CI::A_K_coefs(0, 1, 1, 0, 2, 0);
    double A = 0.0;
    for (std::size_t ip = 0; ip < levels.num_solutions(); ++ip) {
      const auto dE = Psi_0e.energy(0) - levels.energy(ip);
      const auto d_bp = CI::ReducedME(Psi_0e, 0, levels, ip, d_me, 1, -1);
      const auto d_pa = CI::ReducedME(levels, ip, Psi_0e, 0, d_me, 1, -1);
      A += (c1 + c2) * d_bp * d_pa / dE;
    }
    return A;
  };
  const auto A_sos_own = sos(Psi_1o);
  const auto A_sos_gs = sos(Psi_1o_gs);

  // Mixed states: with the solved 1- given, and without
  const auto A_ms_own = CI::A_K(0, Psi_0e, 0, Psi_0e, 0, &d, d_me, &d, d_me,
                                0.0, 0.0, ints, {}, {Psi_1o}, quiet)
                          .first;
  const auto A_ms_gs = CI::A_K(0, Psi_0e, 0, Psi_0e, 0, &d, d_me, &d, d_me, 0.0,
                               0.0, ints, {}, {}, quiet)
                         .first;

  fmt::print("{:<22} {:>16} {:>16}\n", "A^0", "sum-over-states",
             "mixed states");
  fmt::print("{:<22} {:16.8e} {:16.8e}\n", "1- at own E0", A_sos_own, A_ms_own);
  fmt::print("{:<22} {:16.8e} {:16.8e}\n", "1- at ground E0", A_sos_gs,
             A_ms_gs);

  // The two Hamiltonians must differ, so that the comparison is meaningful
  REQUIRE(std::abs(A_sos_own - A_sos_gs) / std::abs(A_sos_own) > 1.0e-4);
  // The E0 of the solved J/Pi is converged to 1e-8, so the match is not exact
  REQUIRE(std::abs(A_ms_own - A_sos_own) / std::abs(A_sos_own) < 1.0e-6);
  REQUIRE(std::abs(A_ms_gs - A_sos_gs) / std::abs(A_sos_gs) < 1.0e-10);
}

//==============================================================================
// Tests the two contributions from intermediate states with a core hole, for
// which the reference values are the standard sum-over-states expressions of
// the core and core-valence parts of the polarisability (as in the
// polarisability module):
//   alpha_core = -(2/3) sum_{c,m} |<c||d||m>|^2/(e_c - e_m)
//   alpha_cv   = -(2/3) (n_v/[j_v]) sum_c |<v||d||c>|^2/(e_v - e_c)
TEST_CASE("CI: second-order core and core-valence", "[CI][SecondOrder][unit]") {

  std::cout << "\nCI second-order: core and core-valence contributions\n";

  // Mg: the core has a p shell, so the core-valence term is non-zero
  Wavefunction wf({800, 1.0e-5, 40.0, 0.33 * 40.0, "loglinear"},
                  {"Mg", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-8);
  wf.formBasis(SplineBasis::Parameters("15spd", 30, 7, {1.0e-3, 2.9e-2}, 40.0));

  const auto excited =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel()).second;
  REQUIRE(!wf.core().empty());
  REQUIRE(!excited.empty());

  const DiracOperator::E1 d{wf.grid()};

  // A one-orbital CI basis: the only J=0 even CSF is then 3s^2, so the CI
  // state is exactly the closed 3s shell, and the core-valence term must be
  // the Pauli blocking of the c -> 3s excitations, at full weight
  const auto ci_basis =
    CI::basis_subset(wf.basis(), "3s", wf.coreConfiguration());
  REQUIRE(ci_basis.size() == 1);
  const auto &Fv = ci_basis.front();

  CI::Integrals ints;
  ints.ci_basis = ci_basis;
  {
    const Coulomb::YkTable yk(ints.ci_basis);
    ints.qk.fill(ints.ci_basis, yk, 8, false);
  }
  ints.h1 = CI::calculate_h1_table(ints.ci_basis, {}, {}, ints.qk, false);
  REQUIRE(ints.availableQ());

  CI::PsiJPi Psi(0, +1, ints.ci_basis);
  Psi.solve(CI::construct_Hci(Psi, ints));
  REQUIRE(Psi.CSFs().size() == 1);

  //----------------------------------------------------------------------------
  // Core and core-valence: A^0/sqrt(3[J]) must be the core (core-valence)
  // polarisability. Also at omega != 0, where s carries -omega, and each
  // reference sum has the two denominators of the dynamic polarisability
  fmt::print("\n{:<16} {:>14} {:>14} {:>9}\n", "", "expected", "found", "eps");

  for (const auto omega : {0.0, 0.05}) {
    const auto omega_s = -omega;

    // Core
    {
      double expected = 0.0;
      for (const auto &c : wf.core()) {
        for (const auto &m : excited) {
          if (d.isZero(c, m))
            continue;
          const auto d_cm = d.reducedME(c, m);
          const auto de = c.en() - m.en();
          expected += -1.0 / 3.0 * d_cm * d_cm *
                      (1.0 / (de + omega) + 1.0 / (de + omega_s));
        }
      }

      const auto A_core =
        CI::A_K_core(0, Psi.twoJ(), &d, &d, omega, omega_s, wf.core(), excited);
      const auto found = A_core / std::sqrt(3.0 * (Psi.twoJ() + 1));
      const auto eps = std::abs(found - expected) / std::abs(expected);

      fmt::print("{:<16} {:14.8e} {:14.8e} {:9.1e}\n",
                 fmt::format("alpha_core({:.2f})", omega), expected, found,
                 eps);
      REQUIRE(eps < 1.0e-12);
    }

    // Core-valence: the blocking of c -> 3s, with the 3s shell full
    {
      double expected = 0.0;
      for (const auto &c : wf.core()) {
        if (d.isZero(Fv, c))
          continue;
        const auto d_vc = d.reducedME(Fv, c);
        const auto de = Fv.en() - c.en();
        // occupation 2 in a shell of [j] = 2: the blocking weight is one
        expected += -1.0 / 3.0 * d_vc * d_vc *
                    (1.0 / (de - omega) + 1.0 / (de - omega_s));
      }

      const auto A_cv = CI::A_K_cv(0, Psi, 0, Psi, 0, &d, &d, omega, omega_s,
                                   wf.core(), ci_basis);
      const auto found = A_cv / std::sqrt(3.0 * (Psi.twoJ() + 1));
      const auto eps = std::abs(found - expected) / std::abs(expected);

      fmt::print("{:<16} {:14.8e} {:14.8e} {:9.1e}\n",
                 fmt::format("alpha_cv({:.2f})", omega), expected, found, eps);
      REQUIRE(expected < 0.0);
      REQUIRE(eps < 1.0e-12);
    }
  }

  //----------------------------------------------------------------------------
  // The core is a J=0 state: only a scalar, even-parity, amplitude survives
  {
    const auto no_K2 =
      CI::A_K_core(2, Psi.twoJ(), &d, &d, 0.0, 0.0, wf.core(), excited);
    fmt::print("\nA^2_core (must vanish): {:.1e}\n", no_K2);
    REQUIRE(no_K2 == 0.0);

    // sigma is even parity, so the amplitude is odd overall
    const DiracOperator::s spin{};
    const auto no_parity =
      CI::A_K_core(0, Psi.twoJ(), &d, &spin, 0.0, 0.0, wf.core(), excited);
    fmt::print("A^0_core (odd parity, must vanish): {:.1e}\n", no_parity);
    REQUIRE(no_parity == 0.0);
  }
}
