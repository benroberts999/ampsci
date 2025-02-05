#include "DiracODE/DiracODE.hpp"
#include "Angular/Angular.hpp"
#include "DiracODE/AsymptoticSpinor.hpp"
#include "DiracOperator/include.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
//! Unit tests for solving (local) Dirac equation ODE
TEST_CASE("DiracODE: AsymptoticSpinor expansion", "[DiracODE][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "DiracODE: AsymptoticSpinor expansion (large r)\n";

  fmt::print("{:<3s} {:>2s} {:1s} {:1s} {:>5s}  {:16s}  {:17s}   {}\n", "z",
             "k", "l", "n", "r", "f/g (asym)", "f/g (exact)", "eps");

  for (auto z : {0.1, 0.5, 1.0, 10.0, 100.0, 150.0}) {
    for (auto kappa : {-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6}) {
      const auto l = Angular::l_k(kappa);
      double w_eps = -1.0;
      int w_n{0};
      double w_asym{0.0}, w_exact{0.0}, w_r{0.0};

      if (z * PhysConst::alpha > 1.0 && std::abs(kappa) == 1)
        continue;

      const auto r0{10.0 / z};
      const auto rmax{150.0 / z};
      const auto num_grid_points{15ul};
      const auto b{(rmax + r0) / 4};
      const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                     GridType::loglinear, b);
      for (auto n : {1, 2, 3, 4, 5, 6, 7}) {

        if (l >= n)
          continue;

        const auto e = AtomData::diracen(z, n, kappa, PhysConst::alpha);
        const auto F1s = DiracSpinor::exactHlike(n, kappa, grid, z);
        DiracODE::AsymptoticSpinor x{kappa, z, e};
        for (std::size_t i = 0; i < num_grid_points; ++i) {
          const auto r = grid->r(i);
          const auto [f, g] = x.fg(r);
          const auto f0 = F1s.f(i);
          const auto g0 = F1s.g(i);
          if (g == 0.0 || f0 == 0.0 || g0 == 0.0) {
            continue;
          }

          const auto ratio_asym = f / g;
          const auto ratio_exact = f0 / g0;
          const auto eps = std::abs(ratio_asym / ratio_exact - 1.0);

          if (eps > w_eps) {
            w_eps = eps;
            w_n = n;
            w_asym = ratio_asym;
            w_exact = ratio_exact;
            w_r = r;
          }

          auto eps_targ = z > 99.0 ? 1.0e-10 :
                          z > 9.0  ? 1.0e-8 :
                          z > 0.9  ? 1.0e-7 :
                                     1.0e-5;

          if (l >= 5)
            eps_targ *= 100;
          if (l >= 5 && z < 1.0)
            eps_targ *= 100;

          REQUIRE(eps < eps_targ);
        }
      }

      // if (l >= 5)
      fmt::print("{:3g} {:2} {:1} {:1} {:5.1f} {:+.10e} [{:+.10e}]  {:.1e}\n",
                 z, kappa, l, w_n, w_r, w_asym, w_exact, w_eps);
    }
  }
}

//==============================================================================
//! Unit tests for solving (local) Dirac equation ODE
TEST_CASE("DiracODE: Low-r solution", "[DiracODE][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "DiracODE: Low-r solution (pointlike)\n";

  // Set up radial grid:
  const auto r0{1.0e-6};
  const auto rmax{100.0}; // NB: rmax depends on Zeff
  const auto num_grid_points{2000ul};
  const auto b{10.0};
  const auto grid = Grid(r0, rmax, num_grid_points, GridType::loglinear, b);

  std::cout << " Zeff  k n ri r        (f/g)_r          [ Exact           ]  "
               "eps     [targ]\n";
  for (const auto Zeff : {0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 50.0, 100.0}) {
    for (auto kappa : {-1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6}) {
      double w_eps = -1.0;
      int w_n{0};
      std::size_t w_i{0};
      double w_asym{0.0}, w_exact{0.0}, w_r{0.0};
      for (auto n : {1, 2, 3, 4, 5, 6, 7}) {

        const auto l = Angular::l_k(kappa);
        if (l >= n)
          continue;

        const auto e = AtomData::diracen(Zeff, n, kappa, PhysConst::alpha);
        const auto F1s = DiracSpinor::exactHlike(
            n, kappa, std::make_shared<Grid>(grid), Zeff);

        const auto v_nuc =
            Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid.r());

        DiracODE::Internal::DiracDerivative Hd(grid, v_nuc, kappa, e,
                                               PhysConst::alpha);

        std::size_t max_pt = 20;
        std::vector<double> f(max_pt);
        std::vector<double> g(max_pt);
        solve_Dirac_outwards(f, g, Hd, max_pt);

        for (std::size_t i = max_pt / 2; i < max_pt; ++i) {
          const auto r = grid.r(i);
          const auto f0 = F1s.f(i);
          const auto g0 = F1s.g(i);

          const auto ratio_asym = f[i] / g[i];
          const auto ratio_exact = f0 / g0;
          const auto eps = std::abs((ratio_asym - ratio_exact) / ratio_exact);

          if (eps > w_eps) {
            w_eps = eps;
            w_n = n;
            w_asym = ratio_asym;
            w_exact = ratio_exact;
            w_r = r;
            w_i = i;
          }
        }
      }

      const auto magnitude = std::abs(w_exact);
      auto target = //magnitude * 1.0e-4;
          magnitude > 1e2  ? 1.0e-6 :
          magnitude > 1e1  ? 1.0e-5 :
          magnitude > 1e0  ? 1.0e-4 :
          magnitude > 1e-1 ? 1.0e-3 :
          magnitude > 1e-2 ? 1.0e-2 :
                             1.0e-1;
      if (Zeff < 1.0)
        target *= 6.0;

      fmt::print(
          "{:5g} {:2} {} {} {:.1e} {:+.10e} [{:+.10e}]  {:.0e} [{:.0e}]\n",
          Zeff, kappa, w_n, w_i, w_r, w_asym, w_exact, w_eps, target);

      REQUIRE(w_eps < target);
    }
  }
}

//==============================================================================
//! Unit tests for solving (local) Dirac equation ODE
TEST_CASE("DiracODE: Adams-Moulton method", "[DiracODE][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "DiracODE: Adams-Moulton method\n";

  std::cout << "\nTest Hydrogen-like numerical solutions vs. exact Dirac:\n";
  for (const auto Zeff : {0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0}) {
    std::cout << "\nZ = " << Zeff << "\n";

    // Set up radial grid:
    const auto r0{5.0e-7 / Zeff};
    const auto rmax{500.0 / Zeff}; // NB: rmax depends on Zeff
    const auto num_grid_points{2000ul};
    const auto b{10.0};
    const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                   GridType::loglinear, b);

    // States to solve for:
    const std::string states = "10spdfghi";

    // Sperical potential w/ R_nuc = 0.0 is a pointlike potential
    const auto v_nuc = Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid->r());

    // Solve Dirac ODE for each state, store in 'orbitals' vector:
    const auto converge_target = 1.0e-15;
    std::vector<DiracSpinor> orbitals;
    const auto states_list = AtomData::listOfStates_nk(states);
    for (const auto &[n, k, en] : states_list) {
      auto &Fnk = orbitals.emplace_back(n, k, grid);
      // Use non-rel formula for guess (alpha = 0.0 gives non-rel)
      const auto en_guess = -(Zeff * Zeff) / (2.0 * n * n);
      DiracODE::boundState(Fnk, en_guess, v_nuc, {}, PhysConst::alpha,
                           converge_target);
    }

    // In the following, we find the the _worst_ orbital (by means of comparison
    // to the expected exact Dirac equation solution) for a number of properties,
    // and check if it meets the criteria.

    { // Check convergence:
      const auto comp_eps = [](const auto &Fa, const auto &Fb) {
        return Fa.eps() < Fb.eps();
      };
      const auto worst_F =
          std::max_element(cbegin(orbitals), cend(orbitals), comp_eps);

      std::cout << "converge: " << worst_F->shortSymbol() << " "
                << worst_F->eps() << "\n";
      CHECK(std::abs(worst_F->eps()) < converge_target);
    }

    { // Check orthogonality of orbitals:
      const auto [eps, worst] = DiracSpinor::check_ortho(orbitals, orbitals);
      std::cout << "orth: " << worst << " " << eps << "\n";
      const auto orthog_targ = Zeff >= 1.0 ? 1.0e-10 : 1.0e-8;
      CHECK(std::abs(eps) < orthog_targ);
    }

    { // Compare energy to exact (Dirac) value:
      auto comp_eps_en = [Zeff](const auto &Fa, const auto &Fb) {
        const auto exact_a =
            AtomData::diracen(Zeff, Fa.n(), Fa.kappa(), PhysConst::alpha);
        const auto exact_b =
            AtomData::diracen(Zeff, Fb.n(), Fb.kappa(), PhysConst::alpha);
        const auto eps_a = std::abs((Fa.en() - exact_a) / exact_a);
        const auto eps_b = std::abs((Fb.en() - exact_b) / exact_b);
        return eps_a < eps_b;
      };

      const auto worst_F =
          std::max_element(cbegin(orbitals), cend(orbitals), comp_eps_en);

      const auto exact = AtomData::diracen(Zeff, worst_F->n(), worst_F->kappa(),
                                           PhysConst::alpha);
      const auto eps = std::abs((worst_F->en() - exact) / exact);

      std::cout << "en vs. exact (eps): " << worst_F->shortSymbol() << " "
                << eps << "\n";
      const auto en_targ = Zeff >= 1.0 ? 1.0e-9 : 1.0e-7;
      CHECK(std::abs(eps) < en_targ);
    }

    { // Check radial integrals (r, r^2, 1/r, 1/r^2)

      // Define four radial operators. Designed to test wavefunction at low,
      // medium, and large radial distances
      const auto rhat1 = DiracOperator::RadialF(*grid, 1);
      const auto rhat2 = DiracOperator::RadialF(*grid, 2);
      const auto rinv1 = DiracOperator::RadialF(*grid, -1);
      const auto rinv2 = DiracOperator::RadialF(*grid, -2);

      // Lambda: finds worst comparison of <a|o|a> to <A|o|A>
      // |A> is exact orbital, |a> is solution from DiracODE
      // Returns pair
      const auto get_worst = [&orbitals, &grid, Zeff](const auto &o) {
        std::pair<std::string, double> worst{"", 0.0};
        for (const auto &Fa : orbitals) {
          const auto Fexact =
              DiracSpinor::exactHlike(Fa.n(), Fa.kappa(), grid, Zeff);
          const auto aoa = o.radialIntegral(Fa, Fa);
          const auto AoA = o.radialIntegral(Fexact, Fexact);
          const auto eps = std::abs((aoa - AoA) / AoA);
          if (eps > worst.second) {
            worst.first = Fa.shortSymbol();
            worst.second = eps;
          }
        }
        return worst;
      };

      const auto worst1 = get_worst(rhat1);
      const auto worst2 = get_worst(rhat2);
      const auto winv1 = get_worst(rinv1);
      const auto winv2 = get_worst(rinv2);

      std::cout << "<r>: " << worst1.first << " " << worst1.second << "\n";
      std::cout << "<r^2>: " << worst2.first << " " << worst2.second << "\n";
      std::cout << "<r^-1>: " << winv1.first << " " << winv1.second << "\n";
      std::cout << "<r^-2>: " << winv2.first << " " << winv2.second << "\n";

      auto eps_targ = Zeff >= 1.0 ? 1.0e-9 : 1.0e-6;
      if (Zeff >= 100) {
        eps_targ *= 10.0;
      }
      CHECK(std::abs(worst1.second) < eps_targ);
      CHECK(std::abs(worst2.second) < eps_targ);
      CHECK(std::abs(winv1.second) < eps_targ);
      CHECK(std::abs(winv2.second) < eps_targ);
    }
  }
}

//==============================================================================
// Test inhomogenous (Green's) method:
TEST_CASE("DiracODE: inhomogenous (Green's) method", "[DiracODE][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "DiracODE: inhomogenous (Green's) method\n";

  const double Zeff = 5.0;

  // Set up radial grid:
  const auto r0{1.0e-7};
  const auto rmax{100.0}; // NB: rmax depends on Zeff
  const auto num_grid_points{2000ul};
  const auto b{10.0};
  const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                 GridType::loglinear, b);
  // Sperical potential w/ R_nuc = 0.0 is a pointlike potential
  const auto v_nuc = Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid->r());

  {
    // Solve: (Fa and Fb should be equal)
    // (H + v + vp - e)Fa = 0
    // (H + v - e)Fb = -vp*Fa
    std::cout << "Test inhomogenous (Green's) method:   \n"
                 "(H + v + vp - e)Fa = 0     vs.\n"
                 "(H + v - e)Fb = -vp*Fa    (Fa and Fb should be equal)\n";
    auto max_eps_dF = -1.0;
    auto max_eps_orthNorm = -1.0;
    const auto states_new = AtomData::listOfStates_nk("5spdf");

    // This will act as a "non-local" potential
    std::vector<double> vp;
    for (const auto r : grid->r()) {
      vp.push_back(-0.3 / (r * r * r * r + 1.0));
    }
    const auto v_tot = qip::add(v_nuc, vp);

    for (const auto &[n, k, en] : states_new) {

      auto Fa = DiracSpinor(n, k, grid);
      auto Fap1 = DiracSpinor(n + 1, k, grid); // for orthog
      auto Fb = DiracSpinor(n, k, grid);
      const auto en_guess = -(Zeff * Zeff) / (2.0 * n * n);
      const auto en_guess_p1 = -(Zeff * Zeff) / (2.0 * (n + 1) * (n + 1));

      DiracODE::boundState(Fa, en_guess, v_tot, {}, PhysConst::alpha, 1.0e-15);
      DiracODE::boundState(Fap1, en_guess_p1, v_tot, {}, PhysConst::alpha,
                           1.0e-15);

      const auto dvFa = vp * Fa; // "non-local"
      DiracODE::solve_inhomog(Fb, Fa.en(), v_nuc, {}, PhysConst::alpha,
                              -1 * dvFa);

      const auto eps_norm = std::abs(Fb * Fb - 1.0); //<b|b> - norm

      Fb.normalise(); // don't propogate norm error:

      const auto eps_orth = std::abs(Fap1 * Fb);  //<a+1|b> - orthogonality
      const auto eps_1 = std::abs(Fa * Fb - 1.0); // <a|b> - check values
      const auto eps_2 = (Fa - Fb) * (Fa - Fb);   // <a-b>^2 - check values

      const auto eps_dF = std::max(eps_1, eps_2);
      const auto eps_orthNorm = std::max(eps_norm, eps_orth);
      printf("%4s <b|b>-1 = %.0e, <a|b>-1 = %.0e, <a-b>^2 = %.0e, <a+1|b> = "
             "%.0e\n",
             Fa.shortSymbol().c_str(), eps_norm, eps_1, eps_2, eps_orth);

      if (!(eps_dF < max_eps_dF))
        max_eps_dF = eps_dF;
      if (!(eps_orthNorm < max_eps_orthNorm))
        max_eps_orthNorm = eps_orthNorm;
    }
    //"Inhomog (G): orthonorm"
    REQUIRE(std::abs(max_eps_orthNorm) < 1.0e-5);
    //"Inhomog (G): value"
    REQUIRE(std::abs(max_eps_dF) < 1.0e-11);
  }

  { // Test DiracODE HartreeFock method:
    // Solve: (Fa and Fb should be equal)
    // (H + v + vp - e)Fa = 0
    // (H + v - e)Fb + X = 0, where X = VpFa (for now, local Vp)
    std::cout << "Test DiracODE HartreeFock method:   \n"
                 "(H + v + vx - e)Fa = 0     vs.\n"
                 "(H + v - e)Fb = -vx*Fa    (Fa and Fb should be equal)\n";
    auto max_eps_dF = -1.0;
    auto max_eps_en = -1.0;
    auto max_eps_orthNorm = -1.0;
    const auto states_new = AtomData::listOfStates_nk("5spdf");

    // This will act as a "non-local" potential
    std::vector<double> vp;
    for (const auto r : grid->r()) {
      vp.push_back(0.1 / (r * r * r * r + 1.0));
    }
    const auto v_tot = qip::add(v_nuc, vp);

    for (const auto &[n, k, en] : states_new) {

      auto Fa = DiracSpinor(n, k, grid);
      auto Fap1 = DiracSpinor(n + 1, k, grid); // for orthog
      auto Fb = DiracSpinor(n, k, grid);
      const auto en_guess = -(Zeff * Zeff) / (2.0 * n * n);
      const auto en_guess_p1 = -(Zeff * Zeff) / (2.0 * (n + 1) * (n + 1));

      // Solve 'a' version (local)
      DiracODE::boundState(Fa, en_guess, v_tot, {}, PhysConst::alpha, 1.0e-15);
      DiracODE::boundState(Fap1, en_guess_p1, v_tot, {}, PhysConst::alpha,
                           1.0e-15);

      auto dvFa = vp * Fa;
      DiracODE::boundState(Fb, Fa.en(), v_nuc, {}, PhysConst::alpha, 1.0e-15,
                           &dvFa, &Fa, 1);

      const auto eps_norm = std::abs(Fb * Fb - 1.0); //<b|b> - norm

      Fb.normalise(); // don't propogate norm error:

      const auto eps_en = std::abs((Fb.en() - Fa.en()) / (Fa.en()));
      const auto eps_orth = std::abs(Fap1 * Fb);  //<a+1|b> - orthogonality
      const auto eps_1 = std::abs(Fa * Fb - 1.0); // <a|b> - check values
      const auto eps_2 = (Fa - Fb) * (Fa - Fb);   // <a-b>^2 - check values

      const auto eps_dF = std::max(eps_1, eps_2);
      const auto eps_orthNorm = std::max(eps_norm, eps_orth);
      printf("%4s <b|b>-1 = %.0e, <a|b>-1 = %.0e, <a-b>^2 "
             "= %.0e, <a+1|b> = "
             "%.0e, en:%.0e\n",
             Fa.shortSymbol().c_str(), eps_norm, eps_1, eps_2, eps_orth,
             eps_en);

      if (!(eps_dF < max_eps_dF))
        max_eps_dF = eps_dF;
      if (!(eps_orthNorm < max_eps_orthNorm))
        max_eps_orthNorm = eps_orthNorm;
      if (!(eps_en < max_eps_en))
        max_eps_en = eps_en;
    }
    //"Dirac ODE HF: orthonorm"
    REQUIRE(std::abs(max_eps_orthNorm) < 5.0e-9);
    //"Dirac ODE HF: value"
    REQUIRE(std::abs(max_eps_dF) < 1.0e-14);
    //"Dirac ODE HF: energy"
    REQUIRE(std::abs(max_eps_en) < 5.0e-9);
  }
}

//==============================================================================
TEST_CASE("DiracODE: continuum", "[DiracODE][cntm][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "DiracODE: continuum\n";

  // Data from Mathematica: using exact H-like functions (non-relativistic)
  // [z,l,n,e,<Fn|r|Fe>,<Fn|1/r|Fe>]
  const auto data = std::vector{std::tuple{1, 0, 1, 0.01, -1.05444, 0.537766},
                                {5, 0, 1, 0.01, -0.0432611, 0.541197},
                                {10, 0, 1, 0.01, -0.0108239, 0.541305},
                                {1, 1, 2, 0.01, -3.39833, 0.229388},
                                {5, 1, 2, 0.01, -0.152383, 0.238861},
                                {10, 1, 2, 0.01, -0.0382367, 0.239171},
                                {1, 2, 3, 0.01, -5.27605, 0.115292},
                                {5, 2, 3, 0.01, -0.278067, 0.129662},
                                {10, 2, 3, 0.01, -0.0701639, 0.130167},
                                {1, 0, 2, 0.01, -3.70084, 0.302915},
                                {5, 0, 2, 0.01, -0.164999, 0.310496},
                                {10, 0, 2, 0.01, -0.0413950, 0.310744},
                                {1, 1, 3, 0.01, -6.76090, 0.180263},
                                {5, 1, 3, 0.01, -0.345874, 0.193675},
                                {10, 1, 3, 0.01, -0.0871896, 0.194138},
                                {5, 2, 4, 0.01, -0.545582, 0.127082},
                                {10, 2, 4, 0.01, -0.138554, 0.127784},
                                {1, 0, 3, 0.01, -6.96622, 0.211203},
                                {5, 0, 3, 0.01, -0.354696, 0.222719},
                                {10, 0, 3, 0.01, -0.0894003, 0.223119},
                                {5, 1, 4, 0.01, -0.598245, 0.159932},
                                {10, 1, 4, 0.01, -0.151798, 0.160551},
                                {5, 2, 5, 0.01, -0.860366, 0.116040},
                                {10, 2, 5, 0.01, -0.220328, 0.116914},
                                {1, 0, 1, 0.1, -0.847141, 0.508285},
                                {5, 0, 1, 0.1, -0.0428494, 0.539903},
                                {10, 0, 1, 0.1, -0.0107980, 0.540981},
                                {1, 1, 2, 0.1, -1.50201, 0.168976},
                                {5, 1, 2, 0.1, -0.145865, 0.235207},
                                {10, 1, 2, 0.1, -0.0378163, 0.238243},
                                {1, 2, 3, 0.1, -1.05930, 0.0549268},
                                {5, 2, 3, 0.1, -0.249619, 0.123885},
                                {10, 2, 3, 0.1, -0.0682495, 0.128663},
                                {1, 0, 2, 0.1, -1.72179, 0.253816},
                                {5, 0, 2, 0.1, -0.158282, 0.307574},
                                {10, 0, 2, 0.1, -0.0409619, 0.310002},
                                {1, 1, 3, 0.1, -1.68704, 0.115957},
                                {5, 1, 3, 0.1, -0.314038, 0.188336},
                                {10, 1, 3, 0.1, -0.0850546, 0.192756},
                                {5, 2, 4, 0.1, -0.456319, 0.119275},
                                {10, 2, 4, 0.1, -0.132231, 0.125703},
                                {1, 0, 3, 0.1, -1.81834, 0.156790},
                                {5, 0, 3, 0.1, -0.322617, 0.218126},
                                {10, 0, 3, 0.1, -0.0872496, 0.221928},
                                {5, 1, 4, 0.1, -0.505526, 0.153021},
                                {10, 1, 4, 0.1, -0.145244, 0.158714},
                                {5, 2, 5, 0.1, -0.659711, 0.106668},
                                {10, 2, 5, 0.1, -0.205188, 0.114339},
                                {1, 0, 1, 1., -0.231566, 0.347350},
                                {1, 1, 2, 1., -0.0876112, 0.0492813},
                                {5, 1, 2, 1., -0.0991548, 0.204507},
                                {10, 1, 2, 1., -0.0339833, 0.229388},
                                {5, 2, 3, 1., -0.106991, 0.0851966},
                                {10, 2, 3, 1., -0.0527605, 0.115292},
                                {1, 0, 2, 1., -0.142368, 0.135524},
                                {5, 0, 2, 1., -0.109897, 0.282901},
                                {10, 0, 2, 1., -0.0370084, 0.302915},
                                {5, 1, 3, 1., -0.149001, 0.150338},
                                {10, 1, 3, 1., -0.0676090, 0.180263},
                                {5, 2, 4, 1., -0.135334, 0.0747194},
                                {10, 2, 4, 1., -0.0876013, 0.108334},
                                {1, 0, 3, 1., -0.0876221, 0.0753679},
                                {5, 0, 3, 1., -0.155845, 0.185765},
                                {10, 0, 3, 1., -0.0696622, 0.211203},
                                {5, 1, 4, 1., -0.164807, 0.112249},
                                {10, 1, 4, 1., -0.0986910, 0.143264},
                                {5, 2, 5, 1., -0.138528, 0.0614803},
                                {10, 2, 5, 1., -0.113942, 0.0944541},
                                {1, 0, 1, 10., -0.0114138, 0.119845},
                                {5, 0, 1, 10., -0.0193264, 0.434844},
                                {1, 1, 2, 10., -0.00108025, 0.00546877},
                                {5, 1, 2, 10., -0.0141231, 0.0926831},
                                {10, 1, 2, 10., -0.0150201, 0.168976},
                                {5, 2, 3, 10., -0.00457608, 0.0173721},
                                {10, 2, 3, 10., -0.0105930, 0.0549268},
                                {1, 0, 2, 10., -0.00433815, 0.0428903},
                                {5, 0, 2, 10., -0.0186887, 0.185463},
                                {10, 0, 2, 10., -0.0172179, 0.253816},
                                {5, 1, 3, 10., -0.0107625, 0.0581273},
                                {10, 1, 3, 10., -0.0168704, 0.115957},
                                {5, 2, 4, 10., -0.00389425, 0.0137832},
                                {10, 2, 4, 10., -0.0109358, 0.0457548},
                                {1, 0, 3, 10., -0.00239395, 0.0234002},
                                {5, 0, 3, 10., -0.0132774, 0.105927},
                                {10, 0, 3, 10., -0.0181834, 0.156790},
                                {5, 1, 4, 10., -0.00791567, 0.0396239},
                                {10, 1, 4, 10., -0.0151507, 0.0826362},
                                {5, 2, 5, 10., -0.00312208, 0.0106815},
                                {10, 2, 5, 10., -0.00975771, 0.0364069},
                                {1, 0, 1, 100., -0.000267412, 0.0268749},
                                {5, 0, 1, 100., -0.00175216, 0.197118},
                                {10, 0, 1, 100., -0.00231566, 0.347350},
                                {1, 1, 2, 100., -7.74840e-6, 0.000387904},
                                {5, 1, 2, 100., -0.000275495, 0.0142052},
                                {10, 1, 2, 100., -0.000876112, 0.0492813},
                                {5, 2, 3, 100., -0.0000261183, 0.000882703},
                                {10, 2, 3, 100., -0.000171345, 0.00602882},
                                {1, 0, 2, 100., -0.0000952534, 0.00951355},
                                {5, 0, 2, 100., -0.000735639, 0.0717620},
                                {10, 0, 2, 100., -0.00142368, 0.135524},
                                {5, 1, 3, 100., -0.000168313, 0.00847526},
                                {10, 1, 3, 100., -0.000581805, 0.0299574},
                                {5, 2, 4, 100., -0.0000204521, 0.000685889},
                                {10, 2, 4, 100., -0.000138197, 0.00471937},
                                {5, 0, 3, 100., -0.000414090, 0.0392846},
                                {10, 0, 3, 100., -0.000876221, 0.0753679},
                                {5, 1, 4, 100., -0.000113492, 0.00566696},
                                {10, 1, 4, 100., -0.000404630, 0.0201691},
                                {5, 2, 5, 100., -0.0000157475, 0.000526220},
                                {10, 2, 5, 100., -0.000107907, 0.00363341}};

  const auto r0{1.0e-6};
  const auto rmax{100.0}; // NB: rmax depends on Zeff
  const auto num_grid_points{10000ul};
  const auto b{10.0};
  const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                 GridType::loglinear, b);

  // Find 'inital guess' for asymptotic region:

  double worst = 0.0;

  std::cout << "\nTest continuum solutions. Compare against exact "
               "(non-relativistic "
               "formula) from Mathematica (nb: MMA solutions not perfect)\n";
  std::cout << "Compare <F_nl|r|F_el> and <F_nl|1/r|F_el>:\n\n";

  std::cout << " Z  n  l   En    <Fn|r|Fe>             <Fn|1/r|Fe>\n";
  for (auto &[z, l, n, e, e1, e3] : data) {

    int kappa = -(l + 1);

    const auto v0 = Nuclear::sphericalNuclearPotential(z, 0.0, grid->r());

    DiracSpinor Fe{0, kappa, grid};

    DiracODE::solveContinuum(Fe, e, v0, 1.0e-10 * PhysConst::alpha);

    // const auto F1s =
    //     DiracSpinor::exactHlike(n, kappa, grid, z, 1.0e-10 * PhysConst::alpha);
    // Above fails on M1 mac, leads to 'nan'. Probably alpha is too small
    const auto F1s =
        DiracSpinor::exactHlike(n, kappa, grid, z, 1.0e-8 * PhysConst::alpha);

    auto v1 = Fe * (grid->r() * F1s);
    auto v3 = Fe * (grid->rpow(-1) * F1s);

    auto eps1 = std::min(std::abs((v1 - e1) / e1), std::abs(v1 - e1));
    auto eps3 = std::min(std::abs((v3 - e3) / e3), std::abs(v3 - e3));

    printf("%2i %2i %2i %6.2f %+7.5f [%+7.5f]   ", z, n, l, e, v1, e1);
    auto eps = std::max(eps1, eps3);

    printf("%+7.5f [%+7.5f] - %.0e [%.0e]", v3, e3, eps, Fe.eps());

    if (eps > worst) {
      worst = eps;
      if (eps > 0.01)
        std::cout << " *";
    }
    std::cout << "\n";

    if (e > 0.05) {
      REQUIRE(eps < 1.0e-2);
    } else {
      REQUIRE(eps < 1.0e-1);
    }
  }
  std::cout << worst << "\n";
}
