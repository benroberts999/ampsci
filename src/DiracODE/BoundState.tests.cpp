#include "Angular/include.hpp"
#include "DiracODE/AsymptoticSpinor.hpp"
#include "DiracODE/include.hpp"
#include "DiracOperator/include.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/DiracContinuum.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
//! Unit tests for solving (local) Dirac equation ODE
TEST_CASE("DiracODE: Low-r solution", "[DiracODE][bound][unit]") {

  std::cout << "Low-r solution (pointlike)\n";

  // Set up radial grid:
  const auto r0{1.0e-6};
  const auto rmax{100.0}; // NB: rmax depends on Zeff
  const auto num_grid_points{2000ul};
  const auto b{10.0};
  const auto grid = Grid(r0, rmax, num_grid_points, GridType::loglinear, b);

  bool print_table = false;

  if (print_table) {
    std::cout << " Zeff  k n ri r        (f/g)_r          [ Exact           ]  "
                 "eps     [targ]\n";
  }
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
        const auto F1s =
          DiracSpinor::exactHlike(n, kappa, std::make_shared<Grid>(grid), Zeff);

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
        magnitude > 1e2 ? 1.0e-6 : 1.0e-5;
      // magnitude > 1e0  ? 1.0e-4 :
      // magnitude > 1e-1 ? 1.0e-3 :
      // magnitude > 1e-2 ? 1.0e-2 :
      //                    1.0e-1;
      if (Zeff < 1.0)
        target *= 6.0;

      if (print_table) {
        fmt::print(
          "{:5g} {:2} {} {} {:.1e} {:+.10e} [{:+.10e}]  {:.0e} [{:.0e}]\n",
          Zeff, kappa, w_n, w_i, w_r, w_asym, w_exact, w_eps, target);
      }

      REQUIRE(w_eps < target);
    }
  }
}

//==============================================================================
//! Unit tests for solving (local) Dirac equation ODE
TEST_CASE("DiracODE: Adams-Moulton method", "[DiracODE][bound][unit]") {

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
TEST_CASE("DiracODE: exotic atoms - unit", "[Exotic][bound][unit]") {
  // Simple convergence check: DiracODE::boundState converges for exotic
  // (muonic-like) atoms with non-unit mass, for both m=1 and m=m_muon

  const double Z = 1.0;
  const double alpha = PhysConst::alpha;

  for (const double mass : {1.0, PhysConst::m_muon}) {
    const double r0 = 1.0e-5 / Z / mass;
    const double rmax = 1.0e2 / Z / mass;
    const double b = rmax / 10.0;

    auto grid = std::make_shared<const Grid>(
      Grid{r0, rmax, 1000ul, GridType::loglinear, b});

    const auto V0 = Nuclear::sphericalNuclearPotential(Z, 0.0, grid->r());

    for (const auto [n, kappa, x_en] : AtomData::listOfStates_nk("2sp")) {
      const auto e_excat = mass * AtomData::diracen(Z, n, kappa, alpha);
      const auto e0 = 0.85 * e_excat;
      const auto F = DiracODE::boundState(n, kappa, e0, grid, V0, {}, alpha,
                                          1.0e-14, nullptr, nullptr, Z, mass);
      REQUIRE(F.eps() < 1.0e-12);
      REQUIRE(F.en() == Approx(e_excat).epsilon(1.0e-9));
    }
  }
}

//==============================================================================
TEST_CASE("DiracODE: exotic atoms - numerical", "[Exotic][bound]") {
  // Numerical test: energies for pointlike H-like atom match the analytic
  // Dirac formula E = m * diracen(Z, n, kappa, alpha) for a range of Z and
  // masses from m=1 (electron) up to m=m_muon.

  const double alpha = PhysConst::alpha;

  const std::vector<int> Zs = {1, 10, 50, 100};
  const std::vector<double> masses = {0.1, 10.0, PhysConst::m_muon, 500.0};

  int count = 0;
  for (const int Z : Zs) {
    for (const double mass : masses) {
      fmt::print("\nZ = {:3}, mass = {:.4f} au = {:.4f} MeV\n", Z, mass,
                 mass * PhysConst::m_e_MeV);

      const double r0 = std::max(1.0e-5 / Z / mass, 1.0e-9);
      const double rmax = std::max(200.0 / Z / mass, 1.0);
      const double b = rmax / 10.0;
      std::cout << "Grid: [" << r0 << ", " << rmax << "] au\n";

      auto grid = std::make_shared<const Grid>(
        Grid{r0, rmax, 5000ul, GridType::loglinear, b});

      const auto V0 = Nuclear::sphericalNuclearPotential(Z, 0.0, grid->r());

      fmt::print("{:6s}  {:>15s}  Rinf    (conv )  [{:>15s}]  {:<8s}\n",
                 "state", "E (au)", "expected", "eps");
      for (const auto [n, kappa, x_en] : AtomData::listOfStates_nk("2sp3d4f")) {
        const auto e_exact = mass * AtomData::diracen(Z, n, kappa, alpha);
        const auto e0 = e_exact;

        const auto F = DiracODE::boundState(n, kappa, e0, grid, V0, {}, alpha,
                                            1.0e-14, nullptr, nullptr, Z, mass);

        const auto rel_err = std::abs(F.en() / e_exact - 1.0);
        std::string warning = F.eps() > 1.0e-9 ? "**" : "";
        fmt::print("{:6s}  {:15.8e}  {:.1e} ({:.0e})  [{:15.8e}]  {:.1e}  {}\n",
                   F.shortSymbol(), F.en(), F.rinf(), F.eps(), e_exact, rel_err,
                   warning);
        if (F.eps() > 1.0e-9)
          ++count;

        // If Dirac Eq. converged, ensure energy matches to high degree
        const auto target = mass > 1.0 ? 1.0e-10 : 1.0e-6;
        if (F.eps() < 1.0e-9)
          REQUIRE(F.en() == Approx(e_exact).epsilon(target));
      }
    }
  }
  std::cout << "\nNot converged: " << count << "\n";
  REQUIRE(count < 20);
  std::cout << "(It's OK if not all converged, so long as we noticed with conv "
               "parameter)\n";
}

//==============================================================================
TEST_CASE("Wavefunction: H-like Exotic muon and tauon",
          "[DiracODE][bound][unit][Exotic][muon][tauon]") {

  std::cout << "\nSolve H-like muonic and tauonic atoms, for few Z, \n"
               "compare energies, matrix elements, and asymptotics to exact "
               "Dirac functions\n\n";

  for (const int Z : {2, 55}) {
    for (const double m_mu : {PhysConst::m_muon, PhysConst::m_tau}) {

      const std::string ident =
        m_mu == Approx(PhysConst::m_tau) ? "tauonic" : "muonic";
      const auto label =
        fmt::format("Z = {} - {} (m={:.1f} m_e)", Z, ident, m_mu);

      // Grid scaled for exotic orbits (~m times smaller than electronic)
      Wavefunction wf(
        {2000, 1.0e-5 / m_mu, 100.0 / m_mu / Z, 10.0 / m_mu, "loglinear"},
        {std::to_string(Z), -1, "pointlike", -1.0, -1.0});

      wf.solve_core("HartreeFock", "[]");
      wf.solve_exotic("4sp3d4f", m_mu, false);

      fmt::print("\n{}: solve ODE\n", label);
      double temp_worst = 0.0;
      for (const auto &Fv : wf.valence()) {
        temp_worst = std::max(temp_worst, Fv.eps());
        CHECK(Fv.eps() < 1.0e-12);
      }
      std::cout << "Worst eps: " << temp_worst << "\n";

      // Exact exotic H-like orbitals (with exotic-particle mass)
      const auto hlike = DiracSpinor::HlikeBasis("4sp3d4f", wf.grid_sptr(), Z,
                                                 PhysConst::alpha, m_mu);

      const auto &rv = wf.grid().r();

      // Energies
      fmt::print("\n{}: energies (valence vs exact):\n", label);
      fmt::print("{:<4} [{:>15}] {:>15} {:>7}\n", "State", "Exact", "E_val",
                 "eps");
      for (const auto &Fv : wf.valence()) {
        const auto *pFh = DiracSpinor::find(Fv.n(), Fv.kappa(), hlike);
        if (!pFh)
          continue;
        const auto eps = std::abs((Fv.en() - pFh->en()) / pFh->en());
        fmt::print("{:<4} [{:>15.8f}] {:>15.8f} {:>7.0e}\n", Fv.shortSymbol(),
                   pFh->en(), Fv.en(), eps);
        const auto target = Fv.twoj() <= 1 ? 1.0e-9 : 1.0e-12;
        REQUIRE(eps < target);
      }

      // Orthonormality: exact H-like basis
      fmt::print("\n{}: orthonormality <a|b> (exact):\n", label);
      double e_worst = 0.0;
      for (const auto &Fa : hlike) {
        for (const auto &Fb : hlike) {
          if (Fa.kappa() != Fb.kappa() || Fb.n() > Fa.n())
            continue;
          const double expected = (Fa.n() == Fb.n()) ? 1.0 : 0.0;
          const double val = Fa * Fb;
          const double eps = std::abs(val - expected);
          e_worst = std::max(eps, e_worst);
          CHECK(eps < 1.0e-9);
        }
      }
      fmt::print("Worst: {:.2e}\n", e_worst);

      // Orthonormality: exact H-like basis
      fmt::print("\n{}: orthonormality <a|b> (valence):\n", label);
      e_worst = 0.0;
      for (const auto &Fa : wf.valence()) {
        for (const auto &Fb : wf.valence()) {
          if (Fa.kappa() != Fb.kappa() || Fb.n() > Fa.n())
            continue;
          const double expected = (Fa.n() == Fb.n()) ? 1.0 : 0.0;
          const double val = Fa * Fb;
          const double eps = std::abs(val - expected);
          e_worst = std::max(eps, e_worst);
          CHECK(eps < 1.0e-9);
        }
      }
      fmt::print("Worst: {:.2e}\n", e_worst);

      // <r> diagonal
      fmt::print("\n{}: <r> diagonal:\n", label);
      fmt::print("{:<10} [{:>13}] {:>13} {:>7}\n", "State", "Exact", "<r>_val",
                 "eps");
      for (const auto &Fv : wf.valence()) {
        const auto *pFh = DiracSpinor::find(Fv.n(), Fv.kappa(), hlike);
        if (!pFh)
          continue;
        const double r_val = Fv * (rv * Fv);
        const double r_ex = (*pFh) * (rv * (*pFh));
        const double eps = std::abs((r_val - r_ex) / r_ex);
        fmt::print("{:<10} [{:>13.8f}] {:>13.8f} {:>7.0e}\n", Fv.shortSymbol(),
                   r_ex, r_val, eps);
        REQUIRE(eps < 1.0e-9);
      }

      // <n|r|n+1> off-diagonal
      fmt::print("\n{}: <n|r|n+1> off-diagonal:\n", label);
      fmt::print("{:<20} [{:>13}] {:>13} {:>7}\n", "States", "Exact", "<r>_val",
                 "eps");
      for (const auto &Fv : wf.valence()) {
        const auto *pFh = DiracSpinor::find(Fv.n(), Fv.kappa(), hlike);
        const auto *pFvb = wf.getState(Fv.n() + 1, Fv.kappa());
        const auto *pFhb = DiracSpinor::find(Fv.n() + 1, Fv.kappa(), hlike);
        if (!pFh || !pFvb || !pFhb)
          continue;
        const double r_val = Fv * (rv * (*pFvb));
        const double r_ex = (*pFh) * (rv * (*pFhb));
        const double eps = std::abs((r_val - r_ex) / r_ex);
        fmt::print(
          "{:<20} [{:>13.8f}] {:>13.8f} {:>7.0e}\n",
          fmt::format("<{}|r|{}>", Fv.shortSymbol(), pFvb->shortSymbol()), r_ex,
          r_val, eps);
        REQUIRE(eps < 1.0e-8);
      }

      // Small-r sign: f and g at r ~ 0.1/(Z*m_mu)
      const double r_small = 0.1 / (Z * m_mu);
      const auto ir = wf.grid().getIndex(r_small);
      fmt::print("\n{}: f,g at small r = {:.3e} (index {}):\n", label, rv[ir],
                 ir);
      fmt::print("{:<10} [{:>11}] {:>11} [{:>11}] {:>11}\n", "State", "f_ex",
                 "f_val", "g_ex", "g_val");
      for (const auto &Fv : wf.valence()) {
        const auto *pFh = DiracSpinor::find(Fv.n(), Fv.kappa(), hlike);
        if (!pFh)
          continue;
        fmt::print("{:<10} [{:>11.4e}] {:>11.4e} [{:>11.4e}] {:>11.4e}\n",
                   Fv.shortSymbol(), pFh->f(ir), Fv.f(ir), pFh->g(ir),
                   Fv.g(ir));

        REQUIRE(pFh->f(ir) == Approx(Fv.f(ir)).epsilon(1.0e-6));
        REQUIRE(pFh->g(ir) == Approx(Fv.g(ir)).epsilon(1.0e-5));
        REQUIRE(Fv.f(ir) > 0.0);
        REQUIRE(pFh->f(ir) > 0.0);

        // don't both print, but also test large
        const auto i_large = std::size_t(0.9 * double(Fv.max_pt()));
        REQUIRE(Fv.f(i_large) == Approx(pFh->f(i_large)).epsilon(1.0e-7));
      }
    }
  }
}

//==============================================================================
// Tests the large-r tail extension (DiracODE::Internal::extendTail). After the
// bound state is found, the decaying orbital is continued outward past the
// practical infinity until |f| drops below Param::tail_cut * max|f|.
TEST_CASE("DiracODE: tail extension", "[DiracODE][unit]") {
  std::cout << "Tail extension past practical infinity\n";

  using namespace DiracODE;
  const double Z = 20.0;
  const auto alpha = PhysConst::alpha;
  // Grid must extend well past the cALR cutoff so there is room to extend:
  const auto grid = std::make_shared<const Grid>(1.0e-7, 150.0, 4000ul,
                                                 GridType::loglinear, 10.0);
  const auto v = Nuclear::sphericalNuclearPotential(Z, 0.0, grid->r());
  const auto &r = grid->r();
  const auto num_points = grid->num_points();

  bool any_extended = false;
  for (const auto &[n, k, en_lst] : AtomData::listOfStates_nk("4spdf")) {
    (void)en_lst;
    auto F = DiracSpinor(n, k, grid);
    const auto en_guess = -0.5 * Z * Z / double(n * n);
    boundState(F, en_guess, v, {}, alpha, 1.0e-15);

    // Practical infinity *before* extension (same definition boundState uses):
    const auto old_pinf =
      Internal::findPracticalInfinity(F.en(), v, r, Internal::Param::cALR);
    const auto new_pinf = F.max_pt();

    double fmax = 0.0;
    for (const auto x : F.f())
      fmax = std::max(fmax, std::abs(x));
    REQUIRE(fmax > 0.0);

    // Energy still correct: extension must not corrupt the solution:
    const auto e_exact = DiracHydrogen::enk(n, k, Z, alpha);
    REQUIRE(F.en() == Approx(e_exact).epsilon(1.0e-6));
    // ...and still normalised:
    REQUIRE(F.norm() == Approx(1.0).epsilon(1.0e-12));

    // Tail was extended outward (there is room on this grid):
    REQUIRE(new_pinf >= old_pinf);
    if (new_pinf > old_pinf)
      any_extended = true;

    // Extended tail is the physical decaying solution: fixed sign, and
    // monotonically decreasing in |f| (i.e. no spurious growing solution):
    const auto sign0 = (F.f()[old_pinf - 1] >= 0.0) ? 1 : -1;
    for (auto i = old_pinf; i < new_pinf; ++i) {
      const auto fi = F.f(i);
      REQUIRE((fi == 0.0 || (fi > 0.0 ? 1 : -1) == sign0));
      REQUIRE(std::abs(fi) <= std::abs(F.f()[i - 1]));
    }

    // Tail reached the cutoff (unless it ran into the grid edge):
    if (new_pinf < num_points) {
      REQUIRE(std::abs(F.f()[new_pinf - 1]) <
              10.0 * Internal::Param::tail_cut * fmax);
    }

    // Everything beyond the new practical infinity is zeroed:
    for (auto i = new_pinf; i < num_points; ++i) {
      REQUIRE(F.f(i) == 0.0);
      REQUIRE(F.g(i) == 0.0);
    }

    // Reduce case: calling again with a *larger* cutoff must move pinf back
    // inward, to where |f| ~ the new (larger) cutoff:
    {
      const double big_cut = 1.0e-8;
      const double big_fcut = big_cut * fmax;
      const Internal::DiracDerivative Hd(*grid, v, k, F.en(), alpha);
      const auto reduced =
        Internal::extendTail(F.f(), F.g(), Hd, new_pinf, big_cut);
      REQUIRE(reduced < new_pinf);
      REQUIRE(reduced >= 2);
      REQUIRE(std::abs(F.f(reduced - 2)) >=
              0.1 * big_fcut); // last point above cutoff
      REQUIRE(std::abs(F.f(reduced - 1)) < 0.1 * big_fcut); // first point below
      REQUIRE(F.f(reduced) == 0.0); // trimmed tail is zeroed
    }
  }
  REQUIRE(any_extended);
}

//==============================================================================
//! Inhomogeneous solve: a local potential vx entered as the source term
//! VxFa = vx*F, iterated to self-consistency, must give the same solution as
//! including vx in the potential directly. Tests the particular + homogeneous
//! inward solution and its join to the outward solution.
TEST_CASE("DiracODE: inhomogeneous bound state", "[DiracODE][bound][unit]") {

  const auto grid = std::make_shared<const Grid>(1.0e-6, 150.0, 2000,
                                                 GridType::loglinear, 10.0);
  const double Zeff = 5.0;
  const auto v_nuc = Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid->r());
  // Short-range attractive potential, standing in for exchange/correlation
  std::vector<double> vx(grid->num_points());
  for (std::size_t i = 0; i < vx.size(); ++i) {
    vx[i] = -2.0 * std::exp(-grid->r(i));
  }
  const auto v_tot = qip::add(v_nuc, vx);

  double worst_eps = 0.0;
  for (int n = 1; n <= 6; ++n) {
    for (int l = 0; l < n && l <= 3; ++l) {
      for (const int kappa : {l, -l - 1}) {
        if (kappa == 0) {
          continue;
        }
        const auto en_guess = -0.5 * Zeff * Zeff / (n * n);
        const auto F_direct = DiracODE::boundState(
          n, kappa, en_guess, grid, v_tot, {}, PhysConst::alpha, 1.0e-15);

        auto Fa = DiracODE::boundState(n, kappa, en_guess, grid, v_nuc, {},
                                       PhysConst::alpha, 1.0e-15);
        for (int it = 0; it < 100; ++it) {
          const auto Fa_prev = Fa;
          const auto VxFa = vx * Fa_prev;
          DiracODE::boundState(Fa, Fa.en(), v_nuc, {}, PhysConst::alpha,
                               1.0e-15, &VxFa, &Fa_prev);
          if (std::abs((Fa.en() - Fa_prev.en()) / Fa.en()) < 1.0e-15)
            break;
        }

        REQUIRE(Fa.num_nodes() == n - l - 1);
        const auto eps = std::abs((Fa.en() - F_direct.en()) / F_direct.en());
        worst_eps = std::max(worst_eps, eps);
      }
    }
  }
  REQUIRE(worst_eps < 1.0e-12);
}
