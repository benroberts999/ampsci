#pragma once
#include "DiracODE/DiracODE.hpp"
#include "DiracOperator/Operators.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Check.hpp"
#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace UnitTest {

//******************************************************************************
//! Unit tests for solving (local) Dirac equation ODE
bool DiracODE(std::ostream &obuff) {
  bool pass = true;

  const double Zeff = 5.0;

  // Set up radial grid:
  const auto r0{1.0e-7};
  const auto rmax{100.0}; // NB: rmax depends on Zeff
  const auto num_grid_points{2000ul};
  const auto b{10.0};
  const auto grid = std::make_shared<const Grid>(r0, rmax, num_grid_points,
                                                 GridType::loglinear, b);

  // States to solve for:
  const std::string states = "10spdfghi";

  // Sperical potential w/ R_nuc = 0.0 is a pointlike potential
  const auto v_nuc = Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid->r());

  // Solve Dirac ODE for each state, store in 'orbitals' vector:
  std::vector<DiracSpinor> orbitals;
  const auto states_list = AtomData::listOfStates_nk(states);
  for (const auto &[n, k, en] : states_list) {
    auto &Fnk = orbitals.emplace_back(n, k, grid);
    // Use non-rel formula for guess (alpha = 0.0 gives non-rel)
    const auto en_guess = -(Zeff * Zeff) / (2.0 * n * n);
    DiracODE::boundState(Fnk, en_guess, v_nuc, {}, PhysConst::alpha, 15);
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
    pass &= qip::check_value(&obuff, "converge " + worst_F->shortSymbol(),
                             worst_F->eps(), 0.0, 1.0e-14);
  }

  { // Check orthogonality of orbitals:
    const auto [eps, worst] = DiracSpinor::check_ortho(orbitals, orbitals);
    pass &= qip::check_value(&obuff, "orth " + worst, eps, 0.0, 1.0e-10);
  }

  { // Compare energy to exact (Dirac) value:
    auto comp_eps_en = [Zeff](const auto &Fa, const auto &Fb) {
      const auto exact_a =
          AtomData::diracen(Zeff, Fa.n, Fa.k, PhysConst::alpha);
      const auto exact_b =
          AtomData::diracen(Zeff, Fb.n, Fb.k, PhysConst::alpha);
      const auto eps_a = std::abs((Fa.en() - exact_a) / exact_a);
      const auto eps_b = std::abs((Fb.en() - exact_b) / exact_b);
      return eps_a < eps_b;
    };

    const auto worst_F =
        std::max_element(cbegin(orbitals), cend(orbitals), comp_eps_en);

    const auto exact =
        AtomData::diracen(Zeff, worst_F->n, worst_F->k, PhysConst::alpha);
    const auto eps = std::abs((worst_F->en() - exact) / exact);

    pass &=
        qip::check_value(&obuff, "en vs. exact (eps) " + worst_F->shortSymbol(),
                         eps, 0.0, 1.0e-10);
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
        const auto Fexact = DiracSpinor::exactHlike(Fa.n, Fa.k, grid, Zeff);
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

    pass &= qip::check_value(&obuff, "<r>    " + worst1.first, worst1.second,
                             0.0, 1.0e-10);
    pass &= qip::check_value(&obuff, "<r^2>  " + worst2.first, worst2.second,
                             0.0, 1.0e-10);
    pass &= qip::check_value(&obuff, "<r^-1> " + winv1.first, winv1.second, 0.0,
                             1.0e-10);
    pass &= qip::check_value(&obuff, "<r^-2> " + winv2.first, winv2.second, 0.0,
                             1.0e-10);
  }

  { // Test inhomogenous (Green's) method:
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

      DiracODE::boundState(Fa, en_guess, v_tot, {}, PhysConst::alpha, 15);
      DiracODE::boundState(Fap1, en_guess_p1, v_tot, {}, PhysConst::alpha, 15);

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
    pass &= qip::check_value(&obuff, "Inhomog (G): orthonorm", max_eps_orthNorm,
                             0.0, 1.0e-5);
    pass &= qip::check_value(&obuff, "Inhomog (G): value", max_eps_dF, 0.0,
                             1.0e-11);
  }

  return pass;
}

} // namespace UnitTest
