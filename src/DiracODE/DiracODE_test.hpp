#pragma once
#include "DiracODE/DiracODE.hpp"
#include "DiracOperator/Operators.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Check.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace UnitTest {

//******************************************************************************
//! Unit tests for solveing (local) Dirac equation ODE
bool DiracODE(std::ostream &obuff) {
  bool pass = true;

  const double Zeff = 5.0;
  const double rmax = 100.0; // note: rmax depends on Zeff

  const auto grid =
      std::make_shared<const Grid>(1.0e-7, rmax, 2000, GridType::loglinear, 10);

  const auto v_nuc = [&]() {
    std::vector<double> v;
    for (const auto &r : grid->r) {
      v.push_back(-Zeff / r);
    }
    return v;
  }();

  const auto states = AtomData::listOfStates_nk("10spdfghi");

  std::vector<DiracSpinor> orbitals;
  for (const auto &[n, k, en] : states) {
    auto &Fnk = orbitals.emplace_back(n, k, grid);
    const auto en_guess =
        -(Zeff * Zeff) / (2.0 * n * n); // non-rel formula for guess
    DiracODE::boundState(Fnk, en_guess, v_nuc, {}, PhysConst::alpha, 15);
    // printf("%7s en=%11.5f rmax=%5.1f it=%3i eps=%.1e\n",
    // Fnk.symbol().c_str(),
    //        Fnk.en, Fnk.rinf(), Fnk.its, Fnk.eps);
  }

  { // Check convergence:
    const auto comp_eps = [](const auto &Fa, const auto &Fb) {
      return Fa.eps < Fb.eps;
    };
    const auto worst_F =
        std::max_element(cbegin(orbitals), cend(orbitals), comp_eps);
    pass &= qip::check_value(&obuff, "converge " + worst_F->shortSymbol(),
                             worst_F->eps, 0.0, 1.0e-14);
  }

  // Check orthogonality of orbitals:
  {
    const auto [eps, worst] = DiracSpinor::check_ortho(orbitals, orbitals);
    pass &= qip::check_value(&obuff, "orth " + worst, eps, 0.0, 1.0e-10);
  }

  { // Compare energy to exact (Dirac) value:
    auto comp_eps_en = [Zeff](const auto &Fa, const auto &Fb) {
      const auto exact_a =
          AtomData::diracen(Zeff, Fa.n, Fa.k, PhysConst::alpha);
      const auto exact_b =
          AtomData::diracen(Zeff, Fb.n, Fb.k, PhysConst::alpha);
      const auto eps_a = std::abs((Fa.en - exact_a) / exact_a);
      const auto eps_b = std::abs((Fb.en - exact_b) / exact_b);
      return eps_a < eps_b;
    };

    const auto worst_F =
        std::max_element(cbegin(orbitals), cend(orbitals), comp_eps_en);

    const auto exact =
        AtomData::diracen(Zeff, worst_F->n, worst_F->k, PhysConst::alpha);
    const auto eps = std::abs((worst_F->en - exact) / exact);

    pass &=
        qip::check_value(&obuff, "en vs. exact (eps) " + worst_F->shortSymbol(),
                         eps, 0.0, 1.0e-10);
  }

  { // Check radial integrals (r, r^2, 1/r, 1/r^2)
    const auto rhat1 = DiracOperator::RadialF(*grid, 1);
    const auto rhat2 = DiracOperator::RadialF(*grid, 2);
    const auto rinv1 = DiracOperator::RadialF(*grid, -1);
    const auto rinv2 = DiracOperator::RadialF(*grid, -2);

    std::pair<std::string, double> worst1{"", 0.0};
    std::pair<std::string, double> worst2{"", 0.0};
    std::pair<std::string, double> winv1{"", 0.0};
    std::pair<std::string, double> winv2{"", 0.0};

    for (const auto &Fa : orbitals) {
      const auto Fexact = DiracSpinor::exactHlike(Fa.n, Fa.k, grid, Zeff);

      const auto eps1 = std::abs((rhat1.radialIntegral(Fa, Fa) -
                                  rhat1.radialIntegral(Fexact, Fexact)) /
                                 rhat1.radialIntegral(Fa, Fa));
      const auto eps2 = std::abs((rhat2.radialIntegral(Fa, Fa) -
                                  rhat2.radialIntegral(Fexact, Fexact)) /
                                 rhat2.radialIntegral(Fa, Fa));
      const auto inv1 = std::abs((rinv1.radialIntegral(Fa, Fa) -
                                  rinv1.radialIntegral(Fexact, Fexact)) /
                                 rinv1.radialIntegral(Fa, Fa));
      const auto inv2 = std::abs((rinv2.radialIntegral(Fa, Fa) -
                                  rinv2.radialIntegral(Fexact, Fexact)) /
                                 rinv2.radialIntegral(Fa, Fa));

      if (eps1 > worst1.second) {
        worst1.first = Fa.shortSymbol();
        worst1.second = eps1;
      }
      if (eps2 > worst2.second) {
        worst2.first = Fa.shortSymbol();
        worst2.second = eps2;
      }
      if (inv1 > winv1.second) {
        winv1.first = Fa.shortSymbol();
        winv1.second = inv1;
      }
      if (inv2 > winv2.second) {
        winv2.first = Fa.shortSymbol();
        winv2.second = inv2;
      }
    }

    pass &= qip::check_value(&obuff, "<r>    " + worst1.first, worst1.second,
                             0.0, 1.0e-10);
    pass &= qip::check_value(&obuff, "<r^2>  " + worst2.first, worst2.second,
                             0.0, 1.0e-10);
    pass &= qip::check_value(&obuff, "<r^-1> " + winv1.first, winv1.second, 0.0,
                             1.0e-10);
    pass &= qip::check_value(&obuff, "<r^-2> " + winv2.first, winv2.second, 0.0,
                             1.0e-10);
  }

  return pass;
}

} // namespace UnitTest
