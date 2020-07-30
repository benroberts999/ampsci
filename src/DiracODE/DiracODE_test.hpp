#pragma once
#include "DiracODE/DiracODE.hpp"
#include "Maths/Grid.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Check.hpp"
// #include "qip/Maths.hpp"
// #include "qip/Vector.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
//
#include "DiracOperator/Operators.hpp"

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
    DiracODE::boundState(Fnk, en_guess, v_nuc, {}, PhysConst::alpha, 16);
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

  // auto rhat = DiracOperator::RadialF(*grid, 1);
  // auto x = rhat.radialIntegral(orbitals.front(), orbitals.front());
  // std::cout << x << "\n";

  // auto wostd::accumulate(cbegin(orbitals), cend(orbitals), 0.0, eps_en);

  // double diracen(double z, double n, int k, double alpha =
  // 0.00729735256635);

  return pass;
}

} // namespace UnitTest
