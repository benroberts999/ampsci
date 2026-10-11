#include "Wavefunction/DiracSpinor.hpp"
#include "DiracODE/include.hpp"
#include "Maths/Grid.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <vector>

//==============================================================================
TEST_CASE("DiracSpinor: num_nodes", "[DiracSpinor][unit]") {

  const auto grid = std::make_shared<const Grid>(1.0e-6, 100.0, 2000,
                                                 GridType::loglinear, 10.0);
  const double Zeff = 10.0;
  const auto v_nuc = Nuclear::sphericalNuclearPotential(Zeff, 0.0, grid->r());

  for (int n = 1; n <= 7; ++n) {
    for (int l = 0; l < n; ++l) {
      for (const int kappa : {l, -l - 1}) {
        if (kappa == 0) {
          continue;
        }
        DiracSpinor Fnk(n, kappa, grid);
        DiracODE::boundState(Fnk, -0.5 * Zeff * Zeff / (n * n), v_nuc);
        const int expected = n - l - 1;
        REQUIRE(Fnk.num_nodes() == expected);

        // Small-amplitude sign changes in the tail are not physical nodes
        const auto f_max = std::abs(*std::max_element(
          Fnk.f().begin(), Fnk.f().end(),
          [](double a, double b) { return std::abs(a) < std::abs(b); }));
        for (auto i = Fnk.max_pt() - 20; i < Fnk.max_pt(); ++i) {
          Fnk.f(i) = (i % 2 == 0 ? 1.0 : -1.0) * 1.0e-5 * f_max;
        }
        REQUIRE(Fnk.num_nodes() == expected);
        // Without the amplitude cut, they would be counted
        REQUIRE(Fnk.num_nodes(0.0) > expected);

        // Inner lobes are real nodes, even when small (e.g., core-penetrating
        // Rydberg states): shrink everything inside the first node
        if (expected > 0) {
          auto first_node = Fnk.min_pt() + 1;
          while (Fnk.f(first_node - 1) * Fnk.f(first_node) >= 0.0) {
            ++first_node;
          }
          for (auto i = Fnk.min_pt(); i < first_node; ++i) {
            Fnk.f(i) *= 1.0e-6;
          }
          REQUIRE(Fnk.num_nodes() == expected);
        }
      }
    }
  }

  // Failed solutions (zero or non-finite orbital) are reported as -1
  DiracSpinor Fzero(2, -1, grid);
  REQUIRE(Fzero.num_nodes() == -1);
  DiracODE::boundState(Fzero, -0.5 * Zeff * Zeff / 4.0, v_nuc);
  Fzero.f(Fzero.max_pt() / 2) = std::numeric_limits<double>::quiet_NaN();
  REQUIRE(Fzero.num_nodes() == -1);
}
