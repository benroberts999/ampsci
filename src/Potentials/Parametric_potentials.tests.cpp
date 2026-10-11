#include "Parametric_potentials.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <algorithm>

//==============================================================================
//==============================================================================
//! Unit tests for Parametric potentials
TEST_CASE("Parametric potentials", "[prm][unit]") {

  // ensure no invalid numbers
  for (int z = 2; z < 137; ++z) {
    double H{-1.0}, d{-1.0}, g{-1.0}, t{-1.0};
    Parametric::defaultGreenCore(z, H, d);
    REQUIRE((H > 0.05 && H < 100.0 && d > 0.01 && d < 10.0));
    Parametric::defaultGreen(z, H, d);
    REQUIRE((H > 0.05 && H < 100.0 && d > 0.01 && d < 10.0));
    Parametric::defaultTietz(z, g, t);
    REQUIRE((g > 0.05 && g < 100.0 && t > 0.01 && t < 10.0));
  }

  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-6, 75.0, 50.0, GridType::loglinear});
  const double zeff = 15.0;
  const int lmax = 3;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
    }
    orbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
  }

  const auto vg = Parametric::GreenPotential(int(zeff), radial_grid->r());
  const auto vt = Parametric::TietzPotential(int(zeff), radial_grid->r());

  // Z = 1: no screening electrons, so electron potential is exactly zero
  const auto v_H = Parametric::GreenPotential(1, radial_grid->r());
  const auto result =
    std::all_of(v_H.cbegin(), v_H.cend(), [](double v) { return v == 0.0; });
  REQUIRE(result);

  // geneated with 5000 pts
  const auto aga0s =
    std::vector{35.851733939957235, 23.857128488699516, 23.806698412657056,
                15.520681948768228, 15.508963307396499, 10.439854223306053,
                10.436137600652057};
  // updated: were incorrect
  const auto ata0s =
    std::vector{33.386808487089361, 24.565794637034891, 24.524491657740008,
                16.856261639250768, 16.844655540371761, 11.455226123489282,
                11.451146314342445};
  std::size_t count = 0;
  for (const auto &Fa : orbs) {
    auto aga = Fa * (vg * Fa);
    const auto aga0 = aga0s.at(count);
    REQUIRE(std::abs(aga - aga0) < 1.0e-6);
    auto ata = Fa * (vt * Fa);
    const auto ata0 = ata0s.at(count);
    REQUIRE(std::abs(ata - ata0) < 1.0e-6);
    ++count;
  }
}
