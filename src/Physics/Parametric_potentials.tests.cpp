#include "DiracOperator/DiracOperator.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"

//==============================================================================
//==============================================================================
//! Unit tests for Ginges/Flambaum Radiative potential method
TEST_CASE("Physics: Parametric potentials", "[prm][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Physics: Parametric potentials\n";

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

  // geneated with 5000 pts
  const auto aga0s =
      std::vector{35.851733939957235, 23.857128488699516, 23.806698412657056,
                  15.520681948768228, 15.508963307396499, 10.439854223306053,
                  10.436137600652057};
  const auto ata0s =
      std::vector{30.913482805598647, 23.504270400835946, 23.468519915732003,
                  16.557099169686712, 16.546379062750177, 11.415329473402792,
                  11.411383617964162};
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
