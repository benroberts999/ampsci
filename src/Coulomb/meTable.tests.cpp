#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp"
#include "catch2/catch.hpp"
#include <random>

//==============================================================================

//==============================================================================
TEST_CASE("Coulomb: meTable", "[Coulomb][meTable][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: meTable\n";

  const auto radial_grid = std::make_shared<const Grid>(
      GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 6;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
    }
    orbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
  }

  const auto h = DiracOperator::E1(*radial_grid);

  Coulomb::meTable<double> met;

  // store all ME's in table - directly calculate each (no symmetry check)
  for (const auto &a : orbs) {
    for (const auto &b : orbs) {
      if (h.isZero(a, b))
        continue;
      const auto d_ab = h.reducedME(a, b);
      met.add(a, b, d_ab);
    }
  }

  // Now, loop through and check, only for a>=b:
  for (const auto &a : orbs) {
    for (const auto &b : orbs) {
      if (b > a)
        continue;
      const auto dab = h.reducedME(a, b);
      const auto dba = h.symm_sign(a, b) * dab;
      if (h.isZero(a, b)) {
        REQUIRE(!met.contains(a, b));
        REQUIRE(!met.contains(b, a));
      } else {
        REQUIRE(std::abs(dab - *met.get(a, b)) < 1.0e-12);
        REQUIRE(std::abs(dba - *met.get(b, a)) < 1.0e-12);

        REQUIRE(*met.get(a, b) ==
                Approx(*met.get(a.shortSymbol(), b.shortSymbol())));
        REQUIRE(*met.get(b, a) ==
                Approx(*met.get(b.shortSymbol(), a.shortSymbol())));
        REQUIRE(*met.get(a, b) == Approx(*met.get(a.symbol(), b.symbol())));
        REQUIRE(*met.get(b, a) == Approx(*met.get(b.symbol(), a.symbol())));
      }
    }
  }

  // Test converting key to symbols
  for (auto &[key, val] : met) {
    const auto [a, b] = Coulomb::meTable<double>::index_to_symbols(key);
    const auto [a2, b2] = met.index_to_symbols(key);
    REQUIRE(a == a2);
    REQUIRE(b == b2);
    const auto retrieved = met.get(a, b);
    REQUIRE(retrieved != nullptr);
    REQUIRE(*retrieved == val);
  }
  // .. and for const case!
  for (const auto &[key, val] : met) {
    const auto [a, b] = Coulomb::meTable<double>::index_to_symbols(key);
    const auto [a2, b2] = met.index_to_symbols(key);
    REQUIRE(a == a2);
    REQUIRE(b == b2);
    const auto retrieved = met.get(a, b);
    REQUIRE(retrieved != nullptr);
    REQUIRE(*retrieved == val);
  }

  // Test insert:
  // Split table into two
  Coulomb::meTable<double> meta, metb;
  int count = 0;
  for (const auto &key_val : met) {
    ++count;
    if (count % 2 == 0) {
      meta->insert(key_val);
    } else {
      metb->insert(key_val);
    }
  }
  // ensure none of the "b" elements are in the "a" map
  for (const auto &[key, val] : metb) {
    const auto [a, b] = met.index_to_symbols(key);
    REQUIRE(meta.get(a, b) == nullptr);
  }
  // insert b back into a
  meta.add(metb);
  // Now, ensure all of the "b" elements are in the "a" map
  for (const auto &[key, val] : metb) {
    const auto [a, b] = met.index_to_symbols(key);
    REQUIRE(meta.get(a, b) != nullptr);
  }

  // test update
  for (const auto &a : orbs) {
    for (const auto &b : orbs) {
      if (h.isZero(a, b))
        continue;
      met.update(a, b, 2.0);
      REQUIRE(std::abs(*met.get(a, b) - 2.0) < 1.0e-12);
    }
  }
}

//============================================================================
TEST_CASE("Coulomb: meTable<vector>", "[Coulomb][meTable][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Coulomb: meTable<vector>\n";

  // build set of blank orbitals, one n for each kappa up to l=lmax
  const int lmax = 6;
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.emplace_back(n_min, l, nullptr);
    }
    orbs.emplace_back(n_min, -l - 1, nullptr);
  }

  Coulomb::meTable<std::vector<int>> met_v;

  int i = 0;
  for (const auto &a : orbs) {
    for (const auto &b : orbs) {
      ++i;
      auto v = std::vector<int>{i, i + 1, i + 2};
      const auto p1 = v.data();
      met_v.add(a, b, std::move(v));
      const auto gotten = met_v.get(a, b);
      const auto p2 = gotten->data();
      // test move works correctly
      REQUIRE(p1 == p2);
      REQUIRE(gotten->at(0) == i);
      REQUIRE(gotten->at(1) == i + 1);
      REQUIRE(gotten->at(2) == i + 2);
      REQUIRE(gotten->size() == 3);
    }
  }
}