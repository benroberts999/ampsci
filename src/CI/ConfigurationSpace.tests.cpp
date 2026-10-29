#include "ConfigurationSpace.hpp"
#include "Angular/Wigner369j.hpp"
#include "CSF.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <array>
#include <memory>
#include <set>
#include <string>
#include <variant>
#include <vector>

namespace UnitTest {
// Orbitals (no radial functions) of each l, from n_min to n_max: one entry
// {n_min, n_max, l} per l; both j for l > 0. Energies are hydrogen-like with
// a quantum defect that falls with l, so low l lies lower
std::vector<DiracSpinor>
config_space_basis(const std::vector<std::array<int, 3>> &n_ranges);

// nk_index of each orbital {n, kappa}
std::vector<DiracSpinor::Index>
config_space_orbitals(const std::vector<std::pair<int, int>> &n_kappas);
} // namespace UnitTest

//==============================================================================
TEST_CASE("ConfigSpace: parse configurations", "[CI][unit]") {
  const auto s2p = CI::configuration_string(CI::parse_configuration("6s2.6p"));
  REQUIRE(s2p == "6s2.6p1");

  // A shell given twice has its electrons added
  const auto s_twice =
    CI::configuration_string(CI::parse_configuration("6s.6s"));
  REQUIRE(s_twice == "6s2");

  const auto f14 = CI::configuration_string(CI::parse_configuration("5f14"));
  REQUIRE(f14 == "5f14");

  // Invalid: no separator, too many electrons, n <= l, empty
  for (const auto &invalid : {"6s26p", "6s3", "2d1", ""}) {
    const auto configuration = CI::parse_configuration(invalid);
    REQUIRE(configuration.empty());
  }
}

//==============================================================================
// Stepwise: step 1 to 9spd, step 2 to 6spd. A double excitation needs one
// electron within 6spd; the other may go to 9spd
TEST_CASE("ConfigSpace: stepwise excitations", "[CI][unit]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 9, 0}, {4, 9, 1}, {3, 9, 2}});
  const CI::ConfigurationSpace space{basis, 2, {"4s2"}, {"1:9spd", "2:6spd"}};
  REQUIRE(!space.unrestricted());

  using UnitTest::config_space_orbitals;
  const auto leading =
    space.included(config_space_orbitals({{4, -1}, {4, -1}}));
  const auto single = space.included(config_space_orbitals({{4, -1}, {9, -1}}));
  const auto high_low =
    space.included(config_space_orbitals({{6, -1}, {9, -1}}));
  const auto high_high =
    space.included(config_space_orbitals({{9, -1}, {9, -1}}));
  const auto mid_mid =
    space.included(config_space_orbitals({{7, -1}, {8, -1}}));
  // 5p- and 6d+: both j of a shell count as that shell
  const auto low_low = space.included(config_space_orbitals({{5, 1}, {6, -3}}));

  REQUIRE(leading);
  REQUIRE(single);
  REQUIRE(high_low);
  REQUIRE(!high_high);
  REQUIRE(!mid_mid);
  REQUIRE(low_low);
}

//==============================================================================
// Two full steps include every two-electron configuration; one restricted
// step does not
TEST_CASE("ConfigSpace: two electrons, full CI", "[CI][unit]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 8, 0}, {4, 8, 1}, {3, 8, 2}});
  const CI::ConfigurationSpace default_space{basis, 2, {}, {"1", "2"}};
  const CI::ConfigurationSpace given_leading{basis, 2, {"5s.5p"}, {"1", "2"}};
  const CI::ConfigurationSpace restricted{basis, 2, {}, {"1", "2:6spd"}};

  REQUIRE(default_space.unrestricted());
  REQUIRE(given_leading.unrestricted());
  REQUIRE(!restricted.unrestricted());
}

//==============================================================================
// Lowest of each parity, by the sum of shell energies, with the j energies of
// a shell averaged with weights 2j+1
TEST_CASE("ConfigSpace: default leading", "[CI][unit]") {
  auto basis = UnitTest::config_space_basis({{4, 5, 0}, {4, 4, 1}, {3, 3, 2}});
  // 4s, 5s, 4p-, 4p+, 3d-, 3d+. 4p, weighted by 2j+1: (2*-0.25 + 4*-0.35)/6
  // = -0.317 (unweighted, -0.30)
  const std::array<double, 6> energies{-0.5, -0.32, -0.25, -0.35, -0.3, -0.3};
  for (std::size_t i = 0; i < basis.size(); ++i) {
    basis[i].en() = energies[i];
  }
  // No excitations: exactly the leading configurations
  const CI::ConfigurationSpace space{basis, 2, {}, {"0"}};

  // Even, lowest five: 4s2 (-1.0), 4s.5s (-0.82), 3d.4s (-0.8), 5s2 (-0.64),
  // 4p2 (-0.633); not 3d.5s (-0.62) or 3d2 (-0.6). Unweighted, 4p2 (-0.60)
  // would be left out instead of 3d.5s. Odd: all three
  using UnitTest::config_space_orbitals;
  const auto s2_5 = space.included(config_space_orbitals({{5, -1}, {5, -1}}));
  const auto p2 = space.included(config_space_orbitals({{4, 1}, {4, -2}}));
  const auto ds_5 = space.included(config_space_orbitals({{3, 2}, {5, -1}}));
  const auto d2 = space.included(config_space_orbitals({{3, 2}, {3, -3}}));
  const auto dp = space.included(config_space_orbitals({{3, -3}, {4, -2}}));

  REQUIRE(s2_5);
  REQUIRE(p2);
  REQUIRE(!ds_5);
  REQUIRE(!d2);
  REQUIRE(dp);
}

//==============================================================================
// Singles and doubles over the full basis are exactly the configurations with
// at least one electron in a shell of the leading configurations (whichever
// they are)
TEST_CASE("ConfigSpace: three electrons, singles and doubles",
          "[CI][unit][3VCI]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 9, 0}, {4, 9, 1}, {3, 9, 2}});
  const CI::ConfigurationSpace space{
    basis, 3, {"4s2.4p", "3d2.5s", "4s.4p.5d"}, {"1", "2"}};
  REQUIRE(!space.unrestricted());
  const std::set<std::pair<int, int>> leading_shells{
    {4, 0}, {4, 1}, {3, 2}, {5, 0}, {5, 2}};

  int num_wrong = 0;
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = i; j < basis.size(); ++j) {
      for (std::size_t k = j; k < basis.size(); ++k) {
        bool has_leading_shell = false;
        for (const auto *Fn : {&basis[i], &basis[j], &basis[k]}) {
          if (leading_shells.count({Fn->n(), Fn->l()}) != 0) {
            has_leading_shell = true;
          }
        }
        const auto included = space.included(
          {basis[i].nk_index(), basis[j].nk_index(), basis[k].nk_index()});
        if (included != has_leading_shell) {
          ++num_wrong;
        }
      }
    }
  }
  REQUIRE(num_wrong == 0);
}

//==============================================================================
TEST_CASE("ConfigSpace: adding a step only adds configurations",
          "[CI][unit][3VCI]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 9, 0}, {4, 9, 1}, {3, 9, 2}});
  const CI::ConfigurationSpace doubles{basis, 3, {"4s2.4p"}, {"1", "2"}};
  const CI::ConfigurationSpace triples{
    basis, 3, {"4s2.4p"}, {"1", "2", "3:5spd"}};
  // Steps not given are full: '3:5spd' is '1, 2, 3:5spd'
  const CI::ConfigurationSpace triples_short{basis, 3, {"4s2.4p"}, {"3:5spd"}};

  int num_lost = 0;
  int num_added = 0;
  int num_short_different = 0;
  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = i; j < basis.size(); ++j) {
      for (std::size_t k = j; k < basis.size(); ++k) {
        const std::vector<DiracSpinor::Index> orbitals{
          basis[i].nk_index(), basis[j].nk_index(), basis[k].nk_index()};
        const auto in_doubles = doubles.included(orbitals);
        const auto in_triples = triples.included(orbitals);
        num_lost += in_doubles && !in_triples ? 1 : 0;
        num_added += !in_doubles && in_triples ? 1 : 0;
        num_short_different +=
          triples_short.included(orbitals) != in_triples ? 1 : 0;
      }
    }
  }
  REQUIRE(num_lost == 0);
  REQUIRE(num_added > 0);
  REQUIRE(num_short_different == 0);
}

//==============================================================================
// Only CSFs of included configurations are formed; all of them
TEST_CASE("ConfigSpace: three-electron CSFs", "[CI][unit][3VCI]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 6, 0}, {4, 6, 1}, {3, 6, 2}});
  const CI::ConfigurationSpace space{basis, 3, {"4s2.4p"}, {"1", "2:5spd"}};

  for (const auto parity : {+1, -1}) {
    const auto all_CSFs = CI::form_CSF3s(3, parity, basis);
    std::size_t num_expected = 0;
    for (const auto &csf : all_CSFs) {
      if (space.included(csf.orbitals())) {
        ++num_expected;
      }
    }
    const auto num_CSFs = CI::form_CSF3s(3, parity, basis, space).size();
    REQUIRE(num_CSFs == num_expected);
    REQUIRE(num_CSFs < all_CSFs.size());
  }
}

//==============================================================================
TEST_CASE("ConfigSpace: parse hole configurations", "[CI][unit][ciholes]") {
  const auto hole =
    CI::configuration_string(CI::parse_configuration("6p2.6s-1"));
  REQUIRE(hole == "6s-1.6p2");

  // Invalid: no number after '-', a non-integer number, an empty shell
  for (const auto &invalid : {"6s-.6p2", "6s-x.6p2", "6s0.6p"}) {
    const auto configuration = CI::parse_configuration(invalid);
    REQUIRE(configuration.empty());
  }
}

//==============================================================================
// One electron, with hole CSFs (hole in 3s; the two particles within 6spd).
// The hole is created by an excitation from a configuration without one, so
// 3s-1.4p2 and 3s-1.4p.5s are single excitations from 4p; 3s-1.5s.5p is a
// double. From a leading configuration with the hole, the particles are
// excited; the hole is never filled
TEST_CASE("ConfigSpace: hole excitations", "[CI][unit][ciholes]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 9, 0}, {4, 9, 1}, {3, 9, 2}});
  const auto holes = UnitTest::config_space_basis({{3, 3, 0}});
  const auto hole_basis =
    UnitTest::config_space_basis({{4, 6, 0}, {4, 6, 1}, {3, 6, 2}});
  using UnitTest::config_space_orbitals;
  const auto hole_3s = config_space_orbitals({{3, -1}});
  const auto p4 = config_space_orbitals({{4, 1}});
  const auto p4_p4 = config_space_orbitals({{4, 1}, {4, -2}});
  const auto p4_s5 = config_space_orbitals({{4, 1}, {5, -1}});
  const auto s5_p5 = config_space_orbitals({{5, -1}, {5, 1}});

  const CI::ConfigurationSpace singles{basis, 1,     {"4p"},
                                       {"1"}, holes, hole_basis};
  REQUIRE(!singles.unrestricted());
  const auto single_p2 = singles.included(p4_p4, hole_3s);
  const auto single_ps = singles.included(p4_s5, hole_3s);
  const auto double_sp = singles.included(s5_p5, hole_3s);
  REQUIRE(single_p2);
  REQUIRE(single_ps);
  REQUIRE(!double_sp);

  // Two particles and a hole: two full steps include every configuration
  const CI::ConfigurationSpace doubles{basis,      1,     {"4p"},
                                       {"1", "2"}, holes, hole_basis};
  REQUIRE(doubles.unrestricted());

  const CI::ConfigurationSpace from_hole{basis, 1,     {"3s-1.4p2"},
                                         {"1"}, holes, hole_basis};
  REQUIRE(!from_hole.unrestricted());
  const auto hole_kept = from_hole.included(p4_s5, hole_3s);
  const auto hole_filled = from_hole.included(p4);
  const auto hole_double = from_hole.included(s5_p5, hole_3s);
  REQUIRE(hole_kept);
  REQUIRE(!hole_filled);
  REQUIRE(!hole_double);

  // The default leading configurations have no holes
  const CI::ConfigurationSpace leading_only{basis, 1,     {},
                                            {"0"}, holes, hole_basis};
  const auto default_electron = leading_only.included(p4);
  const auto default_hole = leading_only.included(p4_p4, hole_3s);
  REQUIRE(default_electron);
  REQUIRE(!default_hole);
}

//==============================================================================
// Only hole CSFs of included configurations are formed; all of them
TEST_CASE("ConfigSpace: hole CSFs", "[CI][unit][ciholes]") {
  const auto basis =
    UnitTest::config_space_basis({{4, 6, 0}, {4, 6, 1}, {3, 6, 2}});
  const auto holes = UnitTest::config_space_basis({{3, 3, 0}});
  const CI::ConfigurationSpace space{basis, 1, {"4p"}, {"1"}, holes};

  for (const auto parity : {+1, -1}) {
    const auto all_CSFs = CI::form_CSF2p1hs(3, parity, basis, holes);
    std::size_t num_expected = 0;
    for (const auto &csf : all_CSFs) {
      const auto &hole_csf = std::get<CI::CSF2p1h>(csf.variant());
      if (space.included({hole_csf.state(0), hole_csf.state(1)},
                         {hole_csf.hole})) {
        ++num_expected;
      }
    }
    const auto num_CSFs =
      CI::form_CSF2p1hs(3, parity, basis, holes, space).size();
    REQUIRE(num_CSFs == num_expected);
    REQUIRE(num_CSFs > 0);
    REQUIRE(num_CSFs < all_CSFs.size());
  }
}

//==============================================================================
std::vector<DiracSpinor>
UnitTest::config_space_basis(const std::vector<std::array<int, 3>> &n_ranges) {
  const auto grid = std::make_shared<const Grid>(
    GridParameters{10, 1.0e-4, 10.0, 5.0, GridType::loglinear});
  std::vector<DiracSpinor> basis;
  for (const auto &[n_min, n_max, l] : n_ranges) {
    const auto quantum_defect = std::max(0.0, 1.0 - 0.4 * l);
    for (int n = n_min; n <= n_max; ++n) {
      const auto n_effective = n - quantum_defect;
      const auto energy = -0.5 / (n_effective * n_effective);
      if (l > 0) {
        basis.emplace_back(n, l, grid);
        basis.back().en() = energy;
      }
      basis.emplace_back(n, -l - 1, grid);
      basis.back().en() = energy;
    }
  }
  return basis;
}

//==============================================================================
std::vector<DiracSpinor::Index> UnitTest::config_space_orbitals(
  const std::vector<std::pair<int, int>> &n_kappas) {
  std::vector<DiracSpinor::Index> orbitals;
  for (const auto &[n, kappa] : n_kappas) {
    orbitals.push_back(DiracSpinor::Index(Angular::nk_to_index(n, kappa)));
  }
  return orbitals;
}
