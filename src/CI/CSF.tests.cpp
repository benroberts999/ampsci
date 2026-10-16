#include "CSF.hpp"
#include "Angular/include.hpp"
#include "CI_Integrals.hpp"
#include "Maths/Grid.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

//==============================================================================
// Tests the text form of a CI level reference: J{+,-}:index, and e/oJ:index
TEST_CASE("CI: level references", "[CI][unit]") {

  std::cout << "CI level references\n";

  // {text, 2J, parity, index}
  const std::vector<std::tuple<std::string, int, int, std::size_t>> cases{
    {"2+:3", 4, 1, 3},  {"e2:3", 4, 1, 3},    {"2-:3", 4, -1, 3},
    {"o2:3", 4, -1, 3}, {"0+", 0, 1, 0},      {"e0", 0, 1, 0},
    {"1-", 2, -1, 0},   {"o1", 2, -1, 0},     {"10+:12", 20, 1, 12},
    {"1-:2", 2, -1, 2}, {"3/2-:1", 3, -1, 1}, {"o1/2", 1, -1, 0},
    {"5/2+", 5, 1, 0},  {"1.5-:2", 3, -1, 2}};

  for (const auto &[text, twoJ, parity, index] : cases) {
    const auto level = CI::parse_level(text);
    REQUIRE(level.has_value());
    REQUIRE(level->twoJ == twoJ);
    REQUIRE(level->parity == parity);
    REQUIRE(level->index == index);
  }

  // Round trip: the standard form must parse back to itself
  for (const auto &text : {"0+:0", "2-:3", "10+:12", "3/2-:1"}) {
    const auto level = CI::parse_level(text);
    REQUIRE(level.has_value());
    REQUIRE(CI::to_string(*level) == text);
  }

  // Malformed
  for (const auto &text : {"", " ", "2", "+2", "2+:", "2+:x", "x2:3", "2*:3",
                           "-2+", "2+3", "1/3-:0", "/2-:0", "1.25+"}) {
    fmt::print("'{}' -> {}\n", text,
               CI::parse_level(text) ? "parsed" : "rejected");
    REQUIRE(!CI::parse_level(text).has_value());
  }
}

//==============================================================================
// Tests the jj -> LS recoupling amplitudes and <L^2>, <S^2> expectation values
TEST_CASE("CI: LS coupling", "[CI][unit]") {

  std::cout << "CI LS coupling\n";

  // Sum over (L, S) of A^2, for CSF with given quantum numbers
  const auto sum_A2 = [](int n1, int k1, int n2, int k2, int twoJ) {
    const auto l1 = Angular::l_k(k1);
    const auto l2 = Angular::l_k(k2);
    const auto twoj1 = Angular::twoj_k(k1);
    const auto twoj2 = Angular::twoj_k(k2);
    double sum = 0.0;
    for (int S = 0; S <= 1; ++S) {
      for (int L = std::abs(l1 - l2); L <= l1 + l2; ++L) {
        const auto A =
          CI::LS_amplitude(n1, l1, twoj1, n2, l2, twoj2, L, S, twoJ);
        sum += A * A;
      }
    }
    return sum;
  };

  // Unitarity: each CSF carries unit total LS weight.
  // Cases cover: distinct configs, same shell same j, same shell different j,
  // same l different n. {n1, kappa1, n2, kappa2, twoJ}
  const std::vector<std::array<int, 5>> csf_cases{
    {4, -1, 3, 2, 4},  // 4s1/2 3d3/2, J=2
    {4, -1, 3, -3, 4}, // 4s1/2 3d5/2, J=2
    {4, -1, 3, -3, 6}, // 4s1/2 3d5/2, J=3
    {3, 1, 3, -2, 2},  // 3p1/2 3p3/2, J=1
    {3, 1, 3, -2, 4},  // 3p1/2 3p3/2, J=2
    {3, 1, 3, 1, 0},   // (3p1/2)^2,   J=0
    {3, -2, 3, -2, 0}, // (3p3/2)^2,   J=0
    {3, -2, 3, -2, 4}, // (3p3/2)^2,   J=2
    {3, 1, 4, -2, 2},  // 3p1/2 4p3/2, J=1 (different n: distinct config)
    {3, 2, 3, -3, 2},  // 3d3/2 3d5/2, J=1
    {3, 2, 3, -3, 6},  // 3d3/2 3d5/2, J=3
    {3, -3, 3, -3, 8}, // (3d5/2)^2,   J=4
    {3, 2, 4, 1, 2} // 3d3/2 4p1/2, J=1
  };
  for (const auto &[n1, k1, n2, k2, twoJ] : csf_cases) {
    REQUIRE(sum_A2(n1, k1, n2, k2, twoJ) == Approx(1.0).epsilon(1.0e-12));
  }

  // Known values: p^2 (same n), J=0.
  // jj basis: {(p1/2)^2, (p3/2)^2}; LS terms: 1S0, 3P0
  // Weights: (p1/2)^2 = {1/3, 2/3}; (p3/2)^2 = {2/3, 1/3}
  const auto a11 = CI::LS_amplitude(3, 1, 1, 3, 1, 1, 0, 0, 0);
  const auto a12 = CI::LS_amplitude(3, 1, 1, 3, 1, 1, 1, 1, 0);
  const auto a21 = CI::LS_amplitude(3, 1, 3, 3, 1, 3, 0, 0, 0);
  const auto a22 = CI::LS_amplitude(3, 1, 3, 3, 1, 3, 1, 1, 0);
  REQUIRE(a11 * a11 == Approx(1.0 / 3).epsilon(1.0e-12));
  REQUIRE(a12 * a12 == Approx(2.0 / 3).epsilon(1.0e-12));
  REQUIRE(a21 * a21 == Approx(2.0 / 3).epsilon(1.0e-12));
  REQUIRE(a22 * a22 == Approx(1.0 / 3).epsilon(1.0e-12));
  // rows of the transformation matrix are orthogonal:
  REQUIRE(a11 * a21 + a12 * a22 == Approx(0.0).margin(1.0e-12));

  // 3p1/2 3p3/2, J=1: only 3P1 exists (1P1 is Pauli-forbidden: L+S odd)
  REQUIRE(CI::LS_amplitude(3, 1, 1, 3, 1, 3, 1, 0, 2) == 0.0);
  const auto a3P1 = CI::LS_amplitude(3, 1, 1, 3, 1, 3, 1, 1, 2);
  REQUIRE(a3P1 * a3P1 == Approx(1.0).epsilon(1.0e-12));

  //----------------------------------------------------------------------------
  // <L^2>, <S^2> for hand-built states:

  const auto grid = std::make_shared<const Grid>(
    GridParameters{10, 1.0e-4, 10.0, 5.0, GridType::loglinear});
  const DiracSpinor s4{4, -1, grid};
  const DiracSpinor p1{3, 1, grid};
  const DiracSpinor p3{3, -2, grid};
  const DiracSpinor d3{3, 2, grid};
  const DiracSpinor d5{3, -3, grid};

  const auto view = [](const std::vector<double> &c) {
    return LinAlg::View<const double>{c.data(), 0, c.size(), 1};
  };

  {
    // Single CSF {4s, 3d5/2}, J=3: pure 3D3
    const std::vector<CI::CSF2> csfs{{s4, d5}};
    const std::vector<double> c{1.0};
    const auto [L2, S2] = CI::expectation_L2S2(view(c), csfs, 6);
    REQUIRE(L2 == Approx(6.0).epsilon(1.0e-12));
    REQUIRE(S2 == Approx(2.0).epsilon(1.0e-12));
  }

  {
    // Single CSF {3p1/2, 3p3/2}, J=1: pure 3P1 (tests the sqrt(2) factor)
    const std::vector<CI::CSF2> csfs{{p1, p3}};
    const std::vector<double> c{1.0};
    const auto [L2, S2] = CI::expectation_L2S2(view(c), csfs, 2);
    REQUIRE(L2 == Approx(2.0).epsilon(1.0e-12));
    REQUIRE(S2 == Approx(2.0).epsilon(1.0e-12));
  }

  {
    // {4s, 3d3/2} and {4s, 3d5/2}, J=2: mixtures of 1D2 and 3D2.
    // L2 = 6 always; total 3D2 weight across the two CSFs is 1
    const std::vector<CI::CSF2> csfs{{s4, d3}, {s4, d5}};
    const std::vector<double> ca{1.0, 0.0};
    const std::vector<double> cb{0.0, 1.0};
    const auto [L2a, S2a] = CI::expectation_L2S2(view(ca), csfs, 4);
    const auto [L2b, S2b] = CI::expectation_L2S2(view(cb), csfs, 4);
    REQUIRE(L2a == Approx(6.0).epsilon(1.0e-12));
    REQUIRE(L2b == Approx(6.0).epsilon(1.0e-12));
    REQUIRE(S2a + S2b == Approx(2.0).epsilon(1.0e-12));
  }

  {
    // (3p)^2, J=0: CI coefficients equal to the LS amplitudes of 3P0
    // give the pure 3P0 state; those of 1S0 give pure 1S0
    const std::vector<CI::CSF2> csfs{{p1, p1}, {p3, p3}};
    const std::vector<double> c_3P0{a12, a22};
    const std::vector<double> c_1S0{a11, a21};
    const auto [L2a, S2a] = CI::expectation_L2S2(view(c_3P0), csfs, 0);
    const auto [L2b, S2b] = CI::expectation_L2S2(view(c_1S0), csfs, 0);
    REQUIRE(L2a == Approx(2.0).epsilon(1.0e-12));
    REQUIRE(S2a == Approx(2.0).epsilon(1.0e-12));
    REQUIRE(L2b == Approx(0.0).margin(1.0e-12));
    REQUIRE(S2b == Approx(0.0).margin(1.0e-12));
  }
}

//==============================================================================
// Tests the three-electron CSFs: the number formed for each configuration and
// J (against an independent count of m-scheme determinants), and that each is
// normalised, orthogonal, and fully antisymmetric (by explicit expansion over
// single-particle m)
TEST_CASE("CI: three-electron CSFs", "[CI][3VCI][unit]") {

  std::cout << "CI three-electron CSFs\n";

  const auto grid = std::make_shared<const Grid>(
    GridParameters{10, 1.0e-4, 10.0, 5.0, GridType::loglinear});

  // j = 1/2 to 9/2: includes j^3 with unique J, e.g. (5/2)^3, and with
  // repeated J: (9/2)^3 has two J = 9/2 states
  std::vector<DiracSpinor> basis;
  for (const auto &[n, kappa] : std::vector<std::pair<int, int>>{
         {4, -1}, {5, -1}, {4, 1}, {4, -2}, {3, 2}, {4, 3}, {5, -5}}) {
    basis.emplace_back(n, kappa, grid);
  }

  // Number of determinants of the configuration with total 2*M = twoM.
  // Equal orbitals are adjacent (sorted); their m strictly increase (Pauli)
  const auto num_determinants =
    [](const std::array<DiracSpinor::Index, 3> &orbitals, int twoM) {
      const auto tj0 = Angular::nkindex_to_twoj(orbitals[0]);
      const auto tj1 = Angular::nkindex_to_twoj(orbitals[1]);
      const auto tj2 = Angular::nkindex_to_twoj(orbitals[2]);
      int count = 0;
      for (int m0 = -tj0; m0 <= tj0; m0 += 2) {
        for (int m1 = -tj1; m1 <= tj1; m1 += 2) {
          if (orbitals[1] == orbitals[0] && m1 <= m0)
            continue;
          const auto m2 = twoM - m0 - m1;
          if (std::abs(m2) > tj2)
            continue;
          if (orbitals[2] == orbitals[1] && m2 <= m1)
            continue;
          ++count;
        }
      }
      return count;
    };

  // Total number of CSFs for each {2J, parity}, summed over configurations
  std::map<std::pair<int, int>, std::size_t> expected_totals;

  std::size_t worst_count_error = 0;
  double worst_norm = 0.0;
  double worst_antisymmetry = 0.0;
  double worst_overlap = 0.0;

  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = i; j < basis.size(); ++j) {
      for (std::size_t k = j; k < basis.size(); ++k) {
        std::array orbitals{basis[i].nk_index(), basis[j].nk_index(),
                            basis[k].nk_index()};
        std::sort(orbitals.begin(), orbitals.end());
        const auto parity =
          basis[i].parity() * basis[j].parity() * basis[k].parity();
        const auto max_twoJ =
          basis[i].twoj() + basis[j].twoj() + basis[k].twoj();

        for (int twoJ = 1; twoJ <= max_twoJ; twoJ += 2) {
          const auto csfs = CI::form_CSF3s(orbitals, twoJ);

          const auto expected =
            std::size_t(num_determinants(orbitals, twoJ) -
                        num_determinants(orbitals, twoJ + 2));
          expected_totals[{twoJ, parity}] += expected;
          worst_count_error = std::max(
            worst_count_error, csfs.size() > expected ? csfs.size() - expected :
                                                        expected - csfs.size());

          // Expansion of each CSF over single-particle states, at M = J:
          // amplitude of |o1 m1>|o2 m2>|o3 m3>, key {o1, 2m1, o2, 2m2, o3, 2m3}
          std::vector<std::map<std::array<int, 6>, double>> expansions;
          for (const auto &csf : csfs) {
            auto &amplitude = expansions.emplace_back();
            for (const auto &parent : csf.parents()) {
              const auto [x, y] = parent.pair;
              const auto z = parent.spectator;
              const auto tjx = Angular::nkindex_to_twoj(x);
              const auto tjy = Angular::nkindex_to_twoj(y);
              const auto tjz = Angular::nkindex_to_twoj(z);
              // Parent pair, antisymmetric in electrons 1 and 2 (as CSF2):
              // {orbital of 1, orbital of 2, coefficient}
              const auto eta = x == y ? 1.0 : 1.0 / std::sqrt(2.0);
              std::vector<std::tuple<int, int, double>> pair_terms{{x, y, eta}};
              if (x != y) {
                const auto s = Angular::neg1pow_2(tjx + tjy - parent.twoJ2);
                pair_terms.push_back({y, x, -s * eta});
              }
              for (const auto &[o1, o2, c] : pair_terms) {
                const auto tj1 = Angular::nkindex_to_twoj(std::uint64_t(o1));
                const auto tj2 = Angular::nkindex_to_twoj(std::uint64_t(o2));
                for (int m1 = -tj1; m1 <= tj1; m1 += 2) {
                  for (int m2 = -tj2; m2 <= tj2; m2 += 2) {
                    const auto twoM2 = m1 + m2;
                    const auto m3 = twoJ - twoM2;
                    if (std::abs(m3) > tjz)
                      continue;
                    const auto cg =
                      Angular::cg_2(tj1, m1, tj2, m2, parent.twoJ2, twoM2) *
                      Angular::cg_2(parent.twoJ2, twoM2, tjz, m3, twoJ, twoJ);
                    if (cg == 0.0)
                      continue;
                    amplitude[{o1, m1, o2, m2, int(z), m3}] +=
                      parent.coef * c * cg;
                  }
                }
              }
            }
          }

          for (std::size_t a = 0; a < expansions.size(); ++a) {
            const auto &amplitude = expansions[a];
            double norm2 = 0.0;
            for (const auto &[key, value] : amplitude) {
              norm2 += value * value;
              const auto [o1, m1, o2, m2, o3, m3] = key;
              // Exchange of electrons 1-2, 1-3, 2-3: amplitude changes sign
              for (const auto &exchanged :
                   {std::array{o2, m2, o1, m1, o3, m3},
                    std::array{o3, m3, o2, m2, o1, m1},
                    std::array{o1, m1, o3, m3, o2, m2}}) {
                const auto it = amplitude.find(exchanged);
                const auto value_exchanged =
                  it == amplitude.end() ? 0.0 : it->second;
                worst_antisymmetry = std::max(
                  worst_antisymmetry, std::abs(value + value_exchanged));
              }
            }
            worst_norm = std::max(worst_norm, std::abs(norm2 - 1.0));

            for (std::size_t b = a + 1; b < expansions.size(); ++b) {
              double overlap = 0.0;
              for (const auto &[key, value] : expansions[b]) {
                const auto it = amplitude.find(key);
                if (it != amplitude.end()) {
                  overlap += value * it->second;
                }
              }
              worst_overlap = std::max(worst_overlap, std::abs(overlap));
            }
          }
        }
      }
    }
  }

  fmt::print("N_CSF vs N(M=J) - N(M=J+1): {}\n", worst_count_error);
  fmt::print("<A|A> - 1:                  {:.1e}\n", worst_norm);
  fmt::print("P_ij|A> + |A>:              {:.1e}\n", worst_antisymmetry);
  fmt::print("<A|B>:                      {:.1e}\n", worst_overlap);
  REQUIRE(worst_count_error == 0);
  REQUIRE(worst_norm < 1.0e-12);
  REQUIRE(worst_antisymmetry < 1.0e-12);
  REQUIRE(worst_overlap < 1.0e-12);

  // All configurations of the basis together
  for (const auto &[twoJ_parity, expected] : expected_totals) {
    const auto [twoJ, parity] = twoJ_parity;
    REQUIRE(CI::form_CSF3s(twoJ, parity, basis).size() == expected);
  }

  // j^3 with repeated J: two J = 9/2 states of (5g9/2)^3, labelled by the
  // J12 of the coupling each was formed from
  const auto g = basis.back().nk_index();
  const auto g3 = CI::form_CSF3s({g, g, g}, 9);
  REQUIRE(g3.size() == 2);
  REQUIRE(g3[0].twoJ12() != g3[1].twoJ12());
  REQUIRE(g3[0].config(true) == "5g+^3");
  REQUIRE(g3[0].config(false) == "5g^3");
}

//==============================================================================
// The smaller CSF-level functions: J text forms, CSF3 accessors and
// comparisons, the electron-2-3 exchange, PsiJPi accessors, the leading
// configuration, and the two-spectator dSigma/dE correction
TEST_CASE("CI: CSF functions", "[CI][3VCI][unit]") {

  std::cout << "CI CSF functions\n";

  // J as text
  REQUIRE(CI::J_string(0) == "0");
  REQUIRE(CI::J_string(4) == "2");
  REQUIRE(CI::J_string(3) == "3/2");
  REQUIRE(CI::parse_twoJ(" 2 ") == 4);
  REQUIRE(CI::parse_twoJ("5/2") == 5);
  REQUIRE(CI::parse_twoJ("1.5") == 3);
  REQUIRE(CI::parse_twoJ("0.5") == 1);
  REQUIRE(CI::parse_twoJ("2.0") == 4);
  REQUIRE(CI::parse_twoJ("2.50") == 5);
  for (const auto &text : {"", "5/", "/2", "1/3", "1.25", "1.", ".5", "x"}) {
    REQUIRE(!CI::parse_twoJ(text).has_value());
  }

  const auto grid = std::make_shared<const Grid>(
    GridParameters{10, 1.0e-4, 10.0, 5.0, GridType::loglinear});
  const DiracSpinor s{4, -1, grid};
  const DiracSpinor p1{4, 1, grid};
  const DiracSpinor p3{4, -2, grid};
  const DiracSpinor d5{3, -3, grid};
  const auto is = s.nk_index();
  const auto ip1 = p1.nk_index();
  const auto ip3 = p3.nk_index();
  const auto id5 = d5.nk_index();

  // CSF3: sorted storage, accessors, parity, electron number
  const CI::CSF3 spp{{ip3, s.nk_index(), ip1}, 2, {}};
  REQUIRE(spp.state(0) == is);
  REQUIRE(spp.state(1) == ip1);
  REQUIRE(spp.state(2) == ip3);
  REQUIRE(spp.twoJ12() == 2);
  REQUIRE(spp.parity() == +1);
  REQUIRE(spp.num_electrons() == 3);
  REQUIRE(CI::CSF2{s, p1}.num_electrons() == 2);
  REQUIRE(CI::CSF3{{is, is, ip1}, 0, {}}.parity() == -1);

  // Equality: same orbitals and the same J12 label
  REQUIRE(spp == CI::CSF3{{is, ip1, ip3}, 2, {}});
  REQUIRE(spp != CI::CSF3{{is, ip1, ip3}, 4, {}});
  REQUIRE(spp != CI::CSF3{{is, ip1, id5}, 2, {}});

  // Number of different orbitals, counting repeats
  REQUIRE(CI::CSF3::num_different(spp, spp) == 0);
  REQUIRE(CI::CSF3::num_different(spp, {{is, ip1, id5}, 0, {}}) == 1);
  REQUIRE(CI::CSF3::num_different(spp, {{is, id5, id5}, 0, {}}) == 2);
  REQUIRE(CI::CSF3::num_different(spp, {{id5, id5, id5}, 0, {}}) == 3);
  REQUIRE(CI::CSF3::num_different({{ip3, ip3, is}, 0, {}},
                                  {{ip3, is, is}, 0, {}}) == 1);

  // Configuration strings: non-rel shells merged, powers separated
  REQUIRE(spp.config(false) == "4s4p^2");
  REQUIRE(spp.config(true) == "4s+4p-4p+");
  REQUIRE(CI::CSF3({is, is, ip3}, 0, {}).config(false) == "4s^2 4p");
  REQUIRE(CI::CSF3({ip3, ip3, is}, 0, {}).config(true) == "4s+4p+^2");

  // Exchange of electrons 2, 3: on the parent space, A = (1 - 2 P23)/3 is a
  // projector of rank N_CSF, so Tr(P23) = (N_parents - 3 N_CSF)/2.
  // (j)^3, j = 5/2: parents [jj]J2 (J2 even) coupled with j to J
  for (const int twoJ : {3, 5, 9}) {
    std::vector<CI::CSF3::Parent> parents;
    for (int twoJ2 = 0; twoJ2 <= 8; twoJ2 += 4) {
      if (Angular::triangle(twoJ2, d5.twoj(), twoJ) != 0) {
        parents.push_back({{id5, id5}, id5, twoJ2, 0.0});
      }
    }
    double trace = 0.0;
    for (const auto &p : parents) {
      trace += CI::parent_exchange_23(p, p, twoJ);
      for (const auto &q : parents) {
        REQUIRE(CI::parent_exchange_23(p, q, twoJ) ==
                Approx(CI::parent_exchange_23(q, p, twoJ)).margin(1.0e-14));
      }
    }
    const auto num_CSFs = CI::form_CSF3s({id5, id5, id5}, twoJ).size();
    const auto expected =
      0.5 * (double(parents.size()) - 3.0 * double(num_CSFs));
    REQUIRE(trace == Approx(expected).margin(1.0e-14));
  }
  // Distinct orbitals: P23 moves z out of the third position: no diagonal
  const CI::CSF3::Parent ab_c{{is, ip1}, ip3, 2, 0.0};
  REQUIRE(CI::parent_exchange_23(ab_c, ab_c, 3) == 0.0);

  // PsiJPi accessors, for two and three electrons
  const std::vector<DiracSpinor> basis{s, p1, p3, d5};
  const CI::PsiJPi psi3(3, +1, basis, 3);
  REQUIRE(psi3.num_electrons() == 3);
  REQUIRE(psi3.num_CSFs() == CI::form_CSF3s(3, +1, basis).size());
  REQUIRE(psi3.CSF3s().size() == psi3.num_CSFs());
  REQUIRE(psi3.csf(0).num_electrons() == 3);
  REQUIRE(psi3.csf(1).config(true) == psi3.CSF3s()[1].config(true));
  const CI::PsiJPi psi2(2, +1, basis);
  REQUIRE(psi2.num_electrons() == 2);
  REQUIRE(psi2.num_CSFs() == psi2.CSFs().size());
  REQUIRE(psi2.CSF(0) == psi2.CSFs()[0]);
  REQUIRE(psi2.csf(0).num_electrons() == 2);

  // Leading configuration: weights summed over the relativistic CSFs of each
  // non-rel configuration. J = 3/2: two 4s4p-4p+ CSFs (c = 0.5 each) beat
  // the single 4s^2 4p CSF with the largest coefficient (c = 0.6)
  std::vector<CI::CSF3> csfs = CI::form_CSF3s({is, ip1, ip3}, 3);
  REQUIRE(csfs.size() == 2);
  csfs.push_back(CI::form_CSF3s({is, is, ip3}, 3).at(0));
  const std::vector<double> coefs{0.5, 0.5, 0.6};
  const LinAlg::View<const double> view{coefs.data(), 0, coefs.size(), 1};
  const auto [config, weight] = CI::leading_config(view, csfs);
  REQUIRE(config == "4s4p^2");
  REQUIRE(weight == Approx(0.5));

  // dSigma/dE correction with two spectators: the spectator energy is the
  // sum, so it equals the one-spectator form with that energy
  CI::Sigma1Correction s1c;
  s1c.S1.add(is, is, -0.01);
  s1c.dS1.add(is, is, -0.002);
  s1c.e_sigma[s.kappa()] = -0.2;
  s1c.en[ip1] = -0.10;
  s1c.en[ip3] = -0.15;
  s1c.en[id5] = -0.25;
  s1c.E0 = -0.8;
  const auto two_spectators = s1c.delta_h1(is, is, ip1, ip3);
  REQUIRE(two_spectators != 0.0);
  REQUIRE(two_spectators == Approx(s1c.delta_h1(is, is, id5)));
  // Unknown spectator: no correction
  REQUIRE(s1c.delta_h1(is, is, ip1, s.nk_index()) == 0.0);

  // L^2, S^2 between two-electron CSFs (non-rel limit): 4p- 4p+, J=1 is pure
  // 3P1; different non-rel configurations do not connect
  const auto [L2_pp, S2_pp] = CI::L2S2_pair(ip1, ip3, ip1, ip3, 2);
  REQUIRE(L2_pp == Approx(2.0));
  REQUIRE(S2_pp == Approx(2.0));
  const auto [L2_sd, S2_sd] = CI::L2S2_pair(is, ip1, is, id5, 2);
  REQUIRE(L2_sd == 0.0);
  REQUIRE(S2_sd == 0.0);

  // Nearest term {2S, L}: S has the parity of J
  REQUIRE(CI::Term_S_L_from_expectation(0.0, 0.75, 1) == std::pair{1, 0});
  REQUIRE(CI::Term_S_L_from_expectation(2.0, 3.75, 3) == std::pair{3, 1});
  REQUIRE(CI::Term_S_L_from_expectation(2.0, 2.0, 2) == std::pair{2, 1});
  REQUIRE(CI::Term_S_L_from_expectation(6.0, 0.0, 4) == std::pair{0, 2});
}
