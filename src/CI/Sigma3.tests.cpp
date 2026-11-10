#include "CI/Sigma3.hpp"
#include "Angular/include.hpp"
#include "CI/CI_Integrals.hpp"
#include "CI/CSF.hpp"
#include "CI/ConfigurationInteraction.hpp"
#include "Coulomb/include.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/Sigma2.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace UnitTest {
namespace Sigma3Test {

// Slater determinants. The spin-orbitals (orbital, 2m) are numbered in order,
// core first; a determinant is a sorted list of spin-orbital numbers,
// |D> = a+_D[0] a+_D[1] ... |0>, and includes the core. A state is a map from
// determinants to amplitudes
struct SpinOrbital {
  // Index into Atom::all
  std::size_t orbital;
  int twom;
};
using Determinant = std::vector<std::size_t>;
using Expansion = std::map<Determinant, double>;

// Al: [Ne] core, valence orbitals 3s, 3p, 3d, 4s, 4p; Q^k over all of them.
// Some of the three-electron CSFs of J = 3/2, odd parity
struct Atom {
  std::vector<DiracSpinor> core{}, valence{}, all{};
  Coulomb::QkTable qk{};
  std::vector<SpinOrbital> so{};
  // Number of core spin-orbitals (the first ones)
  std::size_t num_core{0};
  // Energy of each spin-orbital
  std::vector<double> en{};
  // <ij|g|kl> between spin-orbitals; see g()
  std::vector<double> g_table{};
  int twoJ{3};
  // Every CSF of the block, and a sample of them
  std::vector<CI::CSF> block{};
  std::vector<CI::CSF3> csfs{};

  double g(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    const auto n = so.size();
    return g_table[((i * n + j) * n + k) * n + l];
  }
};
Atom atom();

// a+_i |D> and a_i |D>: the sign (zero if not possible) and the result
std::pair<double, Determinant> create(std::size_t i, Determinant D);
std::pair<double, Determinant> annihilate(std::size_t i, Determinant D);

// The CSF (M = J) as determinants, from its parent expansion
Expansion expand(const CI::CSF3 &csf, const Atom &atom);

// <A|B>
double dot(const Expansion &A, const Expansion &B);

// Sum of the orbital energies of the determinant
double energy(const Determinant &D, const Atom &atom);

// V|state>, V = (1/2) sum g_pqrs a+_p a+_q a_s a_r, kept only for the
// determinants with exactly one core hole
Expansion V_one_hole(const Expansion &state, const Atom &atom);

// Sigma_3 |state>, with the operator in second-quantised form:
// -sum g_xavw g_yzau / D  a+_x a+_y a+_z a_w a_v a_u.
// Denominators RS (averaged with the bra-ket partner) or BW (with E0)
Expansion Sigma3_on(const Expansion &state, const Atom &atom,
                    MBPT::Denominators denominators, double E0);

} // namespace Sigma3Test
} // namespace UnitTest

//==============================================================================
// When all three orbitals of two CSFs differ, the one- and two-body parts of
// second-order perturbation theory vanish: the whole second-order matrix
// element through the determinants with one core hole is Sigma_3
TEST_CASE("Sigma3: second-order perturbation theory",
          "[CI][Sigma3][3VCI][unit]") {
  using namespace UnitTest::Sigma3Test;
  const auto at = atom();
  const auto orbitals = CI::sigma3_orbitals(at.core, at.valence);

  double worst = 0.0;
  double largest = 0.0;
  int num_pairs = 0;
  for (const auto &A : at.csfs) {
    const auto A_dets = expand(A, at);
    const auto VA = V_one_hole(A_dets, at);
    const auto E_A = energy(A_dets.begin()->first, at);
    for (const auto &B : at.csfs) {
      if (CI::CSF3::num_different(A, B) != 3)
        continue;
      const auto B_dets = expand(B, at);
      const auto VB = V_one_hole(B_dets, at);
      const auto E_B = energy(B_dets.begin()->first, at);

      // (1/2) [1/(E_A - E_Q) + 1/(E_B - E_Q)] <A|V|Q><Q|V|B>
      double expected = 0.0;
      for (const auto &[Q, amp_B] : VB) {
        const auto it = VA.find(Q);
        if (it == VA.end())
          continue;
        const auto E_Q = energy(Q, at);
        expected +=
          0.5 * (1.0 / (E_A - E_Q) + 1.0 / (E_B - E_Q)) * it->second * amp_B;
      }

      const auto found =
        CI::Sigma3_AB(A, B, at.twoJ, at.qk, orbitals, MBPT::Denominators::RS);
      worst = std::max(worst, std::abs(found - expected));
      largest = std::max(largest, std::abs(expected));
      ++num_pairs;
    }
  }
  REQUIRE(num_pairs > 10);
  REQUIRE(largest > 1.0e-5);
  REQUIRE(worst < 1.0e-12);
}

//==============================================================================
// Any pair of CSFs (repeated orbitals, diagonal): against the second-quantised
// operator acting on the determinants
TEST_CASE("Sigma3: operator on determinants", "[CI][Sigma3][3VCI][unit]") {
  using namespace UnitTest::Sigma3Test;
  const auto at = atom();
  const auto orbitals = CI::sigma3_orbitals(at.core, at.valence);

  std::vector<Expansion> dets;
  for (const auto &A : at.csfs) {
    dets.push_back(expand(A, at));
  }

  // Each CSF is normalised
  double worst_norm = 0.0;
  for (const auto &A_dets : dets) {
    worst_norm = std::max(worst_norm, std::abs(dot(A_dets, A_dets) - 1.0));
  }
  REQUIRE(worst_norm < 1.0e-12);

  const double E0 = 3.0 * orbitals.en_min - 0.1;
  for (const auto denominators :
       {MBPT::Denominators::RS, MBPT::Denominators::BW}) {
    INFO(MBPT::parse_Denominators(denominators));
    double worst = 0.0;
    double largest = 0.0;
    for (std::size_t j = 0; j < at.csfs.size(); ++j) {
      const auto S3_B = Sigma3_on(dets[j], at, denominators, E0);
      for (std::size_t i = 0; i < at.csfs.size(); ++i) {
        const auto expected = dot(dets[i], S3_B);
        const auto found = CI::Sigma3_AB(at.csfs[i], at.csfs[j], at.twoJ, at.qk,
                                         orbitals, denominators, E0);
        worst = std::max(worst, std::abs(found - expected));
        largest = std::max(largest, std::abs(expected));
      }
    }
    REQUIRE(largest > 1.0e-5);
    REQUIRE(worst < 1.0e-12);
  }
}

//==============================================================================
// The matrix is symmetric for every choice of denominators
TEST_CASE("Sigma3: symmetric", "[CI][Sigma3][3VCI][unit]") {
  using namespace UnitTest::Sigma3Test;
  const auto at = atom();
  const auto orbitals = CI::sigma3_orbitals(at.core, at.valence);
  const double E0 = 3.0 * orbitals.en_min - 0.1;

  for (const auto denominators :
       {MBPT::Denominators::RS, MBPT::Denominators::DFK, MBPT::Denominators::BW,
        MBPT::Denominators::Fermi, MBPT::Denominators::Fermi0}) {
    INFO(MBPT::parse_Denominators(denominators));
    double worst = 0.0;
    double largest = 0.0;
    for (const auto &A : at.csfs) {
      for (const auto &B : at.csfs) {
        const auto S_AB =
          CI::Sigma3_AB(A, B, at.twoJ, at.qk, orbitals, denominators, E0);
        const auto S_BA =
          CI::Sigma3_AB(B, A, at.twoJ, at.qk, orbitals, denominators, E0);
        worst = std::max(worst, std::abs(S_AB - S_BA));
        largest = std::max(largest, std::abs(S_AB));
      }
    }
    REQUIRE(largest > 1.0e-5);
    REQUIRE(worst < 1.0e-12);
  }
}

//==============================================================================
// The CSFs Sigma_3 is included for: those of the lowest non-relativistic
// configurations
TEST_CASE("Sigma3: lowest configurations", "[CI][Sigma3][3VCI][unit]") {
  using namespace UnitTest::Sigma3Test;
  const auto at = atom();
  const auto orbitals = CI::sigma3_orbitals(at.core, at.valence);

  // Lowest configuration of this block: 3s^2 3p
  const auto lowest = CI::sigma3_selected(at.block, orbitals, 1);
  std::size_t num_lowest = 0;
  bool only_ground_config = true;
  for (std::size_t i = 0; i < at.block.size(); ++i) {
    const auto is_ground_config = at.block[i].config(false) == "3s^2 3p";
    if (lowest[i]) {
      ++num_lowest;
    }
    if (bool(lowest[i]) != is_ground_config) {
      only_ground_config = false;
    }
  }
  REQUIRE(only_ground_config);
  REQUIRE(num_lowest == 1);

  // Both j of each shell are selected together: a CSF is selected if any CSF
  // of the same non-relativistic configuration is
  const auto several = CI::sigma3_selected(at.block, orbitals, 4);
  std::set<std::string> selected_configs;
  for (std::size_t i = 0; i < at.block.size(); ++i) {
    if (several[i]) {
      selected_configs.insert(at.block[i].config(false));
    }
  }
  bool whole_configs = true;
  for (std::size_t i = 0; i < at.block.size(); ++i) {
    const auto in_selected_config =
      selected_configs.count(at.block[i].config(false)) != 0;
    if (in_selected_config != bool(several[i])) {
      whole_configs = false;
    }
  }
  REQUIRE(selected_configs.size() == 4);
  REQUIRE(whole_configs);

  // None, and all
  const auto none = CI::sigma3_selected(at.block, orbitals, 0);
  const auto all = CI::sigma3_selected(at.block, orbitals, 1000);
  const auto num_none = std::count(none.begin(), none.end(), 1);
  const auto num_all = std::count(all.begin(), all.end(), 1);
  REQUIRE(num_none == 0);
  REQUIRE(std::size_t(num_all) == at.block.size());
}

//==============================================================================
// CI matrix with Sigma_3: the matrix without it, plus Sigma_3 for exactly the
// pairs with a CSF of the lowest configurations (also when all three orbitals
// differ); sparse and dense matrices agree
TEST_CASE("Sigma3: in the CI matrix", "[CI][Sigma3][3VCI][unit]") {
  using namespace UnitTest::Sigma3Test;
  const auto at = atom();

  Coulomb::meTable<double> h1;
  for (const auto &v : at.valence) {
    h1.add(v, v, v.en());
  }

  CI::Sigma3Parameters sigma3;
  sigma3.orbitals = CI::sigma3_orbitals(at.core, at.valence);
  sigma3.denominators = MBPT::Denominators::DFK;
  sigma3.num_configs = 3;
  const auto selected =
    CI::sigma3_selected(at.block, sigma3.orbitals, sigma3.num_configs);

  const auto H0 = CI::construct_Hci(at.block, at.twoJ, h1, at.qk);
  const auto H3 =
    CI::construct_Hci(at.block, at.twoJ, h1, at.qk, nullptr, nullptr, nullptr,
                      {}, nullptr, 0.0, &sigma3);

  double worst = 0.0;
  double largest_three_different = 0.0;
  for (std::size_t i = 0; i < at.block.size(); ++i) {
    for (std::size_t j = 0; j < at.block.size(); ++j) {
      const auto expected =
        selected[i] || selected[j] ?
          CI::Sigma3_AB(at.block[i], at.block[j], at.twoJ, at.qk,
                        sigma3.orbitals, sigma3.denominators) :
          0.0;
      const auto found = H3(i, j) - H0(i, j);
      worst = std::max(worst, std::abs(found - expected));
      if (CI::num_different(at.block[i], at.block[j]) == 3) {
        largest_three_different =
          std::max(largest_three_different, std::abs(found));
      }
    }
  }
  REQUIRE(worst < 1.0e-14);
  REQUIRE(largest_three_different > 1.0e-6);

  const auto H3_sparse =
    CI::construct_Hci<LinAlg::SparseSymmetricMatrix<double>>(
      at.block, at.twoJ, h1, at.qk, nullptr, nullptr, nullptr, {}, nullptr, 0.0,
      &sigma3)
      .dense();
  double worst_sparse = 0.0;
  for (std::size_t i = 0; i < at.block.size(); ++i) {
    for (std::size_t j = 0; j < at.block.size(); ++j) {
      worst_sparse =
        std::max(worst_sparse, std::abs(H3_sparse(i, j) - H3(i, j)));
    }
  }
  REQUIRE(worst_sparse == 0.0);
}

//==============================================================================
// Full CI run with the sigma3 option: the levels are the eigenvalues of the
// Hamiltonian with Sigma_3, as re-formed from the stored integrals; with BW
// denominators, of the Hamiltonian at the converged E0
TEST_CASE("Sigma3: CI run", "[CI][Sigma3][3VCI][unit]") {

  // B: [He] core, three valence electrons
  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"B", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("6spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  const std::string common =
    "ci_basis = 3sp; num_electrons = 3; J- = 1/2; num_solutions = 2; "
    "sigma1 = true; sigma2 = true; davidson = false; print_details = false; "
    "qk_file = false; sk_file = false; dsk_file = false; ds1_file = false; "
    "ci_file = false;";

  for (const std::string denominators : {"DFK", "BW"}) {
    INFO(denominators);
    const auto options = common + "denominators = " + denominators + ";";

    const auto without =
      CI::configuration_interaction(IO::InputBlock{"CI", options}, wf);
    const auto with = CI::configuration_interaction(
      IO::InputBlock{"CI", options + "sigma3 = true; sigma3_configs = 2;"}, wf);

    const auto &psi = with.levels.at(0);
    const auto shift = psi.energy(0) - without.levels.at(0).energy(0);
    REQUIRE(std::abs(shift) > 1.0e-7);

    const auto Hci = CI::construct_Hci(psi, with.integrals);
    const auto E = LinAlg::symmhEigensystem(Hci).first;
    REQUIRE(psi.energy(0) == Approx(E[0]).epsilon(1.0e-9));
    REQUIRE(psi.energy(1) == Approx(E[1]).epsilon(1.0e-9));
  }
}

//==============================================================================
//==============================================================================

UnitTest::Sigma3Test::Atom UnitTest::Sigma3Test::atom() {
  Atom at;

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Al", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("4spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  at.core = wf.core();
  for (const auto &Fn : wf.basis()) {
    if (wf.isInCore(Fn.n(), Fn.kappa()))
      continue;
    if (Fn.n() == 3 || (Fn.n() == 4 && Fn.l() <= 1)) {
      at.valence.push_back(Fn);
    }
  }
  at.all = at.core;
  at.all.insert(at.all.end(), at.valence.begin(), at.valence.end());
  at.qk.fill(at.all, Coulomb::YkTable(at.all), -1, false);

  // Spin-orbitals, core first
  for (std::size_t i = 0; i < at.all.size(); ++i) {
    for (int twom = -at.all[i].twoj(); twom <= at.all[i].twoj(); twom += 2) {
      at.so.push_back({i, twom});
      at.en.push_back(at.all[i].en());
    }
    if (i + 1 == at.core.size()) {
      at.num_core = at.so.size();
    }
  }

  const auto n = at.so.size();
  at.g_table.assign(n * n * n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t l = 0; l < n; ++l) {
          // m is conserved
          if (at.so[i].twom + at.so[j].twom != at.so[k].twom + at.so[l].twom)
            continue;
          at.g_table[((i * n + j) * n + k) * n + l] =
            at.qk.g(at.all[at.so[i].orbital], at.all[at.so[j].orbital],
                    at.all[at.so[k].orbital], at.all[at.so[l].orbital],
                    at.so[i].twom, at.so[j].twom, at.so[k].twom, at.so[l].twom);
        }
      }
    }
  }

  // A sample of the CSFs: every seventh, so that many configurations (with
  // and without repeated orbitals) are present
  at.block = CI::form_CSF3s(at.twoJ, -1, at.valence);
  for (std::size_t i = 0; i < at.block.size(); i += 7) {
    at.csfs.push_back(std::get<CI::CSF3>(at.block[i].variant()));
  }

  return at;
}

//==============================================================================
std::pair<double, UnitTest::Sigma3Test::Determinant>
UnitTest::Sigma3Test::create(std::size_t i, Determinant D) {
  const auto pos = std::lower_bound(D.begin(), D.end(), i);
  if (pos != D.end() && *pos == i)
    return {0.0, {}};
  const auto sign = (pos - D.begin()) % 2 == 0 ? 1.0 : -1.0;
  D.insert(pos, i);
  return {sign, D};
}

std::pair<double, UnitTest::Sigma3Test::Determinant>
UnitTest::Sigma3Test::annihilate(std::size_t i, Determinant D) {
  const auto pos = std::lower_bound(D.begin(), D.end(), i);
  if (pos == D.end() || *pos != i)
    return {0.0, {}};
  const auto sign = (pos - D.begin()) % 2 == 0 ? 1.0 : -1.0;
  D.erase(pos);
  return {sign, D};
}

//==============================================================================
UnitTest::Sigma3Test::Expansion
UnitTest::Sigma3Test::expand(const CI::CSF3 &csf, const Atom &atom) {

  // Number of the valence spin-orbital (nk_index, 2m)
  const auto spin_orbital = [&atom](DiracSpinor::Index nk, int twom) {
    for (std::size_t i = atom.num_core; i < atom.so.size(); ++i) {
      if (atom.all[atom.so[i].orbital].nk_index() == nk &&
          atom.so[i].twom == twom)
        return i;
    }
    return atom.so.size();
  };

  // Amplitudes of the products |i>_1 |j>_2 |k>_3
  std::map<std::array<std::size_t, 3>, double> products;
  for (const auto &parent : csf.parents()) {
    const auto [x, y] = parent.pair;
    const auto z = parent.spectator;
    const auto tjx = Angular::nkindex_to_twoj(x);
    const auto tjy = Angular::nkindex_to_twoj(y);
    const auto tjz = Angular::nkindex_to_twoj(z);
    const auto eta = x == y ? 0.5 : 1.0 / std::sqrt(2.0);
    for (int mx = -tjx; mx <= tjx; mx += 2) {
      for (int my = -tjy; my <= tjy; my += 2) {
        const auto twoM2 = mx + my;
        const auto mz = atom.twoJ - twoM2;
        if (std::abs(mz) > tjz)
          continue;
        const auto cg =
          Angular::cg_2(tjx, mx, tjy, my, parent.twoJ2, twoM2) *
          Angular::cg_2(parent.twoJ2, twoM2, tjz, mz, atom.twoJ, atom.twoJ);
        if (cg == 0.0)
          continue;
        const auto i = spin_orbital(x, mx);
        const auto j = spin_orbital(y, my);
        const auto k = spin_orbital(z, mz);
        products[{i, j, k}] += parent.coef * eta * cg;
        products[{j, i, k}] -= parent.coef * eta * cg;
      }
    }
  }

  // Antisymmetric: the determinant i < j < k has amplitude sqrt(3!) times
  // that of the product
  Determinant core(atom.num_core);
  for (std::size_t c = 0; c < atom.num_core; ++c) {
    core[c] = c;
  }
  Expansion out;
  for (const auto &[ijk, amplitude] : products) {
    if (ijk[0] < ijk[1] && ijk[1] < ijk[2] && std::abs(amplitude) > 1.0e-14) {
      auto D = core;
      D.insert(D.end(), ijk.begin(), ijk.end());
      out[D] = std::sqrt(6.0) * amplitude;
    }
  }
  return out;
}

//==============================================================================
double UnitTest::Sigma3Test::dot(const Expansion &A, const Expansion &B) {
  double sum = 0.0;
  for (const auto &[D, amplitude] : A) {
    const auto it = B.find(D);
    if (it != B.end()) {
      sum += amplitude * it->second;
    }
  }
  return sum;
}

double UnitTest::Sigma3Test::energy(const Determinant &D, const Atom &atom) {
  double sum = 0.0;
  for (const auto i : D) {
    sum += atom.en[i];
  }
  return sum;
}

//==============================================================================
UnitTest::Sigma3Test::Expansion
UnitTest::Sigma3Test::V_one_hole(const Expansion &state, const Atom &atom) {
  const auto n = atom.so.size();
  Expansion out;
  for (const auto &[D, amplitude] : state) {
    for (std::size_t ir = 0; ir < D.size(); ++ir) {
      for (std::size_t is = ir + 1; is < D.size(); ++is) {
        const auto r = D[ir];
        const auto s = D[is];
        // Exactly one of the two removed is a core electron
        if (r >= atom.num_core || s < atom.num_core)
          continue;
        const auto [sign_r, D_r] = annihilate(r, D);
        const auto [sign_s, D_rs] = annihilate(s, D_r);
        // Both created are valence
        for (std::size_t p = atom.num_core; p < n; ++p) {
          for (std::size_t q = p + 1; q < n; ++q) {
            const auto g = atom.g(p, q, r, s) - atom.g(p, q, s, r);
            if (g == 0.0)
              continue;
            // a+_p a+_q a_s a_r
            const auto [sign_q, D_q] = create(q, D_rs);
            if (sign_q == 0.0)
              continue;
            const auto [sign_p, D_pq] = create(p, D_q);
            if (sign_p == 0.0)
              continue;
            out[D_pq] += amplitude * sign_r * sign_s * sign_q * sign_p * g;
          }
        }
      }
    }
  }
  return out;
}

//==============================================================================
UnitTest::Sigma3Test::Expansion
UnitTest::Sigma3Test::Sigma3_on(const Expansion &state, const Atom &atom,
                                MBPT::Denominators denominators, double E0) {
  const auto n = atom.so.size();
  const auto &en = atom.en;
  Expansion out;
  for (const auto &[D, amplitude] : state) {
    // The three valence spin-orbitals of the determinant
    const std::array<std::size_t, 3> valence{
      D[atom.num_core], D[atom.num_core + 1], D[atom.num_core + 2]};
    // Every ordered choice (v, w, u)
    for (const auto v : valence) {
      for (const auto w : valence) {
        for (const auto u : valence) {
          if (v == w || v == u || w == u)
            continue;
          // a_w a_v a_u
          const auto [sign_u, D_u] = annihilate(u, D);
          const auto [sign_v, D_uv] = annihilate(v, D_u);
          const auto [sign_w, D_core] = annihilate(w, D_uv);
          const auto sign_vwu = sign_u * sign_v * sign_w;
          const auto twom = atom.so[v].twom + atom.so[w].twom + atom.so[u].twom;

          // a+_x a+_y a+_z
          for (std::size_t z = atom.num_core; z < n; ++z) {
            const auto [sign_z, D_z] = create(z, D_core);
            for (std::size_t y = atom.num_core; y < n; ++y) {
              const auto [sign_y, D_yz] = create(y, D_z);
              if (sign_y == 0.0)
                continue;
              for (std::size_t x = atom.num_core; x < n; ++x) {
                if (atom.so[x].twom + atom.so[y].twom + atom.so[z].twom != twom)
                  continue;
                const auto [sign_x, D_xyz] = create(x, D_yz);
                if (sign_x == 0.0)
                  continue;

                double sum = 0.0;
                for (std::size_t a = 0; a < atom.num_core; ++a) {
                  const auto gg = atom.g(x, a, v, w) * atom.g(y, z, a, u);
                  if (gg == 0.0)
                    continue;
                  const auto e_int = en[v] + en[w] + en[y] + en[z] - en[a];
                  const auto inv_de =
                    denominators == MBPT::Denominators::BW ?
                      1.0 / (E0 - e_int) :
                      0.5 * (1.0 / (en[v] + en[w] + en[u] - e_int) +
                             1.0 / (en[x] + en[y] + en[z] - e_int));
                  sum += gg * inv_de;
                }
                if (sum != 0.0) {
                  out[D_xyz] -=
                    amplitude * sign_vwu * sign_z * sign_y * sign_x * sum;
                }
              }
            }
          }
        }
      }
    }
  }
  return out;
}
