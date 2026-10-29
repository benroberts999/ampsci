#include "CSF.hpp"
#include "Angular/include.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "fmt/format.hpp"
#include "fmt/ostream.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <tuple>
#include <type_traits>
#include <variant>

namespace CI {

//==============================================================================
CSF1::CSF1(const DiracSpinor &a) : states{a.nk_index()}, m_parity(a.parity()) {}

DiracSpinor::Index CSF1::state(std::size_t i) const { return states.at(i); }

bool operator==(const CSF1 &A, const CSF1 &B) { return A.states == B.states; }

bool operator!=(const CSF1 &A, const CSF1 &B) { return !(A == B); }

int CSF1::num_different(const CSF1 &A, const CSF1 &B) { return A == B ? 0 : 1; }

std::string CSF1::config(bool relativistic) const {
  const auto [n, kappa] = Angular::index_to_nk(states[0]);
  auto symbol = DiracSpinor::shortSymbol(n, kappa);
  if (!relativistic) {
    symbol.pop_back();
  }
  return symbol;
}

//==============================================================================
// Forms list of all possible (1-particle) CSFs with given J and parity: the
// orbitals with j = J and that parity
std::vector<CSF> form_CSF1s(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis,
                            const ConfigurationSpace &config_space) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ > 0 && twoJ % 2 == 1 && "2*J must be odd for one electron");

  std::vector<CSF> CSFs;
  for (const auto &v : cisp_basis) {
    if (v.twoj() != twoJ || v.parity() != parity)
      continue;
    if (!config_space.included({v.nk_index()}))
      continue;
    CSFs.emplace_back(CSF1{v});
  }
  return CSFs;
}

//==============================================================================
CSF2::CSF2(const DiracSpinor &a, const DiracSpinor &b)
  : states(a <= b ? std::array{a.nk_index(), b.nk_index()} :
                    std::array{b.nk_index(), a.nk_index()}),
    m_parity(a.parity() * b.parity()) {}

DiracSpinor::Index CSF2::state(std::size_t i) const { return states.at(i); }

bool operator==(const CSF2 &A, const CSF2 &B) {
  // only works because states is sorted
  return A.states == B.states;
}

bool operator!=(const CSF2 &A, const CSF2 &B) {
  // only works because states is sorted
  return !(A.states == B.states);
}

// Returns number of different orbitals between two CSFs
int CSF2::num_different(const CSF2 &A, const CSF2 &B) {
  if (A == B)
    return 0;
  if (A.state(0) == B.state(0) || A.state(1) == B.state(1) || //
      A.state(0) == B.state(1) || A.state(1) == B.state(0))
    return 1;
  return 2;
}

// returns _different_ orbitals, for case where CSFs differ by 1.
// i.e., returns {n,a} where |A> = |B_a^n> (i.e., A has n, but not a)
std::array<DiracSpinor::Index, 2> CSF2::diff_1_na(const CSF2 &A,
                                                  const CSF2 &B) {
  assert(num_different(A, B) == 1); // only valid in this case
  if (A.state(0) == B.state(0))
    return {A.state(1), B.state(1)};
  if (A.state(1) == B.state(1))
    return {A.state(0), B.state(0)};
  if (A.state(0) == B.state(1))
    return {A.state(1), B.state(0)};
  if (A.state(1) == B.state(0))
    return {A.state(0), B.state(1)};
  assert(false); // should be unreachable, for testing
}

DiracSpinor::Index CSF2::same_1_j(const CSF2 &A, const CSF2 &B) {
  assert(num_different(A, B) == 1); // only valid in this case
  if (A.state(0) == B.state(0))
    return A.state(0);
  if (A.state(1) == B.state(1))
    return A.state(1);
  if (A.state(0) == B.state(1))
    return A.state(0);
  if (A.state(1) == B.state(0))
    return A.state(1);
  assert(false); // should be unreachable, for testing
}

std::string CSF2::config(bool relativistic) const {
  // auto s1 = states[0]->shortSymbol();
  // auto s2 = states[1]->shortSymbol();
  const auto [na, ka] = Angular::index_to_nk(states[0]);
  const auto [nb, kb] = Angular::index_to_nk(states[1]);
  auto s1 = DiracSpinor::shortSymbol(na, ka);
  auto s2 = DiracSpinor::shortSymbol(nb, kb);
  if (!relativistic) {
    s1.pop_back();
    s2.pop_back();
  }
  return s1 == s2 ? s1 + "^2" : s1 + s2;
}

//==============================================================================
// Parses "2+:3" or "e2:3" (the index is optional)
std::optional<Level> parse_level(std::string_view str) {

  if (str.empty())
    return std::nullopt;

  // Reads the digits of str[first,last), returns empty if there are none, or
  // if any character is not a digit
  const auto read_int = [&str](std::size_t first,
                               std::size_t last) -> std::optional<int> {
    if (last <= first)
      return std::nullopt;
    int value = 0;
    for (auto i = first; i < last; ++i) {
      if (!std::isdigit(str.at(i)))
        return std::nullopt;
      value = 10 * value + (str.at(i) - '0');
    }
    return value;
  };

  // Leading 'e' or 'o' form
  int parity = 0;
  std::size_t first = 0;
  if (str.front() == 'e' || str.front() == 'E') {
    parity = 1;
    first = 1;
  } else if (str.front() == 'o' || str.front() == 'O') {
    parity = -1;
    first = 1;
  }

  // Optional ':index'
  std::size_t index = 0;
  auto last = str.size();
  const auto colon = str.find(':');
  if (colon != std::string_view::npos) {
    const auto t_index = read_int(colon + 1, str.size());
    if (!t_index)
      return std::nullopt;
    index = std::size_t(*t_index);
    last = colon;
  }

  // Trailing '+' or '-' form
  if (parity == 0) {
    if (last == 0)
      return std::nullopt;
    const auto pi_char = str.at(last - 1);
    if (pi_char == '+') {
      parity = 1;
    } else if (pi_char == '-') {
      parity = -1;
    } else {
      return std::nullopt;
    }
    --last;
  }

  const auto twoJ = parse_twoJ(str.substr(first, last - first));
  if (!twoJ)
    return std::nullopt;

  return Level{*twoJ, parity, index};
}

//==============================================================================
std::string to_string(const Level &level) {
  return J_string(level.twoJ) + (level.parity == 1 ? "+" : "-") + ":" +
         std::to_string(level.index);
}

//==============================================================================
std::string J_string(int twoJ) {
  return twoJ % 2 == 0 ? std::to_string(twoJ / 2) : std::to_string(twoJ) + "/2";
}

//==============================================================================
std::optional<int> parse_twoJ(std::string_view str) {
  // trim whitespace
  while (!str.empty() && std::isspace(str.front())) {
    str.remove_prefix(1);
  }
  while (!str.empty() && std::isspace(str.back())) {
    str.remove_suffix(1);
  }
  if (str.empty())
    return std::nullopt;

  // Digits only, else empty
  const auto read_int = [](std::string_view digits) -> std::optional<int> {
    if (digits.empty())
      return std::nullopt;
    int value = 0;
    for (const auto c : digits) {
      if (!std::isdigit(c))
        return std::nullopt;
      value = 10 * value + (c - '0');
    }
    return value;
  };

  // Decimal: must be an integer or half-integer
  const auto point = str.find('.');
  if (point != std::string_view::npos) {
    const auto whole = read_int(str.substr(0, point));
    const auto fraction = str.substr(point + 1);
    if (!whole || fraction.empty())
      return std::nullopt;
    // .0 or .5 only (any number of trailing zeros)
    const auto first = fraction.front();
    const auto rest = fraction.substr(1);
    if ((first != '0' && first != '5') ||
        rest.find_first_not_of('0') != std::string_view::npos) {
      return std::nullopt;
    }
    return 2 * *whole + (first == '5' ? 1 : 0);
  }

  const auto slash = str.find('/');
  if (slash == std::string_view::npos) {
    const auto J = read_int(str);
    return J ? std::optional{2 * *J} : std::nullopt;
  }
  // Fraction: numerator / 2 only
  const auto numerator = read_int(str.substr(0, slash));
  const auto denominator = read_int(str.substr(slash + 1));
  if (!numerator || !denominator || *denominator != 2)
    return std::nullopt;
  return *numerator;
}

//==============================================================================
// Forms list of all possible (2-particle) CSFs with given J and parity
std::vector<CSF> form_CSF2s(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis,
                            const ConfigurationSpace &config_space) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ >= 0 && "2*J must not be negative");

  std::vector<CSF> CSFs;

  // This is always smaller than actual # of CSFs, but good head-start
  CSFs.reserve(cisp_basis.size());

  for (const auto &v : cisp_basis) {
    for (const auto &w : cisp_basis) {

      // Symmetric: only include unique CSFs once
      if (w < v)
        continue;

      // Parity symmetry:
      if (v.parity() * w.parity() != parity)
        continue;

      // J triangle rule (use M=Jz=J):
      if (v.twoj() + w.twoj() < twoJ || std::abs(v.twoj() - w.twoj()) > twoJ)
        continue;

      // identical particles can only give even J (cannot have mv=mw)
      if (v == w && twoJ % 4 != 0)
        continue;

      if (!config_space.included({v.nk_index(), w.nk_index()}))
        continue;

      CSFs.emplace_back(CSF2{v, w});
    }
  }

  return CSFs;
}

//==============================================================================
CSF3::CSF3(const std::array<DiracSpinor::Index, 3> &orbitals, int twoJ12,
           std::vector<Parent> parents)
  : states(orbitals),
    m_twoJ12(twoJ12),
    m_parents(std::move(parents)),
    m_parity(Angular::neg1pow(Angular::nkindex_to_l(orbitals[0]) +
                              Angular::nkindex_to_l(orbitals[1]) +
                              Angular::nkindex_to_l(orbitals[2]))) {
  std::sort(states.begin(), states.end());
}

DiracSpinor::Index CSF3::state(std::size_t i) const { return states.at(i); }

int CSF3::twoJ12() const { return m_twoJ12; }

const std::vector<CSF3::Parent> &CSF3::parents() const { return m_parents; }

//------------------------------------------------------------------------------
int CSF3::num_different(const CSF3 &A, const CSF3 &B) {
  // Count common orbitals, with repeats; works because states are sorted
  int num_common = 0;
  std::size_t i = 0;
  std::size_t j = 0;
  while (i < A.states.size() && j < B.states.size()) {
    if (A.states[i] == B.states[j]) {
      ++num_common;
      ++i;
      ++j;
    } else if (A.states[i] < B.states[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  return 3 - num_common;
}

bool operator==(const CSF3 &A, const CSF3 &B) {
  // only works because states is sorted
  return A.states == B.states && A.m_twoJ12 == B.m_twoJ12;
}

bool operator!=(const CSF3 &A, const CSF3 &B) { return !(A == B); }

std::string CSF3::config(bool relativistic) const {
  std::array<std::string, 3> symbols;
  for (std::size_t i = 0; i < states.size(); ++i) {
    const auto [n, kappa] = Angular::index_to_nk(states[i]);
    symbols[i] = DiracSpinor::shortSymbol(n, kappa);
    if (!relativistic) {
      symbols[i].pop_back();
    }
  }

  // Orbitals of one shell are adjacent, since states is sorted
  std::string out;
  for (std::size_t i = 0; i < symbols.size();) {
    std::size_t count = 1;
    while (i + count < symbols.size() && symbols[i + count] == symbols[i]) {
      ++count;
    }
    out += symbols[i];
    if (count > 1) {
      out += "^" + std::to_string(count);
      // a power would otherwise run into the n of the next orbital
      if (i + count < symbols.size()) {
        out += " ";
      }
    }
    i += count;
  }
  return out;
}

//==============================================================================
double parent_exchange_23(const CSF3::Parent &bra, const CSF3::Parent &ket,
                          int twoJ) {

  // Each parent in terms of states with electrons 1, 2, 3 in definite
  // orbitals {o1, o2, o3}, electrons 1 and 2 coupled to J2:
  //   |[xy]J2, z> = (|(xy)J2, z> - (-1)^(jx+jy-J2) |(yx)J2, z>) / sqrt(2),
  // or |(xx)J2, z> if x = y (J2 even).
  // Returns the terms {o1, o2, o3, coefficient}, and how many there are
  using Term = std::tuple<DiracSpinor::Index, DiracSpinor::Index,
                          DiracSpinor::Index, double>;
  const auto expand = [](const CSF3::Parent &p) {
    const auto [x, y] = p.pair;
    if (x == y) {
      return std::pair{std::array{Term{x, x, p.spectator, 1.0}, Term{}}, 1};
    }
    const auto eta = 1.0 / std::sqrt(2.0);
    const auto s = Angular::neg1pow_2(Angular::nkindex_to_twoj(x) +
                                      Angular::nkindex_to_twoj(y) - p.twoJ2);
    return std::pair{std::array{Term{x, y, p.spectator, eta},
                                Term{y, x, p.spectator, -s * eta}},
                     2};
  };

  const auto [bra_terms, num_bra] = expand(bra);
  const auto [ket_terms, num_ket] = expand(ket);

  double out = 0.0;
  for (int ib = 0; ib < num_bra; ++ib) {
    const auto &[b1, b2, b3, c_bra] = bra_terms.at(std::size_t(ib));
    for (int ik = 0; ik < num_ket; ++ik) {
      const auto &[k1, k2, k3, c_ket] = ket_terms.at(std::size_t(ik));
      // P23 |(k1 k2)J12, k3> = sum_K R_K |(k1 k3)K, k2>: only the K equal to
      // the bra coupling survives, and only if the orbitals match
      if (b1 != k1 || b2 != k3 || b3 != k2)
        continue;
      const auto tj1 = Angular::nkindex_to_twoj(k1);
      const auto tj2 = Angular::nkindex_to_twoj(k2);
      const auto tj3 = Angular::nkindex_to_twoj(k3);
      const auto twoK = bra.twoJ2;
      const auto twoJ12 = ket.twoJ2;
      const auto recouple = Angular::neg1pow_2(tj3 + tj2 + twoK + twoJ12) *
                            std::sqrt(double((twoK + 1) * (twoJ12 + 1))) *
                            Angular::sixj_2(tj3, tj1, twoK, tj2, twoJ, twoJ12);
      out += c_bra * c_ket * recouple;
    }
  }
  return out;
}

//==============================================================================
std::vector<CSF3> form_CSF3s(std::array<DiracSpinor::Index, 3> orbitals,
                             int twoJ) {

  std::sort(orbitals.begin(), orbitals.end());
  const auto [a, b, c] = orbitals;

  // Order in which orbitals are assigned to the third electron. The standard
  // couplings, ((ab)J12, c), or ((bc)J12, a) if only b = c (equivalent
  // electrons coupled first), come first: these are then the CSFs formed,
  // except for j^3, where Gram-Schmidt removes the linear dependence
  const auto spectator_order =
    (b == c && a != b) ? std::array{0, 2, 1} : std::array{2, 0, 1};

  // Parent basis: states antisymmetric in electrons 1 and 2
  std::vector<CSF3::Parent> parent_basis;
  for (const auto s : spectator_order) {
    const auto z = orbitals.at(std::size_t(s));
    const auto pair = s == 0 ? std::array{b, c} :
                      s == 1 ? std::array{a, c} :
                               std::array{a, b};
    // If orbitals coincide, the same (pair, z) arises more than once
    const auto seen = std::any_of(
      parent_basis.cbegin(), parent_basis.cend(),
      [&pair, z](const auto &p) { return p.pair == pair && p.spectator == z; });
    if (seen)
      continue;

    const auto tjx = Angular::nkindex_to_twoj(pair[0]);
    const auto tjy = Angular::nkindex_to_twoj(pair[1]);
    const auto tjz = Angular::nkindex_to_twoj(z);
    for (int twoJ2 = std::abs(tjx - tjy); twoJ2 <= tjx + tjy; twoJ2 += 2) {
      // identical particles can only give even J2
      if (pair[0] == pair[1] && twoJ2 % 4 != 0)
        continue;
      if (Angular::triangle(twoJ2, tjz, twoJ) == 0)
        continue;
      parent_basis.push_back({pair, z, twoJ2, 0.0});
    }
  }

  // Antisymmetriser, on the space antisymmetric in electrons 1 and 2:
  // A = (1 - P13 - P23)/3 = (1 - 2 P23)/3 there. It is a projector
  const auto num_parents = parent_basis.size();
  if (num_parents == 0)
    return {};
  LinAlg::Matrix<double> A(num_parents, num_parents);
  for (std::size_t i = 0; i < num_parents; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      const auto delta = i == j ? 1.0 : 0.0;
      const auto Aij =
        (delta -
         2.0 * parent_exchange_23(parent_basis[i], parent_basis[j], twoJ)) /
        3.0;
      A(i, j) = Aij;
      A(j, i) = Aij;
    }
  }

  // Gram-Schmidt on the columns A|p>, which are all antisymmetric
  std::vector<CSF3> CSFs;
  std::vector<std::vector<double>> orthonormal;
  for (std::size_t s = 0; s < num_parents; ++s) {
    std::vector<double> v(num_parents);
    for (std::size_t i = 0; i < num_parents; ++i) {
      v[i] = A(i, s);
    }
    for (const auto &u : orthonormal) {
      double overlap = 0.0;
      for (std::size_t i = 0; i < num_parents; ++i) {
        overlap += u[i] * v[i];
      }
      for (std::size_t i = 0; i < num_parents; ++i) {
        v[i] -= overlap * u[i];
      }
    }
    double norm2 = 0.0;
    for (const auto x : v) {
      norm2 += x * x;
    }
    // |A|p>|^2 = <p|A|p>, which is of order 1 unless A|p> depends linearly on
    // those already found
    if (norm2 < 1.0e-10)
      continue;
    const auto norm = std::sqrt(norm2);
    for (auto &x : v) {
      x /= norm;
    }

    std::vector<CSF3::Parent> parents;
    for (std::size_t i = 0; i < num_parents; ++i) {
      if (std::abs(v[i]) > 1.0e-12) {
        parents.push_back(parent_basis[i]);
        parents.back().coef = v[i];
      }
    }
    CSFs.emplace_back(orbitals, parent_basis[s].twoJ2, std::move(parents));
    orthonormal.push_back(std::move(v));
  }

  return CSFs;
}

//==============================================================================
std::vector<CSF> form_CSF3s(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis,
                            const ConfigurationSpace &config_space) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ > 0 && twoJ % 2 == 1 && "2*J must be odd for three electrons");

  std::vector<CSF> CSFs;
  for (std::size_t i = 0; i < cisp_basis.size(); ++i) {
    const auto &u = cisp_basis[i];
    for (std::size_t j = i; j < cisp_basis.size(); ++j) {
      const auto &v = cisp_basis[j];
      for (std::size_t k = j; k < cisp_basis.size(); ++k) {
        const auto &w = cisp_basis[k];

        if (u.parity() * v.parity() * w.parity() != parity)
          continue;
        if (u.twoj() + v.twoj() + w.twoj() < twoJ)
          continue;
        if (!config_space.included({u.nk_index(), v.nk_index(), w.nk_index()}))
          continue;

        for (auto &csf :
             form_CSF3s({u.nk_index(), v.nk_index(), w.nk_index()}, twoJ)) {
          CSFs.emplace_back(std::move(csf));
        }
      }
    }
  }

  return CSFs;
}

//==============================================================================
CSF2p1h::CSF2p1h(DiracSpinor::Index v, DiracSpinor::Index w,
                 DiracSpinor::Index a, int twoJ12)
  : states(v <= w ? std::array{v, w} : std::array{w, v}),
    hole(a),
    m_twoJ12(twoJ12),
    m_parity(Angular::neg1pow(Angular::nkindex_to_l(v) +
                              Angular::nkindex_to_l(w) +
                              Angular::nkindex_to_l(a))) {}

DiracSpinor::Index CSF2p1h::state(std::size_t i) const { return states.at(i); }

int CSF2p1h::twoJ12() const { return m_twoJ12; }

bool operator==(const CSF2p1h &A, const CSF2p1h &B) {
  // only works because states is sorted
  return A.states == B.states && A.hole == B.hole && A.m_twoJ12 == B.m_twoJ12;
}

bool operator!=(const CSF2p1h &A, const CSF2p1h &B) { return !(A == B); }

int CSF2p1h::num_different(const CSF2p1h &A, const CSF2p1h &B) {
  // Particles of A not in B (counting repeats; states are sorted), plus the
  // hole of B if A does not have it
  const auto common_particles =
    A.states == B.states ? 2 :
    (A.state(0) == B.state(0) || A.state(0) == B.state(1) ||
     A.state(1) == B.state(0) || A.state(1) == B.state(1)) ?
                           1 :
                           0;
  return 2 - common_particles + (A.hole == B.hole ? 0 : 1);
}

std::string CSF2p1h::config(bool relativistic) const {
  const auto [na, ka] = Angular::index_to_nk(states[0]);
  const auto [nb, kb] = Angular::index_to_nk(states[1]);
  const auto [nh, kh] = Angular::index_to_nk(hole);
  auto s1 = DiracSpinor::shortSymbol(na, ka);
  auto s2 = DiracSpinor::shortSymbol(nb, kb);
  auto sh = DiracSpinor::shortSymbol(nh, kh);
  if (!relativistic) {
    s1.pop_back();
    s2.pop_back();
    sh.pop_back();
  }
  return (s1 == s2 ? s1 + "^2" : s1 + s2) + " " + sh + "^-1";
}

//==============================================================================
std::vector<CSF> form_CSF2p1hs(int twoJ, int parity,
                               const std::vector<DiracSpinor> &particle_basis,
                               const std::vector<DiracSpinor> &holes,
                               const ConfigurationSpace &config_space) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ > 0 && twoJ % 2 == 1 &&
         "2*J must be odd for two particles and one hole");

  std::vector<CSF> CSFs;
  for (const auto &a : holes) {
    for (const auto &v : particle_basis) {
      for (const auto &w : particle_basis) {

        // Symmetric: only include unique pairs once
        if (w < v)
          continue;

        if (v.parity() * w.parity() * a.parity() != parity)
          continue;

        if (!config_space.included({v.nk_index(), w.nk_index()},
                                   {a.nk_index()}))
          continue;

        // One CSF per pair angular momentum J12 that couples with the hole
        // to J; identical particles can only give even J12
        for (int twoJ12 = std::abs(v.twoj() - w.twoj());
             twoJ12 <= v.twoj() + w.twoj(); twoJ12 += 2) {
          if (v == w && twoJ12 % 4 != 0)
            continue;
          if (Angular::triangle(twoJ12, a.twoj(), twoJ) == 0)
            continue;
          CSFs.emplace_back(
            CSF2p1h{v.nk_index(), w.nk_index(), a.nk_index(), twoJ12});
        }
      }
    }
  }
  return CSFs;
}

//==============================================================================
std::vector<CSF> form_CSFs(int twoJ, int parity,
                           const std::vector<DiracSpinor> &cisp_basis,
                           int num_electrons,
                           const ConfigurationSpace &config_space,
                           const std::vector<DiracSpinor> &holes,
                           const std::vector<DiracSpinor> &holes_basis) {
  std::vector<CSF> CSFs;
  if (num_electrons == 1) {
    CSFs = form_CSF1s(twoJ, parity, cisp_basis, config_space);
  } else if (num_electrons == 2) {
    CSFs = form_CSF2s(twoJ, parity, cisp_basis, config_space);
  } else if (num_electrons == 3) {
    CSFs = form_CSF3s(twoJ, parity, cisp_basis, config_space);
  } else {
    assert(false &&
           "CI is implemented for one, two, or three valence electrons");
  }

  // Hole CSFs: two particles and one hole, with one valence electron
  if (!holes.empty()) {
    assert(num_electrons == 1 &&
           "Hole CSFs are implemented for one valence electron only");
    const auto &particle_basis = holes_basis.empty() ? cisp_basis : holes_basis;
    for (auto &csf :
         form_CSF2p1hs(twoJ, parity, particle_basis, holes, config_space)) {
      CSFs.emplace_back(std::move(csf));
    }
  }
  return CSFs;
}

//==============================================================================
int CSF::parity() const {
  return std::visit([](const auto &csf) { return csf.parity(); }, m_csf);
}

std::size_t CSF::num_electrons() const {
  return std::visit([](const auto &csf) { return csf.num_electrons(); }, m_csf);
}

std::size_t CSF::num_particles() const {
  return std::visit([](const auto &csf) { return csf.num_particles(); }, m_csf);
}

std::size_t CSF::num_holes() const {
  return std::visit([](const auto &csf) { return csf.num_holes(); }, m_csf);
}

std::vector<DiracSpinor::Index> CSF::holes() const {
  return std::visit(
    [](const auto &csf) -> std::vector<DiracSpinor::Index> {
      if constexpr (std::is_same_v<std::decay_t<decltype(csf)>, CSF2p1h>) {
        return {csf.hole};
      } else {
        return {};
      }
    },
    m_csf);
}

std::string CSF::config(bool relativistic) const {
  return std::visit(
    [relativistic](const auto &csf) { return csf.config(relativistic); },
    m_csf);
}

std::vector<DiracSpinor::Index> CSF::orbitals() const {
  return std::visit(
    [](const auto &csf) {
      return std::vector<DiracSpinor::Index>(csf.states.cbegin(),
                                             csf.states.cend());
    },
    m_csf);
}

DiracSpinor::Index CSF::state(std::size_t i) const {
  return std::visit([i](const auto &csf) { return csf.state(i); }, m_csf);
}

bool operator==(const CSF &A, const CSF &B) {
  return std::visit(
    [](const auto &a, const auto &b) {
      if constexpr (std::is_same_v<decltype(a), decltype(b)>) {
        return a == b;
      } else {
        return false;
      }
    },
    A.m_csf, B.m_csf);
}

bool operator!=(const CSF &A, const CSF &B) { return !(A == B); }

//------------------------------------------------------------------------------
int num_different(const CSF &A, const CSF &B) {
  // One electron x, and two particles (v, w) with a hole: the determinants
  // differ by one orbital (the hole replaced by the other particle) if x is
  // one of the particles, else by two
  const auto one_electron_vs_hole = [](const CSF1 &x, const CSF2p1h &vw) {
    return (x.state(0) == vw.state(0) || x.state(0) == vw.state(1)) ? 1 : 2;
  };
  return std::visit(
    [&one_electron_vs_hole](const auto &a, const auto &b) -> int {
      using A_type = std::decay_t<decltype(a)>;
      using B_type = std::decay_t<decltype(b)>;
      if constexpr (std::is_same_v<A_type, B_type>) {
        return A_type::num_different(a, b);
      } else if constexpr (std::is_same_v<A_type, CSF1> &&
                           std::is_same_v<B_type, CSF2p1h>) {
        return one_electron_vs_hole(a, b);
      } else if constexpr (std::is_same_v<A_type, CSF2p1h> &&
                           std::is_same_v<B_type, CSF1>) {
        return one_electron_vs_hole(b, a);
      } else {
        assert(false && "num_different between CSFs of different kinds is "
                        "not implemented");
        return 99;
      }
    },
    A.variant(), B.variant());
}

//==============================================================================
std::size_t num_interacting_pairs(const std::vector<CSF> &CSFs) {
  // Configurations (particle orbitals, and holes): index of one of its CSFs,
  // and number of CSFs. Two configurations interact for all of their CSFs,
  // or for none
  using Configuration =
    std::pair<std::vector<DiracSpinor::Index>, std::vector<DiracSpinor::Index>>;
  std::map<Configuration, std::pair<std::size_t, std::size_t>> config_map;
  for (std::size_t i = 0; i < CSFs.size(); ++i) {
    auto &[first_CSF, num_CSFs] =
      config_map
        .try_emplace(Configuration{CSFs[i].orbitals(), CSFs[i].holes()}, i, 0)
        .first->second;
    ++num_CSFs;
  }
  std::vector<std::pair<std::size_t, std::size_t>> configs;
  for (const auto &config : config_map) {
    configs.push_back(config.second);
  }

  std::size_t count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : count)
  for (std::size_t a = 0; a < configs.size(); ++a) {
    const auto &[A_index, num_A] = configs[a];
    // Pairs within the configuration, including each CSF with itself
    count += num_A * (num_A + 1) / 2;
    for (std::size_t b = 0; b < a; ++b) {
      const auto &[B_index, num_B] = configs[b];
      if (num_different(CSFs[A_index], CSFs[B_index]) <= 2) {
        count += num_A * num_B;
      }
    }
  }
  return count;
}

//==============================================================================
std::vector<OneBodyConnected>
one_body_connected(const std::vector<CSF> &CSFAs,
                   const std::vector<CSF> &CSFBs) {
  // A configuration: its particle orbitals (sorted) and its holes (sorted)
  using Orbitals = std::vector<DiracSpinor::Index>;
  using Config = std::pair<Orbitals, Orbitals>;

  // Sorted orbitals with the one at 'position' removed (stays sorted)
  const auto remove_orbital = [](const Orbitals &orbitals,
                                 std::size_t position) {
    Orbitals sub_orbitals;
    sub_orbitals.reserve(orbitals.size() - 1);
    for (std::size_t i = 0; i < orbitals.size(); ++i) {
      if (i != position) {
        sub_orbitals.push_back(orbitals[i]);
      }
    }
    return sub_orbitals;
  };

  // The configuration with one orbital (particle or hole) removed, in every
  // way: two configurations of the same size differing by one orbital share
  // one of these
  const auto remove_one = [&remove_orbital](const Config &config) {
    const auto &[particles, holes] = config;
    std::vector<Config> sub_configs;
    for (std::size_t position = 0; position < particles.size(); ++position) {
      sub_configs.push_back({remove_orbital(particles, position), holes});
    }
    for (std::size_t position = 0; position < holes.size(); ++position) {
      sub_configs.push_back({particles, remove_orbital(holes, position)});
    }
    return sub_configs;
  };

  // The configuration with one particle and one hole removed, in every way:
  // a configuration is connected to those it is a hole excitation of
  const auto remove_particle_hole = [&remove_orbital](const Config &config) {
    const auto &[particles, holes] = config;
    std::vector<Config> sub_configs;
    for (std::size_t ip = 0; ip < particles.size(); ++ip) {
      for (std::size_t ih = 0; ih < holes.size(); ++ih) {
        sub_configs.push_back(
          {remove_orbital(particles, ip), remove_orbital(holes, ih)});
      }
    }
    return sub_configs;
  };

  // Configurations of B (CSF indices); then, for each configuration with one
  // orbital removed, and with one particle and one hole removed, the B
  // configurations it came from; and the B configurations as they are
  std::map<Config, std::vector<std::size_t>> B_config_map;
  for (std::size_t j = 0; j < CSFBs.size(); ++j) {
    B_config_map[{CSFBs[j].orbitals(), CSFBs[j].holes()}].push_back(j);
  }
  std::vector<std::vector<std::size_t>> B_configs;
  std::map<Config, std::vector<std::size_t>> B_containing;
  std::map<Config, std::vector<std::size_t>> B_excitation_of;
  std::map<Config, std::size_t> B_exact;
  for (const auto &[config, CSF_indices] : B_config_map) {
    for (const auto &sub_config : remove_one(config)) {
      B_containing[sub_config].push_back(B_configs.size());
    }
    for (const auto &sub_config : remove_particle_hole(config)) {
      B_excitation_of[sub_config].push_back(B_configs.size());
    }
    B_exact[config] = B_configs.size();
    B_configs.push_back(CSF_indices);
  }

  std::map<Config, std::vector<std::size_t>> A_config_map;
  for (std::size_t i = 0; i < CSFAs.size(); ++i) {
    A_config_map[{CSFAs[i].orbitals(), CSFAs[i].holes()}].push_back(i);
  }
  std::vector<std::pair<Config, std::vector<std::size_t>>> A_configs(
    std::make_move_iterator(A_config_map.begin()),
    std::make_move_iterator(A_config_map.end()));
  std::vector<OneBodyConnected> connections(A_configs.size());
#pragma omp parallel for schedule(dynamic, 64)
  for (std::size_t a = 0; a < A_configs.size(); ++a) {
    const auto &[config, CSF_indices] = A_configs[a];
    // B configurations of the same size sharing all but one orbital
    // (repeated orbitals, or identical configurations, give duplicates);
    // those this one is a hole excitation of; and the hole excitations of it
    std::vector<std::size_t> partners;
    for (const auto &sub_config : remove_one(config)) {
      const auto it = B_containing.find(sub_config);
      if (it != B_containing.end()) {
        partners.insert(partners.end(), it->second.begin(), it->second.end());
      }
    }
    for (const auto &sub_config : remove_particle_hole(config)) {
      const auto it = B_exact.find(sub_config);
      if (it != B_exact.end()) {
        partners.push_back(it->second);
      }
    }
    {
      const auto it = B_excitation_of.find(config);
      if (it != B_excitation_of.end()) {
        partners.insert(partners.end(), it->second.begin(), it->second.end());
      }
    }
    std::sort(partners.begin(), partners.end());
    partners.erase(std::unique(partners.begin(), partners.end()),
                   partners.end());

    std::vector<std::size_t> connected_CSFs;
    for (const auto b : partners) {
      connected_CSFs.insert(connected_CSFs.end(), B_configs[b].begin(),
                            B_configs[b].end());
    }
    connections[a] = {CSF_indices, std::move(connected_CSFs)};
  }
  return connections;
}

//==============================================================================
double particle_hole_sum(const CSF2p1h &A, const CSF2p1h &B, int twoJ,
                         const PairElement &ph_element) {
  const auto [x, y] = A.states;
  const auto [v, w] = B.states;
  const auto b = A.hole;
  const auto a = B.hole;
  const auto twoJ12_A = A.twoJ12();
  const auto twoJ12_B = B.twoJ12();

  const auto tjx = Angular::nkindex_to_twoj(x);
  const auto tjy = Angular::nkindex_to_twoj(y);
  const auto tjv = Angular::nkindex_to_twoj(v);
  const auto tjw = Angular::nkindex_to_twoj(w);
  const auto tja = Angular::nkindex_to_twoj(a);
  const auto tjb = Angular::nkindex_to_twoj(b);

  // Recoupling of the pair to a particle-hole pair, with the first particle
  // as spectator: |(vw)J12, a~; J> = sum_K R_K(vwa) |v, (w a~)K; J>
  const auto R = [twoJ](int tj1, int tj2, int twoJ12, int tjh, int twoK) {
    return Angular::neg1pow_2(tj1 + tj2 + tjh + twoJ) *
           std::sqrt(double((twoJ12 + 1) * (twoK + 1))) *
           Angular::sixj_2(tj1, tj2, twoJ12, tjh, twoJ, twoK);
  };
  // Phase for exchanging the two particles of a pair
  const auto s_vw = double(Angular::neg1pow_2(tjv + tjw - twoJ12_B + 2));
  const auto s_xy = double(Angular::neg1pow_2(tjx + tjy - twoJ12_A + 2));
  const auto eta_vw = v == w ? 1.0 / std::sqrt(2.0) : 1.0;
  const auto eta_xy = x == y ? 1.0 / std::sqrt(2.0) : 1.0;

  // Spectator (first particle of bra and ket), then the particles of the
  // interacting particle-hole pairs; with the exchange phases
  struct Term {
    DiracSpinor::Index spectator_A, spectator_B;
    DiracSpinor::Index particle_A, particle_B;
    int tj_spectator_A, tj_spectator_B, tj_particle_A, tj_particle_B;
    double phase;
  };
  const std::array<Term, 4> terms{
    Term{x, v, y, w, tjx, tjv, tjy, tjw, 1.0},
    Term{x, w, y, v, tjx, tjw, tjy, tjv, s_vw},
    Term{y, v, x, w, tjy, tjv, tjx, tjw, s_xy},
    Term{y, w, x, v, tjy, tjw, tjx, tjv, s_vw * s_xy}};

  double sum = 0.0;
  for (const auto &t : terms) {
    if (t.spectator_A != t.spectator_B)
      continue;
    // K: triangle with the spectator and J, and with each particle-hole pair
    const auto twoK_min = std::max({std::abs(t.tj_spectator_B - twoJ),
                                    std::abs(t.tj_particle_B - tja),
                                    std::abs(t.tj_particle_A - tjb)});
    const auto twoK_max = std::min(
      {t.tj_spectator_B + twoJ, t.tj_particle_B + tja, t.tj_particle_A + tjb});
    for (int twoK = twoK_min; twoK <= twoK_max; twoK += 2) {
      const auto R_B =
        R(t.tj_spectator_B, t.tj_particle_B, twoJ12_B, tja, twoK);
      const auto R_A =
        R(t.tj_spectator_A, t.tj_particle_A, twoJ12_A, tjb, twoK);
      if (R_A == 0.0 || R_B == 0.0)
        continue;
      sum += t.phase * R_A * R_B *
             ph_element(t.particle_A, b, t.particle_B, a, twoK);
    }
  }
  return eta_vw * eta_xy * sum;
}

//==============================================================================
double LS_amplitude(int n1, int l1, int twoj1, int n2, int l2, int twoj2, int L,
                    int S, int twoJ) {

  const bool same_nl = (n1 == n2 && l1 == l2);

  // Same non-rel shell: only L+S even terms exist (Pauli)
  if (same_nl && (L + S) % 2 != 0)
    return 0.0;

  // Same shell, different j: L+S odd components cancel in the
  // antisymmetrisation; the surviving even part carries sqrt(2)
  const double eta = (same_nl && twoj1 != twoj2) ? std::sqrt(2.0) : 1.0;

  const auto factor =
    std::sqrt(double((twoj1 + 1) * (twoj2 + 1) * (2 * L + 1) * (2 * S + 1)));

  return eta * factor *
         Angular::ninej_2(2 * l1, 2 * l2, 2 * L, 1, 1, 2 * S, twoj1, twoj2,
                          twoJ);
}

//==============================================================================
std::pair<double, double> L2S2_pair(DiracSpinor::Index x, DiracSpinor::Index y,
                                    DiracSpinor::Index xp,
                                    DiracSpinor::Index yp, int twoJ) {
  const auto [n1, k1] = Angular::index_to_nk(x);
  const auto [n2, k2] = Angular::index_to_nk(y);
  const auto [n1p, k1p] = Angular::index_to_nk(xp);
  const auto [n2p, k2p] = Angular::index_to_nk(yp);
  const auto l1 = Angular::l_k(k1);
  const auto l2 = Angular::l_k(k2);
  // Same non-rel configuration, in the same (stored) order
  if (n1 != n1p || n2 != n2p || l1 != Angular::l_k(k1p) ||
      l2 != Angular::l_k(k2p)) {
    return {0.0, 0.0};
  }
  const auto twoj1 = Angular::twoj_k(k1);
  const auto twoj2 = Angular::twoj_k(k2);
  const auto twoj1p = Angular::twoj_k(k1p);
  const auto twoj2p = Angular::twoj_k(k2p);

  double L2 = 0.0;
  double S2 = 0.0;
  for (int S = 0; S <= 1; ++S) {
    for (int L = std::abs(l1 - l2); L <= l1 + l2; ++L) {
      if (twoJ > 2 * (L + S) || twoJ < 2 * std::abs(L - S))
        continue;
      const auto A = LS_amplitude(n1, l1, twoj1, n2, l2, twoj2, L, S, twoJ);
      const auto Ap = LS_amplitude(n1, l1, twoj1p, n2, l2, twoj2p, L, S, twoJ);
      L2 += A * Ap * L * (L + 1.0);
      S2 += A * Ap * S * (S + 1.0);
    }
  }
  return {L2, S2};
}

//==============================================================================
std::pair<double, double> L2S2_particle_hole(DiracSpinor::Index y,
                                             DiracSpinor::Index b,
                                             DiracSpinor::Index w,
                                             DiracSpinor::Index a, int twoK) {
  // (l_1 + l_2)^2 between the pairs (w a) and (y b), in product coupling:
  // a particle and a hole are never in the same shell
  auto [L2, S2] = L2S2_pair(w, a, y, b, twoK);
  // Less the one-body terms, which are diagonal
  if (y == w && a == b) {
    const auto lw = Angular::nkindex_to_l(w);
    const auto la = Angular::nkindex_to_l(a);
    L2 -= lw * (lw + 1.0) + la * (la + 1.0);
    S2 -= 1.5;
  }
  return {L2, S2};
}

//==============================================================================
std::pair<double, double> L2S2_2p1h(const CSF2p1h &A, const CSF2p1h &B,
                                    int twoJ) {
  const auto [x, y] = A.states;
  const auto [v, w] = B.states;
  const auto same_hole = A.hole == B.hole;
  const auto same_J12 = A.twoJ12() == B.twoJ12();

  double L2 = 0.0;
  double S2 = 0.0;

  // The pair, with the hole as spectator
  if (same_hole && same_J12) {
    const auto [L2_pair, S2_pair] = L2S2_pair(v, w, x, y, B.twoJ12());
    L2 += L2_pair;
    S2 += S2_pair;
  }

  // The hole, with the pair as spectator
  if (same_hole && same_J12 && A.states == B.states) {
    const auto la = Angular::nkindex_to_l(A.hole);
    L2 += la * (la + 1.0);
    S2 += 0.75;
  }

  // Cross terms between each particle and the hole
  const auto L_cross = [](DiracSpinor::Index yp, DiracSpinor::Index bp,
                          DiracSpinor::Index wp, DiracSpinor::Index ap,
                          int twoK) {
    return L2S2_particle_hole(yp, bp, wp, ap, twoK).first;
  };
  const auto S_cross = [](DiracSpinor::Index yp, DiracSpinor::Index bp,
                          DiracSpinor::Index wp, DiracSpinor::Index ap,
                          int twoK) {
    return L2S2_particle_hole(yp, bp, wp, ap, twoK).second;
  };
  L2 += particle_hole_sum(A, B, twoJ, L_cross);
  S2 += particle_hole_sum(A, B, twoJ, S_cross);

  return {L2, S2};
}

//==============================================================================
std::pair<double, double> L2S2_element(const CSF &A, const CSF &B, int twoJ) {
  return std::visit(
    [twoJ](const auto &a, const auto &b) -> std::pair<double, double> {
      using CSFType = std::decay_t<decltype(a)>;
      using B_type = std::decay_t<decltype(b)>;
      if constexpr ((std::is_same_v<CSFType, CSF1> &&
                     std::is_same_v<B_type, CSF2p1h>) ||
                    (std::is_same_v<CSFType, CSF2p1h> &&
                     std::is_same_v<B_type, CSF1>)) {
        // L and S conserve the numbers of particles and holes: no element
        // between a one-electron and a hole CSF
        return {0.0, 0.0};
      } else if constexpr (!std::is_same_v<CSFType, B_type>) {
        assert(false && "L^2 and S^2 between CSFs of different kinds are not "
                        "implemented");
        return {0.0, 0.0};
      } else if constexpr (std::is_same_v<CSFType, CSF1>) {
        // Each orbital has definite l, and s = 1/2: both operators are diagonal
        if (a != b) {
          return {0.0, 0.0};
        }
        const auto l = Angular::nkindex_to_l(a.state(0));
        return {l * (l + 1.0), 0.75};
      } else if constexpr (std::is_same_v<CSFType, CSF2>) {
        return L2S2_pair(b.state(0), b.state(1), a.state(0), a.state(1), twoJ);
      } else if constexpr (std::is_same_v<CSFType, CSF2p1h>) {
        return L2S2_2p1h(a, b, twoJ);
      } else {
        // One-body terms: diagonal in the CSFs
        double L2 = 0.0;
        double S2 = 0.0;
        if (a == b) {
          for (const auto orbital : a.states) {
            const auto l = Angular::nkindex_to_l(orbital);
            L2 -= l * (l + 1.0);
          }
          S2 -= 9.0 / 4.0;
        }
        // Pair terms, via parents with the same third orbital and pair J2
        for (const auto &p : a.parents()) {
          for (const auto &q : b.parents()) {
            if (p.spectator != q.spectator || p.twoJ2 != q.twoJ2)
              continue;
            const auto [L2_pq, S2_pq] =
              L2S2_pair(q.pair[0], q.pair[1], p.pair[0], p.pair[1], p.twoJ2);
            L2 += 3.0 * p.coef * q.coef * L2_pq;
            S2 += 3.0 * p.coef * q.coef * S2_pq;
          }
        }
        return {L2, S2};
      }
    },
    A.variant(), B.variant());
}

//==============================================================================
std::pair<double, double>
expectation_L2S2(const LinAlg::View<const double> &coefs,
                 const std::vector<CSF> &csfs, int twoJ) {
  assert(coefs.size() == csfs.size());

  // L^2 and S^2 connect only CSFs of the same non-rel configuration: group
  // the CSFs by configuration, and sum over the pairs within each group
  std::map<std::string, std::vector<std::size_t>> groups;
  for (std::size_t i = 0; i < csfs.size(); ++i) {
    if (coefs[i] != 0.0) {
      groups[csfs[i].config(false)].push_back(i);
    }
  }

  double L2 = 0.0;
  double S2 = 0.0;
  for (const auto &[config, members] : groups) {
    for (const auto a : members) {
      for (const auto b : members) {
        const auto [L2_ab, S2_ab] = L2S2_element(csfs[a], csfs[b], twoJ);
        L2 += coefs[a] * coefs[b] * L2_ab;
        S2 += coefs[a] * coefs[b] * S2_ab;
      }
    }
  }
  return {L2, S2};
}

//==============================================================================
PsiJPi::PsiJPi(int twoJ, int pi, const std::vector<DiracSpinor> &cisp_basis,
               int num_electrons, const ConfigurationSpace &config_space,
               const std::vector<DiracSpinor> &holes,
               const std::vector<DiracSpinor> &holes_basis)
  : m_twoj(twoJ),
    m_pi(pi),
    m_num_electrons(num_electrons),
    m_CSFs(form_CSFs(twoJ, pi, cisp_basis, num_electrons, config_space, holes,
                     holes_basis)) {}

int PsiJPi::num_electrons() const { return m_num_electrons; }

bool PsiJPi::has_holes() const {
  return std::any_of(m_CSFs.cbegin(), m_CSFs.cend(),
                     [](const auto &csf) { return csf.num_holes() > 0; });
}

std::size_t PsiJPi::num_CSFs() const { return m_CSFs.size(); }

const CSF &PsiJPi::csf(std::size_t i) const { return m_CSFs.at(i); }

const std::vector<CSF> &PsiJPi::CSFs() const { return m_CSFs; }

//==============================================================================
// Solves the CI equation for Hamiltonian matrix Hci;
// finds first num_solutions solutions.
// Doesn't set the Config info: have to call update_config_info() manually.
void PsiJPi::solve(const LinAlg::Matrix<double> &Hci, int num_solutions,
                   std::optional<double> all_below_cm) {
  assert(Hci.rows() == Hci.cols());
  assert(Hci.rows() == num_CSFs());

  if (all_below_cm) {
    const auto all_below_au = *all_below_cm / PhysConst::Hartree_invcm;
    const auto [t_num, t_evals, t_evecs] =
      LinAlg::symmhEigensystem(Hci, all_below_au);
    m_num_solutions = std::size_t(t_num);
    m_Solution.first = std::move(t_evals);
    m_Solution.second = std::move(t_evecs);
  } else if (num_solutions <= 0 || num_solutions >= int(num_CSFs())) {
    m_Solution = LinAlg::symmhEigensystem(Hci);
    m_num_solutions = num_CSFs();
  } else {
    m_Solution = LinAlg::symmhEigensystem(Hci, num_solutions);
    m_num_solutions = std::size_t(num_solutions);
  }

  if (m_num_solutions >= m_Info.size()) {
    m_Info.resize(m_num_solutions);
  }
}

//==============================================================================
// Finds the lowest num_solutions solutions (or all below all_below_cm)
// iteratively (Davidson), for a sparse Hci. Doesn't set the Config info (as
// solve(), above)
LinAlg::DavidsonConvergence
PsiJPi::solve(const LinAlg::SparseSymmetricMatrix<double> &Hci,
              int num_solutions, std::optional<double> all_below_cm,
              std::ostream &outstream) {
  assert(Hci.rows() == num_CSFs());
  const auto N = int(num_CSFs());

  // A single CSF is its own solution (Davidson needs fewer solutions than
  // CSFs)
  if (N == 1) {
    set_solution(Hci.diagonal()[0], LinAlg::Vector<double>{1.0});
    return {0, 0.0, true};
  }

  // Number of solutions to find in each Davidson solve: as requested, or,
  // for all below an energy, a few to start with, doubled until enough.
  // Davidson finds at most N - 1: if more are requested, fewer are returned
  constexpr int initial_below = 8;
  auto number = all_below_cm       ? initial_below :
                num_solutions <= 0 ? N :
                                     num_solutions;
  number = std::min(number, N - 1);
  const auto all_below_au =
    all_below_cm ? *all_below_cm / PhysConst::Hartree_invcm : 0.0;
  // Beyond this many solutions, Davidson is slower, and takes more memory,
  // than the dense solver
  const auto many_solutions = std::max(1, N / 20);
  bool warned = false;

  while (true) {
    if (number > many_solutions && !warned) {
      fmt::print(outstream,
                 "Warning: Davidson for {}/{} solutions: dense may be better\n",
                 number, N);
      warned = true;
    }
    auto [values, vectors, convergence] = LinAlg::Davidson(Hci, number);

    // All below: enough once the highest found is above the threshold (or
    // no more can be found)
    if (all_below_cm && convergence.converged && number < N - 1 &&
        values[std::size_t(number) - 1] <= all_below_au) {
      number = std::min(2 * number, N - 1);
      continue;
    }

    m_Solution.first = std::move(values);
    m_Solution.second = std::move(vectors);
    m_num_solutions = std::size_t(number);
    if (all_below_cm) {
      m_num_solutions = 0;
      while (m_num_solutions < std::size_t(number) &&
             m_Solution.first[m_num_solutions] <= all_below_au) {
        ++m_num_solutions;
      }
    }
    if (m_num_solutions >= m_Info.size()) {
      m_Info.resize(m_num_solutions);
    }
    return convergence;
  }
}

//==============================================================================
// Stores a single solution, that was not found by diagonalising Hci
void PsiJPi::set_solution(double energy, const LinAlg::Vector<double> &coefs) {
  assert(coefs.size() == num_CSFs());

  m_num_solutions = 1;
  m_Solution.first = LinAlg::Vector<double>{energy};
  m_Solution.second = LinAlg::Matrix<double>(1, coefs.size());
  for (std::size_t i = 0; i < coefs.size(); ++i) {
    m_Solution.second(0, i) = coefs[i];
  }
  m_Info.clear();
  m_Info.resize(m_num_solutions);
}

void PsiJPi::replace_solution(std::size_t i, double energy,
                              const LinAlg::View<const double> &coefs) {
  assert(i < m_num_solutions);
  assert(coefs.size() == num_CSFs());
  m_Solution.first[i] = energy;
  for (std::size_t j = 0; j < coefs.size(); ++j) {
    m_Solution.second(i, j) = coefs[j];
  }
}

// You must manually update the config. info for each solution (if required)
void PsiJPi::update_config_info(std::size_t i, const ConfigInfo &info) {
  assert(m_Info.size() == m_num_solutions);
  m_Info.at(i) = info;
}

// Energy of the ith CI solution
double PsiJPi::energy(std::size_t i) const {
  assert(i < m_num_solutions);
  return m_Solution.first.at(i);
}

// List of CI expansion coefs for the ith CI solution
LinAlg::View<const double> PsiJPi::coefs(std::size_t i) const {
  assert(i < m_num_solutions);
  return m_Solution.second.row_view(i);
}

// The CI coeficient for the ith CI solution, corresponding to the jth CSF
double PsiJPi::coef(std::size_t i, std::size_t j) const {
  assert(i < m_num_solutions);
  return m_Solution.second.at(i, j);
}

// Parity for the CI solutions (+/-1)
int PsiJPi::parity() const { return m_pi; }

// 2J for the CI solutions
int PsiJPi::twoJ() const { return m_twoj; }

// Number of CI solutions stored
std::size_t PsiJPi::num_solutions() const { return m_num_solutions; }

// Configuration info for the ith CI solution, if it has been set
const ConfigInfo &PsiJPi::info(std::size_t i) const {
  assert(m_Info.size() == m_num_solutions);
  return m_Info.at(i);
}

//==============================================================================
bool PsiJPi::read_write(const std::string &fname, IO::FRW::RoW rw,
                        std::ostream &outstream) {
  using IO::FRW::rw_binary;
  const bool readQ = (rw == IO::FRW::read);

  const std::size_t nc = num_CSFs();

  // Marks the file layout: a file without it is never read as this one
  constexpr std::array<char, 8> format_tag{'a', 'm', 'p', 's',
                                           'c', 'i', '0', '2'};
  // Each block: twoJ, pi, num_csfs, num_solutions; then the payload
  constexpr std::size_t header_bytes =
    2 * sizeof(int) + 2 * sizeof(std::size_t);
  // Payload: energies E[num_solutions], coefficients M[num_solutions x num_csfs]
  const auto payload_bytes = [](std::size_t num_csfs, std::size_t num_stored) {
    return (num_stored + num_stored * num_csfs) * sizeof(double);
  };

  // Opens fname for reading, and checks the format tag
  const auto open_tagged = [&fname, &format_tag](std::fstream &f) {
    if (!IO::FRW::file_exists(fname))
      return false;
    IO::FRW::open_binary(f, fname, IO::FRW::read);
    auto tag = format_tag;
    tag.fill('\0');
    rw_binary(f, IO::FRW::read, tag);
    return bool(f) && tag == format_tag;
  };

  // ===== READ =====
  if (readQ) {
    std::fstream f;
    if (!open_tagged(f)) {
      if (IO::FRW::file_exists(fname)) {
        outstream << "Cannot read from " << fname
                  << ": not a CI solutions file of the current format\n";
      }
      return false;
    }
    while (true) {
      int twoj = 0, pi = 0;
      std::size_t num_csfs = 0, num_stored = 0;
      rw_binary(f, IO::FRW::read, twoj, pi, num_csfs, num_stored);
      if (!f) {
        return false;
      }
      if (twoj != m_twoj || pi != m_pi) {
        f.seekg(
          static_cast<std::streamoff>(payload_bytes(num_csfs, num_stored)),
          std::ios::cur);
        continue;
      }
      if (num_csfs != nc) {
        // in theory, should never happen: filename includes basis string..
        outstream << "\nCannot read from " << fname << ": CSF count mismatch "
                  << "(file has " << num_csfs << ", expected " << nc << ").\n";
        return false;
      }
      std::vector<double> E(num_stored), M(num_stored * nc);
      f.read(reinterpret_cast<char *>(E.data()),
             static_cast<std::streamsize>(E.size() * sizeof(double)));
      f.read(reinterpret_cast<char *>(M.data()),
             static_cast<std::streamsize>(M.size() * sizeof(double)));
      if (!f) {
        outstream << "\nCannot read from " << fname << ": file is truncated\n";
        return false;
      }
      outstream << "Reading J=" << J_string(m_twoj)
                << " pi=" << (m_pi > 0 ? "+" : "-") << " (" << num_stored
                << " solutions) from " << fname << "\n";
      m_num_solutions = num_stored;
      m_Solution.first = LinAlg::Vector<double>(std::move(E));
      m_Solution.second = LinAlg::Matrix<double>(num_stored, nc, std::move(M));
      m_Info.resize(m_num_solutions);
      return true;
    }
  }

  // ===== WRITE =====
  // A block new to the file is appended. If this J/Pi is in the file already
  // (its size may differ), or the last block is incomplete, the file is
  // re-written instead: the other complete blocks, as they are, then this
  // one. A file without the format tag is replaced
  bool tagged = false;
  bool rewrite = false;
  // End of the last complete block
  std::streamoff blocks_end = 0;
  {
    std::fstream f;
    tagged = open_tagged(f);
    if (tagged) {
      blocks_end = f.tellg();
      f.seekg(0, std::ios::end);
      const std::streamoff file_end = f.tellg();
      f.seekg(blocks_end);
      while (blocks_end < file_end) {
        int twoj = 0, pi = 0;
        std::size_t num_csfs = 0, num_stored = 0;
        rw_binary(f, IO::FRW::read, twoj, pi, num_csfs, num_stored);
        const auto block_end =
          blocks_end + static_cast<std::streamoff>(
                         header_bytes + payload_bytes(num_csfs, num_stored));
        if (!f || block_end > file_end) {
          rewrite = true;
          break;
        }
        if (twoj == m_twoj && pi == m_pi) {
          if (num_csfs != nc) {
            outstream << "\nCannot write to " << fname
                      << ": CSF count mismatch " << "(file has " << num_csfs
                      << ", expected " << nc << ").\n";
            return false;
          }
          rewrite = true;
        }
        f.seekg(block_end);
        blocks_end = block_end;
      }
    }
  }

  // The other complete blocks, if re-writing
  std::vector<char> other_blocks;
  if (tagged && rewrite) {
    std::fstream f;
    open_tagged(f);
    while (f && f.tellg() < blocks_end) {
      const auto block_start = f.tellg();
      int twoj = 0, pi = 0;
      std::size_t num_csfs = 0, num_stored = 0;
      rw_binary(f, IO::FRW::read, twoj, pi, num_csfs, num_stored);
      const auto block_bytes =
        header_bytes + payload_bytes(num_csfs, num_stored);
      if (twoj != m_twoj || pi != m_pi) {
        const auto old_size = other_blocks.size();
        other_blocks.resize(old_size + block_bytes);
        f.seekg(block_start);
        f.read(other_blocks.data() + old_size,
               static_cast<std::streamsize>(block_bytes));
      }
      f.seekg(block_start + static_cast<std::streamoff>(block_bytes));
    }
  }

  outstream << "Writing J=" << J_string(m_twoj)
            << " pi=" << (m_pi > 0 ? "+" : "-") << " (" << m_num_solutions
            << " solutions) to " << fname << "\n";

  // Solutions are rows of the coefficient matrix (which may hold more rows
  // than solutions, e.g., if all eigenvectors were found); the energy vector
  // may hold more entries than solutions
  assert(m_num_solutions == 0 || (m_Solution.second.cols() == nc &&
                                  m_Solution.second.rows() >= m_num_solutions &&
                                  m_Solution.first.size() >= m_num_solutions));

  std::fstream f;
  if (tagged && !rewrite) {
    IO::FRW::open_binary(f, fname, IO::FRW::update);
    f.seekp(0, std::ios::end);
  } else {
    IO::FRW::open_binary(f, fname, IO::FRW::write);
    auto tag = format_tag;
    rw_binary(f, IO::FRW::write, tag);
    f.write(other_blocks.data(),
            static_cast<std::streamsize>(other_blocks.size()));
  }
  auto num_csfs = nc;
  auto num_stored = m_num_solutions;
  rw_binary(f, IO::FRW::write, m_twoj, m_pi, num_csfs, num_stored);
  f.write(reinterpret_cast<const char *>(m_Solution.first.data()),
          static_cast<std::streamsize>(num_stored * sizeof(double)));
  f.write(reinterpret_cast<const char *>(m_Solution.second.data()),
          static_cast<std::streamsize>(num_stored * nc * sizeof(double)));

  return bool(f);
}

} // namespace CI