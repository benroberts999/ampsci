#include "CSF.hpp"
#include "Angular/include.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <tuple>

namespace CI {

//==============================================================================
CSF2::CSF2(const DiracSpinor &a, const DiracSpinor &b)
  : CSF(a.parity() * b.parity()),
    states(a <= b ? std::array{a.nk_index(), b.nk_index()} :
                    std::array{b.nk_index(), a.nk_index()}) {}

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

  const auto J = read_int(first, last);
  if (!J)
    return std::nullopt;

  return Level{2 * *J, parity, index};
}

//==============================================================================
std::string to_string(const Level &level) {
  return std::to_string(level.twoJ / 2) + (level.parity == 1 ? "+" : "-") +
         ":" + std::to_string(level.index);
}

//==============================================================================
// Forms list of all possible (2-particle) CSFs with given J and parity
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ >= 0 && "2*J must not be negative");

  std::vector<CSF2> CSFs;

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

      CSFs.emplace_back(v, w);
    }
  }

  return CSFs;
}

//==============================================================================
CSF3::CSF3(const std::array<DiracSpinor::Index, 3> &orbitals, int twoJ12,
           std::vector<Parent> parents)
  : CSF(Angular::neg1pow(Angular::nkindex_to_l(orbitals[0]) +
                         Angular::nkindex_to_l(orbitals[1]) +
                         Angular::nkindex_to_l(orbitals[2]))),
    states(orbitals),
    m_twoJ12(twoJ12),
    m_parents(std::move(parents)) {
  std::sort(states.begin(), states.end());
}

DiracSpinor::Index CSF3::state(std::size_t i) const { return states.at(i); }

int CSF3::twoJ12() const { return m_twoJ12; }

const std::vector<CSF3::Parent> &CSF3::parents() const { return m_parents; }

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
std::vector<CSF3> form_CSF3s(int twoJ, int parity,
                             const std::vector<DiracSpinor> &cisp_basis) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ > 0 && twoJ % 2 == 1 && "2*J must be odd for three electrons");

  std::vector<CSF3> CSFs;
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

        auto config_CSFs =
          form_CSF3s({u.nk_index(), v.nk_index(), w.nk_index()}, twoJ);
        CSFs.insert(CSFs.end(), std::make_move_iterator(config_CSFs.begin()),
                    std::make_move_iterator(config_CSFs.end()));
      }
    }
  }

  return CSFs;
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
std::pair<double, double>
expectation_L2S2(const LinAlg::View<const double> &coefs,
                 const std::vector<CSF2> &csfs, int twoJ) {
  assert(coefs.size() == csfs.size());

  // L^2 and S^2 are diagonal in the LS-coupled basis, and (in the non-rel
  // limit) only connect CSFs of the same non-rel configuration {n1l1, n2l2}.
  // Accumulate B(L,S) = sum_I c_I * A_I(L,S) per configuration, keyed by
  // {n1, l1, n2, l2, L, S}. Stored (sorted) orbital order is consistent
  // across the CSFs of one configuration, so signs are consistent.
  std::map<std::array<int, 6>, double> B{};

  for (std::size_t i = 0; i < csfs.size(); ++i) {
    const auto ci = coefs[i];
    if (ci == 0.0)
      continue;
    const auto [n1, k1] = Angular::index_to_nk(csfs[i].state(0));
    const auto [n2, k2] = Angular::index_to_nk(csfs[i].state(1));
    const auto l1 = Angular::l_k(k1);
    const auto l2 = Angular::l_k(k2);
    const auto twoj1 = Angular::twoj_k(k1);
    const auto twoj2 = Angular::twoj_k(k2);

    for (int S = 0; S <= 1; ++S) {
      for (int L = std::abs(l1 - l2); L <= l1 + l2; ++L) {
        if (twoJ > 2 * (L + S) || twoJ < 2 * std::abs(L - S))
          continue;
        const auto A = LS_amplitude(n1, l1, twoj1, n2, l2, twoj2, L, S, twoJ);
        if (A != 0.0) {
          B[{n1, l1, n2, l2, L, S}] += ci * A;
        }
      }
    }
  }

  double L2 = 0.0;
  double S2 = 0.0;
  for (const auto &[key, b] : B) {
    const auto L = key[4];
    const auto S = key[5];
    L2 += b * b * L * (L + 1.0);
    S2 += b * b * S * (S + 1.0);
  }

  return {L2, S2};
}

//==============================================================================

// Solves the CI equation for Hamiltonian matrix Hci;
// finds first num_solutions solutions.
// Doesn't set the Config info: have to call update_config_info() manually.
void PsiJPi::solve(const LinAlg::Matrix<double> &Hci, int num_solutions,
                   std::optional<double> all_below_cm) {
  assert(Hci.rows() == Hci.cols());
  assert(Hci.rows() == m_CSFs.size());

  if (all_below_cm) {
    const auto all_below_au = *all_below_cm / PhysConst::Hartree_invcm;
    const auto [t_num, t_evals, t_evecs] =
      LinAlg::symmhEigensystem(Hci, all_below_au);
    m_num_solutions = std::size_t(t_num);
    m_Solution.first = std::move(t_evals);
    m_Solution.second = std::move(t_evecs);
  } else if (num_solutions <= 0 || num_solutions >= int(m_CSFs.size())) {
    m_Solution = LinAlg::symmhEigensystem(Hci);
    m_num_solutions = m_CSFs.size();
  } else {
    m_Solution = LinAlg::symmhEigensystem(Hci, num_solutions);
    m_num_solutions = std::size_t(num_solutions);
  }

  if (m_num_solutions >= m_Info.size()) {
    m_Info.resize(m_num_solutions);
  }
}

// Stores a single solution, that was not found by diagonalising Hci
void PsiJPi::set_solution(double energy, const LinAlg::Vector<double> &coefs) {
  assert(coefs.size() == m_CSFs.size());

  m_num_solutions = 1;
  m_Solution.first = LinAlg::Vector<double>{energy};
  m_Solution.second = LinAlg::Matrix<double>(1, coefs.size());
  for (std::size_t i = 0; i < coefs.size(); ++i) {
    m_Solution.second(0, i) = coefs[i];
  }
  m_Info.clear();
  m_Info.resize(m_num_solutions);
}

// You must manually update the config. info for each solution (if required)
void PsiJPi::update_config_info(std::size_t i, const ConfigInfo &info) {
  assert(m_Info.size() == m_num_solutions);
  m_Info.at(i) = info;
}

// Full list of CSFs
const std::vector<CSF2> &PsiJPi::CSFs() const { return m_CSFs; }

// The ith CSF
const CSF2 &PsiJPi::CSF(std::size_t i) const { return m_CSFs.at(i); }

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

  std::size_t nc = m_CSFs.size();

  // Bytes of one J/Pi block: num_solutions + E[nc] + M[nc*nc]
  const auto block_bytes = [](std::size_t n) {
    return sizeof(std::size_t) + (n + n * n) * sizeof(double);
  };
  // Bytes of one block's prefix: twoJ + pi + num_csfs
  constexpr std::size_t prefix_bytes = 2 * sizeof(int) + sizeof(std::size_t);

  // Write this block's solutions (num_solutions, E, M) to an open stream
  const auto write_jpi_block = [&](std::fstream &f) {
    rw_binary(f, IO::FRW::write, m_num_solutions);
    // Energies: always nc doubles (zero-pad if partial solve)
    std::vector<double> E(nc, 0.0);
    std::copy_n(m_Solution.first.data(), std::min(m_Solution.first.size(), nc),
                E.data());
    f.write(reinterpret_cast<const char *>(E.data()),
            static_cast<std::streamsize>(nc * sizeof(double)));
    // Coefficient matrix: always nc*nc doubles (zero-pad extra rows)
    std::vector<double> M(nc * nc, 0.0);
    const auto nr = std::min(m_Solution.second.rows(), nc);
    const auto nc2 = std::min(m_Solution.second.cols(), nc);
    for (std::size_t i = 0; i < nr; ++i) {
      std::copy_n(m_Solution.second.data() + i * m_Solution.second.cols(), nc2,
                  M.data() + i * nc);
    }
    f.write(reinterpret_cast<const char *>(M.data()),
            static_cast<std::streamsize>(nc * nc * sizeof(double)));
  };

  if (readQ && !IO::FRW::file_exists(fname)) {
    return false;
  }

  // Scan file for matching (twoJ, pi) block; record its offset if found
  std::streamoff match_offset = -1;
  std::size_t match_num_csfs = 0;
  if (IO::FRW::file_exists(fname)) {
    std::fstream f;
    IO::FRW::open_binary(f, fname, IO::FRW::read);
    while (f) {
      const auto pos = f.tellg();
      int twoj = 0, pi = 0;
      std::size_t num_csfs = 0;
      rw_binary(f, IO::FRW::read, twoj, pi, num_csfs);
      if (!f) {
        break;
      }
      if (twoj == m_twoj && pi == m_pi) {
        match_offset = pos;
        match_num_csfs = num_csfs;
        break;
      }
      f.seekg(static_cast<std::streamoff>(block_bytes(num_csfs)),
              std::ios::cur);
    }
  }

  // ===== READ =====
  if (readQ) {
    if (match_offset < 0) {
      return false;
    }
    if (match_num_csfs != nc) {
      // in theory, should never happen: filename includes basis string..
      outstream << "\nCannot read from " << fname << ": CSF count mismatch "
                << "(file has " << match_num_csfs << ", expected " << nc
                << ").\n";
      return false;
    }
    std::fstream f;
    IO::FRW::open_binary(f, fname, IO::FRW::read);
    f.seekg(match_offset + static_cast<std::streamoff>(prefix_bytes));
    rw_binary(f, IO::FRW::read, m_num_solutions);
    const auto Jstr = (m_twoj % 2 == 0) ? std::to_string(m_twoj / 2) :
                                          std::to_string(m_twoj) + "/2";
    outstream << "Reading J=" << Jstr << " pi=" << (m_pi > 0 ? "+" : "-")
              << " (" << m_num_solutions << " solutions) from " << fname
              << "\n";
    std::vector<double> E(nc), M(nc * nc);
    f.read(reinterpret_cast<char *>(E.data()),
           static_cast<std::streamsize>(nc * sizeof(double)));
    f.read(reinterpret_cast<char *>(M.data()),
           static_cast<std::streamsize>(nc * nc * sizeof(double)));
    m_Solution.first = LinAlg::Vector<double>(std::move(E));
    m_Solution.second = LinAlg::Matrix<double>(nc, nc, std::move(M));
    m_Info.resize(m_num_solutions);
    return true;
  }

  // ===== WRITE =====
  {
    const auto Jstr = (m_twoj % 2 == 0) ? std::to_string(m_twoj / 2) :
                                          std::to_string(m_twoj) + "/2";
    outstream << "Writing J=" << Jstr << " pi=" << (m_pi > 0 ? "+" : "-")
              << " (" << m_num_solutions << " solutions) to " << fname << "\n";
  }

  if (match_offset >= 0) {
    // Block exists: overwrite in-place (payload size fixed by num_csfs)
    if (match_num_csfs != nc) {
      outstream << "\nCannot write to " << fname << ": CSF count mismatch "
                << "(file has " << match_num_csfs << ", expected " << nc
                << ").\n";
      return false;
    }
    std::fstream f;
    IO::FRW::open_binary(f, fname, IO::FRW::update);
    f.seekp(match_offset + static_cast<std::streamoff>(prefix_bytes));
    write_jpi_block(f);
  } else {
    // New block: append prefix + payload at end of file
    std::fstream f;
    IO::FRW::open_binary(f, fname, IO::FRW::update);
    if (!f) {
      IO::FRW::open_binary(f, fname, IO::FRW::write);
    }
    f.seekp(0, std::ios::end);
    rw_binary(f, IO::FRW::write, m_twoj, m_pi, nc);
    write_jpi_block(f);
  }

  return true;
}

} // namespace CI