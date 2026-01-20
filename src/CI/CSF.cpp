#include "CSF.hpp"
#include "Angular/include.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <algorithm>
#include <array>
#include <cassert>

namespace CI {

//==============================================================================
CSF2::CSF2(const DiracSpinor &a, const DiracSpinor &b)
    : m_parity(a.parity() * b.parity()),
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

int CSF2::parity() const { return m_parity; }

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

// Solves the CI equation for Hamiltonian matrix Hci;
// finds first num_solutions solutions.
// Doesn't set the Config info: have to call update_config_info() manually.
void PsiJPi::solve(const LinAlg::Matrix<double> &Hci, int num_solutions,
                   std::optional<double> all_below) {
  assert(Hci.rows() == Hci.cols());
  assert(Hci.rows() == m_CSFs.size());

  if (all_below) {
    const auto [t_num, t_evals, t_evecs] =
        LinAlg::symmhEigensystem(Hci, *all_below / PhysConst::Hartree_invcm);
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

} // namespace CI