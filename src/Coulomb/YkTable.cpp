#include "YkTable.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

namespace Coulomb {

//******************************************************************************
YkTable::YkTable(const Grid *const in_grid,
                 const std::vector<DiracSpinor> *const in_a_orbs,
                 const std::vector<DiracSpinor> *const in_b_orbs)
    : m_a_orbs(in_a_orbs),
      m_b_orbs(in_b_orbs == nullptr ? in_a_orbs : in_b_orbs),
      m_grid(in_grid),
      m_aisb([&]() {
        return (in_b_orbs == nullptr || in_a_orbs == in_b_orbs);
      }()) {
  if (!m_a_orbs->empty() && !m_b_orbs->empty())
    update_y_ints();
}

//******************************************************************************
int YkTable::max_tj() const {
  if (m_a_orbs->empty()) {
    return 0;
  }
  const auto maxtj_a = std::max_element(m_a_orbs->cbegin(), m_a_orbs->cend(),
                                        DiracSpinor::comp_j);
  if (m_aisb || m_b_orbs->empty()) {
    return maxtj_a->twoj();
  }
  const auto maxtj_b = std::max_element(m_b_orbs->cbegin(), m_b_orbs->cend(),
                                        DiracSpinor::comp_j);
  return std::max(maxtj_a->twoj(), maxtj_b->twoj());
}

//******************************************************************************
std::pair<int, int> YkTable::k_minmax(const DiracSpinor &a,
                                      const DiracSpinor &b) {
  return std::make_pair(std::abs(a.twoj() - b.twoj()) / 2,
                        (a.twoj() + b.twoj()) / 2);
}

//******************************************************************************
const std::vector<double> &YkTable::get_yk_ab(const int k,
                                              const DiracSpinor &Fa,
                                              const DiracSpinor &Fb) const {
  const auto &[kmin, kmax] = k_minmax(Fa, Fb);
  const auto ik = std::size_t(k - kmin);
  const auto &yab = get_y_ab(Fa, Fb);
  assert(k >= kmin);
  assert(k <= kmax);
  assert(ik < yab.size());
  return yab[ik];
}
//------------------------------------------------------------------------------
const std::vector<std::vector<double>> &
YkTable::get_y_ab(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  const auto a = std::find(m_a_orbs->begin(), m_a_orbs->end(), Fa);
  const auto b =
      Fa == Fb ? a : std::find(m_b_orbs->begin(), m_b_orbs->end(), Fb);
  assert(a != m_a_orbs->end());
  assert(b != m_b_orbs->end());
  const auto ia = std::size_t(a - m_a_orbs->begin());
  const auto ib = std::size_t(b - m_b_orbs->begin());
  return (m_aisb && ib > ia) ? m_y_abkr[ib][ia] : m_y_abkr[ia][ib];
}

//******************************************************************************
void YkTable::update_y_ints() {
  resize_y();
  const auto tj_max = max_tj();
  const auto k_max = tj_max; // k_max = |j + j'| = 2*j_max
  m_Ck.fill_maxK_twojmax(k_max, tj_max);

  a_size = m_a_orbs->size();
  b_size = m_b_orbs->size();

#pragma omp parallel for
  for (std::size_t ia = 0; ia < a_size; ia++) {
    const auto &Fa = (*m_a_orbs)[ia];
    const auto b_max = m_aisb ? ia : b_size - 1;
    for (std::size_t ib = 0; ib <= b_max; ib++) {
      const auto &Fb = (*m_b_orbs)[ib];
      const auto &[kmin, kmax] = k_minmax(Fa, Fb); // weird that this works?
      for (auto k = kmin; k <= kmax; k++) {
        const auto Lk = m_Ck.get_Lambdakab(k, Fa.k, Fb.k);
        if (Lk == 0)
          continue;
        const auto ik = std::size_t(k - kmin);
        Coulomb::yk_ab(Fa, Fb, k, m_y_abkr[ia][ib][ik]);
      } // k
    }   // b
  }     // a
}
//******************************************************************************
void YkTable::update_y_ints(const DiracSpinor &Fn) {
  //
  bool nisa = true;
  auto n = std::find(m_a_orbs->begin(), m_a_orbs->end(), Fn);
  if (n == m_a_orbs->end()) {
    nisa = false;
    n = std::find(m_b_orbs->begin(), m_b_orbs->end(), Fn);
  }
  assert(n != m_b_orbs->end());

  const auto in = nisa ? std::size_t(n - m_a_orbs->begin())
                       : std::size_t(n - m_b_orbs->begin());
  const auto &m_orbs = nisa ? *m_b_orbs : *m_a_orbs;
  const auto m_size = m_orbs.size();

#pragma omp parallel for
  for (std::size_t im = 0; im < m_size; im++) {
    const auto &Fm = m_orbs[im];
    const auto &[kmin, kmax] = k_minmax(Fm, Fn);
    for (auto k = kmin; k <= kmax; k++) {
      const auto Lk = m_Ck.get_Lambdakab(k, Fn.k, Fm.k);
      if (Lk == 0)
        continue;
      const auto ik = std::size_t(k - kmin);
      if (m_aisb) {
        if (in > im)
          Coulomb::yk_ab(Fm, Fn, k, m_y_abkr[in][im][ik]);
        else
          Coulomb::yk_ab(Fm, Fn, k, m_y_abkr[im][in][ik]);
      } else {
        if (nisa)
          Coulomb::yk_ab(Fm, Fn, k, m_y_abkr[in][im][ik]);
        else
          Coulomb::yk_ab(Fm, Fn, k, m_y_abkr[im][in][ik]);
      } // misa
    }   // k
  }     // m
}

//******************************************************************************
void YkTable::resize_y() {
  if (m_a_orbs->size() == a_size && m_b_orbs->size() == b_size)
    return;
  a_size = m_a_orbs->size();
  b_size = m_b_orbs->size();

  m_y_abkr.resize(a_size);
  for (std::size_t ia = 0; ia < a_size; ia++) {
    const auto &Fa = (*m_a_orbs)[ia];
    const auto b_max = m_aisb ? ia : b_size - 1;
    m_y_abkr[ia].resize(b_max + 1);
    for (std::size_t ib = 0; ib <= b_max; ib++) {
      const auto &Fb = (*m_b_orbs)[ib];
      const auto &[kmin, kmax] = k_minmax(Fa, Fb);
      const auto num_k = std::size_t(kmax - kmin + 1);
      m_y_abkr[ia][ib].resize(num_k);
      for (auto &yk_ab_r : m_y_abkr[ia][ib]) {
        yk_ab_r.resize(m_grid->num_points);
      } // k
    }   // b
  }     // a
}

} // namespace Coulomb
