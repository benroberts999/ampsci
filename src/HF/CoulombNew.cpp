#include "CoulombNew.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "HF/CoulombInts.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

//******************************************************************************
YkTable::YkTable(const Grid *const in_grid,
                 const std::vector<DiracSpinor> *const in_a_orbs,
                 const std::vector<DiracSpinor> *const in_b_orbs)
    : m_a_orbs(in_a_orbs),                                    //
      m_b_orbs(in_b_orbs == nullptr ? in_a_orbs : in_b_orbs), //
      m_grid(in_grid),                                        //
      m_aisb([&]() {
        return (in_b_orbs == nullptr || in_a_orbs == in_b_orbs) ? true : false;
      }()) //
{
  update_y_ints();
}

//******************************************************************************
int YkTable::max_tj() const {
  auto comp2j = [](const auto &a, const auto &b) {
    return a.twoj() < b.twoj();
  };
  if (m_a_orbs->size() == 0)
    return 0;
  const auto maxtj_a =
      std::max_element(m_a_orbs->cbegin(), m_a_orbs->cend(), comp2j)->twoj();
  if (m_aisb || m_b_orbs->size() == 0)
    return maxtj_a;
  const auto maxtj_b =
      std::max_element(m_b_orbs->cbegin(), m_b_orbs->cend(), comp2j)->twoj();
  return std::max(maxtj_a, maxtj_b);
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
  const auto [kmin, kmax] = k_minmax(Fa, Fb);
  const auto ik = std::size_t(k - kmin);
  const auto &yab = get_y_ab(Fa, Fb);
  if constexpr (check_bounds) {
    if (k < kmin || k > kmax || ik > yab.size()) {
      std::cerr << "Fail 35 in Coulomb: k too big/small: " << k << ": " << kmin
                << "/" << kmax << " " << yab.size() << "\n";
      std::abort();
    }
  }
  return yab[ik];
}
//------------------------------------------------------------------------------
const std::vector<std::vector<double>> &
YkTable::get_y_ab(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  const auto a = std::find(m_a_orbs->begin(), m_a_orbs->end(), Fa);
  const auto b = std::find(m_b_orbs->begin(), m_b_orbs->end(), Fb);
  if constexpr (check_bounds) {
    if (a == m_a_orbs->end()) {
      std::cerr << "Fail 40 in Coulomb: Fa not in a: " << Fa.symbol() << "\n";
      std::abort();
    }
    if (b == m_b_orbs->end()) {
      std::cerr << "Fail 44 in Coulomb: Fb not in b: " << Fb.symbol() << "\n";
      std::abort();
    }
  }
  const std::size_t ia = a - m_a_orbs->begin();
  const std::size_t ib = b - m_b_orbs->begin();
  return (m_aisb && ib > ia) ? m_y_abkr[ib][ia] : m_y_abkr[ia][ib];
}

//******************************************************************************
void YkTable::update_y_ints() { //
  resize_y();
  const auto tj_max = max_tj();
  const auto k_max = tj_max;
  // k_min = |j - j'|; k_max = |j + j'|
  m_Ck.fill_maxK_twojmax(k_max, tj_max);

  a_size = m_a_orbs->size();
  b_size = m_b_orbs->size();

  // #pragma omp parallel for
  for (std::size_t ia = 0; ia < a_size; ia++) {
    std::cerr << __LINE__ << "\n";
    const auto &Fa = (*m_a_orbs)[ia];
    std::cerr << __LINE__ << "\n";
    std::cerr << Fa.symbol() << " -- ";
    const auto b_max = m_aisb ? ia : b_size - 1;
    for (std::size_t ib = 0; ib <= b_max; ib++) {
      const auto &Fb = (*m_b_orbs)[ib];
      std::cerr << Fb.symbol() << "\n";
      auto rmaxi = (Fb == Fa) ? 0 : std::min(Fa.pinf, Fb.pinf); // XXX check?
      const auto [kmin, kmax] = k_minmax(Fa, Fb);
      for (auto k = kmin; k <= kmax; k++) {
        std::cerr << k << " " << __LINE__ << "\n";
        const auto Lk = m_Ck.get_Lambdakab(k, Fa.k, Fb.k);
        std::cerr << __LINE__ << "\n";
        if (Lk == 0)
          continue;
        const auto ik = std::size_t(k - kmin);
        std::cerr << "\n" << __LINE__ << "\n";
        CoulombInts::yk_ab(Fa, Fb, k, m_y_abkr[ia][ib][ik], rmaxi);
        std::cerr << __LINE__ << "\n";
      }
    }
  }
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
  if constexpr (check_bounds) {
    if (n == m_b_orbs->end()) {
      std::cerr << "Fail 108 in Coulomb: Fn not in a or b: " << Fn.symbol()
                << "\n";
      std::abort();
    }
  }

  const auto in = nisa ? std::size_t(n - m_a_orbs->begin())
                       : std::size_t(n - m_b_orbs->begin());
  const auto &m_orbs = nisa ? *m_b_orbs : *m_a_orbs;
  const auto m_size = m_orbs.size();

#pragma omp parallel for
  for (std::size_t im = 0; im < m_size; im++) {
    const auto &Fm = m_orbs[im];
    auto rmaxi = (Fm == Fn) ? 0 : std::min(Fm.pinf, Fn.pinf); // XXX check?
    const auto [kmin, kmax] = k_minmax(Fm, Fn);
    //
    for (auto k = kmin; k <= kmax; k++) {
      const auto Lk = m_Ck.get_Lambdakab(k, Fn.k, Fm.k);
      if (Lk == 0)
        continue;
      const auto ik = std::size_t(k - kmin);
      if (m_aisb) {
        if (in > im)
          CoulombInts::yk_ab(Fm, Fn, k, m_y_abkr[in][im][ik], rmaxi);
        else
          CoulombInts::yk_ab(Fm, Fn, k, m_y_abkr[im][in][ik], rmaxi);
      } else {
        if (nisa)
          CoulombInts::yk_ab(Fm, Fn, k, m_y_abkr[in][im][ik], rmaxi);
        else
          CoulombInts::yk_ab(Fm, Fn, k, m_y_abkr[im][in][ik], rmaxi);
      } // misa
    }   // k
  }     // m
}

//******************************************************************************
void YkTable::resize_y() { //
  if (m_a_orbs->size() == a_size && m_b_orbs->size() == b_size)
    return; // XXX
  a_size = m_a_orbs->size();
  b_size = m_b_orbs->size();

  m_y_abkr.resize(a_size);
  for (std::size_t ia = 0; ia < a_size; ia++) {
    const auto &Fa = (*m_a_orbs)[ia];
    const auto b_max = m_aisb ? ia : b_size - 1;
    m_y_abkr[ia].resize(b_max + 1);
    for (std::size_t ib = 0; ib <= b_max; ib++) {
      const auto &Fb = (*m_b_orbs)[ib];
      const auto [kmin, kmax] = k_minmax(Fa, Fb);
      const auto num_k = kmax - kmin + 1;
      m_y_abkr[ia][ib].resize(num_k);
      for (auto &yakb_r : m_y_abkr[ia][ib]) {
        yakb_r.resize(m_grid->num_points);
      }
    }
  }
}
