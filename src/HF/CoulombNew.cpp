#include "CoulombNew.hpp"
#include "Angular/Angular_tables.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

//******************************************************************************
int CoulombNew::max_tj() const {
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
std::pair<int, int> CoulombNew::k_minmax(const DiracSpinor &a,
                                         const DiracSpinor &b) {
  return std::make_pair(std::abs(a.twoj() - b.twoj()) / 2,
                        (a.twoj() + b.twoj()) / 2);
}

//******************************************************************************
const std::vector<double> &CoulombNew::get_yk_ab(const int k,
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
CoulombNew::get_y_ab(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
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
void CoulombNew::update_y_ints() { //
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
        calculate_y_ijk(Fa, Fb, k, m_y_abkr[ia][ib][ik], rmaxi);
        std::cerr << __LINE__ << "\n";
      }
    }
  }
}
//******************************************************************************
void CoulombNew::update_y_ints(const DiracSpinor &Fn) {
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
          calculate_y_ijk(Fm, Fn, k, m_y_abkr[in][im][ik], rmaxi);
        else
          calculate_y_ijk(Fm, Fn, k, m_y_abkr[im][in][ik], rmaxi);
      } else {
        if (nisa)
          calculate_y_ijk(Fm, Fn, k, m_y_abkr[in][im][ik], rmaxi);
        else
          calculate_y_ijk(Fm, Fn, k, m_y_abkr[im][in][ik], rmaxi);
      } // misa
    }   // k
  }     // m
}

//******************************************************************************
void CoulombNew::resize_y() { //
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

//******************************************************************************
//******************************************************************************
// Templates for faster method to calculate r^k
template <int k> static inline double powk_new(const double x) {
  return x * powk_new<k - 1>(x);
}
template <> inline double powk_new<0>(const double) { return 1.0; }

//******************************************************************************
template <int k>
static inline void yk_ijk(const int l, const DiracSpinor &Fa,
                          const DiracSpinor &Fb, std::vector<double> &vabk,
                          const std::size_t maxi)
// Calculalates y^k_ab screening function.
// Note: is symmetric: y_ab = y_ba
//
// Stores in vabk (in/out parameter, reference to whatever)
//
// r_min := min(r,r')
// rho(r') := fa(r')*fb(r') + ga(r')gb(r')
// y^k_ab(r) = Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
//           = Int_0^r [r'^k/r^(k+1)]*rho(r') dr'
//             + Int_r^inf [r^k/r'^(k+1)]*rho(r') dr'
//          := A(r)/r^(k+1) + B(r)*r^k
// A(r0)  = 0
// B(r0)  = Int_0^inf [r^k/r'^(k+1)]*rho(r') dr'
// A(r_n) = A(r_{n-1}) + (rho(r_{n-1})*r_{n-1}^k)*dr
// B(r_n) = A(r_{n-1}) + (rho(r_{n-1})/r_{n-1}^(k+1))*dr
// y^k_ab(rn) = A(rn)/rn^(k+1) + B(rn)*rn^k
//
// Also uses Quadrature integer rules! (Defined in NumCalc)
{
  const auto &gr = Fa.p_rgrid; // just save typing
  const auto du = gr->du;
  const auto num_points = gr->num_points;
  vabk.resize(num_points); // for safety
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;

  // faster method to calculate r^k
  auto powk = [=](double x) {
    if constexpr (k < 0)
      return std::pow(x, l);
    else
      return powk_new<k>(x);
  };

  // Quadrature integration weights:
  auto w = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  double Ax = 0.0, Bx = 0.0;

  const auto bmax = std::min(Fa.pinf, Fb.pinf);
  for (std::size_t i = 0; i < bmax; i++) {
    Bx += w(i) * gr->drduor[i] * (Fa.f[i] * Fb.f[i] + Fa.g[i] * Fb.g[i]) /
          powk(gr->r[i]);
  }

  vabk[0] = Bx * du * powk(gr->r[0]);
  for (std::size_t i = 1; i < irmax; i++) {
    const auto rm1_to_k = powk(gr->r[i - 1]);
    const auto inv_rm1_to_kp1 = 1.0 / (rm1_to_k * gr->r[i - 1]);
    const auto r_to_k = powk(gr->r[i]);
    const auto inv_r_to_kp1 = 1.0 / (r_to_k * gr->r[i]);
    const auto Fdr = gr->drdu[i - 1] *
                     (Fa.f[i - 1] * Fb.f[i - 1] + Fa.g[i - 1] * Fb.g[i - 1]);
    const auto wi = w(i - 1);
    Ax += wi * Fdr * rm1_to_k;
    Bx -= wi * Fdr * inv_rm1_to_kp1;
    vabk[i] = du * (Ax * inv_r_to_kp1 + Bx * r_to_k);
  }
  for (std::size_t i = irmax; i < num_points; i++) {
    vabk[i] = 0.0;
  }
}

//******************************************************************************
void CoulombNew::calculate_y_ijk(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                 const int k, std::vector<double> &vabk,
                                 const std::size_t maxi)
// This is static
{
  auto sp1 = SafeProfiler::profile(__func__);

  if (k == 0)
    yk_ijk<0>(k, Fa, Fb, vabk, maxi);
  else if (k == 1)
    yk_ijk<1>(k, Fa, Fb, vabk, maxi);
  else if (k == 2)
    yk_ijk<2>(k, Fa, Fb, vabk, maxi);
  else if (k == 3)
    yk_ijk<3>(k, Fa, Fb, vabk, maxi);
  else if (k == 4)
    yk_ijk<4>(k, Fa, Fb, vabk, maxi);
  else if (k == 5)
    yk_ijk<5>(k, Fa, Fb, vabk, maxi);
  else if (k == 6)
    yk_ijk<6>(k, Fa, Fb, vabk, maxi);
  else if (k == 7)
    yk_ijk<7>(k, Fa, Fb, vabk, maxi);
  else if (k == 8)
    yk_ijk<8>(k, Fa, Fb, vabk, maxi);
  else
    yk_ijk<-1>(k, Fa, Fb, vabk, maxi);
}
