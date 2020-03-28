#include "Coulomb.hpp"
#include "Angular/Angular_tables.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

namespace Coulomb {

//******************************************************************************
// Templates for faster method to calculate r^k
template <int k> static inline double powk_new(const double x) {
  return x * powk_new<k - 1>(x);
}
template <> inline double powk_new<0>(const double) { return 1.0; }

//******************************************************************************
template <int k>
static inline void yk_ijk_impl(const int l, const DiracSpinor &Fa,
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
  const auto bmin = std::max(Fa.p0, Fb.p0);
  for (std::size_t i = bmin; i < bmax; i++) {
    Bx += w(i) * gr->drduor[i] * (Fa.f[i] * Fb.f[i] + Fa.g[i] * Fb.g[i]) /
          powk(gr->r[i]);
  }

  vabk[0] = Bx * du * powk(gr->r[0]);
  for (std::size_t i = 1; i < irmax; i++) {
    // XXX below: don't do if Fa[i-1] or Fb.f[i - 1] = 0 !
    // but DO need to calc vabk for "all"
    const auto rm1_to_k = powk(gr->r[i - 1]);
    const auto inv_rm1_to_kp1 = 1.0 / (rm1_to_k * gr->r[i - 1]);
    const auto r_to_k = powk(gr->r[i]);
    const auto inv_r_to_kp1 = 1.0 / (r_to_k * gr->r[i]);
    const auto Fdr = gr->drdu[i - 1] *
                     (Fa.f[i - 1] * Fb.f[i - 1] + Fa.g[i - 1] * Fb.g[i - 1]);
    const auto wi = w(i - 1);
    Ax += wi * Fdr * rm1_to_k;
    Bx -= wi * Fdr * inv_rm1_to_kp1;
    // above
    vabk[i] = du * (Ax * inv_r_to_kp1 + Bx * r_to_k);
  }
  for (std::size_t i = irmax; i < num_points; i++) {
    vabk[i] = 0.0;
  }
}

//******************************************************************************
std::vector<double> yk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb,
                          const int k, const std::size_t maxi) {
  std::vector<double> ykab; //
  yk_ab(Fa, Fb, k, ykab, maxi);
  return ykab;
}

//------------------------------------------------------------------------------
void yk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &vabk, const std::size_t maxi) {
  auto sp1 = SafeProfiler::profile(__func__);

  if (k == 0)
    yk_ijk_impl<0>(k, Fa, Fb, vabk, maxi);
  else if (k == 1)
    yk_ijk_impl<1>(k, Fa, Fb, vabk, maxi);
  else if (k == 2)
    yk_ijk_impl<2>(k, Fa, Fb, vabk, maxi);
  else if (k == 3)
    yk_ijk_impl<3>(k, Fa, Fb, vabk, maxi);
  else if (k == 4)
    yk_ijk_impl<4>(k, Fa, Fb, vabk, maxi);
  else if (k == 5)
    yk_ijk_impl<5>(k, Fa, Fb, vabk, maxi);
  else if (k == 6)
    yk_ijk_impl<6>(k, Fa, Fb, vabk, maxi);
  else if (k == 7)
    yk_ijk_impl<7>(k, Fa, Fb, vabk, maxi);
  else if (k == 8)
    yk_ijk_impl<8>(k, Fa, Fb, vabk, maxi);
  else
    yk_ijk_impl<-1>(k, Fa, Fb, vabk, maxi);
}

//******************************************************************************
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd,
               const int k) //
{
  auto sp1 = SafeProfiler::profile(__func__);
  const auto yk_bd = yk_ab(Fb, Fd, k, std::min(Fa.pinf, Fc.pinf));
  return Rk_abcd(Fa, Fc, yk_bd);
}
//------------------------------------------------------------------------------
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fc,
               const std::vector<double> &yk_bd) {
  auto sp1 = SafeProfiler::profile(__func__, "yk");
  const auto &drdu = Fa.p_rgrid->drdu;
  const auto i0 = std::max(Fa.p0, Fc.p0);
  const auto imax = std::min(Fa.pinf, Fc.pinf);
  const auto Rff = NumCalc::integrate(1.0, i0, imax, Fa.f, Fc.f, yk_bd, drdu);
  const auto Rgg = NumCalc::integrate(1.0, i0, imax, Fa.g, Fc.g, yk_bd, drdu);
  return (Rff + Rgg) * Fa.p_rgrid->du;
}

//******************************************************************************
DiracSpinor Rk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                        const DiracSpinor &Fc, const DiracSpinor &Fd,
                        const int k) {
  auto sp1 = SafeProfiler::profile(__func__);
  const auto ykbd = yk_ab(Fb, Fd, k, Fc.pinf);
  return Rk_abcd_rhs(Fa.k, Fc, ykbd);
}
//------------------------------------------------------------------------------
DiracSpinor Rk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fc,
                        const std::vector<double> &ykbd) {
  // auto sp1 = SafeProfiler::profile(__func__);
  // auto out = DiracSpinor(0, Fa.k, *(Fa.p_rgrid));
  // out.p0 = Fc.p0;
  // out.pinf = Fc.pinf;
  // // out.f = NumCalc::mult_vectors(Fc.f, ykbd);
  // // out.g = NumCalc::mult_vectors(Fc.g, ykbd);
  // return ykbd * out;
  return Rk_abcd_rhs(Fa.k, Fc, ykbd);
}
//------------------------------------------------------------------------------
DiracSpinor Rk_abcd_rhs(const int kappa_a, const DiracSpinor &Fc,
                        const std::vector<double> &ykbd) {
  auto sp1 = SafeProfiler::profile(__func__);
  auto out = DiracSpinor(0, kappa_a, *(Fc.p_rgrid));
  out.p0 = Fc.p0;
  out.pinf = Fc.pinf;
  out.f = NumCalc::mult_vectors(Fc.f, ykbd);
  out.g = NumCalc::mult_vectors(Fc.g, ykbd);
  return out;
}
//------------------------------------------------------------------------------
void Rk_abcd_rhs(DiracSpinor *const Rkv, const DiracSpinor &Fc,
                 const std::vector<double> &ykbd) {
  Rkv->p0 = Fc.p0;
  Rkv->pinf = Fc.pinf;
  for (auto i = 0ul; i < Rkv->p0; ++i) {
    Rkv->f[i] = 0.0;
    Rkv->g[i] = 0.0;
  }
  for (auto i = Rkv->p0; i < Rkv->pinf; ++i) {
    Rkv->f[i] = Fc.f[i] * ykbd[i];
    Rkv->g[i] = Fc.g[i] * ykbd[i];
  }
  for (auto i = Rkv->pinf; i < Rkv->p_rgrid->num_points; ++i) {
    Rkv->f[i] = 0.0;
    Rkv->g[i] = 0.0;
  }
}

//******************************************************************************
double Zk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  // Z^k_abcd = s ( Q^k_abcd + sum_l [k] 6j * Q^l_abdc)
  auto sp1 = SafeProfiler::profile(__func__);
  const auto s = Angular::evenQ_2(Fa.twoj() + Fb.twoj() + 2) ? 1 : -1;
  return s * Wk_abcd(Fa, Fb, Fc, Fd, k);
}

//******************************************************************************
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  auto sp1 = SafeProfiler::profile(__func__);
  const auto Qkabcd = Qk_abcd(Fa, Fb, Fc, Fd, k);
  const auto Pkabcd = Pk_abcd(Fa, Fb, Fc, Fd, k);
  return (Qkabcd + Pkabcd);
}
//******************************************************************************
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  auto sp1 = SafeProfiler::profile(__func__);

  const auto tkp1 = 2 * k + 1;
  const auto min_twol = std::max(std::abs(Fd.twoj() - Fa.twoj()),
                                 std::abs(Fc.twoj() - Fb.twoj()));
  const auto max_twol = std::min(Fd.twoj() + Fa.twoj(), Fc.twoj() + Fb.twoj());
  double sum = 0.0;
  for (int tl = min_twol; tl <= max_twol; tl += 2) {
    const auto sixj = Angular::sixj_2(Fc.twoj(), Fa.twoj(), 2 * k, //
                                      Fd.twoj(), Fb.twoj(), tl);
    if (sixj == 0)
      continue;
    const auto Qlabdc = Qk_abcd(Fa, Fb, Fd, Fc, tl / 2);
    sum += sixj * Qlabdc;
  }
  return (tkp1 * sum);
}

//******************************************************************************
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<double> &ykbd,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  auto sp1 = SafeProfiler::profile(__func__, "yk");

  const auto Qkabcd = Qk_abcd(Fa, Fb, Fc, Fd, k, ykbd, Ck);
  const auto Pkabcd = Pk_abcd(Fa, Fb, Fc, Fd, k, ybc, Ck, sixj);

  return (Qkabcd + Pkabcd);
}
//******************************************************************************
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               // const std::vector<double> &ykbd,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj) {
  // W^k_abcd = Q^k_abcd +
  // W = Q + P
  // P_abcd = sum_l [k] 6j * Q^l_abdc //not index order!
  auto sp1 = SafeProfiler::profile(__func__, "yk");

  const auto tkp1 = 2 * k + 1;
  double sum = 0.0;
  const auto min_l = std::abs(Fb.twoj() - Fc.twoj()) / 2;
  auto count = 0;
  for (const auto &ybc_l : ybc) {
    const auto l = min_l + count;
    ++count;
    const auto sj =
        sixj.get_6j(Fc.twoj(), Fa.twoj(), Fd.twoj(), Fb.twoj(), k, l);
    if (sj == 0.0)
      continue;
    sum += sj * Qk_abcd(Fa, Fb, Fd, Fc, l, ybc_l, Ck); // a,b,d,c [exch.]
  }

  return (tkp1 * sum);
}
//------------------------------------------------------------------------------
void Pkv_bcd(DiracSpinor *Pkv, const DiracSpinor &Fb, const DiracSpinor &Fc,
             const DiracSpinor &Fd, const int k,
             const std::vector<std::vector<double>> &ybc,
             const Angular::Ck_ab &Ck, const Angular::SixJ &sixj) {
  // W^k_abcd = Q^k_abcd +
  // W = Q + P
  // P_abcd = sum_l [k] 6j * Q^l_abdc //not index order!
  auto sp1 = SafeProfiler::profile(__func__, "yk");

  *Pkv *= 0.0;
  // Pkv->p0 = 0;
  // Pkv->pinf=Pkv->p0=0;

  const auto kappa = Pkv->k;

  const auto tkp1 = 2 * k + 1;
  // double sum = 0.0;
  const auto min_l = std::abs(Fb.twoj() - Fc.twoj()) / 2;
  auto count = 0;
  for (const auto &ybc_l : ybc) {
    const auto l = min_l + count;
    ++count;
    const auto sj = sixj.get_6j(Fc.twoj(), Angular::twoj_k(kappa), Fd.twoj(),
                                Fb.twoj(), k, l);
    if (sj == 0.0)
      continue;
    *Pkv += sj * Qkv_bcd(kappa, Fb, Fd, Fc, l, ybc_l, Ck); // a,b,d,c [exch.]
  }
  *Pkv *= tkp1;
}

//******************************************************************************
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  auto sp1 = SafeProfiler::profile(__func__);
  const auto tCac = Angular::tildeCk_kk(k, Fa.k, Fc.k);
  if (tCac == 0.0)
    return 0.0;
  const auto tCbd = Angular::tildeCk_kk(k, Fb.k, Fd.k);
  if (tCbd == 0.0)
    return 0.0;
  const auto Rkabcd = Rk_abcd(Fa, Fb, Fc, Fd, k);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rkabcd;
}

//******************************************************************************
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<double> &ykbd, const Angular::Ck_ab &Ck) {
  auto sp1 = SafeProfiler::profile(__func__, "yk");

  const auto tCac = Ck.get_tildeCkab(k, Fa.k, Fc.k);
  if (tCac == 0.0)
    return 0.0;
  const auto tCbd = Ck.get_tildeCkab(k, Fb.k, Fd.k);
  if (tCbd == 0.0)
    return 0.0;
  const auto Rkabcd = Rk_abcd(Fa, Fc, ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rkabcd;
}

//******************************************************************************
DiracSpinor Qkv_bcd(const int kappa_v, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
                    const std::vector<double> &ykbd, const Angular::Ck_ab &Ck) {
  auto sp1 = SafeProfiler::profile(__func__, "yk");
  const auto tCac = Ck.get_tildeCkab(k, kappa_v, Fc.k);
  const auto tCbd = Ck.get_tildeCkab(k, Fb.k, Fd.k);
  const auto tCC = tCbd * tCac;
  auto Rkv = tCC == 0 ? DiracSpinor(0, kappa_v, *(Fb.p_rgrid))
                      : Rk_abcd_rhs(kappa_v, Fc, ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  Rkv.scale(m1tk * tCC);
  return Rkv;
}

//******************************************************************************
void Qkv_bcd(DiracSpinor *const Qkv, const DiracSpinor &Fb,
             const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
             const std::vector<double> &ykbd, const Angular::Ck_ab &Ck) {
  auto sp1 = SafeProfiler::profile(__func__, "yk");
  const auto tCac = Ck.get_tildeCkab(k, Qkv->k, Fc.k);
  const auto tCbd = Ck.get_tildeCkab(k, Fb.k, Fd.k);
  const auto tCC = tCbd * tCac;
  if (tCC == 0.0) {
    Qkv->scale(0.0);
    return;
  }
  Rk_abcd_rhs(Qkv, Fc, ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  Qkv->scale(m1tk * tCC);
  return;
}

// //******************************************************************************
// double Xk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
//                const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
//   // implement in terms of Q instead ..but check!
//   auto sp1 = SafeProfiler::profile(__func__);
//   const auto tCac = Angular::Ck_kk(k, Fa.k, Fc.k);
//   if (tCac == 0.0)
//     return 0.0;
//   const auto tCbd = Angular::Ck_kk(k, Fb.k, Fd.k);
//   if (tCbd == 0.0)
//     return 0.0;
//   const auto Rkabcd = Rk_abcd(Fa, Fb, Fc, Fd, k);
//   const auto m1tk = Angular::evenQ(k) ? 1 : -1;
//   return m1tk * tCac * tCbd * Rkabcd;
// }

} // namespace Coulomb
