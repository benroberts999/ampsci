#include "Coulomb.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

namespace Coulomb {

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
// --> nb: slightly different to above description..
//
// Also uses Quadrature integer rules! (Defined in NumCalc)
{
  const auto &gr = Fa.rgrid; // just save typing
  const auto du = gr->du();
  const auto num_points = gr->num_points();
  vabk.resize(num_points); // for safety
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;

  // faster method to calculate r^k
  const auto powk = [l]() {
    if constexpr (k < 0) {
      return [l](double x) { return std::pow(x, l); };
    } else {
      (void)l; // l not used
      return qip::pow<k, double>;
    }
  }();

  // Quadrature integration weights:
  const auto w = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto ff = [&](std::size_t i) {
    return (Fa.f(i) * Fb.f(i) + Fa.g(i) * Fb.g(i)) * w(i) * gr->drduor()[i];
  };
  const auto &r = gr->r();

  double Ax = 0.0, Bx = 0.0;

  vabk[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {
    const auto rat = r[i - 1] / r[i];
    Ax = (Ax + ff(i - 1)) * (rat * powk(rat));
    vabk[i] = Ax * du;
  }

  const auto bmax =
      std::min(std::min(Fa.max_pt(), Fb.max_pt()), num_points - 1);
  for (auto i = bmax; i >= 1; --i) {
    Bx = Bx * powk(r[i - 1] / r[i]) + ff(i - 1);
    vabk[i - 1] += Bx * du;
  }

  for (std::size_t i = irmax; i < num_points; i++) {
    vabk[i] = 0.0;
  }
}

//------------------------------------------------------------------------------
std::vector<double> yk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb,
                          const int k, const std::size_t maxi) {
  std::vector<double> ykab; //
  yk_ab(Fa, Fb, k, ykab, maxi);
  return ykab;
}

//------------------------------------------------------------------------------
void yk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &vabk, const std::size_t maxi) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // faster method to calculate r^k
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
  else if (k == 9)
    yk_ijk_impl<9>(k, Fa, Fb, vabk, maxi);
  else if (k == 10)
    yk_ijk_impl<10>(k, Fa, Fb, vabk, maxi);
  else
    yk_ijk_impl<-1>(k, Fa, Fb, vabk, maxi);
}

//******************************************************************************
template <int k, int pm>
static inline void Breit_abk_impl(const int l, const DiracSpinor &Fa,
                                  const DiracSpinor &Fb, //
                                  std::vector<double> &b0,
                                  std::vector<double> &binf,
                                  const std::size_t maxi) {
  static_assert(pm == 1 || pm == -1,
                "Breit_abk_impl must be called with pm=+/-1 only\n");
  const auto &gr = Fa.rgrid; // just save typing
  const auto du = gr->du();
  const auto num_points = gr->num_points();
  b0.resize(num_points);   // needed
  binf.resize(num_points); // needed
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;

  // faster method to calculate r^k
  const auto powk = [l]() {
    if constexpr (k < 0) {
      return [l](double x) { return std::pow(x, l); };
    } else {
      (void)l; // don't use l..
      return qip::pow<k, double>;
    }
  }();

  // Quadrature integration weights:
  auto w = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  double Ax = 0.0, Bx = 0.0;

  auto fgfg = [&](std::size_t i) {
    if constexpr (pm == -1)
      return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
    else
      return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  const auto bmax = std::min(Fa.max_pt(), Fb.max_pt());
  const auto bmin = std::max(Fa.min_pt(), Fb.min_pt());
  for (std::size_t i = bmin; i < bmax; i++) {
    Bx += gr->drduor()[i] * w(i) * fgfg(i) / powk(gr->r()[i]);
  }

  b0[0] = 0.0;
  binf[0] = Bx * du * powk(gr->r()[0]);
  for (std::size_t i = 1; i < irmax; i++) {
    const auto rm1_to_k = powk(gr->r()[i - 1]);
    const auto inv_rm1_to_kp1 = 1.0 / (rm1_to_k * gr->r()[i - 1]);
    const auto r_to_k = powk(gr->r()[i]);
    const auto inv_r_to_kp1 = 1.0 / (r_to_k * gr->r()[i]);
    const auto Fdr = gr->drdu()[i - 1] * fgfg(i - 1) * w(i - 1);
    Ax += Fdr * rm1_to_k;
    Bx -= Fdr * inv_rm1_to_kp1;
    b0[i] = du * Ax * inv_r_to_kp1;
    binf[i] = du * Bx * r_to_k;
  }
  // std::cout << "\n";
  // for (std::size_t i = irmax; i < num_points; i++) {
  //   b0[i] = 0.0;
  //   binf[i] = 0.0;
  // }
}
//------------------------------------------------------------------------------
void bk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &b0, std::vector<double> &binf,
           const std::size_t maxi) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);

  constexpr int pm = -1; // for fg - gf
  // faster method to calculate r^k
  if (k == 0)
    Breit_abk_impl<0, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 1)
    Breit_abk_impl<1, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 2)
    Breit_abk_impl<2, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 3)
    Breit_abk_impl<3, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 4)
    Breit_abk_impl<4, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 5)
    Breit_abk_impl<5, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 6)
    Breit_abk_impl<6, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 7)
    Breit_abk_impl<7, pm>(k, Fa, Fb, b0, binf, maxi);
  else if (k == 8)
    Breit_abk_impl<8, pm>(k, Fa, Fb, b0, binf, maxi);
  else
    Breit_abk_impl<-1, pm>(k, Fa, Fb, b0, binf, maxi);
}
//------------------------------------------------------------------------------
void gk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &g0, std::vector<double> &ginf,
           const std::size_t maxi) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);

  constexpr int pm = +1; // for fg + gf
  if (k == 0)
    Breit_abk_impl<0, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 1)
    Breit_abk_impl<1, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 2)
    Breit_abk_impl<2, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 3)
    Breit_abk_impl<3, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 4)
    Breit_abk_impl<4, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 5)
    Breit_abk_impl<5, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 6)
    Breit_abk_impl<6, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 7)
    Breit_abk_impl<7, pm>(k, Fa, Fb, g0, ginf, maxi);
  else if (k == 8)
    Breit_abk_impl<8, pm>(k, Fa, Fb, g0, ginf, maxi);
  else
    Breit_abk_impl<-1, pm>(k, Fa, Fb, g0, ginf, maxi);
}

//******************************************************************************
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd,
               const int k) //
{
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  const auto yk_bd = yk_ab(Fb, Fd, k, std::min(Fa.max_pt(), Fc.max_pt()));
  return Rk_abcd(Fa, Fc, yk_bd);
}
//------------------------------------------------------------------------------
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fc,
               const std::vector<double> &yk_bd) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");
  const auto &drdu = Fa.rgrid->drdu();
  const auto i0 = std::max(Fa.min_pt(), Fc.min_pt());
  const auto imax = std::min(Fa.max_pt(), Fc.max_pt());
  const auto Rff =
      NumCalc::integrate(1.0, i0, imax, Fa.f(), Fc.f(), yk_bd, drdu);
  const auto Rgg =
      NumCalc::integrate(1.0, i0, imax, Fa.g(), Fc.g(), yk_bd, drdu);
  return (Rff + Rgg) * Fa.rgrid->du();
}

//******************************************************************************
DiracSpinor Rkv_bcd(const int kappa_a, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  return Rkv_bcd(kappa_a, Fc, yk_ab(Fb, Fd, k, Fc.max_pt()));
}
//------------------------------------------------------------------------------
DiracSpinor Rkv_bcd(const int kappa_a, const DiracSpinor &Fc,
                    const std::vector<double> &ykbd) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  auto out = DiracSpinor(0, kappa_a, Fc.rgrid);
  out.set_min_pt() = Fc.min_pt();
  out.set_max_pt() = Fc.max_pt();
  out.set_f() = qip::multiply(Fc.f(), ykbd);
  out.set_g() = qip::multiply(Fc.g(), ykbd);
  return out;
}
//------------------------------------------------------------------------------
void Rkv_bcd(DiracSpinor *const Rkv, const DiracSpinor &Fc,
             const std::vector<double> &ykbd) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  Rkv->set_min_pt() = Fc.min_pt();
  Rkv->set_max_pt() = Fc.max_pt();
  for (auto i = 0ul; i < Rkv->min_pt(); ++i) {
    Rkv->set_f(i) = 0.0;
    Rkv->set_g(i) = 0.0;
  }
  for (auto i = Rkv->min_pt(); i < Rkv->max_pt(); ++i) {
    Rkv->set_f(i) = Fc.f(i) * ykbd[i];
    Rkv->set_g(i) = Fc.g(i) * ykbd[i];
  }
  for (auto i = Rkv->max_pt(); i < Rkv->rgrid->num_points(); ++i) {
    Rkv->set_f(i) = 0.0;
    Rkv->set_g(i) = 0.0;
  }
}

//******************************************************************************
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  const auto tCac = Angular::tildeCk_kk(k, Fa.k, Fc.k);
  if (Angular::zeroQ(tCac))
    return 0.0;
  const auto tCbd = Angular::tildeCk_kk(k, Fb.k, Fd.k);
  if (Angular::zeroQ(tCbd))
    return 0.0;
  const auto Rkabcd = Rk_abcd(Fa, Fb, Fc, Fd, k);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rkabcd;
}

//------------------------------------------------------------------------------
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<double> &ykbd, const Angular::Ck_ab &Ck) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");

  const auto tCac = Ck.get_tildeCkab(k, Fa.k, Fc.k);
  if (Angular::zeroQ(tCac))
    return 0.0;
  const auto tCbd = Ck.get_tildeCkab(k, Fb.k, Fd.k);
  if (Angular::zeroQ(tCbd))
    return 0.0;
  const auto Rkabcd = Rk_abcd(Fa, Fc, ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rkabcd;
}

//------------------------------------------------------------------------------
DiracSpinor Qkv_bcd(const int kappa_a, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  const auto tCac = Angular::tildeCk_kk(k, kappa_a, Fc.k);
  const auto tCbd = Angular::tildeCk_kk(k, Fb.k, Fd.k);
  if (Angular::zeroQ(tCbd) || Angular::zeroQ(tCac))
    return DiracSpinor(0, kappa_a, Fc.rgrid);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return (m1tk * tCac * tCbd) * Rkv_bcd(kappa_a, Fb, Fc, Fd, k);
}

//------------------------------------------------------------------------------
DiracSpinor Qkv_bcd(const int kappa_v, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
                    const std::vector<double> &ykbd, const Angular::Ck_ab &Ck) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");
  const auto tCac = Ck.get_tildeCkab(k, kappa_v, Fc.k);
  const auto tCbd = Ck.get_tildeCkab(k, Fb.k, Fd.k);
  const auto tCC = tCbd * tCac;
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  if (tCC == 0) //?????
    return DiracSpinor(0, kappa_v, Fb.rgrid);
  return (m1tk * tCC) * Rkv_bcd(kappa_v, Fc, ykbd);
}

//------------------------------------------------------------------------------
void Qkv_bcd(DiracSpinor *const Qkv, const DiracSpinor &Fb,
             const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
             const std::vector<double> &ykbd, const Angular::Ck_ab &Ck) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");
  const auto tCac = Ck.get_tildeCkab(k, Qkv->k, Fc.k);
  const auto tCbd = Ck.get_tildeCkab(k, Fb.k, Fd.k);
  const auto tCC = tCbd * tCac;
  if (tCC == 0.0) {
    Qkv->scale(0.0);
    return;
  }
  Rkv_bcd(Qkv, Fc, ykbd);
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  Qkv->scale(m1tk * tCC);
  return;
}

//******************************************************************************
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);

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
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               // const std::vector<double> &ykbd,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj) {
  // W^k_abcd = Q^k_abcd +
  // W = Q + P
  // P_abcd = sum_l [k] 6j * Q^l_abdc //not index order!
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");

  const auto tkp1 = 2 * k + 1;
  double sum = 0.0;
  const auto min_l = std::abs(Fb.twoj() - Fc.twoj()) / 2;
  auto count = 0;
  for (const auto &ybc_l : ybc) {
    const auto l = min_l + count;
    ++count;

    if (!Angular::Ck_kk_SR(l, Fb.k, Fc.k) || !Angular::Ck_kk_SR(l, Fa.k, Fd.k))
      continue;

    const auto sj =
        sixj.get_6j(Fc.twoj(), Fa.twoj(), Fd.twoj(), Fb.twoj(), k, l);
    if (Angular::zeroQ(sj))
      continue;
    sum += sj * Qk_abcd(Fa, Fb, Fd, Fc, l, ybc_l, Ck); // a,b,d,c [exch.]
  }

  return (tkp1 * sum);
}

//------------------------------------------------------------------------------
DiracSpinor Pkv_bcd(int kappa_a, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k) {
  auto out = DiracSpinor(0, kappa_a, Fc.rgrid);
  const auto tkp1 = 2 * k + 1;
  const auto tja = Angular::twoj_k(kappa_a);
  const auto min_twol =
      std::max(std::abs(Fd.twoj() - tja), std::abs(Fc.twoj() - Fb.twoj()));
  const auto max_twol = std::min(Fd.twoj() + tja, Fc.twoj() + Fb.twoj());
  for (int tl = min_twol; tl <= max_twol; tl += 2) {
    if (!Angular::Ck_kk_SR(tl / 2, Fb.k, Fc.k) ||
        !Angular::Ck_kk_SR(tl / 2, kappa_a, Fd.k))
      continue;
    const auto sixj =
        Angular::sixj_2(Fc.twoj(), tja, 2 * k, Fd.twoj(), Fb.twoj(), tl);
    if (sixj == 0)
      continue;
    out += sixj * Qkv_bcd(kappa_a, Fb, Fd, Fc, tl / 2);
  }
  return (tkp1 * out);
}

//------------------------------------------------------------------------------
void Pkv_bcd(DiracSpinor *Pkv, const DiracSpinor &Fb, const DiracSpinor &Fc,
             const DiracSpinor &Fd, const int k,
             const std::vector<std::vector<double>> &ybc,
             const Angular::Ck_ab &Ck, const Angular::SixJ &sixj) {
  // W^k_abcd = Q^k_abcd +
  // W = Q + P
  // P_abcd = sum_l [k] 6j * Q^l_abdc //not index order!
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");

  *Pkv *= 0.0;

  const auto kappa = Pkv->k;

  const auto tkp1 = 2 * k + 1;
  const auto min_l = std::abs(Fb.twoj() - Fc.twoj()) / 2;
  auto count = 0;
  for (const auto &ybc_l : ybc) {
    const auto l = min_l + count;
    ++count;

    if (!Angular::Ck_kk_SR(l, Fb.k, Fc.k))
      continue;

    const auto sj = sixj.get_6j(Fc.twoj(), Angular::twoj_k(kappa), Fd.twoj(),
                                Fb.twoj(), k, l);
    if (Angular::zeroQ(sj))
      continue;
    *Pkv += sj * Qkv_bcd(kappa, Fb, Fd, Fc, l, ybc_l, Ck); // a,b,d,c [exch.]
  }
  *Pkv *= tkp1;
}
//------------------------------------------------------------------------------
void Pkv_bcd_2(DiracSpinor *Pkv, const DiracSpinor &Fb, const DiracSpinor &Fc,
               const DiracSpinor &Fd, const int k,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj,
               const std::vector<double> &f2k) {
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");

  *Pkv *= 0.0;

  const auto fk = [&f2k](int l) {
    // nb: only screens l, k assumed done outside...
    if (l < int(f2k.size())) {
      return f2k[std::size_t(l)];
    }
    return 1.0;
  };

  const auto kappa = Pkv->k;

  const auto tkp1 = 2 * k + 1;
  const auto min_l = std::abs(Fb.twoj() - Fc.twoj()) / 2;
  auto count = 0;
  for (const auto &ybc_l : ybc) {
    const auto l = min_l + count;
    ++count;
    // Include screening factor here (early escape if zero)
    const auto sj = fk(l) * sixj.get_6j(Fc.twoj(), Angular::twoj_k(kappa),
                                        Fd.twoj(), Fb.twoj(), k, l);
    if (Angular::zeroQ(sj))
      continue;
    *Pkv += sj * Qkv_bcd(kappa, Fb, Fd, Fc, l, ybc_l, Ck); // a,b,d,c [exch.]
  }
  *Pkv *= tkp1;
}

//******************************************************************************
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<double> &ykbd,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__, "yk");

  const auto Qkabcd = Qk_abcd(Fa, Fb, Fc, Fd, k, ykbd, Ck);
  const auto Pkabcd = Pk_abcd(Fa, Fb, Fc, Fd, k, ybc, Ck, sixj);

  return (Qkabcd + Pkabcd);
}

//------------------------------------------------------------------------------
DiracSpinor Wkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k) {
  auto out = Pkv_bcd(kappa_v, Fb, Fc, Fd, k);
  if (Angular::Ck_kk_SR(k, kappa_v, Fc.k) && Angular::Ck_kk_SR(k, Fb.k, Fd.k)) {
    out += Qkv_bcd(kappa_v, Fb, Fc, Fd, k);
  }
  return out;
}

//------------------------------------------------------------------------------
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  const auto Qkabcd = Qk_abcd(Fa, Fb, Fc, Fd, k);
  const auto Pkabcd = Pk_abcd(Fa, Fb, Fc, Fd, k);
  return (Qkabcd + Pkabcd);
}

//------------------------------------------------------------------------------
double Zk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k) {
  // Z^k_abcd = s ( Q^k_abcd + sum_l [k] 6j * Q^l_abdc)
  [[maybe_unused]] auto sp1 = IO::Profile::safeProfiler(__func__);
  const auto s = Angular::evenQ_2(Fa.twoj() + Fb.twoj() + 2) ? 1 : -1;
  return s * Wk_abcd(Fa, Fb, Fc, Fd, k);
}

//******************************************************************************
std::pair<int, int> k_minmax(const DiracSpinor &a, const DiracSpinor &b) {
  return k_minmax_tj(a.twoj(), b.twoj());
}
std::pair<int, int> k_minmax_tj(int tja, int tjb) {
  return std::make_pair(std::abs(tja - tjb) / 2, (tja + tjb) / 2);
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_Q(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {

  // Determine if K needs to be even/odd (parity selection rule)
  const auto k_even_ac = (a.l() + c.l()) % 2 == 0;
  const auto k_even_bd = (b.l() + d.l()) % 2 == 0;
  if (k_even_ac != k_even_bd) {
    // no K satisfies selection rule!
    return {1, 0};
  }

  // Find min/max k from triangle rule:
  const auto [l1, u1] = k_minmax(a, c);
  const auto [l2, u2] = k_minmax(b, d);
  auto min_k = std::max(l1, l2);
  auto max_k = std::min(u1, u2);

  // Adjust min/max k due to parity selectrion rule:
  // Allows one to safely use only evey second k to calculate Q
  if ((b.l() + d.l() + min_k) % 2 != 0) {
    ++min_k;
  }
  if ((b.l() + d.l() + max_k) % 2 != 0) {
    --max_k;
  }

  return {min_k, max_k};
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_P(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {
  // P^k_abcd = sum_l {a, c, k \\ b, d, l} * Q^l_abdc
  //  |b-d| <= k <=|b+d|
  //  |a-c| <= k <=|a+c|
  // min/max k Comes from 6j symbol ONLY
  const auto [l1, u1] = k_minmax(a, c);
  const auto [l2, u2] = k_minmax(b, d);
  return {std::max({l1, l2}), std::min({u1, u2})};
}
std::pair<int, int> k_minmax_P(int kappa_a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {
  // P^k_abcd = sum_l {a, c, k \\ b, d, l} * Q^l_abdc
  //  |b-d| <= k <=|b+d|
  //  |a-c| <= k <=|a+c|
  // min/max k Comes from 6j symbol ONLY
  const auto [l1, u1] = k_minmax_tj(Angular::twoj_k(kappa_a), c.twoj());
  const auto [l2, u2] = k_minmax(b, d);
  return {std::max({l1, l2}), std::min({u1, u2})};
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_W(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {
  // NOTE: Cannot safely k++2, since parity rules may be opposite for P and Q
  // parts!
  const auto [l1, u1] = k_minmax_Q(a, b, c, d);
  const auto [l2, u2] = k_minmax_P(a, b, c, d);
  // nb: min/max swapped, since W = Q+P, so only 1 needs to survive!
  return {std::min(l1, l2), std::max(u1, u2)};
}

} // namespace Coulomb
