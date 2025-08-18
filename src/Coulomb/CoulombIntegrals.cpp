#include "CoulombIntegrals.hpp"
#include "Angular/Angular.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace Coulomb {

//==============================================================================
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
  const auto &gr = Fa.grid(); // just save typing
  const auto du = gr.du();
  const auto num_points = gr.num_points();
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
    return (Fa.f(i) * Fb.f(i) + Fa.g(i) * Fb.g(i)) * w(i) * gr.drduor(i);
  };
  const auto &r = gr.r();

  double Ax = 0.0, Bx = 0.0;

  vabk[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {
    const auto rat = r[i - 1] / r[i];
    Ax = (Ax + ff(i - 1)) * (rat * powk(rat));
    vabk[i] = Ax * du;
  }

  const auto bmax = std::min({Fa.max_pt(), Fb.max_pt(), num_points});
  // nb bmax may be num_points
  const auto rbmax =
      bmax == num_points ? r.back() + gr.drdu().back() * du : r[bmax];
  Bx = Bx * powk(r[bmax - 1] / rbmax) + ff(bmax - 1);
  vabk[bmax - 1] += Bx * du;
  for (auto i = bmax - 1; i >= 1; --i) {
    Bx = Bx * powk(r[i - 1] / r[i]) + ff(i - 1);
    vabk[i - 1] += Bx * du;
  }

  for (std::size_t i = irmax; i < num_points; i++) {
    vabk[i] = 0.0;
  }
}

//------------------------------------------------------------------------------
// Used for Breit
template <int k, typename Function>
static inline void yk_ijk_gen_impl(const int l, const Function &ff,
                                   const Grid &gr, std::vector<double> &v0,
                                   std::vector<double> &vi,
                                   const std::size_t maxi)
// Calculalates genergic v^k_ab-like screening function.
// Stores "0" and "infinity" parts seperately.
//
// Stores in v0 and vi (in/out parameter, reference to whatever)
//
// r_min := min(r,r')
// rho(r') := fa(r')*fb(r') + ga(r')gb(r')
// v^k_ab(r) = Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
//           = Int_0^r [r'^k/r^(k+1)]*rho(r') dr'
//             + Int_r^inf [r^k/r'^(k+1)]*rho(r') dr'
//          := A(r)/r^(k+1) + B(r)*r^k
//           = v0(r)        + vi(r)
// A(r0)  = 0
// B(r0)  = Int_0^inf [r^k/r'^(k+1)]*rho(r') dr'
// A(r_n) = A(r_{n-1}) + (rho(r_{n-1})*r_{n-1}^k)*dr
// B(r_n) = A(r_{n-1}) + (rho(r_{n-1})/r_{n-1}^(k+1))*dr
// y^k_ab(rn) = A(rn)/rn^(k+1) + B(rn)*rn^k
//
// Also uses Quadrature integer rules! (Defined in NumCalc)
{
  // const auto &gr = Fa.grid(); // just save typing
  const auto du = gr.du();
  const auto num_points = gr.num_points();
  v0.resize(num_points); // for safety
  vi.resize(num_points); // for safety
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

  const auto &r = gr.r();

  double Ax = 0.0, Bx = 0.0;

  v0[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {
    const auto rat = r[i - 1] / r[i];
    Ax = (Ax + ff(i - 1) * w(i - 1) * gr.drduor(i - 1)) * (rat * powk(rat));
    v0[i] = Ax * du;
  }
  for (std::size_t i = irmax; i < num_points; i++) {
    v0[i] = 0.0;
  }

  // nb bmax may be num_points
  const auto bmax = irmax;
  for (std::size_t i = 0; i <= bmax; i++) {
    vi[i] = 0.0;
  }
  const auto rbmax =
      bmax == num_points ? r.back() + gr.drdu().back() * du : r[bmax];
  Bx = Bx * powk(r[bmax - 1] / rbmax) + ff(bmax - 1);
  vi[bmax - 1] += Bx * du;
  for (auto i = bmax - 1; i >= 1; --i) {
    Bx = Bx * powk(r[i - 1] / r[i]) + ff(i - 1) * w(i - 1) * gr.drduor(i - 1);
    vi[i - 1] = Bx * du;
  }
}

//------------------------------------------------------------------------------

// Used for frequency-dependent Breit
template <int k, typename Function>
static inline void
yk_ijk_gen_impl_freq(const int l, const Function &ff, const Grid &gr,
                     std::vector<double> &v0, std::vector<double> &vi,
                     const std::size_t maxi, const double w) {
  const auto du = gr.du();
  const auto num_points = gr.num_points();
  v0.resize(num_points); // for safety
  vi.resize(num_points); // for safety
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;

  // Quadrature integration weights:
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto &r = gr.r();

  auto jkwr_func = [](const int &l, const double &w, const double &r) {
    return SphericalBessel::exactGSL_JL_alt(l, w * r);
  };
  auto ykwr_func = [](const int &l, const double &w, const double &r) {
    return SphericalBessel::exactGSL_YL_alt(l, w * r);
  };

  // fills vector with spherical Bessel functions of first kind j_k(\alpha*\omega*r)
  //const auto jkwr = SphericalBessel::fillBesselVec_kr_alt(l, w, r);

  // fills vector with spherical Bessel functions of second kind y_k(\alpha*\omega*r)
  //const auto ykwr = SphericalBessel::fillSecondBesselVec_kr_alt(l, w, r);

  // values to keep running track of numerical integrals
  double Ax = 0.0;
  double Bx = 0.0;

  // performs numerical integral from r'=0 to r'=r using single for loop trick
  v0[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {
    Ax = Ax + jkwr_func(k, w, r[i - 1]) * ff(i - 1) * weights(i - 1) *
                  gr.drdu(i - 1);
    v0[i] = -w * (2 * k + 1.0) * ykwr_func(k, w, r[i]) * Ax * du;
  }
  for (std::size_t i = irmax; i < num_points; i++) {
    v0[i] = 0.0;
  }

  // nb bmax may be num_points
  const auto bmax = irmax;
  for (std::size_t i = 0; i <= bmax; i++) {
    vi[i] = 0.0;
  }
  const auto rbmax =
      bmax == num_points ? r.back() + gr.drdu().back() * du : r[bmax];

  Bx = Bx + ykwr_func(k, w, r[bmax - 1]) * ff(bmax - 1) * gr.drdu(bmax - 1);

  vi[bmax - 1] = -w * (2.0 * k + 1.0) * jkwr_func(k, w, r[bmax - 1]) * Bx * du;

  // loop for the integral that goes from r'=r up to r<infinity
  // in this loop I have to shift the index up by 1 so that i goes down to i>=1 rather than i>=0, since it freaks out when I do the latter, for some reason
  for (auto i = bmax - 1; i >= 1; i--) {
    Bx = Bx + ykwr_func(k, w, r[i - 1]) * ff(i - 1) * weights(i - 1) *
                  gr.drdu(i - 1);
    vi[i - 1] = -w * (2 * k + 1.0) * jkwr_func(k, w, r[i - 1]) * Bx * du;
  }
}

//---------------------------------------------------------------------------------

// function for calculating the Breit frequency-dependent screening factor vkabcd(w) - stores first term v1 and second term v2 separately
template <int k>
static inline void
vkabcd_freqw(const int l, const std::vector<double> &Pkbd,
             const std::vector<double> &Qkbd, const Grid &gr,
             std::vector<double> &v1, std::vector<double> &v2,
             std::vector<double> &v3, std::vector<double> &v4,
             const std::size_t maxi, const double w) {

  const auto du = gr.du();
  const auto num_points = gr.num_points();
  v1.resize(num_points); // for safety
  v2.resize(num_points); // for safety
  v3.resize(num_points); // for safety
  v4.resize(num_points); // for safety
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
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto &r = gr.r();
  // // fills vector with spherical Bessel functions of first kind j_k(\alpha*\omega*r)
  // const auto jkpluswr = SphericalBessel::fillBesselVec_kr_alt(l + 1, w, r);
  // const auto jkminuswr = SphericalBessel::fillBesselVec_kr_alt(l - 1, w, r);

  // // fills vector with spherical Bessel functions of second kind y_k(\alpha*\omega*r)
  // const auto ykpluswr =
  //     SphericalBessel::fillSecondBesselVec_kr_alt(l + 1, w, r);
  // const auto ykminuswr =
  //     SphericalBessel::fillSecondBesselVec_kr_alt(l - 1, w, r);

  auto jkwr = [](const int &l, const double &w, const double &r) {
    return SphericalBessel::exactGSL_JL_alt(l, w * r);
  };
  auto ykwr = [](const int &l, const double &w, const double &r) {
    return SphericalBessel::exactGSL_YL_alt(l, w * r);
  };

  // values to keep running track of numerical integrals
  double A1 = 0.0;
  double A2 = 0.0;
  double A3 = 0.0;

  // performs numerical integrals for v1 and v2
  v1[0] = 0.0;
  v2[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {
    // first term in first part of v
    //A1 = A1 + jkminuswr[i - 1] * Pkbd[i - 1] * weights(i - 1) * gr.drdu(i - 1);
    A1 = A1 + jkwr(k - 1, w, r[i - 1]) * Pkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);
    //A1 = 0.0;

    // second term in first part of v
    double ratio = r[i - 1] / r[i];
    // A2 = (A2 + Pkbd[i - 1] * weights(i - 1) * gr.drduor(i - 1) /
    //                (r[i - 1] * r[i - 1])) *
    //      pow(ratio, k + 2);

    A2 = A2 + pow(r[i - 1], double(k - 1)) * Pkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);
    //A2 = 0.0;

    //v1[i] = w * ykpluswr[i] * A1 + ((2.0 * k + 1.0) / (w * w)) * A2;
    //v1[i] = w * ykpluswr[i] * A1 + ((2.0 * k + 1.0) / (w * w)) * (1.0 / pow(r[i], double(k + 2))) * A2;
    v1[i] = w * ykwr(k + 1, w, r[i]) * A1 +
            ((2.0 * k + 1.0) / (w * w)) * (1.0 / pow(r[i], double(k + 2))) * A2;
    v1[i] = -2.0 * v1[i] * du;

    // integral for second term in v
    //A3 = A3 + jkpluswr[i - 1] * Qkbd[i - 1] * weights(i - 1) * gr.drdu(i - 1);
    A3 = A3 + jkwr(k + 1, w, r[i - 1]) * Qkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);
    //v2[i] = -2.0 * w * ykminuswr[i] * A3 * du;
    v2[i] = -2.0 * w * ykwr(k - 1, w, r[i]) * A3 * du;
  }

  for (std::size_t i = irmax; i < num_points; i++) {
    v1[i] = 0.0;
    v2[i] = 0.0;
  }

  double B1 = 0.0;
  double B2 = 0.0;
  double B3 = 0.0;

  // nb bmax may be num_points
  const auto bmax = irmax;
  for (std::size_t i = 0; i <= bmax; i++) {
    v3[i] = 0.0;
    v4[i] = 0.0;
  }
  const auto rbmax =
      bmax == num_points ? r.back() + gr.drdu().back() * du : r[bmax];

  // calculating screening functions at end point
  //B1 = B1 + ykpluswr[bmax - 1] * Qkbd[bmax - 1] * gr.drdu(bmax - 1);
  B1 = B1 + ykwr(k + 1, w, r[bmax - 1]) * Qkbd[bmax - 1] * gr.drdu(bmax - 1);
  B2 = B2 + ((2.0 * k + 1.0) / (w * w)) * (1.0 / pow(r[bmax - 1], k + 2)) *
                Qkbd[bmax - 1] * gr.drdu(bmax - 1);

  //v3[bmax - 1] = -2.0 * w * jkminuswr[bmax - 1] * B1 * du + pow(r[bmax - 1], k - 1) * B2 * du;
  v3[bmax - 1] = -2.0 * w * jkwr(k - 1, w, r[bmax - 1]) * B1 * du +
                 pow(r[bmax - 1], k - 1) * B2 * du;

  //B3 = B3 + ykminuswr[bmax - 1] * Pkbd[bmax - 1] * gr.drdu(bmax - 1);
  B3 = B3 + ykwr(k - 1, w, r[bmax - 1]) * Pkbd[bmax - 1] * gr.drdu(bmax - 1);

  //v4[bmax - 1] = -2.0 * w * jkpluswr[bmax - 1] * B3 * du;
  v4[bmax - 1] = -2.0 * w * jkwr(k + 1, w, r[bmax - 1]) * B3 * du;

  for (auto i = bmax - 1; i >= 1; i--) {
    // j_{k-1}(wr1)y_{k+1}(wr2) term in integral that is integrated from r1 <= r2 <=infty
    //B1 = B1 + ykpluswr[i - 1] * Qkbd[i - 1] * weights(i - 1) * gr.drdu(i - 1);
    B1 = B1 + ykwr(k + 1, w, r[i - 1]) * Qkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);

    //B1 = 0.0;

    // r2^{k-1}/r1^{k+2} term in integral that is integrated from r1 <= r2 <=infty
    B2 = B2 + ((2.0 * k + 1.0) / (w * w)) * (1.0 / pow(r[i - 1], k + 2)) *
                  Qkbd[i - 1] * weights(i - 1) * gr.drdu(i - 1);
    //B2 = 0.0;

    //v3[i - 1] = -2.0 * w * jkminuswr[i - 1] * B1 * du - 2.0 * pow(r[i - 1], k - 1) * B2 * du;
    v3[i - 1] = -2.0 * w * jkwr(k - 1, w, r[i - 1]) * B1 * du -
                2.0 * pow(r[i - 1], k - 1) * B2 * du;

    // j_{k+1}(wr1)y_{k-1}(wr2) term in integral that is integrated from r1 <= r2 <=infty
    //B3 = B3 + ykminuswr[i - 1] * Pkbd[i - 1] * weights(i - 1) * gr.drdu(i - 1);
    B3 = B3 + ykwr(k - 1, w, r[i - 1]) * Pkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);
    //v4[i - 1] = -2.0 * w * jkpluswr[i - 1] * B3 * du;
    v4[i - 1] = -2.0 * w * jkwr(k + 1, w, r[i - 1]) * B3 * du;
  }
}

//------------------------------------------------------------------------------
std::vector<double> yk_ab(const int k, const DiracSpinor &Fa,
                          const DiracSpinor &Fb, const std::size_t maxi) {
  std::vector<double> ykab; //
  yk_ab(k, Fa, Fb, ykab, maxi);
  return ykab;
}

//------------------------------------------------------------------------------
void yk_ab(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
           std::vector<double> &vabk, const std::size_t maxi) {

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

//==============================================================================
template <int k, int pm>
static inline void Breit_abk_impl(const int l, const DiracSpinor &Fa,
                                  const DiracSpinor &Fb, //
                                  std::vector<double> &b0,
                                  std::vector<double> &binf,
                                  const std::size_t maxi) {
  static_assert(pm == 1 || pm == -1,
                "Breit_abk_impl must be called with pm=+/-1 only\n");
  auto fgfg = [&](std::size_t i) {
    if constexpr (pm == -1)
      return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
    else
      return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };
  yk_ijk_gen_impl<k>(l, fgfg, Fa.grid(), b0, binf, maxi);
  return;
}
//------------------------------------------------------------------------------
void bk_ab(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
           std::vector<double> &b0, std::vector<double> &binf,
           std::size_t maxi) {

  auto fgfg = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  if (k == 0)
    yk_ijk_gen_impl<0>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 1)
    yk_ijk_gen_impl<1>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 2)
    yk_ijk_gen_impl<2>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 3)
    yk_ijk_gen_impl<3>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 4)
    yk_ijk_gen_impl<4>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 5)
    yk_ijk_gen_impl<5>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 6)
    yk_ijk_gen_impl<6>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 7)
    yk_ijk_gen_impl<7>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else if (k == 8)
    yk_ijk_gen_impl<8>(k, fgfg, Fa.grid(), b0, binf, maxi);
  else
    yk_ijk_gen_impl<-1>(k, fgfg, Fa.grid(), b0, binf, maxi);
}
//------------------------------------------------------------------------------
void gk_ab(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
           std::vector<double> &g0, std::vector<double> &ginf,
           const std::size_t maxi) {

  auto fgfg = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  if (k == 0)
    yk_ijk_gen_impl<0>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 1)
    yk_ijk_gen_impl<1>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 2)
    yk_ijk_gen_impl<2>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 3)
    yk_ijk_gen_impl<3>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 4)
    yk_ijk_gen_impl<4>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 5)
    yk_ijk_gen_impl<5>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 6)
    yk_ijk_gen_impl<6>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 7)
    yk_ijk_gen_impl<7>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else if (k == 8)
    yk_ijk_gen_impl<8>(k, fgfg, Fa.grid(), g0, ginf, maxi);
  else
    yk_ijk_gen_impl<-1>(k, fgfg, Fa.grid(), g0, ginf, maxi);
}

//==============================================================================

// frequency-dependent Breit integrals which i call g^k_{ab}(r,w) and h^k_{ab}(r,w) for X_{ij} and Y_{ij}, respectively
// also need to define v^{k,1}_{ab} and v^k{k,2}_{ab} for the two respective terms in v^k_{ab} which will be multiplying different parts of the
void gk_ab_freqw(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 std::vector<double> &g0, std::vector<double> &ginf,
                 const std::size_t maxi, const double w) {

  auto fgfg = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  if (k == 0)
    yk_ijk_gen_impl_freq<0>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 1)
    yk_ijk_gen_impl_freq<1>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 2)
    yk_ijk_gen_impl_freq<2>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 3)
    yk_ijk_gen_impl_freq<3>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 4)
    yk_ijk_gen_impl_freq<4>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 5)
    yk_ijk_gen_impl_freq<5>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 6)
    yk_ijk_gen_impl_freq<6>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 7)
    yk_ijk_gen_impl_freq<7>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else if (k == 8)
    yk_ijk_gen_impl_freq<8>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
  else
    yk_ijk_gen_impl_freq<-1>(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
}

void hk_ab_freqw(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 std::vector<double> &b0, std::vector<double> &binf,
                 std::size_t maxi, const double w) {

  auto fgfg = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  if (k == 0)
    yk_ijk_gen_impl_freq<0>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 1)
    yk_ijk_gen_impl_freq<1>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 2)
    yk_ijk_gen_impl_freq<2>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 3)
    yk_ijk_gen_impl_freq<3>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 4)
    yk_ijk_gen_impl_freq<4>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 5)
    yk_ijk_gen_impl_freq<5>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 6)
    yk_ijk_gen_impl_freq<6>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 7)
    yk_ijk_gen_impl_freq<7>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else if (k == 8)
    yk_ijk_gen_impl_freq<8>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
  else
    yk_ijk_gen_impl_freq<-1>(k, fgfg, Fa.grid(), b0, binf, maxi, w);
}

void vk_ab_freqw(const int k, const DiracSpinor &Fi, const DiracSpinor &Fj,
                 const Grid &gr, std::vector<double> &v1,
                 std::vector<double> &v2, std::vector<double> &v3,
                 std::vector<double> &v4, std::size_t maxi, const double w) {

  auto Xij = [&](std::size_t i) {
    return (Fi.f(i) * Fj.g(i) + Fi.g(i) * Fj.f(i));
  };

  auto Yij = [&](std::size_t i) {
    return (Fi.f(i) * Fj.g(i) - Fi.g(i) * Fj.f(i));
  };

  std::vector<double> Pkbd(gr.size(), 0.0);
  std::vector<double> Qkbd(gr.size(), 0.0);

  for (int i = 0; i < gr.size(); i++) {
    Pkbd[i] = ((Fi.kappa() - Fj.kappa()) / k) * Xij(i) - Yij(i);
    Qkbd[i] = ((Fi.kappa() - Fj.kappa()) / (k + 1.0)) * Xij(i) + Yij(i);
  }

  if (k == 0)
    vkabcd_freqw<0>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 1)
    vkabcd_freqw<1>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 2)
    vkabcd_freqw<2>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 3)
    vkabcd_freqw<3>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 4)
    vkabcd_freqw<4>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 5)
    vkabcd_freqw<5>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 6)
    vkabcd_freqw<6>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 7)
    vkabcd_freqw<7>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else if (k == 8)
    vkabcd_freqw<8>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
  else
    vkabcd_freqw<-1>(k, Pkbd, Qkbd, gr, v1, v2, v3, v4, maxi, w);
}

//==============================================================================
//==============================================================================

//==============================================================================
double Rk_abcd(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd) {
  // note: Returns non-zero value, even if angular part would be zero
  const auto yk_bd = yk_ab(k, Fb, Fd, std::min(Fa.max_pt(), Fc.max_pt()));
  return Rk_abcd(Fa, Fc, yk_bd);
}
//------------------------------------------------------------------------------
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fc,
               const std::vector<double> &yk_bd) {
  const auto &drdu = Fa.grid().drdu();
  const auto i0 = std::max(Fa.min_pt(), Fc.min_pt());
  const auto imax = std::min(Fa.max_pt(), Fc.max_pt());
  const auto Rff =
      NumCalc::integrate(1.0, i0, imax, Fa.f(), Fc.f(), yk_bd, drdu);
  const auto Rgg =
      NumCalc::integrate(1.0, i0, imax, Fa.g(), Fc.g(), yk_bd, drdu);
  return (Rff + Rgg) * Fa.grid().du();
}

//==============================================================================
DiracSpinor Rkv_bcd(const int k, const int kappa_a, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) {

  return Rkv_bcd(kappa_a, Fc, yk_ab(k, Fb, Fd, Fc.max_pt()));
}
//------------------------------------------------------------------------------
DiracSpinor Rkv_bcd(const int kappa_a, const DiracSpinor &Fc,
                    const std::vector<double> &ykbd) {

  auto out = DiracSpinor(0, kappa_a, Fc.grid_sptr());
  out.min_pt() = Fc.min_pt();
  out.max_pt() = Fc.max_pt();
  out.f() = qip::multiply(Fc.f(), ykbd);
  out.g() = qip::multiply(Fc.g(), ykbd);
  return out;
}
//------------------------------------------------------------------------------
void Rkv_bcd(DiracSpinor *const Rkv, const DiracSpinor &Fc,
             const std::vector<double> &ykbd) {

  Rkv->min_pt() = Fc.min_pt();
  Rkv->max_pt() = Fc.max_pt();
  for (auto i = 0ul; i < Rkv->min_pt(); ++i) {
    Rkv->f(i) = 0.0;
    Rkv->g(i) = 0.0;
  }
  for (auto i = Rkv->min_pt(); i < Rkv->max_pt(); ++i) {
    Rkv->f(i) = Fc.f(i) * ykbd[i];
    Rkv->g(i) = Fc.g(i) * ykbd[i];
  }
  for (auto i = Rkv->max_pt(); i < Rkv->grid().num_points(); ++i) {
    Rkv->f(i) = 0.0;
    Rkv->g(i) = 0.0;
  }
}

//==============================================================================
double Qk_abcd(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd) {

  // const auto tCac = Angular::tildeCk_kk(k, Fa.kappa(), Fc.kappa());
  // if (Angular::zeroQ(tCac))
  //   return 0.0;
  // const auto tCbd = Angular::tildeCk_kk(k, Fb.kappa(), Fd.kappa());
  // if (Angular::zeroQ(tCbd))
  //   return 0.0;
  const auto tCac = Angular::Ck_kk(k, Fa.kappa(), Fc.kappa());
  if (tCac == 0.0)
    return 0.0;
  const auto tCbd = Angular::Ck_kk(k, Fb.kappa(), Fd.kappa());
  if (tCbd == 0.0)
    return 0.0;
  const auto Rkabcd = Rk_abcd(k, Fa, Fb, Fc, Fd);
  // const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  const auto m1tk2 = Angular::neg1pow_2(2 * k + Fa.twoj() - Fb.twoj());
  return m1tk2 * tCac * tCbd * Rkabcd;
}

//------------------------------------------------------------------------------
bool Qk_abcd_SR(int k, int ka, int kb, int kc, int kd) {
  return Angular::Ck_kk_SR(k, ka, kc) && Angular::Ck_kk_SR(k, kb, kd);
}

bool Qk_abcd_SR(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d) {
  return Qk_abcd_SR(k, a.kappa(), b.kappa(), c.kappa(), d.kappa());
}
//------------------------------------------------------------------------------
bool Pk_abcd_SR(int k, int ka, int kb, int kc, int kd) {
  return Angular::Ck_kk_SR(k, ka, kd) && Angular::Ck_kk_SR(k, kb, kc);
}

bool Pk_abcd_SR(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d) {
  return Pk_abcd_SR(k, a.kappa(), b.kappa(), d.kappa(), c.kappa());
}

//------------------------------------------------------------------------------
DiracSpinor Qkv_bcd(const int k, const int kappa_a, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) {
  assert(k >= 0 && "Check k and kappa order");

  const auto tCac = Angular::tildeCk_kk(k, kappa_a, Fc.kappa());
  const auto tCbd = Angular::tildeCk_kk(k, Fb.kappa(), Fd.kappa());
  if (Angular::zeroQ(tCbd) || Angular::zeroQ(tCac))
    return DiracSpinor(0, kappa_a, Fc.grid_sptr());
  const auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return (m1tk * tCac * tCbd) * Rkv_bcd(k, kappa_a, Fb, Fc, Fd);
}

//------------------------------------------------------------------------------
// void Qkv_bcd(DiracSpinor *const Qkv, const DiracSpinor &Fb,
//              const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
//              const std::vector<double> &ykbd, const Angular::CkTable &Ck) {
//   const auto tCac = Ck.get_tildeCkab(k, Qkv->kappa(), Fc.kappa());
//   const auto tCbd = Ck.get_tildeCkab(k, Fb.kappa(), Fd.kappa());
//   const auto tCC = tCbd * tCac;
//   if (tCC == 0.0) {
//     Qkv->scale(0.0);
//     return;
//   }
//   Rkv_bcd(Qkv, Fc, ykbd);
//   const auto m1tk = Angular::evenQ(k) ? 1 : -1;
//   Qkv->scale(m1tk * tCC);
//   return;
// }

//==============================================================================
double Pk_abcd(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc

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
    const auto Qlabdc = Qk_abcd(tl / 2, Fa, Fb, Fd, Fc);
    sum += sixj * Qlabdc;
  }
  return (tkp1 * sum);
}

//------------------------------------------------------------------------------
DiracSpinor Pkv_bcd(const int k, int kappa_a, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) {
  assert(k >= 0 && "Check k and kappa order");
  auto out = DiracSpinor(0, kappa_a, Fc.grid_sptr());
  const auto tkp1 = 2 * k + 1;
  const auto tja = Angular::twoj_k(kappa_a);
  const auto min_twol =
      std::max(std::abs(Fd.twoj() - tja), std::abs(Fc.twoj() - Fb.twoj()));
  const auto max_twol = std::min(Fd.twoj() + tja, Fc.twoj() + Fb.twoj());
  for (int tl = min_twol; tl <= max_twol; tl += 2) {
    if (!Angular::Ck_kk_SR(tl / 2, Fb.kappa(), Fc.kappa()) ||
        !Angular::Ck_kk_SR(tl / 2, kappa_a, Fd.kappa()))
      continue;
    const auto sixj =
        Angular::sixj_2(Fc.twoj(), tja, 2 * k, Fd.twoj(), Fb.twoj(), tl);
    if (sixj == 0)
      continue;
    out += sixj * Qkv_bcd(tl / 2, kappa_a, Fb, Fd, Fc);
  }
  return (tkp1 * out);
}

//------------------------------------------------------------------------------
DiracSpinor Wkv_bcd(const int k, int kappa_v, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd) {
  assert(k >= 0 && "Check k and kappa order");
  auto out = Pkv_bcd(k, kappa_v, Fb, Fc, Fd);
  if (Angular::Ck_kk_SR(k, kappa_v, Fc.kappa()) &&
      Angular::Ck_kk_SR(k, Fb.kappa(), Fd.kappa())) {
    out += Qkv_bcd(k, kappa_v, Fb, Fc, Fd);
  }
  return out;
}

//------------------------------------------------------------------------------
double Wk_abcd(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd) {
  // W^k_abcd = Q^k_abcd + sum_l [k] 6j * Q^l_abdc
  return Qk_abcd(k, Fa, Fb, Fc, Fd) + Pk_abcd(k, Fa, Fb, Fc, Fd);
}

//==============================================================================
double g_abcd(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
              const DiracSpinor &d, int tma, int tmb, int tmc, int tmd) {

  if (tmc - tma != tmb - tmd)
    return 0.0;
  const int twoq = tmc - tma;

  //
  double g = 0.0;
  const auto [k0, ki] = k_minmax_Q(a, b, c, d);
  for (int k = k0; k <= ki; k += 2) {
    // for (int q = -k; q <= k; q++) {
    //   int twoq = 2 * q;
    if (std::abs(twoq) > 2 * k)
      continue;
    const auto A =
        Angular::neg1pow_2(twoq) *
        Angular::Ck_kk_mmq(k, a.kappa(), c.kappa(), tma, tmc, -twoq) *
        Angular::Ck_kk_mmq(k, b.kappa(), d.kappa(), tmb, tmd, twoq);
    if (A != 0.0) {
      g += A * Rk_abcd(k, a, b, c, d);
    }
  }
  // }
  return g;
}

//==============================================================================
std::pair<int, int> k_minmax_Ck(const DiracSpinor &a, const DiracSpinor &b) {
  // return k_minmax_tj(a.twoj(), b.twoj());
  auto minmax = k_minmax_tj(a.twoj(), b.twoj());
  auto &min_k = minmax.first;
  auto &max_k = minmax.second;
  if ((a.l() + b.l() + min_k) % 2 != 0) {
    ++min_k;
  }
  if ((a.l() + b.l() + max_k) % 2 != 0) {
    --max_k;
  }
  return minmax;
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_Ck(int kappa_a, int kappa_b) {

  auto minmax = k_minmax_tj(Angular::twoj_k(kappa_a), Angular::twoj_k(kappa_b));
  auto &min_k = minmax.first;
  auto &max_k = minmax.second;
  const auto la = Angular::l_k(kappa_a);
  const auto lb = Angular::l_k(kappa_b);
  if ((la + lb + min_k) % 2 != 0) {
    ++min_k;
  }
  if ((la + lb + max_k) % 2 != 0) {
    --max_k;
  }
  return minmax;
}

//------------------------------------------------------------------------------
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
  const auto [l1, u1] = k_minmax_Ck(a, c);
  const auto [l2, u2] = k_minmax_Ck(b, d);
  return {std::max(l1, l2), std::min(u1, u2)};
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_Q(int kap_a, int kap_b, int kap_c, int kap_d) {

  // Determine if K needs to be even/odd (parity selection rule)
  const auto k_even_ac = (Angular::l_k(kap_a) + Angular::l_k(kap_c)) % 2 == 0;
  const auto k_even_bd = (Angular::l_k(kap_b) + Angular::l_k(kap_d)) % 2 == 0;
  if (k_even_ac != k_even_bd) {
    // no K satisfies selection rule!
    return {1, 0};
  }

  // Find min/max k from triangle rule:
  const auto [l1, u1] = k_minmax_Ck(kap_a, kap_c);
  const auto [l2, u2] = k_minmax_Ck(kap_b, kap_d);
  return {std::max(l1, l2), std::min(u1, u2)};
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_P(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {
  // P^k_abcd = sum_l {a, c, k \\ b, d, l} * Q^l_abdc

  // From the 6j part:
  // |b-d| <= k <=|b+d|
  // |a-c| <= k <=|a+c|
  const auto [lk1, uk1] = k_minmax_tj(a.twoj(), c.twoj());
  const auto [lk2, uk2] = k_minmax_tj(b.twoj(), d.twoj());

  // From the Qk part: nb: this is for l, internal sum!
  const auto [lq0, uq0] = k_minmax_Q(a, b, d, c);
  // but if there are _no_ Q^l, then it is zero
  if (lq0 > uq0)
    return {1, 0};

  return {std::max({lk1, lk2}), std::min({uk1, uk2})};
}

//------------------------------------------------------------------------------
std::pair<int, int> k_minmax_W(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {
  const auto [l1, u1] = k_minmax_Q(a, b, c, d);
  const auto [l2, u2] = k_minmax_P(a, b, c, d);
  // if (l1 > u1 && l2 > u2)
  //   return {1, 0};
  if (l1 > u1)
    return {l2, u2};
  if (l2 > u2)
    return {l1, u1};
  // nb: min/max swapped, since W = Q+P, so only 1 needs to survive!
  return {std::min(l1, l2), std::max(u1, u2)};
}

} // namespace Coulomb
