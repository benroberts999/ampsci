#include "CoulombBreit.hpp"
#include "Angular/include.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <vector>

namespace Coulomb {

//==============================================================================
// Used for Breit
template <int KT, typename Function>
static inline void yk_ijk_gen_impl(const int k, const Function &ff,
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
  const auto powk = [k]() {
    if constexpr (KT < 0) {
      return [k](double x) { return std::pow(x, k); };
    } else {
      (void)k;
      return qip::pow<KT, double>;
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
  Bx = Bx * powk(r[bmax - 1] / rbmax) +
       ff(bmax - 1) * w(bmax - 1) * gr.drduor(bmax - 1);
  vi[bmax - 1] += Bx * du;
  for (auto i = bmax - 1; i >= 1; --i) {
    Bx = Bx * powk(r[i - 1] / r[i]) + ff(i - 1) * w(i - 1) * gr.drduor(i - 1);
    vi[i - 1] = Bx * du;
  }
}

//==============================================================================

// Implements the freq-dep Breit screening function (arXiv:2602.17129).
// Applies the kernel replacement:
//   r_<^k / r_>^{k+1}  ->  -w*(2k+1)*j_k(w*r_<)*y_k(w*r_>)
// Used by gk_ab_freqw (X_ab density) and hk_ab_freqw (Y_ab density).
template <typename Function>
static inline void
yk_ijk_gen_impl_freq(const int k, const Function &ff, const Grid &gr,
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
  const auto jL = SphericalBessel::exactGSL_JL_alt<double>;
  const auto yL = SphericalBessel::exactGSL_YL_alt<double>;

  // For small w the Bessel kernel -w(2k+1) j_k(wr_<) y_k(wr_>) causes
  // bad cancellation (y_k diverges, j_k vanishes). Use O(w^2) Taylor
  // expansion of the kernel instead:
  //   -w(2k+1) j_k(wr_<) y_k(wr_>) = r_<^k/r_>^{k+1}
  //     * [1 - w^2*r_<^2/(2(2k+3)) + w^2*r_>^2/(2(2k-1))] + O(w^4)
  // This requires two running sums each for v0 and vi.
  const bool cut_off = w <= PhysConst::alpha;
  const double w2 = w * w;
  const double c_lo = 1.0 / (2.0 * (2 * k + 3)); // coefficient of r_<^2 term
  const double c_hi = 1.0 / (2.0 * (2 * k - 1)); // coefficient of r_>^2 term

  // Ax = A_k(r)/r^{k+1}  where A_k(r) = int_0^r r'^k ff dr'  (static part)
  // Cx = A_{k+2}(r)/r^{k+1}                                   (w^2 correction)
  // Bx = r^k * int_r^inf ff/r'^{k+1} dr'                      (static part)
  // Dx = r^k * int_r^inf ff*r'^{1-k} dr' = r^k * B_{k-1}(r)  (w^2 correction)
  double Ax = 0.0, Cx = 0.0;
  double Bx = 0.0, Dx = 0.0;

  // performs numerical integral from r'=0 to r'=r using single for loop trick
  v0[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {
    if (cut_off) {
      const auto rat = r[i - 1] / r[i];
      const auto ratp1 = rat * std::pow(rat, k); // (r_{i-1}/r_i)^{k+1}
      Ax = (Ax + ff(i - 1) * weights(i - 1) * gr.drduor(i - 1)) * ratp1;
      Cx = (Cx + ff(i - 1) * weights(i - 1) * r[i - 1] * r[i - 1] *
                   gr.drduor(i - 1)) *
           ratp1;
      v0[i] = (Ax * (1.0 + w2 * r[i] * r[i] * c_hi) - w2 * c_lo * Cx) * du;
    } else {
      Ax =
        Ax + jL(k, w * r[i - 1]) * ff(i - 1) * weights(i - 1) * gr.drdu(i - 1);
      v0[i] = -w * (2 * k + 1.0) * yL(k, w * r[i]) * Ax * du;
    }
  }
  for (std::size_t i = irmax; i < num_points; i++) {
    v0[i] = 0.0;
  }

  // nb bmax may be num_points
  const auto bmax = irmax;
  for (std::size_t i = 0; i <= bmax; i++) {
    vi[i] = 0.0;
  }

  if (cut_off) {
    Bx = ff(bmax - 1) * weights(bmax - 1) * gr.drduor(bmax - 1);
    Dx = ff(bmax - 1) * weights(bmax - 1) * r[bmax - 1] * r[bmax - 1] *
         gr.drduor(bmax - 1);
    vi[bmax - 1] =
      (Bx * (1.0 - w2 * r[bmax - 1] * r[bmax - 1] * c_lo) + w2 * c_hi * Dx) *
      du;
  } else {
    Bx = Bx + yL(k, w * r[bmax - 1]) * ff(bmax - 1) * weights(bmax - 1) *
                gr.drdu(bmax - 1);
    vi[bmax - 1] = -w * (2.0 * k + 1.0) * jL(k, w * r[bmax - 1]) * Bx * du;
  }

  // loop for the integral that goes from r'=r up to r<infinity
  // in this loop I have to shift the index up by 1 so that i goes down to i>=1
  // rather than i>=0 (unsigned integer)
  for (auto i = bmax - 1; i >= 1; i--) {
    if (cut_off) {
      const auto rat = r[i - 1] / r[i];
      const auto ratk = std::pow(rat, k); // (r_{i-1}/r_i)^k
      Bx = Bx * ratk + ff(i - 1) * weights(i - 1) * gr.drduor(i - 1);
      Dx = Dx * ratk +
           ff(i - 1) * weights(i - 1) * r[i - 1] * r[i - 1] * gr.drduor(i - 1);
      vi[i - 1] =
        (Bx * (1.0 - w2 * r[i - 1] * r[i - 1] * c_lo) + w2 * c_hi * Dx) * du;
    } else {
      Bx =
        Bx + yL(k, w * r[i - 1]) * ff(i - 1) * weights(i - 1) * gr.drdu(i - 1);
      vi[i - 1] = -w * (2 * k + 1.0) * jL(k, w * r[i - 1]) * Bx * du;
    }
  }
}

//---------------------------------------------------------------------------------

// Computes the freq-dep v^k_abcd screening factor (arXiv:2602.17129).
// Pkbd = P^k_ij, Qkbd = Q^k_ij.
// v1,v2: (0,r) parts; v3,v4: (r,inf) parts; see vk_ab_freqw docs for full breakdown.
template <int KT>
static inline void
vkabcd_freqw(const int k, const std::vector<double> &Pkbd,
             const std::vector<double> &Qkbd, const Grid &gr,
             std::vector<double> &v1, std::vector<double> &v2,
             std::vector<double> &v3, std::vector<double> &v4,
             const std::size_t maxi, const double w) {

  bool cut_off = w <= PhysConst::alpha;
  // bool cut_off = true;
  // bool cut_off = false;

  const auto du = gr.du();
  const auto num_points = gr.num_points();
  v1.resize(num_points); // for safety
  v2.resize(num_points); // for safety
  v3.resize(num_points); // for safety
  v4.resize(num_points); // for safety
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;

  // faster method to calculate r^k
  const auto powk = [k]() {
    if constexpr (KT < 0) {
      return [k](double x) { return std::pow(x, k); };
    } else {
      (void)k;
      return qip::pow<KT, double>;
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
  const auto jL = SphericalBessel::exactGSL_JL_alt<double>;
  const auto yL = SphericalBessel::exactGSL_YL_alt<double>;

  // values to keep running track of numerical integrals
  double A1 = 0.0;
  double A2 = 0.0;
  double C1 = 0.0;
  double A3 = 0.0;

  const double wsd2 = w * w / 2.0;
  const double wcubed = w * w * w;
  const double odw2 = 1.0 / (w * w);

  // performs numerical integrals for v1 and v2
  v1[0] = 0.0;
  v2[0] = 0.0;
  for (std::size_t i = 1; i < irmax; ++i) {

    double ratio = r[i - 1] / r[i];

    //! This term numerically unstable for small omega
    if (cut_off) {
      A1 = (A1 + Pkbd[i - 1] * weights(i - 1) * gr.drduor(i - 1)) * powk(ratio);
      A2 = (A2 + Pkbd[i - 1] * weights(i - 1) * gr.drduor(i - 1)) *
           powk(ratio) * ratio * ratio;
      C1 = (C1 + Pkbd[i - 1] * weights(i - 1) * r[i - 1] * gr.drdu(i - 1)) *
           powk(ratio);
      v1[i] = (A1 - A2 + wsd2 * C1) * du;
    } else {
      // first term in first part of v
      A1 = A1 + jL(k - 1, w * r[i - 1]) * Pkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);

      // second term in first part of v
      A2 = (A2 + Pkbd[i - 1] * weights(i - 1) * gr.drduor(i - 1) * 1.0 /
                   (r[i - 1] * r[i - 1])) *
           powk(ratio) * ratio * ratio;

      v1[i] = wcubed * yL(k + 1, w * r[i]) * A1 + (2.0 * k + 1.0) * A2;
      v1[i] = -2.0 * odw2 * v1[i] * du;
    }

    // integral for second term in v
    A3 = A3 + jL(k + 1, w * r[i - 1]) * Qkbd[i - 1] * weights(i - 1) *
                gr.drdu(i - 1);
    v2[i] = -2.0 * w * yL(k - 1, w * r[i]) * A3 * du;
  }

  for (std::size_t i = irmax; i < num_points; i++) {
    v1[i] = 0.0;
    v2[i] = 0.0;
  }

  double B1 = 0.0;
  double B2 = 0.0;
  double D1 = 0.0;
  double B3 = 0.0;

  // nb bmax may be num_points
  const auto bmax = irmax;
  for (std::size_t i = 0; i <= bmax; i++) {
    v3[i] = 0.0;
    v4[i] = 0.0;
  }
  // const auto rbmax =
  //   bmax == num_points ? r.back() + gr.drdu().back() * du : r[bmax];

  // calculating screening functions at end point
  if (cut_off) {
    // B1 = Qkbd[bmax - 1] * weights(bmax - 1) * gr.drduor(bmax - 1);
    // B2 = Qkbd[bmax - 1] * weights(bmax - 1) * gr.drduor(bmax - 1);
    B1 = Qkbd[bmax - 1];
    B2 = Qkbd[bmax - 1];
    D1 = Qkbd[bmax - 1] * weights(bmax - 1) * r[bmax - 1] * gr.drdu(bmax - 1);
    v3[bmax - 1] = (B1 - B2 + wsd2 * D1) * du;
  } else {
    B1 = B1 + yL(k + 1, w * r[bmax - 1]) * Qkbd[bmax - 1] * gr.drdu(bmax - 1);
    B2 = B2 + (pow(r[bmax - 1], k - 1) / pow(r[bmax - 1], k + 2)) *
                Qkbd[bmax - 1] * gr.drdu(bmax - 1);
    v3[bmax - 1] = -2.0 * w * jL(k - 1, w * r[bmax - 1]) * B1 * du -
                   2.0 * ((2.0 * k + 1.0) / (w * w)) * B2 * du;
  }

  B3 = B3 + yL(k - 1, w * r[bmax - 1]) * Pkbd[bmax - 1] * gr.drdu(bmax - 1);

  v4[bmax - 1] = -2.0 * w * jL(k + 1, w * r[bmax - 1]) * B3 * du;

  for (auto i = bmax - 1; i >= 1; i--) {

    double ratio = r[i - 1] / r[i];

    if (cut_off) {
      B1 = B1 * powk(ratio) * (1.0 / ratio) +
           Qkbd[i - 1] * weights(i - 1) * gr.drduor(i - 1);
      B2 = B2 * powk(ratio) * ratio +
           Qkbd[i - 1] * weights(i - 1) * gr.drduor(i - 1);
      D1 = D1 * powk(ratio) * ratio +
           Qkbd[i - 1] * weights(i - 1) * r[i - 1] * gr.drdu(i - 1);
      v3[i - 1] = (B1 - B2 + wsd2 * D1) * du;
    } else {
      // j_{l-1}(wr1)y_{l+1}(wr2) term in integral that is integrated from r1 <= r2 <=infty
      B1 = B1 + yL(k + 1, w * r[i - 1]) * Qkbd[i - 1] * weights(i - 1) *
                  gr.drdu(i - 1);

      // r2^{l-1}/r1^{l+2} term in integral that is integrated from r1 <= r2 <=infty
      B2 = powk(ratio) * (1.0 / ratio) * B2 + Qkbd[i - 1] * weights(i - 1) *
                                                gr.drduor(i - 1) /
                                                (r[i - 1] * r[i - 1]);

      v3[i - 1] = wcubed * jL(k - 1, w * r[i - 1]) * B1 + (2.0 * k + 1.0) * B2;

      v3[i - 1] = -2.0 * odw2 * v3[i - 1] * du;
    }

    // j_{l+1}(wr1)y_{l-1}(wr2) term in integral that is integrated from r1 <= r2 <=infty
    B3 = B3 + yL(k - 1, w * r[i - 1]) * Pkbd[i - 1] * weights(i - 1) *
                gr.drdu(i - 1);
    v4[i - 1] = -2.0 * w * jL(k + 1, w * r[i - 1]) * B3 * du;
  }
}

//==============================================================================

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

//==============================================================================

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

void gk_ab_freqw(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 std::vector<double> &g0, std::vector<double> &ginf,
                 const std::size_t maxi, const double w) {

  auto fgfg = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  yk_ijk_gen_impl_freq(k, fgfg, Fa.grid(), g0, ginf, maxi, w);
}

//==============================================================================

void hk_ab_freqw(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 std::vector<double> &b0, std::vector<double> &binf,
                 std::size_t maxi, const double w) {

  auto fgfg = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  yk_ijk_gen_impl_freq(k, fgfg, Fa.grid(), b0, binf, maxi, w);
}

//==============================================================================

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

  // P^k_ij and Q^k_ij of arXiv:2602.17129:
  for (std::size_t i = 0; i < gr.size(); i++) {
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

} // namespace Coulomb
