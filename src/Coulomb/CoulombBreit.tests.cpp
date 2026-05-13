#include "Coulomb/CoulombBreit.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/SixJTable.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

namespace UnitTest {
//==============================================================================
// This is a Naive (slow, but simple + correct) implementation of vk1_ab
// This forms the "baseline" to compare against
inline std::vector<double> vk1_naive(int k, const DiracSpinor &Fa,
                                     const DiracSpinor &Fb, const double &w);

// This is a Naive (slow, but simple + correct) implementation of vk2_ab
// This forms the "baseline" to compare against
inline std::vector<double> vk2_naive(int k, const DiracSpinor &Fa,
                                     const DiracSpinor &Fb, const double &w);

// This is a Naive (slow, but simple + correct) implementation of vk3_ab
// This forms the "baseline" to compare against
inline std::vector<double> vk3_naive(int k, const DiracSpinor &Fa,
                                     const DiracSpinor &Fb, const double &w);

// This is a Naive (slow, but simple + correct) implementation of vk4_ab
// This forms the "baseline" to compare against
inline std::vector<double> vk4_naive(int k, const DiracSpinor &Fa,
                                     const DiracSpinor &Fb, const double &w);

} // namespace UnitTest

//==============================================================================
//! Test that the gk and vk formulas match the naive double integration result
TEST_CASE("fBreit: gk, vk formulas", "[fBreit][Breit][unit][k4]") {

  const auto radial_grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 6;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
    }
    orbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
  }

  std::vector<double> v1, v2, v3, v4;

  // Compare vk formulas to naive formulas
  using namespace qip::overloads;
  for (const auto &Fa : orbs) {
    for (const auto &Fb : orbs) {
      if (Fb < Fa)
        continue;
      auto [k0, ki] = Angular::kminmax_Ck(Fa.kappa(), Fb.kappa());
      for (int k = k0; k <= ki; k += 2) {
        Coulomb::vk_ab_freqw(k, Fa, Fb, Fa.grid(), v1, v2, v3, v4,
                             Fa.grid().num_points(), 0.1);
        std::vector<double> v1_naive = UnitTest::vk1_naive(k, Fa, Fb, 0.1);
        std::vector<double> v2_naive = UnitTest::vk2_naive(k, Fa, Fb, 0.1);
        std::vector<double> v3_naive = UnitTest::vk3_naive(k, Fa, Fb, 0.1);
        std::vector<double> v4_naive = UnitTest::vk4_naive(k, Fa, Fb, 0.1);

        const auto delv1 = std::abs(qip::compare(v1, v1_naive).first);
        CHECK(delv1 < 1.0e-14);
        const auto delv2 = std::abs(qip::compare(v2, v2_naive).first);
        CHECK(delv2 < 1.0e-14);
        const auto delv3 = std::abs(qip::compare(v3, v3_naive).first);
        CHECK(delv3 < 1.0e-14);
        const auto delv4 = std::abs(qip::compare(v4, v4_naive).first);
        CHECK(delv4 < 1.0e-14);
      }
    }
  }
}

//==============================================================================
TEST_CASE("FreqwBreit: symmetry", "[Coulomb][unit][fBreit][Breit]") {
  // X_ab = fa*gb + ga*fb is symmetric  -> gk(a,b) = gk(b,a)
  // Y_ab = fa*gb - ga*fb is antisymm.  -> hk(a,b) = -hk(b,a), hk(a,a) = 0
  // P^k_ii = Q^k_ii = 0 (kappa_i-kappa_i=0, Y_ii=0) -> vk(a,a) = 0

  const auto grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const auto Fa = DiracSpinor::exactHlike(1, -1, grid, 1.0); // 1s
  const auto Fb = DiracSpinor::exactHlike(2, 1, grid, 1.0);  // 2p_{1/2}

  using namespace qip::overloads;
  std::vector<double> g0_ab, ginf_ab, g0_ba, ginf_ba;
  std::vector<double> h0_ab, hinf_ab, h0_ba, hinf_ba;
  std::vector<double> h0_aa, hinf_aa;

  for (const int k : {0, 1, 2, 3}) {
    for (const double w : {0.01, 0.1, 1.0}) {

      Coulomb::gk_ab_freqw(k, Fa, Fb, g0_ab, ginf_ab, 0, w);
      Coulomb::gk_ab_freqw(k, Fb, Fa, g0_ba, ginf_ba, 0, w);
      REQUIRE(std::abs(qip::compare(g0_ab, g0_ba).first) < 1.0e-14);
      REQUIRE(std::abs(qip::compare(ginf_ab, ginf_ba).first) < 1.0e-14);

      Coulomb::hk_ab_freqw(k, Fa, Fb, h0_ab, hinf_ab, 0, w);
      Coulomb::hk_ab_freqw(k, Fb, Fa, h0_ba, hinf_ba, 0, w);
      REQUIRE(std::abs(qip::compare(h0_ab, -1.0 * h0_ba).first) < 1.0e-14);
      REQUIRE(std::abs(qip::compare(hinf_ab, -1.0 * hinf_ba).first) < 1.0e-14);
    }

    // Y_aa = 0 exactly
    Coulomb::hk_ab_freqw(k, Fa, Fa, h0_aa, hinf_aa, 0, 0.1);
    for (const auto v : h0_aa)
      REQUIRE(v == 0.0);
    for (const auto v : hinf_aa)
      REQUIRE(v == 0.0);
  }

  // P^k_ii = Q^k_ii = 0, so all vk outputs zero for equal states
  std::vector<double> v1, v2, v3, v4;
  Coulomb::vk_ab_freqw(1, Fa, Fa, Fa.grid(), v1, v2, v3, v4, 0, 0.1);
  for (const auto v : v1)
    REQUIRE(std::abs(v) < 1.0e-15);
  for (const auto v : v2)
    REQUIRE(std::abs(v) < 1.0e-15);
  for (const auto v : v3)
    REQUIRE(std::abs(v) < 1.0e-15);
  for (const auto v : v4)
    REQUIRE(std::abs(v) < 1.0e-15);
}

//==============================================================================
TEST_CASE("FreqwBreit: w->0 limit", "[Coulomb][unit][fBreit][Breit][b7]") {
  // gk_ab_freqw(w) -> gk_ab (static) as w -> 0.
  // Kernel replacement: -w(2k+1) j_k(wr) y_k(wr') -> r^k/r'^{k+1} as w->0.
  // Deviation from static is O((wr)^2) (??)

  using namespace qip::overloads;

  const auto grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const auto Fa = DiracSpinor::exactHlike(1, -1, grid, 10.0); // 1s
  const auto Fb = DiracSpinor::exactHlike(2, 1, grid, 10.0);  // 2p_{1/2}
  const auto &gr = Fa.grid();

  const std::vector<double> r_tgts = {1.0e-6, 0.01, 0.1, 5.0, 50.0};
  const std::vector<double> omegas = {1, 0.05, 1.0e-6, 1.0e-9};

  //-------------------------------------------------------------------
  fmt::print("\nFreq. dep. Breit: freqw(r,w) vs static\n");
  using FreqwFn =
    void (*)(int, const DiracSpinor &, const DiracSpinor &,
             std::vector<double> &, std::vector<double> &, std::size_t, double);
  using StaticFn =
    void (*)(int, const DiracSpinor &, const DiracSpinor &,
             std::vector<double> &, std::vector<double> &, std::size_t);
  const std::pair<FreqwFn, StaticFn> fn_pairs[] = {
    {Coulomb::gk_ab_freqw, Coulomb::gk_ab},
    {Coulomb::hk_ab_freqw, Coulomb::bk_ab},
  };
  const char *fn_names[] = {"gk_ab", "hk_ab/bk_ab"};

  for (std::size_t ifn = 0; ifn < 2; ++ifn) {
    const auto [freqw_fn, static_fn] = fn_pairs[ifn];
    fmt::print("\n{}:\n", fn_names[ifn]);
    for (int k = 0; k <= 3; ++k) {

      // Static reference
      std::vector<double> g0_s, ginf_s;
      static_fn(k, Fa, Fb, g0_s, ginf_s, 0);

      // Pre-compute freqw results for all omegas
      std::vector<std::vector<double>> g0_w(omegas.size());
      std::vector<std::vector<double>> gI_w(omegas.size());
      for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
        std::vector<double> g0_w_temp, gI_w_temp;
        freqw_fn(k, Fa, Fb, g0_w_temp, gI_w_temp, 0, omegas[iw]);

        g0_w[iw] = g0_w_temp;
        gI_w[iw] = gI_w_temp;
      }

      const auto g_stat = g0_s + ginf_s;

      // Print table
      const int width = 13;

      // Table headers
      std::cout << std::endl << "k = " << k << " : " << fn_names[ifn] << "\n";
      fmt::print("{:>5}", "r");
      for (const auto w : omegas)
        fmt::print("      w={:.0e}", w);
      fmt::print("{:>{}}\n", "static", width);

      // Table values
      // 0 -> r1 gk/hk integral and comparison to static case
      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        fmt::print("{:>5.0e}", r_tgts.at(ir));
        const auto ir2 = Fa.grid().getIndex(r_tgts.at(ir));
        for (std::size_t iw = 0; iw < omegas.size(); ++iw) {

          fmt::print("{:>{}.4e}", g0_w[iw][ir2], width);
        }
        fmt::print("{:>{}.4e}", g0_s[ir2], width);
        fmt::print(" {:+.1e}\n",
                   g0_w[omegas.size() - 2][ir2] / g0_s[ir2] - 1.0);
      }

      // r1 -> infty gk/hk integral and comparison to static case

      std::cout << std::endl << "r1 -> infty" << std::endl;
      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        fmt::print("{:>5.0e}", r_tgts.at(ir));
        const auto ir2 = Fa.grid().getIndex(r_tgts.at(ir));
        for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
          fmt::print("{:>{}.4e}", gI_w[iw][ir2], width);
        }
        fmt::print("{:>{}.4e}", ginf_s[ir2], width);
        fmt::print(" {:+.1e}\n",
                   gI_w[omegas.size() - 2][ir2] / ginf_s[ir2] - 1.0);
      }

      const double eps1 = 1.0e-6;
      const double eps2 = 1.0e-9;

      // 1. Ensure the w~1.0e-3 case matches static case for both 0 and Inf integral
      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        const auto ir2 = Fa.grid().getIndex(r_tgts.at(ir));
        REQUIRE(g0_w[omegas.size() - 2][ir2] + gI_w[omegas.size() - 2][ir2] ==
                Approx(g_stat[ir2]).epsilon(eps1));
        REQUIRE(gI_w[omegas.size() - 2][ir2] ==
                Approx(ginf_s[ir2]).epsilon(eps1));
      }

      // 2. w=1e-9 (last) agrees with static to eps for both 0 and Inf integral
      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        const auto ir2 = Fa.grid().getIndex(r_tgts.at(ir));
        REQUIRE(g0_w[omegas.size() - 1][ir2] ==
                Approx(g0_s[ir2]).epsilon(eps2));
        REQUIRE(gI_w[omegas.size() - 1][ir2] ==
                Approx(ginf_s[ir2]).epsilon(eps2));
      }

      // // 3. Ensure sign and order-of-magnitude does not change
      // for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
      //   for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
      //     // Ensure sign is same
      //     REQUIRE(g_w[iw][ir] / g_stat[ir] > 0.0);

      //     // This FAILS sometimes:
      //     // Ensure All within factior of ~2 from static case
      //     // if (iw != 0)
      //     // Might not be wrong, but should be tested!
      //     // Also w=5 and w=10 breaks the trend? Again, maybe correct, but maybe not!
      //     CHECK(g_w[iw][ir] == Approx(g_stat[ir]).epsilon(0.1));
      //   }
      // }

    } // k loop
  } // fn_pairs loop

  //--------------------------------------
  // Test that

  for (int k = 1; k <= 3; ++k) {

    // Fill static g0, ginf, h0 and hinf
    std::vector<double> g0_km1_s, g0_kp1_s, ginf_km1_s, ginf_kp1_s, h0_km1_s,
      h0_kp1_s, hinf_km1_s, hinf_kp1_s;
    Coulomb::gk_ab_freqw(k - 1, Fa, Fb, g0_km1_s, ginf_km1_s);
    Coulomb::gk_ab_freqw(k + 1, Fa, Fb, g0_kp1_s, ginf_kp1_s);
    Coulomb::hk_ab_freqw(k - 1, Fa, Fb, h0_km1_s, hinf_km1_s);
    Coulomb::hk_ab_freqw(k + 1, Fa, Fb, h0_kp1_s, hinf_kp1_s);

    // Pre-compute total V0_w = v1 + v2 and Vinf_w = v3 + v4
    // Need to be careful about what these should match in small w limit:
    // V0_w should match ((kappa_a - kappa_b) / k) * (g0_{k-1} - g0_{k+1})  - (h0_{k-1} - h0_{k+1})
    // Vinf_w should match ((kappa_a - kappa_b) / (k + 1)) * (gI_{k-1} - gI_{k+1})  + (hI_{k-1} - hI_{k+1})
    // v2 and v4 should go to zero in small omega limit
    std::vector<std::vector<double>> V0_w(omegas.size());
    std::vector<std::vector<double>> Vinf_w(omegas.size());
    std::vector<std::vector<double>> V2_w(omegas.size());
    std::vector<std::vector<double>> V4_w(omegas.size());
    std::vector<std::vector<double>> V1_w(omegas.size());
    std::vector<std::vector<double>> V3_w(omegas.size());
    for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
      std::vector<double> v1, v2, v3, v4;
      Coulomb::vk_ab_freqw(k, Fa, Fb, gr, v1, v2, v3, v4, 0, omegas[iw]);
      V0_w[iw] = v1 + v2;
      Vinf_w[iw] = v3 + v4;
      V2_w[iw] = v2;
      V4_w[iw] = v4;
      V1_w[iw] = v1;
      V3_w[iw] = v3;
    }

    // Require that v2 and v4 go to zero in w -> 0 limit
    for (const auto r_t : r_tgts) {
      const auto i = gr.getIndex(r_t);
      // std::cout << V2_w[omegas.size() - 1][i] << std::endl;
      REQUIRE(abs(V2_w[omegas.size() - 1][i]) <=
              1.0e-12 * abs(V1_w[omegas.size() - 1][i]));
      REQUIRE(abs(V4_w[omegas.size() - 1][i]) <=
              1.0e-12 * abs(V3_w[omegas.size() - 1][i]));
    }

    // Static limits to compare to
    const std::vector<double> V0_s =
      (double(Fa.kappa() - Fb.kappa()) / k) * (g0_km1_s - g0_kp1_s) -
      (h0_km1_s - h0_kp1_s);
    const std::vector<double> Vinf_s =
      (double(Fa.kappa() - Fb.kappa()) / (k + 1)) * (ginf_km1_s - ginf_kp1_s) +
      (hinf_km1_s - hinf_kp1_s);

    const int width = 13;
    std::cout << std::endl << "k = " << k << " [vk_ab]\n";
    fmt::print("{:>5}", "r");
    for (const auto w : omegas)
      fmt::print("      w={:.0e}", w);
    fmt::print("\n");

    //!!! Should write definitive tests here. I am currently only printing tables

    // V0 (r2: 0 -> r1) integral and comparison to static case
    for (const auto r_t : r_tgts) {
      const auto i = gr.getIndex(r_t);
      fmt::print("{:>5.0e}", gr.r(i));
      for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
        fmt::print("{:>{}.4e}", V0_w[iw][i], width);
      }
      fmt::print("{:>{}.4e}", V0_s[i], width);
      fmt::print(" {:+.1e}\n", V0_w[omegas.size() - 2][i] / V0_s[i] - 1.0);
    }

    // Vinf (r2: r1 -> infty) integral and comparison to static case
    std::cout << std::endl << "r1 -> infty" << std::endl;
    for (const auto r_t : r_tgts) {
      const auto i = gr.getIndex(r_t);
      fmt::print("{:>5.0e}", gr.r(i));
      for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
        fmt::print("{:>{}.4e}", Vinf_w[iw][i], width);
      }
      fmt::print("{:>{}.4e}", Vinf_s[i], width);
      fmt::print(" {:+.1e}\n", Vinf_w[omegas.size() - 2][i] / Vinf_s[i] - 1.0);
    }
  }
}

//============================================================================
//============================================================================

//============================================================================
inline std::vector<double> UnitTest::vk1_naive(int k, const DiracSpinor &Fa,
                                               const DiracSpinor &Fb,
                                               const double &w) {

  const auto &gr = Fa.grid();
  const auto &r = gr.r();
  const auto &du = gr.du();
  const auto &num_points = gr.num_points();
  std::vector<double> v1(num_points);
  const auto odw2 = 1.0 / (w * w);

  // Quadrature integration weights:
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto phiL = SphericalBessel::PhiL;
  const auto psiL = SphericalBessel::PsiL;

  auto Xij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  auto Yij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  std::vector<double> Pkab(num_points);
  for (int i = 0; i < gr.num_points(); i++) {
    Pkab[i] = (double(Fa.kappa() - Fb.kappa()) / k) * Xij(i) - Yij(i);
  }

  // double integral
  for (auto n = 1ul; n < v1.size(); n++) {
    double A1 = 0.0; // keep track of running inner integral
    double A2 = 0.0;
    for (auto j = 0; j < n; j++) {
      double ratiokm1 = qip::pow(r[j] / r[n], k - 1);
      // first term
      A1 += ratiokm1 * phiL(k - 1, w * r[j], true) * Pkab[j] * weights(j) *
            gr.drdu(j);
      A2 += ratiokm1 * Pkab[j] * weights(j) * gr.drdu(j);
    } // integral over "r2" (j)
    double odr3 = 1.0 / (r[n] * r[n] * r[n]);
    v1[n] =
      2 * (2 * k + 1.0) * odw2 * odr3 *
      (A1 * psiL(k + 1, w * r[n], false) + A2 * psiL(k + 1, w * r[n], true)) *
      du;
  } // integral over "r1" (n)

  return v1;
}

//============================================================================
inline std::vector<double> UnitTest::vk2_naive(int k, const DiracSpinor &Fa,
                                               const DiracSpinor &Fb,
                                               const double &w) {

  const auto &gr = Fa.grid();
  const auto &r = gr.r();
  const auto &du = gr.du();
  const auto &num_points = gr.num_points();
  std::vector<double> v2(num_points);
  const auto w2dfact =
    2.0 * w * w / ((2 * k + 3.0) * (2 * k + 1.0) * (2 * k - 1.0));

  // Quadrature integration weights:
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto phiL = SphericalBessel::PhiL;
  const auto psiL = SphericalBessel::PsiL;

  auto Xij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  auto Yij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  std::vector<double> Qkab(num_points);
  for (int i = 0; i < gr.num_points(); i++) {
    Qkab[i] = (double(Fa.kappa() - Fb.kappa()) / (k + 1.0)) * Xij(i) + Yij(i);
  }

  // double integral
  for (auto n = 1ul; n < v2.size(); n++) {
    for (auto j = 0; j < n; j++) {
      double ratiok = qip::pow(r[j] / r[n], k);
      v2[n] += ratiok * r[j] * phiL(k + 1, w * r[j], false) * Qkab[j] *
               weights(j) * gr.drdu(j);
    } // integral over "r2" (j)
    v2[n] = w2dfact * v2[n] * psiL(k - 1, w * r[n], false) * du;
  } // integral over "r1" (n)

  return v2;
}

//============================================================================
inline std::vector<double> UnitTest::vk3_naive(int k, const DiracSpinor &Fa,
                                               const DiracSpinor &Fb,
                                               const double &w) {

  const auto &gr = Fa.grid();
  const auto &r = gr.r();
  const auto &du = gr.du();
  const auto &num_points = gr.num_points();
  std::vector<double> v3(num_points);
  const auto odw2 = 1.0 / (w * w);

  // Quadrature integration weights:
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto phiL = SphericalBessel::PhiL;
  const auto psiL = SphericalBessel::PsiL;

  auto Xij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  auto Yij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  std::vector<double> Qkab(num_points);
  for (int i = 0; i < gr.num_points(); i++) {
    Qkab[i] = (double(Fa.kappa() - Fb.kappa()) / (k + 1.0)) * Xij(i) + Yij(i);
  }

  // double integral
  for (auto n = 0ul; n < v3.size(); n++) {
    double A1 = 0.0;
    double A2 = 0.0;
    for (auto j = n; j < v3.size(); j++) {
      double ratiokm1 = qip::pow(r[n] / r[j], k - 1);
      double odr2 = 1.0 / (r[j] * r[j]);
      // first term
      A1 += ratiokm1 * odr2 * psiL(k + 1, w * r[j], true) * Qkab[j] *
            weights(j) * gr.drduor(j);
      A2 += ratiokm1 * odr2 * Qkab[j] * weights(j) * gr.drduor(j);
    } // integral over "r2" (j)
    v3[n] =
      2 * (2 * k + 1.0) * odw2 *
      (A1 * phiL(k - 1, w * r[n], false) + A2 * phiL(k - 1, w * r[n], true)) *
      du;
  } // integral over "r1" (n)

  return v3;
}

//============================================================================
inline std::vector<double> UnitTest::vk4_naive(int k, const DiracSpinor &Fa,
                                               const DiracSpinor &Fb,
                                               const double &w) {

  const auto &gr = Fa.grid();
  const auto &r = gr.r();
  const auto &du = gr.du();
  const auto &num_points = gr.num_points();
  std::vector<double> v4(num_points);
  const double w2dfact =
    2.0 * w * w / ((2 * k + 3.0) * (2 * k + 1.0) * (2 * k - 1.0));

  // Quadrature integration weights:
  const auto weights = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  const auto phiL = SphericalBessel::PhiL;
  const auto psiL = SphericalBessel::PsiL;

  auto Xij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) + Fa.g(i) * Fb.f(i));
  };

  auto Yij = [&](std::size_t i) {
    return (Fa.f(i) * Fb.g(i) - Fa.g(i) * Fb.f(i));
  };

  std::vector<double> Pkab(num_points);
  for (int i = 0; i < gr.num_points(); i++) {
    Pkab[i] = (double(Fa.kappa() - Fb.kappa()) / k) * Xij(i) - Yij(i);
  }

  // double integral
  for (auto n = 0ul; n < v4.size(); n++) {
    double A1 = 0.0;
    double A2 = 0.0;
    for (auto j = n; j < v4.size(); j++) {
      double ratiok = qip::pow(r[n] / r[j], k);
      // first term
      v4[n] += ratiok * psiL(k - 1, w * r[j], false) * Pkab[j] * weights(j) *
               gr.drdu(j);
    } // integral over "r2" (j)
    v4[n] = w2dfact * r[n] * v4[n] * phiL(k + 1, w * r[n], false) * du;
  } // integral over "r1" (n)

  return v4;
}
