#include "Coulomb/CoulombBreit.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

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
  // Test the vk_ab_freqw function
  // For now: just check 1.0e-3 matches 1.0e-9 (not a great test)
  // Test might not mean anything
  // Ideally: Real test?

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
    std::vector<std::vector<double>> V0_w(omegas.size());
    std::vector<std::vector<double>> Vinf_w(omegas.size());
    for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
      std::vector<double> v1, v2, v3, v4;
      Coulomb::vk_ab_freqw(k, Fa, Fb, gr, v1, v2, v3, v4, 0, omegas[iw]);
      V0_w[iw] = v1 + v2;
      Vinf_w[iw] = v3 + v4;
    }

    // Static limits to compare to
    const std::vector<double> V0_s =
      ((Fa.kappa() - Fb.kappa()) / k) * (g0_km1_s - g0_kp1_s) -
      (h0_km1_s - h0_kp1_s);
    const std::vector<double> Vinf_s =
      ((Fa.kappa() - Fb.kappa()) / (k + 1)) * (ginf_km1_s - ginf_kp1_s) +
      (hinf_km1_s - hinf_kp1_s);

    const int width = 13;
    std::cout << "k = " << k << " [vk_ab]\n";
    fmt::print("{:>5}", "r");
    for (const auto w : omegas)
      fmt::print("      w={:.0e}", w);
    fmt::print("\n");

    // V0 (0 -> r1) integral and comparison to static case
    for (const auto r_t : r_tgts) {
      const auto i = gr.getIndex(r_t);
      fmt::print("{:>5.0e}", gr.r(i));
      for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
        fmt::print("{:>{}.4e}", V0_w[iw][i], width);
      }
      fmt::print("{:>{}.4e}", V0_s[i], width);
      fmt::print(" {:+.1e}\n", V0_w[omegas.size() - 2][i] / V0_s[i] - 1.0);
    }

    // Vinf (r -> infty) integral and comparison to static case
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

    // const auto &v_ref = v_w.back();

    // // ensure the 1.0e-3 case matches the 1.0e-9 one
    // const double eps = 1.0e-6;
    // for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
    //   const auto i = gr.getIndex(r_tgts[ir]);
    //   REQUIRE(v_w[omegas.size() - 2][i] == Approx(v_ref[i]).epsilon(eps));
    // }
    // std::cout << "\n";
  }
}
