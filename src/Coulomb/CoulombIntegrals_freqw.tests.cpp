#include "Coulomb/CoulombIntegrals.hpp"
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
TEST_CASE("FreqwBreit: w->0 limit", "[Coulomb][unit][fBreit][Breit]") {
  // gk_ab_freqw(w) -> gk_ab (static) as w -> 0.
  // Kernel replacement: -w(2k+1) j_k(wr) y_k(wr') -> r^k/r'^{k+1} as w->0.
  // Deviation from static is O((wr)^2).

  const auto grid = std::make_shared<const Grid>(
    GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const auto Fa = DiracSpinor::exactHlike(1, -1, grid, 10.0);
  const auto Fb = DiracSpinor::exactHlike(2, 1, grid, 10.0);
  const auto &gr = Fa.grid();

  fmt::print("\nFreq. dep. Breit: freqw(r,w) vs static\n");
  const std::vector<double> r_tgts = {1.0e-6, 0.01, 0.1, 5.0, 100.0};
  const std::vector<double> omegas = {10.0, 0.5, 0.05, 0.001, 1.0e-9};

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

  const int k = 3;
  for (std::size_t ifn = 0; ifn < 2; ++ifn) {
    const auto [freqw_fn, static_fn] = fn_pairs[ifn];
    fmt::print("\n{}:\n", fn_names[ifn]);
    for (int k = 0; k <= 3; ++k) {

      // Static reference
      std::vector<double> g0_s, ginf_s;
      static_fn(k, Fa, Fb, g0_s, ginf_s, 0);

      using namespace qip::overloads;

      // Pre-compute freqw results for all omegas
      std::vector<std::vector<double>> g_w(omegas.size());
      for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
        std::vector<double> g0_w, ginf_w;
        freqw_fn(k, Fa, Fb, g0_w, ginf_w, 0, omegas[iw]);

        g_w[iw] = g0_w + ginf_w;
      }

      const auto g_stat = g0_s + ginf_s;

      // Print table
      const int width = 13;

      // Table headers
      std::cout << "k = " << k << " : " << fn_names[ifn] << "\n";
      fmt::print("{:>5}", "r");
      for (const auto w : omegas)
        fmt::print("      w={:.0e}", w);
      fmt::print("{:>{}}\n", "static", width);

      // Table values
      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        fmt::print("{:>5.0e}", r_tgts.at(ir));
        for (std::size_t iw = 0; iw < omegas.size(); ++iw) {
          fmt::print("{:>{}.4e}", g_w[iw][ir], width);
        }
        fmt::print("{:>{}.4e}", g_stat[ir], width);
        fmt::print(" {:+.1e}\n", g_w[omegas.size() - 2][ir] / g_stat[ir] - 1.0);
      }

      const double eps = 1.0e-6;

      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        REQUIRE(g_w[omegas.size() - 2][ir] == Approx(g_stat[ir]).epsilon(eps));
      }

      // 2. w=1e-9 (last) agrees with static to eps
      for (std::size_t ir = 0; ir < r_tgts.size(); ++ir) {
        CHECK(g_w[omegas.size() - 1][ir] == Approx(g_stat[ir]).epsilon(eps));
      }

    } // k loop
  }   // fn_pairs loop
  std::cout << "\n";
}
