#include "SphericalBessel.hpp"
#include "SphericalBessel.testdata.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <iostream>
#include <vector>

TEST_CASE("Maths::SphericalBessel", "[jL][unit]") {

  const auto r_list = qip::loglinear_range(1.0e-5, 10.0, 1.0, 100);

  for (int L = 0; L <= 15; ++L) {
    for (auto r : r_list) {
      REQUIRE(SphericalBessel::jL(L, r) ==
              Approx(SphericalBessel::exactGSL_JL(L, r)).margin(1.0e-15));

      REQUIRE(SphericalBessel::jL(L, r) ==
              Approx(SphericalBessel::exactGSL_JL(L, r)).epsilon(1.0e-3));
    }

    if (L < 4) {
      const auto jl = SphericalBessel::fillBesselVec(L, r_list);
      const auto jlkr = SphericalBessel::fillBesselVec_kr(L, 0.15, r_list);
      for (auto i = 0ul; i < r_list.size(); ++i) {
        const auto r = r_list.at(i);
        REQUIRE(jl.at(i) == Approx(SphericalBessel::jL(L, r)));
        REQUIRE(jlkr.at(i) == Approx(SphericalBessel::jL(L, 0.15 * r)));
      }
    }
  }
}

//==============================================================================
TEST_CASE("Maths::JL_table", "[JL_table][jL][unit]") {

  const std::vector<double> r = qip::logarithmic_range(1.0e-5, 100.0, 20);
  // dense q grid, for interp etc.
  const std::vector<double> q = qip::logarithmic_range(1.0e-2, 1000.0, 5000);

  const std::vector<std::size_t> iq_exact = {
    0ul,   1ul,   2ul,          15ul,         279ul,
    556ul, 999ul, q.size() - 3, q.size() - 2, q.size() - 1};

  const int max_L = 4;
  // cell_average = false: this test asserts pointwise equality of stored
  // values with j_L(q*r). The default cell-averaging behaviour intentionally
  // breaks pointwise equality on coarse cells (see the
  // "fillBesselVec: cell-averaged quadrature" test below).
  SphericalBessel::JL_table table(max_L, q, r, /*cell_average=*/false);

  SECTION("at() returns correct values") {
    for (int L = 0; L <= max_L; ++L) {
      for (std::size_t iq = 0; iq < q.size(); iq += 150) {
        const auto &v = table.at(std::size_t(L), iq);
        REQUIRE(v.size() == r.size());
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::jL(L, q[iq] * r[ir])).epsilon(1e-12));
        }
      }
    }
  }

  SECTION("jL() returns first q_i >= q") {
    for (int L = 0; L <= max_L; ++L) {
      for (double qq :
           {0.99e-2, 1.0e-2, 0.1, 25.0, 164.7, 999.0, 1000.0, 1001.0}) {
        const auto &v = table.jL(L, qq);

        for (std::size_t ir = 0; ir < r.size(); ir++) {
          // This one may be quite rough, particularly for small jL values
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).margin(1e-2));
        }
      }

      // test "exact" cases
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto &v = table.jL(L, qq);
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      // test "exact" cases, with small error
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto eps = 1.0e-5;
        const auto &v1 = table.jL(L, qq * (1.0 + eps));
        const auto &v2 = table.jL(L, qq * (1.0 - eps));
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          // less exact: went to far!
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).margin(1e-2));
          // very exact
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::jL(
                           L, q[std::min(iq + 1, q.size() - 1)] * r[ir]))
                    .epsilon(1e-12));
          // very exact
          REQUIRE(v2[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      //test edge cases:
      const auto &v0 = table.jL(L, 0.1 * q.front());
      const auto &v1 = table.jL(L, 2.0 * q.back());
      for (std::size_t ir = 0; ir < r.size(); ir++) {
        REQUIRE(
          v0[ir] ==
          Approx(SphericalBessel::jL(L, q.front() * r[ir])).epsilon(1e-12));
        REQUIRE(
          v1[ir] ==
          Approx(SphericalBessel::jL(L, q.back() * r[ir])).epsilon(1e-12));
      }
    }
  }

  SECTION("jL_nearest() returns for closest q") {
    for (int L = 0; L <= max_L; ++L) {
      for (double qq :
           {0.99e-2, 1.0e-2, 0.1, 25.0, 164.7, 999.0, 1000.0, 1001.0}) {
        const auto &v = table.jL_nearest(L, qq);

        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).margin(1e-2));
        }
      }

      // test "exact" cases
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto &v = table.jL_nearest(L, qq);
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      // test "exact" cases, with small error
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto eps = 1.0e-5;
        const auto &v1 = table.jL_nearest(L, qq * (1.0 + eps));
        const auto &v2 = table.jL_nearest(L, qq * (1.0 - eps));
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-12));
          REQUIRE(v2[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      //test edge cases:
      const auto &v0 = table.jL_nearest(L, 0.1 * q.front());
      const auto &v1 = table.jL_nearest(L, 2.0 * q.back());
      for (std::size_t ir = 0; ir < r.size(); ir++) {
        REQUIRE(
          v0[ir] ==
          Approx(SphericalBessel::jL(L, q.front() * r[ir])).epsilon(1e-12));
        REQUIRE(
          v1[ir] ==
          Approx(SphericalBessel::jL(L, q.back() * r[ir])).epsilon(1e-12));
      }
    }
  }

  SECTION("jL_interp() interpolates correctly") {
    for (int L = 0; L <= max_L; ++L) {
      for (double qq : {1.0e-2, 0.1, 25.0, 164.7, 999.0, 1000.0}) {
        const auto &v = table.jL_interp(L, qq);
        REQUIRE(v.size() == r.size());
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          if (qq * r[ir] < 1.0) {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-4));
          } else if (qq * r[ir] < 10.0) {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-3));
          } else if (qq * r[ir] < 50.0) {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::jL(L, qq * r[ir])).margin(1e-3));
          } else {
            REQUIRE(v[ir] ==
                    Approx(SphericalBessel::jL(L, qq * r[ir])).margin(1e-2));
          }
        }
      }

      // test "exact" cases
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto &v = table.jL_interp(L, qq);
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          REQUIRE(v[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-12));
        }
      }

      // test "exact" cases, with small error
      for (std::size_t iq : iq_exact) {
        double qq = q[iq];
        const auto eps = 1.0e-8;
        const auto &v1 = table.jL_interp(L, qq * (1.0 + eps));
        const auto &v2 = table.jL_interp(L, qq * (1.0 - eps));
        for (std::size_t ir = 0; ir < r.size(); ir++) {
          // interp still not expected to be exact
          REQUIRE(v1[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-3));
          REQUIRE(v2[ir] ==
                  Approx(SphericalBessel::jL(L, qq * r[ir])).epsilon(1e-3));
        }
      }

      //test edge cases:
      const auto &v0 = table.jL_interp(L, 0.1 * q.front());
      const auto &v1 = table.jL_interp(L, 2.0 * q.back());
      for (std::size_t ir = 0; ir < r.size(); ir++) {
        REQUIRE(
          v0[ir] ==
          Approx(SphericalBessel::jL(L, q.front() * r[ir])).epsilon(1e-12));
        REQUIRE(
          v1[ir] ==
          Approx(SphericalBessel::jL(L, q.back() * r[ir])).epsilon(1e-12));
      }
    }
  }
}

//==============================================================================
TEST_CASE("Maths::SphericalBessel cell-averaged quadrature",
          "[jL][Bessel][unit][xyz]") {
  // Verifies the cell-averaging mode of fillBesselVec_kr:
  // integrating j_L(q*r) * psi(r) on a coarse log grid with cell-averaged
  // j_L values should match a high-resolution pointwise reference, while
  // the same coarse grid with pure point sampling should be much worse
  // once q is large enough for the cells to not resolve j_L's oscillation.
  // psi(r) is a smooth decaying envelope, which is the typical situation
  // (smooth wavefunction times oscillatory Bessel factor).

  auto psi = [](double r) { return r * std::exp(-0.5 * r); };

  // Dense linear reference grid: resolves j_L for all q values tested.
  const std::size_t N_fine = 5000;
  const double r_max = 25.0;
  std::vector<double> r_fine(N_fine);
  for (std::size_t i = 0; i < N_fine; ++i)
    r_fine[i] = double(i + 1) * r_max / double(N_fine);

  // Coarse log grid: cells widen with r and stop resolving j_L for q >> 1.
  const auto r_coarse = qip::logarithmic_range(1.0e-3, r_max, 60);

  auto integrate_trap = [](const std::vector<double> &r,
                           const std::vector<double> &fr) {
    double s = 0.0;
    for (std::size_t i = 1; i < r.size(); ++i)
      s += 0.5 * (fr[i - 1] + fr[i]) * (r[i] - r[i - 1]);
    return s;
  };

  auto integrate_jl_psi = [&](int L, double q, const std::vector<double> &r,
                              bool cell_average) {
    const auto jl = SphericalBessel::fillBesselVec_kr(L, q, r, cell_average);
    std::vector<double> integrand(r.size());
    for (std::size_t i = 0; i < r.size(); ++i)
      integrand[i] = jl[i] * psi(r[i]);
    return integrate_trap(r, integrand);
  };

  for (int L : {0, 1, 2, 3}) {
    for (double q : {1.0, 5.0, 20.0, 100.0}) {
      const double I_ref = integrate_jl_psi(L, q, r_fine, false);
      const double I_avg = integrate_jl_psi(L, q, r_coarse, true);
      const double I_pt = integrate_jl_psi(L, q, r_coarse, false);

      const double err_avg = std::abs(I_avg - I_ref);
      const double err_pt = std::abs(I_pt - I_ref);

      // Averaged coarse integral matches the dense reference reasonably well.
      REQUIRE(err_avg < 0.05 * std::abs(I_ref) + 1.0e-4);

      // For q large enough that cells do not resolve j_L, averaging is a
      // strict improvement over point sampling.
      if (q >= 20.0) {
        REQUIRE(err_avg < err_pt);
      }
    }
  }
}

//==============================================================================
TEST_CASE("Maths::SphericalBessel jL testdata", "[jL][Bessel][unit]") {
  for (const auto &data : UnitTest::jk_DATA) {
    const auto value = SphericalBessel::jL(data.k, data.x);
    // fmt::print("{:2} {:.1e} {:.6e} {:.6e} {:.1e}\n", data.k, data.x, value,
    //            data.value, (value - data.value) / data.value);

    if (std::abs(data.value) > 1.0e-10) {
      REQUIRE(value == Approx(data.value).epsilon(1e-10));
    } else {
      REQUIRE(value == Approx(data.value).epsilon(1e-5));
      REQUIRE(value == Approx(data.value).margin(1e-12));
    }
  }
}

// //==============================================================================
TEST_CASE("Maths::SphericalBessel yL testdata", "[yL][Bessel][unit]") {
  for (const auto &data : UnitTest::yk_DATA) {
    const auto value = SphericalBessel::yL(data.k, data.x);
    // fmt::print("{:2} {:.1e} {:.6e} {:.6e} {:.1e}\n", data.k, data.x, value,
    //            data.value, (value - data.value) / data.value);

    if (std::abs(data.value) > 1.0e-10) {
      REQUIRE(value == Approx(data.value).epsilon(1e-10));
    } else {
      REQUIRE(value == Approx(data.value).epsilon(1e-5));
      REQUIRE(value == Approx(data.value).margin(1e-12));
    }
  }
}

// //==============================================================================
TEST_CASE("Maths::SphericalBessel PhiL testdata", "[jL][Bessel][unit]") {
  // phi(x) = [(2k+1)!! / x^k] * j_k(x); tilde data stores 1 - phi
  for (const auto &d : UnitTest::PhikTilde_DATA) {
    // non-tilde: phi = 1 - (tilde value)

    // const auto value = SphericalBessel::PhiL(d.k, d.x);
    // const auto valuetilde = SphericalBessel::PhiL(d.k, d.x, true);
    // fmt::print("{:2} {:.1e} {:.6e} {:.6e} {:.1e}    {:.6e} {:.6e} {:.1e}\n",
    //            d.k, d.x, value, 1.0 - d.value,
    //            (value - (1.0 - d.value)) / (1.0 - d.value), valuetilde,
    //            -d.value, (valuetilde - (-d.value)) / -d.value);

    REQUIRE(SphericalBessel::PhiL(d.k, d.x) ==
            Approx(1.0 - d.value).epsilon(1e-10));

    // tilde: 1 - phi
    if (std::abs(d.value) > 1.0e-9) {
      REQUIRE(SphericalBessel::PhiL(d.k, d.x, true) ==
              Approx(-d.value).epsilon(1e-6));
      REQUIRE(SphericalBessel::PhiL(d.k, d.x, true) ==
              Approx(-d.value).margin(1e-12));
    } else if (d.k < 10) {
      REQUIRE(SphericalBessel::PhiL(d.k, d.x, true) ==
              Approx(-d.value).epsilon(1e-2));
      REQUIRE(SphericalBessel::PhiL(d.k, d.x, true) ==
              Approx(-d.value).margin(1e-12));
    } else {
      REQUIRE(SphericalBessel::PhiL(d.k, d.x, true) ==
              Approx(-d.value).epsilon(1e-1));
      REQUIRE(SphericalBessel::PhiL(d.k, d.x, true) ==
              Approx(-d.value).margin(1e-12));
    }
  }
}

//==============================================================================
TEST_CASE("Maths::SphericalBessel PsiL testdata", "[yL][Bessel][unit]") {
  // psi(x) = [x^{k+1} / (2k-1)!!] * y_k(x); tilde data stores 1 - psi
  for (const auto &d : UnitTest::PsikTilde_DATA) {

    // const auto value = SphericalBessel::PsiL(d.k, d.x);
    // const auto valuetilde = SphericalBessel::PsiL(d.k, d.x, true);
    // fmt::print("{:2} {:.1e} {:.6e} {:.6e} {:.1e}    {:.6e} {:.6e} {:.1e}\n",
    //            d.k, d.x, value, 1.0 - d.value,
    //            (value - (1.0 - d.value)) / (1.0 - d.value), valuetilde,
    //            -d.value, (valuetilde - (-d.value)) / -d.value);

    // non-tilde: psi = 1 - (tilde value)
    const auto value = SphericalBessel::PsiL(d.k, d.x);
    const auto valuetilde = SphericalBessel::PsiL(d.k, d.x, true);

    fmt::print("{:2} {:.1e} {:.6e} {:.6e} {:.1e}    {:.6e} {:.6e} {:.1e}\n",
               d.k, d.x, value, 1.0 - d.value,
               (value - (1.0 - d.value)) / (1.0 - d.value), valuetilde,
               -d.value, (valuetilde - (-d.value)) / -d.value);

    REQUIRE(SphericalBessel::PsiL(d.k, d.x) ==
            Approx(1.0 - d.value).epsilon(1e-10));

    // tilde: 1 - psi
    if (std::abs(d.value) > 1.0e-9) {
      REQUIRE(SphericalBessel::PsiL(d.k, d.x, true) ==
              Approx(-d.value).epsilon(1e-6));
      REQUIRE(SphericalBessel::PsiL(d.k, d.x, true) ==
              Approx(-d.value).margin(1e-12));
    } else if (d.k < 9) {
      REQUIRE(SphericalBessel::PsiL(d.k, d.x, true) ==
              Approx(-d.value).epsilon(1e-2));
      REQUIRE(SphericalBessel::PsiL(d.k, d.x, true) ==
              Approx(-d.value).margin(1e-12));
    } else {
      REQUIRE(SphericalBessel::PsiL(d.k, d.x, true) ==
              Approx(-d.value).epsilon(1e-1));
      REQUIRE(SphericalBessel::PsiL(d.k, d.x, true) ==
              Approx(-d.value).margin(1e-12));
    }
  }
}