#include "MixedStatesContinuum.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <cmath>
#include <memory>
#include <vector>

//==============================================================================
//! Anderson mixing coefficients: the constrained least-squares
//! solution for known residuals, the history-dropping on an ill-conditioned
//! Gram matrix, and exact convergence of the extrapolation for a small
//! linear fixed-point problem whose damped iteration diverges.
TEST_CASE("External Field: Anderson coefficients",
          "[ExternalField][MixedStatesCntm][TDHFcntm][unit]") {
  using ExternalField::anderson_coefficients;

  // Single entry: the plain step
  {
    const std::vector<std::vector<double>> gram{{2.0}};
    const auto c = anderson_coefficients(gram);
    REQUIRE(c.size() == 1);
    REQUIRE(c[0] == Approx(1.0));
  }

  // Orthogonal residuals: minimising |c1 r1 + c2 r2|^2 with c1 + c2 = 1
  // gives c_i proportional to 1/|r_i|^2
  {
    const std::vector<std::vector<double>> gram{{4.0, 0.0}, {0.0, 1.0}};
    const auto c = anderson_coefficients(gram);
    REQUIRE(c.size() == 2);
    REQUIRE(c[0] == Approx(0.2));
    REQUIRE(c[1] == Approx(0.8));
    REQUIRE(c[0] + c[1] == Approx(1.0));
  }

  // Ill-conditioned Gram matrix (one negligible residual, oldest): the
  // oldest entry is dropped and coefficients are returned for the rest
  {
    const std::vector<std::vector<double>> gram{{1.0e-20, 0.0, 0.0, 0.0},
                                                {0.0, 1.0, 0.0, 0.0},
                                                {0.0, 0.0, 1.0, 0.0},
                                                {0.0, 0.0, 0.0, 1.0}};
    const auto c = anderson_coefficients(gram);
    REQUIRE(c.size() == 3);
    for (const auto ci : c) {
      REQUIRE(ci == Approx(1.0 / 3.0));
    }
  }

  // Linear fixed point x = A x + b in two dimensions, with eigenvalues of A
  // outside the unit circle (damped iteration diverges for any damping):
  // the Anderson extrapolation reaches the exact solution once the residual
  // history spans the space (three iterates), as used in the mixed-state
  // solvers
  {
    const double A[2][2] = {{1.5, 0.3}, {0.2, -1.2}};
    const double b[2] = {1.0, 1.0};
    // Exact: (1 - A) x = b
    const double M[2][2] = {{1.0 - A[0][0], -A[0][1]},
                            {-A[1][0], 1.0 - A[1][1]}};
    const double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
    const double x_exact[2] = {(M[1][1] * b[0] - M[0][1] * b[1]) / det,
                               (M[0][0] * b[1] - M[1][0] * b[0]) / det};

    std::vector<double> x{0.0, 0.0};
    std::vector<std::vector<double>> g_hist;
    std::vector<std::vector<double>> r_hist;
    int its = 0;
    double residual = 1.0;
    for (; its < 10; ++its) {
      const std::vector<double> g{A[0][0] * x[0] + A[0][1] * x[1] + b[0],
                                  A[1][0] * x[0] + A[1][1] * x[1] + b[1]};
      const std::vector<double> r{g[0] - x[0], g[1] - x[1]};
      residual = std::sqrt(r[0] * r[0] + r[1] * r[1]);
      if (residual < 1.0e-10)
        break;
      g_hist.push_back(g);
      r_hist.push_back(r);
      const auto m = r_hist.size();
      std::vector<std::vector<double>> gram(m, std::vector<double>(m, 0.0));
      for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
          gram[i][j] =
            r_hist[i][0] * r_hist[j][0] + r_hist[i][1] * r_hist[j][1];
        }
      }
      const auto c = anderson_coefficients(gram);
      const auto n_drop = long(m - c.size());
      g_hist.erase(g_hist.begin(), g_hist.begin() + n_drop);
      r_hist.erase(r_hist.begin(), r_hist.begin() + n_drop);
      REQUIRE(!c.empty());
      x = {0.0, 0.0};
      for (std::size_t i = 0; i < c.size(); ++i) {
        x[0] += c[i] * g_hist[i][0];
        x[1] += c[i] * g_hist[i][1];
      }
    }
    fmt::print("\nAnderson mixing, 2D linear map with |eigenvalues| > 1: "
               "converged in {} iterations, residual {:.1e}\n",
               its, residual);
    REQUIRE(its <= 4);
    REQUIRE(x[0] == Approx(x_exact[0]).epsilon(1.0e-8));
    REQUIRE(x[1] == Approx(x_exact[1]).epsilon(1.0e-8));
  }
}

//==============================================================================
//! AndersonHistory: its incrementally kept Gram matrix, with the drops of the
//! oldest entries on overflow and on ill-conditioning, must reproduce the
//! extrapolation of the full recomputation (anderson_coefficients over the
//! residuals), including the carried amplitude. Also DiracSpinor::add_scaled
//! against the out-of-place form.
TEST_CASE("External Field: Anderson history",
          "[ExternalField][MixedStatesCntm][TDHFcntm][unit]") {
  using ExternalField::anderson_coefficients;
  const auto grid =
    std::make_shared<const Grid>(1.0e-5, 20.0, 300, GridType::loglinear, 3.0);

  // Smooth functions of r of varying shape, standing in for the iterates
  const auto make = [&grid](int seed, double scale) {
    DiracSpinor F{0, -1, grid};
    for (std::size_t i = 0; i < grid->num_points(); ++i) {
      const auto r = grid->r(i);
      F.f(i) = scale * r * std::exp(-r / (1.0 + 0.3 * seed)) *
               std::cos(0.5 * seed * r);
      F.g(i) = 0.1 * scale * r * std::exp(-r / (1.5 + 0.2 * seed)) *
               std::sin(0.3 * seed * r);
    }
    return F;
  };

  {
    auto a = make(1, 1.0);
    const auto b = make(2, 1.0);
    const auto sum = a + 0.3 * b;
    a.add_scaled(0.3, b);
    REQUIRE((a - sum).norm2() <= 1.0e-28 * sum.norm2());
  }

  // Depth 4, so that the history overflows. Entry 2 has a negligible
  // residual, which makes the Gram matrix ill-conditioned while it is kept
  // (as the newest, then a middle, then the oldest entry), so entries are
  // dropped for conditioning at three steps
  constexpr std::size_t depth = 4;
  ExternalField::AndersonHistory history(depth);
  std::vector<DiracSpinor> g_ref;
  std::vector<DiracSpinor> r_ref;
  std::vector<double> k_ref;
  int n_conditioning_drops = 0;
  for (int n = 0; n < 10; ++n) {
    const auto g = make(n + 1, 1.0);
    const auto r = make(n + 3, n == 2 ? 1.0e-9 : 0.1);
    const auto k = 0.7 * n - 1.0;

    // Reference: the Gram matrix rebuilt from all the kept residuals
    g_ref.push_back(g);
    r_ref.push_back(r);
    k_ref.push_back(k);
    if (g_ref.size() > depth) {
      g_ref.erase(g_ref.begin());
      r_ref.erase(r_ref.begin());
      k_ref.erase(k_ref.begin());
    }
    const auto c = anderson_coefficients(r_ref);
    REQUIRE(!c.empty());
    const auto n_drop = long(r_ref.size() - c.size());
    if (n_drop > 0) {
      ++n_conditioning_drops;
    }
    g_ref.erase(g_ref.begin(), g_ref.begin() + n_drop);
    r_ref.erase(r_ref.begin(), r_ref.begin() + n_drop);
    k_ref.erase(k_ref.begin(), k_ref.begin() + n_drop);
    auto x_ref = c[0] * g_ref[0];
    auto k_x_ref = c[0] * k_ref[0];
    for (std::size_t i = 1; i < c.size(); ++i) {
      x_ref += c[i] * g_ref[i];
      k_x_ref += c[i] * k_ref[i];
    }

    history.push(g, r, r * r, k);
    DiracSpinor x{0, -1, grid};
    double k_x = 0.0;
    history.mix(&x, &k_x);
    REQUIRE(history.size() == c.size());
    REQUIRE((x - x_ref).norm2() <= 1.0e-24 * x_ref.norm2());
    REQUIRE(k_x == Approx(k_x_ref).epsilon(1.0e-12).margin(1.0e-12));
  }
  REQUIRE(n_conditioning_drops == 3);
}
