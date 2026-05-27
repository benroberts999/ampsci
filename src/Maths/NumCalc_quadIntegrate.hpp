#pragma once
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_coeficients.hpp"
#include "qip/Vector.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

//! Numerical integration and differentiation routines.
namespace NumCalc {

//! Quadrature integration order; must be odd and in [1, 13]
constexpr std::size_t Nquad = 13;
static_assert(
  Nquad >= 1 && Nquad <= 13 && Nquad % 2 != 0,
  "\nFAIL10 in NumCalc: Nquad must be in [1,13], and must be odd\n");

// Instantiate quadrature coefficients for the chosen order:
constexpr QintCoefs<Nquad> quintcoef;
constexpr auto cq = quintcoef.cq;
constexpr auto dq_inv = quintcoef.dq_inv;

//==============================================================================
/*!
  @brief Quadrature integration of one or more vectors over a regular grid in t.
  @details
  Integrates the element-wise product `f1[i] * rest[i] * ...` from index `beg`
  to `end-1` (i.e., `end` is exclusive), multiplied by the step size `dt`.

  The grid must be uniform in the parametric variable t, but not necessarily
  in the physical variable x. For a non-uniform x grid, fold the Jacobian
  (dx/dt) into the integrand: pass `f(x) * (dx/dt)` rather than `f(x)` alone.

  End-point corrections (using `Nquad`-point quadrature weights) are applied
  at the start and end of the range. The additivity property is preserved:
  `int(a->b) + int(b->c) == int(a->c)`.

  No safety checks are performed. Requirements:
  - `end - beg > 2 * Nquad`
  - `end >= Nquad`
  - `beg + 2 * Nquad <= f1.size()`

  @param dt   Step size in the parametric variable t.
  @param beg  First grid index (inclusive).
  @param end  Last grid index (exclusive); if 0, defaults to `f1.size()`.
  @param f1   First vector to integrate (include Jacobian if grid is non-uniform in x).
  @param rest Additional vectors; all multiplied element-wise with f1.
  @returns    Integral of the element-wise product times dt.
*/
template <typename C, typename... Args>
inline double integrate(const double dt, std::size_t beg, std::size_t end,
                        const C &f1, const Args &...rest) {

  const auto max_grid = f1.size();
  if (end == 0)
    end = max_grid;
  const auto end_mid = std::min(max_grid - Nquad, end);
  const auto start_mid = std::max(Nquad, beg);

  double Rint_ends = qip::inner_product_sub(beg, Nquad, cq, f1, rest...);

  const double Rint_mid =
    qip::inner_product_sub(start_mid, end_mid, f1, rest...);

  for (auto i = end_mid; i < end; ++i) {
    Rint_ends += cq[end_mid + Nquad - i - 1] * qip::multiply_at(i, f1, rest...);
  }
  return (Rint_mid + dq_inv * Rint_ends) * dt;
}

//==============================================================================
/*!
  @brief Computes the derivative df/dr on a non-uniform grid.
  @details
  Given f sampled on a non-uniform grid r(t) with uniform step dt,
  computes df/dr using finite difference coefficients:

  `df/dr = (df/di) / (dt * dr/dt)`

  where i is the grid index. Coefficients from:
  http://en.wikipedia.org/wiki/Finite_difference_coefficient

  Uses a 7-point stencil in the interior; lower-order one-sided differences
  near the endpoints.

  For `order > 1`, applies the derivative recursively.

  @param f     Function values on the grid.
  @param drdt  Jacobian dr/dt at each grid point.
  @param dt    Uniform step size in the parametric variable t.
  @param order Derivative order (default 1); higher orders applied recursively.
  @returns     df/dr at each grid point.
*/
template <typename T>
inline std::vector<T> derivative(const std::vector<T> &f,
                                 const std::vector<T> &drdt, const T dt,
                                 const int order = 1) {
  // df/dr = df/dt * dt/dr = (df/dt) / (dr/dt)
  //       = (df/di) * (di/dt) / (dr/dt)
  //       = (df/di) / (dt * dr/dt)

  std::size_t num_points = f.size();
  std::vector<T> df(num_points);
  if (num_points < 4)
    return df;

  df[0] = (f[1] - f[0]) / (dt * drdt[0]);
  df[num_points - 1] =
    (f[num_points - 1] - f[num_points - 2]) / (dt * drdt[num_points - 1]);

  df[1] = (f[2] - f[0]) / (2 * dt * drdt[1]);
  df[num_points - 2] =
    (f[num_points - 1] - f[num_points - 3]) / (2 * dt * drdt[num_points - 2]);

  df[2] = (f[0] - 8 * f[1] + 8 * f[3] - f[4]) / (12 * dt * drdt[2]);
  df[num_points - 3] = (f[num_points - 5] - 8 * f[num_points - 4] +
                        8 * f[num_points - 2] - f[num_points - 1]) /
                       (12 * dt * drdt[num_points - 3]);

  df[3] = (-1 * f[0] + 9 * f[1] - 45 * f[2] + 45 * f[4] - 9 * f[5] + 1 * f[6]) /
          (60 * dt * drdt[3]);
  df[num_points - 4] =
    (-1 * f[num_points - 7] + 9 * f[num_points - 6] - 45 * f[num_points - 5] +
     45 * f[num_points - 3] - 9 * f[num_points - 2] + 1 * f[num_points - 1]) /
    (60 * dt * drdt[num_points - 4]);

  for (std::size_t i = 4; i < (num_points - 4); i++) {
    df[i] = ((1.0 / 8) * f[i - 4] - (4.0 / 3) * f[i - 3] + 7 * f[i - 2] -
             28 * f[i - 1] - (1.0 / 8) * f[i + 4] + (4.0 / 3) * f[i + 3] -
             7 * f[i + 2] + 28 * f[i + 1]) /
            (35 * dt * drdt[i]);
  }

  if (order > 1)
    df = derivative(df, drdt, dt, order - 1);

  return df;
}

//==============================================================================

//! Integration direction for additivePIntegral
enum Direction {
  //! Integrate from 0 outward to r
  zero_to_r,
  //! Integrate from infinity inward to r
  r_to_inf
};

/*!
  @brief Additively accumulates a partial integral into `answer`.
  @details
  Computes and adds to `answer`:

  `answer(r) += f(r) * Int[g(r') * h(r'), {r', 0, r}]`

  (or from infinity, depending on `direction`). Uses the trapezoidal rule.

  @note This is additive (`+=`). `answer` must already exist and be
  initialised (typically to zero).

  @param answer  Accumulation vector (modified in place).
  @param f       Prefactor at each grid point.
  @param g       First integrand factor.
  @param h       Second integrand factor.
  @param gr      Radial grid (provides dr/du and du).
  @param pinf    Upper index limit; 0 means use full grid.
*/
template <Direction direction, typename Real>
inline void
additivePIntegral(std::vector<Real> &answer, const std::vector<Real> &f,
                  const std::vector<Real> &g, const std::vector<Real> &h,
                  const Grid &gr, std::size_t pinf = 0) {
  const auto size = g.size();
  if (pinf == 0 || pinf >= size)
    pinf = size;
  const auto max = static_cast<int>(pinf - 1); // must be signed

  constexpr const bool forward = (direction == zero_to_r);
  constexpr const int inc = forward ? +1 : -1;
  const auto init = forward ? 0ul : std::size_t(max);
  const auto fin = forward ? std::size_t(max) : 0ul;

  Real x = int(init) < max ? 0.5 * g[init] * h[init] * gr.drdu()[init] : 0.0;
  answer[init] += f[init] * x * gr.du();
  for (auto i = int(init) + inc; i != int(fin) + inc; i += inc) {
    const auto im = std::size_t(i - inc);
    const auto i2 = std::size_t(i);
    x += 0.5 * (g[im] * h[im] * gr.drdu()[im] + g[i2] * h[i2] * gr.drdu()[i2]);
    answer[i2] += f[i2] * x * gr.du();
  }
  if (int(fin) < max)
    answer[fin] += 0.5 * f[fin] * g[fin] * h[fin] * gr.drdu()[fin] * gr.du();
}

/*!
  @brief Additively accumulates a partial integral into `answer` (two-function overload).
  @details
  As above, but integrand has only one factor g (no h):

  `answer(r) += f(r) * Int[g(r'), {r', 0, r}]`

  @note This is additive (`+=`). `answer` must already exist and be
  initialised (typically to zero).
*/
template <Direction direction, typename Real>
inline void additivePIntegral(std::vector<Real> &answer,
                              const std::vector<Real> &f,
                              const std::vector<Real> &g, const Grid &gr,
                              std::size_t pinf = 0) {
  const auto size = g.size();
  if (pinf == 0 || pinf >= size)
    pinf = size;
  const auto max = static_cast<int>(pinf - 1); // must be signed

  constexpr const bool forward = (direction == zero_to_r);
  constexpr const int inc = forward ? +1 : -1;
  const auto init = forward ? 0ul : std::size_t(max);
  const auto fin = forward ? std::size_t(max) : 0ul;

  Real x = int(init) < max ? 0.5 * g[init] * gr.drdu()[init] : 0.0;
  answer[init] += f[init] * x * gr.du();
  for (auto i = int(init) + inc; i != int(fin) + inc; i += inc) {
    const auto im = std::size_t(i - inc);
    const auto i2 = std::size_t(i);
    x += 0.5 * (g[im] * gr.drdu()[im] + g[i2] * gr.drdu()[i2]);
    answer[i2] += f[i2] * x * gr.du();
  }
  if (int(fin) < max)
    answer[fin] += 0.5 * f[fin] * g[fin] * gr.drdu()[fin] * gr.du();
}

//------------------------------------------------------------------------------
/*!
  @brief Returns the partial integral `f(r) * Int[g(r')*h(r'), {r', 0, r}]`.
  @details
  Allocates and returns the result vector; wrapper around additivePIntegral().
*/
template <Direction direction, typename Real>
std::vector<Real> partialIntegral(const std::vector<Real> &f,
                                  const std::vector<Real> &g,
                                  const std::vector<Real> &h, const Grid &gr,
                                  std::size_t pinf = 0) {
  std::vector<Real> answer(f.size(), 0.0);
  additivePIntegral<direction>(answer, f, g, h, gr, pinf);
  return answer;
}

//==============================================================================
//==============================================================================

//! Grid type for num_integrate
enum t_grid { linear, logarithmic };

//! Returns a function giving linearly-spaced x values: x(i) = a + i*dt
inline std::function<double(long unsigned)> linx(double a, double dt) {
  return [=](long unsigned i) { return a + double(i) * dt; };
}

//! Returns a function that always returns 1.0 (uniform Jacobian)
inline std::function<double(long unsigned)> one() {
  return [=](long unsigned) { return 1.0; };
}

//! Returns a function giving logarithmically-spaced x values: x(i) = a*exp(i*dt)
inline std::function<double(long unsigned)> logx(double a, double dt) {
  return [=](long unsigned i) { return a * std::exp(double(i) * dt); };
}

//==============================================================================
/*!
  @brief Numerically integrates a function f(x) over [a, b].
  @details
  Uses the same quadrature scheme as `integrate()`, on either a linear or
  logarithmic grid of `n_pts` points.

  Returns 0 if `a >= b` or `n_pts <= 1`.

  @param f      Function to integrate.
  @param a      Lower bound.
  @param b      Upper bound.
  @param n_pts  Number of quadrature points.
  @param type   Grid spacing: linear (default) or logarithmic.
  @returns      Approximate integral of f over [a, b].
*/
inline double num_integrate(const std::function<double(double)> &f, double a,
                            double b, long unsigned n_pts,
                            t_grid type = linear) {
  if (a >= b || n_pts <= 1)
    return 0.0;
  const auto dt = (type == linear) ? (b - a) / double(n_pts - 1) :
                                     std::log(b / a) / double(n_pts - 1);

  std::function<double(long unsigned)> x =
    (type == linear) ? linx(a, dt) : logx(a, dt);
  std::function<double(long unsigned)> dxdt =
    (type == linear) ? one() : logx(a, dt);

  double Rint_s = 0.0;
  for (long unsigned i = 0; i < Nquad; i++) {
    Rint_s += cq[i] * f(x(i)) * dxdt(i);
  }

  double Rint_m = 0.0;
  for (auto i = Nquad; i < n_pts - Nquad; i++) {
    Rint_m += f(x(i)) * dxdt(i);
  }

  double Rint_e = 0;
  for (auto i = n_pts - Nquad; i < n_pts; i++) {
    Rint_e += cq[n_pts - i - 1] * f(x(i)) * dxdt(i);
  }

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;
}

} // namespace NumCalc
