#pragma once
#include <cmath>
#include <functional>
#include <utility>

namespace qip {

/*!
  @brief Numerical derivative of y(x) at point x; returns {dy/dx, error}.

  @details
  Iteratively halves the step size until
  |dy/dx_n - dy/dx_{n-1}| < delta_target.

  @tparam Function Callable as `y(x)`, returning a value convertible to Real.
  @tparam Real     Floating-point type.

  @param y            Function to differentiate.
  @param x            Point at which to evaluate the derivative.
  @param delta_target Convergence target (default 1e-6).
  @param dx           Initial step size (default 0.01).
  @param it_limit     Maximum number of iterations (default 250).
*/
template <typename Function, typename Real>
std::pair<Real, Real>
derivative(Function y, Real x, Real delta_target = Real{1.0e-6},
           Real dx = Real{0.01}, unsigned it_limit = 250) {
  Real dydx{0.0};
  Real delta{1.0};
  for (auto i = 0u;; ++i) {
    const auto dydx0 = dydx;
    const auto dy = 0.5 * (y(x + dx) - y(x - dx));
    dydx = dy / dx;
    delta = std::abs(dydx - dydx0);
    if (1.1 * delta < delta_target || i > it_limit)
      break;
    dx *= 0.5;
  }
  return {dydx, 1.1 * delta};
}

/*!
  @brief Solves f(x) = 0 using Newton's method; returns {root, error}.

  @tparam Function Callable as `f(x)`, returning a value convertible to Real.
  @tparam Real     Floating-point type.

  @param f            Function to solve.
  @param x            Initial guess.
  @param delta_target Convergence target for |x_n - x_{n+1}| (default 1e-6).
  @param dx           Initial step size for derivative (default 0.01).
  @param it_limit     Maximum number of iterations (default 250).
*/
template <typename Function, typename Real>
std::pair<Real, Real> Newtons(Function f, Real x,
                              Real delta_target = Real{1.0e-6},
                              Real dx = Real{0.01}, unsigned it_limit = 250) {
  Real delta{1.0};
  for (auto i = 0u;; ++i) {
    // const auto df = (f(x + dx) - f(x)) / dx;
    const auto df = derivative(f, x, delta_target, dx, 4).first;
    delta = f(x) / df;
    x = x - delta;
    if (std::abs(delta) < delta_target || i > it_limit)
      break;
  }
  return {x, 1.1 * delta};
}

/*!
  @brief Solves f(x) = 0 using Newton's method with bounds; returns {root, error}.

  @details Solution is clamped to [bounds.first, bounds.second].

  @param f            Function to solve.
  @param x            Initial guess.
  @param bounds       Allowed range [lower, upper].
  @param delta_target Convergence target for |x_n - x_{n+1}| (default 1e-6).
  @param dx           Initial step size for derivative (default 0.01).
  @param it_limit     Maximum number of iterations (default 250).
*/
template <typename Function, typename Real>
std::pair<Real, Real> Newtons(Function f, Real x, std::pair<Real, Real> bounds,
                              Real delta_target = Real{1.0e-6},
                              Real dx = Real{0.01}, unsigned it_limit = 250) {
  Real delta{1.0};
  for (auto i = 0u;; ++i) {
    // const auto df = (f(x + dx) - f(x)) / dx;
    const auto df = derivative(f, x, delta_target, dx, 4).first;
    delta = f(x) / df;
    x = x - delta;
    if (x < bounds.first) {
      x = bounds.first;
      break;
    }
    if (x > bounds.second) {
      x = bounds.second;
      break;
    }
    if (std::abs(delta) < delta_target || i > it_limit)
      break;
  }
  return {x, 1.1 * delta};
}

} // namespace qip
