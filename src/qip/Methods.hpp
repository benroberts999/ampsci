#pragma once
#include <cmath>
#include <functional>
#include <utility>

namespace qip {

//! Slow, but accurate, method of finding derivative of function (y) at a point
//! (x). Returns derivative + error estimate
/*! @details
 delta_target is target for |dy/dx_n - dy/dx_{n+1}| < delta_target (1.0e-6);
 dx is initial step-size used to find derivative of f (0.01);
 it_limit is maximum number of iterations.
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

//! Solve f(x) = 0 for x using Newtons method. Returns root + error estimate/
/*! @details
 x is be initial guess for the root;
 delta_target is target for |x_n - x_{n+1}| < delta_target (1.0e-6);
 dx is initial step-size used to find derivative of f (0.01);
 it_limit is maximum number of iterations.
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
//! Solve f(x) = 0 for x using Newtons method. Returns root + error estimate.
//! Enforced to be between bounds.
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
