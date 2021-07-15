#pragma once
#include <cmath>
#include <functional>
#include <utility>

namespace qip {

//! Solve f(x) = 0 for x using Newtons method.
/*! @details
 dx is used to find derivative of f;
 x is be initial guess for the root;
 delta_target is target for |x_n - x_{n+1}| < delta_target;
 it_limit is maximum number of iterations.
*/
template <typename Function, typename Real>
std::pair<Real, Real> Newtons(Function f, Real x, Real dx = Real{0.01},
                              Real delta_target = Real{1.0e-6},
                              unsigned it_limit = 50) {
  Real delta{1.0};
  for (auto i = 0u;; ++i) {
    const auto df = (f(x + dx) - f(x)) / dx;
    delta = f(x) / df;
    x = x - delta;
    if (std::abs(delta) < delta_target || i > it_limit)
      break;
  }
  return {x, delta};
}

} // namespace qip
