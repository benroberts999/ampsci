#pragma once
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_coeficients.hpp"
#include "qip/Vector.hpp"
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

//! Numerical integration and differentiation. Bit of a mess right now..
namespace NumCalc {

// Quadrature integration order [1,13], only odd
constexpr std::size_t Nquad = 13;
static_assert(
    Nquad >= 1 && Nquad <= 13 && Nquad % 2 != 0,
    "\nFAIL10 in NumCalc: Nquad must be in [1,13], and must be odd\n");

// instantiate coefs for correct Nquad order:
constexpr QintCoefs<Nquad> quintcoef;
constexpr auto cq = quintcoef.cq;
constexpr auto dq_inv = quintcoef.dq_inv;

//******************************************************************************
template <typename C, typename... Args>
inline double integrate(const double dt, std::size_t beg, std::size_t end,
                        const C &f1, const Args &... rest)
// // Note: includes no safety checks!
// // Integrates from (point) beg to end-1 (i.e., not including end)
// // Require:
// //   * (beg-end) > 2*Nquad
// //   * end - Nquad > Nquad
// //   * beg + 2*Nquad
// // NB: end-point corrections only applied if beg < Nquad, end > max - nquad
// // Not sure if this is best choice - but it ensures that:
// // int_{a->b} + int_{b->c} = int_{a->c}
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "4");

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

//******************************************************************************
template <typename T>
inline std::vector<T> derivative(const std::vector<T> &f,
                                 const std::vector<T> &drdt, const T dt,
                                 const int order = 1)
// df/dr = df/dt * dt/dr = (df/dt) / (dr/dt)
//       = (df/di) * (di/dt) / (dr/dt)
//       = (df/di) / (dt * dr/dt)
// coeficients from: http://en.wikipedia.org/wiki/Finite_difference_coefficient
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

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

//******************************************************************************
enum Direction { zero_to_r, r_to_inf };
template <Direction direction, typename Real>
inline void
additivePIntegral(std::vector<Real> &answer, const std::vector<Real> &f,
                  const std::vector<Real> &g, const std::vector<Real> &h,
                  const Grid &gr, std::size_t pinf = 0)
//  answer(r) += f(r) * Int[g(r')*h(r') , {r',0,r}]
// Note '+=' - this is additive!
// Note: vector answer must ALREADY exist, and be correctly initialised!
{
  const auto size = g.size();
  if (pinf == 0 || pinf >= size)
    pinf = size;
  const auto max = static_cast<int>(pinf - 1); // must be signed

  constexpr const bool forward = (direction == zero_to_r);
  constexpr const int inc = forward ? +1 : -1;
  const auto init = forward ? 0ul : std::size_t(max);
  const auto fin = forward ? std::size_t(max) : 0ul;

  // Quadrature weights; doesn't seem to help
  // auto Wquad = [&](std::size_t i) {
  //   return i < Nquad ? cq[i] * dq_inv
  //                    : i >= size - Nquad ? cq[size - i - 1] * dq_inv : 1.0;
  // };

  Real x = int(init) < max ? 0.5 * g[init] * h[init] * gr.drdu()[init] : 0.0;
  answer[init] += f[init] * x * gr.du();
  for (auto i = int(init) + inc; i != int(fin) + inc; i += inc) {
    const auto im = std::size_t(i - inc);
    const auto i2 = std::size_t(i);
    x += 0.5 * (g[im] * h[im] * gr.drdu()[im] /** Wquad(im)*/ +
                g[i2] * h[i2] * gr.drdu()[i2] /** Wquad(i2)*/);
    answer[i2] += f[i2] * x * gr.du();
  }
  if (int(fin) < max)
    answer[fin] += 0.5 * f[fin] * g[fin] * h[fin] * gr.drdu()[fin] *
                   gr.du(); // * Wquad(fin);
}
//------------------------
template <Direction direction, typename Real>
inline void additivePIntegral(std::vector<Real> &answer,
                              const std::vector<Real> &f,
                              const std::vector<Real> &g, const Grid &gr,
                              std::size_t pinf = 0)
//  answer(r) += f(r) * Int[g(r') , {r',0,r}]
// Note '+=' - this is additive!
// Note: vector answer must ALREADY exist, and be correctly initialised!
{
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
    answer[fin] +=
        0.5 * f[fin] * g[fin] * gr.drdu()[fin] * gr.du(); // * Wquad(fin);
}

//------------------------------------------------------------------------------
template <Direction direction, typename Real>
std::vector<Real> partialIntegral(const std::vector<Real> &f,
                                  const std::vector<Real> &g,
                                  const std::vector<Real> &h, const Grid &gr,
                                  std::size_t pinf = 0) {
  std::vector<Real> answer(f.size(), 0.0);
  additivePIntegral(answer, f, g, h, gr, pinf);
  return answer;
}

//******************************************************************************
//******************************************************************************

enum t_grid { linear, logarithmic };

inline std::function<double(long unsigned)> linx(double a, double dt) {
  return [=](long unsigned i) { return a + double(i) * dt; };
}

inline std::function<double(long unsigned)> one() {
  return [=](long unsigned) { return 1.0; };
}

inline std::function<double(long unsigned)> logx(double a, double dt) {
  return [=](long unsigned i) { return a * std::exp(double(i) * dt); };
}

//******************************************************************************
inline double num_integrate(const std::function<double(double)> &f, double a,
                            double b, long unsigned n_pts,
                            t_grid type = linear) {
  //

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
