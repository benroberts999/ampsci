#pragma once
#include <array>
#include <iostream>
#include <vector>

namespace NumCalc {

// Quadrature integration order [1,13], only odd
static const std::size_t Nquad = 3;
static_assert(
    Nquad >= 1 && Nquad <= 13 && Nquad % 2 != 0,
    "\nFAIL10 in NumCalc: Nquad must be in [1,13], and must be odd\n");

//******************************************************************************
template <typename T>
inline std::vector<T> derivative(const std::vector<T> &f,
                                 const std::vector<T> &drdt, const T dt,
                                 const int order = 1)
/*
df/dr = df/dt * dt/dr = (df/dt) / (dr/dt) = (df/di) * (di/dt) / (dr/dt) =
= (df/di)  / (dt * dr/dt)
coeficients from: http://en.wikipedia.org/wiki/Finite_difference_coefficient
*/
{

  std::size_t ngp = f.size();
  std::vector<T> df(ngp);

  df[0] = (f[1] - f[0]) / (dt * drdt[0]);
  df[ngp - 1] = (f[ngp - 1] - f[ngp - 2]) / (dt * drdt[ngp - 1]);

  df[1] = (f[2] - f[0]) / (2 * dt * drdt[1]);
  df[ngp - 2] = (f[ngp - 1] - f[ngp - 3]) / (2 * dt * drdt[ngp - 2]);

  df[2] = (f[0] - 8 * f[1] + 8 * f[3] - f[4]) / (12 * dt * drdt[2]);
  df[ngp - 3] = (f[ngp - 5] - 8 * f[ngp - 4] + 8 * f[ngp - 2] - f[ngp - 1]) /
                (12 * dt * drdt[ngp - 3]);

  df[3] = (-1 * f[0] + 9 * f[1] - 45 * f[2] + 45 * f[4] - 9 * f[5] + 1 * f[6]) /
          (60 * dt * drdt[3]);
  df[ngp - 4] = (-1 * f[ngp - 7] + 9 * f[ngp - 6] - 45 * f[ngp - 5] +
                 45 * f[ngp - 3] - 9 * f[ngp - 2] + 1 * f[ngp - 1]) /
                (60 * dt * drdt[ngp - 4]);

  for (std::size_t i = 4; i < (ngp - 4); i++) {
    df[i] = ((1. / 8) * f[i - 4] - (4. / 3) * f[i - 3] + 7 * f[i - 2] -
             28 * f[i - 1] - (1. / 8) * f[i + 4] + (4. / 3) * f[i + 3] -
             7 * f[i + 2] + 28 * f[i + 1]) /
            (35 * dt * drdt[i]);
  }

  if (order > 1)
    df = derivative(df, drdt, dt, order - 1);

  return df;
}

//******************************************************************************
// Define the coeficients for quadrature integration:
template <std::size_t N> struct QintCoefs {};

template <> struct QintCoefs<13> {
  static const std::size_t N = 13;
  const std::array<double, N> cq{
      {1382741929621, 9535909891802, -5605325192308, 28323664941310,
       -32865015189975, 53315213499588, -41078125154304, 39022895874876,
       -13155015007785, 12465244770050, 3283609164916, 5551687979302,
       5206230892907}};
  const double dq_inv = 1. / 5230697472000;
};
template <> struct QintCoefs<11> {
  static const std::size_t N = 11;
  const std::array<double, N> cq{
      {262747265, 1637546484, -454944189, 3373884696, -2145575886, 3897945600,
       -1065220914, 1942518504, 636547389, 1021256716, 952327935}};
  const double dq_inv = 1. / 958003200;
};
template <> struct QintCoefs<9> {
  static const std::size_t N = 9;
  const std::array<double, N> cq{{2082753, 11532470, 261166, 16263486, -1020160,
                                  12489922, 5095890, 7783754, 7200319}};
  const double dq_inv = 1. / 7257600;
};
template <> struct QintCoefs<7> {
  static const std::size_t N = 7;
  const std::array<double, N> cq{
      {36799, 176648, 54851, 177984, 89437, 130936, 119585}};
  const double dq_inv = 1. / 120960;
};
template <> struct QintCoefs<5> {
  static const std::size_t N = 5;
  const std::array<double, N> cq{{475, 1902, 1104, 1586, 1413}};
  const double dq_inv = 1. / 1440;
};
template <> struct QintCoefs<3> {
  static const std::size_t N = 3;
  const std::array<double, N> cq{{9, 28, 23}};
  const double dq_inv = 1. / 24;
};
template <> struct QintCoefs<1> {
  static const std::size_t N = 1;
  const std::array<double, N> cq{{1}};
  const double dq_inv = 1. / 2;
};

// instantiate coefs for correct Nquad order:
static const QintCoefs<Nquad> quintcoef;
static const auto &cq = quintcoef.cq;
static const auto dq_inv = quintcoef.dq_inv;
//******************************************************************************
template <typename C>
inline double integrate(const C &f1, const double dt = 1., std::size_t beg = 0,
                        std::size_t end = 0)
/*
Note: includes no safety checks!
Integrates from (point) beg to end-1 (i.e., not including end)
Require:
  * (beg-end) > 2*Nquad
  * end - Nquad > Nquad
  * beg + 2*Nquad
*/
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i];

  double Rint_e = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate1

//******************************************************************************
template <typename C>
inline double integrate(const C &f1, const C &f2, const double dt = 1.,
                        std::size_t beg = 0, std::size_t end = 0)
// Copy-paste from above
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i] * f2[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i] * f2[i];

  double Rint_e = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1] * f2[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate 2

//******************************************************************************
template <typename C>
inline double integrate(const C &f1, const C &f2, const C &f3,
                        const double dt = 1., std::size_t beg = 0,
                        std::size_t end = 0)
// Copy-paste from above
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i] * f2[beg + i] * f3[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i] * f2[i] * f3[i];

  double Rint_e = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1] * f2[end - i - 1] * f3[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate 3

//******************************************************************************
template <typename C>
inline double integrate(const C &f1, const C &f2, const C &f3, const C &f4,
                        const double dt = 1., std::size_t beg = 0,
                        std::size_t end = 0)
// Copy-paste from above
{

  if (end == 0)
    end = f1.size();

  // if (end - beg < 2 * Nquad)
  //   std::cerr << "\nFAIL 71 in INT: interval too small\n";

  double Rint_s = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_s += cq[i] * f1[beg + i] * f2[beg + i] * f3[beg + i] * f4[beg + i];

  double Rint_m = 0;
  for (auto i = beg + Nquad; i < end - Nquad; i++)
    Rint_m += f1[i] * f2[i] * f3[i] * f4[i];

  double Rint_e = 0;
  for (std::size_t i = 0; i < Nquad; i++)
    Rint_e += cq[i] * f1[end - i - 1] * f2[end - i - 1] * f3[end - i - 1] *
              f4[end - i - 1];

  return (Rint_m + dq_inv * (Rint_s + Rint_e)) * dt;

} // END integrate 4

} // namespace NumCalc
