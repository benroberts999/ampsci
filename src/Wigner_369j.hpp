#pragma once
#include <cmath>
#include <gsl/gsl_sf_coupling.h>

namespace Wigner {

/*
Wrapper functions to calculate wigner 3,6,9-J symbols.
Uses GSL:
https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=3j#coupling-coefficients
NOTE:
Since j always integer or half-integer, GSL inputs are always integer.
Three versions of each symbol:
 * 'regular', takes in double. Converts to integer safely. Slower (marginally),
    but easier
 * '_1' version takes integers as-is. Only work for
integer angular momentum (l).
 * and '_2' version, takes in 2*j (as an integer). Works for l and j
*/

//******************************************************************************
inline int l_k(int ka) { return (ka > 0) ? ka : -ka - 1; }
inline int twoj_k(int ka) { return (ka > 0) ? 2 * ka - 1 : -2 * ka - 1; }
inline double j_k(int ka) {
  return (ka > 0) ? double(ka) - 0.5 : double(-ka) - 0.5;
}

//******************************************************************************
inline int parity(int la, int lb, int k)
// Parity rule. Returns 1 only if la+lb+k is even
{
  if ((la + lb + k) % 2 == 0)
    return 1;
  return 0;
}

//******************************************************************************
inline int triangle(double j1, double j2, double J)
// Triangle rule.
{
  if ((j1 + j2 < J) || (fabs(j1 - j2) > J))
    return 0;
  return 1;
}

inline int triangle(int j1, int j2, int J) {
  // nb: can be called with wither j or twoj!
  if ((j1 + j2 < J) || (abs(j1 - j2) > J))
    return 0;
  return 1;
}

inline int sumsToZero(int m1, int m2, int m3) {
  if (m1 + m2 + m3 != 0)
    return 0;
  return 1;
}

inline int sumsToZero(double m1, double m2, double m3) {
  if (fabs(m1 + m2 + m3) > 0.0001)
    return 0;
  return 1;
}

//******************************************************************************
inline double threej(double j1, double j2, double j3, double m1, double m2,
                     double m3)
// Calculates wigner 3j symbol:
//   (j1 j2 j3)
//   (m1 m2 m3)
// Note: this function takes DOUBLE values.
// Works for l and j (integer and half-integer)
{
  if (triangle(j1, j2, j3) * sumsToZero(m1, m2, m3) == 0)
    return 0;
  int two_j1 = (int)round(2 * j1);
  int two_j2 = (int)round(2 * j2);
  int two_j3 = (int)round(2 * j3);
  int two_m1 = (int)round(2 * m1);
  int two_m2 = (int)round(2 * m2);
  int two_m3 = (int)round(2 * m3);
  return gsl_sf_coupling_3j(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

//------------------------------------------------------------------------------
inline double threej_1(int j1, int j2, int j3, int m1, int m2, int m3)
// Calculates wigner 3j symbol:
//   (j1 j2 j3)
//   (m1 m2 m3)
// Note: this function takes INTEGER values, only works for l (not half-integer
// j)!
{
  if (triangle(j1, j2, j3) * sumsToZero(m1, m2, m3) == 0)
    return 0;
  return gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3);
}

//------------------------------------------------------------------------------
inline double threej_2(int two_j1, int two_j2, int two_j3, int two_m1,
                       int two_m2, int two_m3)
// Calculates wigner 3j symbol:
//   (j1 j2 j3)
//   (m1 m2 m3)
// Note: this function takes INTEGER values, that have already multiplied by 2!
// Works for l and j (integer and half-integer)
{
  if (triangle(two_j1, two_j2, two_j3) * sumsToZero(two_m1, two_m2, two_m3) ==
      0)
    return 0;
  return gsl_sf_coupling_3j(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

//******************************************************************************
inline double cg(double j1, double m1, double j2, double m2, double J, double M)
// <j1 m1, j2 m2 | J M> = (-1)^(j1-j2+M) * Sqrt(2J+1) * (j1 j2  J)
// .                                                    (m1 m2 -M)
// (Last term is 3j symbol)
// Note: this function takes DOUBLE values.
// Works for l and j (integer and half-integer)
{
  if (triangle(j1, j2, J) * sumsToZero(m1, m2, -M) == 0)
    return 0;
  int two_j1 = (int)round(2 * j1);
  int two_j2 = (int)round(2 * j2);
  int two_m1 = (int)round(2 * m1);
  int two_m2 = (int)round(2 * m2);
  int two_J = (int)round(2 * J);
  int two_M = (int)round(2 * M);
  int sign = -1;
  if ((two_j1 - two_j2 + two_M) % 4 == 0)
    sign = 1; // mod 4 (instead 2), since x2
  return sign * sqrt(two_J + 1.) *
         gsl_sf_coupling_3j(two_j1, two_j2, two_J, two_m1, two_m2, -two_M);
}

//------------------------------------------------------------------------------
inline double cg_1(int j1, int m1, int j2, int m2, int J, int M)
// Calculates Clebsh-Gordon coeficient:
// <j1 m1, j2 m2 | J M> = (-1)^(j1-j2+M) * Sqrt(2J+1) * (j1 j2  J)
// .                                                    (m1 m2 -M)
// (Last term is 3j symbol)
// Note: this function takes INTEGER values, only works for l (not half-integer
// j)!
{
  if (triangle(j1, j2, J) * sumsToZero(m1, m2, -M) == 0)
    return 0;
  int sign = -1;
  if ((j1 - j2 + M) % 2 == 0)
    sign = 1;
  return sign * sqrt(2. * J + 1.) *
         gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * J, 2 * m1, 2 * m2, -2 * M);
}

//------------------------------------------------------------------------------
inline double cg_2(int two_j1, int two_m1, int two_j2, int two_m2, int two_J,
                   int two_M)
// <j1 m1, j2 m2 | J M> = (-1)^(j1-j2+M) * Sqrt(2J+1) * (j1 j2  J)
// .                                                    (m1 m2 -M)
// (Last term is 3j symbol)
// Note: this function takes INTEGER values, that have already multiplied by 2!
// Works for l and j (integer and half-integer)
{
  if (triangle(two_j1, two_j2, two_J) * sumsToZero(two_m1, two_m2, -two_M) == 0)
    return 0;
  int sign = -1;
  if ((two_j1 - two_j2 + two_M) % 4 == 0)
    sign = 1; // mod 4 (instead 2), since x2
  return sign * sqrt(two_J + 1.) *
         gsl_sf_coupling_3j(two_j1, two_j2, two_J, two_m1, two_m2, -two_M);
}

//******************************************************************************
inline double sixj(double j1, double j2, double j3, double j4, double j5,
                   double j6)
// Calculates wigner 6j symbol:
//   {j1 j2 j3}
//   {j4 j5 j6}
// Note: this function takes DOUBLE values.
// Works for l and j (integer and half-integer)
{

  if (triangle(j1, j2, j3) * triangle(j1, j5, j6) * triangle(j4, j2, j6) *
          triangle(j4, j5, j3) ==
      0)
    return 0;

  int two_j1 = (int)round(2 * j1);
  int two_j2 = (int)round(2 * j2);
  int two_j3 = (int)round(2 * j3);
  int two_j4 = (int)round(2 * j4);
  int two_j5 = (int)round(2 * j5);
  int two_j6 = (int)round(2 * j6);
  return gsl_sf_coupling_6j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6);
}

//------------------------------------------------------------------------------
inline double sixj_1(int j1, int j2, int j3, int j4, int j5, int j6)
// Calculates wigner 6j symbol:
//   {j1 j2 j3}
//   {j4 j5 j6}
// Note: this function takes INTEGER values, only works for l (not half-integer
// j)!
{
  if (triangle(j1, j2, j3) * triangle(j1, j5, j6) * triangle(j4, j2, j6) *
          triangle(j4, j5, j3) ==
      0)
    return 0;
  return gsl_sf_coupling_6j(2 * j1, 2 * j2, 2 * j3, 2 * j4, 2 * j5, 2 * j6);
}

//------------------------------------------------------------------------------
inline double sixj_2(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5,
                     int two_j6)
// Calculates wigner 6j symbol:
//   {j1 j2 j3}
//   {j4 j5 j6}
// Note: this function takes INTEGER values, that have already multiplied by 2!
// Works for l and j (integer and half-integer)
{
  if (triangle(two_j1, two_j2, two_j3) * triangle(two_j1, two_j5, two_j6) *
          triangle(two_j4, two_j2, two_j6) * triangle(two_j4, two_j5, two_j3) ==
      0)
    return 0;
  return gsl_sf_coupling_6j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6);
}

//******************************************************************************
inline double ninej(double j1, double j2, double j3, double j4, double j5,
                    double j6, double j7, double j8, double j9)
// Calculates wigner 9j symbol:
//   {j1 j2 j3}
//   {j4 j5 j6}
//   {j7 j8 j9}
// Note: this function takes DOUBLE values.
// Works for l and j (integer and half-integer)
{
  int two_j1 = (int)round(2 * j1);
  int two_j2 = (int)round(2 * j2);
  int two_j3 = (int)round(2 * j3);
  int two_j4 = (int)round(2 * j4);
  int two_j5 = (int)round(2 * j5);
  int two_j6 = (int)round(2 * j6);
  int two_j7 = (int)round(2 * j7);
  int two_j8 = (int)round(2 * j8);
  int two_j9 = (int)round(2 * j9);
  return gsl_sf_coupling_9j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6,
                            two_j7, two_j8, two_j9);
}

//------------------------------------------------------------------------------
inline double ninej_1(int j1, int j2, int j3, int j4, int j5, int j6, int j7,
                      int j8, int j9)
// Calculates wigner 9j symbol:
//   {j1 j2 j3}
//   {j4 j5 j6}
//   {j7 j8 j9}
// Note: this function takes INTEGER values, only works for l (not half-integer
// j)!
{
  return gsl_sf_coupling_9j(2 * j1, 2 * j2, 2 * j3, 2 * j4, 2 * j5, 2 * j6,
                            2 * j7, 2 * j8, 2 * j9);
}

//------------------------------------------------------------------------------
inline double ninej_2(int two_j1, int two_j2, int two_j3, int two_j4,
                      int two_j5, int two_j6, int two_j7, int two_j8,
                      int two_j9)
// Calculates wigner 9j symbol:
//   {j1 j2 j3}
//   {j4 j5 j6}
//   {j7 j8 j9}
// Note: this function takes INTEGER values, that have already multiplied by 2!
// Works for l and j (integer and half-integer)
{
  return gsl_sf_coupling_9j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6,
                            two_j7, two_j8, two_j9);
}

//******************************************************************************
inline double Ck_kk(int k, int ka, int kb)
// Reduced (relativistic) angular ME:
// <ka||C^k||kb> = (-1)^(ja+1/2) * 3js(ja jb k, -1/2 1/2 0) * Pi
// Note: takes in kappa! (not j!)
{
  if (parity(l_k(ka), l_k(kb), k) == 0) {
    return 0;
  }
  auto two_ja = twoj_k(ka);
  auto two_jb = twoj_k(kb);
  auto sign = ((two_ja + 1) / 2 % 2 == 0) ? 1 : -1;
  auto f = sqrt((two_ja + 1) * (two_jb + 1));
  auto g = gsl_sf_coupling_3j(two_ja, two_jb, 2 * k, -1, 1, 0);
  return sign * f * g;
}

} // namespace Wigner
