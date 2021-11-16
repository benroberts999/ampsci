#pragma once
#include <algorithm> //std::min!
#include <cmath>
#include <gsl/gsl_sf_coupling.h>
#include <optional>
#include <utility>

/*!
@brief
Calculate wigner 3,6,9-J symbols + Clebsh-Gordon coefs etc..
@details
Wrapper functions to calculate wigner 3,6,9-J symbols.
Uses GSL:
https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=3j#coupling-coefficients
NOTE:
Since j always integer or half-integer, GSL inputs are always integer.
Three versions of each symbol:
 - 'regular', takes in double. Converts to integer safely. Slower (marginally),
    but easier
 - '_1' version takes integers as-is. Only work for integer angular momentum
(l).
 - and '_2' version, takes in 2*j (as an integer). Works for l and j
*/
namespace Angular {

//******************************************************************************
//! @brief returns l given kappa (all folliwing do similar)
constexpr int l_k(int ka) { return (ka > 0) ? ka : -ka - 1; }
constexpr int twoj_k(int ka) { return (ka > 0) ? 2 * ka - 1 : -2 * ka - 1; }
constexpr double j_k(int ka) {
  return (ka > 0) ? double(ka) - 0.5 : double(-ka) - 0.5;
}
constexpr int parity_k(int ka) {
  return (ka % 2 == 0) ? ((ka > 0) ? 1 : -1) : ((ka < 0) ? 1 : -1);
}
constexpr int parity_l(int l) { return (l % 2 == 0) ? 1 : -1; }
//! @brief "Complimentary" l (l for lower component). l-tilde = (2j-l) = l +/-
//! 1, for j = l +/- 1/2
constexpr int l_tilde_k(int ka) { return (ka > 0) ? ka - 1 : -ka; }
constexpr int kappa_twojl(int twoj, int l) {
  return ((2 * l - twoj) * (twoj + 1)) / 2;
}
//! Given 2*j and parity (+/-1), returns kappa
constexpr int kappa_twojpi(const int twoj, const int pi) {
  const auto l_minus = (twoj - 1) / 2;
  const auto l_minus_pi = (l_minus % 2 == 0) ? 1 : -1;
  const auto l = l_minus_pi == pi ? l_minus : l_minus + 1;
  return kappa_twojl(twoj, l); // must be simpler way?
}

inline constexpr bool zeroQ(double x, double eps = 1.0e-10) {
  return x >= 0 ? x < eps : x > -eps;
}

//******************************************************************************
//    Kappa Index:
// For easy array access, define 1-to-1 index for each kappa:
// kappa: -1  1 -2  2 -3  3 -4  4 ...
// index:  0  1  2  3  4  5  6  7 ...
// kappa(i) = (-1,i+1)*(int(i/2)+1)
//! @brief kappa index. kappa=-1,1,-2,...; index=0,1,2.
constexpr int indexFromKappa(int ka) {
  return (ka < 0) ? -2 * ka - 2 : 2 * ka - 1;
}
constexpr int kappaFromIndex(int i) {
  return (i % 2 == 0) ? -(i + 2) / 2 : (i + 1) / 2;
}
constexpr int twojFromIndex(int i) { return (i % 2 == 0) ? i + 1 : i; }
constexpr int lFromIndex(int i) { return (i % 2 == 0) ? i / 2 : (i + 1) / 2; }

//******************************************************************************

//! @brief Minimum/Maximum 'l' allowed in {a,b,k \ c,d,l} 6j symbol
inline int min_lambda_tj(int tja, int tjb, int tjc, int tjd) {
  return std::max(std::abs(tja - tjd), std::abs(tjb - tjc)) / 2;
}
//! @brief Minimum/Maximum 'l' allowed in {a,b,k \ c,d,l} 6j symbol
inline int max_lambda_tj(int tja, int tjb, int tjc, int tjd) {
  return std::min((tja + tjd), (tjb + tjc)) / 2;
}
//******************************************************************************
//! @brief Returns true if a is even
constexpr bool evenQ(int a) { return (a % 2 == 0); }
//! @brief Returns true if a is even, given 2*a (i.e., true if two_a/2 is even)
constexpr bool evenQ_2(int two_a) { return (two_a % 4 == 0); }
//! Evaluates (-1)^{a}
constexpr int neg1pow(int a) { return evenQ(a) ? 1 : -1; }
//! Evaluates (-1)^{two_a/2}
constexpr int neg1pow_2(int two_a) { return evenQ_2(two_a) ? 1 : -1; }

//******************************************************************************
//! @brief Parity rule. Returns 1 only if la+lb+k is even
constexpr int parity(int la, int lb, int k) {
  return ((la + lb + k) % 2 == 0) ? 1 : 0;
}

//******************************************************************************
//! @brief Returns 1 if triangle rule is satisfied
constexpr int triangle(double j1, double j2, double J) {
  return ((j1 + j2 < J) || (std::fabs(j1 - j2) > J)) ? 0 : 1;
}
//! @brief Returns 1 if triangle rule is satisfied. nb: works with j OR twoj!
constexpr int triangle(int j1, int j2, int J) {
  // nb: can be called with wither j or twoj!
  return ((j1 + j2 < J) || (std::abs(j1 - j2) > J)) ? 0 : 1;
}

constexpr int sumsToZero(int m1, int m2, int m3) {
  return (m1 + m2 + m3 != 0) ? 0 : 1;
}

constexpr int sumsToZero(double m1, double m2, double m3) {
  return ((m1 + m2 + m3) > 0.00001 || (m1 + m2 + m3) < -0.00001) ? 0 : 1;
}

//******************************************************************************
//! @brief Calculates wigner 3j symbol: Works for l and j (integer and
//! half-integer)
/*! @details
   \f[ \begin{pmatrix}j1&j2&j3\\m1&m2&m3\end{pmatrix} \f]
*/
inline double threej(double j1, double j2, double j3, double m1, double m2,
                     double m3) {
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
//! @brief Calculates wigner 3j symbol: only works for l (not half-integer j)!
{
  if (triangle(j1, j2, j3) * sumsToZero(m1, m2, m3) == 0)
    return 0;
  return gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3);
}

//------------------------------------------------------------------------------
//! @brief Calculates wigner 3j symbol: takes INTEGER values, that have already
//! multiplied by 2. Works for l and j (integer and half-integer).
inline double threej_2(int two_j1, int two_j2, int two_j3, int two_m1,
                       int two_m2, int two_m3) {
  if (triangle(two_j1, two_j2, two_j3) * sumsToZero(two_m1, two_m2, two_m3) ==
      0)
    return 0;
  return gsl_sf_coupling_3j(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

//******************************************************************************
//!@brief  Special (common) 3js case:  (j1 j2 k, -0.5, 0.5, 0)
/*! @details
   \f[ \begin{pmatrix}j1&j2&k\\-1/2&1/2&0\end{pmatrix} \f]
*/
inline double special_threej_2(int two_j1, int two_j2, int two_k) {
  if (triangle(two_j1, two_j2, two_k) == 0)
    return 0.0;
  if (two_k == 0) {
    // auto s = ((two_j1 + 1) % 4 == 0) ? 1.0 : -1.0;
    auto s = evenQ_2(two_j1 + 1) ? 1.0 : -1.0;
    return s / std::sqrt(two_j1 + 1);
  }
  // else if(two_k == 1){
  // XXX Simple formula??
  // }
  return gsl_sf_coupling_3j(two_j1, two_j2, two_k, -1, 1, 0);
}

//******************************************************************************
//!@brief Clebsh-Gordon coeficient <j1 m1, j2 m2 | J M>
inline double cg(double j1, double m1, double j2, double m2, double J, double M)
// <j1 m1, j2 m2 | J M> = (-1)^(j1-j2+M) * std::sqrt(2J+1) * (j1 j2  J)
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
  return sign * std::sqrt(two_J + 1.) *
         gsl_sf_coupling_3j(two_j1, two_j2, two_J, two_m1, two_m2, -two_M);
}

//------------------------------------------------------------------------------
inline double cg_1(int j1, int m1, int j2, int m2, int J, int M)
// Calculates Clebsh-Gordon coeficient:
// <j1 m1, j2 m2 | J M> = (-1)^(j1-j2+M) * std::sqrt(2J+1) * (j1 j2  J)
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
  return sign * std::sqrt(2. * J + 1.) *
         gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * J, 2 * m1, 2 * m2, -2 * M);
}

//------------------------------------------------------------------------------
inline double cg_2(int two_j1, int two_m1, int two_j2, int two_m2, int two_J,
                   int two_M)
// <j1 m1, j2 m2 | J M> = (-1)^(j1-j2+M) * std::sqrt(2J+1) * (j1 j2  J)
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
  return sign * std::sqrt(two_J + 1.) *
         gsl_sf_coupling_3j(two_j1, two_j2, two_J, two_m1, two_m2, -two_M);
}

//******************************************************************************
//! Checks triangle conditions for 6j symbols (for 2*j)
inline bool sixj_zeroQ(int a, int b, int c, int d, int e, int f) {
  // nb: works for 2*integers ONLY
  // check triangle consitions
  // Note: there are some  6j symbols that pass this test, though are still zero
  // GSL calculates these to have extremely small (but non-zero) value
  if (triangle(a, b, c) == 0)
    return true;
  if (triangle(c, d, e) == 0)
    return true;
  if (triangle(b, d, f) == 0)
    return true;
  if (triangle(e, f, a) == 0)
    return true;
  if (!evenQ(a + b + c))
    return true;
  if (!evenQ(c + d + e))
    return true;
  if (!evenQ(b + d + f))
    return true;
  if (!evenQ(e + f + a))
    return true;
  return false;
}

inline bool sixjTriads(std::optional<int> a, std::optional<int> b,
                       std::optional<int> c, std::optional<int> d,
                       std::optional<int> e, std::optional<int> f) {
  // nb: works for 2*integers ONLY
  // checks triangle consitions, with optional parameters
  if (a && b && c && triangle(*a, *b, *c) == 0)
    return false;
  if (c && d && e && triangle(*c, *d, *e) == 0)
    return false;
  if (a && e && f && triangle(*a, *e, *f) == 0)
    return false;
  if (b && d && f && triangle(*b, *d, *f) == 0)
    return false;

  return true;
}

//******************************************************************************
//!@brief 6j symbol {j1 j2 j3 \\ j4 j5 j6}
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
  if (sixj_zeroQ(j1, j2, j3, j4, j5, j6))
    return 0.0;
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
  if (sixj_zeroQ(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6))
    return 0.0;
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
//!@brief Reduced (relativistic) angular ME: <ka||C^k||kb>
inline double Ck_kk(int k, int ka, int kb)
// Reduced (relativistic) angular ME:
// <ka||C^k||kb> = (-1)^(ja+1/2) * srt([ja][jb]) * 3js(ja jb k, -1/2 1/2 0) * Pi
// Note: takes in kappa! (not j!)
{
  if (parity(l_k(ka), l_k(kb), k) == 0) {
    return 0;
  }
  auto two_ja = twoj_k(ka);
  auto two_jb = twoj_k(kb);
  // auto sign = ((two_ja + 1) / 2 % 2 == 0) ? 1 : -1;
  auto sign = evenQ_2(two_ja + 1) ? 1 : -1;
  auto f = std::sqrt((two_ja + 1) * (two_jb + 1));
  auto g = special_threej_2(two_ja, two_jb, 2 * k);
  // auto g = gsl_sf_coupling_3j(two_ja, two_jb, 2 * k, -1, 1, 0);
  return sign * f * g;
}

//! Ck selection rule only. Returns false if C^k=0, true if non-zero.
inline bool Ck_kk_SR(int k, int ka, int kb) {
  if (parity(l_k(ka), l_k(kb), k) == 0)
    return false;
  if (triangle(twoj_k(ka), twoj_k(kb), 2 * k) == 0)
    return false;
  return true;
}

//! tilde version: symmetric
inline double tildeCk_kk(int k, int ka, int kb) {
  // tildeCk_kk = (-1)^{ja+1/2}*Ck_kk
  auto m1tjph = evenQ_2(twoj_k(ka) + 1) ? 1 : -1;
  return m1tjph * Ck_kk(k, ka, kb);
}

//! Returns [k_min, k_kmax] for C^k (min/max non-zero k, given kappa_a, kappa_b)
inline std::pair<int, int> kminmax_Ck(int ka, int kb) {
  // j = |k|-0.5
  // kmin = |ja-jb| = | |ka| - |kb| |
  // kmax = ja+jb   = |ka| + |kb| - 1
  const auto aka = std::abs(ka);
  const auto akb = std::abs(kb);
  return {std::abs(aka - akb), aka + akb - 1};

  // return std::make_pair(std::abs(a.twoj() - b.twoj()) / 2,
  //                       (a.twoj() + b.twoj()) / 2);
}
//******************************************************************************
inline double Ck_2j2j(int k, int two_ja, int two_jb)
// Reduced (relativistic) angular ME:
// <ka||C^k||kb> = (-1)^(ja+1/2) * srt([ja][jb]) * 3js(ja jb k, -1/2 1/2 0) * Pi
// Note: takes in two*j!
// NOTE: DOESNT check parity! Only use if that's already known to be true
{
  // auto sign = (((two_ja + 1) / 2) % 2 == 0) ? 1 : -1;
  auto sign = evenQ_2(two_ja + 1) ? 1 : -1;
  auto f = std::sqrt((two_ja + 1) * (two_jb + 1));
  // auto g = gsl_sf_coupling_3j(two_ja, two_jb, 2 * k, -1, 1, 0);
  auto g = special_threej_2(two_ja, two_jb, 2 * k);
  return sign * f * g;
}
//******************************************************************************
//!@brief Reduced spin angular ME: (for spin 1/2): <ka||S||kb>
inline double S_kk(int ka, int kb)
// Reduced spin angular ME: (for spin 1/2!)
// <ka||S||kb> = d(la,lb) * (-1)^{ja+la+3/2} * Sqrt([ja][jb](3/2)) *
//             *  sjs{ja, 1, jb,  1/2, la, 1/2}
// Special 6j case:
// sjs{ja, 1, jb,  1/2, la, 1/2}
//    = 0.5 (-1)^{ka+kb} * Sqrt(abs[{(ka-1)^2 - kb^2}/{3ka(1+2ka)}])
//     * triangle rule for j!
// At least ~least 20% faster
{
  auto la = l_k(ka);
  if (la != l_k(kb))
    return 0.0;
  auto tja = twoj_k(ka);
  auto tjb = twoj_k(kb);
  if (triangle(tja, 2, tjb) == 0)
    return 0;
  auto sign = (((tja + 2 * la + 3) / 2) % 2 == 0) ? 1 : -1;
  auto f = std::sqrt((tja + 1) * (tjb + 1) * 1.5);
  // auto sixj = gsl_sf_coupling_6j(tja, 2, tjb, 1, 2 * la, 1);
  auto sixj_sign = ((ka + kb) % 2 == 0) ? 1.0 : -1.0;
  auto sixj = 0.5 * sixj_sign *
              std::sqrt(std::fabs(((ka - 1) * (ka - 1) - kb * kb) /
                                  (3.0 * ka * (1.0 + 2.0 * ka))));
  return sign * f * sixj;
}

} // namespace Angular
