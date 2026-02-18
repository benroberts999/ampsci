#pragma once
#include <algorithm> //std::min!
#include <cmath>
#include <gsl/gsl_sf_coupling.h>
#include <optional>
#include <utility>

/*!
@brief
Angular provides functions and classes for calculating and storing angular
factors (3,6,9-J symbols etc.)

@details
Provides functions to:
 - Calculate wigner 3,6,9-J symbols + Clebsh-Gordon coefs etc..
 - Also provied look-up tables, which are faster than calculating symbols
on-the-fly
 - Wrapper functions to calculate wigner 3,6,9-J symbols.
 - Uses GSL:
https://www.gnu.org/software/gsl/doc/html/specfunc.html?highlight=3j#coupling-coefficients

The general equations are defined:
\f{align}{
    \frac{1}{r_{12}} &= \sum_{kq} \frac{r_<^k}{r_>^{k+1}}(-1)^q
                        C^k_{-q}(\hat{n}_1)C^k_{q}(\hat{n}_2)\\
    C^k_{q} &\equiv \sqrt{\frac{4\pi}{2k+1}} Y_{kq}(\hat{n}),
\f}
and
\f{align}{
  C^k_{ab} &\equiv \langle{\kappa_a}||C^k||{\kappa_b}\rangle 
            \equiv  (-1)^{j_a+1/2} \widetilde C^k_{ab} ,\\
            &=  (-1)^{j_a+1/2}\sqrt{[j_a][j_b]}
            \begin{pmatrix}
              {j_a} & {j_b} & {k} \\ {-1/2} & {1/2} &{0}
            \end{pmatrix}
            \pi(l_a+l_b+k). 
\f}
*/
namespace Angular {

//==============================================================================
//! returns l given kappa
constexpr int l_k(int ka) { return (ka > 0) ? ka : -ka - 1; }
//! returns 2j given kappa
constexpr int twoj_k(int ka) { return (ka > 0) ? 2 * ka - 1 : -2 * ka - 1; }
//! returns parity [(-1)^l] given kappa
constexpr int parity_k(int ka) {
  return (ka % 2 == 0) ? ((ka > 0) ? 1 : -1) : ((ka < 0) ? 1 : -1);
}
//! returns parity [(-1)^l] given l
constexpr int parity_l(int l) { return (l % 2 == 0) ? 1 : -1; }
//! "Complimentary" l := (2j-l) = l+/-1, for  j=l+/-1/2
constexpr int l_tilde_k(int ka) { return (ka > 0) ? ka - 1 : -ka; }
//! returns kappa, given 2j and l
constexpr int kappa_twojl(int twoj, int l) {
  return ((2 * l - twoj) * (twoj + 1)) / 2;
}
//! returns kappa, given 2*j and parity (+/-1),
constexpr int kappa_twojpi(const int twoj, const int pi) {
  const auto l_minus = (twoj - 1) / 2;
  const auto l_minus_pi = (l_minus % 2 == 0) ? 1 : -1;
  const auto l = l_minus_pi == pi ? l_minus : l_minus + 1;
  return kappa_twojl(twoj, l); // must be simpler way?
}
//! Checks if a double is within eps of 0.0 [eps=1.0e-10 by default]
inline constexpr bool zeroQ(double x, double eps = 1.0e-10) {
  return x >= 0 ? x < eps : x > -eps;
}

//==============================================================================

//! @brief returns kappa_index, given kappa
/*! @details
 Kappa Index:
 For easy array access, define 1-to-1 index for each kappa:
 kappa: -1  1 -2  2 -3  3 -4  4 ...
 index:  0  1  2  3  4  5  6  7 ...
 kappa(i) = (-1,i+1)*(int(i/2)+1)
*/
constexpr int indexFromKappa(int ka) {
  return (ka < 0) ? -2 * ka - 2 : 2 * ka - 1;
}
//! Returns kappa, given kappa_index
constexpr int kappaFromIndex(int i) {
  return (i % 2 == 0) ? -(i + 2) / 2 : (i + 1) / 2;
}
//! returns 2j, given kappa_index
constexpr int twojFromIndex(int i) { return (i % 2 == 0) ? i + 1 : i; }
//! returns l, given kappa_index
constexpr int lFromIndex(int i) { return (i % 2 == 0) ? i / 2 : (i + 1) / 2; }

//==============================================================================
//! Returns number of possible states _below_ given n
constexpr int states_below_n(int n) { return n * n - 2 * n + 1; }

//! return nk_index given {n, kappa}: nk_index(n,k) := n^2 - 2n + 1 +
//! kappa_index
/*! @details   nk_index:
 For easy array access, define 1-to-1 index for each {n, kappa}:
 nk_index(n,k) := n^2 - 2n + 1 + kappa_index.
 nb: n^2 - 2n + 1 = states_below_n - number of possible states with n'<n.
 Note: ONLY valid for n >= 1 (i.e., cannot be used for general basis states)
*/
constexpr int nk_to_index(int n, int k) {
  return states_below_n(n) + Angular::indexFromKappa(k);
}

//! return {n, kappa} given nk_index:
inline std::pair<int, int> index_to_nk(int index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(index + 0.01));
  // int n = 1 + int_sqrt(index);
  const auto kappa_index = index - states_below_n(n);
  return {n, Angular::kappaFromIndex(kappa_index)};
}

inline std::pair<int, int> nkindex_to_n_kindex(int index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(index + 0.01));
  // int n = 1 + int_sqrt(index);
  const auto kappa_index = index - states_below_n(n);
  return {n, kappa_index};
}

//! Returns kappa, given nk_index
inline int nkindex_to_kappa(int index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(index + 0.01));
  // int n = 1 + int_sqrt(index);
  const auto kappa_index = index - states_below_n(n);
  return kappaFromIndex(kappa_index);
}

//! Returns 2*j, given nk_index
inline int nkindex_to_twoj(int index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(index + 0.01));
  // int n = 1 + int_sqrt(index);
  const auto kappa_index = index - states_below_n(n);
  return twojFromIndex(kappa_index);
}

//! Returns l, given nk_index
inline int nkindex_to_l(int index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(index + 0.01));
  // int n = 1 + int_sqrt(index);
  const auto kappa_index = index - states_below_n(n);
  return l_k(kappaFromIndex(kappa_index));
}

//==============================================================================
//! Returns true if a is even - for integer values
constexpr bool evenQ(int a) { return (a % 2 == 0); }
//! Returns true if a (an int) is even, given 2*a (true if two_a/2 is even)
constexpr bool evenQ_2(int two_a) { return (two_a % 4 == 0); }
//! Evaluates (-1)^{a} (for integer a)
constexpr int neg1pow(int a) { return evenQ(a) ? 1 : -1; }
//! Evaluates (-1)^{two_a/2} (for integer a; two_a is even)
constexpr int neg1pow_2(int two_a) { return evenQ_2(two_a) ? 1 : -1; }

//==============================================================================
//! Parity rule. Returns 1 only if la+lb+k is even
constexpr int parity(int la, int lb, int k) {
  return evenQ(la + lb + k) ? 1 : 0;
}

//==============================================================================
//! Returns 1 if triangle rule is satisfied. nb: works with j OR twoj!
constexpr int triangle(int j1, int j2, int J) {
  // nb: can be called with wither j or twoj!
  return ((j1 + j2 < J) || (std::abs(j1 - j2) > J)) ? 0 : 1;
}
//! Checks if three ints sum to zero, returns 1 if they do
constexpr int sumsToZero(int m1, int m2, int m3) {
  return (m1 + m2 + m3 != 0) ? 0 : 1;
}

//==============================================================================
//! @brief Calculates wigner 3j symbol. Takes in 2*j (or 2*l) - intput is
//! integer
/*! @details
   \f[ \begin{pmatrix}j1&j2&j3\\m1&m2&m3\end{pmatrix} \f]
*/
inline double threej_2(int two_j1, int two_j2, int two_j3, int two_m1,
                       int two_m2, int two_m3) {
  if (triangle(two_j1, two_j2, two_j3) * sumsToZero(two_m1, two_m2, two_m3) ==
      0)
    return 0;
  return gsl_sf_coupling_3j(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

//==============================================================================
//!@brief  Special (common) 3js case:  (j1 j2 k, -0.5, 0.5, 0)
/*! @details
   \f[ \begin{pmatrix}j1&j2&k\\-1/2&1/2&0\end{pmatrix} \f]
*/
inline double special_threej_2(int two_j1, int two_j2, int two_k) {
  if (triangle(two_j1, two_j2, two_k) == 0)
    return 0.0;
  if (two_k == 0) {
    auto s = evenQ_2(two_j1 + 1) ? 1.0 : -1.0;
    return s / std::sqrt(two_j1 + 1);
  }
  return gsl_sf_coupling_3j(two_j1, two_j2, two_k, -1, 1, 0);
}

//==============================================================================
//! Clebsh-Gordon coeficient <j1 m1, j2 m2 | J M> [takes 2*j, as int]
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

//==============================================================================
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

//! Checks if a 6j symbol is valid - each input is optional
//! @details i.e., sixjTriads(a,b,c,{},{},{}) will check if _any_ 6J symbol of
//! form {a b c \ * * *} is valid (* can be any value)
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

//==============================================================================
//! 6j symbol {j1 j2 j3 \\ j4 j5 j6} - [takes 2*j as int]
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

//------------------------------------------------------------------------------
//! 9j symbol {j1 j2 j3 \\ j4 j5 j6 \\ j7 j8 j9} [takes 2*j as int]
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

//==============================================================================
//! Reduced (relativistic) angular ME: <ka||C^k||kb> [takes k and kappa]
inline double Ck_kk(int k, int ka, int kb)
// Reduced (relativistic) angular ME:
// <ka||C^k||kb> = (-1)^(ja+1/2) * srt([ja][jb]) * 3js(ja jb k, -1/2 1/2 0) * Pi
// Note: takes in kappa! (not j!)
{
  if (parity(l_k(ka), l_k(kb), k) == 0) {
    return 0.0;
  }
  auto two_ja = twoj_k(ka);
  auto two_jb = twoj_k(kb);
  // auto sign = ((two_ja + 1) / 2 % 2 == 0) ? 1 : -1;
  auto sign = evenQ_2(two_ja + 1) ? 1 : -1;
  auto f = std::sqrt((two_ja + 1) * (two_jb + 1));
  auto g = special_threej_2(two_ja, two_jb, 2 * k);
  return sign * f * g;
}

//! Full C^k_q matrix element - rarely used
inline double Ck_kk_mmq(int k, int ka, int kb, int twoma, int twomb, int twoq) {
  const auto tja = twoj_k(ka);
  const auto tjb = twoj_k(kb);
  return neg1pow_2(tja - twoma) *
         threej_2(tja, 2 * k, tjb, -twoma, twoq, twomb) * Ck_kk(k, ka, kb);
}

//! Ck selection rule only. Returns false if C^k=0, true if non-zero.
inline bool Ck_kk_SR(int k, int ka, int kb) {
  if (parity(l_k(ka), l_k(kb), k) == 0)
    return false;
  if (triangle(twoj_k(ka), twoj_k(kb), 2 * k) == 0)
    return false;
  return true;
}

//! tildeCk_kk = (-1)^{ja+1/2}*Ck_kk
inline double tildeCk_kk(int k, int ka, int kb) {
  // tildeCk_kk = (-1)^{ja+1/2}*Ck_kk
  auto m1tjph = evenQ_2(twoj_k(ka) + 1) ? 1 : -1;
  return m1tjph * Ck_kk(k, ka, kb);
}

//! Returns [k_min, k_kmax] for C^k (min/max non-zero k, given kappa_a, kappa_b)
//! Accounts for parity
inline std::pair<int, int> kminmax_Ck(int ka, int kb) {
  // j = |k|-0.5
  // kmin = |ja-jb| = | |ka| - |kb| |
  // kmax = ja+jb   = |ka| + |kb| - 1
  const auto aka = std::abs(ka);
  const auto akb = std::abs(kb);

  auto kmin = std::abs(aka - akb);
  auto kmax = aka + akb - 1;
  // account for parity
  if (parity(l_k(ka), l_k(kb), kmin) == 0)
    ++kmin;
  if (parity(l_k(ka), l_k(kb), kmax) == 0)
    --kmax;

  return {kmin, kmax};
}

//==============================================================================
//! @brief Reduced spin angular ME: (for spin 1/2): <ka||S||kb>
/*! @details
  Reduced spin angular ME: (for spin 1/2!)
  <ka||S||kb> = d(la,lb) * (-1)^{ja+la+3/2} * Sqrt([ja][jb](3/2)) *
              *  sjs{ja, 1, jb,  1/2, la, 1/2}
  Special 6j case:
  sjs{ja, 1, jb,  1/2, la, 1/2}
    = 0.5 (-1)^{ka+kb} * Sqrt(abs[{(ka-1)^2 - kb^2}/{3ka(1+2ka)}])
      * triangle rule for j!
  At least ~least 20% faster
*/
inline double S_kk(int ka, int kb) {
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
