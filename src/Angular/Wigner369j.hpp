#pragma once
#include <algorithm> //std::min!
#include <cmath>
#include <cstdint>
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
  https://www.gnu.org/software/gsl/doc/html/specfunc.html#coupling-coefficients

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

  | kappa | -1 | 1 | -2 | 2 | -3 | 3 | -4 | 4 | ... |
  |-------|----|---|----|---|----|---|----|---|-----|
  | index | 0  | 1 | 2  | 3 | 4  | 5 | 6  | 7 | ... |
 
  \f[
    i_k(\kappa)=
    \begin{cases}
      -2\kappa-2, & \kappa<0, \\
      2\kappa-1,  & \kappa>0 .
    \end{cases}
  \f]

  \f[
    \kappa(i_k)=
    \begin{cases}
      -\left(\frac{i_k}{2}+1\right), & i_k\ \text{even},\\
      \frac{i_k+1}{2}, & i_k\ \text{odd}.
    \end{cases}
  \f]
*/
constexpr std::uint64_t kappa_to_kindex(int ka) {
  return (ka < 0) ? std::uint64_t(-2 * ka - 2) : std::uint64_t(2 * ka - 1);
}

//! Returns kappa, given kappa_index
//! @details see \ref Angular::kappa_to_kindex()
constexpr int kindex_to_kappa(std::uint64_t i) {
  return (i % 2 == 0) ? int(-(i + 2) / 2) : int((i + 1) / 2);
}

//! returns 2j, given kappa_index
//! @details see \ref Angular::kappa_to_kindex()
constexpr int kindex_to_twoj(std::uint64_t i) {
  return (i % 2 == 0) ? int(i + 1) : int(i);
}

//! returns l, given kappa_index
//! @details see \ref Angular::kappa_to_kindex()
constexpr int kindex_to_l(std::uint64_t i) {
  return (i % 2 == 0) ? int(i / 2) : int((i + 1) / 2);
}

//==============================================================================
//! Returns number of possible states _below_ given n (i.e., n' < n)
/*! @details

  This could be made more "compact" with a maximum l.
  In practice, we go to large n, but never very large l.

  @warning
   - n must be >1 (no non-standard states)

  Cannot overflow; converts to std::uint64_t.
*/
constexpr std::uint64_t states_below_n(int n) {
  const std::uint64_t nn = static_cast<std::uint64_t>(n);
  return nn * nn - 2 * nn + 1;
}

//! @brief Returns nk_index, given {n, kappa}
/*! @details
  For convenient array access we define a one-to-one mapping between the pair
  \f$(n,\kappa)\f$ and a non-negative index:

  \f[
    i_{nk}(n,\kappa) = (n-1)^2 + i_k(\kappa),
  \f]

  where \f$i_k(\kappa)\f$ is the kappa index defined by \ref Angular::kappa_to_kindex()

  Equivalently,

  \f[
    i_{nk}(n,\kappa) = n^2 - 2n + 1 + i_k(\kappa).
  \f]

  Since

  \f[
    n^2 - 2n + 1 = (n-1)^2 = \texttt{states\_below\_n}(n),
  \f]
  
  this counts the number of states with principal quantum number \f$n' < n\f$,
  plus the index within the \f$n\f$ shell.

  @note This indexing is only valid for physical bound-state quantum numbers, n >= 1,
  and cannot in general be used for arbitrary basis states.

  @warning To safely convert the returned index (`uint64_t`) to `int`, one must
  have `n <= 46340`.

  @warning To safely convert the returned index to \ref Coulomb::nkIndex
  (`uint16_t`), one must have `n <= 256`.
*/
constexpr std::uint64_t nk_to_index(int n, int k) {
  return states_below_n(n) + Angular::kappa_to_kindex(k);
}

//! Returns {n, kappa} given nk_index
//! @details Inverse of \ref Angular::nk_to_index().
inline std::pair<int, int> index_to_nk(std::uint64_t nk_index) {
  // Better way? isqrt?
  const auto n = 1 + int(std::sqrt(nk_index));
  const auto kappa_index = nk_index - states_below_n(n);
  return {n, Angular::kindex_to_kappa(kappa_index)};
}

//! Returns {n, kappa_index} given nk_index
//! @details Inverse of \ref Angular::nk_to_index().
inline std::pair<int, std::uint64_t>
nkindex_to_n_kindex(std::uint64_t nk_index) {
  const auto n = 1 + int(std::sqrt(nk_index));
  const auto kappa_index = nk_index - states_below_n(n);
  return {n, kappa_index};
}

//! Returns kappa, given nk_index
//! @details See \ref Angular::nk_to_index().
inline int nkindex_to_kappa(std::uint64_t nk_index) {
  const auto n = 1 + int(std::sqrt(nk_index));
  const auto kappa_index = nk_index - states_below_n(n);
  return kindex_to_kappa(kappa_index);
}

//! Returns 2*j, given nk_index
//! @details See \ref Angular::nk_to_index().
inline int nkindex_to_twoj(std::uint64_t nk_index) {
  const auto n = 1 + int(std::sqrt(nk_index));
  const auto kappa_index = nk_index - states_below_n(n);
  return kindex_to_twoj(kappa_index);
}

//! Returns l, given nk_index
//! @details See \ref Angular::nk_to_index().
inline int nkindex_to_l(std::uint64_t nk_index) {
  const auto n = 1 + int(std::sqrt(double(nk_index) + 0.01));
  const auto kappa_index = nk_index - states_below_n(n);
  return l_k(kindex_to_kappa(kappa_index));
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
/*!
  @brief Wigner 3j symbol. Inputs are 2*j (or 2*l) as integers.
  @details
  Computes the Wigner 3j symbol:
  \f[
    \begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}
  \f]
  Returns zero if the triangle rule or \f$ m_1+m_2+m_3=0 \f$ is not satisfied.
  Wraps \c gsl_sf_coupling_3j.

  @param two_j1  2*j1
  @param two_j2  2*j2
  @param two_j3  2*j3
  @param two_m1  2*m1
  @param two_m2  2*m2
  @param two_m3  2*m3
  @return Value of the 3j symbol; 0 if selection rules are not met.
*/
inline double threej_2(int two_j1, int two_j2, int two_j3, int two_m1,
                       int two_m2, int two_m3) {
  if (triangle(two_j1, two_j2, two_j3) * sumsToZero(two_m1, two_m2, two_m3) ==
      0)
    return 0;
  return gsl_sf_coupling_3j(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

//==============================================================================
/*!
  @brief Special 3j symbol: (j1 j2 k; -1/2 1/2 0). Inputs are 2*j as integers.
  @details
  Evaluates the commonly occurring 3j symbol:
  \f[
    \begin{pmatrix} j_1 & j_2 & k \\ -\tfrac{1}{2} & \tfrac{1}{2} & 0 \end{pmatrix}
  \f]
  Handles the \f$ k=0 \f$ case analytically; otherwise wraps \c gsl_sf_coupling_3j.
  Returns zero if the triangle rule is not satisfied.

  @param two_j1  2*j1
  @param two_j2  2*j2
  @param two_k   2*k
  @return Value of the 3j symbol; 0 if selection rules are not met.
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
/*!
  @brief Clebsch-Gordon coefficient <j1 m1, j2 m2 | J M>. Takes 2*j as integers.
  @details
  Computes the Clebsch-Gordon coefficient:
  \f[
    \langle j_1 m_1,\, j_2 m_2 \,|\, J M \rangle
    = (-1)^{j_1-j_2+M}\sqrt{2J+1}
      \begin{pmatrix} j_1 & j_2 & J \\ m_1 & m_2 & -M \end{pmatrix}
  \f]
  Returns zero if the triangle rule or \f$ m_1+m_2=M \f$ is not satisfied.

  @param two_j1  2*j1
  @param two_m1  2*m1
  @param two_j2  2*j2
  @param two_m2  2*m2
  @param two_J   2*J
  @param two_M   2*M
  @return Value of the CG coefficient; 0 if selection rules are not met.
*/
inline double cg_2(int two_j1, int two_m1, int two_j2, int two_m2, int two_J,
                   int two_M) {
  if (triangle(two_j1, two_j2, two_J) * sumsToZero(two_m1, two_m2, -two_M) == 0)
    return 0;
  int sign = -1;
  if ((two_j1 - two_j2 + two_M) % 4 == 0)
    sign = 1; // mod 4 (instead 2), since x2
  return sign * std::sqrt(two_J + 1.) *
         gsl_sf_coupling_3j(two_j1, two_j2, two_J, two_m1, two_m2, -two_M);
}

//==============================================================================
/*!
  @brief Returns true if the 6j symbol is zero by triangle/parity rules. Inputs are 2*j.
  @details
  Checks all four triangle triads and integer-sum (parity) conditions for:
  \f[
    \sixj{a}{b}{c}{d}{e}{f}
  \f]
  The triads (abc), (cde), (bdf), (efa) must each satisfy the triangle rule
  and sum to an even integer (inputs are 2*j).

  @note Some symbols that pass these checks are still numerically zero.
        GSL may return a very small non-zero value in such cases.
  @warning Only valid for inputs that are 2*(integer or half-integer).
*/
inline bool sixj_zeroQ(int a, int b, int c, int d, int e, int f) {
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

/*!
  @brief Checks 6j triangle conditions with optional arguments (wildcards). Inputs are 2*j.
  @details
  Checks only the triads formed by the provided (non-null) arguments.
  Missing (empty) optionals act as wildcards.
  For example, \c sixjTriads(a,b,c,{},{},{}) checks whether any 6j symbol
  \f$ \{ a\, b\, c \;|\; *\, *\, * \} \f$ could be non-zero.

  The four triads of \f$ \{a\,b\,c\;|\;d\,e\,f\} \f$ are: (abc), (cde), (aef), (bdf).

  @param a  2*j1 (optional)
  @param b  2*j2 (optional)
  @param c  2*j3 (optional)
  @param d  2*j4 (optional)
  @param e  2*j5 (optional)
  @param f  2*j6 (optional)
  @return false if any complete triad fails the triangle rule; true otherwise.
*/
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
/*!
  @brief Wigner 6j symbol {j1 j2 j3 | j4 j5 j6}. Inputs are 2*j as integers.
  @details
  Computes the Wigner 6j symbol:
  \f[
    \sixj{j_1}{j_2}{j_3}{j_4}{j_5}{j_6}
  \f]
  Returns zero if any triangle or parity condition is violated (via \ref sixj_zeroQ).
  Wraps \c gsl_sf_coupling_6j.

  @param two_j1  2*j1
  @param two_j2  2*j2
  @param two_j3  2*j3
  @param two_j4  2*j4
  @param two_j5  2*j5
  @param two_j6  2*j6
  @return Value of the 6j symbol; 0 if selection rules are not met.
*/
inline double sixj_2(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5,
                     int two_j6) {
  if (sixj_zeroQ(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6))
    return 0.0;
  return gsl_sf_coupling_6j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6);
}

//------------------------------------------------------------------------------
/*!
  @brief Wigner 9j symbol. Inputs are 2*j as integers.
  @details
  Computes the Wigner 9j symbol:
  \f[
    \begin{Bmatrix} j_1 & j_2 & j_3 \\ j_4 & j_5 & j_6 \\ j_7 & j_8 & j_9 \end{Bmatrix}
  \f]
  Wraps \c gsl_sf_coupling_9j.

  @param two_j1  2*j1
  @param two_j2  2*j2
  @param two_j3  2*j3
  @param two_j4  2*j4
  @param two_j5  2*j5
  @param two_j6  2*j6
  @param two_j7  2*j7
  @param two_j8  2*j8
  @param two_j9  2*j9
  @return Value of the 9j symbol.
*/
inline double ninej_2(int two_j1, int two_j2, int two_j3, int two_j4,
                      int two_j5, int two_j6, int two_j7, int two_j8,
                      int two_j9) {
  return gsl_sf_coupling_9j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6,
                            two_j7, two_j8, two_j9);
}

//==============================================================================
/*!
  @brief Reduced relativistic angular ME: <ka||C^k||kb>. Takes kappa values.
  @details
  Computes the reduced matrix element of the spherical tensor \f$ C^k \f$
  between relativistic states labelled by Dirac quantum numbers:
  \f[
    \langle \kappa_a \| C^k \| \kappa_b \rangle
    = (-1)^{j_a+\tfrac{1}{2}} \sqrt{[j_a][j_b]}
      \begin{pmatrix} j_a & j_b & k \\ -\tfrac{1}{2} & \tfrac{1}{2} & 0 \end{pmatrix}
      \pi(l_a+l_b+k),
  \f]
  where \f$ [j] \equiv 2j+1 \f$ and \f$ \pi(l_a+l_b+k) \f$ is 1 if \f$ l_a+l_b+k \f$
  is even, 0 otherwise.

  @param k   Multipole rank.
  @param ka  Dirac quantum number \f$ \kappa_a \f$ (bra).
  @param kb  Dirac quantum number \f$ \kappa_b \f$ (ket).
  @return \f$ \langle \kappa_a \| C^k \| \kappa_b \rangle \f$; 0 if selection rules are not met.
*/
inline double Ck_kk(int k, int ka, int kb) {
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

/*!
  @brief Full (non-reduced) C^k_q matrix element. Takes kappa and 2*m values.
  @details
  Computes the full matrix element including magnetic quantum numbers via
  the Wigner-Eckart theorem:
  \f[
    \langle \kappa_a m_a | C^k_q | \kappa_b m_b \rangle
    = (-1)^{j_a-m_a}
      \begin{pmatrix} j_a & k & j_b \\ -m_a & q & m_b \end{pmatrix}
      \langle \kappa_a \| C^k \| \kappa_b \rangle
  \f]

  @param k      Multipole rank.
  @param ka     Dirac quantum number \f$ \kappa_a \f$ (bra).
  @param kb     Dirac quantum number \f$ \kappa_b \f$ (ket).
  @param twoma  2*m_a
  @param twomb  2*m_b
  @param twoq   2*q
  @return \f$ \langle \kappa_a m_a | C^k_q | \kappa_b m_b \rangle \f$
*/
inline double Ck_kk_mmq(int k, int ka, int kb, int twoma, int twomb, int twoq) {
  const auto tja = twoj_k(ka);
  const auto tjb = twoj_k(kb);
  return neg1pow_2(tja - twoma) *
         threej_2(tja, 2 * k, tjb, -twoma, twoq, twomb) * Ck_kk(k, ka, kb);
}

/*!
  @brief Selection rule check for C^k: returns true if <ka||C^k||kb> may be non-zero.
  @details
  Checks parity and triangle conditions only; does not compute the value.
  Equivalent to \c Ck_kk(k,ka,kb)!=0 but avoids the full computation.

  @param k   Multipole rank.
  @param ka  Dirac quantum number \f$ \kappa_a \f$.
  @param kb  Dirac quantum number \f$ \kappa_b \f$.
  @return false if \f$ \langle \kappa_a \| C^k \| \kappa_b \rangle = 0 \f$ by selection rules; true otherwise.
*/
inline bool Ck_kk_SR(int k, int ka, int kb) {
  if (parity(l_k(ka), l_k(kb), k) == 0)
    return false;
  if (triangle(twoj_k(ka), twoj_k(kb), 2 * k) == 0)
    return false;
  return true;
}

/*!
  @brief Tilde variant of the C^k reduced ME: (-1)^{ja+1/2} * <ka||C^k||kb>
  @details
  Computes \f$ \widetilde{C}^k_{ab} \equiv (-1)^{j_a+1/2} \langle \kappa_a \| C^k \| \kappa_b \rangle \f$.
  This removes the \f$ (-1)^{j_a+1/2} \f$ phase from \ref Ck_kk, leaving
  the purely geometric 3j piece times \f$ \sqrt{[j_a][j_b]} \f$.

  @param k   Multipole rank.
  @param ka  Dirac quantum number \f$ \kappa_a \f$.
  @param kb  Dirac quantum number \f$ \kappa_b \f$.
  @return \f$ \widetilde{C}^k_{ab} \f$
*/
inline double tildeCk_kk(int k, int ka, int kb) {
  auto m1tjph = evenQ_2(twoj_k(ka) + 1) ? 1 : -1;
  return m1tjph * Ck_kk(k, ka, kb);
}

/*!
  @brief Min and max non-zero multipole rank k for C^k, given kappa_a and kappa_b.
  @details
  Returns \f$ (k_\text{min}, k_\text{max}) \f$ such that
  \f$ \langle \kappa_a \| C^k \| \kappa_b \rangle \neq 0 \f$ only for
  \f$ k_\text{min} \leq k \leq k_\text{max} \f$.
  Steps are in increments of 2, enforced by parity.

  Derived from:
  \f[
    k_\text{min} = |j_a - j_b|, \quad k_\text{max} = j_a + j_b,
  \f]
  then adjusted so that \f$ l_a + l_b + k \f$ is even.

  @param ka  Dirac quantum number \f$ \kappa_a \f$.
  @param kb  Dirac quantum number \f$ \kappa_b \f$.
  @return \f$ \{k_\text{min},\, k_\text{max}\} \f$
*/
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
/*!
  @brief Reduced spin-1/2 angular ME: <ka||S||kb>. Takes kappa values.
  @details
  Computes the reduced matrix element of the spin operator \f$ \mathbf{S} \f$
  (for spin-1/2) between relativistic states:
  \f[
    \langle \kappa_a \| S \| \kappa_b \rangle
    = \delta_{l_a l_b} \, (-1)^{j_a+l_a+3/2} \sqrt{[j_a][j_b]\tfrac{3}{2}}
      \sixj{j_a}{1}{j_b}{\tfrac{1}{2}}{l_a}{\tfrac{1}{2}}
  \f]
  The 6j symbol is evaluated analytically:
  \f[
    \sixj{j_a}{1}{j_b}{\tfrac{1}{2}}{l_a}{\tfrac{1}{2}}
    = \frac{(-1)^{\kappa_a+\kappa_b}}{2}
      \sqrt{\left|\frac{(\kappa_a-1)^2 - \kappa_b^2}{3\kappa_a(1+2\kappa_a)}\right|}
  \f]
  This special-case formula is ~20% faster than the general 6j routine.
  Returns 0 if \f$ l_a \neq l_b \f$ or the triangle rule is violated.

  @param ka  Dirac quantum number \f$ \kappa_a \f$.
  @param kb  Dirac quantum number \f$ \kappa_b \f$.
  @return \f$ \langle \kappa_a \| S \| \kappa_b \rangle \f$; 0 if selection rules are not met.
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
