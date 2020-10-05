#pragma once
#include "Angular/Angular_369j.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

namespace Angular {

//******************************************************************************
// "Helper" functions
constexpr int twoj(int jindex) { return 2 * jindex + 1; }
constexpr int jindex(int twoj) { return (twoj - 1) / 2; }
constexpr int jindex_kappa(int ka) { return (ka > 0) ? ka - 1 : -ka - 1; }

template <typename T> inline T max4(T a, T b, T c, T d) {
  return std::max(std::max(a, b), std::max(c, d));
  // Actually faster than c++17s: std::max({a,b,c,d})? Probs not.
}
template <typename T> inline T min4(T a, T b, T c, T d) {
  return std::min(std::min(a, b), std::min(c, d));
}

//******************************************************************************
/*!
@brief Lookup table for C^k and 3j symbols (special m=1/2, q=0 case)
@details
  - Builds 3j symbol lookup table for given maximum k and maximum j (2j)
  - 3j symbols, special case: (ja jb k \\ -1/2, 1/2, 0)
  - Ck_ab      := \f$C^k_{ab} = <ka||C^k||kb>\f$ [symmetric up to +/- sign]
  - TildeCk_ab := \f$(-1)^{ja+1/2} C^k_{ab}\f$ [symmetric]
  - Slightly faster than calculating on-the-fly, but adds some overhead
\par Construction
  - Needs maximum two*j values. Will build look-up tables for all possible
symbols.
\par Usage
  - Note: Functions take k and kappa_a, kappa_b as input!
*/
class Ck_ab {

public:
  Ck_ab(const int in_max_twoj = 0) { fill(in_max_twoj); }

public:
  //! @brief Extends existing look-up table to new twoj.
  //! @details nb: called on construction automatically, you only need to call
  //! this if you need to extend the table after original construction (rare)
  void fill(const int in_max_twoj);

  //! @brief 'mutable' getters will calculate values if they don't exist yet
  //! @details Safer, but NOT thread safe
  double get_Ckab_mutable(int k, int ka, int kb);
  double get_tildeCkab_mutable(int k, int ka, int kb);
  double get_3jkab_mutable(int k, int ka, int kb);

  //! @brief Use const versions if sure value already exists
  //! @details Note: undefined behaviour to call if the value doesn't exist yet
  double get_Ckab(int k, int ka, int kb) const;
  double get_tildeCkab(int k, int ka, int kb) const;
  double get_3jkab(int k, int ka, int kb) const;

  double operator()(int k, int ka, int kb) const { return get_Ckab(k, ka, kb); }

  //! Lambda^k_ij := 3js((ji,jj,k),(-1/2,1/2,0))^2 * parity(li+lj+k)
  double get_Lambdakab(int k, int ka, int kb) const;

  int max_tj() const { return twoj(m_max_jindex_sofar); }
  int max_k() const { return m_max_k_sofar; }

private:
  std::vector<std::vector<std::vector<double>>> m_3j_k_a_b = {};
  std::vector<std::vector<double>> m_Rjab_a_b = {}; // Sqrt([ja][jb])
  int m_max_jindex_sofar = -1;
  int m_max_k_sofar = -1;
};

//******************************************************************************

//******************************************************************************
class SixJTable_constk
//  { ja, jb, k}
//  { jc, jd, l}
// k is fixed
{
public:
  SixJTable_constk(int in_k, int tj_max = -1) : m_k(in_k) { fill(tj_max); }

private:
  const int m_k;
  int max_ji_sofar = -1;
  std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>
      m_k_a_bcdl = {};

public:
  double get_6j(int tja, int tjb, int tjc, int tjd, int l) const;
  double get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const;
  //"safe" version that will allocate new 6j if not yet (NOT thresh-safe!)
  double get_6j_mutable(int tja, int tjb, int tjc, int tjd, int l);

  int get_k() const { return m_k; }
  int get_maxtj() const { return twoj(max_ji_sofar); }

public:
  void fill(const int tj_max);
};

//******************************************************************************
/*!
@brief Lookup table for 6j symbols
@details
  - Lookup table for 6j symbols: { ja, jb, k \\ jc, jd, l}
  - j's half-integer (called using integer 2j). Integer k and l
  - Much faster than calculating on the fly.
  - Stores 6j symbols up to + including given max_K and max_2j (stores all
allowed l)
\par Construction
  - Needs max k and two*j values. Will build look-up tables for all possible
symbols.
\par Usage
  - Note: Functions take k and {two*j} as input! Also: note input order (strange
order..)
*/
class SixJ
//  { ja, jb, k}
//  { jc, jd, l}
// any k
{
public:
  SixJ(int in_max_twoj = 1) { fill(in_max_twoj); }
  //! @brief Extends existing look-up table to new maximum K and twoj.
  void fill(int in_max_twoj);

  //! @brief Thread-safe, but will seg-fault if 6j doesn't exist XXX CHANGE
  //! ORDER ?? XXX
  double get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const;

  double operator()(int tja, int tjb, int tjc, int tjd, int k, int l) const {
    return get_6j(tja, tjb, tjc, tjd, k, l);
  }

  //! @brief Will calculate + store 6j if it doesn't exist, but not thread-safe
  double get_6j_mutable(int tja, int tjb, int tjc, int tjd, int k, int l);

  int max_tj() const { return max_tj_sofar; }
  int max_k() const { return max_k_sofar; }

private:
  int max_k_sofar = -1;
  int max_tj_sofar = -1;
  std::vector<SixJTable_constk> m_sixj_k = {};
};

} // namespace Angular
