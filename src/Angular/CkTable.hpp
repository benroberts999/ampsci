#pragma once
#include "Angular/Wigner369j.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

namespace Angular {

//******************************************************************************
// "Helper" functions
constexpr int twoj(int jindex) { return 2 * jindex + 1; }
constexpr int jindex(int twoj) { return (twoj - 1) / 2; }
constexpr int jindex_kappa(int ka) { return (ka > 0) ? ka - 1 : -ka - 1; }

//******************************************************************************
/*!
@brief Lookup table for C^k and 3j symbols (special m=1/2, q=0 case)
@details
  - Builds 3j symbol lookup table for given maximum k and maximum j (2j)
  - 3j symbols, special case: (ja jb k \\ -1/2, 1/2, 0)
  - Ck(ab)      := \f$C^k_{ab} = <ka||C^k||kb>\f$ [symmetric up to +/- sign]
  - TildeCk_ab := \f$(-1)^{ja+1/2} C^k_{ab}\f$ [symmetric]
  - Slightly faster than calculating on-the-fly, but adds some overhead
\par Construction
  - Needs maximum two*j values. Will build look-up tables for all possible
symbols.
\par Usage
  - Note: Functions take k and kappa_a, kappa_b as input!
*/
class CkTable {

public:
  CkTable(const int in_max_twoj = 0) { fill(in_max_twoj); }

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

} // namespace Angular
