#pragma once
#include "Angular/Wigner369j.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

namespace Angular {

//==============================================================================
// "Helper" functions
//! Converts jindex to 2*j [helper function]
constexpr int twoj(int jindex) { return 2 * jindex + 1; }
//! Converts 2*j to jindex {1/2, 3/2, 5/2} -> {0, 1, 2} [helper function]
constexpr int jindex(int twoj) { return (twoj - 1) / 2; }
//! Converts kappa to jindex {-1, 1, -2} -> {0, 0, 1} [helper function]
constexpr int jindex_kappa(int ka) { return (ka > 0) ? ka - 1 : -ka - 1; }

//==============================================================================
/*!
@brief Lookup table for C^k and 3j symbols (special m=1/2, q=0 case)
@details
  - Builds 3j symbol lookup table for given maximum k and maximum j (2j)
  - 3j symbols, special case: (ja jb k \\ -1/2, 1/2, 0)
  - Ckab : \f$C^k_{ab} = \langle k_a||C^k||k_b \rangle\f$ 
    [symmetric up to +/- sign]
  - TildeCkab : \f$\widetilde C^k_{ab} = (-1)^{ja+1/2} C^k_{ab}\f$ [symmetric]
  - Slightly faster than calculating on-the-fly, but adds some overhead
\par Construction
  - Needs maximum two*j values. Will build look-up tables for all possible
symbols.
\par Mutable versions
  - Some functions come with '_mutable' versions
  - These will calculate (+ store) values if they don't exist yet
  - '_mutable' versions are NOT thread safe (non-mutable are)
  - The non-mutable versions: must not be called if symbol hasn't been 
    calculated. This is undefined behaviour
  - You can check which symbols exist by calling max_tj() and max_k()
  - All "lower" symbols are calculated [max_tj, max_k defined all]
\par Usage
  - Note: Functions take k and kappa_a, kappa_b as input!
*/
class CkTable {

public:
  //! Calculates and stored all Ck/3j symbols up to given maximum 2j
  CkTable(const int in_max_twoj = 0) { fill(in_max_twoj); }

public:
  //! Extends existing look-up table to new twoj.
  /*! 
  @details 
  nb: called on construction automatically, you only need to call this if you 
  need to extend the table after original construction (rare)
  */
  void fill(const int in_max_twoj);

  //! Ckab. mutable: will calculate if needed. Not thread safe
  double get_Ckab_mutable(int k, int ka, int kb);
  //! tildeCkab. mutable: will calculate if needed. Not thread safe
  double get_tildeCkab_mutable(int k, int ka, int kb);
  //! special 3j(k, ka, kb). mutable: will calculate if needed. Not thread safe
  double get_3jkab_mutable(int k, int ka, int kb);

  //! Ckab. Undefined if k, ka, or kb are out-of-bounds [check with max_tj()].
  double get_Ckab(int k, int ka, int kb) const;
  //! tildeCkab. Undefined if k, ka, or kb are out-of-bounds
  double get_tildeCkab(int k, int ka, int kb) const;
  //! special 3j(k, ka, kb). Undefined if k, ka, or kb are out-of-bounds
  double get_3jkab(int k, int ka, int kb) const;

  //! Operator overload: returns Ckab
  double operator()(int k, int ka, int kb) const { return get_Ckab(k, ka, kb); }

  //! Lambda^k_ij := 3js((ji,jj,k),(-1/2,1/2,0))^2 * parity(li+lj+k)
  double get_Lambdakab(int k, int ka, int kb) const;

  //! Maximum value for 2j currently stored in tables
  int max_tj() const { return twoj(m_max_jindex_sofar); }
  //! Maximum value for k currently stored in tables
  int max_k() const { return m_max_k_sofar; }

private:
  std::vector<std::vector<std::vector<double>>> m_3j_k_a_b = {};
  std::vector<std::vector<double>> m_Rjab_a_b = {}; // Sqrt([ja][jb])
  int m_max_jindex_sofar = -1;
  int m_max_k_sofar = -1;
};

} // namespace Angular
