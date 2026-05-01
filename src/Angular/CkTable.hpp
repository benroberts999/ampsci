#pragma once
#include "Angular/Wigner369j.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

namespace Angular {

//==============================================================================
// Helper index-conversion functions
//! Converts jindex to 2j; jindex = 0,1,2,... maps to 2j = 1,3,5,...
constexpr int jindex_to_twoj(int jindex) { return 2 * jindex + 1; }
//! Converts 2j to jindex; 2j = 1,3,5,... maps to jindex = 0,1,2,...
constexpr int twoj_to_jindex(int twoj) { return (twoj - 1) / 2; }
//! Converts kappa to jindex; e.g. {-1, 1, -2, 2} -> {0, 0, 1, 1}
constexpr int kappa_to_jindex(int ka) { return (ka > 0) ? ka - 1 : -ka - 1; }

//==============================================================================
/*!
  @brief Lookup table for Ck_ab reduced matrix elements and the
         special 3j symbol (j_a, j_b, k; -1/2, 1/2, 0).
  @details
  Pre-computes and caches the angular symbols needed for radial matrix-element
  evaluation. Three quantities are stored for each triple \f$(k, \kappa_a,
  \kappa_b)\f$:

  - The reduced matrix element of the normalised spherical harmonic:
    \f[
      C^k_{ab} \equiv \redmatel{\kappa_a}{C^k}{\kappa_b}
      = (-1)^{j_a+\frac{1}{2}}\sqrt{[j_a][j_b]}
        \threej{j_a}{j_b}{k}{-\tfrac{1}{2}}{\tfrac{1}{2}}{0}
        \,\pi(l_a+l_b+k),
    \f]
    where \f$[x]\equiv 2x+1\f$ and \f$\pi\f$ encodes the parity selection
    rule. \f$C^k_{ab}\f$ is antisymmetric under interchange up to a phase.

  - The symmetrised (tilde) variant:
    \f[
      \widetilde C^k_{ab} = (-1)^{j_a+\frac{1}{2}} C^k_{ab},
    \f]
    which is fully symmetric: \f$\widetilde C^k_{ab} = \widetilde C^k_{ba}\f$.

  - The underlying special 3j symbol:
    \f[
      \threej{j_a}{j_b}{k}{-\tfrac{1}{2}}{\tfrac{1}{2}}{0}.
    \f]

  All table lookups take \f$k\f$ and the Dirac quantum numbers
  \f$\kappa_a, \kappa_b\f$ as input.

  @note
  - Table is filled up to a maximum \f$2j\f$ at construction; call fill()
    to extend it afterwards.
  - All symbols with \f$2j \le \f$ max_tj() and \f$k \le \f$ max_k() are
    guaranteed to be present; there are no gaps.
  - Non-mutable accessors are thread-safe. The `_mutable` variants extend the
    table on demand but are **not** thread-safe.

  @warning Calling a non-mutable accessor with out-of-range indices is
           undefined behaviour; check against max_tj() and max_k() first.
*/
class CkTable {

public:
  /*!
    @brief Constructs the table, pre-filling all symbols up to @p in_max_twoj.
    @details Calls fill() internally; passing zero (the default) creates an
             empty table that can be extended later with fill() or the mutable
             accessors.
    @param in_max_twoj  Maximum value of \f$2j\f$ to pre-compute.
  */
  CkTable(const int in_max_twoj = 0) { fill(in_max_twoj); }

public:
  /*!
    @brief Extends the lookup table to cover all symbols up to @p in_max_twoj.
    @details Called automatically on construction. Only needed explicitly when
             the required \f$2j\f$ grows beyond the original maximum.
    @param in_max_twoj  New maximum value of \f$2j\f$; no-op if already met.
  */
  void fill(const int in_max_twoj);

  //! Returns Ck_ab; extends table if needed. Not thread-safe.
  double get_Ckab_mutable(int k, int ka, int kb);
  //! Returns tilde-Ck_ab; extends table if needed. Not thread-safe.
  double get_tildeCkab_mutable(int k, int ka, int kb);
  //! Returns 3j(j_a,j_b,k; -1/2,1/2,0); extends table if needed. Not thread-safe.
  double get_3jkab_mutable(int k, int ka, int kb);

  //! Returns Ck_ab. Undefined behaviour if indices exceed max_tj() or max_k().
  double get_Ckab(int k, int ka, int kb) const;
  //! Returns tilde-Ck_ab. Undefined behaviour if indices exceed max_tj() or max_k().
  double get_tildeCkab(int k, int ka, int kb) const;
  //! Returns 3j(j_a,j_b,k; -1/2,1/2,0). Undefined behaviour if indices exceed max_tj() or max_k().
  double get_3jkab(int k, int ka, int kb) const;

  //! Equivalent to get_Ckab(@p k, @p ka, @p kb).
  double operator()(int k, int ka, int kb) const { return get_Ckab(k, ka, kb); }

  /*!
    @brief Returns Lambda^k_ab := 3j(j_a,j_b,k; -1/2,1/2,0)^2 * parity(l_a+l_b+k).
    @details
    This is the angular factor that appears in the Coulomb interaction after
    angular reduction, and equals \f$\left(C^k_{ab}\right)^2 / ([j_a][j_b])\f$
    up to a phase.
  */
  double get_Lambdakab(int k, int ka, int kb) const;

  //! Maximum value of 2j currently stored in the table.
  int max_tj() const { return jindex_to_twoj(m_max_jindex_sofar); }
  //! Maximum value of k currently stored in the table.
  int max_k() const { return m_max_k_sofar; }

private:
  std::vector<std::vector<std::vector<double>>> m_3j_k_a_b = {};
  std::vector<std::vector<double>> m_Rjab_a_b = {}; // Sqrt([ja][jb])
  int m_max_jindex_sofar = -1;
  int m_max_k_sofar = -1;
};

} // namespace Angular
