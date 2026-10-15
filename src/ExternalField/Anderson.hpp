#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <cstddef>
#include <vector>

namespace ExternalField {

/*!
  @brief Anderson mixing coefficients for a fixed-point
  iteration x -> G(x).
  @details
  Given the residuals \f$ r_k = G(x_k) - x_k \f$ of the stored iterates
  (oldest first) through their Gram matrix \f$ B_{kl} = \braket{r_k}{r_l} \f$,
  returns the \f$ c_k \f$ minimising
  \f[ \Big\|\sum_k c_k r_k\Big\|^2 \quad\text{subject to}\quad \sum_k c_k = 1, \f]
  via the bordered system [B 1; 1^T 0][c; lambda] = [0; 1]. The next iterate
  is \f$ x = \sum_k c_k G(x_k) \f$. For a linear map this is equivalent to
  (truncated) GMRES: it converges wherever the linear problem is
  non-singular, where damped iteration may not.

  A saturated history makes B ill-conditioned: the oldest entries are then
  dropped, and the returned vector holds one coefficient per KEPT entry, the
  last c.size() entries of the history (the caller drops the same entries
  from its own history). Empty if no entry is usable (take a plain step).
*/
std::vector<double>
anderson_coefficients(const std::vector<std::vector<double>> &gram);

//! As above, from the residual spinors directly (Gram matrix of their inner
//! products).
std::vector<double>
anderson_coefficients(const std::vector<DiracSpinor> &residuals);

/*!
  @brief Anderson mixing history of a fixed-point iteration x -> G(x) on
  spinors, as used by the mixed-state solvers.
  @details
  Holds the map outputs \f$ g_k = G(x_k) \f$, the residuals
  \f$ r_k = g_k - x_k \f$ and their Gram matrix, and a scalar carried with
  each entry (an amplitude linear in the state, mixed with the same
  coefficients; zero if unused). The Gram matrix is kept incrementally: a
  new entry costs one inner product per kept entry, not a rebuild of the
  whole matrix. At most @p depth entries are kept: the oldest is dropped
  when the history is full, and (see @ref anderson_coefficients) while the
  Gram matrix is ill-conditioned.
*/
class AndersonHistory {
public:
  explicit AndersonHistory(std::size_t depth) : m_depth(depth) {
    assert(depth > 0 && "AndersonHistory needs a depth of at least one");
  }

  //! Adds the entry (g, r), with |r|^2 = r2 (which the caller has from its
  //! convergence test) and the carried scalar k; drops the oldest entry if
  //! the history is full
  void push(DiracSpinor g, DiracSpinor r, double r2, double k = 0.0);

  /*!
    @brief The Anderson extrapolation over the kept entries,
    x = sum_k c_k g_k and (if @p k is non-null) k = sum_k c_k k_k, with the
    coefficients of @ref anderson_coefficients; entries that had to be
    dropped for conditioning are dropped from the history too.
    @details If no entry is usable (the coefficients could not be solved
    for), x is the plain step: the latest map output, and the history is
    cleared.
  */
  void mix(DiracSpinor *x, double *k = nullptr);

  //! Number of kept entries
  std::size_t size() const { return m_g.size(); }

private:
  std::size_t m_depth;
  std::vector<DiracSpinor> m_g{};
  std::vector<DiracSpinor> m_r{};
  std::vector<double> m_k{};
  // Gram matrix of the residuals, <r_i|r_j>, over the kept entries
  std::vector<std::vector<double>> m_gram{};

  // Drops the n oldest entries (and their Gram rows and columns)
  void drop_oldest(std::size_t n);
};

/*!
  @brief Anderson mixing history of a fixed-point iteration x -> G(x) on a
  state held as one real vector, as used by the TDHF solvers (the whole set
  of core corrections, flattened).
  @details
  As @ref AndersonHistory: holds the map outputs \f$ g_k = G(x_k) \f$, the
  residuals \f$ r_k = g_k - x_k \f$ and their Gram matrix (Euclidean inner
  product), kept incrementally. At most @p depth entries are kept: the
  oldest is dropped when the history is full, and (see
  @ref anderson_coefficients) while the Gram matrix is ill-conditioned.
*/
class AndersonVectorHistory {
public:
  explicit AndersonVectorHistory(std::size_t depth) : m_depth(depth) {
    assert(depth > 0 && "AndersonVectorHistory needs a depth of at least one");
  }

  //! Adds the entry for iterate x with map output g = G(x) (residual
  //! g - x); drops the oldest entry if the history is full
  void push(const std::vector<double> &x, std::vector<double> g);

  /*!
    @brief The Anderson extrapolation over the kept entries,
    x = sum_k c_k g_k, with the coefficients of @ref anderson_coefficients;
    entries that had to be dropped for conditioning are dropped from the
    history too.
    @details If no entry is usable (the coefficients could not be solved
    for), returns the plain step: the latest map output, and the history is
    cleared.
  */
  std::vector<double> mix();

  //! Number of kept entries
  std::size_t size() const { return m_g.size(); }

private:
  std::size_t m_depth;
  std::vector<std::vector<double>> m_g{};
  std::vector<std::vector<double>> m_r{};
  // Gram matrix of the residuals, <r_i|r_j>, over the kept entries
  std::vector<std::vector<double>> m_gram{};

  // Drops the n oldest entries (and their Gram rows and columns)
  void drop_oldest(std::size_t n);
};

} // namespace ExternalField
