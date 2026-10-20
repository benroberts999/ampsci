#pragma once
#include "Matrix.hpp"
#include "Vector.hpp"
#include <cstdint>
#include <vector>

namespace LinAlg {

//==============================================================================
/*!
  @brief Real symmetric sparse matrix, stored as the non-zero elements of its
  lower triangle.
  @details
  For large symmetric matrices whose elements are mostly zero, where the dense
  matrix may not fit in memory. Its main use is the matrix-vector product,
  e.g., for iterative eigensolvers such as the Davidson method.

  Storage:

  * Each row i holds two arrays: the column indices \f$ j \le i \f$
    (ascending) and the matching values. Only the lower triangle is stored;
    the upper triangle follows by symmetry.
  * Each stored element takes `sizeof(T) + 4` bytes: the value plus a 4-byte
    column index.
  * The product \f$ y = A x \f$ streams through each row's arrays once, in
    order: each stored \f$ A_{ij} \f$ adds \f$ A_{ij} x_j \f$ to \f$ y_i \f$
    and \f$ A_{ij} x_i \f$ to \f$ y_j \f$.
  * There is no single-element access. The layout suits sweeps over all
    elements (sequential, cache-friendly), not random lookup of individual
    elements.

  Rows are set independently with set_row(), so different rows may be filled
  in parallel.

  @tparam T Element type: `double` or `float`.
*/
template <typename T>
class SparseSymmetricMatrix {
  // For each row i: column indices (j <= i, ascending) and values of the
  // stored elements
  std::vector<std::vector<std::uint32_t>> m_columns{};
  std::vector<std::vector<T>> m_values{};

public:
  //! Default construct: empty (0 x 0) matrix
  SparseSymmetricMatrix() {}

  //! Construct a zero matrix of size `dimension x dimension`
  explicit SparseSymmetricMatrix(std::size_t dimension)
    : m_columns(dimension), m_values(dimension) {}

  //! Number of rows (equal to the number of columns)
  std::size_t rows() const { return m_columns.size(); }
  //! Number of columns (equal to the number of rows)
  std::size_t cols() const { return m_columns.size(); }

  //! Number of stored elements: the non-zero elements of the lower triangle,
  //! including the diagonal
  std::size_t num_stored() const;

  /*!
    @brief Sets row @p row of the lower triangle.
    @details
    Element \f$ A_{\rm row, columns[k]} = \f$ values[k] (and its symmetric
    partner). Any previous contents of the row are replaced. Different rows
    may be set in parallel; the same row must not be.

    The row is copied, not moved, from @p columns and @p values. Arrays built
    up element by element (push_back) may have up to twice the capacity they
    need; moving them in would keep that spare capacity, up to doubling the
    memory of the matrix. The copy allocates exactly the size needed. The
    cost of the copy is small next to computing the elements.

    @param row      Row index.
    @param columns  Column indices of the elements: ascending, and each
                    <= @p row.
    @param values   Element values; same length as @p columns.
  */
  void set_row(std::size_t row, const std::vector<std::uint32_t> &columns,
               const std::vector<T> &values);

  //! The diagonal elements, \f$ A_{ii} \f$
  [[nodiscard]] Vector<T> diagonal() const;

  /*!
    @brief Dense sub-matrix: the rows and columns given by @p indices.
    @details
    Returns the dense matrix \f$ S_{ab} = A_{i_a i_b} \f$, both triangles
    filled, with \f$ i_a \f$ = indices[a]. The indices must be distinct, but
    need not be ordered. Costs one pass over the stored elements of the rows
    in @p indices.
  */
  [[nodiscard]] Matrix<T>
  dense_submatrix(const std::vector<std::size_t> &indices) const;

  //! The full dense matrix (both triangles filled)
  [[nodiscard]] Matrix<T> dense() const;

  /*!
    @brief Product of each row of @p x with the matrix: \f$ Y = X A \f$.
    @details
    Since A is symmetric, row r of the result is A applied to row r of @p x:
    \f$ y_r = A x_r \f$. Applying A to several vectors at once reads the
    stored elements only once. @p x must have a.rows() columns. Parallelised
    over rows (OpenMP); the result does not depend on the scheduling, only
    (at the level of rounding) on the number of threads.
  */
  template <typename U>
  friend Matrix<U> operator*(const Matrix<U> &x,
                             const SparseSymmetricMatrix<U> &a);

  //! Matrix-vector product: returns \f$ A b \f$
  template <typename U>
  friend Vector<U> operator*(const SparseSymmetricMatrix<U> &a,
                             const Vector<U> &b);
};

} // namespace LinAlg

#include "SparseSymmetricMatrix.ipp"
