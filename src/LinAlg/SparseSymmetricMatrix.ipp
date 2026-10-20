#pragma once
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <limits>

namespace LinAlg {

//==============================================================================

template <typename T>
std::size_t SparseSymmetricMatrix<T>::num_stored() const {
  std::size_t count = 0;
  for (const auto &row_columns : m_columns) {
    count += row_columns.size();
  }
  return count;
}

//==============================================================================

template <typename T>
void SparseSymmetricMatrix<T>::set_row(
  std::size_t row, const std::vector<std::uint32_t> &columns,
  const std::vector<T> &values) {
  assert(row < rows());
  assert(columns.size() == values.size());
  assert(std::is_sorted(columns.begin(), columns.end()));
  assert(columns.empty() || columns.back() <= row);
  // Copies (not moves): the copy has exactly the size needed, without the
  // spare capacity the input may have from being built by push_back
  m_columns[row] = std::vector<std::uint32_t>(columns.begin(), columns.end());
  m_values[row] = std::vector<T>(values.begin(), values.end());
}

//==============================================================================

template <typename T>
Vector<T> SparseSymmetricMatrix<T>::diagonal() const {
  Vector<T> diag(rows());
  for (std::size_t i = 0; i < rows(); ++i) {
    // Columns are ascending and <= i, so the diagonal is the last element
    if (!m_columns[i].empty() && m_columns[i].back() == i) {
      diag[i] = m_values[i].back();
    }
  }
  return diag;
}

//==============================================================================

template <typename T>
Matrix<T> SparseSymmetricMatrix<T>::dense_submatrix(
  const std::vector<std::size_t> &indices) const {
  // Position of each full-matrix index in the sub-matrix (or none)
  constexpr auto not_included = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> position(rows(), not_included);
  for (std::size_t a = 0; a < indices.size(); ++a) {
    assert(indices[a] < rows());
    assert(position[indices[a]] == not_included);
    position[indices[a]] = a;
  }

  // Both indices of each stored element are included, so every element of
  // the sub-matrix is stored in one of its rows
  Matrix<T> sub(indices.size(), indices.size());
  for (std::size_t a = 0; a < indices.size(); ++a) {
    const auto i = indices[a];
    for (std::size_t k = 0; k < m_columns[i].size(); ++k) {
      const auto b = position[m_columns[i][k]];
      if (b == not_included)
        continue;
      sub(a, b) = m_values[i][k];
      sub(b, a) = m_values[i][k];
    }
  }
  return sub;
}

//==============================================================================

template <typename T>
Matrix<T> SparseSymmetricMatrix<T>::dense() const {
  std::vector<std::size_t> all_indices(rows());
  for (std::size_t i = 0; i < rows(); ++i) {
    all_indices[i] = i;
  }
  return dense_submatrix(all_indices);
}

//==============================================================================

template <typename U>
Matrix<U> operator*(const Matrix<U> &x, const SparseSymmetricMatrix<U> &a) {
  assert(x.cols() == a.rows());
  const auto dimension = a.rows();
  const auto num_vectors = x.rows();

  // Vectors as columns: element j of every vector is then contiguous
  const auto x_columns = x.transpose();

  // A stored element (i, j < i) also stands for (j, i), which adds to row j
  // of the product: rows written by any thread. So each thread accumulates
  // its own product, and these are summed in thread order. The static
  // schedule fixes which thread does which row (round-robin chunks balance
  // the triangular rows), so the result does not change between runs.
  const auto num_threads = std::size_t(omp_get_max_threads());
  std::vector<Matrix<U>> thread_products(num_threads,
                                         Matrix<U>(dimension, num_vectors));
#pragma omp parallel for schedule(static, 64)
  for (std::size_t i = 0; i < dimension; ++i) {
    auto &product = thread_products[std::size_t(omp_get_thread_num())];
    const auto &columns = a.m_columns[i];
    const auto &values = a.m_values[i];
    const U *x_i = x_columns[i];
    U *y_i = product[i];
    for (std::size_t k = 0; k < columns.size(); ++k) {
      const std::size_t j = columns[k];
      const U a_ij = values[k];
      const U *x_j = x_columns[j];
      for (std::size_t r = 0; r < num_vectors; ++r) {
        y_i[r] += a_ij * x_j[r];
      }
      if (j != i) {
        U *y_j = product[j];
        for (std::size_t r = 0; r < num_vectors; ++r) {
          y_j[r] += a_ij * x_i[r];
        }
      }
    }
  }

  auto y_columns = std::move(thread_products.front());
  for (std::size_t t = 1; t < num_threads; ++t) {
    y_columns += thread_products[t];
  }
  return y_columns.transpose();
}

//==============================================================================

template <typename U>
Vector<U> operator*(const SparseSymmetricMatrix<U> &a, const Vector<U> &b) {
  assert(b.size() == a.cols());
  // b as a single-row matrix
  const Matrix<U> b_row(1, b.size(), std::vector<U>(b.cbegin(), b.cend()));
  const auto product = b_row * a;
  return Vector<U>(std::vector<U>(product.cbegin(), product.cend()));
}

} // namespace LinAlg
