#pragma once
#include "Template.hpp"
#include <array>
#include <cassert>
#include <numeric>
#include <vector>

/*!
  @brief General-purpose utility library.

  @details
  Provides lightweight, header-only utilities used throughout ampsci, including:
  - array and vector operations (@ref Array.hpp, @ref Vector.hpp),
  - numerical methods (@ref Maths.hpp, @ref Methods.hpp),
  - string helpers (@ref String.hpp),
  - strong typing (@ref StrongType.hpp),
  - template metaprogramming helpers (@ref Template.hpp),
  - random number utilities (@ref Random.hpp),
  - OpenMP helpers (@ref omp.hpp),
    - include this instead of `<omp.h>` to allow compilation with+without OpenMP
  - output/display widgets (@ref Widgets.hpp).
*/
namespace qip {

//==============================================================================
/*!
  @brief Iterator with a configurable stride.

  @details
  Not fully compliant with the standard forward iterator concept,
  but works with most standard library algorithms.
  Reverse iterators can be formed with a negative stride.
  See also @ref ConstStrideIterator.
*/
template <typename T>
class StrideIterator : public qip::Comparison<StrideIterator<T>> {
protected:
  T *m_ptr;
  long m_stride;

public:
  StrideIterator(T *ptr, long stride) : m_ptr(ptr), m_stride(stride) {}

  T &operator*() { return *m_ptr; }

  const T &operator*() const { return *m_ptr; }

  bool operator==(const StrideIterator &other) const {
    return m_ptr == other.m_ptr;
  }
  bool operator<(const StrideIterator &other) const {
    return m_ptr < other.m_ptr;
  }

  StrideIterator &operator++() {
    m_ptr += m_stride;
    return *this;
  }
  StrideIterator &operator--() {
    m_ptr -= m_stride;
    return *this;
  }
  StrideIterator operator++(int) {
    auto temp = *this;
    ++*this;
    return temp;
  }
  StrideIterator operator--(int) {
    auto temp = *this;
    --*this;
    return temp;
  }

  StrideIterator &operator+=(long n) {
    m_ptr += n * m_stride;
    return *this;
  }

  StrideIterator &operator-=(long n) {
    m_ptr -= n * m_stride;
    return *this;
  }

  StrideIterator operator+(long n) const {
    auto out_iter = *this;
    return out_iter += n;
  }

  StrideIterator operator-(long n) const {
    auto out_iter = *this;
    return out_iter -= n;
  }
};

//! Const version of @ref StrideIterator.
template <typename T>
class ConstStrideIterator : public qip::Comparison<ConstStrideIterator<T>> {
protected:
  T *m_ptr;
  long m_stride;

public:
  ConstStrideIterator(T *ptr, long stride) : m_ptr(ptr), m_stride(stride) {}

  const T &operator*() const { return *m_ptr; }

  bool operator==(const ConstStrideIterator &other) const {
    return m_ptr == other.m_ptr;
  }
  bool operator<(const ConstStrideIterator &other) const {
    return m_ptr < other.m_ptr;
  }

  ConstStrideIterator &operator++() {
    m_ptr += m_stride;
    return *this;
  }
  ConstStrideIterator &operator--() {
    m_ptr -= m_stride;
    return *this;
  }
  ConstStrideIterator operator++(int) {
    auto temp = *this;
    ++*this;
    return temp;
  }
  ConstStrideIterator operator--(int) {
    auto temp = *this;
    --*this;
    return temp;
  }

  ConstStrideIterator &operator+=(long n) {
    m_ptr += n * m_stride;
    return *this;
  }

  ConstStrideIterator &operator-=(long n) {
    m_ptr -= n * m_stride;
    return *this;
  }

  ConstStrideIterator operator+(long n) const {
    auto out_iter = *this;
    return out_iter += n;
  }

  ConstStrideIterator operator-(long n) const {
    auto out_iter = *this;
    return out_iter -= n;
  }
};

//==============================================================================
/*!
  @brief Non-owning view onto a 1D contiguous or strided array segment.

  @details
  - Size is the total number of elements.
  - Can have zero size, but accessing data is undefined in that case.
  - Interface mostly mirrors std::vector; no surprises.
*/
template <typename T = double>
class ArrayView {

private:
  std::size_t m_size; // number of elements
  std::size_t m_stride;
  T *m_data;

public:
  ArrayView(T *data, std::size_t size, std::size_t stride = 1)
    : m_size(size), m_stride(stride), m_data(data) {}

  std::size_t size() const { return m_size; }

  T &operator[](std::size_t i) { return m_data[i * m_stride]; }

  T operator[](std::size_t i) const { return m_data[i * m_stride]; }

  T &at(std::size_t i) {
    assert(i < m_size);
    return m_data[i * m_stride];
  }

  T at(std::size_t i) const {
    assert(i < m_size);
    return m_data[i * m_stride];
  }

  T &operator()(std::size_t i) { return at(i); }
  T operator()(std::size_t i) const { return at(i); }

  T &front() { return at(0); }
  T front() const { return at(0); }
  T &back() { return at(m_size - 1); }
  T back() const { return at(m_size - 1); }

  T *data() { return m_data; }

  //! Iterator to the beginning
  auto begin() { return StrideIterator(m_data, long(m_stride)); }
  auto cbegin() const { return ConstStrideIterator(m_data, long(m_stride)); }

  auto end() {
    return StrideIterator(m_data + long(m_size * m_stride), long(m_stride));
  }
  auto cend() const {
    return ConstStrideIterator(m_data + long(m_size * m_stride),
                               long(m_stride));
  }

  auto rbegin() {
    return StrideIterator(m_data + long(m_size * m_stride) - long(m_stride),
                          -long(m_stride));
  }
  auto crbegin() const {
    return ConstStrideIterator(
      m_data + long(m_size * m_stride) - long(m_stride), -long(m_stride));
  }

  auto rend() {
    return StrideIterator(m_data - long(m_stride), -long(m_stride));
  }
  auto crend() const {
    return ConstStrideIterator(m_data - long(m_stride), -long(m_stride));
  }

  //! Returns a copy of the view as a std::vector.
  std::vector<T> vector() {
    std::vector<T> out;
    for (std::size_t i = 0; i < m_size; ++i) {
      out.push_back(m_data[i * m_stride]);
    }
    return out;
  }
};

//==============================================================================

//! Variadic product - helper function
template <typename T, typename... Args>
T product(T first, Args... rest) {
  if constexpr (sizeof...(rest) == 0) {
    return first;
  } else {
    return first * product(rest...);
  }
}

template <std::size_t N>
void NDrange_impl(std::vector<std::array<std::size_t, N>> &result,
                  std::array<std::size_t, N> &current,
                  const std::array<std::size_t, N> &maxValues,
                  std::size_t index) {
  if (index == N) {
    result.push_back(current);
    return;
  }

  for (std::size_t i = 0; i < maxValues[index]; ++i) {
    current[index] = i;
    NDrange_impl<N>(result, current, maxValues, index + 1);
  }
}

/*!
  @brief Returns all index tuples for an N-dimensional range.

  @details
  Enables iteration over every index combination, e.g.:
  @code
  for (auto [i,j,k] : qip::NDrange(3,4,2)) { array.at(i,j,k); }
  @endcode
  @note Not memory efficient; do not use for large arrays.
*/
template <typename... Args>
auto NDrange(std::size_t first, Args... rest) {
  constexpr std::size_t N = sizeof...(rest) + 1;

  const std::array<std::size_t, N> maxValues = {
    first, static_cast<std::size_t>(rest)...};

  std::vector<std::array<std::size_t, N>> result;
  result.reserve(product(first, static_cast<std::size_t>(rest)...));

  std::array<std::size_t, N> current = {0};

  NDrange_impl<N>(result, current, maxValues, 0);
  return result;
}

//==============================================================================
/*!
  @brief N-dimensional array with arithmetic operators.

  @details
  Stores data in a flat std::vector; indices are computed from the shape.
  Element-wise arithmetic (+, -, *, /) is provided between arrays of identical
  shape, and between an array and its element type (see operator+=, etc.).
  Row and column views (@ref ArrayView) are available for 2D arrays.
*/
template <typename T = double>
class Array : public Arithmetic<Array<T>>, Arithmetic2<Array<T>, T> {

private:
  // List of sizes for each array dimension
  std::vector<std::size_t> m_sizes;
  // Number of array dimensions
  std::size_t m_Ndim;
  // Cumulative sizes (used to index into data)
  std::vector<std::size_t> m_cumulative_sizes;
  // Total number of elements
  std::size_t m_total_size;
  // Raw data
  std::vector<T> m_data;

public:
  /*!
    @brief Constructs an N-dimensional array with the given dimension sizes.

    @details e.g., Array(3,6,2) creates a 362 array.
  */
  template <typename... Args>
  Array(std::size_t first, Args... rest);

  //! Resizes the array. Note: invalidates all underlying data.
  template <typename... Args>
  void resize(std::size_t first, Args... rest);

  //! Returns the total number of elements.
  std::size_t size() const { return m_total_size; }

  //! Returns the size of a specific dimension.
  std::size_t size(std::size_t dim) const { return m_sizes.at(dim); }

  //! Returns the number of dimensions.
  std::size_t dimensions() const { return m_Ndim; }

  //! Returns the shape (sizes of all dimensions) of the array.
  const std::vector<std::size_t> &shape() const { return m_sizes; }

  //! Access element with bounds checking.
  template <typename... Args>
  T &at(std::size_t first, Args... rest);

  //! Const access to element with bounds checking.
  template <typename... Args>
  T at(std::size_t first, Args... rest) const;

  //! Access element without bounds checking.
  template <typename... Args>
  T &operator()(std::size_t first, Args... rest);

  //! Const access to element without bounds checking.
  template <typename... Args>
  T operator()(std::size_t first, Args... rest) const;

  //! Pointer to the first element of the underlying contiguous storage.
  T *data() { return m_data.data(); }
  const T *data() const { return m_data.data(); }

  //! Const reference to the underlying flat data vector.
  const std::vector<T> &vector() const { return m_data; }

  //! Number of rows (equivalent to size(0)).
  std::size_t rows() const { return size(0); }
  //! Number of columns (equivalent to size(1)).
  std::size_t cols() const { return size(1); }

  //! A view onto the ith row. Only defined for 2D arrays.
  ArrayView<T> row(std::size_t i);
  ArrayView<const T> row(std::size_t i) const;
  //! A view onto the jth column. Only defined for 2D arrays.
  ArrayView<T> col(std::size_t j);
  ArrayView<const T> col(std::size_t j) const;

  //! Iterator to the beginning
  auto begin() { return m_data.begin(); }
  //! Constant iterator to the beginning
  auto cbegin() const { return m_data.cbegin(); }
  //! Iterator to the end
  auto end() { return m_data.end(); }
  //! Constant iterator to the end
  auto cend() const { return m_data.cend(); }

  //! Reverse iterator to the beginning of the data
  auto rbegin() { return m_data.rbegin(); }
  //! Constant reverse iterator to the beginning of the data
  auto crbegin() const { return m_data.crbegin(); }
  //! Reverse iterator to the end of the data
  auto rend() { return m_data.rend(); }
  //! Constant reverse iterator to the end of the data
  auto crend() const { return m_data.crend(); }

  //! Element-wise arithmetic between arrays of identical shape (+, -, *, /).
  Array<T> &operator+=(const Array<T> &other);
  Array<T> &operator-=(const Array<T> &other);
  Array<T> &operator*=(const Array<T> &other);
  Array<T> &operator/=(const Array<T> &other);

  /*!
    @brief Element-wise scalar arithmetic (+, -, *, /).

    @details
    Array op T is provided for +, -, *, /.
    T op Array is provided only for * and +.
  */
  Array<T> &operator+=(const T &t);
  Array<T> &operator-=(const T &t);
  Array<T> &operator*=(const T &t);
  Array<T> &operator/=(const T &t);

private:
  // Calculates the cumulative_sizes array
  std::vector<std::size_t> calc_cumulative_size() const;

  // Unchecked index calculation
  template <typename... Args>
  std::size_t unchecked_index(std::size_t first, Args... rest) const;

  // Helper for unchecked index calculation
  template <typename... Args>
  std::size_t unchecked_index_impl(std::size_t dim, std::size_t first,
                                   Args... rest) const;

  // Checked index calculation
  template <typename... Args>
  std::size_t checked_index(std::size_t first, Args... rest) const;

  // Helper for checked index calculation
  template <typename... Args>
  std::size_t checked_index_impl(std::size_t dim, std::size_t first,
                                 Args... rest) const;
};

//==============================================================================
// Implementations
//==============================================================================
template <typename T>
std::vector<std::size_t> Array<T>::calc_cumulative_size() const {
  std::vector<std::size_t> cumulative_size;
  cumulative_size.reserve(m_sizes.size());
  for (std::size_t i = 0; i < m_Ndim; ++i) {
    cumulative_size.push_back(std::accumulate(m_sizes.cbegin() + long(i) + 1,
                                              m_sizes.cend(), 1ul,
                                              std::multiplies<std::size_t>()));
  }
  return cumulative_size;
}

template <typename T>
template <typename... Args>
Array<T>::Array(std::size_t first, Args... rest)
  : m_sizes({first, static_cast<std::size_t>(rest)...}),
    m_Ndim(m_sizes.size()),
    m_cumulative_sizes(calc_cumulative_size()),
    m_total_size(std::accumulate(m_sizes.cbegin(), m_sizes.cend(), 1ul,
                                 std::multiplies<std::size_t>())),
    m_data(m_total_size) {}

template <typename T>
template <typename... Args>
void Array<T>::resize(std::size_t first, Args... rest) {
  m_sizes = std::vector{first, static_cast<std::size_t>(rest)...};
  m_Ndim = m_sizes.size();
  m_cumulative_sizes = calc_cumulative_size();
  m_total_size = std::accumulate(m_sizes.cbegin(), m_sizes.cend(), 1ul,
                                 std::multiplies<std::size_t>());
  m_data.resize(m_total_size);
}

template <typename T>
template <typename... Args>
std::size_t Array<T>::unchecked_index(std::size_t first, Args... rest) const {
  return unchecked_index_impl(0, first, rest...);
}

template <typename T>
template <typename... Args>
std::size_t Array<T>::unchecked_index_impl(std::size_t dim, std::size_t first,
                                           Args... rest) const {
  if constexpr (sizeof...(rest) == 0) {
    return first * m_cumulative_sizes[dim];
  } else {
    return first * m_cumulative_sizes[dim] +
           unchecked_index_impl(dim + 1, rest...);
  }
}

template <typename T>
template <typename... Args>
std::size_t Array<T>::checked_index(std::size_t first, Args... rest) const {
  assert(sizeof...(rest) + 1 == m_Ndim &&
         "Number of arguments must match number of dimensions");
  return checked_index_impl(0, first, rest...);
}

template <typename T>
template <typename... Args>
std::size_t Array<T>::checked_index_impl(std::size_t dim, std::size_t first,
                                         Args... rest) const {
  assert(first < m_sizes[dim]);
  if constexpr (sizeof...(rest) == 0) {
    return first * m_cumulative_sizes[dim];
  } else {
    return first * m_cumulative_sizes[dim] +
           checked_index_impl(dim + 1, rest...);
  }
}

template <typename T>
template <typename... Args>
T &Array<T>::at(std::size_t first, Args... rest) {
  return m_data.at(checked_index(first, rest...));
}

template <typename T>
template <typename... Args>
T Array<T>::at(std::size_t first, Args... rest) const {
  return m_data.at(checked_index(first, rest...));
}

template <typename T>
template <typename... Args>
T &Array<T>::operator()(std::size_t first, Args... rest) {
  return m_data[unchecked_index(first, rest...)];
}

template <typename T>
template <typename... Args>
T Array<T>::operator()(std::size_t first, Args... rest) const {
  return m_data[unchecked_index(first, rest...)];
}

template <typename T>
Array<T> &Array<T>::operator+=(const Array<T> &other) {
  assert(m_sizes == other.m_sizes &&
         "Arithmetic only defined for equal-dimension arrays");
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] += other.m_data[i];
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator-=(const Array<T> &other) {
  assert(m_sizes == other.m_sizes &&
         "Arithmetic only defined for equal-dimension arrays");
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] -= other.m_data[i];
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator*=(const Array<T> &other) {
  assert(m_sizes == other.m_sizes &&
         "Arithmetic only defined for equal-dimension arrays");
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] *= other.m_data[i];
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator/=(const Array<T> &other) {
  assert(m_sizes == other.m_sizes &&
         "Arithmetic only defined for equal-dimension arrays");
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] /= other.m_data[i];
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator+=(const T &t) {
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] += t;
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator-=(const T &t) {
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] -= t;
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator*=(const T &t) {
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] *= t;
  }
  return *this;
}

template <typename T>
Array<T> &Array<T>::operator/=(const T &t) {
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    this->m_data[i] /= t;
  }
  return *this;
}

template <typename T>
ArrayView<T> Array<T>::row(std::size_t i) {
  assert(m_Ndim == 2 && "Row only defined for 2D array");
  return ArrayView(m_data.data() + i * m_sizes[1], m_sizes[1]);
}
template <typename T>
ArrayView<T> Array<T>::col(std::size_t j) {
  assert(m_Ndim == 2 && "Col only defined for 2D array");
  return ArrayView(m_data.data() + j, m_sizes[0], m_sizes[1]);
}
template <typename T>
ArrayView<const T> Array<T>::row(std::size_t i) const {
  assert(m_Ndim == 2 && "Row only defined for 2D array");
  return ArrayView(m_data.data() + i * m_sizes[1], m_sizes[1]);
}
template <typename T>
ArrayView<const T> Array<T>::col(std::size_t j) const {
  assert(m_Ndim == 2 && "Col only defined for 2D array");
  return ArrayView<const T>(m_data.data() + j, m_sizes[0], m_sizes[1]);
}

} // namespace qip
