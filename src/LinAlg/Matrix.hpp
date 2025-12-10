#pragma once
#include "qip/Vector.hpp" // for std::vector overloads
#include <array>
#include <cassert>
#include <complex>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

template <typename T>
struct is_complex : std::false_type {};
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};
template <typename T>
constexpr bool is_complex_v = is_complex<T>::value;

//==============================================================================
//! Defines Matrix, Vector classes, and linear some algebra functions
namespace LinAlg {

//! Proved a "view" onto an array
template <typename T>
class View;

//==============================================================================
//! Matrix class; row-major
template <typename T = double>
class Matrix {

protected:
  std::size_t m_rows;
  std::size_t m_cols;
  std::vector<T> m_data{};

public:
  //! Default initialiser
  Matrix() : m_rows(0), m_cols(0) {}

  //! Initialise a blank matrix rows*cols, filled with 0
  Matrix(std::size_t rows, std::size_t cols)
      : m_rows(rows), m_cols(cols), m_data(rows * cols) {}

  //! Initialise a matrix rows*cols, filled with 'value'
  Matrix(std::size_t rows, std::size_t cols, const T &value)
      : m_rows(rows), m_cols(cols), m_data(rows * cols, value) {}

  //! Initialise a blank square matrix dimension*dimension, filled with 0
  // excplicit, since don't alow flot->int converions
  explicit Matrix(std::size_t dimension)
      : m_rows(dimension), m_cols(dimension), m_data(dimension * dimension) {}

  //! Initialise a matrix from initialiser list. {{},{},{}}. Each row must be
  //! same length
  Matrix(std::initializer_list<std::initializer_list<T>> ll)
      : m_rows(ll.size()), m_cols(ll.begin()->size()), m_data{} {
    // way to avoid copy?
    m_data.reserve(m_rows * m_cols);
    for (auto &l : ll) {
      m_data.insert(m_data.end(), l.begin(), l.end());
    }
  }

  //! Initialise a matrix from single initialiser list. {...}.
  Matrix(std::size_t rows, std::size_t cols, std::initializer_list<T> l)
      : m_rows(rows), m_cols(cols), m_data{l} {
    assert(m_data.size() == rows * cols &&
           "initializer_list must be rows*cols");
  }

  //! Initialise a matrix from single initialiser list. {...}.
  Matrix(std::size_t rows, std::size_t cols, std::vector<T> &&v)
      : m_rows(rows), m_cols(cols), m_data{std::forward<std::vector<T>>(v)} {
    assert(m_data.size() == rows * cols &&
           "initializer_list must be rows*cols");
  }
  //! Initialise a matrix from single initialiser list. {...}.
  Matrix(std::size_t rows, std::size_t cols, const std::vector<T> &v)
      : m_rows(rows), m_cols(cols), m_data{v} {
    assert(m_data.size() == rows * cols &&
           "initializer_list must be rows*cols");
  }

  //============================================================================
  //! Resizes matrix to new dimension; all values reset to default
  void resize(std::size_t rows, std::size_t cols) {
    m_rows = rows;
    m_cols = cols;
    m_data.assign(rows * cols, T{});
  }

  //! Resizes matrix to new dimension; all values reset to 'value'
  void resize(std::size_t rows, std::size_t cols, const T &value) {
    m_rows = rows;
    m_cols = cols;
    m_data.assign(rows * cols, value);
  }

  //============================================================================

  //! Return rows [major index size]
  std::size_t rows() const { return m_rows; }
  //! Return columns [minor index size]
  std::size_t cols() const { return m_cols; }
  //! Return rows*columns [total array size]
  std::size_t size() const { return m_data.size(); }

  //! Returns pointer to first element. Note: for std::complex<T>, this is a
  //! pointer to complex<T>, not T
  T *data() { return m_data.data(); }
  //! As above, but const
  const T *data() const { return m_data.data(); }

  //============================================================================

  //! [] index access (with no range checking). [i][j] returns ith row, jth col
  const T *operator[](std::size_t i) const { return &(m_data[i * m_cols]); }
  //! As above, but const
  T *operator[](std::size_t i) { return &(m_data[i * m_cols]); }

  //! () index access (with range checking). (i,j) returns ith row, jth col
  T &at(std::size_t row_i, std::size_t col_j) {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }
  //! const () index access (with range checking). (i,j) ith row, jth col
  T at(std::size_t row_i, std::size_t col_j) const {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }

  //! const ref () index access (with range checking). (i,j) ith row, jth col
  const T &atc(std::size_t row_i, std::size_t col_j) const {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }

  //! () index access (with range checking). (i,j) returns ith row, jth col
  T &operator()(std::size_t i, std::size_t j) { return at(i, j); }
  //! As above, but const
  T operator()(std::size_t i, std::size_t j) const { return at(i, j); }

  //============================================================================

  //! iterators for underlying std::vector (entire data)
  auto begin() { return m_data.begin(); }
  auto cbegin() const { return m_data.cbegin(); }
  auto end() { return m_data.end(); }
  auto cend() const { return m_data.cend(); }

  //! Returns raw c pointer to start of a row
  // [[deprecated]]
  const T *row(std::size_t row) const {
    return m_data.data() + long(row * m_cols);
  }

  //! Returns a mutable 'View' of a row
  [[nodiscard]] View<T> row_view(std::size_t row) {
    return View<T>(this->data(), row * m_cols, m_cols, 1ul);
  }
  //! Returns an immutable 'View' of a row
  [[nodiscard]] View<const T> row_view(std::size_t row) const {
    return View<const T>(this->data(), row * m_cols, m_cols, 1ul);
  }
  //! Returns a mutable 'View' of a column
  [[nodiscard]] View<T> column_view(std::size_t col) {
    return View<T>(this->data(), col, m_rows, m_rows);
  }
  //! Returns an immutable 'View' of a column
  [[nodiscard]] View<const T> column_view(std::size_t col) const {
    return View<const T>(this->data(), col, m_rows, m_rows);
  }

  //============================================================================
  //! Returns gsl_matrix_view (or _float_view, _complex_view,
  //! _complex_float_view). Call .matrix to use as a GSL matrix (no copy is
  //! involved). Allows one to use all GSL built-in functions. Note: non-owning
  //! pointer - matrix AND gsl_view must remain in scope.
  [[nodiscard]] auto as_gsl_view();

  //! As above, but const
  [[nodiscard]] auto as_gsl_view() const;

  //============================================================================
  // Basic matrix operations:
  //============================================================================
  //============================================================================

  //! Returns the determinant. Uses GSL; via LU decomposition. Only works for
  //! double/complex<double>
  [[nodiscard]] T determinant() const;

  //! Inverts the matrix, in place. Uses GSL; via LU decomposition. Only works
  //! for double/complex<double>.
  Matrix<T> &invert_in_place();

  //! Returns inverse of the matrix. Leaves original matrix intact. Uses GSL;
  //! via LU decomposition. Only works for double/complex<double>
  [[nodiscard]] Matrix<T> inverse() const;
  //! Returns transpose of matrix
  [[nodiscard]] Matrix<T> transpose() const;

  //============================================================================

  //! Constructs a diagonal unit matrix (identity), in place; only for square
  Matrix<T> &make_identity();
  //! Sets all elements to zero, in place
  Matrix<T> &zero();
  // //! M -> M + aI, for I=identity (adds a to diag elements), in place
  // Matrix<T> &plusIdent(T a = T(1));
  //============================================================================

  //! Returns conjugate of matrix
  [[nodiscard]] Matrix<T> conj() const;
  //! Returns real part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] auto real() const;
  //! Returns imag part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] auto imag() const;
  //! Converts a real to complex matrix (changes type; returns a complex matrix)
  [[nodiscard]] auto complex() const;

  //! Conjugates matrix, in place
  Matrix<T> &conj_in_place();

  //============================================================================
  //! Muplitplies all the elements by those of matrix a, in place: M_ij *= a_ij
  Matrix<T> &mult_elements_by(const Matrix<T> &a);

  //! Returns new matrix, C_ij = A_ij*B_ij
  [[nodiscard]] friend Matrix<T> mult_elements(Matrix<T> a,
                                               const Matrix<T> &b) {
    return a.mult_elements_by(b);
  }

  //============================================================================
  // Operator overloads: +,-, scalar */
  //! Overload standard operators: do what expected
  Matrix<T> &operator+=(const Matrix<T> &rhs);
  Matrix<T> &operator-=(const Matrix<T> &rhs);
  Matrix<T> &operator*=(const T x);
  Matrix<T> &operator/=(const T x);

  //============================================================================
  // nb: these are defined inline here to avoid ambiguous overload?
  [[nodiscard]] friend Matrix<T> operator+(Matrix<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs += rhs);
  }
  [[nodiscard]] friend Matrix<T> operator-(Matrix<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs -= rhs);
  }
  [[nodiscard]] friend Matrix<T> operator*(const T x, Matrix<T> rhs) {
    return (rhs *= x);
  }
  [[nodiscard]] friend Matrix<T> operator*(Matrix<T> lhs, const T x) {
    return (lhs *= x);
  }
  [[nodiscard]] friend Matrix<T> operator/(Matrix<T> lhs, const T x) {
    return (lhs /= x);
  }

  //============================================================================
  //! Matrix<T> += T : T assumed to be *Identity!
  Matrix<T> &operator+=(T aI);
  //! Matrix<T> -= T : T assumed to be *Identity!
  Matrix<T> &operator-=(T aI);

  [[nodiscard]] friend Matrix<T> operator+(Matrix<T> M, T aI) {
    return (M += aI);
  }
  [[nodiscard]] friend Matrix<T> operator-(Matrix<T> M, T aI) {
    return (M -= aI);
  }

  //============================================================================
  //! Matrix multiplication: C_ij = sum_k A_ik*B_kj
  template <typename U>
  friend Matrix<U> operator*(const Matrix<U> &a, const Matrix<U> &b);

  //============================================================================
  template <typename U>
  friend std::ostream &operator<<(std::ostream &os, const Matrix<U> &a);
};

//==============================================================================
//==============================================================================
//==============================================================================

//==============================================================================
//==============================================================================
//==============================================================================

template <typename T>
constexpr auto myEps();
//! Compares two matrices; returns true iff all elements compare relatively to
//! better than eps
template <typename T>
bool equal(const Matrix<T> &lhs, const Matrix<T> &rhs, T eps = myEps<T>());

} // namespace LinAlg

#include "Matrix.ipp"
