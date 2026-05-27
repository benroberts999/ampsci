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

//! Type trait: true iff T is std::complex<U> for some U
template <typename T>
struct is_complex : std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};

//! Convenience variable template: is_complex_v<T> == is_complex<T>::value
template <typename T>
constexpr bool is_complex_v = is_complex<T>::value;

//==============================================================================
/*!
  @brief Linear algebra: matrices, vectors, views, and solvers.
  @details
  Provides dense matrix and vector types, non-owning views, and a set of
  linear algebra solvers wrapping LAPACK and GSL.

  **Core types** (Matrix.hpp, Vector.hpp):
  - `Matrix<T>` — owning row-major dense matrix; supports real and complex
    element types. Provides arithmetic, GSL interop, and element-wise access.
  - `Vector<T>` — owning 1D array; used for eigenvectors, right-hand sides,
    and eigenvalue arrays.

  **Views** (Matrix.hpp):
  - `View<T>` — non-owning strided view over a 1D segment of an array.
    Used to provide row and column access into a Matrix without copying.
    Obtained via `Matrix::row_view()` and `Matrix::column_view()`.
  - `Matrix_view<T>` — non-owning view of a 2D subblock of a Matrix.
    Obtained via `Matrix::submatrix_view()`.

  **Solvers** (Solvers.hpp):
  - `solve_Axeqb(A, b)` — solves Ax = b via LU decomposition (GSL).
  - `symmhEigensystem(A)` — all eigenvalues/vectors of a symmetric or
    Hermitian matrix (LAPACK `dsyev`/`zheev`).
  - `symmhEigensystem(A, n)` — first n eigenvalues/vectors (LAPACK `dsyevx`).
  - `symmhEigensystem(A, threshold)` — all eigenvalues below a threshold
    (LAPACK `dsyevr`).
  - `symmhEigensystem(A, B)` — generalised problem Av = eBv (LAPACK
    `dsygv`/`zhegv`).
  - `genEigensystem(A, sort)` — general non-symmetric real matrix (GSL).
*/
namespace LinAlg {

// Forward declaration (Matrix_view constructors accept Matrix<T> by reference)
template <typename T>
class Matrix;

/*!
  @brief Non-owning strided view onto a 1D segment of an array.
  @details
  Provides indexed read/write access into an existing array without ownership.
  The stride allows views over non-contiguous memory (e.g., a column of a
  row-major matrix).

  Used by `Matrix::row_view()` and `Matrix::column_view()`.

  @note Does not perform bounds checking in `operator[]`; use `at()` for
        checked access.
*/
template <typename T>
class View {
  std::size_t m_size;
  std::size_t m_stride;
  T *m_data;

public:
  //! Construct view over @p size elements starting at offset @p start, with @p stride
  View(T *data, std::size_t start, std::size_t size, std::size_t stride)
    : m_size(size), m_stride(stride), m_data(data + long(start)) {}

  //! Number of elements in the view
  std::size_t size() const { return m_size; }

  //! [] index access (no range checking), mutable
  T &operator[](std::size_t i) { return m_data[i * m_stride]; }
  //! [] index access (no range checking), const
  T operator[](std::size_t i) const { return m_data[i * m_stride]; }

  //! at(i): element access with range checking, mutable
  T &at(std::size_t i) {
    assert(i < m_size);
    return m_data[i * m_stride];
  }
  //! at(i): element access with range checking, const
  T at(std::size_t i) const {
    assert(i < m_size);
    return m_data[i * m_stride];
  }
  //! () index access with range checking, mutable
  T &operator()(std::size_t i) { return at(i); }
  //! () index access with range checking, const
  T operator()(std::size_t i) const { return at(i); }

  //! Raw pointer to first viewed element
  T *data() { return m_data; }
};

//==============================================================================
/*!
  @brief Non-owning 2D view onto a Matrix.
  @details
  Provides row/column element access to an existing matrix buffer without
  ownership or resize capability. Supports both mutable (`Matrix_view<T>`)
  and read-only (`Matrix_view<const T>`) access.

  Implicitly constructible from `Matrix<T>` (mutable view) and from
  `const Matrix<T>` (const view only, via SFINAE).

  @note Iterators (`begin()`/`end()`) yield raw pointers into the flat buffer.
*/
template <typename T>
class Matrix_view {
  std::size_t m_rows;
  std::size_t m_cols;
  T *m_data;

public:
  //! Construct from raw pointer and dimensions (no ownership)
  Matrix_view(T *data, std::size_t rows, std::size_t cols)
    : m_rows(rows), m_cols(cols), m_data(data) {}

  //! Implicit conversion from mutable Matrix (works for both mutable and const view)
  Matrix_view(Matrix<std::remove_const_t<T>> &m)
    : m_rows(m.rows()), m_cols(m.cols()), m_data(m.data()) {}

  //! Implicit conversion from const Matrix (only for Matrix_view<const T>)
  template <typename U = T, typename = std::enable_if_t<std::is_const_v<U>>>
  Matrix_view(const Matrix<std::remove_const_t<T>> &m)
    : m_rows(m.rows()), m_cols(m.cols()), m_data(m.data()) {}

  //! Return rows [major index size]
  std::size_t rows() const { return m_rows; }
  //! Return columns [minor index size]
  std::size_t cols() const { return m_cols; }
  //! Return rows*columns [total array size]
  std::size_t size() const { return m_rows * m_cols; }
  //! Bool - is empty? (same as size==0)
  bool empty() const { return m_rows == 0 || m_cols == 0; }

  //! Raw pointer to first element, mutable
  T *data() { return m_data; }
  //! Raw pointer to first element, const
  const T *data() const { return m_data; }

  //! [] index access (no range checking). [i] returns pointer to ith row, mutable
  T *operator[](std::size_t i) { return m_data + i * m_cols; }
  //! [] index access (no range checking). [i] returns const pointer to ith row
  const T *operator[](std::size_t i) const { return m_data + i * m_cols; }

  //! at(i,j): element access with range checking, mutable
  T &at(std::size_t i, std::size_t j) {
    assert(i < m_rows && j < m_cols);
    return m_data[i * m_cols + j];
  }
  //! at(i,j): element access with range checking, const
  T at(std::size_t i, std::size_t j) const {
    assert(i < m_rows && j < m_cols);
    return m_data[i * m_cols + j];
  }
  //! at(i,j): element access with range checking, const ref
  const T &atc(std::size_t i, std::size_t j) const {
    assert(i < m_rows && j < m_cols);
    return m_data[i * m_cols + j];
  }

  //! () index access with range checking, mutable
  T &operator()(std::size_t i, std::size_t j) { return at(i, j); }
  //! () index access with range checking, const
  T operator()(std::size_t i, std::size_t j) const { return at(i, j); }

  //! Iterators over flat data buffer
  auto begin() { return m_data; }
  auto end() { return m_data + m_rows * m_cols; }
  auto cbegin() const { return m_data; }
  auto cend() const { return m_data + m_rows * m_cols; }
};

//==============================================================================
/*!
  @brief Row-major dense matrix with arithmetic and linear algebra support.
  @details
  Owns its data as a contiguous `std::vector<T>` stored in row-major order.
  Supports scalar types `T = float`, `double`, `std::complex<float>`, and
  `std::complex<double>`.

  Provides:
  - Element access via `operator[]` (no bounds checking) and `at()`/`operator()` (bounds-checked).
  - Basic arithmetic: addition, subtraction, scalar multiply/divide, and matrix
    multiplication via BLAS.
  - Linear algebra via GSL: `determinant()`, `inverse()`, `transpose()`.
  - Conversion utilities: `conj()`, `real()`, `imag()`, `complex()`.
  - Non-owning views: `row_view()`, `column_view()`, `matrix_view()`.
  - GSL interop: `as_gsl_view()`.

  @note The scalar `+=` and `-=` operators treat the scalar as a multiple of
        the identity: `M += a` adds `a` to all diagonal elements. Only valid
        for square matrices.
*/
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

  //! Bool - is empty? (same as size==0)
  bool empty() const { return m_data.empty(); }

  //! Pointer to first element; for std::complex<T> this is complex<T>*, not T*
  T *data() { return m_data.data(); }
  //! Const pointer to first element; for std::complex<T> this is const complex<T>*, not const T*
  const T *data() const { return m_data.data(); }

  //============================================================================

  //! [] index access (no range checking). [i] returns pointer to ith row, mutable
  T *operator[](std::size_t i) { return &(m_data[i * m_cols]); }
  //! [] index access (no range checking). [i] returns const pointer to ith row
  const T *operator[](std::size_t i) const { return &(m_data[i * m_cols]); }

  //! at(i,j): element access with range checking, mutable
  T &at(std::size_t row_i, std::size_t col_j) {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }
  //! at(i,j): element access with range checking, const
  T at(std::size_t row_i, std::size_t col_j) const {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }

  //! at(i,j): element access with range checking, const ref
  const T &atc(std::size_t row_i, std::size_t col_j) const {
    assert(row_i < m_rows && col_j < m_cols);
    return m_data[row_i * m_cols + col_j];
  }

  //! () index access with range checking, mutable
  T &operator()(std::size_t i, std::size_t j) { return at(i, j); }
  //! () index access with range checking, const
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

  //! Returns a mutable view of the entire matrix (no ownership, no resize)
  [[nodiscard]] Matrix_view<T> matrix_view() {
    return Matrix_view<T>(m_data.data(), m_rows, m_cols);
  }
  //! Returns an immutable view of the entire matrix (no ownership, no resize)
  [[nodiscard]] Matrix_view<const T> matrix_view() const {
    return Matrix_view<const T>(m_data.data(), m_rows, m_cols);
  }

  //============================================================================
  /*!
    @brief Returns a GSL matrix view for use with GSL functions (no copy).
    @details
    Returns a `gsl_matrix_view` (or `_float_view`, `_complex_view`,
    `_complex_float_view`) wrapping the matrix data in place. Access the
    underlying GSL matrix via `.matrix` on the returned object. Allows
    use of all GSL/BLAS built-in functions without any data copy.

    Supported types: `double`, `float`, `std::complex<double>`,
    `std::complex<float>`.

    @note Both the `Matrix` and the returned view must remain in scope
          during use; the view is non-owning.
    @warning Fails at runtime (assert) for unsupported scalar types.
  */
  [[nodiscard]] auto as_gsl_view();

  //! Const GSL matrix view; see as_gsl_view()
  [[nodiscard]] auto as_gsl_view() const;

  //============================================================================
  // Basic matrix operations:
  //============================================================================
  //============================================================================

  /*!
    @brief Returns the determinant via LU decomposition (GSL).
    @details
    Makes an internal copy of the matrix, performs LU decomposition, and
    returns the determinant. The original matrix is not modified.

    @return Determinant of the matrix.
    @note Only available for `T = double` or `T = std::complex<double>`.
    @warning Only defined for square matrices; asserts otherwise.
  */
  [[nodiscard]] T determinant() const;

  /*!
    @brief Inverts the matrix in place via LU decomposition (GSL).
    @details
    Performs LU decomposition on a copy, then overwrites `*this` with its
    inverse using `gsl_linalg_LU_invert` (or the complex equivalent).

    @return Reference to `*this` (the inverted matrix).
    @note Only available for `T = double` or `T = std::complex<double>`.
    @warning Only defined for square matrices; asserts otherwise.
  */
  Matrix<T> &invert_in_place();

  /*!
    @brief Returns inverse of the matrix; original is unchanged.
    @details
    Copies `*this` and calls `invert_in_place()` on the copy.

    @note Only available for `T = double` or `T = std::complex<double>`.
  */
  [[nodiscard]] Matrix<T> inverse() const {
    auto inv = *this;
    return inv.invert_in_place();
  }

  /*!
    @brief Returns the transpose of the matrix.
    @details
    Allocates and returns a new `Matrix<T>` of size `cols() x rows()` with
    elements transposed. Uses GSL `gsl_matrix_transpose_memcpy` for `double`,
    `float`, and their complex variants; falls back to a naive loop for other
    types.

    @return New matrix of dimensions `cols() x rows()`.
  */
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
  /*!
    @brief Elementwise multiply in place: M_ij *= a_ij
    @details
    Each element of `*this` is multiplied by the corresponding element of @p a.
    Matrices must have the same dimensions.

    @param a  Matrix whose elements multiply `*this`.
    @return   Reference to `*this`.
    @warning  Dimensions must match; asserts otherwise.
  */
  Matrix<T> &mult_elements_by(const Matrix<T> &a);

  /*!
    @brief Returns new matrix C with C_ij = A_ij * B_ij (elementwise product).
    @details Calls `a.mult_elements_by(b)` on a copy of @p a.
  */
  [[nodiscard]] friend Matrix<T> mult_elements(Matrix<T> a,
                                               const Matrix<T> &b) {
    return a.mult_elements_by(b);
  }

  //============================================================================
  //! In-place elementwise addition; dimensions must match
  Matrix<T> &operator+=(const Matrix<T> &rhs);
  //! In-place elementwise subtraction; dimensions must match
  Matrix<T> &operator-=(const Matrix<T> &rhs);
  //! In-place scalar multiply: M_ij *= x
  Matrix<T> &operator*=(const T x);
  //! In-place scalar divide: M_ij /= x
  Matrix<T> &operator/=(const T x);

  //============================================================================
  //! Elementwise addition: returns M + N
  [[nodiscard]] friend Matrix<T> operator+(Matrix<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs += rhs);
  }
  //! Elementwise subtraction: returns M - N
  [[nodiscard]] friend Matrix<T> operator-(Matrix<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs -= rhs);
  }
  //! Scalar multiply: returns x * M
  [[nodiscard]] friend Matrix<T> operator*(const T x, Matrix<T> rhs) {
    return (rhs *= x);
  }
  //! Scalar multiply: returns M * x
  [[nodiscard]] friend Matrix<T> operator*(Matrix<T> lhs, const T x) {
    return (lhs *= x);
  }
  //! Scalar divide: returns M / x
  [[nodiscard]] friend Matrix<T> operator/(Matrix<T> lhs, const T x) {
    return (lhs /= x);
  }

  //============================================================================
  //! M += aI: adds a to diagonal elements (scalar treated as a*Identity)
  Matrix<T> &operator+=(T aI);
  //! M -= aI: subtracts a from diagonal elements (scalar treated as a*Identity)
  Matrix<T> &operator-=(T aI);

  //! Returns M + aI (adds a to diagonal)
  [[nodiscard]] friend Matrix<T> operator+(Matrix<T> M, T aI) {
    return (M += aI);
  }
  //! Returns M - aI (subtracts a from diagonal)
  [[nodiscard]] friend Matrix<T> operator-(Matrix<T> M, T aI) {
    return (M -= aI);
  }

  //============================================================================
  //! Matrix multiplication: C_ij = sum_k A_ik*B_kj
  template <typename U>
  friend Matrix<U> operator*(const Matrix<U> &a, const Matrix<U> &b);

  //! Prints matrix elements row by row to output stream
  template <typename U>
  friend std::ostream &operator<<(std::ostream &os, const Matrix<U> &a);
};

//==============================================================================
//==============================================================================
//==============================================================================

//! Converts bool to CBLAS_TRANSPOSE enum (CblasTrans if true, CblasNoTrans if false)
inline CBLAS_TRANSPOSE to_cblas_trans(bool trans) {
  return trans ? CblasTrans : CblasNoTrans;
}

/*!
  @brief Matrix multiplication C = op(A) * op(B) via CBLAS GEMM (row-major).
  @details
  Computes \f$ C = \mathrm{op}(A)\,\mathrm{op}(B) \f$ where
  \f$ \mathrm{op}(X) = X \f$ or \f$ X^T \f$ depending on @p trans_A /
  @p trans_B. @p c is overwritten with the result.

  Wraps CBLAS `gemm` for `double`, `float`, `std::complex<double>`, and
  `std::complex<float>`.

  @param a        Left-hand matrix.
  @param b        Right-hand matrix.
  @param[out] c   Output matrix; must be pre-allocated to the correct size.
  @param trans_A  If true, use \f$A^T\f$ in the product.
  @param trans_B  If true, use \f$B^T\f$ in the product.
  @warning @p c must not alias @p a or @p b.
*/
template <typename T>
void GEMM(const Matrix<T> &a, const Matrix<T> &b, Matrix<T> *c,
          bool trans_A = false, bool trans_B = false);

//------------------------------------------------------------------------------

//! 5-matrix contraction for N*N matrices: M_ab = A_ai B_aj C_ij D_ib E_jb, with BLAS
/*!
  @details
  All input matrices are assumed square with the same dimension N.
  The contraction sums over indices i and j:
  \f[
    M_{ab} = \sum_{i,j} A_{ai}\,B_{aj}\,C_{ij}\,D_{ib}\,E_{jb}.
  \f]

  This implementation evaluates the contraction by looping over the
  shared index i and using a BLAS GEMM for the inner j-summation:
  \f[
    X^{(i)}_{jb} = C_{ij}\,E_{jb}, \qquad
    Y^{(i)}_{ab} = \sum_j B_{aj}\,X^{(i)}_{jb},
  \f]
  followed by the accumulation
  \f[
    M_{ab} \mathrel{+}= A_{ai}\,Y^{(i)}_{ab}\,D_{ib}.
  \f]

  Thus each i-iteration performs one dense matrix multiplication
  (B * X^{(i)}) via BLAS, while the remaining operations are
  elementwise scalings and accumulations.

  The routine allocates temporary work matrices X and Y of size N×N
  (reused across i) and writes the result into @p M, which must be
  preallocated and is overwritten on entry (i.e., initialised to zero).

  @note
  ~20% faster than Naiive implementation, BUT not easily parallelisable.

  @tparam T  Scalar type (float, double, std::complex<...>)
  @param A,B,C,D,E  Input N×N matrices
  @param[out] M     Output N×N matrix receiving the contraction. Must be previously sized
*/
template <typename T>
void PENTA_GEMM(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C,
                const Matrix<T> &D, const Matrix<T> &E, Matrix<T> *M);

//! M_ab = A_ai B_aj C_ij D_ib E_jb. Assuming all square matrices. Uses blas.
/*!
  @details
  Calls:
  \ref PENTA_GEMM(const Matrix<T>&, const Matrix<T>&,
                  const Matrix<T>&, const Matrix<T>&,
                  const Matrix<T>&, Matrix<T>*).
*/
template <typename T>
Matrix<T> PENTA_GEMM(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C,
                     const Matrix<T> &D, const Matrix<T> &E) {
  Matrix<T> M(A.rows(), A.cols());
  PENTA_GEMM<T>(A, B, C, D, E, &M);
  return M;
}

//-----------

//! 5-matrix contraction for N*N matrices: M_ab = A_ai B_aj C_ij D_ib E_jb, without BLAS, but with optional OpenMP parallelisation
/*! @details
  Computes the tensor contraction
  \f[
    M_{ab} = \sum_{i,j} A_{ai}\,B_{aj}\,C_{ij}\,D_{ib}\,E_{jb},
  \f]
  assuming all input matrices are square with the same dimension.

  This implementation uses explicit nested loops (no BLAS).
  If @tparam PARALLEL is true, OpenMP pragmas are enabled to
  parallelise the outer loops.

  For a BLAS-accelerated implementation, see
  \ref PENTA_GEMM(const Matrix<T>&, const Matrix<T>&,
                  const Matrix<T>&, const Matrix<T>&,
                  const Matrix<T>&, Matrix<T>*).
  
  Usually, this version will be significantly faster than the BLAS one.
  However, if it iself is being called in a parallel region, the BLAS one will be faster.

  @tparam T         Scalar type: double, float, std::complex
  @tparam PARALLEL  Enable OpenMP parallelisation when true (compile-time)
  @param A,B,C,D,E  Input N×N matrices
  @param[out] M     Output N×N matrix receiving the contraction. Must be previously sized
*/
template <typename T, bool PARALLEL = false>
void PENTA(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C,
           const Matrix<T> &D, const Matrix<T> &E, Matrix<T> *M);

//! 5-index contraction for N*N matrices: M_ab = A_ai B_aj C_ij D_ib E_jb, without BLAS, but with optional OpenMP parallelisation
//! @details calls above
template <typename T, bool PARALLEL = false>
Matrix<T> PENTA(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C,
                const Matrix<T> &D, const Matrix<T> &E) {
  Matrix<T> M(A.rows(), A.cols());
  PENTA<T, PARALLEL>(A, B, C, D, E, &M);
  return M;
}

//==============================================================================
//==============================================================================
//==============================================================================

//! Default relative tolerance for `equal()`: 1e-6 for float, 1e-12 for double
template <typename T>
constexpr auto myEps();

/*!
  @brief Compares two matrices element-wise to within a relative tolerance.
  @details
  Returns true iff all elements satisfy
  \f[
    |A_{ij} - B_{ij}| \leq |\varepsilon\,(A_{ij} + B_{ij})|.
  \f]
  Matrices of different dimensions compare as not equal.

  @param lhs  Left matrix.
  @param rhs  Right matrix.
  @param eps  Relative tolerance (default: `myEps<T>()`).
  @return `true` if all elements agree within @p eps; `false` otherwise.
*/
template <typename T>
bool equal(const Matrix<T> &lhs, const Matrix<T> &rhs, T eps = myEps<T>());

} // namespace LinAlg

#include "Matrix.ipp"
