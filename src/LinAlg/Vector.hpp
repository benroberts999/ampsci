#pragma once
#include "Matrix.hpp"

namespace LinAlg {

//==============================================================================
/*!
  @brief Owning 1D array; inherits from Matrix<T> with a single column.
  @details
  Stores elements contiguously. Provides 1D index access via `[]`, `()`,
  and `at()`, plus arithmetic operators and GSL interop inherited from Matrix.

  Supports inner product (`a * b`), matrix-vector product (`A * v`), and
  outer product (`outer_product(a, b)`).

  **Vector vs View<T>:**
  Use `Vector<T>` when you need to own the data (e.g. storing a result,
  passing to a solver). Use `View<T>` (@ref LinAlg::View) for non-owning,
  zero-copy access into an existing array — e.g. a row or column of a Matrix,
  obtained via `Matrix::row_view()` or `Matrix::column_view()`.

  @tparam T Element type (default: `double`).
*/
template <typename T = double>
class Vector : public Matrix<T> {
public:
  //! Default construct: empty vector
  Vector() : Matrix<T>() {}

  //! Construct zero-initialised vector of length `dimension`
  Vector(std::size_t dimension) : Matrix<T>(dimension, 1) {}

  //! Construct from initialiser list: `Vector<double> v = {1.0, 2.0, 3.0};`
  Vector(std::initializer_list<T> l) : Matrix<T>(l.size(), 1, l) {}

  //! Construct from std::vector by move
  Vector(std::vector<T> &&v)
    : Matrix<T>(v.size(), 1, std::forward<std::vector<T>>(v)) {}

  //! Construct from std::vector by copy
  Vector(const std::vector<T> &v) : Matrix<T>(v.size(), 1, v) {}

  //! Construct from single-column Matrix by move
  Vector(const Matrix<T> &&m) : Matrix<T>(std::move(m)) {
    assert(m.cols() == 1 && "Can only convert Matrix to Vector if matrix has 1 "
                            "column. Traspose first?");
  }

  //! Construct from single-column Matrix by copy
  Vector(const Matrix<T> &m) : Matrix<T>(m) {
    assert(m.cols() == 1 && "Can only convert Matrix to Vector if matrix has 1 "
                            "column. Traspose first?");
  }

  //============================================================================

  //! Element access by index, no range checking, mutable
  T &operator[](std::size_t i) { return this->data()[i]; }
  //! Element access by index, no range checking, const
  T operator[](std::size_t i) const { return this->data()[i]; }
  //! Element access by index, with range checking, mutable
  T &at(std::size_t i) {
    assert(i < this->size());
    return this->data()[i];
  }
  //! Element access by index, with range checking, const
  T at(std::size_t i) const {
    assert(i < this->size());
    return this->data()[i];
  }
  //! Element access by index, with range checking, mutable
  T &operator()(std::size_t i) { return at(i); }
  //! Element access by index, with range checking, const
  T operator()(std::size_t i) const { return at(i); }

  //============================================================================

  //! Returns element-wise complex conjugate
  [[nodiscard]] Vector<T> conj() const;
  //! Returns real part of each element as a new Vector
  [[nodiscard]] auto real() const;
  //! Returns imaginary part of each element as a new Vector
  [[nodiscard]] auto imag() const;
  //! Returns a complex-valued copy of this Vector
  [[nodiscard]] auto complex() const;

  //! Transpose is not defined for Vector (deleted)
  Vector<T> transpose() const = delete;

  //============================================================================

  //! Addition assignment: element-wise `*this += rhs`
  Vector<T> &operator+=(const Vector<T> &rhs);
  //! Subtraction assignment: element-wise `*this -= rhs`
  Vector<T> &operator-=(const Vector<T> rhs);
  //! Scalar multiplication assignment: `*this *= x`
  Vector<T> &operator*=(const T x);
  //! Scalar division assignment: `*this /= x`
  Vector<T> &operator/=(const T x);

  //! Element-wise addition: `lhs + rhs`
  [[nodiscard]] friend Vector<T> operator+(Vector<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs += rhs);
  }
  //! Element-wise subtraction: `lhs - rhs`
  [[nodiscard]] friend Vector<T> operator-(Vector<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs -= rhs);
  }
  //! Scalar multiplication: `x * v`
  [[nodiscard]] friend Vector<T> operator*(const T x, Vector<T> rhs) {
    return (rhs *= x);
  }
  //! Scalar multiplication: `v * x`
  [[nodiscard]] friend Vector<T> operator*(Vector<T> lhs, const T x) {
    return (lhs *= x);
  }
  //! Scalar division: `v / x`
  [[nodiscard]] friend Vector<T> operator/(Vector<T> lhs, const T x) {
    return (lhs /= x);
  }

  //============================================================================

  //! Matrix-vector product: returns `A * b`, i.e. `v_i = sum_j A_ij * b_j`
  [[nodiscard]] friend Vector<T> operator*(const Matrix<T> &a,
                                           const Vector<T> &b) {
    // https://www.gnu.org/software/gsl/doc/html/blas.html
    assert(a.cols() == b.rows());
    Vector<T> product(b.rows());
    const auto a_gsl = a.as_gsl_view();
    const auto b_gsl = b.as_gsl_view();
    auto product_gsl = product.as_gsl_view();
    if constexpr (std::is_same_v<T, double>) {
      gsl_blas_dgemv(CblasNoTrans, 1.0, &a_gsl.matrix, &b_gsl.vector, 0.0,
                     &product_gsl.vector);
    } else if constexpr (std::is_same_v<T, float>) {
      gsl_blas_sgemv(CblasNoTrans, 1.0f, &a_gsl.matrix, &b_gsl.vector, 0.0f,
                     &product_gsl.vector);
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
      gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, &a_gsl.matrix,
                     &b_gsl.vector, GSL_COMPLEX_ZERO, &product_gsl.vector);
    } else if constexpr (std::is_same_v<T, std::complex<float>>) {
      const gsl_complex_float one{1.0f, 0.0f};
      const gsl_complex_float zero{0.0f, 0.0f};
      gsl_blas_cgemv(CblasNoTrans, one, &a_gsl.matrix, &b_gsl.vector, zero,
                     &product_gsl.vector);
    }

    return product;
  }

  //! Inner (dot) product: returns `sum_i a_i * b_i`
  [[nodiscard]] friend T operator*(const Vector<T> &a, const Vector<T> &b) {
    // https://www.gnu.org/software/gsl/doc/html/blas.html
    assert(a.rows() == b.rows());
    const auto a_gsl = a.as_gsl_view();
    const auto b_gsl = b.as_gsl_view();
    T product = T(0);
    if constexpr (std::is_same_v<T, double>) {
      gsl_blas_ddot(&a_gsl.vector, &b_gsl.vector, &product);
    } else if constexpr (std::is_same_v<T, float>) {
      gsl_blas_sdot(&a_gsl.vector, &b_gsl.vector, &product);
    } else if constexpr (std::is_same_v<T, std::complex<double>>) {
      gsl_blas_zdotu(&a_gsl.vector, &b_gsl.vector,
                     reinterpret_cast<gsl_complex *>(&product));
    } else if constexpr (std::is_same_v<T, std::complex<float>>) {
      gsl_blas_cdotu(&a_gsl.vector, &b_gsl.vector,
                     reinterpret_cast<gsl_complex_float *>(&product));
    }

    return product;
  }

  //! Outer product: returns matrix M with `M_ij = a_i * b_j`
  [[nodiscard]] friend Matrix<T> outer_product(const Vector<T> &a,
                                               const Vector<T> &b) {
    Matrix<T> op(a.rows(), b.rows());
    for (std::size_t i = 0; i < op.rows(); ++i) {
      for (std::size_t j = 0; j < op.cols(); ++j) {
        op[i][j] = a[i] * b[j];
      }
    }
    return op;
  }

  //============================================================================

  //! Returns a GSL vector view of the underlying data (no copy).
  //! The Vector must remain in scope for the lifetime of the view.
  [[nodiscard]] auto as_gsl_view();

  //! Returns a const GSL vector view of the underlying data (no copy).
  //! The Vector must remain in scope for the lifetime of the view.
  [[nodiscard]] auto as_gsl_view() const;
};

} // namespace LinAlg

#include "Vector.ipp"
