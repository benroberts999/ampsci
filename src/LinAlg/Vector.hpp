#pragma once
#include "Matrix.hpp"

namespace LinAlg {

//******************************************************************************
//! Vector class (inherits from Matrix)
template <typename T = double> class Vector : public Matrix<T> {
public:
  //! Initialise a blank square matrix dimension*dimension, filled with 0
  Vector(std::size_t dimension) : Matrix<T>(dimension, 1) {}

  //! Initialise a matrix from initialiser list. {{},{},{}}. Each row must be
  //! same length
  Vector(std::initializer_list<T> l) : Matrix<T>(l.size(), 1, l) {}

  //! Initialise from std::vector using move()
  Vector(std::vector<T> &&v)
      : Matrix<T>(v.size(), 1, std::forward<std::vector<T>>(v)) {}

  //! Initialise from std::vector by copy
  Vector(const std::vector<T> &v) : Matrix<T>(v.size(), 1, v) {}

  //! Initialise from Matrix by move
  Vector(const Matrix<T> &&m) : Matrix<T>(std::move(m)) {
    assert(m.cols() == 1 && "Can only convert Matrix to Vector if matrix has 1 "
                            "column. Traspose first?");
  }

  //! Initialise from Matrix by copy
  Vector(const Matrix<T> &m) : Matrix<T>(m) {
    assert(m.cols() == 1 && "Can only convert Matrix to Vector if matrix has 1 "
                            "column. Traspose first?");
  }

  //****************************************************************************

  //! [] index access (with no range checking). [i][j] returns ith row, jth col
  T operator[](std::size_t i) const { return this->data()[i]; }
  //! As above, but const
  T &operator[](std::size_t i) { return this->data()[i]; }
  //! () index access (with range checking). (i,j) returns ith row, jth col
  T &at(std::size_t i) {
    assert(i < this->size());
    return this->data()[i];
  }
  //! As above, but const
  T at(std::size_t i) const {
    assert(i < this->size());
    return this->data()[i];
  }
  //! () index access (with range checking). (i,j) returns ith row, jth col
  T &operator()(std::size_t i) { return at(i); }
  //! As above, but const
  T operator()(std::size_t i) const { return at(i); }

  //****************************************************************************
  [[nodiscard]] Vector<T> conj() const;
  [[nodiscard]] auto real() const;
  [[nodiscard]] auto imag() const;
  [[nodiscard]] auto complex() const;

  Vector<T> transpose() const = delete;

  //****************************************************************************
  Vector<T> &operator+=(const Vector<T> &rhs);
  Vector<T> &operator-=(const Vector<T> rhs);
  Vector<T> &operator*=(const T x);
  Vector<T> &operator/=(const T x);

  [[nodiscard]] friend Vector<T> operator+(Vector<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs += rhs);
  }
  [[nodiscard]] friend Vector<T> operator-(Vector<T> lhs,
                                           const Matrix<T> &rhs) {
    return (lhs -= rhs);
  }
  [[nodiscard]] friend Vector<T> operator*(const T x, Vector<T> rhs) {
    return (rhs *= x);
  }
  [[nodiscard]] friend Vector<T> operator*(Vector<T> lhs, const T x) {
    return (lhs *= x);
  }
  [[nodiscard]] friend Vector<T> operator/(Vector<T> lhs, const T x) {
    return (lhs /= x);
  }

  //****************************************************************************
  //! Matrix*Vector multiplication: v_i = sum_j A_ij*B_j
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

  //! Inner product: = sum_i a_i*b_i
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

  //! Outer product:M_ij = a_i*b_j
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

  //****************************************************************************
  //! Returns gsl_vector_view (or _float_view, _complex_view,
  //! _complex_float_view). Call .matrix to use as a GSL matrix (no copy is
  //! involved). Allows one to use all GSL built-in functions. Note: non -
  //! owning pointer - matrix must remain in scope.
  [[nodiscard]] auto as_gsl_view();

  //! As above, but const
  [[nodiscard]] auto as_gsl_view() const;
};

} // namespace LinAlg

#include "Vector.ipp"
