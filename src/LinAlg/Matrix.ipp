#pragma once

namespace LinAlg {

//******************************************************************************
// Returns the determinant. Uses GSL; via LU decomposition. Only works for
// double/complex<double>
template <typename T> T Matrix<T>::determinant() const {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, std::complex<double>>,
                "Determinant only works for double");

  assert(rows() == cols() && "Determinant only defined for square matrix");
  // Make a copy, since this is destructive. (Performs LU decomp)
  auto LU = *this; // will become LU decomposed version
  int sLU = 0;
  auto gsl_view = LU.as_gsl_view();
  gsl_permutation *permutn = gsl_permutation_alloc(rows());
  if constexpr (std::is_same_v<T, double>) {
    gsl_linalg_LU_decomp(&gsl_view.matrix, permutn, &sLU);
    gsl_permutation_free(permutn);
    return gsl_linalg_LU_det(&gsl_view.matrix, sLU);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    gsl_linalg_complex_LU_decomp(&gsl_view.matrix, permutn, &sLU);
    gsl_permutation_free(permutn);
    const auto gsl_cmplx = gsl_linalg_complex_LU_det(&gsl_view.matrix, sLU);
    // Can probably avoid this copy? doesn't really matter.
    return {GSL_REAL(gsl_cmplx), GSL_IMAG(gsl_cmplx)};
  }
}

//******************************************************************************
// Inverts the matrix, in place. Uses GSL; via LU decomposition. Only works
// for double/complex<double>.
template <typename T> Matrix<T> &Matrix<T>::invert() {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, std::complex<double>>,
                "invert only works for double");

  assert(rows() == cols() && "Inverse only defined for square matrix");
  int sLU = 0;
  // gsl_linalg_LU_decomp(m, permutn, &sLU);
  // gsl_linalg_LU_invx(m, permutn);
  // In-place inversion gsl_linalg_LU_invx added sometime after GSL v:2.1
  // Getafix only has 2.1 installed, so can't use this for now
  auto LU = *this; // copy! to be LU decomposed
  auto LU_gsl = LU.as_gsl_view();
  auto iverse_gsl = this->as_gsl_view();
  gsl_permutation *permutn = gsl_permutation_alloc(m_rows);
  if constexpr (std::is_same_v<T, double>) {
    gsl_linalg_LU_decomp(&LU_gsl.matrix, permutn, &sLU);
    gsl_linalg_LU_invert(&LU_gsl.matrix, permutn, &iverse_gsl.matrix);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    gsl_linalg_complex_LU_decomp(&LU_gsl.matrix, permutn, &sLU);
    gsl_linalg_complex_LU_invert(&LU_gsl.matrix, permutn, &iverse_gsl.matrix);
  }
  gsl_permutation_free(permutn);
  return *this;
}

template <typename T> Matrix<T> Matrix<T>::inverse() const {
  auto inverse = *this; // copy
  return inverse.invert();
}

//******************************************************************************
template <typename T> Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> Tr(m_cols, m_rows);
  if constexpr (std::is_same_v<T, double>) {
    auto Tr_gsl = Tr.as_gsl_view();
    const auto this_gsl = as_gsl_view();
    gsl_matrix_transpose_memcpy(&Tr_gsl.matrix, &this_gsl.matrix);
  } else if constexpr (std::is_same_v<T, float>) {
    auto Tr_gsl = Tr.as_gsl_view();
    const auto this_gsl = as_gsl_view();
    gsl_matrix_float_transpose_memcpy(&Tr_gsl.matrix, &this_gsl.matrix);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    auto Tr_gsl = Tr.as_gsl_view();
    const auto this_gsl = as_gsl_view();
    gsl_matrix_complex_transpose_memcpy(&Tr_gsl.matrix, &this_gsl.matrix);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    auto Tr_gsl = Tr.as_gsl_view();
    const auto this_gsl = as_gsl_view();
    gsl_matrix_complex_float_transpose_memcpy(&Tr_gsl.matrix, &this_gsl.matrix);
  } else {
    // backup, works for any type
    for (auto i = 0ul; i < Tr.rows(); ++i) {
      for (auto j = 0ul; j < Tr.cols(); ++j) {
        Tr[i][j] = (*this)[j][i];
      }
    }
  }
  return Tr;
}

//******************************************************************************
// Constructs a diagonal unit matrix (identity)
template <typename T> Matrix<T> &Matrix<T>::make_identity() {
  for (auto i = 0ul; i < m_rows; ++i) {
    for (auto j = 0ul; j < m_cols; ++j) {
      at(i, j) = i == j ? T(1) : T(0);
    }
  }
  return *this;
}
// Sets all elements to zero
template <typename T> Matrix<T> &Matrix<T>::zero() {
  for (std::size_t i = 0; i < size(); ++i) {
    m_data[i] = T(0);
  }
  return *this;
}
// M -> M + aI, for I=identity (add a to diag elements)
template <typename T> Matrix<T> &Matrix<T>::plusIdent(T a) {
  for (auto i = 0ul; i < std::min(m_rows, m_cols); ++i) {
    at(i, i) += a;
  }
  return *this;
}

//******************************************************************************
template <typename T> Matrix<T> Matrix<T>::conj() const {
  static_assert(is_complex_v<T>, "conj() only available for complex Matrix");
  std::vector<T> conj_data;
  conj_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    conj_data.push_back(std::conj(m_data[i]));
  }
  return Matrix<T>{m_rows, m_cols, std::move(conj_data)};
}
//------------------------------------------------------------------------------
template <typename T> auto Matrix<T>::real() const {
  static_assert(is_complex_v<T>, "real() only available for complex Matrix");
  std::vector<typename T::value_type> real_data;
  real_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    real_data.push_back(std::real(m_data[i]));
  }
  return Matrix<typename T::value_type>{m_rows, m_cols, std::move(real_data)};
}
//------------------------------------------------------------------------------
template <typename T> auto Matrix<T>::imag() const {
  static_assert(is_complex_v<T>, "imag() only available for complex Matrix");
  std::vector<typename T::value_type> imag_data;
  imag_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    imag_data.push_back(std::imag(m_data[i]));
  }
  return Matrix<typename T::value_type>{m_rows, m_cols, std::move(imag_data)};
}
//------------------------------------------------------------------------------
template <typename T> auto Matrix<T>::complex() const {
  static_assert(!is_complex_v<T>, "complex() only available for real Matrix");
  // use move constructor to avoid default Matrix construction
  std::vector<std::complex<T>> new_data;
  new_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    new_data.push_back(m_data[i]);
  }
  return Matrix<std::complex<T>>{m_rows, m_cols, std::move(new_data)};
}

//******************************************************************************
template <typename T> Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &rhs) {
  assert(rows() == rhs.rows() && cols() == rhs.cols());
  using namespace qip::overloads;
  this->m_data += rhs.m_data;
  return *this;
}
template <typename T> Matrix<T> &Matrix<T>::operator-=(const Matrix<T> rhs) {
  assert(rows() == rhs.rows() && cols() == rhs.cols());
  using namespace qip::overloads;
  this->m_data -= rhs.m_data;
  return *this;
}
template <typename T> Matrix<T> &Matrix<T>::operator*=(const T x) {
  using namespace qip::overloads;
  this->m_data *= x;
  return *this;
}
template <typename T> Matrix<T> &Matrix<T>::operator/=(const T x) {
  using namespace qip::overloads;
  this->m_data /= x;
  return *this;
}

//******************************************************************************
template <typename T>
Matrix<T> &Matrix<T>::mult_elements_by(const Matrix<T> &a) {
  assert(rows() == a.rows() && cols() == a.cols());
  for (auto i = 0ul; i < m_data.size(); ++i) {
    m_data[i] *= a.m_data[i];
  }
  return *this;
}

//******************************************************************************
template <typename T>
[[nodiscard]] Matrix<T> operator*(const Matrix<T> &a, const Matrix<T> &b) {
  // https://www.gnu.org/software/gsl/doc/html/blas.html
  assert(a.cols() == b.rows());
  Matrix<T> product(a.rows(), b.cols());
  const auto a_gsl = a.as_gsl_view();
  const auto b_gsl = b.as_gsl_view();
  auto product_gsl = product.as_gsl_view();
  if constexpr (std::is_same_v<T, double>) {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &a_gsl.matrix,
                   &b_gsl.matrix, 0.0, &product_gsl.matrix);
  } else if constexpr (std::is_same_v<T, float>) {
    gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1.0f, &a_gsl.matrix,
                   &b_gsl.matrix, 0.0f, &product_gsl.matrix);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, &a_gsl.matrix,
                   &b_gsl.matrix, GSL_COMPLEX_ZERO, &product_gsl.matrix);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    const gsl_complex_float one{1.0f, 0.0f};
    const gsl_complex_float zero{0.0f, 0.0f};
    gsl_blas_cgemm(CblasNoTrans, CblasNoTrans, one, &a_gsl.matrix,
                   &b_gsl.matrix, zero, &product_gsl.matrix);
  }

  return product;
}

//******************************************************************************
template <typename T> auto Matrix<T>::as_gsl_view() {
  if constexpr (std::is_same_v<T, double>) {
    return gsl_matrix_view_array(m_data.data(), m_rows, m_cols);
  } else if constexpr (std::is_same_v<T, float>) {
    return gsl_matrix_float_view_array(m_data.data(), m_rows, m_cols);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    // reinterpret_cast OK: cppreference.com/w/cpp/numeric/complex
    return gsl_matrix_complex_view_array(
        reinterpret_cast<double *>(m_data.data()), m_rows, m_cols);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    return gsl_matrix_complex_float_view_array(
        reinterpret_cast<float *>(m_data.data()), m_rows, m_cols);
  } else {
    assert(false &&
           "as_gsl_view only for double/float (or complex double/float)");
  }
}

template <typename T> auto Matrix<T>::as_gsl_view() const {
  if constexpr (std::is_same_v<T, double>) {
    return gsl_matrix_const_view_array(m_data.data(), m_rows, m_cols);
  } else if constexpr (std::is_same_v<T, float>) {
    return gsl_matrix_float_const_view_array(m_data.data(), m_rows, m_cols);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    return gsl_matrix_complex_const_view_array(
        reinterpret_cast<const double *>(m_data.data()), m_rows, m_cols);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    return gsl_matrix_complex_float_const_view_array(
        reinterpret_cast<const float *>(m_data.data()), m_rows, m_cols);
  } else {
    assert(false &&
           "as_gsl_view only for double/float (or complex double/float)");
  }
}

//******************************************************************************
template <typename T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &a) {
  for (auto i = 0ul; i < a.rows(); ++i) {
    for (auto j = 0ul; j < a.cols(); ++j) {
      os << a(i, j) << " ";
    }
    os << "\n";
  }
  os << "\n";
  return os;
}

//******************************************************************************
//******************************************************************************
//******************************************************************************
template <typename T> struct Matrix<T>::CollumnIterator {
  using iterator_category = std::forward_iterator_tag; //?
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = T *;   // or also value_type*
  using reference = T &; // or also value_type&
  CollumnIterator(pointer ptr, std::size_t rows) : m_ptr(ptr), m_rows(rows) {}

  reference operator*() const { return *m_ptr; }
  pointer operator->() { return m_ptr; }

  // Prefix increment
  CollumnIterator &operator++() {
    m_ptr += m_rows;
    return *this;
  }

  friend bool operator==(const CollumnIterator &a, const CollumnIterator &b) {
    return a.m_ptr == b.m_ptr;
  };
  friend bool operator!=(const CollumnIterator &a, const CollumnIterator &b) {
    return a.m_ptr != b.m_ptr;
  };

private:
  pointer m_ptr;
  std::size_t m_rows;
};

//******************************************************************************
// Helper for equal()
template <typename T> constexpr auto myEps() {
  if constexpr (std::is_same_v<T, float> ||
                std::is_same_v<T, std::complex<float>>) {
    return 1.0e-6f;
  } else if constexpr (std::is_same_v<T, double> ||
                       std::is_same_v<T, std::complex<double>>) {
    return 1.0e-12;
  } else {
    return 0;
  }
}

// Compares two matrices; returns true iff all elements compare relatively to
// better than eps
template <typename T>
bool equal(const Matrix<T> &lhs, const Matrix<T> &rhs, T eps) {
  if (lhs.rows() != rhs.rows())
    return false;
  if (lhs.cols() != rhs.cols())
    return false;
  for (auto i = 0ul; i < lhs.rows(); ++i) {
    for (auto j = 0ul; j < lhs.cols(); ++j) {
      // need abs on eps in case of complex
      if (std::abs(lhs(i, j) - rhs(i, j)) >
          std::abs(eps * (lhs(i, j) + rhs(i, j))))
        return false;
    }
  }
  return true;
}

} // namespace LinAlg
