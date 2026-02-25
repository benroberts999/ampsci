#pragma once

namespace LinAlg {

//==============================================================================
template <typename T>
class View {
  std::size_t m_size;
  std::size_t m_stride;
  T *m_data;

public:
  View(T *data, std::size_t start, std::size_t size, std::size_t stride)
    : m_size(size), m_stride(stride), m_data(data + long(start)) {}

  std::size_t size() const { return m_size; }

  //! [] index access (with no range checking). [i][j] returns ith row, jth col
  T &operator[](std::size_t i) { return m_data[i * m_stride]; }
  //! As above, but const
  T operator[](std::size_t i) const { return m_data[i * m_stride]; }

  //! () index access (with range checking). (i,j) returns ith row, jth col
  T &at(std::size_t i) {
    assert(i < m_size);
    return m_data[i * m_stride];
  }
  //! As above, but const
  T at(std::size_t i) const {
    assert(i < m_size);
    return m_data[i * m_stride];
  }
  //! () index access (with range checking). (i,j) returns ith row, jth col
  T &operator()(std::size_t i) { return at(i); }
  //! As above, but const
  T operator()(std::size_t i) const { return at(i); }

  T *data() { return m_data; }
};

//==============================================================================
//==============================================================================
//==============================================================================

//==============================================================================
// Returns the determinant. Uses GSL; via LU decomposition. Only works for
// double/complex<double>
template <typename T>
T Matrix<T>::determinant() const {
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

//==============================================================================
// Inverts the matrix, in place. Uses GSL; via LU decomposition. Only works
// for double/complex<double>.
template <typename T>
Matrix<T> &Matrix<T>::invert_in_place() {
  static_assert(
    std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>,
    "invert only works for Matrix<double> or Matrix<complex<double>>");

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

template <typename T>
Matrix<T> Matrix<T>::inverse() const {
  auto inverse = *this; // copy
  return inverse.invert_in_place();
}

//==============================================================================
template <typename T>
Matrix<T> Matrix<T>::transpose() const {
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

//==============================================================================
// Constructs a diagonal unit matrix (identity)
template <typename T>
Matrix<T> &Matrix<T>::make_identity() {
  assert(m_rows == m_cols && "Can only call make_identity() for square matrix");
  for (auto i = 0ul; i < m_rows; ++i) {
    for (auto j = 0ul; j < m_cols; ++j) {
      at(i, j) = i == j ? T(1) : T(0);
    }
  }
  return *this;
}
// Sets all elements to zero
template <typename T>
Matrix<T> &Matrix<T>::zero() {
  for (std::size_t i = 0; i < size(); ++i) {
    m_data[i] = T(0);
  }
  return *this;
}

//==============================================================================
template <typename T>
Matrix<T> Matrix<T>::conj() const {
  static_assert(is_complex_v<T>, "conj() only available for complex Matrix");
  std::vector<T> conj_data;
  conj_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    conj_data.push_back(std::conj(m_data[i]));
  }
  return Matrix<T>{m_rows, m_cols, std::move(conj_data)};
}

template <typename T>
Matrix<T> &Matrix<T>::conj_in_place() {
  static_assert(is_complex_v<T>, "conj() only available for complex Matrix");
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    m_data[i] = std::conj(m_data[i]);
  }
  return *this;
}
//------------------------------------------------------------------------------
template <typename T>
auto Matrix<T>::real() const {
  static_assert(is_complex_v<T>, "real() only available for complex Matrix");
  std::vector<typename T::value_type> real_data;
  real_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    real_data.push_back(std::real(m_data[i]));
  }
  return Matrix<typename T::value_type>{m_rows, m_cols, std::move(real_data)};
}
//------------------------------------------------------------------------------
template <typename T>
auto Matrix<T>::imag() const {
  static_assert(is_complex_v<T>, "imag() only available for complex Matrix");
  std::vector<typename T::value_type> imag_data;
  imag_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    imag_data.push_back(std::imag(m_data[i]));
  }
  return Matrix<typename T::value_type>{m_rows, m_cols, std::move(imag_data)};
}
//------------------------------------------------------------------------------
template <typename T>
auto Matrix<T>::complex() const {
  static_assert(!is_complex_v<T>, "complex() only available for real Matrix");
  // use move constructor to avoid default Matrix construction
  std::vector<std::complex<T>> new_data;
  new_data.reserve(m_data.size());
  for (std::size_t i = 0; i < m_data.size(); ++i) {
    new_data.push_back(m_data[i]);
  }
  return Matrix<std::complex<T>>{m_rows, m_cols, std::move(new_data)};
}

//==============================================================================
template <typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &rhs) {
  assert(rows() == rhs.rows() && cols() == rhs.cols() &&
         "Matrices must have same dimensions for addition");
  using namespace qip::overloads;
  this->m_data += rhs.m_data;
  return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &rhs) {
  assert(rows() == rhs.rows() && cols() == rhs.cols() &&
         "Matrices must have same dimensions for subtraction");
  using namespace qip::overloads;
  this->m_data -= rhs.m_data;
  return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator*=(const T x) {
  using namespace qip::overloads;
  this->m_data *= x;
  return *this;
}
template <typename T>
Matrix<T> &Matrix<T>::operator/=(const T x) {
  using namespace qip::overloads;
  this->m_data /= x;
  return *this;
}

//==============================================================================
// Matrix<T> += T : T assumed to be *Identity!
template <typename T>
Matrix<T> &Matrix<T>::operator+=(T aI) {
  // Adds 'a' to diagonal elements (Assume a*Ident)
  assert(m_rows == m_cols && "Can only call M+a for square matrix");
  for (auto i = 0ul; i < m_rows; ++i) {
    at(i, i) += aI;
  }
  return *this;
}
// Matrix<T> -= T : T assumed to be *Identity!
template <typename T>
Matrix<T> &Matrix<T>::operator-=(T aI) {
  // Adds 'a' to diagonal elements (Assume a*Ident)
  assert(m_rows == m_cols && "Can only call M-a for square matrix");
  for (auto i = 0ul; i < m_rows; ++i) {
    at(i, i) -= aI;
  }
  return *this;
}

//==============================================================================
template <typename T>
Matrix<T> &Matrix<T>::mult_elements_by(const Matrix<T> &a) {
  assert(rows() == a.rows() && cols() == a.cols() &&
         "Matrices must have same dimensions for mult_elements_by");
  for (auto i = 0ul; i < m_data.size(); ++i) {
    m_data[i] *= a.m_data[i];
  }
  return *this;
}

//==============================================================================
template <typename T>
[[nodiscard]] Matrix<T> operator*(const Matrix<T> &a, const Matrix<T> &b) {
  // https://www.gnu.org/software/gsl/doc/html/blas.html
  assert(a.cols() == b.rows() &&
         "Matrices a and b must have correct dimension for multiplication");
  Matrix<T> product(a.rows(), b.cols());

  GEMM(a, b, &product);

  return product;
}

//==============================================================================

// // Matrix multiplication C = A*B. C is overwritten with product, and must be correct size already
// template <typename T>
// void GEMM(const Matrix<T> &a, const Matrix<T> &b, Matrix<T> *c, bool trans_A,
//            bool trans_B) {
//   //
//   assert(a.cols() == b.rows() &&
//          "Matrices a and b must have correct dimension for multiplication");
//   assert(c->rows() == a.rows() && c->cols() == b.cols() &&
//          "Matrices c (output) must have correct dimension for multiplication "
//          "of a and b");

//   const auto a_gsl = a.as_gsl_view();
//   const auto b_gsl = b.as_gsl_view();
//   auto c_gsl = c->as_gsl_view();

//   const auto a_trans = trans_A ? CblasTrans : CblasNoTrans;
//   const auto b_trans = trans_B ? CblasTrans : CblasNoTrans;

//   if constexpr (std::is_same_v<T, double>) {
//     gsl_blas_dgemm(a_trans, b_trans, 1.0, &a_gsl.matrix, &b_gsl.matrix, 0.0,
//                    &c_gsl.matrix);
//   } else if constexpr (std::is_same_v<T, float>) {
//     gsl_blas_sgemm(a_trans, b_trans, 1.0f, &a_gsl.matrix, &b_gsl.matrix, 0.0f,
//                    &c_gsl.matrix);
//   } else if constexpr (std::is_same_v<T, std::complex<double>>) {
//     gsl_blas_zgemm(a_trans, b_trans, GSL_COMPLEX_ONE, &a_gsl.matrix,
//                    &b_gsl.matrix, GSL_COMPLEX_ZERO, &c_gsl.matrix);
//   } else if constexpr (std::is_same_v<T, std::complex<float>>) {
//     const gsl_complex_float one{1.0f, 0.0f};
//     const gsl_complex_float zero{0.0f, 0.0f};
//     gsl_blas_cgemm(a_trans, b_trans, one, &a_gsl.matrix, &b_gsl.matrix, zero,
//                    &c_gsl.matrix);
//   }
// }

//==============================================================================
inline CBLAS_TRANSPOSE to_cblas_trans(bool trans) {
  return trans ? CblasTrans : CblasNoTrans;
}

//------------------------------------------------------------------------------
template <typename T>
void GEMM(const Matrix<T> &a, const Matrix<T> &b, Matrix<T> *c, bool trans_A,
          bool trans_B) {
  assert(c);

  const auto ta = to_cblas_trans(trans_A);
  const auto tb = to_cblas_trans(trans_B);

  // Effective dimensions:
  // op(A): (trans_A ? a.cols x a.rows : a.rows x a.cols)
  // op(B): (trans_B ? b.cols x b.rows : b.rows x b.cols)
  const int A_rows = static_cast<int>(trans_A ? a.cols() : a.rows());
  const int A_cols = static_cast<int>(trans_A ? a.rows() : a.cols());
  const int B_rows = static_cast<int>(trans_B ? b.cols() : b.rows());
  const int B_cols = static_cast<int>(trans_B ? b.rows() : b.cols());

  // GEMM sizes: C = op(A) * op(B), where
  // M = rows(op(A)), N = cols(op(B)), K = cols(op(A)) = rows(op(B))
  const int M = A_rows;
  const int N = B_cols;
  const int K = A_cols;

  assert(A_cols == B_rows && "op(A) cols must equal op(B) rows");
  assert(static_cast<int>(c->rows()) == M && static_cast<int>(c->cols()) == N &&
         "Output matrix c must be sized MxN");

  // Row-major leading dimensions:
  // lda = number of columns in A's *storage* (i.e., a.cols()) regardless of trans
  // same for b, c
  const int lda = static_cast<int>(a.cols());
  const int ldb = static_cast<int>(b.cols());
  const int ldc = static_cast<int>(c->cols());

  if constexpr (std::is_same_v<T, double>) {
    cblas_dgemm(CblasRowMajor, ta, tb, M, N, K, 1.0, a.data(), lda, b.data(),
                ldb, 0.0, c->data(), ldc);

  } else if constexpr (std::is_same_v<T, float>) {
    cblas_sgemm(CblasRowMajor, ta, tb, M, N, K, 1.0f, a.data(), lda, b.data(),
                ldb, 0.0f, c->data(), ldc);

  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    const std::complex<double> alpha{1.0, 0.0};
    const std::complex<double> beta{0.0, 0.0};
    cblas_zgemm(CblasRowMajor, ta, tb, M, N, K, &alpha, a.data(), lda, b.data(),
                ldb, &beta, c->data(), ldc);

  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    const std::complex<float> alpha{1.0f, 0.0f};
    const std::complex<float> beta{0.0f, 0.0f};
    cblas_cgemm(CblasRowMajor, ta, tb, M, N, K, &alpha, a.data(), lda, b.data(),
                ldb, &beta, c->data(), ldc);

  } else {
    static_assert(!sizeof(T), "GEMM: unsupported scalar type");
  }
}

//==============================================================================
// M_ab = A_ai B_aj C_ij D_ib E_jb, using BLAS
template <typename T>
void PENTA_GEMM(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C,
                const Matrix<T> &D, const Matrix<T> &E, Matrix<T> *pM) {
  //
  const auto N = A.rows(); // assume all square
  assert(A.cols() == A.rows() && "Must be square");

  Matrix<T> X(N, N);
  Matrix<T> Y(N, N);
  auto &M = *pM;

  // M_ab = A_ai B_aj C_ij D_ib E_jb
  //      = A_ai B_aj X(i)_jb D_ib
  //      = A_ai Y(i)_ab D_ib
  // X(i)_jb = C_ij * E_j2;
  // Y(i)_aj = B_ij * X(i)_jb

  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      const auto cij = C[i][j];
      for (std::size_t b = 0; b < N; ++b) {
        X[j][b] = cij * E[j][b];
      }
    }
    GEMM(B, X, &Y);
    for (std::size_t a = 0; a < N; ++a) {
      for (std::size_t b = 0; b < N; ++b) {
        M[a][b] += A[a][i] * Y[a][b] * D[i][b];
      }
    }
  }
}

// M_ab = A_ai B_aj C_ij D_ib E_jb
template <typename T, bool PARALLEL>
void PENTA(const Matrix<T> &A, const Matrix<T> &B, const Matrix<T> &C,
           const Matrix<T> &D, const Matrix<T> &E, Matrix<T> *pM) {
  //
  const auto N = A.rows(); // assume all square
  assert(A.cols() == A.rows() && "Must be square");

  auto &M = *pM;

  // M_ab = A_ai B_aj C_ij D_ib E_jb
  if constexpr (PARALLEL) {

#pragma omp parallel for collapse(2)
    for (std::size_t a = 0; a < N; ++a) {
      for (std::size_t b = 0; b < N; ++b) {
        const T *Ba = &B[a][0];
        T Mab = T(0);
        for (std::size_t i = 0; i < N; ++i) {
          const auto AaiDib = A[a][i] * D[i][b];
          const T *Ci = &C[i][0];
          for (std::size_t j = 0; j < N; ++j) {
            Mab += AaiDib * Ba[j] * Ci[j] * E[j][b];
          }
        }
        M[a][b] = Mab;
      }
    }

  } else {

    for (std::size_t a = 0; a < N; ++a) {
      const T *Ba = &B[a][0];
      for (std::size_t b = 0; b < N; ++b) {
        T Mab = T(0);
        for (std::size_t i = 0; i < N; ++i) {
          const auto AaiDib = A[a][i] * D[i][b];
          const T *Ci = &C[i][0];
          for (std::size_t j = 0; j < N; ++j) {
            Mab += AaiDib * Ba[j] * Ci[j] * E[j][b];
          }
        }
        M[a][b] = Mab;
      }
    }
  }
}

//==============================================================================
template <typename T>
auto Matrix<T>::as_gsl_view() {
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
    assert(false && "as_gsl_view() only available for double/float (or complex "
                    "double/float)");
  }
}

template <typename T>
auto Matrix<T>::as_gsl_view() const {
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
    assert(false && "as_gsl_view() only for available double/float (or complex "
                    "double/float)");
  }
}

//==============================================================================
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

//==============================================================================
//==============================================================================
//==============================================================================

//==============================================================================
// Helper for equal()
template <typename T>
constexpr auto myEps() {
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
