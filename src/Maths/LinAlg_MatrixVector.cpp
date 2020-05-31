#include "LinAlg_MatrixVector.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace LinAlg {

//******************************************************************************
// class SqMatrix:
//******************************************************************************
SqMatrix::SqMatrix(std::size_t in_n)
    : n(in_n), m(gsl_matrix_alloc(in_n, in_n)) {}

// SqMatrix::SqMatrix(const std::initializer_list<double> &l)
//     : n(std::sqrt(l.size())), m(gsl_matrix_alloc(n, n)) {
//   auto i = 0;
//   for (auto el : l)
//     m->data[i++] = el;
// }

SqMatrix::~SqMatrix() { gsl_matrix_free(m); }

SqMatrix::SqMatrix(const SqMatrix &matrix) // copy constructor
    : n(matrix.n), m(gsl_matrix_alloc(matrix.n, matrix.n)) {
  gsl_matrix_memcpy(m, matrix.m);
}

SqMatrix &SqMatrix::operator=(const SqMatrix &other) // copy assignment
{
  // Check sizes!
  if (this != &other)
    gsl_matrix_memcpy(m, other.m);
  return *this;
}

//------------------------------------------------------------------------------
void SqMatrix::make_identity() { gsl_matrix_set_identity(this->m); }

void SqMatrix::zero() { gsl_matrix_set_zero(this->m); }

void SqMatrix::clip_low(double value) {
  const auto n2 = n * n;
  for (auto i = 0ul; i < n2; i++) {
    if (std::abs(m->data[i]) < value)
      m->data[i] = 0.0;
  }
}
void SqMatrix::clip_high(double value) {
  const auto n2 = n * n;
  for (auto i = 0ul; i < n2; i++) {
    if (std::abs(m->data[i]) > value) {
      m->data[i] = m->data[i] > 0.0 ? value : -value;
    }
  }
}

void SqMatrix::enforce_symmetric() {
  *this += this->transpose();
  (*this) *= 0.5;
}

double SqMatrix::check_symmetric() {
  double worst = 0.0;
  const auto AmATr = *this - this->transpose();
  for (auto i = 0ul; i < n * n; i++) {
    const auto val = std::abs(AmATr.m->data[i]);
    worst = (val > worst) ? val : worst;
  }
  return worst;
}

void SqMatrix::print() {
  for (auto i = 0ul; i < n; ++i) {
    for (auto j = 0ul; j < n; ++j) {
      printf("%8.1e ", (*this)[i][j]);
    }
    std::cout << "\n";
  }
}

//------------------------------------------------------------------------------
SqMatrix SqMatrix::transpose() const {
  SqMatrix mTr(this->n);
  gsl_matrix_transpose_memcpy(mTr.m, this->m);
  return mTr;
}

double SqMatrix::determinant() const {
  // expensive for this to not be destructive
  gsl_matrix *mLU = gsl_matrix_alloc(n, n);
  gsl_permutation *permutn = gsl_permutation_alloc(n);
  int sLU = 0;
  gsl_matrix_memcpy(mLU, m);
  gsl_linalg_LU_decomp(mLU, permutn, &sLU);
  auto det = gsl_linalg_LU_det(mLU, sLU);
  gsl_matrix_free(mLU);
  gsl_permutation_free(permutn);
  return det;
}

SqMatrix &SqMatrix::invert() {
  // note: this is destuctive: matrix will be inverted
  // uses LU decomposition
  gsl_permutation *permutn = gsl_permutation_alloc(n);
  int sLU = 0;
  // gsl_linalg_LU_decomp(m, permutn, &sLU);
  // gsl_linalg_LU_invx(m, permutn);
  // In-place inversion gsl_linalg_LU_invx added sometime after GSL v:2.1
  // Getafix only has 2.1 installed, so can't use this for now
  auto mLU = *this; // copy!
  gsl_linalg_LU_decomp(mLU.m, permutn, &sLU);
  gsl_linalg_LU_invert(mLU.m, permutn, m);
  gsl_permutation_free(permutn);
  return *this;
}

SqMatrix SqMatrix::inverse() const {
  auto inverse = *this;
  inverse.invert();
  return inverse;
}

//------------------------------------------------------------------------------
double *SqMatrix::operator[](int i) const {
  return &(m->data[std::size_t(i) * n]);
}
double *SqMatrix::operator[](std::size_t i) const { return &(m->data[i * n]); }

SqMatrix operator*(const SqMatrix &lhs, const SqMatrix &rhs) {
  // check matrix sizes?
  SqMatrix product(lhs.n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, lhs.m, rhs.m, 0.0, product.m);
  return product;
}

SqMatrix &SqMatrix::operator+=(const SqMatrix rhs) {
  gsl_matrix_add(this->m, rhs.m);
  return *this;
}
SqMatrix operator+(SqMatrix lhs, const SqMatrix &rhs) {
  lhs += rhs;
  return lhs;
}
SqMatrix &SqMatrix::operator-=(const SqMatrix rhs) {
  gsl_matrix_sub(this->m, rhs.m);
  return *this;
}
SqMatrix operator-(SqMatrix lhs, const SqMatrix &rhs) {
  lhs -= rhs;
  return lhs;
}
SqMatrix &SqMatrix::operator*=(const double x) {
  gsl_matrix_scale(this->m, x);
  return *this;
}
SqMatrix operator*(const double x, SqMatrix rhs) {
  rhs *= x;
  return rhs;
}

void SqMatrix::mult_elements_by(const SqMatrix &rhs) {
  gsl_matrix_mul_elements(this->m, rhs.m);
}
SqMatrix mult_elements(SqMatrix lhs, const SqMatrix &rhs) {
  gsl_matrix_mul_elements(lhs.m, rhs.m);
  return lhs;
}

//******************************************************************************
//******************************************************************************
// class Vector
//******************************************************************************

Vector::Vector(const std::size_t in_n) : n(in_n), vec(gsl_vector_alloc(n)) {}

template <typename T>
Vector::Vector(const std::initializer_list<T> &l)
    : n(l.size()), vec(gsl_vector_alloc(n)) {
  auto i = 0;
  for (auto el : l)
    vec->data[i++] = el;
}

Vector::Vector(const Vector &other) : n(other.n), vec(gsl_vector_alloc(n)) {
  gsl_vector_memcpy(vec, other.vec);
}

Vector &Vector::operator=(const Vector &other) {
  // copy assignment
  // Check dimensions?
  if (this != &other)
    gsl_vector_memcpy(vec, other.vec);
  return *this;
}

Vector::~Vector() { gsl_vector_free(vec); }

//------------------------------------------------------------------------------
void Vector::clip_low(const double value) {
  for (std::size_t i = 0; i < n; ++i) {
    if (std::abs((*this)[i]) < value)
      (*this)[i] = 0.0;
  }
}
void Vector::clip_high(const double value) {
  for (std::size_t i = 0; i < n; ++i) {
    if (std::abs((*this)[i]) > value) {
      auto s = (*this)[i] > 0.0 ? 1 : -1;
      (*this)[i] = s * value;
    }
  }
}
void Vector::print() {
  for (std::size_t i = 0; i < n; ++i) {
    std::cout << (*this)[i] << "\n";
  }
}

//------------------------------------------------------------------------------
double &Vector::operator[](int i) const { return (vec->data[i]); }
double &Vector::operator[](std::size_t i) const { return (vec->data[i]); }

Vector &Vector::operator+=(const Vector rhs) {
  gsl_vector_add(vec, rhs.vec);
  return *this;
}
Vector operator+(Vector lhs, const Vector &rhs) {
  lhs += rhs;
  return lhs;
}
Vector &Vector::operator-=(const Vector rhs) {
  gsl_vector_sub(vec, rhs.vec);
  return *this;
}
Vector operator-(Vector lhs, const Vector &rhs) {
  lhs -= rhs;
  return lhs;
}
Vector &Vector::operator*=(const double x) {
  gsl_vector_scale(vec, x);
  return *this;
}
Vector operator*(const double x, Vector rhs) {
  rhs *= x;
  return rhs;
}
Vector operator*(const SqMatrix &Aij, const Vector &bj) {
  Vector ci(bj.n);
  if (bj.n != Aij.n) {
    std::cerr << "\n Fail 283 in Vector: " << bj.n << "!=" << Aij.n << "\n";
    std::abort();
  }
  for (std::size_t i = 0; i < bj.n; ++i) {
    ci[i] = 0.0;
    for (std::size_t j = 0; j < bj.n; ++j) {
      ci[i] += Aij[i][j] * bj[j];
    }
  }
  return ci;
}

double inner_product(const Vector &a, const Vector &b) {
  auto ip = 0.0;
  for (auto i = 0ul; i < a.n; ++i) {
    ip += a[i] * b[i];
  }
  return ip;
}

double operator*(const Vector &a, const Vector &b) {
  return inner_product(a, b);
}

SqMatrix outer_product(const Vector &a, const Vector &b) {
  SqMatrix op(a.n);
  // assert that a.n = b.n !
  for (auto i = 0ul; i < a.n; ++i) {
    for (auto j = 0ul; j < a.n; ++j) {
      op[i][j] = a[i] * b[j];
    }
  }
  return op;
}

//******************************************************************************
//******************************************************************************
// class ComplexSqMatrix
//******************************************************************************

ComplexSqMatrix::ComplexSqMatrix(std::size_t in_n)
    : n(in_n), m(gsl_matrix_complex_alloc(n, n)) {}

ComplexSqMatrix::~ComplexSqMatrix() { gsl_matrix_complex_free(m); }

ComplexSqMatrix::ComplexSqMatrix(const ComplexSqMatrix &other)
    : n(other.n), m(gsl_matrix_complex_alloc(n, n)) {
  gsl_matrix_complex_memcpy(m, other.m);
}

ComplexSqMatrix &ComplexSqMatrix::operator=(const ComplexSqMatrix &other) {
  if (this != &other)
    gsl_matrix_complex_memcpy(m, other.m);
  return *this;
}

// Constructs a diagonal unit matrix
void ComplexSqMatrix::make_identity() {
  gsl_matrix_complex_set_identity(this->m);
}
// Sets all elements to zero
void ComplexSqMatrix::zero() { gsl_matrix_complex_set_zero(this->m); }

// Returns the transpose of matrix: not destructive
ComplexSqMatrix ComplexSqMatrix::transpose() const {
  ComplexSqMatrix mTr(this->n);
  gsl_matrix_complex_transpose_memcpy(mTr.m, this->m);
  return mTr;
}
// Inverts the matrix: nb: destructive
ComplexSqMatrix &ComplexSqMatrix::invert() {
  // note: this is destuctive: matrix will be inverted
  // uses LU decomposition
  gsl_permutation *permutn = gsl_permutation_alloc(n);
  int sLU = 0;
  // gsl_linalg_complex_LU_decomp(m, permutn, &sLU);
  // gsl_linalg_complex_LU_invx(m, permutn);
  // In-place inversion gsl_linalg_LU_invx added sometime after GSL v:2.1
  // Getafix only has 2.1 installed, so can't use this for now
  auto mLU = *this;
  gsl_linalg_complex_LU_decomp(mLU.m, permutn, &sLU);
  gsl_linalg_complex_LU_invert(mLU.m, permutn, m);
  gsl_permutation_free(permutn);
  return *this;
}
// Returns the inverce of matrix: not destructive
ComplexSqMatrix ComplexSqMatrix::inverse() const {
  auto inverse = *this;
  inverse.invert();
  return inverse;
}

// make a ComplexSqMatrix from a SqMatrix: C = x*mR, x is complex
ComplexSqMatrix ComplexSqMatrix::make_complex(const Complex<double> &x,
                                              const SqMatrix &mR) {
  ComplexSqMatrix mC(mR.n);
  for (auto i = 0ul; i < mR.n; i++) {
    for (auto j = 0ul; j < mR.n; j++) {
      gsl_matrix_complex_set(
          mC.m, i, j, gsl_complex_rect(x.re * mR[i][j], x.im * mR[i][j]));
    }
  }
  return mC;
}

Complex<double> ComplexSqMatrix::get_copy(std::size_t i, std::size_t j) const {
  // really just for testing..?
  gsl_complex val = gsl_matrix_complex_get(m, i, j);
  return {GSL_REAL(val), GSL_IMAG(val)};
}

// Get the real part (copy) of the complex matrix
SqMatrix ComplexSqMatrix::real() const {
  SqMatrix re(n);
  for (auto i = 0ul; i < n; i++) {
    for (auto j = 0ul; j < n; j++) {
      re[i][j] = GSL_REAL(gsl_matrix_complex_get(m, i, j));
    }
  }
  return re;
}
// // Get the imaginary part (copy) of the complex matrix
SqMatrix ComplexSqMatrix::imaginary() const {
  SqMatrix im(n);
  for (auto i = 0ul; i < n; i++) {
    for (auto j = 0ul; j < n; j++) {
      im[i][j] = GSL_IMAG(gsl_matrix_complex_get(m, i, j));
    }
  }
  return im;
}

void ComplexSqMatrix::mult_elements_by(const ComplexSqMatrix &rhs) {
  gsl_matrix_complex_mul_elements(this->m, rhs.m);
}
ComplexSqMatrix mult_elements(ComplexSqMatrix lhs, const ComplexSqMatrix &rhs) {
  gsl_matrix_complex_mul_elements(lhs.m, rhs.m);
  return lhs;
}

//  Multiply elements by constant
ComplexSqMatrix &ComplexSqMatrix::operator*=(const Complex<double> &x) {
  gsl_matrix_complex_scale(this->m, gsl_complex_rect(x.re, x.im));
  return *this;
}

ComplexSqMatrix operator*(const Complex<double> &x, ComplexSqMatrix rhs) {
  return (rhs *= x);
}
ComplexSqMatrix operator*(ComplexSqMatrix rhs, const Complex<double> &x) {
  return (rhs *= x);
}

ComplexSqMatrix operator*(const ComplexSqMatrix &x, const ComplexSqMatrix &y) {
  // These functions compute the matrix-matrix product and sum
  // C = \alpha op(A) op(B) + \beta C
  // where op(A) = A, A^T, A^H
  // for TransA = CblasNoTrans, CblasTrans, CblasConjTrans
  // and similarly for the parameter TransB
  ComplexSqMatrix result(x.n);
  // check for errors?
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, x.m, y.m,
                 GSL_COMPLEX_ZERO, result.m);
  return result;
}

// Add + subtract Complex matrices
ComplexSqMatrix &ComplexSqMatrix::operator+=(const ComplexSqMatrix &rhs) {
  gsl_matrix_complex_add(this->m, rhs.m);
  return *this;
}
ComplexSqMatrix &ComplexSqMatrix::operator-=(const ComplexSqMatrix &rhs) {
  gsl_matrix_complex_sub(this->m, rhs.m);
  return *this;
}
ComplexSqMatrix operator+(ComplexSqMatrix lhs, const ComplexSqMatrix &rhs) {
  return lhs += rhs;
}
ComplexSqMatrix operator-(ComplexSqMatrix lhs, const ComplexSqMatrix &rhs) {
  return lhs -= rhs;
}

//******************************************************************************
//******************************************************************************
// Solve LinAlg equations:

//******************************************************************************
Vector solve_Axeqb(const SqMatrix &Am, const Vector &b) {
  Vector x(b.n);

  gsl_matrix *Am_LU = gsl_matrix_alloc(Am.n, Am.n);
  gsl_permutation *Am_perm = gsl_permutation_alloc(Am.n);
  int sLU = 0;
  gsl_matrix_memcpy(Am_LU, Am.m);
  gsl_linalg_LU_decomp(Am_LU, Am_perm, &sLU);

  gsl_linalg_LU_solve(Am_LU, Am_perm, b.vec, x.vec);
  if constexpr (false) { //??
    // These functions apply an iterative improvement to x, the solution of A x
    // = b, from the precomputed LU decomposition of A into (LU, p). Additional
    // workspace of length N is required in work.
    gsl_vector *work = gsl_vector_alloc(Am.n);
    gsl_linalg_LU_refine(Am.m, Am_LU, Am_perm, b.vec, x.vec, work);
    gsl_vector_free(work);
  }

  gsl_matrix_free(Am_LU);
  gsl_permutation_free(Am_perm);
  return x;
}

//*****************************************************************************
// Eigensystems using GSL. NOTE: e-vectors are stored in COLUMNS (not rows) of
// matrix! Therefore, we transpose the matrix (duuumb)
//*****************************************************************************
std::pair<Vector, SqMatrix> realSymmetricEigensystem(SqMatrix &A, bool sort) {
  // Solves Av = ev for eigenvalues e and eigenvectors v
  // for Real Symmetric Matrices using GSL:
  // https://www.gnu.org/software/gsl/doc/html/eigen.html#real-symmetric-matrices
  // Note: This destroys A and B matrices (see below for details).

  const auto n = A.n;
  std::pair<Vector, SqMatrix> eigen_vv = std::make_pair(n, n);
  auto &[e_values, e_vectors] = eigen_vv;

  // This function computes the eigenvalues and eigenvectors of the real
  // generalized symmetric-definite matrix pair (A, B), and stores them in
  // eval and evec respectively. The computed eigenvectors are normalized to
  // have unit magnitude. On output, B contains its Cholesky decomposition and
  // A is destroyed.
  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv(A.m, e_values.vec, e_vectors.m, work);
  gsl_eigen_symmv_free(work);

  if (sort)
    gsl_eigen_symmv_sort(e_values.vec, e_vectors.m, GSL_EIGEN_SORT_VAL_ASC);

  auto tmp = e_vectors.transpose();
  e_vectors = tmp;

  return eigen_vv;
}

//------------------------------------------------------------------------------
std::pair<Vector, SqMatrix> realSymmetricEigensystem(SqMatrix &A, SqMatrix &B,
                                                     bool sort) {
  // Solves Av = eBv for eigenvalues e and eigenvectors v
  // for Real Generalized Symmetric-Definite Eigensystems using GSL:
  // https://www.gnu.org/software/gsl/doc/html/eigen.html#real-generalized-symmetric-definite-eigensystems
  // Note: This destroys A and B matrices (see below for details).

  const auto n = A.n;
  std::pair<Vector, SqMatrix> eigen_vv = std::make_pair(n, n);
  auto &[e_values, e_vectors] = eigen_vv;

  // This function computes the eigenvalues and eigenvectors of the real
  // generalized symmetric-definite matrix pair (A, B), and stores them in
  // eval and evec respectively. The computed eigenvectors are normalized to
  // have unit magnitude. On output, B contains its Cholesky decomposition and
  // A is destroyed.
  gsl_eigen_gensymmv_workspace *work = gsl_eigen_gensymmv_alloc(n);
  gsl_eigen_gensymmv(A.m, B.m, e_values.vec, e_vectors.m, work);
  gsl_eigen_gensymmv_free(work);

  if (sort)
    gsl_eigen_symmv_sort(e_values.vec, e_vectors.m, GSL_EIGEN_SORT_VAL_ASC);

  auto tmp = e_vectors.transpose();
  e_vectors = tmp;

  return eigen_vv;
}

//*****************************************************************************
std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, bool sort) {
  // Solves for Av = ev
  // for Real Nonsymmetric Matrices, using GSL:
  // https://www.gnu.org/software/gsl/doc/html/eigen.html#real-nonsymmetric-matrices
  // In general, e-values will be complex
  // Note: A is destroyed and should not be used afterwards!

  const auto n = A.n;

  gsl_eigen_nonsymmv_workspace *work = gsl_eigen_nonsymmv_alloc(n);
  gsl_vector_complex *eval = gsl_vector_complex_alloc(n);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n, n);

  gsl_eigen_nonsymmv_params(0, work); // I think this is not needed
  gsl_eigen_nonsymmv(A.m, eval, evec, work);

  if (sort)
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  //(can only sort by ABS here, unfortunately)

  std::pair<Vector, SqMatrix> eigen_vvR = std::make_pair(n, n);
  std::pair<Vector, SqMatrix> eigen_vvI = std::make_pair(n, n);

  std::tuple<Vector, Vector, SqMatrix, SqMatrix> eigen_vv =
      std::make_tuple(n, n, n, n);
  auto &[eval_R, eval_I, evec_R, evec_I] = eigen_vv;

  for (std::size_t i = 0; i < n; ++i) {
    gsl_complex evali = gsl_vector_complex_get(eval, i);
    eval_R[i] = (GSL_REAL(evali));
    eval_I[i] = (GSL_IMAG(evali));
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
    for (std::size_t j = 0; j < n; ++j) {
      gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
      evec_R[i][j] = GSL_REAL(z);
      evec_I[i][j] = GSL_IMAG(z); // already ok? check!
    }
  }

  gsl_eigen_nonsymmv_free(work);
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  return eigen_vv;
}

//------------------------------------------------------------------------------
std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort) {
  // Solves for Av = eBv
  // for Real Generalized Nonsymmetric Eigensystems, using GSL:
  // https://www.gnu.org/software/gsl/doc/html/eigen.html#real-generalized-nonsymmetric-eigensystems
  // In general, e-values will be complex
  // Note: A and B are destroyed and should not be used afterwards!

  const auto n = A.n;

  gsl_eigen_genv_workspace *work = gsl_eigen_genv_alloc(n);
  gsl_vector_complex *alpha = gsl_vector_complex_alloc(n);
  gsl_vector *beta = gsl_vector_alloc(n);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(n, n);

  // This function computes eigenvalues and right eigenvectors of the n-by-n
  // real generalized nonsymmetric matrix pair (A, B). The eigenvalues are
  // stored in (alpha, beta) and the eigenvectors are stored in evec. The
  // computed eigenvectors are normalized to have unit magnitude. On output,
  // (A, B) contains the generalized Schur form (S, T).
  // gsl_eigen_gen_params(0, 0, 0, work); // I think this is not needed
  gsl_eigen_genv(A.m, B.m, alpha, beta, evec, work);
  // eigen value is = alpha/beta

  if (sort)
    gsl_eigen_genv_sort(alpha, beta, evec, GSL_EIGEN_SORT_ABS_ASC);
  //(can only sort by ABS here, unfortunately)

  std::pair<Vector, SqMatrix> eigen_vvR = std::make_pair(n, n);
  std::pair<Vector, SqMatrix> eigen_vvI = std::make_pair(n, n);

  std::tuple<Vector, Vector, SqMatrix, SqMatrix> eigen_vv =
      std::make_tuple(n, n, n, n);
  auto &[eval_R, eval_I, evec_R, evec_I] = eigen_vv;

  for (std::size_t i = 0; i < n; ++i) {
    gsl_complex alphai = gsl_vector_complex_get(alpha, i);
    double betai = gsl_vector_get(beta, i);
    eval_R[i] = (GSL_REAL(alphai) / betai);
    eval_I[i] = (GSL_IMAG(alphai) / betai);
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
    for (std::size_t j = 0; j < n; ++j) {
      gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
      evec_R[i][j] = GSL_REAL(z);
      evec_I[i][j] = GSL_IMAG(z); // already ok? check!
    }
  }

  gsl_eigen_genv_free(work);
  gsl_vector_complex_free(alpha);
  gsl_vector_free(beta);
  gsl_matrix_complex_free(evec);

  return eigen_vv;
}

} // namespace LinAlg
