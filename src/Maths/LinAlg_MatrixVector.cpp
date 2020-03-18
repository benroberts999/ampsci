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

#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace LinAlg {

//******************************************************************************
//******************************************************************************
// class SqMatrix:
SqMatrix::SqMatrix(int in_n) : n(in_n), m(gsl_matrix_alloc(in_n, in_n)) {}

// template <typename T>
SqMatrix::SqMatrix(const std::initializer_list<double> &l)
    : n((int)std::sqrt((int)l.size())), m(gsl_matrix_alloc(n, n)) {
  auto i = 0;
  for (auto el : l)
    m->data[i++] = el;
}

SqMatrix::~SqMatrix() { //
  gsl_matrix_free(m);
  if (m_LU != nullptr)
    gsl_matrix_free(m_LU);
  if (perm != nullptr)
    gsl_permutation_free(perm);
}

SqMatrix::SqMatrix(const SqMatrix &matrix) // copy constructor
    : n(matrix.n), m(gsl_matrix_alloc(matrix.n, matrix.n)) {
  gsl_matrix_memcpy(m, matrix.m);
}

SqMatrix &SqMatrix::operator=(const SqMatrix &other) // copy assignment
{
  if (this->n != other.n) {
    std::cerr << "FAIL 35 in SqMatrix. Cant re-assign matrices of different "
                 "dimension: N_lhs = "
              << this->n << ", N_rhs = " << other.n << "\n\n";
    std::abort();
  }
  gsl_matrix_memcpy(m, other.m);
  return *this;
}

//------------------------------------------------------------------------------
void SqMatrix::make_diag(double value) {
  // value = 1.0 default
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      (*this)[i][j] = (i == j) ? value : 0.0;
    }
  }
}

void SqMatrix::scale(double value) {
  const int n2 = n * n;
  for (int i = 0; i < n2; i++)
    m->data[i] *= value;
}

void SqMatrix::clip_low(double value) {
  const int n2 = n * n;
  for (int i = 0; i < n2; i++) {
    if (std::abs(m->data[i]) < value)
      m->data[i] = 0.0;
  }
}
void SqMatrix::clip_high(double value) {
  const int n2 = n * n;
  for (int i = 0; i < n2; i++) {
    if (std::abs(m->data[i]) > value) {
      auto s = m->data[i] > 0 ? 1 : -1;
      m->data[i] = s * value;
    }
  }
}

void SqMatrix::make_symmetric() {
  *this += this->transpose();
  this->scale(0.5);
}

double SqMatrix::check_symmetric() {
  double worst = 0.0;
  auto AmATr = *this - this->transpose();
  for (int i = 0; i < n * n; i++) {
    auto val = std::abs(AmATr.m->data[i]);
    worst = (val > worst) ? val : worst;
  }
  return worst;
}

void SqMatrix::print() {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%8.1e ", (*this)[i][j]);
      // std::cout << (long)&((*this)[i][j]) << " ";
    }
    std::cout << "\n";
  }
}

//------------------------------------------------------------------------------
SqMatrix SqMatrix::transpose() const {
  SqMatrix mTr(this->n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      mTr[i][j] = (*this)[j][i];
    }
  }
  return mTr;
}

void SqMatrix::LU_decompose() {
  // this may do extra work..
  if (m_LU == nullptr)
    m_LU = gsl_matrix_alloc(n, n);
  if (perm == nullptr)
    perm = gsl_permutation_alloc(n);
  gsl_matrix_memcpy(m_LU, m);
  gsl_linalg_LU_decomp(m_LU, perm, &s_LU);
}

double SqMatrix::determinant() {
  LU_decompose();
  return gsl_linalg_LU_det(m_LU, s_LU);
}

void SqMatrix::invert() {
  // note: this is destuctive: matrix will be inverted
  // Usuing this method, hard not to be: LU decomp changes Matrix
  //(So, would require copy to avoid this)
  LU_decompose();
  gsl_linalg_LU_invert(m_LU, perm, m);
}

SqMatrix SqMatrix::inverse() const {
  auto inverse = *this;
  inverse.invert();
  return inverse;
}

//------------------------------------------------------------------------------
double *SqMatrix::operator[](int i) const { return &(m->data[i * n]); }

SqMatrix operator*(const SqMatrix &lhs, const SqMatrix &rhs) {
  auto n_min = std::min(lhs.n, rhs.n);
  SqMatrix product(n_min);
  for (int i = 0; i < n_min; ++i) {
    for (int j = 0; j < n_min; ++j) {
      double cij = 0.0;
      for (int k = 0; k < n_min; ++k) {
        cij += lhs[i][k] * rhs[k][j];
      }
      product[i][j] = cij;
    }
  }
  return product;
}

SqMatrix &SqMatrix::operator+=(const SqMatrix rhs) {
  int n2 = n * n;
  for (int i = 0; i < n2; i++)
    m->data[i] += rhs.m->data[i];
  return *this;
}
SqMatrix operator+(SqMatrix lhs, const SqMatrix &rhs) {
  lhs += rhs;
  return lhs;
}
SqMatrix &SqMatrix::operator-=(const SqMatrix rhs) {
  int n2 = n * n;
  for (int i = 0; i < n2; i++)
    m->data[i] -= rhs.m->data[i];
  return *this;
}
SqMatrix operator-(SqMatrix lhs, const SqMatrix &rhs) {
  lhs -= rhs;
  return lhs;
}
SqMatrix &SqMatrix::operator*=(const double x) {
  scale(x);
  return *this;
}
SqMatrix operator*(const double x, SqMatrix rhs) {
  rhs *= x;
  return rhs;
}

//******************************************************************************
//******************************************************************************
// class Vector
Vector::Vector(const int in_n) : n(in_n), vec(gsl_vector_alloc(n)) {}

template <typename T>
Vector::Vector(const std::initializer_list<T> &l)
    : n((int)l.size()), vec(gsl_vector_alloc(n)) {
  auto i = 0;
  for (auto el : l)
    vec->data[i++] = el;
}

Vector::Vector(const Vector &other) : n(other.n), vec(gsl_vector_alloc(n)) {
  gsl_vector_memcpy(vec, other.vec);
}

Vector &Vector::operator=(const Vector &other) // copy assignment
{
  if (this->n != other.n) {
    std::cerr << "FAIL 32 in Vector. Cant re-assign Vectors of different "
                 "dimension: N_lhs = "
              << this->n << ", N_rhs = " << other.n << "\n\n";
    std::abort();
  }
  gsl_vector_memcpy(vec, other.vec);
  return *this;
}

Vector::~Vector() { gsl_vector_free(vec); }

//------------------------------------------------------------------------------
void Vector::clip_low(const double value) {
  for (int i = 0; i < n; ++i) {
    if (std::abs((*this)[i]) < value)
      (*this)[i] = 0.0;
  }
}
void Vector::clip_high(const double value) {
  for (int i = 0; i < n; ++i) {
    if (std::abs((*this)[i]) > value) {
      auto s = (*this)[i] > 0.0 ? 1 : -1;
      (*this)[i] = s * value;
    }
  }
}
void Vector::print() {
  for (int i = 0; i < n; ++i) {
    std::cout << (*this)[i] << "\n";
  }
}

//------------------------------------------------------------------------------
double &Vector::operator[](int i) const { return (vec->data[i]); }

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
  for (int i = 0; i < bj.n; ++i) {
    ci[i] = 0.0;
    for (int j = 0; j < bj.n; ++j) {
      ci[i] += Aij[i][j] * bj[j];
    }
  }
  return ci;
}

double inner_produce(const Vector &a, const Vector &b) {
  auto ip = 0.0;
  for (auto i = 0; i < a.n; ++i) {
    ip += a[i] * b[i];
  }
  return ip;
}

double operator*(const Vector &a, const Vector &b) {
  return inner_produce(a, b);
}

SqMatrix outer_produce(const Vector &a, const Vector &b) {
  SqMatrix op(a.n);
  for (auto i = 0; i < a.n; ++i) {
    for (auto j = 0; j < a.n; ++j) {
      op[i][j] = a[i] * b[j];
    }
  }
  return op;
}

//******************************************************************************
//******************************************************************************
// Solve LinAlg equations:

//******************************************************************************
Vector solve_Axeqb(SqMatrix &Am, const Vector &b) {
  Vector x(b.n);
  Am.LU_decompose();
  // gsl_linalg_LU_decomp(Am.m, Am.perm, &s);          // XXX do twice?
  // gsl_linalg_LU_solve(Am.m, Am.perm, b.vec, x.vec); // use Am.m_LU ?
  gsl_linalg_LU_solve(Am.m_LU, Am.perm, b.vec, x.vec);
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

  for (int i = 0; i < n; ++i) {
    gsl_complex evali = gsl_vector_complex_get(eval, i);
    eval_R[i] = (GSL_REAL(evali));
    eval_I[i] = (GSL_IMAG(evali));
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
    for (int j = 0; j < n; ++j) {
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

  for (int i = 0; i < n; ++i) {
    gsl_complex alphai = gsl_vector_complex_get(alpha, i);
    double betai = gsl_vector_get(beta, i);
    eval_R[i] = (GSL_REAL(alphai) / betai);
    eval_I[i] = (GSL_IMAG(alphai) / betai);
    gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);
    for (int j = 0; j < n; ++j) {
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
