#pragma once
#include <algorithm>
#include <array>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>
//
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <tuple>
#include <utility>

/*
Eigen probably better option!

To do:
write asign-from-vector/array classes (can I use move? Not for 2x2 vector!)
*/

namespace LinAlg {

//******************************************************************************
class SqMatrix {

public:
  const int n;

  // private: //make private, declare friends
  gsl_matrix *m = nullptr;
  gsl_matrix *m_LU = nullptr;
  gsl_permutation *perm = nullptr;
  int s_LU = 0;

public:
  SqMatrix(int in_n);
  template <typename T> SqMatrix(const std::initializer_list<T> &l);
  ~SqMatrix();
  // copy+/assignment: doesn't copy LU etc!
  SqMatrix(const SqMatrix &matrix);           // copy constructor;
  SqMatrix &operator=(const SqMatrix &other); // copy assignment

  // private: //make private, declare friends
  void LU_decompose();

public:
  void make_diag(double value = 1.0);
  void scale(double value);
  void clip_low(double value);
  void clip_high(double value);

  void make_symmetric();
  double check_symmetric();
  void print();

  [[nodiscard]] SqMatrix transpose() const;
  double determinant(); // changes m_LU
  void invert();        // nb: destructive!
  [[nodiscard]] SqMatrix inverse() const;

  double *operator[](int i) const;
  friend SqMatrix operator*(const SqMatrix &lhs, const SqMatrix &rhs);
  SqMatrix &operator+=(const SqMatrix rhs);
  friend SqMatrix operator+(SqMatrix lhs, const SqMatrix &rhs);
  SqMatrix &operator-=(const SqMatrix rhs);
  friend SqMatrix operator-(SqMatrix lhs, const SqMatrix &rhs);
  SqMatrix &operator*=(const double x);
  friend SqMatrix operator*(const double x, SqMatrix rhs);
};

//******************************************************************************
class Vector {
public:
  const int n;

  // private:
  gsl_vector *vec;

public:
  Vector(const int in_n) : n(in_n), vec(gsl_vector_alloc(n)) {}

  template <typename T>
  Vector(const std::initializer_list<T> &l)
      : n((int)l.size()), vec(gsl_vector_alloc(n)) {
    auto i = 0;
    for (auto el : l)
      vec->data[i++] = el;
  }

  ~Vector() { gsl_vector_free(vec); }

  Vector(const Vector &vector) // copy constructor
      : n(vector.n), vec(gsl_vector_alloc(vector.n)) {
    gsl_vector_memcpy(vec, vector.vec);
  }

  void clip_low(const double value) {
    for (int i = 0; i < n; ++i) {
      if (std::abs((*this)[i]) < value)
        (*this)[i] = 0.0;
    }
  }
  void clip_high(const double value) {
    for (int i = 0; i < n; ++i) {
      if (std::abs((*this)[i]) > value) {
        auto s = (*this)[i] > 0.0 ? 1 : -1;
        (*this)[i] = s * value;
      }
    }
  }
  void print() {
    for (int i = 0; i < n; ++i) {
      std::cout << (*this)[i] << "\n";
    }
  }

  Vector &operator=(const Vector &other) // copy assignment
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

  Vector &operator+=(const Vector rhs) {
    gsl_vector_add(vec, rhs.vec);
    return *this;
  }
  friend Vector operator+(Vector lhs, const Vector &rhs) {
    lhs += rhs;
    return lhs;
  }
  Vector &operator-=(const Vector rhs) {
    gsl_vector_sub(vec, rhs.vec);
    return *this;
  }
  friend Vector operator-(Vector lhs, const Vector &rhs) {
    lhs -= rhs;
    return lhs;
  }
  Vector &operator*=(const double x) {
    gsl_vector_scale(vec, x);
    return *this;
  }
  friend Vector operator*(const double x, Vector rhs) {
    rhs *= x;
    return rhs;
  }

  double &operator[](int i) const { return (vec->data[i]); }

  friend Vector operator*(const SqMatrix &Aij, const Vector &bj) {
    // ci = Aij bj
    Vector ci(bj.n);
    if (bj.n != Aij.n)
      std::abort(); // XXX write message!
    for (int i = 0; i < bj.n; ++i) {
      ci[i] = 0.0;
      for (int j = 0; j < bj.n; ++j) {
        ci[i] += Aij[i][j] * bj[j];
      }
    }
    return ci;
  }
}; // namespace LinAlg

inline Vector solve_Axeqb(SqMatrix &Am, const Vector &b) {
  Vector x(b.n);
  int s;
  Am.LU_decompose();
  gsl_linalg_LU_decomp(Am.m, Am.perm, &s);
  gsl_linalg_LU_solve(Am.m, Am.perm, b.vec, x.vec);
  return x;
}

//*****************************************************************************
//*****************************************************************************
[[nodiscard]] inline std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, bool sort = true) {
  // Solves Av = ev for eigenvalues e and eigenvectors v
  // for Real Symmetric Matrices using GSL:
  // https://www.gnu.org/software/gsl/doc/html/eigen.html#real-symmetric-matrices
  // Note: This destroys A and B matrices (see below for details).

  const auto n = A.n;
  std::pair<Vector, SqMatrix> eigen_vv = std::make_pair(n, n);
  auto &[e_values, e_vectors] = eigen_vv;

  // This function computes the eigenvalues and eigenvectors of the real
  // generalized symmetric-definite matrix pair (A, B), and stores them in eval
  // and evec respectively. The computed eigenvectors are normalized to have
  // unit magnitude. On output, B contains its Cholesky decomposition and A is
  // destroyed.
  gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(n);
  gsl_eigen_symmv(A.m, e_values.vec, e_vectors.m, work);
  gsl_eigen_symmv_free(work);

  if (sort)
    gsl_eigen_symmv_sort(e_values.vec, e_vectors.m, GSL_EIGEN_SORT_VAL_ASC);

  return eigen_vv;
}

//*****************************************************************************
//------------------------------------------------------------------------------
[[nodiscard]] inline std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, bool sort = true) {
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
      evec_I[i][j] = GSL_IMAG(z);
    }
  }

  gsl_eigen_nonsymmv_free(work);
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  return eigen_vv;
}

//------------------------------------------------------------------------------
[[nodiscard]] inline std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true) {
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
      evec_I[i][j] = GSL_IMAG(z);
    }
  }

  gsl_eigen_genv_free(work);
  gsl_vector_complex_free(alpha);
  gsl_vector_free(beta);
  gsl_matrix_complex_free(evec);

  return eigen_vv;
}

//------------------------------------------------------------------------------
[[nodiscard]] inline std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true) {
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

  return eigen_vv;
}

} // namespace LinAlg
