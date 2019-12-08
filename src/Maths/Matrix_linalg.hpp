#pragma once
#include <algorithm>
#include <array>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>
//
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

/*
Eigen probably better option!

To do:
write asign-from-vector/array classes (can I use move? Not for 2x2 vector!)
*/

namespace LinAlg {

//******************************************************************************
class SqMatrix {
  // Suddenly gets very slow around n=170.. why? When other method doesn't

public:
  const int n;

  // private:
  gsl_matrix *m = nullptr;
  gsl_matrix *m_LU = nullptr;
  gsl_permutation *perm = nullptr;
  int s_LU = 0;

public:
  SqMatrix(int in_n) : n(in_n), m(gsl_matrix_alloc(in_n, in_n)) {}

  template <typename T>
  SqMatrix(const std::initializer_list<T> &l)
      : n((int)std::sqrt((int)l.size())), m(gsl_matrix_alloc(n, n)) {
    auto i = 0;
    for (auto el : l)
      m->data[i++] = el;
  }

  ~SqMatrix() { //
    gsl_matrix_free(m);
    if (m_LU != nullptr)
      gsl_matrix_free(m_LU);
    if (perm != nullptr)
      gsl_permutation_free(perm);
  }

  SqMatrix(const SqMatrix &matrix) // copy constructor
      : n(matrix.n), m(gsl_matrix_alloc(matrix.n, matrix.n)) {
    // Have to add copy constructor for deep copy of m ptr
    gsl_matrix_memcpy(m, matrix.m);
    // XXX doesn't copy LU etc!
  }

  SqMatrix &operator=(const SqMatrix &other) // copy assignment
  {
    if (this->n != other.n) {
      std::cerr << "FAIL 35 in SqMatrix. Cant re-assign matrices of different "
                   "dimension: N_lhs = "
                << this->n << ", N_rhs = " << other.n << "\n\n";
      std::abort();
    }
    gsl_matrix_memcpy(m, other.m);
    return *this;
    // XXX doesn't copy LU etc!
  }

  //** Public functions

  void make_diag(double value = 1.0) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        (*this)[i][j] = (i == j) ? value : 0.0;
      }
    }
  }

  void scale(double value) {
    const int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] *= value;
  }

  void clip_low(double value) {
    const int n2 = n * n;
    for (int i = 0; i < n2; i++) {
      if (std::abs(m->data[i]) < value)
        m->data[i] = 0.0;
    }
  }
  void clip_high(double value) {
    const int n2 = n * n;
    for (int i = 0; i < n2; i++) {
      if (std::abs(m->data[i]) > value) {
        auto s = m->data[i] > 0 ? 1 : -1;
        m->data[i] = s * value;
      }
    }
  }

  void make_symmetric() {
    *this += this->transpose();
    this->scale(0.5);
  }

  double check_symmetric() {
    double worst = 0.0;
    auto AmATr = *this - this->transpose();
    for (int i = 0; i < n * n; i++) {
      auto val = std::abs(AmATr.m->data[i]);
      worst = (val > worst) ? val : worst;
    }
    return worst;
  }

  void print() {
    // mostly for testing
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        // std::cout << (*this)[i][j] << " ";
        printf("%8.1e ", (*this)[i][j]);
      }
      std::cout << "\n";
    }
  }

  [[nodiscard]] SqMatrix transpose() {
    SqMatrix mTr(this->n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        mTr[i][j] = (*this)[j][i];
      }
    }
    return mTr;
  }

  void LU_decompose() {
    if (m_LU != nullptr)
      return; //? what if matrix changed!?
    m_LU = gsl_matrix_alloc(n, n);
    perm = gsl_permutation_alloc(n);
    gsl_matrix_memcpy(m_LU, m);
    gsl_linalg_LU_decomp(m_LU, perm, &s_LU);
  }

  double determinant() {
    LU_decompose();
    return gsl_linalg_LU_det(m_LU, s_LU);
  }

  void invert() {
    // note: this is destuctive: matrix will be inverted
    // Usuing this method, hard not to be: LU decomp changes Matrix
    //(So, would require copy to avoid this)
    LU_decompose();
    gsl_linalg_LU_invert(m_LU, perm, m);
  }

  [[nodiscard]] SqMatrix inverse() const {
    auto inverse = *this;
    inverse.invert();
    return inverse;
  }

  double *operator[](int i) const { return &(m->data[i * n]); }

  friend SqMatrix operator*(const SqMatrix &lhs, const SqMatrix &rhs) {
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

  SqMatrix &operator+=(const SqMatrix rhs) {
    int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] += rhs.m->data[i];
    return *this;
  }
  friend SqMatrix operator+(SqMatrix lhs, const SqMatrix &rhs) {
    lhs += rhs;
    return lhs;
  }
  SqMatrix &operator-=(const SqMatrix rhs) {
    int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] -= rhs.m->data[i];
    return *this;
  }
  friend SqMatrix operator-(SqMatrix lhs, const SqMatrix &rhs) {
    lhs -= rhs;
    return lhs;
  }
  SqMatrix &operator*=(const double x) {
    scale(x);
    return *this;
  }
  friend SqMatrix operator*(const double x, SqMatrix rhs) {
    rhs *= x;
    return rhs;
  }
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

inline void test() {
  double data[] = {1.0,     1 / 2.0, 1 / 3.0, 1 / 4.0, 1 / 2.0, 1 / 3.0,
                   1 / 4.0, 1 / 5.0, 1 / 3.0, 1 / 4.0, 1 / 5.0, 1 / 6.0,
                   1 / 4.0, 1 / 5.0, 1 / 6.0, 1 / 7.0};

  gsl_matrix_view m = gsl_matrix_view_array(data, 4, 4);

  gsl_vector *eval = gsl_vector_alloc(4);
  gsl_matrix *evec = gsl_matrix_alloc(4, 4);

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(4);

  gsl_eigen_symmv(&m.matrix, eval, evec, w);

  gsl_eigen_symmv_free(w);

  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  {
    int i;

    for (i = 0; i < 4; i++) {
      double eval_i = gsl_vector_get(eval, i);
      gsl_vector_view evec_i = gsl_matrix_column(evec, i);

      printf("eigenvalue = %g\n", eval_i);
      printf("eigenvector = \n");
      gsl_vector_fprintf(stdout, &evec_i.vector, "%g");
    }
  }

  gsl_vector_free(eval);
  gsl_matrix_free(evec);
}

inline void test2(const SqMatrix &B) {
  // No. Use: Real Generalized Nonsymmetric Eigensystems

  // double data[] = {-1.0, 1.0, -1.0, 1.0, -8.0, 4.0,  -2.0, 1.0,
  // 27.0, 9.0, 3.0,  1.0, 64.0, 16.0, 4.0,  1.0};
  //   eigenvalue = -6.41391 + 0i
  // eigenvalue = 5.54555 + 3.08545i
  // eigenvalue = 5.54555 + -3.08545i
  // eigenvalue = 2.3228 + 0i

  // gsl_matrix_view m = gsl_matrix_view_array(B.m, B.n, B.n);

  gsl_vector_complex *eval = gsl_vector_complex_alloc(B.n);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc(B.n, B.n);

  gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(B.n);

  gsl_eigen_nonsymmv(B.m, eval, evec, w);

  gsl_eigen_nonsymmv_free(w);

  gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  {

    for (int i = 0; i < B.n; i++) {
      gsl_complex eval_i = gsl_vector_complex_get(eval, i);
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);

      printf("eigenvalue %i = %g + %gi\n", i, GSL_REAL(eval_i),
             GSL_IMAG(eval_i));
      // printf("eigenvector = \n");
      // for (j = 0; j < 4; ++j) {
      //   gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
      //   printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
      // }
    }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}

inline void test3(const SqMatrix &B, const SqMatrix &S) {
  // No. Use: Real Generalized Nonsymmetric Eigensystems

  // double data[] = {-1.0, 1.0, -1.0, 1.0, -8.0, 4.0,  -2.0, 1.0,
  // 27.0, 9.0, 3.0,  1.0, 64.0, 16.0, 4.0,  1.0};
  //   eigenvalue = -6.41391 + 0i
  // eigenvalue = 5.54555 + 3.08545i
  // eigenvalue = 5.54555 + -3.08545i
  // eigenvalue = 2.3228 + 0i

  // gsl_matrix_view m = gsl_matrix_view_array(B.m, B.n, B.n);
  const auto n = B.n;

  gsl_eigen_gen_workspace *work = gsl_eigen_gen_alloc(n);

  // I think this is not needed:
  gsl_eigen_gen_params(0, 0, 0, work);

  gsl_vector_complex *alpha = gsl_vector_complex_alloc(B.n);
  gsl_vector *beta = gsl_vector_alloc(B.n);
  gsl_eigen_gen(B.m, S.m, alpha, beta, work);

  std::vector<double> evals;
  for (int i = 0; i < B.n; i++) {
    gsl_complex eval_ai = gsl_vector_complex_get(alpha, i);
    double eval_bi = gsl_vector_get(beta, i);
    // gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec, i);

    auto evr = GSL_REAL(eval_ai) / eval_bi;
    auto evi = GSL_IMAG(eval_ai) / eval_bi;
    evals.push_back(evr);
    // printf("eigenvalue %i = %g + %gi\n", i, evr, evi);
  }
  std::sort(evals.begin(), evals.end());
  for (const auto &ev : evals) {
    std::cout << ev << "\n";
  }

  gsl_eigen_gen_free(work);
  gsl_vector_complex_free(alpha);
  gsl_vector_free(beta);
}

} // namespace LinAlg
