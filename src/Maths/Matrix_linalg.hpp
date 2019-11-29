#pragma once
#include <algorithm>
#include <array>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>

/*
Eigen probably better option!
*/

namespace Matrix {

class SqMatrix {
  // Suddenly gets very slow around n=170.. why? When other method doesn't
private:
  gsl_matrix *m;

public:
  const int n;

public:
  SqMatrix(int in_n) : m(gsl_matrix_alloc(in_n, in_n)), n(in_n) {}
  ~SqMatrix() { gsl_matrix_free(m); }
  SqMatrix(const SqMatrix &matrix) // copy constructor
      : m(gsl_matrix_alloc(matrix.n, matrix.n)), n(matrix.n) {
    // Have to add copy constructor for deep copy of m ptr
    int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] = matrix.m->data[i];
  }
  SqMatrix &operator=(const SqMatrix &other) // copy assignment
  {
    if (this->n != other.n) {
      std::cerr << "FAIL 35 in SqMatrix. Cant re-assign matrices of different "
                   "dimension: N_lhs = "
                << this->n << ", N_rhs = " << other.n << "\n\n";
      std::abort();
    }
    int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] = other.m->data[i];
    return *this;
  }

  //** Public functions
  void clear(double value = 0.0) {
    const int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] = value;
  }
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

  void print() {
    // mostly for testing
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        std::cout << (*this)[i][j] << " ";
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

  // SqMatrix make_copy() const {
  //   SqMatrix copiedM(n);
  //   int n2 = n * n;
  //   for (int i = 0; i < n2; i++)
  //     copiedM.m->data[i] = m->data[i];
  //   return copiedM;
  // }

  void invert() {
    // note: this is destuctive: matrix will be inverted
    // Usuing this method, hard not to be: LU decomp changes Matrix
    //(So, would require copy to avoid this)
    gsl_matrix *inverse = gsl_matrix_alloc(n, n);
    gsl_permutation *perm = gsl_permutation_alloc(n);
    int s;
    gsl_linalg_LU_decomp(m, perm, &s);
    gsl_linalg_LU_invert(m, perm, inverse);
    int n2 = n * n;
    for (int i = 0; i < n2; i++)
      m->data[i] = inverse->data[i];
    gsl_permutation_free(perm);
    gsl_matrix_free(inverse);
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

/*
auto m = Matrix::SqMatrix(2);
for (int i = 0; i < 2; ++i)
  for (int j = 0; j < 2; ++j)
    m[i][j] = i + j + 6.001;
for (int i = 0; i < 2; ++i)
  for (int j = 0; j < 2; ++j)
    std::cout << m[i][j] << "\n";

auto m2 = m.make_copy();

auto minv = m.inverse();
auto m3 = m * minv;
m3.clip_low(1.e-6);
auto m4 = 6 * m3;
for (int i = 0; i < 2; ++i)
  for (int j = 0; j < 2; ++j)
    std::cout << m3[i][j] << " " << m4[i][j] << "\n";
*/

//
//
//
//******************************************************************************
// template <typename T>
// std::vector<std::vector<T>> invert(const std::vector<std::vector<T>> &M) {
//
//   // size of matrix:
//   auto n = M.size();
//   if (M[0].size() != n)
//     std::cerr << "\nCant invert non-square matrix, silly.\n";
//
//   // Define all the used matrices (for GSL)
//   gsl_matrix *m = gsl_matrix_alloc(n, n);
//   gsl_matrix *inverse = gsl_matrix_alloc(n, n);
//   gsl_permutation *perm = gsl_permutation_alloc(n);
//
//   // fill matrix:
//   for (std::size_t i = 0; i < n; i++) {
//     for (std::size_t j = 0; j < n; j++)
//       gsl_matrix_set(m, i, j, M[i][j]);
//   }
//
//   // peform LU decomposition (using GSL)
//   // and inversion (if non-singular)
//   int s;
//   gsl_linalg_LU_decomp(m, perm, &s);
//   double det = gsl_linalg_LU_det(m, s);
//
//   //"output" inverted matrix
//   std::vector<std::vector<T>> W(n, std::vector<T>(n));
//
//   // Fill the output matrix:
//   if (det != 0) {
//     gsl_linalg_LU_invert(m, perm, inverse);
//     for (std::size_t i = 0; i < n; i++) {
//       for (std::size_t j = 0; j < n; j++)
//         W[i][j] = gsl_matrix_get(inverse, i, j);
//     }
//   }
//
//   // clear memory
//   gsl_permutation_free(perm);
//   gsl_matrix_free(m);
//   gsl_matrix_free(inverse);
//
//   return W;
// }
//
// //******************************************************************************
// template <typename T, std::size_t n>
// std::array<std::array<T, n>, n>
// invert(const std::array<std::array<T, n>, n> &M) {
//
//   // Define all the used matrices (for GSL)
//   gsl_matrix *m = gsl_matrix_alloc(n, n);
//   gsl_matrix *inverse = gsl_matrix_alloc(n, n);
//   gsl_permutation *perm = gsl_permutation_alloc(n);
//
//   // fill matrix:
//   for (std::size_t i = 0; i < n; i++) {
//     for (std::size_t j = 0; j < n; j++)
//       gsl_matrix_set(m, i, j, M[i][j]);
//   }
//
//   // peform LU decomposition (using GSL)
//   // and inversion (if non-singular)
//   int s;
//   gsl_linalg_LU_decomp(m, perm, &s);
//   double det = gsl_linalg_LU_det(m, s);
//
//   // Fill the output matrix:
//   //"output" inverted matrix
//   std::array<std::array<T, n>, n> W;
//   if (det != 0) {
//     gsl_linalg_LU_invert(m, perm, inverse);
//     for (std::size_t i = 0; i < n; i++) {
//       for (std::size_t j = 0; j < n; j++)
//         W[i][j] = gsl_matrix_get(inverse, i, j);
//     }
//   }
//
//   // clear memory
//   gsl_permutation_free(perm);
//   gsl_matrix_free(m);
//   gsl_matrix_free(inverse);
//
//   return W;
// }

} // namespace Matrix

/*
  //CODE TO TEST:
  ChronoTimer sw; // start the overall timer

  int num = 1000;
  const int dim = 50;
  // std::vector<std::array<std::array<double, dim>, dim>> M(num);
  std::array<std::array<double, dim>, dim> m;
  Matrix::SqMatrix m2(dim);

  // std::vector<Matrix::SqMatrix<double>> M3;
  // for (int i = 0; i < num; i++) {
  //   Matrix::SqMatrix<double> M_tmp(dim);
  //   M3.push_back(M_tmp);
  // }

  // for (auto &m : M) {
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      double x = Rand(0.3);
      m[i][j] = (i == j) ? double(i + 1) + x : x;
      m2[i][j] = m[i][j];
      // M2[i][j] = m[i][j];
      // std::cout << m[i][j] << " ";
    }
    // std::cout << "\n";
  }
  // std::cout << "\n";
  // }

  sw.start();
  for (int i = 0; i < num; i++)
    m2.invert();
  std::cout << "\n Invert2: " << sw.lap_reading_str() << "\n";

  sw.start();
  for (int i = 0; i < num; i++)
    m = Matrix::invert(m);
  std::cout << "\n Invert1: " << sw.lap_reading_str() << "\n";

  sw.start();
  for (int i = 0; i < num; i++)
    m = Matrix::invert(m);
  std::cout << "\n Invert1: " << sw.lap_reading_str() << "\n";

  sw.start();
  for (int i = 0; i < num; i++)
    m2.invert();
  std::cout << "\n Invert2: " << sw.lap_reading_str() << "\n";

  // std::cin.get();

  // for (auto &m : M) {
  //   for (int i = 0; i < dim; i++) {
  //     for (int j = 0; j < dim; j++) {
  //       std::cout << m[i][j] << " ";
  //     }
  //     std::cout << "\n";
  //   }
  // }
  // std::cout << "\n";
  //
  // M2.invert();
  // for (int i = 0; i < dim; i++) {
  //   for (int j = 0; j < dim; j++) {
  //     std::cout << M2[i][j] << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "\n";

  // for (int i = 0; i < dim; i++) {
  //   for (int j = 0; j < dim; j++) {
  //     double x = Rand(0.001);
  //     // M2.tempGetSet(i, j) = (i == j) ? double(i + 1) + x : x;
  //     M2[i][j] = (i == j) ? double(i + 1) + x : x;
  //     std::cout << M2[i][j] << " ";
  //     // std::cout << M2.tempGetSet(i, j) << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "\n";
  //
  // for (int i = 0; i < dim; i++) {
  //   for (int j = 0; j < dim; j++) {
  //     // double x = Rand(0.3);
  //     // M2[i][j] = (i == j) ? double(i + 1) + x : x;
  //     // std::cout << M2[i][j] << " ";
  //     std::cout << M2.tempGetSet(i, j) << " ";
  //   }
  //   std::cout << "\n";
  // }
  // std::cout << "\n";

  return 1;

*/
