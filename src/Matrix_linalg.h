#pragma once
#include <array>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>

/*

Note: this is a pretty ineficient way of doing things..
Involves lots of copies into/out of gsl

There's probably a better way to template as well.

Eigen probably better option!

*/

namespace Matrix {

// template <typename T>
class SqMatrix {
  /*
    Suddenly gets very slow around n=170.. why? When other method doesn't
  */
private:
  gsl_matrix *m;

public:
  const int n;

public:
  SqMatrix(int in_n) : m(gsl_matrix_alloc(in_n, in_n)), n(in_n) {}
  ~SqMatrix() { gsl_matrix_free(m); }

  double *operator[](int i) { return &(m->data[i * n]); }

  void invert() {
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
};

//******************************************************************************
template <typename T>
std::vector<std::vector<T>> invert(const std::vector<std::vector<T>> &M) {

  // size of matrix:
  auto n = M.size();
  if (M[0].size() != n)
    std::cerr << "\nCant invert non-square matrix, silly.\n";

  // Define all the used matrices (for GSL)
  gsl_matrix *m = gsl_matrix_alloc(n, n);
  gsl_matrix *inverse = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);

  // fill matrix:
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++)
      gsl_matrix_set(m, i, j, M[i][j]);
  }

  // peform LU decomposition (using GSL)
  // and inversion (if non-singular)
  int s;
  gsl_linalg_LU_decomp(m, perm, &s);
  double det = gsl_linalg_LU_det(m, s);

  //"output" inverted matrix
  std::vector<std::vector<T>> W(n, std::vector<T>(n));

  // Fill the output matrix:
  if (det != 0) {
    gsl_linalg_LU_invert(m, perm, inverse);
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t j = 0; j < n; j++)
        W[i][j] = gsl_matrix_get(inverse, i, j);
    }
  }

  // clear memory
  gsl_permutation_free(perm);
  gsl_matrix_free(m);
  gsl_matrix_free(inverse);

  return W;
}

//******************************************************************************
template <typename T, std::size_t n>
std::array<std::array<T, n>, n>
invert(const std::array<std::array<T, n>, n> &M) {

  // Define all the used matrices (for GSL)
  gsl_matrix *m = gsl_matrix_alloc(n, n);
  gsl_matrix *inverse = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);

  // fill matrix:
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++)
      gsl_matrix_set(m, i, j, M[i][j]);
  }

  // peform LU decomposition (using GSL)
  // and inversion (if non-singular)
  int s;
  gsl_linalg_LU_decomp(m, perm, &s);
  double det = gsl_linalg_LU_det(m, s);

  // Fill the output matrix:
  //"output" inverted matrix
  std::array<std::array<T, n>, n> W;
  if (det != 0) {
    gsl_linalg_LU_invert(m, perm, inverse);
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t j = 0; j < n; j++)
        W[i][j] = gsl_matrix_get(inverse, i, j);
    }
  }

  // clear memory
  gsl_permutation_free(perm);
  gsl_matrix_free(m);
  gsl_matrix_free(inverse);

  return W;
}

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
