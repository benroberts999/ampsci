#pragma once
#include <array>
#include <gsl/gsl_linalg.h>
#include <vector>

/*

Note: this is a pretty ineficient way of doing things..
Involves lots of copies into/out of gsl

There's probably a better way to template as well.

Eigen probably better option!

*/

namespace Matrix {

//******************************************************************************
template <typename T>
std::vector<std::vector<T>> invert(const std::vector<std::vector<T>> &M) {

  // size of matrix:
  size_t n = M.size();
  if (M[0].size() != n)
    std::cerr << "\nCant invert non-square matrix, silly.\n";
  // return 0; // matrix not square!

  //"output" inverted matrix
  std::vector<std::vector<T>> W(n, std::vector<T>(n));
  // W.reserve(n);

  // Define all the used matrices (for GSL)
  gsl_matrix *m = gsl_matrix_alloc(n, n);
  gsl_matrix *inverse = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);

  // fill matrix:
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++)
      gsl_matrix_set(m, i, j, M[i][j]);
  }

  // peform LU decomposition (using GSL)
  // and inversion (if non-singular)
  int s;
  gsl_linalg_LU_decomp(m, perm, &s);
  double det = gsl_linalg_LU_det(m, s);
  if (det != 0)
    gsl_linalg_LU_invert(m, perm, inverse);
  if (det == 0)
    return W; // ERROR!

  // Fill the output matrix:
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++)
      W[i][j] = gsl_matrix_get(inverse, i, j);
  }

  // clear memory
  gsl_permutation_free(perm);
  gsl_matrix_free(m);
  gsl_matrix_free(inverse);

  return W;
}

//******************************************************************************
template <typename T, size_t n>
std::array<std::array<T, n>, n>
invert(const std::array<std::array<T, n>, n> &M) {

  //"output" inverted matrix
  std::array<std::array<T, n>, n> W;

  // Define all the used matrices (for GSL)
  gsl_matrix *m = gsl_matrix_alloc(n, n);
  gsl_matrix *inverse = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);

  // fill matrix:
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++)
      gsl_matrix_set(m, i, j, M[i][j]);
  }

  // peform LU decomposition (using GSL)
  // and inversion (if non-singular)
  int s;
  gsl_linalg_LU_decomp(m, perm, &s);
  double det = gsl_linalg_LU_det(m, s);
  if (det != 0)
    gsl_linalg_LU_invert(m, perm, inverse);
  if (det == 0)
    return W; // ERROR!

  // Fill the output matrix:
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++)
      W[i][j] = gsl_matrix_get(inverse, i, j);
  }

  // clear memory
  gsl_permutation_free(perm);
  gsl_matrix_free(m);
  gsl_matrix_free(inverse);

  return W;
}

} // namespace Matrix

//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
//******************************************************************************
//******************************************************************************
//******************************************************************************
// namespace MAlg {
//
// //***************************************************************************
// double calcDeterminant(const std::vector<std::vector<double>> &inmat)
//
// // 170622.
// // Calculates the determinant of any real square matrix of dimension n.
// // Uses the GNU 'GSL' libraries:
// //
// https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
// // https://www.gnu.org/software/gsl/manual/html_node/LU-Decomposition.html
// // Requires #include <gsl/gsl_linalg.h>
// //
// // INPUT:
// //   inmat     :: double matix of dimension n*n [from std::vector]
// //   n         :: integer, dimension of matrices
// //
// // ### Change Log ###
// // 170702- Uses std::vector input, avoid variable arrays!
//
// {
//
//   // size of array:
//   int n = (int)inmat.size();
//   if (inmat[0].size() != (size_t)n)
//     return 0; // matrix not square!
//
//   // Define all the used matrices (for GSL)
//   gsl_matrix *m = gsl_matrix_alloc(n, n);
//   gsl_permutation *perm = gsl_permutation_alloc(n);
//   // fill matrix:
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       // gsl_matrix_set(m,i,j,inmat[i*n+j]);
//       gsl_matrix_set(m, i, j, inmat[i][j]);
//     }
//   }
//   // peform LU decomposition (using GSL)
//   int s;
//   gsl_linalg_LU_decomp(m, perm, &s);
//   double det = gsl_linalg_LU_det(m, s); // XXX ok as double?
//   // clear memory
//   gsl_permutation_free(perm);
//   gsl_matrix_free(m);
//   return det;
// }
//
// //---- Overloaded:
// ------------------------------------------------------------- double
// calcDeterminant(const std::vector<std::vector<float>> &inmat) {
//
//   int n = (int)inmat.size();
//   int m = (int)inmat[0].size();
//   std::vector<std::vector<double>> dbl_inmat(n, std::vector<double>(m));
//
//   for (int i = 0; i < n; i++) {
//     for (int j = 0; j < m; j++) {
//       dbl_inmat[i][j] = (double)inmat[i][j];
//     }
//   }
//
//   return calcDeterminant(dbl_inmat);
// }
//
// //***************************************************************************
// int linsolve(const std::vector<std::vector<double>> &inmat,
//              const std::vector<double> &invec, std::vector<double> &outvec)
//
// // 170321.
// // Solves the linear matrix equation A.x=b for x
// // Where:
// //   A=inmat    is an n*n square matrix
// //   x=outvec   is an n dim vector (the answer/output!)
// //   b=invec    is an n dim vector (the input)
// //
// // Uses the GNU 'GSL' libraries:
// //
// https://www.gnu.org/software/gsl/manual/html_node/Linear-Algebra-Examples.html
// // https://www.gnu.org/software/gsl/manual/html_node/LU-Decomposition.html
// // Requires #include <gsl/gsl_linalg.h>
// //
// // INPUT:
// //   inmat  :: double 'flat' matix of dimension n*n [call as "(double
// *)inmat"]
// //   invec  :: double vector of dimension n
// //   n      :: integer, dimension of matrices
// // OUTPUT:
// //   outvec :: solution. output n dimensional vector
// //
// //
//
// {
//
//   int n = (int)inmat.size();
//
//   if (inmat[0].size() != (size_t)n)
//     return 1; // matrix not square!
//   if (invec.size() != (size_t)n)
//     return 1; // invec incorrect dimension!
//
//   // ensure outvec has correct dimension:
//   outvec.resize(n);
//
//   int iRet = 0;
//   // Define all the used matrices/vectors:
//   gsl_matrix *A = gsl_matrix_alloc(n, n);
//   gsl_vector *b = gsl_vector_alloc(n);
//   gsl_vector *x = gsl_vector_alloc(n);
//   gsl_permutation *p = gsl_permutation_alloc(n);
//   // fill matrix/vector:
//   for (int i = 0; i < n; i++) {
//     gsl_vector_set(b, i, invec[i]);
//     for (int j = 0; j < n; j++) {
//       // gsl_matrix_set(A,i,j,inmat[i*n+j]);
//       gsl_matrix_set(A, i, j, inmat[i][j]);
//     }
//   }
//   // peform LU decomposition (using GSL)
//   // and solve linear equation A.x=b for x (if non-singular)
//   int s;
//   gsl_linalg_LU_decomp(A, p, &s);
//   double det = gsl_linalg_LU_det(A, s);
//   if (det != 0)
//     gsl_linalg_LU_solve(A, p, b, x);
//   if (det == 0)
//     iRet = 1;
//   // Fill the output vector:
//   for (int i = 0; i < n; i++) {
//     outvec[i] = gsl_vector_get(x, i);
//   }
//   // clear memory
//   gsl_permutation_free(p);
//   gsl_matrix_free(A);
//   gsl_vector_free(x);
//   gsl_vector_free(b);
//   return iRet;
// }
//
// } // namespace MAlg
