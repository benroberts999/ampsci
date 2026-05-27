#pragma once
#include "Matrix.hpp"
#include "Vector.hpp"

namespace LinAlg {

/*!
  @brief Solves the linear system Ax = b for x, given square matrix A and
  vector b.
  @details
  Uses LU decomposition (GSL). A is taken by value and overwritten during
  decomposition.

  @tparam T Element type: `double` or `std::complex<double>` only.
  @param Am Square matrix A (copied; overwritten internally).
  @param b  Right-hand side vector; must satisfy `b.size() == A.rows()`.
  @returns  Solution vector x.
*/
template <typename T>
Vector<T> solve_Axeqb(Matrix<T> Am, const Vector<T> &b);

/*!
  @brief Solves Av = ev for all eigenvalues and eigenvectors of a real
  symmetric or complex Hermitian matrix A.
  @details
  Uses LAPACK `dsyev` (real) or `zheev` (complex Hermitian).
  Eigenvalues are returned in ascending order.
  A is taken by value and overwritten during the computation.

  @tparam T Element type: `double` (symmetric) or `std::complex<double>`
  (Hermitian).
  @param A  Square symmetric/Hermitian matrix (taken by value; overwritten).
  @returns  Pair `{e, V}`, where:
            - `e(i)` is the i-th eigenvalue (always real, ascending order),
            - `V(i,j)` is the j-th element of the i-th eigenvector.
*/
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A);

/*!
  @brief Solves Av = ev for the first `number` eigenvalues and eigenvectors
  of a real symmetric matrix A.
  @details
  Uses LAPACK `dsyevx`. Only available for `T = double`.

  @tparam T Element type: `double` only.
  @param A      Square real symmetric matrix (taken by value; overwritten).
  @param number Number of lowest eigenvalues/eigenvectors to compute.
                Must satisfy `number < A.rows()`.
  @returns  Pair `{e, V}`, where only the first `number` entries of `e` and
            rows of `V` are valid.
*/
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, int number);

/*!
  @brief Solves Av = ev for all eigenvalues below a given threshold of a real
  symmetric matrix A.
  @details
  Uses LAPACK `dsyevr`. Only available for `T = double`.

  @tparam T Element type: `double` only.
  @param A         Square real symmetric matrix (taken by value; overwritten).
  @param all_below Upper bound: only eigenvalues ≤ `all_below` are returned.
  @returns  Tuple `{N, e, V}`, where:
            - `N` is the number of eigenvalues found,
            - `e(i)` is the i-th eigenvalue (ascending),
            - `V(i,j)` is the j-th element of the i-th eigenvector.
            Only the first `N` entries of `e` and rows of `V` are valid.
*/
template <typename T>
std::tuple<int, Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A,
                                                            double all_below);

/*!
  @brief Solves the generalised eigensystem Av = eBv for a real symmetric or
  complex Hermitian matrix pair A, B.
  @details
  Uses LAPACK `dsygv` (real) or `zhegv` (complex Hermitian).
  B must be positive definite. Eigenvalues are returned in ascending order.
  Both A and B are taken by value and overwritten during computation.

  @tparam T Element type: `double` or `std::complex<double>`.
  @param A  Square symmetric/Hermitian matrix (taken by value; overwritten).
  @param B  Square symmetric/Hermitian positive-definite matrix (same size as
            A; taken by value; overwritten).
  @returns  Pair `{e, V}`, where:
            - `e(i)` is the i-th eigenvalue (always real, ascending order),
            - `V(i,j)` is the j-th element of the i-th eigenvector.
*/
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, Matrix<T> B);

/*!
  @brief Solves Av = ev for all eigenvalues and right eigenvectors of a
  general (non-symmetric) real matrix A.
  @details
  Uses GSL `gsl_eigen_nonsymmv`. A must be real (`T = double`); eigenvalues
  and eigenvectors are always complex.

  @tparam T Element type: `double` only.
  @param A    Square real matrix (taken by value; overwritten).
  @param sort If `true`, eigenvalues (and corresponding eigenvectors) are
              sorted by ascending absolute value.
  @returns  Pair `{e, V}`, where:
            - `e(i)` is the i-th complex eigenvalue,
            - `V(i,j)` is the j-th element of the i-th complex eigenvector.
*/
template <typename T>
std::pair<Vector<std::complex<double>>, Matrix<std::complex<double>>>
genEigensystem(Matrix<T> A, bool sort);

//==============================================================================
//==============================================================================
} // namespace LinAlg

#include "Solvers.ipp"
