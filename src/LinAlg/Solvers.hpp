#pragma once
#include "Matrix.hpp"
#include "Vector.hpp"

namespace LinAlg {

//! Solves matrix equation Ax=b for x, for known square matrix A and vector b.
template <typename T> Vector<T> solve_Axeqb(Matrix<T> Am, const Vector<T> &b);

//! Solves Av = ev for eigenvalues e and eigenvectors v of symmetric/Hermetian
//! matrix A. Returns [e,v], where v(i,j) is the jth element of the ith
//! eigenvector corresponding to ith eigenvalue, e(i). e is always real.
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, bool sort);

//! Solves Av = eBv for eigenvalues e and eigenvectors v of symmetric/Hermetian
//! matrix pair A,B. Returns [e,v], where v(i,j) is the jth element of the ith
//! eigenvector corresponding to ith eigenvalue, e(i). e is always real.
template <typename T>
std::pair<Vector<double>, Matrix<T>> symmhEigensystem(Matrix<T> A, Matrix<T> B,
                                                      bool sort);

//! Solves Av = ev for eigenvalues e and eigenvectors v of non-symmetric real
//! matrix A. Returns [e,v], where v(i,j) is the jth element of the ith
//! eigenvector corresponding to ith eigenvalue, e(i). A must be real, while e
//! and v are complex.
template <typename T>
std::pair<Vector<std::complex<double>>, Matrix<std::complex<double>>>
genEigensystem(Matrix<T> A, bool sort);

//******************************************************************************
//******************************************************************************
} // namespace LinAlg

#include "Solvers.ipp"
