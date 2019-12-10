#include "LinAlg_MatrixVector.hpp"
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

namespace LinAlg {

//******************************************************************************
//******************************************************************************
// class SqMatrix:

SqMatrix::SqMatrix(int in_n) : n(in_n), m(gsl_matrix_alloc(in_n, in_n)) {}

template <typename T>
SqMatrix::SqMatrix(const std::initializer_list<T> &l)
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

//******************************************************************************
//******************************************************************************
// Solve LinAlg equations:

} // namespace LinAlg
