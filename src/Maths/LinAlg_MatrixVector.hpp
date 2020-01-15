#pragma once
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <tuple>
#include <utility>
#include <vector>
// #include <iostream>
// #include <algorithm>
// #include <array>

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
  // template <typename T>
  SqMatrix(const std::initializer_list<double> &l);
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

  void make_symmetric(); // change to "symmetrise"
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
  Vector(const int in_n);
  template <typename T> Vector(const std::initializer_list<T> &l);
  Vector &operator=(const Vector &other); // copy assignment

  Vector(const Vector &other); // copy constructor;

  ~Vector();

  // Vector(const Vector &vector); // copy constructor

  void clip_low(const double value);
  void clip_high(const double value);
  void print();

  // add inner product, outerproduct

  double &operator[](int i) const;
  Vector &operator+=(const Vector rhs);
  friend Vector operator+(Vector lhs, const Vector &rhs);
  Vector &operator-=(const Vector rhs);
  friend Vector operator-(Vector lhs, const Vector &rhs);
  Vector &operator*=(const double x);
  friend Vector operator*(const double x, Vector rhs);
  friend Vector operator*(const SqMatrix &Aij, const Vector &bj);

  friend double inner_produce(const Vector &a, const Vector &b);
  friend double operator*(const Vector &a, const Vector &b);
  friend SqMatrix outer_produce(const Vector &a, const Vector &b);
};

//******************************************************************************
//******************************************************************************

Vector solve_Axeqb(SqMatrix &Am, const Vector &b);

//------------------------------------------------------------------------------
// Eigensystems: Note: A and B are destroyed! Don't use afterwards
// (Can't avoid this without needless copy)
[[nodiscard]] std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, bool sort = true);

[[nodiscard]] std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true);

[[nodiscard]] std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, bool sort = true);

[[nodiscard]] std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true);

} // namespace LinAlg
