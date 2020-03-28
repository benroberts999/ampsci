#pragma once
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <tuple>
#include <utility>
#include <vector>

//! Defines SqMatrix, Vector, and linear-algebra solvers (incl Eigensystems)
namespace LinAlg {

//******************************************************************************
//! Basic Square matrix class of constant construct-time size
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
  //! Constructs a diagonal matrix with values value
  void make_diag(double value = 1.0);
  //! Sets all elements to zero
  void zero();
  //! Scale all elements by constant value
  void scale(double value);
  //! All values with |Mij|<value are set to zero
  void clip_low(double value);
  //! All values with |Mij|>value are set to +/-value
  void clip_high(double value);

  //! Forces Matrix to be symmetric by: Mij -> (Mij + Mji)/2
  void make_symmetric(); // change to "symmetrise"?
  //! Returns largest elemt of matrix [Mij - Mji]; zero if symmetric
  double check_symmetric();
  //! Prints Matrix to screen
  void print();

  //! Returns the transpose of matrix: not destructive
  [[nodiscard]] SqMatrix transpose() const;
  double determinant(); // changes m_LU
  //! Inverts the matrix: nb: destructive!
  void invert();
  //! Returns the inverce of matrix: not destructive
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
//! Basic vector class of constant construct-time size
class Vector {
public:
  const int n;

  // private:
  gsl_vector *vec;

public:
  Vector(const int in_n);
  template <typename T> Vector(const std::initializer_list<T> &l);
  Vector &operator=(const Vector &other); // copy assignment
  Vector(const Vector &other);            // copy constructor;
  ~Vector();

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

//! Solves Matrix equationL A*x = b for x
Vector solve_Axeqb(SqMatrix &Am, const Vector &b);

//------------------------------------------------------------------------------
//! @brief Solves Av = ev for eigenvalues e and eigenvectors v
//! for Real Symmetric Matrices
//! @details Eigensystems: Note: A and B are destroyed! Don't use afterwards
//! (Can't avoid this without needless copy). Optionally sorts e and v by e
[[nodiscard]] std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, bool sort = true);

//! @brief Solves Av = eBv for eigenvalues e and eigenvectors v for Real
//! Generalized Symmetric-Definite Eigensystems
[[nodiscard]] std::pair<Vector, SqMatrix>
realSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true);

//! @briefSolves for Av = ev for Real Nonsymmetric Matrices.
//! @details e and v are complex; returned as {real, imag} (seperate
//! vectors/Matrix)
[[nodiscard]] std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, bool sort = true);

//! @brief Solves Av = eBv for eigenvalues e and eigenvectors v for Real
//! Generalized Non-Symmetric-Definite Eigensystems.
[[nodiscard]] std::tuple<Vector, Vector, SqMatrix, SqMatrix>
realNonSymmetricEigensystem(SqMatrix &A, SqMatrix &B, bool sort = true);

} // namespace LinAlg
