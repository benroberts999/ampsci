#pragma once
#include "Maths/LinAlg_MatrixVector.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <random>
#include <string>

namespace UnitTest {

namespace helper {

// returns matrix element that has largest absolute value
inline double max_matrixel(const LinAlg::SqMatrix &m);
// returns vector element that has largest absolute value
inline double max_vecel(const LinAlg::Vector &b);
// returns complex matrix element that has largest absolute value
inline std::pair<double, double> max_matrixel(const LinAlg::ComplexSqMatrix &m);

} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests for Linear algebra classes/functions (Matrix, eigensystems)
bool LinAlg(std::ostream &obuff) {
  bool pass = true;

  { // Determinant: compare against known value
    // This 6x6 matrix has det = 2619387/4096
    LinAlg::SqMatrix m6x6(6);
    for (std::size_t i = 0; i < m6x6.n; ++i) {
      for (std::size_t j = 0; j < m6x6.n; ++j) {
        m6x6[i][j] =
            (i == j) ? double(i + 1) : (0.25 * double(i + 1)) / double(j + 1);
      }
    }
    const auto expected = 2619387.0 / 4096.0;
    pass &= qip::check_value(&obuff, "determinant",
                             (m6x6.determinant() - expected), 0.0, 1.0e-16);
  }

  { // Inverse:
    // also test: plusIdent, make_identity, mult* and -
    LinAlg::SqMatrix m(100);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-0.1, 0.1);
    for (std::size_t i = 0; i < m.n; ++i) {
      for (std::size_t j = 0; j < m.n; ++j) {
        m[i][j] = dis(gen);
      }
    }
    m.plusIdent(1.0); // "ensure" invertable
    const auto minv = m.inverse();
    const auto minv_m = minv * m;
    LinAlg::SqMatrix ident(m.n);
    ident.make_identity();
    const auto del = helper::max_matrixel(minv_m - ident);
    pass &=
        qip::check_value(&obuff, "inverse + matrix-mult", del, 0.0, 1.0e-14);
  }

  { // Transpose
    LinAlg::SqMatrix m(75);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-1.0e5, 1.0e5);
    for (std::size_t i = 0; i < m.n; ++i) {
      for (std::size_t j = 0; j < m.n; ++j) {
        m[i][j] = dis(gen);
      }
    }

    const auto mtr = m.transpose();
    auto sym_sum = 0.0;
    for (std::size_t i = 0; i < m.n; ++i) {
      for (std::size_t j = 0; j < m.n; ++j) {
        sym_sum += m[i][j] - mtr[j][i]; // should all be exactly zero
      }
    }
    pass &= qip::check_value(&obuff, "transpose", sym_sum, 0.0, 1.0e-16);
  }

  { // Multiply elements
    LinAlg::SqMatrix m1(50);
    LinAlg::SqMatrix m2(m1.n);
    LinAlg::SqMatrix m3(m1.n); // m3_ij = m1_ij * m2_ij
    LinAlg::SqMatrix m4(m1.n); // m4 = x * m1;
    const auto x = 3.1415926;
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-1.0e5, 1.0e5);
    for (std::size_t i = 0; i < m1.n; ++i) {
      for (std::size_t j = 0; j < m1.n; ++j) {
        const auto v1 = dis(gen);
        const auto v2 = dis(gen);
        m1[i][j] = v1;
        m2[i][j] = v2;
        m3[i][j] = v1 * v2;
        m4[i][j] = x * v1;
      }
    }

    // Expect: LinAlg::mult_elements(m1, m2) == m3
    const auto m5 = LinAlg::SqMatrix::mult_elements(m1, m2) - m3;

    pass &= qip::check_value(&obuff, "Multiply els", helper::max_matrixel(m5),
                             0.0, 1.0e-16);

    // Also test scalar mult: expect m6 = m4
    const auto m6 = x * m1;
    pass &= qip::check_value(&obuff, "Scalar mult",
                             helper::max_matrixel(m6 - m4), 0.0, 1.0e-16);
  }

  //****************************************************************************
  // Real, symmetric Eigensystems

  { // Ax = b
    LinAlg::SqMatrix m(50);
    LinAlg::Vector b(m.n);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (std::size_t i = 0; i < m.n; ++i) {
      b[i] = dis(gen);
      for (std::size_t j = 0; j < m.n; ++j) {
        m[i][j] = dis(gen);
      }
    }
    m.plusIdent(20.0);

    const auto x1 = LinAlg::solve_Axeqb(m, b);
    const auto x2 = m.inverse() * b;

    pass &= qip::check_value(&obuff, "Solve Ax=b", helper::max_vecel(x1 - x2),
                             0.0, 1.0e-14);
    pass &= qip::check_value(&obuff, "Solve Ax=b",
                             helper::max_vecel(m * x1 - b), 0.0, 1.0e-14);
  }

  { // RS Eigensystems Av = ev
    LinAlg::SqMatrix a(75);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (std::size_t i = 0; i < a.n; ++i) {
      for (std::size_t j = 0; j < a.n; ++j) {
        a[i][j] = dis(gen);
        if (i == j)
          a[i][j] += double(10 + 2 * i);
      }
    }
    a.enforce_symmetric();

    // Note: a is killed by realSymmetricEigensystem!
    auto a_copy = a;
    const auto [evals, evecs] = LinAlg::realSymmetricEigensystem(&a_copy);

    // Av = ev => A*v - ev = 0.0
    std::vector<double> each;
    for (std::size_t i = 0; i < evals.n; ++i) {
      const auto v = evecs.get_row(i);
      each.push_back(helper::max_vecel((a * v) - (evals[i] * v)));
    }
    const auto worst =
        *std::max_element(cbegin(each), cend(each), qip::comp_abs);
    pass &= qip::check_value(&obuff, "Eigen(RS) Av = ev", worst, 0.0, 1.0e-12);
  }

  { // RS Eigensystems, Av = eBv
    LinAlg::SqMatrix a(50);
    LinAlg::SqMatrix b(a.n);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-2.0, 2.0);
    for (std::size_t i = 0; i < a.n; ++i) {
      for (std::size_t j = 0; j < a.n; ++j) {
        a[i][j] = dis(gen);
        b[i][j] = 0.01 * dis(gen);
        if (i == j) {
          a[i][j] += double(10 + 2 * i);
          b[i][j] += 1.5;
        }
      }
    }
    a.enforce_symmetric();
    b.enforce_symmetric();

    // Note: a,b is killed by realSymmetricEigensystem!
    auto a_copy = a;
    auto b_copy = b;
    const auto [evals, evecs] =
        LinAlg::realSymmetricEigensystem(&a_copy, &b_copy);

    // Av = eBv => B^-1*A*v = ev
    std::vector<double> each;
    auto binv = b.inverse();
    for (std::size_t i = 0; i < evals.n; ++i) {
      const auto v = evecs.get_row(i);
      each.push_back(helper::max_vecel((binv * a * v) - (evals[i] * v)));
    }
    const auto worst =
        *std::max_element(cbegin(each), cend(each), qip::comp_abs);
    pass &= qip::check_value(&obuff, "Eigen(RS) Av = eBv", worst, 0.0, 1.0e-12);
  }

  //****************************************************************************
  // Complex matrix:
  {
      // No det method!?
  }

  { // Inverse:
    // also test: plusIdent, make_identity, mult* and -
    LinAlg::ComplexSqMatrix m(20);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-0.1, 0.1);
    for (std::size_t i = 0; i < m.n; ++i) {
      for (std::size_t j = 0; j < m.n; ++j) {
        m[i][j] = LinAlg::ComplexDouble{dis(gen), dis(gen)}.val;
      }
    }
    m.plusIdent(1.0, 0.5); // "ensure" invertable  m += (1.0 + 0.5i)I
    auto minv = m.inverse();
    auto minv_m = minv * m;
    LinAlg::ComplexSqMatrix ident(m.n);
    ident.make_identity();
    const auto [delr, deli] = helper::max_matrixel(minv_m - ident);
    pass &= qip::check_value(&obuff, "Complex inverse + matrix-mult",
                             std::max(std::abs(delr), std::abs(deli)), 0.0,
                             1.0e-14);
  }

  { // Transpose
    LinAlg::ComplexSqMatrix m(75);
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-1.0e5, 1.0e5);
    for (std::size_t i = 0; i < m.n; ++i) {
      for (std::size_t j = 0; j < m.n; ++j) {
        m[i][j] = LinAlg::ComplexDouble{dis(gen), dis(gen)}.val;
      }
    }

    const auto mtr = m.transpose();
    LinAlg::ComplexDouble sym_sum = {0.0, 0.0};
    for (std::size_t i = 0; i < m.n; ++i) {
      for (std::size_t j = 0; j < m.n; ++j) {
        // should all be exactly zero
        sym_sum += m.get_copy(i, j) - mtr.get_copy(j, i);
      }
    }
    const auto [re, im] = sym_sum.unpack();
    pass &=
        qip::check_value(&obuff, "Complex transpose",
                         std::max(std::abs(re), std::abs(im)), 0.0, 1.0e-16);
  }

  { // Multiply elements
    LinAlg::ComplexSqMatrix m1(50);
    LinAlg::ComplexSqMatrix m2(m1.n);
    LinAlg::ComplexSqMatrix m3(m1.n); // m3_ij = m1_ij * m2_ij
    LinAlg::ComplexSqMatrix m4(m1.n); // m4 = x * m1;
    const LinAlg::ComplexDouble x = {3.1415926, 2.71828};
    std::mt19937 gen(0.0); // seeded w/ constant; same each run
    std::uniform_real_distribution<> dis(-1.0e5, 1.0e5);
    for (std::size_t i = 0; i < m1.n; ++i) {
      for (std::size_t j = 0; j < m1.n; ++j) {
        const LinAlg::ComplexDouble v1 = {dis(gen), dis(gen)};
        const LinAlg::ComplexDouble v2 = {dis(gen), dis(gen)};
        m1[i][j] = v1.val;
        m2[i][j] = v2.val;
        m3[i][j] = (v1 * v2).val;
        m4[i][j] = (x * v1).val;
      }
    }

    // Expect: mult_elements(m1, m2) == m3
    const auto m5 = LinAlg::ComplexSqMatrix::mult_elements(m1, m2) - m3;

    const auto [re, im] = helper::max_matrixel(m5);
    pass &=
        qip::check_value(&obuff, "Complex Multiply els",
                         std::max(std::abs(re), std::abs(im)), 0.0, 1.0e-16);

    // Also test scalar mult: expect m6 = m4
    const auto m6 = x * m1;
    const auto [re2, im2] = helper::max_matrixel(m6 - m4);
    pass &=
        qip::check_value(&obuff, "Complex Scalar mult",
                         std::max(std::abs(re2), std::abs(im2)), 0.0, 1.0e-16);

    const auto m1r = m1.real();
    const auto m1i = m1.imaginary();
    const auto m1C = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, m1r) +
                     LinAlg::ComplexSqMatrix::make_complex({0.0, 1.0}, m1i);
    const auto [re3, im3] = helper::max_matrixel(m1C - m1);
    pass &=
        qip::check_value(&obuff, "Complex decompose",
                         std::max(std::abs(re3), std::abs(im3)), 0.0, 1.0e-16);
  }

  return pass;
}

} // namespace UnitTest

//****************************************************************************

inline double UnitTest::helper::max_matrixel(const LinAlg::SqMatrix &m) {
  double max = 0.0;
  for (std::size_t i = 0; i < m.n; ++i) {
    for (std::size_t j = 0; j < m.n; ++j) {
      if (std::abs(m[i][j]) > std::abs(max))
        max = m[i][j];
    }
  }
  return max;
}

inline double UnitTest::helper::max_vecel(const LinAlg::Vector &b) {
  double max = 0.0;
  for (std::size_t i = 0; i < b.n; ++i) {
    if (std::abs(b[i]) > std::abs(max))
      max = b[i];
  }
  return max;
}

inline std::pair<double, double>
UnitTest::helper::max_matrixel(const LinAlg::ComplexSqMatrix &m) {
  double maxre = 0.0;
  double maxim = 0.0;
  for (std::size_t i = 0; i < m.n; ++i) {
    for (std::size_t j = 0; j < m.n; ++j) {
      const auto [re, im] = LinAlg::ComplexDouble(m[i][j]).unpack();
      if (std::abs(re) > std::abs(maxre))
        maxre = re;
      if (std::abs(im) > std::abs(maxim))
        maxim = im;
    }
  }
  return std::pair{maxre, maxim};
}
