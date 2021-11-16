#pragma once
#include "LinAlg/LinAlg.hpp"
#include "qip/Check.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>

namespace UnitTest {

//******************************************************************************

bool LinAlg(std::ostream & /*obuff*/) {
  using namespace LinAlg;
  bool pass = true;

  //**************************************************************************
  // Test access and memory layout (only with double)
  {
    // Check that access is what is expected:
    Matrix<double> a{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    Matrix<double> a2(2, 3);
    a2(0, 0) = 1.0;
    a2(0, 1) = 2.0;
    a2(0, 2) = 3.0;
    a2(1, 0) = 4.0;
    a2(1, 1) = 5.0;
    a2(1, 2) = 6.0;
    assert(equal(a, a2));

    // test access operators () and [][]
    for (std::size_t i = 0; i < a.rows(); ++i) {
      for (std::size_t j = 0; j < a.cols(); ++j) {
        assert(a(i, j) == a[i][j]);
      }
    }

    // Prove memory layout is in correct row-major form
    for (std::size_t i = 0; i < a.rows(); ++i) {
      for (std::size_t j = 0; j < a.cols() - 1; ++j) {
        assert(&a(i, j + 1) == &a(i, j) + 1);
      }
    }

    // Test gsl_view.matix is not doing a data copy
    const auto gv = a.as_gsl_view();
    assert(&(&gv.matrix)->data[0] == &(a[0][0]));

    // Test constructor that takes vector by move:
    std::vector<double> v1{1.0, 1.0, 1.0, 1.0};
    auto mem1 = &v1[0];
    Matrix<double> Ma(2, 2, std::move(v1));
    assert(&(Ma[0][0]) == mem1);
    // Test move constructore
    const auto Mb = std::move(Ma);
    assert(&(Mb[0][0]) == mem1);
  }

  {
    Matrix<int> a{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    assert(std::accumulate(a.row_begin(0), a.row_end(0), 0) == 6);
    assert(std::accumulate(a.row_begin(1), a.row_end(1), 0) == 15);
    assert(std::accumulate(a.row_begin(2), a.row_end(2), 0) == 24);
    assert(std::accumulate(a.col_begin(0), a.col_end(0), 0) == 12);
    assert(std::accumulate(a.col_begin(1), a.col_end(1), 0) == 15);
    assert(std::accumulate(a.col_begin(2), a.col_end(2), 0) == 18);
  }

  //**************************************************************************
  // Test operators, transpose, and "mult elements" for double
  {
    Matrix<double> a{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};

    // test addition + multipcation
    assert(equal(a + a, 2.0 * a));
    assert(equal(2.0 * a, a * 2.0));
    assert(equal(a + 0.5 * a, 1.5 * a));
    assert(equal(a - 0.5 * a, 0.5 * a));
    assert(equal(a / 4.5, (1.0 / 4.5) * a));

    auto b = a;
    const auto b2 = a;
    auto b3 = b2;
    assert(equal(b, a));
    assert(equal(b2, a));
    assert(equal(b3, a));
    b3 = 2.0 * a;
    assert(equal(b3, 2.0 * a));

    const Matrix m{{{1.0, 0.5, 0.1}, {0.5, 1.0, 0.5}, {0.1, 0.5, 1.0}}};
    assert(m.determinant() == 0.54);

    assert(equal(m.inverse(), Matrix{{25.0 / 18, -5.0 / 6, 5.0 / 18},
                                     {-5.0 / 6, 11.0 / 6, -5.0 / 6},
                                     {5.0 / 18, -5.0 / 6, 25.0 / 18}}));

    Matrix<double> c{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
    assert(equal(c.transpose(), {{1.0, 4.0}, {2.0, 5.0}, {3.0, 6.0}}));

    auto d = mult_elements(a, a);
    assert(equal(d, {{1.0, 4.0, 9.0}, {16.0, 25.0, 36.0}}));
    a.mult_elements_by(a);
    assert(equal(a, d));

    assert(equal(a.make_identity(), {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}));
    assert(equal(a.plusIdent(), {{2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}}));
    assert(equal(a.plusIdent(2.0), {{4.0, 0.0, 0.0}, {0.0, 4.0, 0.0}}));
    assert(equal(a.zero(), {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}));
  }

  //**************************************************************************
  // Test operators, transpose, and "mult elements" for int
  {

    Matrix<int> a{{1, 2, 3}, {4, 5, 6}};
    // test addition + multipcation
    assert(equal(a + a, 2 * a));
    assert(equal(2 * a, a * 2));
    assert(equal(a + 2 * a, 3 * a));
    assert(equal(a - 2 * a, -1 * a));
    assert(equal(4 * a / 2, 2 * a)); // careful with int devision, as always

    auto b = a;
    const auto b2 = a;
    auto b3 = b2;
    assert(equal(b, a));
    assert(equal(b2, a));
    assert(equal(b3, a));
    b3 = 2 * a;
    assert(equal(b3, 2 * a));

    Matrix<int> c{{1, 2, 3}, {4, 5, 6}};
    assert(equal(c.transpose(), {{1, 4}, {2, 5}, {3, 6}}));

    auto d = mult_elements(a, a);
    assert(equal(d, {{1, 4, 9}, {16, 25, 36}}));
    a.mult_elements_by(a);
    assert(equal(a, d));
  }

  //**************************************************************************
  // Test operators, transpose for float
  {

    Matrix<float> a{{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}};

    assert(equal(a + a, 2.0f * a));
    assert(equal(a + 0.5f * a, 1.5f * a));
    assert(equal(a - 0.5f * a, 0.5f * a));
    assert(equal(a / 4.5f, (1.0f / 4.5f) * a));
    auto b = a;
    const auto b2 = a;
    auto b3 = b2;
    assert(equal(b, a));
    assert(equal(b2, a));
    assert(equal(b3, a));
    b3 = 2.0f * a;
    assert(equal(b3, 2.0f * a));

    Matrix<float> c{{1.0f, 2.0f, 3.0f}, {4.0f, 5.0f, 6.0f}};
    assert(equal(c.transpose(), {{1.0f, 4.0f}, {2.0f, 5.0f}, {3.0f, 6.0f}}));

    const auto c_tr = c.transpose();
    assert(equal(c_tr, c.transpose()));

    assert(equal(a.make_identity(), {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}}));
    assert(equal(a.plusIdent(), {{2.0f, 0.0f, 0.0f}, {0.0f, 2.0f, 0.0f}}));
    assert(equal(a.plusIdent(2.0f), {{4.0f, 0.0f, 0.0f}, {0.0f, 4.0f, 0.0f}}));
    assert(equal(a.zero(), {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}}));

    // Invert and det only implemented for doubles
  }

  //**************************************************************************
  // Test operators, transpose, and GSL matrix_view for complex double
  {

    using namespace std::complex_literals;
    const Matrix<std::complex<double>> a{{1.5 - 1.i, 2.5 + 2.5i},
                                         {3.5 + 3.5i, 4.5 - 4.5i}};

    // test addition and multiplication
    assert(equal(a + a, 2.0 * a));
    assert(equal(a + 0.5 * a, 1.5 * a));
    assert(equal(a - 0.5 * a, 0.5 * a));
    assert(equal(a / 4.5, (1.0 / 4.5) * a));

    assert(equal(a * (4.0 + 3.0i),
                 {{9. + 0.5i, 2.5 + 17.5i}, {3.5 + 24.5i, 31.5 - 4.5i}}));

    auto b = a;
    const auto b2 = a;
    auto b3 = b2;
    assert(equal(b, a));
    assert(equal(b2, a));
    assert(equal(b3, a));
    b3 = (2.0 + 3.0i) * a;
    assert(equal(b3, (2.0 + 3.0i) * a));

    // Test gsl_view.matix is not doing a data copy
    const auto gv = a.as_gsl_view();
    // have to cast pointer to complex<double>* -> double*
    // Garuenteed to be pointer to real part
    assert(&(&gv.matrix)->data[0] ==
           reinterpret_cast<const double *>(&(a[0][0])));

    auto a_tr = a.transpose();
    assert(equal(a_tr, {{1.5 - 1.i, 3.5 + 3.5i}, {2.5 + 2.5i, 4.5 - 4.5i}}));
    const auto ca_tr = a_tr;
    assert(equal(ca_tr, a.transpose()));

    assert(equal(a.inverse(), {{0.167743874943634 + 0.143393957613107i,
                                0.0796633097850594 - 0.0931910416353525i},
                               {0.111528633699083 - 0.130467458289493i,
                                0.0386291898391703 + 0.0491507590560649i}}));
    const auto ca_inv = a.inverse();
    assert(equal(ca_inv, a.inverse()));

    assert(a.determinant() == 2.25 - 28.75i);
    assert(std::abs(ca_inv.determinant() -
                    (0.0027055463700586 + 0.0345708702840824i)) < 1.0e-12);

    assert(
        equal(a.conj(), {{1.5 + 1.i, 2.5 - 2.5i}, {3.5 - 3.5i, 4.5 + 4.5i}}));
    assert(equal(a.real(), {{1.5, 2.5}, {3.5, 4.5}}));
    assert(equal(a.imag(), {{-1.0, 2.5}, {3.5, -4.5}}));

    auto a2 = a;
    assert(equal(a2.make_identity(),
                 {{1.0 + 0.0i, 0.0 + 0.0i}, {0.0 + 0.0i, 1.0 + 0.0i}}));
    assert(equal(a2.plusIdent(), {{2.0, 0.0}, {0.0, 2.0}}));
    assert(equal(a2.plusIdent(2.0i), {{2.0 + 2.0i, 0.0}, {0.0, 2.0 + 2.0i}}));
    assert(equal(a2.zero(), {{0.0, 0.0}, {0.0, 0.0}}));
  }

  //**************************************************************************
  // Test operators and transpose for complex float
  {
    using namespace std::complex_literals;
    const Matrix<std::complex<float>> a{{1.5f - 1.if, 2.5f + 2.5if},
                                        {3.5f + 3.5if, 4.5f - 4.5if}};

    // test addition and multiplication
    assert(equal(a + a, 2.0f * a));
    assert(equal(a + 0.5f * a, 1.5f * a));
    assert(equal(a - 0.5f * a, 0.5f * a));
    assert(equal(a / 4.5f, (1.0f / 4.5f) * a));

    assert(equal(a * (4.0f + 3.0if), {{9.f + 0.5if, 2.5f + 17.5if},
                                      {3.5f + 24.5if, 31.5f - 4.5if}}));

    auto b = a;
    const auto b2 = a;
    auto b3 = b2;
    assert(equal(b, a));
    assert(equal(b2, a));
    assert(equal(b3, a));
    b3 = (2.0f + 3.0if) * a;
    assert(equal(b3, (2.0f + 3.0if) * a));

    auto a_tr = a.transpose();
    assert(equal(a_tr,
                 {{1.5f - 1.if, 3.5f + 3.5if}, {2.5f + 2.5if, 4.5f - 4.5if}}));
    const auto ca_tr = a_tr;

    assert(equal(ca_tr, a.transpose()));

    // Invert and det only implemented for doubles
  }

  //**************************************************************************
  {
    // Test matrix multiplication, with double, const and non-const
    Matrix a{{-1.0, 6.5, 0.1, 1.0, -4.0},
             {0.5, 3.0, 0.5, 2.0, 14.0},
             {0.1, 0.5, -1.0, 3.0, 8.0},
             {14.1, 0.5, 1.0, 7.0, 9.0}};

    Matrix b{{6.0, 77.0, 2.1, 1.7, 6.0},
             {2.0, 4.1, -2.2, 2.7, -6.0},
             {9.0, 15.8, 2.3, 3.7, 6.0},
             {6.1, -1.0, 2.4, 4.7, -6.0},
             {6.1, 0.4, 2.5, 5.7, 1.0}};

    assert(equal(a * b, {{-10.4, -51.37, -23.77, -1.88, -54.4},
                         {111.1, 62.3, 35.4, 100., -10.},
                         {59.7, -5.85, 24.01, 57.52, -18.4},
                         {192.2, 1100.15, 70.11, 113.22, 54.6}}));
    const auto a2 = a;
    const auto b2 = b;
    assert(equal(a2 * b2, a * b));
  }
  {
    // Test matrix multiplication, with float, const and non-const
    Matrix<float> a{{-1.0f, 6.5f, 0.1f, 1.0f, -4.0f},
                    {0.5f, 3.0f, 0.5f, 2.0f, 14.0f},
                    {0.1f, 0.5f, -1.0f, 3.0f, 8.0f},
                    {14.1f, 0.5f, 1.0f, 7.0f, 9.0f}};

    Matrix<float> b{{6.0f, 77.0f, 2.1f, 1.7f, 6.0f},
                    {2.0f, 4.1f, -2.2f, 2.7f, -6.0f},
                    {9.0f, 15.8f, 2.3f, 3.7f, 6.0f},
                    {6.1f, -1.0f, 2.4f, 4.7f, -6.0f},
                    {6.1f, 0.4f, 2.5f, 5.7f, 1.0f}};

    assert(equal(a * b,
                 Matrix<float>{{-10.4f, -51.37f, -23.77f, -1.88f, -54.4f},
                               {111.1f, 62.3f, 35.4f, 100.f, -10.0f},
                               {59.7f, -5.85f, 24.01f, 57.52f, -18.4f},
                               {192.2f, 1100.15f, 70.11f, 113.22f, 54.6f}}));

    const auto a2 = a;
    const auto b2 = b;
    assert(equal(a2 * b2, a * b));
  }
  {
    // Test matrix multiplication, with complex double, const and non-const
    using namespace std::complex_literals;
    Matrix a{
        {-1.0 + 2.0i + 2.0i, 6.5 + 1.0i, 0.1 + 2.0i, 1.0 + 2.0i, -4.0 + 6.0i},
        {0.5 + 3.0i, 3.0 - 2.0i, 0.5 + 2.0i, 2.0 - 2.0i, 14.0 + 4.0i},
        {0.1 + 4.0i, 0.5 + 2.0i, -1.0 + 3.0i, 3.0 + 2.0i, 8.0 - 2.0i},
        {14.1 + 5.0i, 0.5 - 2.0i, 1.0 + 4.0i, 7.0 - 2.0i, 9.0 + 1.0i}};

    Matrix b{{6.0 + 1.0i, 77.0 + 6.0i, 2.1 - 1.0i, 1.7 + 2.0i, 6.0 + 1.0i},
             {2.0 + 2.0i, 4.1 + 7.0i, -2.2 - 2.0i, 2.7 + 2.0i, -6.0 + 1.0i},
             {9.0 + 3.0i, 15.8 + 8.0i, 2.3 - 3.0i, 3.7 + 2.0i, 6.0 + 1.0i},
             {6.1 + 4.0i, -1.0 + 9.0i, 2.4 - 4.0i, 4.7 + 2.0i, -6.0 + 1.0i},
             {6.1 + 5.0i, 0.4 + 0.0i, 2.5 - 5.0i, 5.7 + 2.0i, 1.0 + 1.0i}};

    assert(equal(a * b, {{-60.4 + 89.1i, -116.37 + 393.4i, 26.23 + 34.3i,
                          -31.88 + 65.7i, -69.4 + 26.6i},
                         {94.1 + 130.2i, 60.3 + 304.i, 52.4 - 65.5i,
                          90. + 60.5i, -15. + 78.i},
                         {44.7 + 105.1i, -85.85 + 383.9i, 39.01 - 39.4i,
                          39.52 + 42.5i, -27.4 + 26.6i},
                         {182.2 + 147.i, 1070.15 + 601.5i, 80.11 - 69.3i,
                          101.22 + 77.4i, 48.6 + 110.6i}}));

    const auto a2 = a;
    const auto b2 = b;
    assert(equal(a2 * b2, a * b));
  }
  {
    // Test matrix multiplication, with complex float, const and non-const
    using namespace std::complex_literals;
    Matrix<std::complex<float>> a{
        {-1.0f + 2.0if + 2.0if, 6.5f + 1.0if, 0.1f + 2.0if, 1.0f + 2.0if,
         -4.0f + 6.0if},
        {0.5f + 3.0if, 3.0f - 2.0if, 0.5f + 2.0if, 2.0f - 2.0if, 14.0f + 4.0if},
        {0.1f + 4.0if, 0.5f + 2.0if, -1.0f + 3.0if, 3.0f + 2.0if, 8.0f - 2.0if},
        {14.1f + 5.0if, 0.5f - 2.0if, 1.0f + 4.0if, 7.0f - 2.0if,
         9.0f + 1.0if}};

    Matrix b{
        {6.0f + 1.0if, 77.0f + 6.0if, 2.1f - 1.0if, 1.7f + 2.0if, 6.0f + 1.0if},
        {2.0f + 2.0if, 4.1f + 7.0if, -2.2f - 2.0if, 2.7f + 2.0if,
         -6.0f + 1.0if},
        {9.0f + 3.0if, 15.8f + 8.0if, 2.3f - 3.0if, 3.7f + 2.0if, 6.0f + 1.0if},
        {6.1f + 4.0if, -1.0f + 9.0if, 2.4f - 4.0if, 4.7f + 2.0if,
         -6.0f + 1.0if},
        {6.1f + 5.0if, 0.4f + 0.0if, 2.5f - 5.0if, 5.7f + 2.0if, 1.0f + 1.0if}};

    assert(equal(a * b, {{-60.4f + 89.1if, -116.37f + 393.4if, 26.23f + 34.3if,
                          -31.88f + 65.7if, -69.4f + 26.6if},
                         {94.1f + 130.2if, 60.3f + 304.if, 52.4f - 65.5if,
                          90.f + 60.5if, -15.f + 78.if},
                         {44.7f + 105.1if, -85.85f + 383.9if, 39.01f - 39.4if,
                          39.52f + 42.5if, -27.4f + 26.6if},
                         {182.2f + 147.if, 1070.15f + 601.5if, 80.11f - 69.3if,
                          101.22f + 77.4if, 48.6f + 110.6if}}));

    const auto a2 = a;
    const auto b2 = b;
    assert(equal(a2 * b2, a * b));
  }

  // Test matrix-vector operations:
  //**************************************************************************
  {
    Vector<double> x{1.0, -1.0};
    Matrix<double> A{{1.0, 2.0}, {3.0, 4.0}};
    const auto b = A * x;
    assert(equal(b, Vector<double>{-1.0, -1.0}));
    assert(std::abs(b * b - 2.0) < 1.0e-15);
  }
  {
    Vector<float> x{1.0f, -1.0f};
    Matrix<float> A{{1.0f, 2.0f}, {3.0f, 4.0f}};
    const auto b = A * x;
    assert(equal(b, Vector<float>{-1.0f, -1.0f}));
    assert(std::abs(b * b - 2.0f) < 1.0e-6f);
  }
  {
    using namespace std::complex_literals;
    Vector<std::complex<double>> x{1.0 + 2.0i, -1.0 + 2.0i};
    Matrix<std::complex<double>> A{{1.0 - 1.0i, 2.0 + 1.0i},
                                   {3.0 + 1.0i, 4.0 - 1.0i}};
    const auto b = A * x;
    assert(equal(b, Vector<std::complex<double>>{-1.0 + 4.0i, -1.0 + 16.0i}));
    assert(std::abs(x * b - (-40.0 - 16.0i)) < 1.0e-15);

    Vector<std::complex<double>> xdag = x.conj();
    auto xdag2 = x.conj();
    assert(equal(xdag2, xdag));

    assert(equal(x.conj(), Vector{1.0 - 2.0i, -1.0 - 2.0i}));
    assert(equal(x.real(), Vector{1.0, -1.0}));
    assert(equal(x.imag(), Vector{2.0, 2.0}));
  }
  {
    using namespace std::complex_literals;
    const Vector<std::complex<float>> x{1.0f + 2.0if, -1.0f + 2.0if};
    const Matrix<std::complex<float>> A{{1.0f - 1.0if, 2.0f + 1.0if},
                                        {3.0f + 1.0if, 4.0f - 1.0if}};
    const auto b = A * x;
    assert(
        equal(b, Vector<std::complex<float>>{-1.0f + 4.0if, -1.0f + 16.0if}));
    assert(std::abs(x * b - (-40.0f - 16.0if)) < 1.0e-6f);
  }

  {
    const Vector<double> a{1.0, 2.0};
    const Vector<double> b{3.0, 4.0};
    const auto Mab = outer_product(a, b);
    assert(equal(Mab, {{3.0, 4.0}, {6.0, 8.0}}));

    // alternate way of performing outer product:
    const Matrix a2 = a;
    const Matrix b2 = b;
    const auto Mab2 = a2 * b2.transpose();
    assert(equal(Mab, Mab2));

    // alternate way of performing inner product:
    const auto inner1 = a2.transpose() * b2; // returns 1*1 matrix
    assert(std::abs(inner1(0, 0) - a * b) < 1.0e-12);
  }

  // Test real/imag/complex conversions:
  //**************************************************************************
  {
    using namespace std::complex_literals;
    Matrix A{{1.0 + 2.0i, 2.0 - 1.0i}, {3.0 + 3.0i, 4.0 + 1.0i}};
    Vector b{-4.0 + 2.0i, -3.0 + 1.0i};

    const auto rA = A.real();
    const auto iA = A.imag();
    const auto A2 = rA.complex() + 1.0i * iA.complex();
    assert(equal(A, A2));

    const auto rb = b.real();
    const auto ib = b.imag();
    const auto b2 = rb.complex() + 1.0i * ib.complex();
    assert(equal(b, b2));
  }

  //**************************************************************************

  //**************************************************************************
  // Solve linear equation:
  {
    Matrix<double> A{{1.0, 2.0}, {3.0, 4.0}};
    Vector<double> b{-1.0, -1.0};
    Vector<double> x_sol{1.0, -1.0};

    auto x = solve_Axeqb(A, b);
    assert(equal(x, x_sol));
  }
  {
    using namespace std::complex_literals;
    Matrix A{{1.0 + 2.0i, 2.0 - 1.0i}, {3.0 + 3.0i, 4.0 + 1.0i}};
    Vector b{-4.0 + 2.0i, -3.0 + 1.0i};
    Vector x_sol{1.0 + 1.0i, -1.0 - 1.0i};

    auto x = solve_Axeqb(A, b);
    assert(equal(x, x_sol));
  }

  //**************************************************************************
  // Eigensystems (symmetric/Hermetian)
  {
    const Matrix A{{1.0, -1.0}, {-1.0, 2.0}};
    const auto [e, v] = symmhEigensystem(A, true);
    assert(equal(e, Vector{0.381966011250105, 2.61803398874989}));
    assert(equal(v, Matrix{{0.850650808352040, 0.525731112119133},
                           {-0.525731112119133, 0.850650808352040}}));
  }
  {
    using namespace std::complex_literals;
    const Matrix A{{1.0 + 0.0i, 0.0 + 1.0i}, {0.0 - 1.0i, -1.0 + 0.0i}};
    const auto [e, v] = symmhEigensystem(A, true);

    const auto root2 = std::sqrt(2.0);
    assert(equal(e, Vector{-root2, root2}));
    assert(equal(v, Matrix<std::complex<double>>{
                        {0.382683432365090, 0.923879532511287i},
                        {0.923879532511287, -0.382683432365090i}}));
  }

  // Generalised Eigensystems (symmetric/Hermetian)
  {
    const Matrix A{{1.0, -1.0}, {-1.0, 2.0}};
    const Matrix B{{1.0, 0.1}, {0.1, 1.0}};
    const auto [e, v] = symmhEigensystem(A, B, true);
    assert(equal(e, Vector{0.350508678167508, 2.88181455415572}));
    assert(equal(v, Matrix{{0.847046341910562, 0.531519044490351},
                           {-0.564870314672654, 0.825179694128265}}));
  }

  {
    using namespace std::complex_literals;
    const Matrix<std::complex<double>> A{{1.0, -1.0 + 0.5i},
                                         {-1.0 - 0.5i, -1.0}};
    const Matrix<std::complex<double>> B{{1.0, 0.5 + 0.1i}, {0.5 - 0.1i, 1.0}};
    const auto [e, v] = symmhEigensystem(A, B, true);
    assert(equal(e, Vector<double>{-1.238601401177898, 2.45481761739411}));

    // Check Av = eBv => Av - eBv = 0
    std::complex<double> sum = 0;
    for (auto n = 0ul; n < A.rows(); ++n) {
      for (auto i = 0ul; i < A.rows(); ++i) {
        for (auto j = 0ul; j < A.rows(); ++j) {
          sum += A(i, j) * v(n, j);
        }
        for (auto j = 0ul; j < A.rows(); ++j) {
          sum -= e(n) * B(i, j) * v(n, j);
        }
      }
    }
    assert(std::abs(sum) < 1.0e-12);
  }

  // Eigensystems (non-symmetric, Real)
  {
    const Matrix A{{1.0, -1.0}, {-2.0, 3.0}};
    const auto [e, v] = genEigensystem(A, true);
    assert(equal(e.real(), Vector{0.267949192431123, 3.73205080756888}));
    assert(equal(v.real(), Matrix{{-0.806898221355073, -0.590690494568872},
                                  {0.343723769333440, -0.939070801588044}}));
  }

  //**************************************************************************
  //**************************************************************************
  //**************************************************************************

  return pass;
}

} // namespace UnitTest
