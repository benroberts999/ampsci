#include "LinAlg.hpp"
#include "catch2/catch.hpp"
#include "qip/Check.hpp"
#include <algorithm>
#include <complex>
#include <numeric>

inline unsigned int Factorial(unsigned int number) {
  return number <= 1 ? number : Factorial(number - 1) * number;
}

TEST_CASE("Factorials are computed", "[factorial]") {
  REQUIRE(Factorial(1) == 1);
  REQUIRE(Factorial(2) == 2);
  REQUIRE(Factorial(3) == 6);
  REQUIRE(Factorial(10) == 3628800);
}

using namespace LinAlg;

//==========================================================================
// Test access and memory layout (only with double)
TEST_CASE("Element access, memory layout", "[LinAlg]") {
  // Check that access is what is expected:
  Matrix<double> a{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  Matrix<double> a2(2, 3);
  a2(0, 0) = 1.0;
  a2(0, 1) = 2.0;
  a2(0, 2) = 3.0;
  a2(1, 0) = 4.0;
  a2(1, 1) = 5.0;
  a2(1, 2) = 6.0;
  REQUIRE(equal(a, a2));

  // test access operators () and [][]
  for (std::size_t i = 0; i < a.rows(); ++i) {
    for (std::size_t j = 0; j < a.cols(); ++j) {
      REQUIRE(a(i, j) == a[i][j]);
    }
  }

  // Prove memory layout is in correct row-major form
  for (std::size_t i = 0; i < a.rows(); ++i) {
    for (std::size_t j = 0; j < a.cols() - 1; ++j) {
      REQUIRE(&a(i, j + 1) == &a(i, j) + 1);
    }
  }

  // Test gsl_view.matix is not doing a data copy
  const auto gv = a.as_gsl_view();
  REQUIRE(&(&gv.matrix)->data[0] == &(a[0][0]));

  std::cout << "Do print cout?\n";
  INFO("Do print?");
  CAPTURE("Do print2?");

  // Test constructor that takes vector by move:
  std::vector<double> v1{1.0, 1.0, 1.0, 1.0};
  auto mem1 = &v1[0];
  Matrix<double> Ma(2, 2, std::move(v1));
  REQUIRE(&(Ma[0][0]) == mem1);
  // Test move constructore
  const auto Mb = std::move(Ma);
  REQUIRE(&(Mb[0][0]) == mem1);
}