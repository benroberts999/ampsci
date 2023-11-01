#include "Maths/NumCalc_quadIntegrate.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <vector>

TEST_CASE("NumCalc_quadIntegrate", "[num_integrate][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "NumCalc_quadIntegrate\n";

  REQUIRE(NumCalc::num_integrate([](double) { return 1.0; }, 0.0, 10.0, 100,
                                 NumCalc::linear) == 10.0);

  REQUIRE(NumCalc::num_integrate([](double) { return 1.0; }, 1.0e-6, 10.0, 100,
                                 NumCalc::logarithmic) == Approx(10.0));

  REQUIRE(NumCalc::num_integrate([](double x) { return x; }, 0.0, 10.0, 100,
                                 NumCalc::linear) == Approx(50.0));

  REQUIRE(NumCalc::num_integrate([](double x) { return x; }, 1.0e-6, 10.0, 100,
                                 NumCalc::logarithmic) == Approx(50.0));

  REQUIRE(NumCalc::num_integrate([](double x) { return 1.0 / x; }, 1.0,
                                 std::exp(3.5), 500,
                                 NumCalc::linear) == Approx(3.5));

  REQUIRE(NumCalc::num_integrate([](double x) { return 1.0 / x; }, 1.0,
                                 std::exp(3.5), 500,
                                 NumCalc::logarithmic) == Approx(3.5));
}