#include "Maths.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <vector>

TEST_CASE("qip::Maths", "[qip][Maths][unit]") {

  const std::vector v_ints{1, -4, 3, -10, 0, 99, -100};
  REQUIRE(*std::min_element(cbegin(v_ints), cend(v_ints), qip::less_abs{}) ==
          0);
  REQUIRE(*std::max_element(cbegin(v_ints), cend(v_ints), qip::less_abs{}) ==
          -100);

  const std::vector v_doubles{1.0, -4.0, 3.0, -10.0, 0.0, 99.0, -100.0};
  REQUIRE(*std::min_element(cbegin(v_doubles), cend(v_doubles),
                            qip::less_abs{}) == 0.0);
  REQUIRE(*std::max_element(cbegin(v_doubles), cend(v_doubles),
                            qip::less_abs{}) == -100.0);

  REQUIRE(qip::max(1, -4, 3, -10, 0, 99, -100) == 99);
  REQUIRE(qip::max_abs(1, -4, 3, -10, 0, 99, -100) == -100);
  REQUIRE(qip::min(1, -4, 3, -10, 0, 99, -100) == -100);
  REQUIRE(qip::min_abs(1, -4, 3, -10, 0, 99, -100) == 0);

  REQUIRE(qip::max(1.0, -4.0, 3.0, -10.0, 0.0, 99.0, -100.0) == 99.0);
  REQUIRE(qip::max_abs(1.0, -4.0, 3.0, -10.0, 0.0, 99.0, -100.0) == -100.0);
  REQUIRE(qip::min(1.0, -4.0, 3.0, -10.0, 0.0, 99.0, -100.0) == -100.0);
  REQUIRE(qip::min_abs(1.0, -4.0, 3.0, -10.0, 0.0, 99.0, -100.0) == 0.0);

  REQUIRE(qip::max_difference(1, -4, 3, -10, 0, 99, -100) == (99 - -100));
  REQUIRE(qip::max_difference(1.0, -4.0, 3.0, -10.0, 0.0, 99.0, -100.0) ==
          Approx(99.0 - -100.0));

  REQUIRE(qip::pow<3>(2) == 8);
  REQUIRE(qip::pow<3>(-2) == -8);
  REQUIRE(qip::pow<-3>(2) == Approx(0.125));
  REQUIRE(qip::pow<3>(0) == 0);
  REQUIRE(qip::pow<0>(2) == 1);

  REQUIRE(qip::pow<3>(2.0) == Approx(8.0));
  REQUIRE(qip::pow<3>(-2.0) == Approx(-8.0));
  REQUIRE(qip::pow<-3>(2.0) == Approx(0.125));
  REQUIRE(qip::pow<3>(0.0) == Approx(0.0));
  REQUIRE(qip::pow<0>(2.0) == Approx(1.0));

  REQUIRE(qip::pow(2.0, 3) == Approx(8));
  REQUIRE(qip::pow(-2.0, 3) == Approx(-8));
  REQUIRE(qip::pow(2.0, -3) == Approx(0.125));
  REQUIRE(qip::pow(0.0, 3) == Approx(0));
  REQUIRE(qip::pow(2.0, 0) == Approx(1));

  REQUIRE(qip::factorial(0) == Approx(1.0));
  REQUIRE(qip::factorial(0ul) == Approx(1.0));
  REQUIRE(qip::factorial(1) == Approx(1.0));
  REQUIRE(qip::factorial(3) == Approx(6.0));
  REQUIRE(qip::factorial(5) == Approx(120.0));
  REQUIRE(qip::factorial(8) == Approx(40320.0));
  REQUIRE(qip::factorial(12) == Approx(479001600.0));
  REQUIRE(qip::factorial(12u) == Approx(479001600.0));
  REQUIRE(qip::factorial(12l) == Approx(479001600.0));
  REQUIRE(qip::factorial(20ul) == Approx(2432902008176640000.0));

  REQUIRE(qip::double_factorial(0) == Approx(1.0));
  REQUIRE(qip::double_factorial(0ul) == Approx(1.0));
  REQUIRE(qip::double_factorial(1) == Approx(1.0));
  REQUIRE(qip::double_factorial(2) == Approx(2.0));
  REQUIRE(qip::double_factorial(5) == Approx(15.0));
  REQUIRE(qip::double_factorial(15) == Approx(2027025.0));
  REQUIRE(qip::double_factorial(20ul) == Approx(3715891200.0));

  REQUIRE(qip::sign(-0) == 0);
  REQUIRE(qip::sign(0) == 0);
  REQUIRE(qip::sign(1) == 1);
  REQUIRE(qip::sign(-1) == -1);
  REQUIRE(qip::sign(16) == 1);
  REQUIRE(qip::sign(-16) == -1);

  REQUIRE(qip::sign(-0.0) == 0);
  REQUIRE(qip::sign(0.0) == 0);
  REQUIRE(qip::sign(1.0) == 1);
  REQUIRE(qip::sign(-1.0) == -1);
  REQUIRE(qip::sign(16.0) == 1);
  REQUIRE(qip::sign(-16.0) == -1);

  REQUIRE(qip::clamp_abs(6, 10) == 6);
  REQUIRE(qip::clamp_abs(12, 10) == 10);
  REQUIRE(qip::clamp_abs(-6, 10) == -6);
  REQUIRE(qip::clamp_abs(-12, 10) == -10);

  REQUIRE(qip::clamp_abs(6.0, 10.0) == 6.0);
  REQUIRE(qip::clamp_abs(12.0, 10.0) == 10.0);
  REQUIRE(qip::clamp_abs(-6.0, 10.0) == -6.0);
  REQUIRE(qip::clamp_abs(-12.0, 10.0) == -10.0);

  REQUIRE(qip::chop(2, 2) == 2);
  REQUIRE(qip::chop(3, 2) == 3);
  REQUIRE(qip::chop(1, 2) == 0);
  REQUIRE(qip::chop(-2, 2) == -2);
  REQUIRE(qip::chop(-3, 2) == -3);
  REQUIRE(qip::chop(-1, 2) == 0);

  REQUIRE(qip::chop(2.0, 2.0) == 2.0);
  REQUIRE(qip::chop(3.0, 2.0) == 3.0);
  REQUIRE(qip::chop(1.0, 2.0) == 0.0);
  REQUIRE(qip::chop(-2.0, 2.0) == -2.0);
  REQUIRE(qip::chop(-3.0, 2.0) == -3.0);
  REQUIRE(qip::chop(-1.0, 2.0) == 0.0);
}