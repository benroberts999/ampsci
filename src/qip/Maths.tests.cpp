#include "Maths.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <vector>

TEST_CASE("qip::Maths", "[qip][Maths][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::Maths\n";

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

  REQUIRE(qip::factorial(0) == 1);
  REQUIRE(qip::factorial(0ul) == 1ul);
  REQUIRE(qip::factorial(1) == 1);
  REQUIRE(qip::factorial(3) == 6);
  REQUIRE(qip::factorial(5) == 120);
  REQUIRE(qip::factorial(8) == 40320);
  REQUIRE(qip::factorial(12) == 479001600);
  REQUIRE(qip::factorial(12u) == 479001600u);
  REQUIRE(qip::factorial(12l) == 479001600l);
  REQUIRE(qip::factorial(20ul) == 2432902008176640000ul);

  // This demonstrates the overflow!
  REQUIRE(qip::factorial(21ul) == 51090942171709440000ul);
  REQUIRE(qip::factorial(21ul) == 14197454024290336768ul);

  REQUIRE(qip::double_factorial(0) == 1);
  REQUIRE(qip::double_factorial(0ul) == 1ul);
  REQUIRE(qip::double_factorial(1) == 1);
  REQUIRE(qip::double_factorial(2) == 2);
  REQUIRE(qip::double_factorial(5) == 15);
  REQUIRE(qip::double_factorial(15) == 2027025);
  REQUIRE(qip::double_factorial(20ul) == 3715891200ul);

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

  REQUIRE(qip::clip(6, 10) == 6);
  REQUIRE(qip::clip(12, 10) == 10);
  REQUIRE(qip::clip(-6, 10) == -6);
  REQUIRE(qip::clip(-12, 10) == -10);

  REQUIRE(qip::clip(6.0, 10.0) == 6.0);
  REQUIRE(qip::clip(12.0, 10.0) == 10.0);
  REQUIRE(qip::clip(-6.0, 10.0) == -6.0);
  REQUIRE(qip::clip(-12.0, 10.0) == -10.0);

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