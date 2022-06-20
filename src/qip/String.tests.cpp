#include "catch2/catch.hpp"
#include "qip/String.hpp"
#include <cassert>
#include <iostream>

TEST_CASE("qip::String", "[qip][String][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::String\n";

  REQUIRE(qip::wildcard_compare("helloben", "helloben") == true);
  REQUIRE(qip::wildcard_compare("helloben", "*ben") == true);
  REQUIRE(qip::wildcard_compare("helloben", "hello*") == true);
  REQUIRE(qip::wildcard_compare("helloben", "hello*ben") == true);
  REQUIRE(qip::wildcard_compare("helloben", "hell*ben") == true);
  REQUIRE(qip::wildcard_compare("helloben", "he*en") == true);
  REQUIRE(qip::wildcard_compare("helloben", "*") == true);
  REQUIRE(qip::wildcard_compare("helloben", "*helloben") == true);
  REQUIRE(qip::wildcard_compare("helloben", "helloben*") == true);
  REQUIRE(qip::wildcard_compare("", "*") == true);
  REQUIRE(qip::wildcard_compare("", "") == true);

  REQUIRE(qip::wildcard_compare("helloben", "xhelloben") == false);
  REQUIRE(qip::wildcard_compare("helloben", "x*ben") == false);
  REQUIRE(qip::wildcard_compare("helloben", "*benx") == false);
  REQUIRE(qip::wildcard_compare("helloben", "xhello*") == false);
  REQUIRE(qip::wildcard_compare("helloben", "hello*x") == false);
  REQUIRE(qip::wildcard_compare("helloben", "*x") == false);
  REQUIRE(qip::wildcard_compare("helloben", "x*") == false);
  REQUIRE(qip::wildcard_compare("helloben", "**") == false);
  REQUIRE(qip::wildcard_compare("helloben", "") == false);
  REQUIRE(qip::wildcard_compare("", "helloben") == false);

  REQUIRE(qip::ci_compare("hello", "hello") == true);
  REQUIRE(qip::ci_compare("hello", "hellob") == false);
  REQUIRE(qip::ci_compare("hellob", "hello") == false);
  REQUIRE(qip::ci_compare("hello", "hEllO") == true);
  REQUIRE(qip::ci_compare("hEllO", "hello") == true);
  REQUIRE(qip::ci_compare("hello!", "hello!") == true);
  REQUIRE(qip::ci_compare("hello!", "hello?") == false);
  REQUIRE(qip::ci_compare("!@#$%&*()-_=+~`,./<>?[]{};':'",
                          "!@#$%&*()-_=+~`,./<>?[]{};':'") == true);
  REQUIRE(qip::ci_compare("!@#$%&*()-_=+~`,./<>?[]{};':'",
                          "!@#$%&*()-_=+~`,./<>?[]{};':'x") == false);

  REQUIRE(qip::fstring("%6.2e %3i %5s", 19517.123524, 2, "ben") ==
          std::string("1.95e+04   2   ben"));
}
