#include "qip/String.hpp"
#include "catch2/catch.hpp"
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

  std::vector<std::string> list = {"ben", "writes", "a", "list"};
  REQUIRE(qip::Levenstein("ben", "ben") == 0);
  REQUIRE(qip::ci_Levenstein("ben", "BEN") == 0);
  REQUIRE(*qip::closest_match("wrote", list) == "writes");
  REQUIRE(*qip::ci_closest_match("WrIt", list) == "writes");
  REQUIRE(*qip::closest_match("ben", list) == "ben");
  REQUIRE(*qip::closest_match("bon", list) == "ben");
  REQUIRE(*qip::closest_match("a", list) == "a");
  REQUIRE(*qip::closest_match("b", list) == "a");
  REQUIRE(*qip::closest_match("list", list) == "list");
  REQUIRE(*qip::closest_match("lost", list) == "list");
}
