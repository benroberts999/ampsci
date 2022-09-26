#include "Random.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <string>

TEST_CASE("qip::Random", "[qip][Random][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "qip::Random\n";

  const auto s1 = qip::random_string(6);
  REQUIRE(s1.size() == 6);

  // Very unlikely two 50 lengths should be same:
  const auto s2 = qip::random_string(50);
  const auto s3 = qip::random_string(50);
  REQUIRE_FALSE(s2 == s3);
}