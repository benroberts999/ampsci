#include "Physics/AtomData.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

TEST_CASE("Physics: AtomData", "[AtomData][Physics]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Physics: AtomData, [AtomData][Physics]\n";

  REQUIRE(AtomData::atomicSymbol(55) == "Cs");
  REQUIRE(AtomData::atomicName(55) == "cesium");
  REQUIRE(AtomData::atomic_Z("55") == 55);
  REQUIRE(AtomData::atomic_Z(55) == 55);
  REQUIRE(AtomData::atomic_Z("Cs") == 55);

  std::vector<int> int_l = {0, 1, 2, 3, 4, 5, 6};
  std::vector<std::string> str_l = {"s", "p", "d", "f", "g", "h", "i"};
  for (std::size_t i = 0; i < int_l.size(); ++i) {
    REQUIRE(AtomData::l_symbol(int_l[i]) == str_l[i]);
    REQUIRE(AtomData::symbol_to_l(str_l[i]) == int_l[i]);
  }

  REQUIRE(AtomData::kappa_symbol(-1) == "s_1/2");
  REQUIRE(AtomData::kappa_symbol(1) == "p_1/2");
  REQUIRE(AtomData::kappa_symbol(-2) == "p_3/2");

  REQUIRE(AtomData::parse_symbol("s") == (std::pair<int, int>{0, -1}));
  REQUIRE(AtomData::parse_symbol("s+") == (std::pair<int, int>{0, -1}));
  REQUIRE(AtomData::parse_symbol("6s+") == (std::pair<int, int>{6, -1}));
  REQUIRE(AtomData::parse_symbol("6s") == (std::pair<int, int>{6, -1}));
  REQUIRE(AtomData::parse_symbol("12d-") == (std::pair<int, int>{12, 2}));
  REQUIRE(AtomData::parse_symbol("12d_3/2") == (std::pair<int, int>{12, 2}));
  REQUIRE(AtomData::parse_symbol("d_3/2") == (std::pair<int, int>{0, 2}));

  // note: NOT COMPLETE tests
}
