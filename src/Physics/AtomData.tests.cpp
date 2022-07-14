#include "AtomData.hpp"
#include "catch2/catch.hpp"
#include "periodicTable.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

TEST_CASE("Physics: AtomData", "[AtomData][Physics][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Physics: AtomData\n";

  REQUIRE(AtomData::atomicSymbol(55) == "Cs");
  REQUIRE(AtomData::atomicName(55) == "cesium");
  REQUIRE(AtomData::atomic_Z("55") == 55);
  REQUIRE(AtomData::atomic_Z(55) == 55);
  REQUIRE(AtomData::atomic_Z("Cs") == 55);

  REQUIRE(AtomData::atomicSymbol(137) == "137");
  REQUIRE(AtomData::atomicName(137) == "E137");

  REQUIRE(AtomData::defaultA(55) == 133);

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

  //===========================================

  REQUIRE(AtomData::coreConfig("[He]") == "1s2");
  REQUIRE(AtomData::coreConfig("[Hg]") ==
          AtomData::coreConfig("[Xe]") + ",4f14,5d10,6s2");
  REQUIRE(AtomData::coreConfig("InvalidNobleGas") == "InvalidNobleGas");

  for (const auto &NG :
       {"[He]", "[Ne]", "[Ar]", "[Kr]", "[Xe]", "[Rn]", "[Og]"})
    REQUIRE(AtomData::niceCoreOutput(AtomData::coreConfig(NG)) == NG);

  REQUIRE(std::abs(AtomData::diracen(1.0, 1, -1) -
                   AtomData::diracen(1.0, 1, 1)) < 1.0e-6);

  //===========================================
  {
    const auto core_string = "[Ne],7p1,6d7";
    const auto expected_symbols =
        std::vector{"1s2", "2s2", "2p6", "7p1", "6d7"};
    const auto core_list = AtomData::core_parser(core_string);
    const auto core_list2 = AtomData::state_parser("1s2,2s2,2p6,7p1,6d7");
    REQUIRE(expected_symbols.size() == core_list.size());
    REQUIRE(core_list.size() == core_list2.size());

    for (std::size_t i = 0; i < core_list.size(); ++i) {
      REQUIRE(core_list.at(i).symbol() == expected_symbols.at(i));
      const auto nlm = AtomData::term_parser(expected_symbols.at(i));
      REQUIRE(nlm == core_list.at(i));
      REQUIRE(nlm.num == core_list.at(i).num);
      REQUIRE(core_list2.at(i) == core_list.at(i));
      REQUIRE(core_list2.at(i).num == core_list.at(i).num);
    }

    REQUIRE(AtomData::guessCoreConfigStr(55) == "[Xe],6s1");
    REQUIRE(AtomData::guessCoreConfigStr(2) == "[He]");
    REQUIRE(AtomData::guessCoreConfigStr(1) == "[],1s1");
    REQUIRE(AtomData::guessCoreConfigStr(29) == "[Ar],3d10,4s1");
    REQUIRE(AtomData::guessCoreConfigStr(47) == "[Kr],4d10,5s1");
    REQUIRE(AtomData::guessCoreConfigStr(79) == "[Xe],4f14,5d10,6s1");
    REQUIRE(AtomData::guessCoreConfigStr(119) == "[Og],8s1");
  }

  //===========================================

  {
    const std::string states = "4sp3d";
    const std::vector full_list =
        std::vector{std::pair{1, -1}, {2, -1}, {3, -1}, {4, -1}, //
                    {2, 1},           {3, 1},  {4, 1},           //
                    {2, -2},          {3, -2}, {4, -2},          //
                    {3, 2},           {3, -3}};
    // Generates a list of DiracSEnken from string: full list
    const auto out_list = AtomData::listOfStates_nk(states);
    REQUIRE(out_list.size() == full_list.size());
    for (std::size_t i = 0; i < out_list.size(); ++i) {
      REQUIRE(out_list.at(i).n == full_list.at(i).first);
      REQUIRE(out_list.at(i).k == full_list.at(i).second);
    }

    const std::vector single_list =
        std::vector{std::pair{4, -1}, {4, 1}, {4, -2}, {3, 2}, {3, -3}};
    // Generates a list of DiracSEnken from string: just max n for each kappa
    const auto out_single = AtomData::listOfStates_singlen(states);
    REQUIRE(out_single.size() == single_list.size());
    for (std::size_t i = 0; i < out_single.size(); ++i) {
      REQUIRE(out_single.at(i).n == single_list.at(i).first);
      REQUIRE(out_single.at(i).k == single_list.at(i).second);
    }
  }

  //===========================================
  REQUIRE(AtomData::int_to_roman(3) == "III");
  REQUIRE(AtomData::int_to_roman(12) == "XII");
  REQUIRE(AtomData::int_to_roman(1) == "I");
  REQUIRE(AtomData::int_to_roman(4000) == "4000");

  //===========================================
  REQUIRE(AtomData::states_below_n(1) == 0);
  REQUIRE(AtomData::states_below_n(2) == 1);
  REQUIRE(AtomData::states_below_n(3) == 4);
  REQUIRE(AtomData::states_below_n(4) == 9);
  REQUIRE(AtomData::states_below_n(5) == 16);

  int index = 0;
  for (int n = 1; n < 15; ++n) {
    for (int l = 0; l < n; ++l) {
      const int k1 = l;
      const int k2 = -l - 1;
      if (l > 0) {
        REQUIRE(AtomData::nk_to_index(n, k1) == index);
        const auto [nx, kx] = AtomData::index_to_nk(index);
        REQUIRE(nx == n);
        REQUIRE(kx == k1);
        ++index;
      }
      REQUIRE(AtomData::nk_to_index(n, k2) == index);
      const auto [ny, ky] = AtomData::index_to_nk(index);
      REQUIRE(ny == n);
      REQUIRE(ky == k2);
      ++index;
    }
  }

  //===========================================
  AtomData::NonRelSEConfig a{3, 0, 2};
  AtomData::NonRelSEConfig a2{3, 0, 1};
  REQUIRE(a.n == 3);
  REQUIRE(a.l == 0);
  REQUIRE(a.num == 2);
  REQUIRE(a.ok());
  REQUIRE(a.symbol() == "3s2");
  REQUIRE(std::abs(a.frac() - 1.0) < 1.0e-6);
  REQUIRE(std::abs(a2.frac() - 0.5) < 1.0e-6);
  REQUIRE(a == a2);
  auto a3 = a2;
  a3 += a2;
  REQUIRE(std::abs(a3.frac() - 1.0) < 1.0e-6);

  AtomData::NonRelSEConfig b{3, 0, 3};
  REQUIRE(!b.ok());
  REQUIRE(b.symbol() == "3s3");

  AtomData::NonRelSEConfig c{1, 1, 1};
  REQUIRE(!c.ok());
  REQUIRE(a != c);

  AtomData::DiracSEnken nk{1, -1, -0.5};
  REQUIRE(nk.n == 1);
  REQUIRE(nk.k == -1);
  REQUIRE(nk.en == -0.5);

  //===========================================
  //===========================================
  // Just require that this causes no errors
  AtomData::printTable();
}
