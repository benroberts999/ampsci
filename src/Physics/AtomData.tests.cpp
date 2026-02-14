#include "AtomData.hpp"
#include "catch2/catch.hpp"
#include "periodicTable.hpp"
#include "qip/String.hpp"
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
  REQUIRE(AtomData::atomic_Z("cesium") == 55);

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

  REQUIRE(AtomData::shortSymbol(1, -1) == "1s+");
  REQUIRE(AtomData::shortSymbol(2, 1) == "2p-");
  REQUIRE(AtomData::shortSymbol(3, -2) == "3p+");
  REQUIRE(AtomData::shortSymbol(4, 2) == "4d-");
  REQUIRE(AtomData::shortSymbol(5, -3) == "5d+");
  REQUIRE(AtomData::shortSymbol(6, 3) == "6f-");
  REQUIRE(AtomData::shortSymbol(7, -4) == "7f+");
  REQUIRE(AtomData::shortSymbol(8, 4) == "8g-");

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

  auto s1 = AtomData::configs_to_string(AtomData::core_parser("[Cs],6s-1"));
  REQUIRE(AtomData::niceCoreOutput(s1) == "[Xe]");
  auto s2 = AtomData::configs_to_string(AtomData::core_parser("[Cs],6s-1,6s1"));
  REQUIRE(AtomData::niceCoreOutput(s2) == "[Xe],6s1");
  auto s3 = AtomData::configs_to_string(AtomData::core_parser("[Ho]"));
  REQUIRE(AtomData::niceCoreOutput(s3) == "[Xe],6s2,4f11");
  REQUIRE(AtomData::niceCoreOutput(
            "1s2,2s2,2p6,3s2,3p6,4s2,3d10,4p6,5s2,4d10,5p6,6s1") == "[Xe],6s1");
  REQUIRE(AtomData::niceCoreOutput(
            "1s2,2s2,2p6,3s2,3p6,4s2,3d10,4p6,5s2,4d10,5p6,6s2,4f14,"
            "5d10,6p6,6d1,7s1") == "[Rn],6d1,7s1");
  REQUIRE(AtomData::niceCoreOutput(
            "4s2,3d10,4p6,5s2,4d10,5p6,4f14,5d10,6s1,1s2,2s2,2p6,3s2,3p6") ==
          "[Xe],4f14,5d10,6s1");

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
    // Generates a list of DiracConfig from string: full list
    const auto out_list = AtomData::listOfStates_nk(states);
    REQUIRE(out_list.size() == full_list.size());
    for (std::size_t i = 0; i < out_list.size(); ++i) {
      REQUIRE(out_list.at(i).n == full_list.at(i).first);
      REQUIRE(out_list.at(i).k == full_list.at(i).second);
    }

    const std::vector single_list =
      std::vector{std::pair{4, -1}, {4, 1}, {4, -2}, {3, 2}, {3, -3}};
    // Generates a list of DiracConfig from string: just max n for each kappa
    const auto out_single = AtomData::listOfStates_singlen(states);
    REQUIRE(out_single.size() == single_list.size());
    for (std::size_t i = 0; i < out_single.size(); ++i) {
      REQUIRE(out_single.at(i).n == single_list.at(i).first);
      REQUIRE(out_single.at(i).k == single_list.at(i).second);
    }
  }

  //===========================================
  AtomData::NonRelConfig a{3, 0, 2};
  AtomData::NonRelConfig a2{3, 0, 1};
  REQUIRE(a.n == 3);
  REQUIRE(a.l == 0);
  REQUIRE(a.num == 2);
  REQUIRE(a.ok());
  REQUIRE(a.symbol() == "3s2");
  REQUIRE(a.frac() == Approx(1.0));
  REQUIRE(a2.frac() == Approx(0.5));
  REQUIRE(a == a2);
  auto a3 = a2;
  a3 += a2;
  REQUIRE(a3.frac() == Approx(1.0));

  const auto a4 = a2 + a2;
  REQUIRE(a4 == a3);
  REQUIRE(a4 == a2); // compares n,l only

  AtomData::NonRelConfig b{3, 0, 3};
  REQUIRE(!b.ok());
  REQUIRE(b.symbol() == "3s3");

  AtomData::NonRelConfig c{1, 1, 1};
  REQUIRE_FALSE(c.ok());
  REQUIRE(a != c);
  REQUIRE(a > c);
  REQUIRE(c < a);
  REQUIRE(c <= a);
  REQUIRE_FALSE(c >= a);

  AtomData::DiracConfig nk{1, -1, -0.5};
  REQUIRE(nk.n == 1);
  REQUIRE(nk.k == -1);
  REQUIRE(nk.en == -0.5);

  REQUIRE(AtomData::n_kappa_list("6sp5d") ==
          std::vector{std::pair{6, -1}, {6, 1}, {6, -2}, {5, 2}, {5, -3}});
  REQUIRE(AtomData::n_kappa_list("12s11p5d17f") ==
          std::vector{std::pair{12, -1},
                      {11, 1},
                      {11, -2},
                      {5, 2},
                      {5, -3},
                      {17, 3},
                      {17, -4}});

  //===========================================
  //===========================================
  // Just require that this causes no errors
  AtomData::printTable();
}
