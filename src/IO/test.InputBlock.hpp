#pragma once
#include "InputBlock.hpp"
#include <cassert>
#include <iostream>

inline void run_tests(const IO::InputBlock &ib);

inline void test_InputBlock() {
  // A basic unit test  of IO::InputBlock
  // Tests each function, but not all that thoroughly

  using namespace IO;

  InputBlock ib("name1", {{"k1", "1"}, {"k2", "2.5"}});
  ib.add(Option{"k3", "number_3"});
  ib.add(InputBlock("blockA", {{"kA1", "old_val"}, {"kA1", "new_val"}}));
  ib.add(InputBlock("blockB", {{"keyB1", "valB"}}));
  // 'true' means will be merged with existing block (if exists)
  ib.add(InputBlock("blockB", {{"keyB2", "17.3"}}), true);
  std::string input_string = "blockC{ kC1=1; InnerBlock{ kib1=-6; } }";
  ib.add(input_string);
  ib.add(Option{"list", "1,2,3,4,5"});
  ib.add(Option{"bool1", "true"});
  ib.add(Option{"bool2", "false"});

  run_tests(ib);

  // Construct a new InputBlock using the string output of another
  std::stringstream ostr1;
  ib.print(ostr1);
  InputBlock ib2("name2", ostr1.str());

  // Run same tests on this one
  run_tests(ib2);

  // test copy construct
  auto ib3 = ib2;
  run_tests(ib3);

  // test copy assign
  ib2 = ib;
  run_tests(ib2);

  // Test that the two string outputs are identical
  std::stringstream ostr2;
  ib2.print(ostr2);
  assert(ostr1.str() == ostr2.str());

  std::cout << "\nPassed all tests :)\n";
}

//==============================================================================
void run_tests(const IO::InputBlock &ib) {

  assert(ib.get("k1", 0) == 1);

  // There is no k109 option, so should return default value
  assert(ib.get("k109", -17.6) == -17.6);
  // Test it returns 'empty' optional
  assert(ib.get("k109") == std::nullopt);
  assert(!ib.get("k109"));

  assert(ib.get<double>("k2") == 2.5);
  // Returns an std::optional, so test that
  assert(ib.get<double>("k2").value() == 2.5);
  // Should instatiate as std::string by default
  assert(ib.get("k3") == "number_3");

  // Should return the later of two options with same key (overrides earlier)
  assert(ib.getBlock("blockA")->get("kA1") == "new_val");
  assert(ib.getBlock("blockZ") == std::nullopt);

  // Adding two 'blockB' - these should be 'consolidated'
  assert(ib.getBlock("blockB")->get("keyB1") == "valB");
  assert(ib.getBlock("blockB")->get("keyB2", 0.0) == 17.3);

  // Test blockC - added via a string
  assert(ib.getBlock("blockC")->get<int>("kC1") == 1);
  // test nested blocks
  assert(ib.getBlock("blockC")->getBlock("InnerBlock")->get<int>("kib1") == -6);

  // Test the 'nested block' retrival
  assert(ib.get<int>({}, "k1") == 1);
  assert(ib.get<std::string>({"blockA"}, "kA1") == "new_val");
  assert(ib.get<int>({"blockC", "InnerBlock"}, "kib1") == -6);

  // test the input list
  const auto in_list = ib.get<std::vector<int>>("list", {});
  const std::vector<int> expected_list{1, 2, 3, 4, 5};
  assert(in_list.size() == expected_list.size() &&
         std::equal(in_list.cbegin(), in_list.cend(), expected_list.cbegin()));

  assert(ib.get<bool>("bool1").value() == true);
  assert(ib.get<bool>("bool2").value() == false);
}
