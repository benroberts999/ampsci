#include "InputBlock2.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <sstream>

TEST_CASE("InputBlock2: basic construction and value access", "[InputBlock2][unit]") {
  using namespace IO;

  nlohmann::json j = {
    {"k_int",    1},
    {"k_double", 2.5},
    {"k_str",    "hello"},
    {"k_bool",   true},
    {"k_empty",  ""},
    {"k_list",   nlohmann::json::array({1, 2, 3})},
    {"k_csv",    "4,5,6"},
    {"SubBlock", {{"inner", 42}}}
  };

  InputBlock2 ib{"root", j};

  REQUIRE(ib.get("k_int",    0)   == 1);
  REQUIRE(ib.get("k_double", 0.0) == Approx(2.5));
  REQUIRE(ib.get("k_str",    std::string{}) == "hello");
  REQUIRE(ib.get("k_bool",   false) == true);

  // missing key -> default
  REQUIRE(ib.get("missing", 99) == 99);

  // empty string -> nullopt
  REQUIRE_FALSE(ib.get<std::string>("k_empty").has_value());

  // optional form
  REQUIRE(ib.get<int>("k_int") == std::optional<int>{1});
  REQUIRE_FALSE(ib.get<int>("missing").has_value());
}

TEST_CASE("InputBlock2: bool coercion", "[InputBlock2][unit]") {
  using namespace IO;

  nlohmann::json j = {
    {"b_native",  true},
    {"b_str_t",   "true"},
    {"b_str_yes", "yes"},
    {"b_str_y",   "Y"},
    {"b_str_f",   "false"},
    {"b_num",     42}   // number is NOT a bool in JSON -> returns default
  };
  InputBlock2 ib{"root", j};

  REQUIRE(ib.get("b_native",  false) == true);
  REQUIRE(ib.get("b_str_t",   false) == true);
  REQUIRE(ib.get("b_str_yes", false) == true);
  REQUIRE(ib.get("b_str_y",   false) == true);
  REQUIRE(ib.get("b_str_f",   true)  == false);
  REQUIRE(ib.get("b_num",     false) == false);  // int is not bool
}

TEST_CASE("InputBlock2: vector and array", "[InputBlock2][unit]") {
  using namespace IO;

  nlohmann::json j = {
    {"arr_int", nlohmann::json::array({10, 20, 30})},
    {"arr_dbl", nlohmann::json::array({1.1, 2.2, 3.3})},
    {"arr_str", nlohmann::json::array({"a", "b", "c"})},
    {"not_arr", 42}
  };
  InputBlock2 ib{"root", j};

  auto v1 = ib.get<std::vector<int>>("arr_int", {});
  REQUIRE(v1 == std::vector<int>{10, 20, 30});

  auto a1 = ib.get<std::array<double, 3>>("arr_dbl", {});
  REQUIRE(a1[0] == Approx(1.1));
  REQUIRE(a1[1] == Approx(2.2));
  REQUIRE(a1[2] == Approx(3.3));

  auto vs = ib.get<std::vector<std::string>>("arr_str", {});
  REQUIRE(vs == std::vector<std::string>{"a", "b", "c"});

  // non-array -> default
  auto v2 = ib.get<std::vector<int>>("not_arr", {-1});
  REQUIRE(v2 == std::vector<int>{-1});
}

TEST_CASE("InputBlock2: sub-block access", "[InputBlock2][unit]") {
  using namespace IO;

  nlohmann::json j = {
    {"Sub", {{"x", 7}, {"y", 8.0}}}
  };
  InputBlock2 ib{"root", j};

  REQUIRE(ib.has_block("Sub"));
  REQUIRE_FALSE(ib.has_block("Other"));

  auto sub = ib.getBlock("Sub");
  REQUIRE(sub.has_value());
  REQUIRE(sub->get("x", 0) == 7);
  REQUIRE(sub->get("y", 0.0) == Approx(8.0));

  auto empty = ib.get_block("Other");
  REQUIRE(empty.get("x", -1) == -1);
}

TEST_CASE("InputBlock2: has_option and set", "[InputBlock2][unit]") {
  using namespace IO;

  InputBlock2 ib{"block"};

  REQUIRE_FALSE(ib.has_option("key"));
  ib.set("key", 42);
  REQUIRE(ib.has_option("key"));
  REQUIRE(ib.get("key", 0) == 42);

  // Overwrite
  ib.set("key", 99);
  REQUIRE(ib.get("key", 0) == 99);

  // set_block
  InputBlock2 child{"child"};
  child.set("v", 3.14);
  ib.set_block("child", child);
  REQUIRE(ib.has_block("child"));
  REQUIRE(ib.getBlock("child")->get("v", 0.0) == Approx(3.14));
}

TEST_CASE("InputBlock2: check validation", "[InputBlock2][unit]") {
  using namespace IO;

  nlohmann::json j = {
    {"allowed_opt", 1},
    {"Sub",         {{"x", 1}}}
  };
  InputBlock2 ib{"root", j};

  // All known -> returns true
  REQUIRE(ib.check({{"allowed_opt", "an option"},
                    {"Sub{}",       "a sub-block"},
                    {"extra",       "not present, still allowed"}}));

  // Unknown option -> returns false (warning printed)
  nlohmann::json j2 = {{"unknown_key", 5}};
  InputBlock2 ib2{"root", j2};
  std::cout << "Note: following warning is expected as part of test:\n";
  REQUIRE_FALSE(ib2.check({{"known_key", "desc"}}));
}

TEST_CASE("InputBlock2: from_file (parse JSON with comments)", "[InputBlock2][unit]") {
  using namespace IO;

  // Write a temp file with comments
  const std::string tmppath = "/tmp/ib2_test.jsonc";
  {
    std::ofstream f(tmppath);
    f << R"({
  // a comment
  "Z": "Cs",   /* block comment */
  "num": 42
})";
  }

  auto ib = InputBlock2::from_file(tmppath);
  REQUIRE(ib.get("Z",   std::string{}) == "Cs");
  REQUIRE(ib.get("num", 0)             == 42);
}
