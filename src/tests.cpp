// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "IO/ChronoTimer.hpp"
#include "catch2/catch.hpp"
#include "fmt/color.hpp"
#include "version/version.hpp"
#include <iostream>

class MyTestPrinter : public Catch::TestEventListenerBase {
public:
  using TestEventListenerBase::TestEventListenerBase;

  void testCaseStarting(Catch::TestCaseInfo const &testInfo) override {

    fmt2::styled_print(fg(fmt::color::blue), "\n=============================="
                                             "============================="
                                             "====================\n");
    fmt2::styled_print(fg(fmt::color::blue), testInfo.name + '\n');
    fmt2::styled_print(fg(fmt::color::blue), "Tags: ");
    std::cout << testInfo.tagsAsString() << '\n';
    fmt2::styled_print(fg(fmt::color::blue), "------------------------------"
                                             "-----------------------------"
                                             "--------------------\n");
  }
};

CATCH_REGISTER_LISTENER(MyTestPrinter)

int main(int argc, char *argv[]) {

  IO::ChronoTimer timer{"Catch2 tests"};

  std::cout << "ampsci tests (Catch2):\n";
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Libraries:\n" << version::libraries() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';

  int result = Catch::Session().run(argc, argv);

  return result;
}