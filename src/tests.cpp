// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "IO/ChronoTimer.hpp"
#include "catch2/catch.hpp"
#include "fmt/color.hpp"
#include "version/version.hpp"
#include <iostream>

struct FailureInfo {
  std::string test_name;
  std::string file;
  std::size_t line = 0;
  std::string expr;
};

class MyTestPrinter : public Catch::TestEventListenerBase {
public:
  using TestEventListenerBase::TestEventListenerBase;

  // std::vector<std::string> failed_tests{};

  std::string current_test{};
  std::vector<FailureInfo> failures{};
  std::set<std::string> already_seen{};

  // When test case starts:
  void testCaseStarting(Catch::TestCaseInfo const &testInfo) override {

    current_test = testInfo.name;

    fmt2::styled_print(
      fg(fmt::color::blue),
      "\n------------------------------------------------------------\n");
    fmt2::styled_print(fg(fmt::color::blue), testInfo.name + '\n');
    fmt2::styled_print(fg(fmt::color::blue), "Tags: ");
    std::cout << testInfo.tagsAsString() << "\n";
  }

  // When each assertion ends:
  bool assertionEnded(Catch::AssertionStats const &stats) override {
    const auto &result = stats.assertionResult;

    if (!result.isOk()) {
      // already_seen.insert(current_test);

      const auto &loc = result.getSourceInfo();
      failures.push_back({current_test, loc.file ? loc.file : "",
                          static_cast<std::size_t>(loc.line),
                          result.getExpression()});
    }

    // let Catch clear stored messages for this assertion:
    return true;
  }

  // When the whole case ends:
  void testCaseEnded(Catch::TestCaseStats const &stats) override {

    if (!stats.totals.assertions.allOk()) {
      // failed_tests.push_back(stats.testInfo.name);
      fmt2::styled_print(fg(fmt::color::red),
                         "FAILED: " + stats.testInfo.name + '\n');
    } else {
      fmt2::styled_print(fg(fmt::color::green),
                         "PASSED: " + stats.testInfo.name + '\n');
    }
  }

  void testRunEnded(Catch::TestRunStats const &) override {
    std::cout << "\n";
    if (failures.empty())
      return;

    fmt2::styled_print(
      fg(fmt::color::red),
      "============================================================\n");
    fmt2::styled_print(fg(fmt::color::red), "FAILED Tests: \n");
    std::string prev{};
    for (const auto &f : failures) {
      if (f.test_name != prev)
        std::cout << " * " << f.test_name << "\n";
      std::cout << "   - " << f.expr << "  [" << f.file << ':' << f.line
                << "]\n";
      prev = f.test_name;
    }
    std::cout << "\n";
  }
};

CATCH_REGISTER_LISTENER(MyTestPrinter)

int main(int argc, char *argv[]) {

  IO::ChronoTimer timer{"Catch2 tests"};

  std::cout << "ampsci tests (Catch2):\n";
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Libraries:\n" << version::libraries() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';

  // int result = Catch::Session().run(argc, argv);
  Catch::Session session;
  session.applyCommandLine(argc, argv);

  int result = session.run();

  std::cout << "Run 'make remove_junk' to clear junk output files\n";
  return result;
}