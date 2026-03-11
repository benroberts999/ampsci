#define CATCH_CONFIG_RUNNER
#include "IO/ChronoTimer.hpp"
#include "catch2/catch.hpp"
#include "fmt/color.hpp"
#include "version/version.hpp"
#include <chrono>
#include <iostream>
#include <limits>

const auto test_colour = fmt::color::blue;
const auto pass_colour = fmt::color::green;
const auto fail_colour = fmt::color::red;
const auto time_colour = fmt::color::dark_slate_gray;

using steady_clock = std::chrono::steady_clock;

// Store information about each failed test, to summarise at end
struct FailureInfo {
  std::string test_name{};
  std::string file{};
  std::size_t line = 0;
  std::string expr{};
};

// Custom "listner" to print info for each test run,
// and to store overall summary details (failed tests, timing info)
class AmpsciTestListner : public Catch::TestEventListenerBase {
public:
  using TestEventListenerBase::TestEventListenerBase;

  std::string current_test{};
  std::vector<FailureInfo> failures{};
  std::set<std::size_t> already_seen{};

  steady_clock::time_point test_start{};
  double slowest_test_time = 0.0;
  std::string slowest_test_name{};
  std::string slowest_test_location{};

  // When test case starts:
  void testCaseStarting(Catch::TestCaseInfo const &testInfo) override {

    current_test = testInfo.name;
    test_start = steady_clock::now();

    fmt2::styled_print(
      fg(test_colour),
      "\n------------------------------------------------------------\n");
    fmt2::styled_print(fg(test_colour), testInfo.name + '\n');
    fmt2::styled_print(fg(test_colour), "Tags: ");
    std::cout << testInfo.tagsAsString() << "\n";
    // this makes it easy to 'click to open in editor'
    fmt2::styled_print(fg(test_colour), "File: ");
    std::cout << testInfo.lineInfo.file << ":"
              << std::to_string(testInfo.lineInfo.line) << "\n";
  }

  // When each assertion ends:
  bool assertionEnded(Catch::AssertionStats const &stats) override {
    const auto &result = stats.assertionResult;
    const auto &loc = result.getSourceInfo();
    const auto &expr = result.getExpression();

    // Store list of failed assertions, to print summary at end
    // (Only unique assertions, hence the set)
    if (!result.isOk() && already_seen.count(loc.line) == 0) {
      already_seen.insert(loc.line);

      failures.push_back(
        {current_test, loc.file ? loc.file : "", loc.line, expr});
    }

    // let Catch clear stored messages for this assertion:
    return true;
  }

  // When the whole case ends:
  void testCaseEnded(Catch::TestCaseStats const &stats) override {

    const auto test_end = steady_clock::now();
    const auto elapsed =
      std::chrono::duration<double>(test_end - test_start).count();
    fmt2::styled_print(fg(time_colour), "Time: {:.3f} s\n", elapsed);

    if (!stats.totals.assertions.allOk()) {
      // failed_tests.push_back(stats.testInfo.name);
      fmt2::styled_print(fg(fail_colour),
                         "FAILED: " + stats.testInfo.name + '\n');
    } else {
      fmt2::styled_print(fg(pass_colour),
                         "PASSED: " + stats.testInfo.name + '\n');
    }
    already_seen.clear();

    if (elapsed > slowest_test_time) {
      slowest_test_time = elapsed;
      slowest_test_name = stats.testInfo.name;
      slowest_test_location = std::string{stats.testInfo.lineInfo.file} + ":" +
                              std::to_string(stats.testInfo.lineInfo.line);
    }
  }

  void testRunEnded(Catch::TestRunStats const &) override {
    std::cout << "\n";

    if (!failures.empty()) {
      fmt2::styled_print(
        fg(fail_colour),
        "============================================================\n");
      fmt2::styled_print(fg(fail_colour), "FAILED Tests: \n");
      std::string prev{};
      for (const auto &f : failures) {
        if (f.test_name != prev) {
          std::cout << "\n * " << f.test_name << "\n";
        }
        std::cout << "   - " << f.expr << "  [" << f.file << ':' << f.line
                  << "]\n";
        prev = f.test_name;
      }
      std::cout << "\n";
    }

    if (!slowest_test_name.empty()) {
      fmt2::styled_print(fg(time_colour), "Slowest: {} ({:.3f} s)\n",
                         slowest_test_name, slowest_test_time);
      fmt2::styled_print(fg(time_colour), slowest_test_location + "\n\n");
    }
  }
};

CATCH_REGISTER_LISTENER(AmpsciTestListner)

int main(int argc, char *argv[]) {

  IO::ChronoTimer timer{"ampsci tests"};

  std::cout << "ampsci tests:\n";
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Libraries:\n" << version::libraries() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';

  int result = Catch::Session().run(argc, argv);
  // Catch::Session session;
  // session.applyCommandLine(argc, argv);
  // int result = session.run();

  std::cout << "nb: Run 'make remove_junk' to clear junk output files\n";
  return result;
}