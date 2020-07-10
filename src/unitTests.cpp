#include "HF/ExternalField_test.hpp"
#include "HF/HartreeFock_test.hpp"
#include "HF/MixedStates_test.hpp"
#include "Maths/LinAlg_test.hpp"
#include "Physics/RadPot_test.hpp"
#include "Wavefunction/BSplineBasis_test.hpp"
//
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include <cassert>
#include <iostream>
#include <string>

// Vector of available unit tests
static const std::vector<std::pair<std::string, bool (*)(std::ostream &obuff)>>
    test_list{
        //
        {"HartreeFock", &UnitTest::HartreeFock},
        {"MixedStates", &UnitTest::MixedStates},
        {"ExternalField", &UnitTest::ExternalField},
        {"RadPot", &UnitTest::RadPot},
        {"BSplineBasis", &UnitTest::BSplineBasis},
        {"LinAlg", &UnitTest::LinAlg}
        //
    };

int main(int argc, char *argv[]) {
  IO::ChronoTimer timer("\nUnit tests");
  const std::string input_file = (argc > 1) ? argv[1] : "unitTests.in";
  IO::print_line();

  // Read in input options file
  std::cout << "Reading input from: " << input_file << "\n";
  const IO::UserInput input(input_file);

  std::ostringstream out_buff;
  out_buff << "diracSCAS test. git:" << IO::git_info() << "\n";
  out_buff << IO::time_date() << "\n";

  const auto default_tf = input.get("UnitTests", "default", true);

  // Run each unit test, count passed/failed
  int passed = 0;
  int failed = 0;
  int total = 0;
  for (const auto &[name, test] : test_list) {
    const auto run = input.get("UnitTests", name, default_tf);
    if (!run)
      continue;
    out_buff << name << "\n";
    const auto passedQ = test(out_buff);
    ++total;
    if (passedQ) {
      out_buff << "PASSED " << name << "\n";
      ++passed;
    } else {
      ++failed;
      out_buff << "FAILED " << name << " ~~~~~\n";
    }
  }
  out_buff << "\n";
  if (failed > 0) {
    out_buff << "FAILS " << failed << "/" << total << " tests; ";
  }
  out_buff << "passes " << passed << "/" << total << " tests.\n";

  // Output results:
  std::cout << out_buff.str();
  std::ofstream of(IO::date() + "_unitTests.txt");
  of << out_buff.str();

  assert(failed == 0);
  return failed;
}

//******************************************************************************
