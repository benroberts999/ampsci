#include "Angular/Angular_test.hpp"
#include "Coulomb/Coulomb_test.hpp"
#include "DiracODE/DiracODE_test.hpp"
#include "HF/ExternalField_test.hpp"
#include "HF/HartreeFock_test.hpp"
#include "HF/MixedStates_test.hpp"
#include "MBPT/CorrelationPotential_test.hpp"
#include "Maths/LinAlg_test.hpp"
#include "Physics/RadPot_test.hpp"
#include "Wavefunction/BSplineBasis_test.hpp"
//
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include <cassert>
#include <iostream>
#include <string>

/*!
@brief Units tests for various pieces of the code.

@details
Each unit test tests a certain piece of the code.

To run the unit tests, compile unitTests (make unitTests), and run from command
line. It takes an (optional) input file. By default every test is run.

Input options specified as (all are either true or false):

UnitTests{
  setting = option;
}

See example input file in doc/unitTests.in.example

Allowed settings: "default", "DiracODE", "HartreeFock", "MixedStates",
                          "ExternalField", "RadPot", "Angular", "LinAlg",
                          "BSplineBasis", "Coulomb", "CorrelationPotential"

Set default = true to run all tests except ones specified as false.
Set default = false to only run tests specifically specified as true.

*/
namespace UnitTest {

// Vector of available unit tests
static const std::vector<std::pair<std::string, bool (*)(std::ostream &obuff)>>
    test_list{
        //
        {"DiracODE", &DiracODE},
        {"HartreeFock", &HartreeFock},
        {"MixedStates", &MixedStates},
        {"ExternalField", &ExternalField},
        {"RadPot", &RadPot},
        {"Angular", &Angular},
        {"LinAlg", &LinAlg},
        {"BSplineBasis", &BSplineBasis},
        {"Coulomb", &Coulomb},
        {"CorrelationPotential", &CorrelationPotential}
        //
    };
} // namespace UnitTest

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

  input.check("UnitTests", {"default", "DiracODE", "HartreeFock", "MixedStates",
                            "ExternalField", "RadPot", "Angular", "LinAlg",
                            "BSplineBasis", "Coulomb", "CorrelationPotential"});

  const auto default_tf = input.get("UnitTests", "default", true);

  // Run each unit test, count passed/failed
  int passed = 0;
  int failed = 0;
  int total = 0;
  for (const auto &[name, test] : UnitTest::test_list) {
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
