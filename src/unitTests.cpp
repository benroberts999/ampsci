#include "Angular/Angular_test.hpp"
#include "Angular/SixJTable_test.hpp"
#include "Coulomb/Coulomb_test.hpp"
#include "Coulomb/QkTable_test.hpp"
#include "DiracODE/DiracODE_test.hpp"
#include "ExternalField/DiagramRPA_test.hpp"
#include "ExternalField/MixedStates_test.hpp"
#include "ExternalField/TDHF_test.hpp"
#include "ExternalField/TDHFbasis_breit_test.hpp"
#include "ExternalField/TDHFbasis_test.hpp"
#include "HF/Breit_test.hpp"
#include "HF/HartreeFock_test.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp" // for time+date
#include "LinAlg/LinAlg_test.hpp"
#include "MBPT/CorrelationPotential_test.hpp"
#include "MBPT/StructureRad_test.hpp"
#include "Maths/BSpline_test.hpp"
#include "Physics/RadPot_test.hpp"
#include "Wavefunction/BSplineBasis_test.hpp"
#include "git.info"
#include <cassert>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

/*!
@brief Units tests for various pieces of the code.

@details
Each unit test tests a certain piece of the code.

To run the unit tests, compile unitTests (make unitTests), and run from command
line. It takes optional command-line options (the names of which tests to run).
By default (if no extra arguments given), every test is run.
If a test that doesn't exist is requested (or mispelled), test will fail, and
the list of available tests will be printed.

$./unitTests help
  * prints list of available tests

 e.g.:
 $./unitTests DiracODE Angular LinAlg
   * Will run the DiracODE, Angular, and LinAlg tests (in that order)

 $./unitTests
   * Will run all available tests

When a nw test is written, for it to be run:
 1) The hpp file for that test must be included into unitTests.cpp
 2) The name + function must be added to the test_list vector in unitTests.cpp

*/
namespace UnitTest {

// Vector of available unit tests
static const std::vector<std::pair<std::string, bool (*)(std::ostream &obuff)>>
    test_list{
        //
        {"DiracODE", &DiracODE},
        {"HartreeFock", &HartreeFock},
        {"Breit", &Breit},
        {"MixedStates", &MixedStates},
        {"TDHF", &TDHF},
        {"TDHFbasis", &TDHFbasis},
        {"TDHFbasis_breit", &TDHFbasis_breit},
        {"RadPot", &RadPot},
        {"Angular", &Angular},
        {"SixJTable", &SixJTable},
        {"BSplineBasis", &BSplineBasis},
        {"Coulomb", &Coulomb},
        {"MBPT2", &MBPT2},
        {"Sigma2", &Sigma2},
        {"SigmaAO", &SigmaAO},
        {"DiagramRPA", &DiagramRPA},
        {"QkTable", &QkTable},
        {"StructureRad", &StructureRad},
        {"BSplineClass", &BSplineClass},
        {"LinAlg", &LinAlg}
        //
    };

static const std::vector<std::pair<std::string, bool (*)(std::ostream &obuff)>>
    quick_list{{"DiracODE", &DiracODE},
               {"HartreeFock", &HartreeFock},
               {"MixedStates", &MixedStates},
               {"Angular", &Angular},
               {"SixJTable", &SixJTable},
               {"LinAlg", &LinAlg},
               {"BSplineBasis", &BSplineBasis},
               {"Coulomb", &Coulomb},
               {"MBPT2", &MBPT2}};

//------------------------------------------------------------------------------
// Looks up test + returns its function. If test not in list, prints list to
// help users know which are available.
// Note: if ever a bad test name is asked for - unit test will fail!
// We don't want a test to spuriously pass!
auto get_test(std::string_view in_name) {

  for (const auto &[name, test] : UnitTest::test_list) {
    if (name == in_name)
      return test;
  }

  if (in_name != "help") {
    std::cout << "No test named: " << in_name << "\n";
  }
  std::cout << "Available tests:\n";
  for (const auto &[name, test] : UnitTest::test_list) {
    std::cout << name << "\n";
  }
  // nb: we do not want to report a working test if we asked for the list!
  assert(false);
  std::abort(); // so won't fail iff compiled without assertions
}

} // namespace UnitTest

//******************************************************************************
int main(int argc, char *argv[]) {

  std::ostringstream out_buff;
  out_buff << "ampsci test. git:" << GitInfo::gitversion << " ("
           << GitInfo::gitbranch << ")\n";
  out_buff << IO::time_date() << "\n";

  // Make a list of all tests to Run
  // (do this to allow possibility of running all tests)
  std::vector<std::string_view> name_list;
  if (argc <= 1) {
    // run all tests
    for (const auto &[name, test] : UnitTest::test_list)
      name_list.emplace_back(name);
  } else if (argv[1] == std::string("quick")) {
    for (const auto &[name, test] : UnitTest::quick_list)
      name_list.emplace_back(name);
  } else {
    // run only those tests specifically asked
    for (int i = 1; i < argc; ++i)
      name_list.emplace_back(argv[i]);
  }

  // Run each unit test, count passed/failed
  int passed = 0;
  int failed = 0;
  int total = 0;
  {
    IO::ChronoTimer timer("\nUnit tests");

    for (const auto &name : name_list) {
      std::cout << "\nRunning: " << name << "\n";
      const auto test = UnitTest::get_test(name);
      out_buff << name << "\n";
      IO::ChronoTimer t(name);
      const auto passedQ = test(out_buff);
      ++total;
      if (passedQ) {
        out_buff << "PASSED " << name << "\n";
        std::cout << "PASSED " << name << "\n";
        ++passed;
      } else {
        ++failed;
        out_buff << "FAILED " << name << " ~~~~~\n";
        std::cout << "FAILED " << name << " ~~~~~\n";
      }
      out_buff << "T= " << t.reading_str() << "\n";
    }
    out_buff << "\n";
    if (failed > 0) {
      out_buff << "FAILS " << failed << "/" << total << " tests; ";
    }
    out_buff << "passes " << passed << "/" << total << " tests.\n";

    // Output results:
    std::cout << "\n" << out_buff.str();
    std::ofstream of("unitTests_" + IO::date() + "_" + IO::time() + ".txt");
    of << out_buff.str();
  }

  assert(failed == 0);
  // do not allow a test to seemingly "pass" just because it wasn't run!
  assert(passed != 0);
  return failed;
}

//******************************************************************************
