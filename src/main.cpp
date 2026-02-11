
#include "DiracOperator/GenerateOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/modules_list.hpp"
#include "Physics/periodicTable.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "ampsci/ampsci.hpp"
#include "fmt/color.hpp"
#include "qip/String.hpp"
#include "qip/omp.hpp"
#include "version/EasterEgg.hpp"
#include "version/version.hpp"
#include <iostream>
#include <string>

//! Man page info
namespace man {

const std::string name{"ampsci - Atomic Many-body Perturbation theory in the "
                       "Screened Coulomb Interaction."};

const std::string author{"Benjamin M. Roberts (https://broberts.io/), "
                         "University of Queensland, Australia."};

const std::string synopsis{"ampsci [InputFile]\n"
                           "ampsci [Atom] [Core] [Valence]\n"
                           "ampsci -a [InputBlock]\n"
                           "ampsci -m [Module]\n"
                           "ampsci -o [Operator]\n"
                           "ampsci -p [Atom] [Isotope]\n"
                           "ampsci -c [value] [unit]\n"};

const std::string description{
    "ampsci is a C++ program for high-precision atomic structure calculations "
    "of single-valence systems.\n"
    "Full documentation (including of the physics methods) can be found "
    "online: https://ampsci.dev/"};

const std::vector<std::pair<std::string, std::string>> options{
    {"[InputFile]",
     "Runs ampsci taking options specified in file 'InputFile'. "
     "See documentation (or option -a) for input file format.\n"
     "Example:\n"
     "./ampsci input.in\n"
     "    - Runs ampsci taking input options from file 'input.in'"},
    {"[Atom] [Core] [Valence]",
     "For quick and simple HF calculation. If core is not given, guesses "
     "core configuration and runs using V^N approximation. Input strings use "
     "same format as input file.\n"
     "Examples:\n"
     "./ampsci Cs\n"
     "    - Runs ampsci for Cs using Hartree Fock (V^N) approximation\n"
     "./ampsci Cs [Xe] 6sd5d\n"
     "    - Runs ampsci for Cs using Hartree Fock with Xe-like core and "
     "valence states up to n=6 for s,p-states and n=5 for d-states\n"},
    {"-a [BlockName] ..., --ampsci",
     "Prints list of available top-level ampsci options. BlockName is "
     "optional; if given it will print options for given ampsci Block. You may "
     "list any number of blocks (space separated).\nExamples:\n"
     "./ampsci -a\n"
     "    - Prints list of all available top-level ampsci options\n"
     "./ampsci -a Atom HartreeFock\n"
     "    - Prints list of all available options in the 'Atom' and "
     "'HartreeFock' input blocks"},
    {"-c [value] [unit], --constants",
     "Prints some handy physical constants, and does unit conversions.\n"
     "Converts energies between: au (=Eh), eV, cm^-1, Hz, nm (wavelength)\n"
     "Converts lengths between: a0, fm, eV^-1\n"
     "Standard SI prefixes allowed: fempto (f) to Terra (T); use `u` for "
     "micro.\n"
     "Examples:\n"
     "./ampsci -c 123.456 eV\n"
     "    - Converts 123.456 eV to other energy units\n"
     "./ampsci -c 27000 cm^-1\n"
     "    - Converts 27000 cm^-1 to other energy units\n"
     "./ampsci -c 4.561 fm\n"
     "    - Converts 4.561 fm to other length units\n"
     "./ampsci -c 1.05 um\n"
     "    - Converts 1.05 micro-m to other length units\n"
     "./ampsci -c 0.01 MeV^-1\n"
     "    - Converts 0.01 MeV^-1 to other length units\n"},
    {"-h, --help, -?",
     "Prints help info, including some detail on input options"},
    {"-l, --libs, --libraries",
     "Prints version details for libaries used by ampsci"},
    {"-m [ModuleName], --modules",
     "Prints list of available Modules. ModuleName is optional; if given, "
     "will list avaiable options for that Module.\n"
     "Examples:\n"
     "./ampsci -m\n"
     "    - Prints list of available modules\n"
     "./ampsci -m MatrixElements\n"
     "    - Prints list of input options for the 'MatrixElements' module\n"},
    {"-o [OperatorName], --operators",
     "Prints list of available operators (for calculating matrix elements). "
     "OperatorName is optional; if given, will list avaiable options for that "
     "operator (most operators take no options).\n"
     "Examples:\n"
     "./ampsci -o\n"
     "    - Prints list of available operators\n"
     "./ampsci -o E1\n"
     "    - Prints list of input optins for the 'E1' operator\n"},
    {"-p [Atom] [Isotope], --periodicTable",
     "Prints textual periodic table with electronic + nuclear information. "
     "Atom and Isotope are optional; if given, will print info for that "
     "isotope. Atom should be atomic symbol (eg Cs), or Z (55). If Isotope is "
     "blank, will print for 'default' isotope. Can also list 'all' known "
     "isotope info\n"
     "Examples:\n"
     "./ampsci -p\n"
     "    - Prints a text periodic table\n"
     "./ampsci -p Cs\n"
     "    - Prints a text periodic table, with atomic information for default "
     "Cs nucleus\n"
     "./ampsci -p Cs 131\n"
     "    - Prints a text periodic table, with atomic information for Cs-131\n"
     "./ampsci -p Cs all\n"
     "    - Prints a text periodic table, with atomic information for all "
     "(available) Cs iotopes\n"},
    {"-s [input options as string], --string",
     "Takes input options as a string, using same format as input file\n"},
    {"-v, --version", "Prints ampsci version (and git commit) details"} //
};

//! Prints 'man page' style info
void print_manual() {
  const int wrap_at = 80;
  std::string tab = "    ";
  fmt2::styled_print(fmt::emphasis::bold, "NAME\n");
  fmt::print(qip::wrap(name, wrap_at, tab + tab));

  std::cout << "\n";
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Libraries:\n" << version::libraries() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';

  std::cout << "\n\n";

  fmt2::styled_print(fmt::emphasis::bold, "SYNOPSIS\n");
  fmt::print(qip::wrap(synopsis, wrap_at, tab + tab));

  std::cout << "\n\n";

  fmt2::styled_print(fmt::emphasis::bold, "DESCRIPTION\n");
  fmt::print(qip::wrap(description, wrap_at, tab + tab));

  std::cout << "\n\n";

  fmt2::styled_print(fmt::emphasis::bold, "OPTIONS\n");
  for (const auto &[option, text] : options) {
    fmt2::styled_print(fg(fmt::color::steel_blue),
                       qip::wrap(option, wrap_at, tab + tab));
    std::cout << "\n";
    fmt::print(qip::wrap(text, wrap_at, tab + tab + tab));
    std::cout << "\n\n";
  }

  fmt2::styled_print(fmt::emphasis::bold, "AUTHOR\n");
  fmt::print(qip::wrap(author, wrap_at, tab + tab));

  std::cout << "\n";
}

} // namespace man

//==============================================================================

//==============================================================================
//==============================================================================
//! Parses command-line input, then runs ampsci
int main(int argc, char *argv[]) {
  using namespace std::string_literals;

  // Parse input text into strings:
  const std::string in_text_1 = (argc > 1) ? argv[1] : "";
  const std::string in_text_2 = (argc > 2) ? argv[2] : "";
  const std::string in_text_3 = (argc > 3) ? argv[3] : "";

  // check for special commands
  if (in_text_1 == "") {
    man::print_manual();
    return 0;
  } else if (in_text_1 == "-v" || in_text_1 == "--version") {
    std::cout << "AMPSCI v: " << version::version() << '\n';
    std::cout << "Libraries:\n" << version::libraries() << '\n';
    std::cout << "Compiled: " << version::compiled() << '\n';
    std::cout << man::author << "\n";
    return 0;
  } else if (in_text_1 == "-l" || in_text_1.substr(0, 5) == "--lib") {
    std::cout << "Libraries:\n" << version::libraries() << '\n';
    return 0;
  } else if (in_text_1 == "-h" || in_text_1 == "--help" || in_text_1 == "-?") {
    man::print_manual();
    return 0;
  } else if (in_text_1 == "-m" || in_text_1 == "--modules") {
    std::cout << "Available modules: \n";
    Module::list_modules();
    const std::string module_name = (argc > 2) ? argv[2] : "";
    if (!module_name.empty()) {
      // run the module, with option 'help' set. This will trigger the helper
      // to print the details for the available options in that module
      Module::runModule(IO::InputBlock{"Module::"s + module_name, {"help;"}},
                        {});
    }
    return 0;
  } else if (in_text_1 == "-o" || in_text_1 == "--operators") {
    std::cout << "Available operators: \n";
    DiracOperator::list_operators();
    const std::string op_name = (argc > 2) ? argv[2] : "";
    if (!op_name.empty()) {
      Wavefunction wf{{1, 1.0, 1.0}, {1, 1}};
      DiracOperator::generate(op_name, IO::InputBlock{op_name, {"help;"}}, wf);
    }
    return 0;
  } else if (in_text_1 == "-a" || in_text_1 == "--ampsci") {
    auto temp_input = IO::InputBlock{"ampsci", {"help;"}};
    for (int i_in = 2; i_in < argc; ++i_in) {
      const std::string block_name = (argc > i_in) ? argv[i_in] : "";
      if (!block_name.empty()) {
        temp_input.add(IO::InputBlock{block_name, {"help;"}});
      }
    }
    ampsci(temp_input);
    return 0;
  } else if (in_text_1 == "-p" || in_text_1 == "--periodicTable") {
    std::string z_str = (argc > 2) ? argv[2] : "";
    std::string a_str = (argc > 3) ? argv[3] : "";
    AtomData::periodicTable(z_str, a_str);
    return 0;
  } else if (in_text_1 == "-e" || in_text_1 == "--EasterEgg") {
    std::cout << EasterEgg::get_egg();
    return 0;
  } else if (in_text_1 == "-c" || in_text_1 == "--constants") {
    AtomData::printConstants();
    if (!in_text_2.empty() && !in_text_3.empty())
      AtomData::conversions(std::stod(in_text_2), in_text_3);
    return 0;
  } else if (in_text_1 == "-s" || in_text_1 == "--string") {
    std::cout << "Reading input from command-line string\n";
  } else if (!in_text_1.empty() && in_text_1.front() == '-') {
    std::cout << "Unrecognised option: " << in_text_1 << '\n';
    man::print_manual();
    return 0;
  }

  // Print git/version info to screen:
  std::cout << '\n';
  IO::print_line();
  std::cout << "AMPSCI v: " << version::version() << '\n';
  std::cout << "Parallel: " << qip::omp_details() << '\n';
  std::cout << "Compiled: " << version::compiled() << '\n';
  std::cout << "Run time: " << IO::time_date() << '\n';

  // If we are not given a valid input text file, assume input is in form:
  // <At> <core> <valence> (e.g., "Cs [Xe] 6sp")
  // All optional, but must appear in order. ("Core" may be skipped for H)
  const auto atom = AtomData::atomicSymbol(AtomData::atomic_Z(in_text_1));
  const auto core =
      atom == "H" ? "" : (in_text_2 == "" ? "[" + atom + "]" : in_text_2);
  // allow core to be skipped for Hydrogen
  const auto valence = (atom == "H" && argc == 3) ? in_text_2 : in_text_3;
  std::string default_input = "Atom{Z=" + atom + ";}" +
                              "HartreeFock { core = " + core +
                              "; valence = " + valence + ";}";

  if (in_text_1 == "-s" || in_text_1 == "--string") {
    default_input = in_text_2;
  }

  // nb: std::filesystem not available in g++-7 (getafix version)
  const auto fstream = std::fstream(in_text_1);
  const auto input = fstream.good() ? IO::InputBlock("ampsci", fstream) :
                                      IO::InputBlock("ampsci", default_input);

  // Run program. Add option to run multiple times
  ampsci(input);
}
