#include "Angular/Wigner369j.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/*
Quick routine that outputs numeric values for Wigner 3,6,9J symbols.
*/

//==============================================================================
void instructions(const std::string &input = "") {
  std::cout << "\nI don't understand the input:\n";
  std::cout << "  '" << input << "'\n";
  std::cout << "Enter 3,6,8j/CG symbols like:\n\n";
  std::cout << "  3j:  \'(j1,j2,j3,m1,m2,m3)\',\n";
  std::cout << "  6j:  \'{j1,j2,j3,j4,j5,j6}\',\n";
  std::cout << "  9j:  \'{j1,j2,j3,j4,j5,j6,j7,j8,j9}\',\n or \n";
  std::cout << "  Clebschscgh Gordon:   \'<j1,m1, j2,m2 | J,M>\'\n\n";
  std::cout
      << "(If inputting from command line, must include the quote marks!)\n"
      << "(You can list symbols in an input file, run like: "
      << "./wigner -f inputfile.in [-f optional if '.in' extension is used])\n"
      << "You can use any symbol in place of the commas (incl space).\n\n";
}

//==============================================================================
void calculateSingleTerm(const std::string &input) {

  // read each number into vector:
  std::istringstream stream(input);
  std::vector<double> v;
  while (true) {
    double j;
    if (stream >> j) {
      v.push_back(j);
    } else if (stream.eof()) {
      break; // We're done
    } else {
      // We've hit a non-double
      stream.clear();  // clear the error flags
      stream.ignore(); // ignore everything till the delimeter
    }
  }

  if (input.substr(0, 1) == "!" || input.substr(0, 1) == "#")
    return;

  if (input.substr(0, 1) == "(") {
    if (v.size() != 6) {
      instructions(input);
      return;
    }
    double cgc = Angular::threej(v[0], v[1], v[2], v[3], v[4], v[5]);
    std::cout << "3j symbol:\n";
    printf("(%4.1f%4.1f%4.1f)\n", v[0], v[1], v[2]);
    printf("(%4.1f%4.1f%4.1f) = %.6f\n", v[3], v[4], v[5], cgc);
  } else if (input.substr(0, 1) == "<") {
    if (v.size() != 6) {
      instructions(input);
      return;
    }
    std::cout << "Clebschscgh Gordon coeficient:\n";
    double cgc = Angular::cg(v[0], v[1], v[2], v[3], v[4], v[5]);
    printf("<%.1f %.1f, %.1f %.1f| %.1f %.1f> = %f\n", v[0], v[1], v[2], v[3],
           v[4], v[5], cgc);
  } else if (input.substr(0, 1) == "{") {
    if (v.size() == 6) {
      std::cout << "6j symbol:\n";
      double cgc = Angular::sixj(v[0], v[1], v[2], v[3], v[4], v[5]);
      printf("{%4.1f%4.1f%4.1f}\n", v[0], v[1], v[2]);
      printf("{%4.1f%4.1f%4.1f} = %.6f\n", v[3], v[4], v[5], cgc);
    } else if (v.size() == 9) {
      std::cout << "9j symbol:\n";
      double cgc =
          Angular::ninej(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
      printf("{%4.1f%4.1f%4.1f}\n", v[0], v[1], v[2]);
      printf("{%4.1f%4.1f%4.1f} = %.6f\n", v[3], v[4], v[5], cgc);
      printf("{%4.1f%4.1f%4.1f}\n", v[6], v[7], v[8]);
    } else {
      instructions(input);
      return;
    }
  } else {
    instructions(input);
    return;
  }
}

//==============================================================================
int main(int num_in, char *argv[]) {

  if (num_in <= 1) {
    instructions();
    return 1;
  }

  auto first_arg = (std::string)argv[1];
  auto f_ext =
      (first_arg.size() > 3) ? first_arg.substr(first_arg.size() - 3) : "";

  std::vector<std::string> input_list;
  if (first_arg == "-f" || f_ext == ".in") {
    auto fname = (first_arg == "-f" && num_in > 2) ? argv[2] : argv[1];
    std::cout << "Reading from input file: " << fname << "\n\n";
    input_list = IO::FRW::readInputFile_byLine(fname);
  } else {
    for (int i = 1; i < num_in; i++) {
      input_list.push_back(argv[i]);
    }
  }

  for (const auto &input : input_list) {
    calculateSingleTerm(input);
    std::cout << "\n";
  }

  return 0;
}
