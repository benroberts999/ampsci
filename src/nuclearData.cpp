#include "Physics/AtomInfo.hpp"
#include "Physics/Nuclear.hpp"
#include <iostream>

void instructions() {
  std::cout << "Input is Z A; Z may be int or string. e.g.: \n "
               "$./nuclearData Cs 133\n"
               "Leave A blank (or put 0) to get dafault A value.\n"
               "Put 'all' to list all available A values.\n"
               "Note: numbers come from online database, and have some errors, "
               "so should be checked if needed.\n";
}

void printData(const Nuclear::Isotope &nuc) {

  std::cout << AtomInfo::atomicSymbol(nuc.Z) << "-" << nuc.A << " (Z=" << nuc.Z
            << ", A=" << nuc.A << ")\n";

  std::cout << "r_rms = ";
  nuc.r_ok() ? std::cout << nuc.r_rms : std::cout << "??";
  std::cout << ", c = ";
  nuc.r_ok() ? std::cout << Nuclear::c_hdr_formula_rrms_t(nuc.r_rms)
             : std::cout << "??";
  std::cout << ", mu = ";
  nuc.mu_ok() ? std::cout << nuc.mu : std::cout << "??";

  std::cout << ", I = ";
  nuc.I_ok() ? std::cout << nuc.I_N : std::cout << "??";
  std::cout << ", parity = ";
  nuc.parity_ok() ? std::cout << nuc.parity : std::cout << "??";
  std::cout << "\n";
}

int parse_A(const std::string &A_str, int z = 0) {
  int a = 0;
  try {
    a = std::stoi(A_str);
  } catch (...) {
    std::cout << "Invalid A: " << A_str << "\n";
    instructions();
    std::abort();
  }

  if (a <= 0)
    a = AtomInfo::defaultA(z);
  return a;
}

//******************************************************************************
int main(int num_in, char *argv[]) {

  if (num_in <= 1) {
    instructions();
    return 1;
  }

  std::string z_str = argv[1];
  const auto z = AtomInfo::get_z(z_str);
  std::string a_str = (num_in > 2) ? argv[2] : "0";

  std::vector<Nuclear::Isotope> isotopes;
  if (a_str == "all" || a_str == "list") {
    isotopes = Nuclear::findIsotopeList(z);
  } else {
    int a = parse_A(a_str, z);
    isotopes.push_back(Nuclear::findIsotopeData(z, a));
  }

  std::cout << AtomInfo::guessCoreConfigStr(z) << "\n";
  for (const auto &nuc : isotopes) {
    std::cout << "\n";
    printData(nuc);
  }
}
