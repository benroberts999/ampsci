#include "Physics/AtomData.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <iostream>
#include <vector>

void instructions() {
  std::cout
      << "Usage:\n"
         "$./periodicTable            Prints periodic table\n"
         "$./periodicTable Cs         Info for Cs with default A\n"
         "$./periodicTable Cs 137     Info for Cs-137\n"
         "$./periodicTable Cs all     Info for all available Cs isotopes\n"
         "(Note: numbers come from online database, and should be checked)\n";
  std::cout << "(Or, enter 'c' to print list of physics constants)\n";
}

void printData(const Nuclear::Isotope &nuc) {

  std::cout << AtomData::atomicSymbol(nuc.Z) << "-" << nuc.A << " (Z=" << nuc.Z
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
    a = AtomData::defaultA(z);
  return a;
}
//******************************************************************************
void printConstants() //
{
  std::cout << "\n";
  printf("1/alpha = %.9f\n", PhysConst::c);
  printf("alpha   = %.14f\n", PhysConst::alpha);
  printf("c       = %.0f m/s\n", PhysConst::c_SI);
  printf("hbar    = %.10e Js\n", PhysConst::hbar_SI);
  printf("mp/me   = %.8f\n", PhysConst::m_p);
  printf("me      = %.10f MeV\n", PhysConst::m_e_MeV);
  printf("aB      = %.12e m\n", PhysConst::aB_m);
  printf("Hy      = %.12f eV = %.12e Hz\n", PhysConst::Hartree_eV,
         PhysConst::Hartree_Hz);
  std::cout << "\n";
  return;
}

//******************************************************************************
int main(int num_in, char *argv[]) {

  if (num_in <= 1) {
    instructions();
    AtomData::printTable();
    return 1;
  }

  std::string z_str = argv[1];

  if (z_str.substr(0, 1) == "c") {
    printConstants();
    return 0;
  }

  const auto z = AtomData::get_z(z_str);
  if (z == 0) {
    instructions();
    return 1;
  }
  z_str = AtomData::atomicSymbol(z);
  std::string a_str = (num_in > 2) ? argv[2] : "0";
  auto name = AtomData::atomicName(z);

  std::vector<Nuclear::Isotope> isotopes;
  int a_default = parse_A("0", z);
  if (a_str == "all" || a_str == "list") {
    isotopes = Nuclear::findIsotopeList(z);
  } else {
    int a = parse_A(a_str, z);
    isotopes.push_back(Nuclear::findIsotopeData(z, a));
  }

  auto core_str = AtomData::guessCoreConfigStr(z);
  auto core_vec = AtomData::core_parser(core_str);

  std::cout << "\n"
            << z_str << ",  " << name << ".\n"
            << "Z = " << z << ";  A = " << a_default << " (default)\n\n";
  std::cout << "Electron config: " << core_str << "   (guess)\n"
            << " = ";
  for (const auto &term : core_vec) {
    if (term.frac() < 1)
      std::cout << "| ";
    std::cout << term.symbol() << " ";
  }
  std::cout << "\n";

  std::cout << "\nIsotpe data:";
  if (isotopes.size() == 0)
    std::cout << " none known\n";
  for (const auto &nuc : isotopes) {
    std::cout << "\n";
    printData(nuc);
  }
}
