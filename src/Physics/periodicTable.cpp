#include "Physics/periodicTable.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/AtomData_PeriodicTable.hpp"
#include "Physics/NuclearData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

namespace AtomData {

void instructions() {
  std::cout
      << "Usage example:\n"
         "Cs         Info for Cs with default A\n"
         "55         Same as above\n"
         "Cs 137     Info for Cs-137\n"
         "Cs all     Info for all available Cs isotopes\n"
         "(Note: numbers come from online database, and should be checked)\n";
}

void printData(const Nuclear::Isotope &nuc) {

  std::cout << AtomData::atomicSymbol(nuc.Z) << "-" << nuc.A << " (Z=" << nuc.Z
            << ", A=" << nuc.A << ")\n";

  std::cout << "r_rms = ";
  nuc.r_rms ? std::cout << *nuc.r_rms : std::cout << "??";
  std::cout << ", c = ";
  nuc.r_rms ? std::cout << Nuclear::c_hdr_formula_rrms_t(*nuc.r_rms) :
              std::cout << "??";
  std::cout << ", mu = ";
  nuc.mu ? std::cout << *nuc.mu : std::cout << "??";
  std::cout << ", Q = ";
  nuc.q ? std::cout << *nuc.q : std::cout << "??";

  std::cout << ", I = ";
  nuc.I_N ? std::cout << *nuc.I_N : std::cout << "??";
  std::cout << ", parity = ";
  nuc.parity ? std::cout << *nuc.parity : std::cout << "??";
  std::cout << "\n";
}

int parse_A(const std::string &A_str, int z) {

  int a = 0;
  if (qip::string_is_integer(A_str)) {
    a = std::stoi(A_str);
  } else {
    std::cout << "Invalid A: " << A_str << "\n";
    instructions();
    std::abort();
  }

  if (a <= 0)
    a = AtomData::defaultA(z);
  return a;
}
//==============================================================================
void printConstants() //
{
  std::cout << "\n";
  printf("1/α     = %.9f\n", PhysConst::c);
  printf("α       = %.14f\n", PhysConst::alpha);
  printf("c       = %.0f m/s\n", PhysConst::c_SI);
  printf("ℏ       = %.10e Js\n", PhysConst::hbar_SI);
  printf("ℏc      = %.7f MeV.fm\n", PhysConst::hbarc_MeVfm);
  printf("mp/me   = %.8f\n", PhysConst::m_p);
  printf("me      = %.11f MeV\n", PhysConst::m_e_MeV);
  printf("aB      = %.12e m\n", PhysConst::aB_m);
  printf("        = %.8f fm\n", PhysConst::aB_fm);
  printf("        = %.12e c/MeV\n", PhysConst::aB_fm / PhysConst::hbarc_MeVfm);
  printf("E_h     = %.12f eV\n", PhysConst::Hartree_eV);
  printf("        = %.12e Hz\n", PhysConst::Hartree_Hz);
  printf("        = %.8f /cm\n", PhysConst::Hartree_invcm);
  printf("λ(E_h)  = %.8f nm\n", 1e9 * PhysConst::c_SI / PhysConst::Hartree_Hz);
  printf("ℏ/E_h   = %.10e s\n", PhysConst::hbar_on_EH);
  std::cout << "\n";
  return;
}

//==============================================================================
void periodicTable(std::string z_str, std::string a_str) {
  if (a_str == "")
    a_str = "0";

  instructions();
  AtomData::printTable();

  const auto z = AtomData::atomic_Z(z_str);
  if (z != 0) {

    z_str = AtomData::atomicSymbol(z);

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
    if (isotopes.empty()) {
      std::cout << " none known\n";
    }
    for (const auto &nuc : isotopes) {
      std::cout << "\n";
      printData(nuc);
    }
  }
}

//==============================================================================

// Energy: au, cm^-1, Hz, nm, eV
// Length: a0, aB, fm
void convert_energies_au(double Eau) {

  const auto Ecm = Eau * PhysConst::Hartree_invcm;
  const auto EHz = Eau * PhysConst::Hartree_Hz;
  const auto EeV = Eau * PhysConst::Hartree_eV;
  const auto wl_nm = PhysConst::HartreeWL_nm / Eau;

  fmt::print("  {:<20.15g} Eh\n", Eau);
  fmt::print("  {:<20.15g} cm^-1\n", Ecm);
  if (EHz < 1.0e3)
    fmt::print("  {:<20.15g} Hz\n", EHz);
  else if (EHz < 1.0e6)
    fmt::print("  {:<20.15g} kHz\n", EHz * 1.0e-3);
  else if (EHz < 1.0e9)
    fmt::print("  {:<20.15g} MHz\n", EHz * 1.0e-6);
  else if (EHz < 1.0e12)
    fmt::print("  {:<20.15g} GHz\n", EHz * 1.0e-9);
  else
    fmt::print("  {:<20.15g} THz\n", EHz * 1.0e-12);

  if (EeV < 1.0e3)
    fmt::print("  {:<20.15g} eV\n", EeV);
  else if (EeV < 1.0e6)
    fmt::print("  {:<20.15g} keV\n", EeV * 1.0e-3);
  else if (EeV < 1.0e9)
    fmt::print("  {:<20.15g} MeV\n", EeV * 1.0e-6);
  else if (EeV < 1.0e12)
    fmt::print("  {:<20.15g} GeV\n", EeV * 1.0e-9);
  else
    fmt::print("  {:<20.15g} TeV\n", EeV * 1.0e-12);

  fmt::print("  {:<20.15g} nm\n", wl_nm);
}

//==============================================================================
void convert_length_au(double La0) {

  const auto Em = La0 * PhysConst::aB_m;

  fmt::print("  {:<20.15g} a0\n", La0);
  fmt::print("  {:<20.15g} fm\n", Em * 1.0e15);
  // always print fm, also others if better
  if (Em > 1.0e-12) {
    if (Em < 1.0e-9)
      fmt::print("  {:<20.15g} pm\n", Em * 1.0e12);
    else if (Em < 1.0e-6)
      fmt::print("  {:<20.15g} nm\n", Em * 1.0e9);
    else if (Em < 1.0e-3)
      fmt::print("  {:<20.15g} um\n", Em * 1.0e6);
    else if (Em > 1.0e-12)
      fmt::print("  {:<20.15g} mm\n", Em * 1.0e3);
  }
  std::cout << "\n";

  // in inverse eV^{-1};
  const auto L_inv_MeV = (Em * 1.0e15) / PhysConst::hbarc_MeVfm;

  // fmt::print("  {:<20.15g} MeV^-1\n", L_inv_MeV);

  if (L_inv_MeV < 1.0e-3)
    fmt::print("  {:<20.15g} TeV^-1\n", L_inv_MeV * 1.0e6);
  else if (L_inv_MeV < 1.0e-0)
    fmt::print("  {:<20.15g} GeV^-1\n", L_inv_MeV * 1.0e3);
  else if (L_inv_MeV < 1.0e3)
    fmt::print("  {:<20.15g} MeV^-1\n", L_inv_MeV);
  else if (L_inv_MeV < 1.0e6)
    fmt::print("  {:<20.15g} keV^-1\n", L_inv_MeV * 1e-3);
  else
    fmt::print("  {:<20.15g} eV^-1\n", L_inv_MeV * 1e-6);
}

//==============================================================================
void conversions(double number, const std::string &unit) {
  assert(unit.size() > 0);

  fmt::print("{:<.15g} {} = \n", number, unit);

  // this part: NOT case insensitive

  if (!qip::ci_contains(unit, "^-1")) {
    if (unit[0] == 'f')
      number *= 1.0e-15;
    if (unit[0] == 'p')
      number *= 1.0e-12;
    if (unit[0] == 'n')
      number *= 1.0e-9;
    if (unit[0] == 'u')
      number *= 1.0e-6;
    if (unit[0] == 'm' && unit.size() > 1)
      number *= 1.0e-3;
    if (unit[0] == 'k')
      number *= 1.0e3;
    if (unit[0] == 'M')
      number *= 1.0e6;
    if (unit[0] == 'G')
      number *= 1.0e9;
    if (unit[0] == 'T')
      number *= 1.0e12;
  } else {
    if (unit[0] == 'f')
      number /= 1.0e-15;
    if (unit[0] == 'p')
      number /= 1.0e-12;
    if (unit[0] == 'n')
      number /= 1.0e-9;
    if (unit[0] == 'u')
      number /= 1.0e-6;
    if (unit[0] == 'm' && unit.size() > 1)
      number /= 1.0e-3;
    if (unit[0] == 'k')
      number /= 1.0e3;
    if (unit[0] == 'M')
      number /= 1.0e6;
    if (unit[0] == 'G')
      number /= 1.0e9;
    if (unit[0] == 'T')
      number /= 1.0e12;
  }

  if (qip::ci_contains(unit, {"au", "a.u.", "Hartree", "Eh"})) {
    return convert_energies_au(number);
  }

  if (qip::ci_contains(unit, "Ry")) {
    return convert_energies_au(number / 2);
  }

  if (qip::ci_contains(unit, {"cm^-1", "/cm", "inversecm", "icm", "cm^{-1}"})) {
    return convert_energies_au(number / PhysConst::Hartree_invcm);
  }

  if (qip::ci_contains(unit, "Hz")) {
    return convert_energies_au(number / PhysConst::Hartree_Hz);
  }

  // Convert inverse energy to length:
  if (qip::ci_contains(unit, std::vector<std::string>{"eV^-1", "eV^{-1}"})) {
    number *= 1.0e6; //back to MeV:
    return convert_length_au(number * PhysConst::hbarc_MeVfm /
                             PhysConst::aB_fm);
  }

  //must come after the eV^-1
  if (qip::ci_contains(unit, "eV")) {
    return convert_energies_au(number / PhysConst::Hartree_eV);
  }

  // This must come _after_ inverse cm check
  // interpret as wavelength (convert to energy), then as a length:
  if ((unit.size() > 1 && unit[1] == 'm') || unit == "m") {
    std::cout << "Wavelength:\n";
    convert_energies_au(PhysConst::HartreeWL_nm / (1.0e9 * number));
    std::cout << "\n";
    convert_length_au((1.0e9 * number) / PhysConst::aB_nm);
  }

  if (qip::ci_contains(unit, std::vector<std::string>{"a0", "ab"})) {
    return convert_length_au(number);
  }
}

} // namespace AtomData