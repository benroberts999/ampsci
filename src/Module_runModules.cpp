#include "Module_runModules.hpp"
#include "./DMionisation/Module_atomicKernal.hpp"
#include "DiracOperator.hpp"
#include "HartreeFockClass.hpp"
#include "Module_matrixElements.hpp"
#include "Operators.hpp"
#include "UserInput.hpp"
#include "Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <string>

namespace Module {

//******************************************************************************
void runModule(const UserInputBlock &module_input, const Wavefunction &wf,
               const HartreeFock &hf) //
{
  // ds
  auto module_name = module_input.name();
  if (module_name.substr(0, 14) == "MatrixElements") {
    matrixElements(module_input, wf);
  } else if (module_name == "Module::Tests") {
    Module_tests(module_input, wf, hf);
  } else if (module_name == "Module::WriteOrbitals") {
    Module_WriteOrbitals(module_input, wf, hf);
  } else if (module_name == "Module::AtomicKernal") {
    atomicKernal(module_input, wf);
  } else {
    std::cerr << "\nWARNING: Module `" << module_name
              << "' not known. Spelling mistake?\n";
  }
}

//******************************************************************************
// Below: some basic modules:

//******************************************************************************
void Module_tests(const UserInputBlock &input, const Wavefunction &wf,
                  const HartreeFock &hf) {
  std::string ThisModule = "Module::Tests";
  if (input.get("orthonormal", true))
    Module_Tests_orthonormality(wf);
  if (input.get("Hamiltonian", false))
    Module_Tests_Hamiltonian(wf, hf);
}

//------------------------------------------------------------------------------
void Module_Tests_orthonormality(const Wavefunction &wf) {
  std::cout << "Test orthonormality: ";
  std::cout << "log10(|1 - <a|a>|) or log10(|<a|b>|)\n";
  std::cout << "(should all read zero).\n";

  for (int i = 0; i < 3; i++) {
    const auto &tmp_b = (i == 2) ? wf.valence_orbitals : wf.core_orbitals;
    const auto &tmp_a = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;

    if (tmp_b.size() == 0 || tmp_a.size() == 0)
      continue;

    // Core-Core:
    if (i == 0)
      std::cout << "\nCore-Core\n    ";
    else if (i == 1)
      std::cout << "\nValence-Core\n    ";
    else
      std::cout << "\nValence-Valence\n    ";
    for (auto &psi_b : tmp_b)
      printf("%2i%2i", psi_b.n, psi_b.k);
    std::cout << "\n";
    for (auto &psi_a : tmp_a) {
      printf("%2i%2i", psi_a.n, psi_a.k);
      for (auto &psi_b : tmp_b) {
        if (psi_b > psi_a) {
          std::cout << "    ";
          continue;
        }
        if (psi_a.k != psi_b.k) {
          std::cout << "    ";
          continue;
        }
        double xo = (psi_a * psi_b);
        if (psi_a.n == psi_b.n)
          xo -= 1;
        if (xo == 0)
          printf("   0");
        else
          printf(" %+3.0f", log10(fabs(xo)));
      }
      std::cout << "\n";
    }
  }
}

//------------------------------------------------------------------------------
void Module_Tests_Hamiltonian(const Wavefunction &wf, const HartreeFock &hf) {
  std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";
  double c = 1. / wf.get_alpha();
  DiracOperator w(c, GammaMatrix::g5, 1, true);
  RadialOperator x_a(wf.rgrid, -1);
  DiracOperator y(c * c, DiracMatrix(0, 0, 0, -2));
  DiracOperator z1(wf.vnuc);
  DiracOperator z2(wf.vdir);
  for (const auto tmp_orbs : {&wf.core_orbitals, &wf.valence_orbitals}) {
    for (auto &psi : *tmp_orbs) {
      auto k = psi.k;
      // DiracOperator z3(hf.get_vex(psi));
      auto vexPsi = hf.vex_psia(psi);
      DiracOperator x_b(c, DiracMatrix(0, 1 - k, 1 + k, 0), 0, true);
      auto rhs = (w * psi) + (x_a * (x_b * psi)) + (y * psi) + (z1 * psi) +
                 (z2 * psi) + vexPsi;
      double R = psi * rhs;
      double ens = psi.en;
      double fracdiff = (R - ens) / ens;
      printf("<%i% i|H|%i% i> = %17.11f, E = "
             "%17.11f; % .0e\n",
             psi.n, psi.k, psi.n, psi.k, R, ens, fracdiff);
    }
  }
}

//******************************************************************************
void Module_WriteOrbitals(const UserInputBlock &input, const Wavefunction &wf,
                          const HartreeFock &) {
  const std::string ThisModule = "Module::WriteOrbitals";

  auto label = input.get<std::string>("label", "");
  std::string oname = wf.atomicSymbol() + "-orbitals";
  if (label != "")
    oname += "_" + label;
  oname += ".txt";

  std::ofstream of(oname);
  of << "r ";
  for (auto &psi : wf.core_orbitals)
    of << "\"" << psi.symbol(true) << "\" ";
  for (auto &psi : wf.valence_orbitals)
    of << "\"" << psi.symbol(true) << "\" ";
  of << "\n";
  for (std::size_t i = 0; i < wf.rgrid.ngp; i++) {
    of << wf.rgrid.r[i] << " ";
    for (auto &psi : wf.core_orbitals)
      of << psi.f[i] << " ";
    for (auto &psi : wf.valence_orbitals)
      of << psi.f[i] << " ";
    of << "\n";
  }
  std::cout << "Orbitals written to file: " << oname << "\n";
}

} // namespace Module
