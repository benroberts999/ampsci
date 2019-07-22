#pragma once
#include "AtomInfo.hpp"
#include "ChronoTimer.hpp"
#include "DiracOperator.hpp"
#include "HartreeFockClass.hpp"
#include "Module_atomicKernal.hpp"
#include "Nuclear.hpp"
#include "Operators.hpp"
#include "PhysConst_constants.hpp"
#include "UserInput.hpp"
#include "Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <string>

inline void Module_orthonormality(const Wavefunction &wf);
inline void Module_Hamiltonian(const Wavefunction &wf, const HartreeFock &hf);
inline void Module_WriteOrbitals(const Wavefunction &wf,
                                 const std::string &label);

// Use enums to avoid looping through by specialising template?
//******************************************************************************
inline void RunModule(const std::string &module, const UserInput &input,
                      const Wavefunction &wf, const HartreeFock &hf) //
{
  //
  if (module == "Module::Tests") {
    if (input.get(module, "orthonormality", true))
      Module_orthonormality(wf);
    if (input.get(module, "Hamiltonian", false))
      Module_Hamiltonian(wf, hf);
  } else if (module == "Module::WriteOrbitals") {
    Module_WriteOrbitals(wf, input.get<std::string>(module, "label", ""));
  } else if (module == "Module::AtomicKernal") {
    Module::atomicKernal(input, wf);
  }
}

//******************************************************************************
inline void Module_orthonormality(const Wavefunction &wf) {
  std::cout << "Test orthonormality [log-scale, should all read 0]:\n";
  for (int i = 0; i < 3; i++) {
    const auto &tmp_b = (i == 2) ? wf.valence_orbitals : wf.core_orbitals;
    const auto &tmp_a = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
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
  std::cout << "\n(Note: Core orbitals are orthogonalised explicitely, as is "
               "each valence state [with respect to the core]. However, "
               "valence states are not explicitely orthogonalised wrt each "
               "other, since there's no self-consistent way to do this with "
               "a finite set of valence orbitals).\n";
}

//******************************************************************************
inline void Module_Hamiltonian(const Wavefunction &wf, const HartreeFock &hf) {
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
      DiracOperator z3(hf.get_vex(psi));
      auto vexPsi = hf.vex_psia(psi);
      DiracOperator x_b(c, DiracMatrix(0, 1 - k, 1 + k, 0), 0, true);
      auto rhs = (w * psi) + (x_a * (x_b * psi)) + (y * psi) + (z1 * psi) +
                 (z2 * psi) + vexPsi;
      double R = psi * rhs;
      double ens = psi.en;
      double fracdiff = (R - ens) / ens;
      printf("<%i% i|H|%i% i> = %17.11f, E = %17.11f; % .0e\n", psi.n, psi.k,
             psi.n, psi.k, R, ens, fracdiff);
    }
  }
}

//******************************************************************************
inline void Module_WriteOrbitals(const Wavefunction &wf,
                                 const std::string &label) {
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
