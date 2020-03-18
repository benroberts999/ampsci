#include "Modules/Module_runModules.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/UserInput.hpp"
#include "Modules/Module_fitParametric.hpp"
#include "Modules/Module_matrixElements.hpp"
#include "Modules/Module_pnc.hpp"
#include "Modules/Module_tests.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Hamiltonian.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace Module {

//******************************************************************************
void runModules(const UserInput &input, const Wavefunction &wf) {
  auto modules = input.module_list();
  for (const auto &module : modules) {
    runModule(module, wf);
  }
}

//******************************************************************************
void runModule(const UserInputBlock &module_input, const Wavefunction &wf) //
{
  const auto &module_name = module_input.name();
  if (module_name.substr(0, 14) == "MatrixElements") {
    matrixElements(module_input, wf);
  } else if (module_name == "Module::Tests") {
    Module_tests(module_input, wf);
  } else if (module_name == "Module::WriteOrbitals") {
    writeOrbitals(module_input, wf);
  } else if (module_name == "Module::AtomicKernal") {
    atomicKernal(module_input, wf);
  } else if (module_name == "Module::FitParametric") {
    fitParametric(module_input, wf);
  } else if (module_name == "Module::BohrWeisskopf") {
    calculateBohrWeisskopf(module_input, wf);
  } else if (module_name == "Module::pnc") {
    calculatePNC(module_input, wf);
  } else if (module_name == "Module::lifetimes") {
    calculateLifetimes(module_input, wf);
  } else if (module_name == "Module::SecondOrder") {
    SecondOrder(module_input, wf);
  } else {
    std::cerr << "\nWARNING: Module " << module_name
              << " not known. Spelling mistake?\n";
  }
}

//******************************************************************************
void writeOrbitals(const UserInputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::WriteOrbitals";
  input.checkBlock({"label"});

  std::cout << "\n Running: " << ThisModule << "\n";
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
  of << "# f block\n";
  for (std::size_t i = 0; i < wf.rgrid.num_points; i++) {
    of << wf.rgrid.r[i] << " ";
    for (auto &psi : wf.core_orbitals)
      of << psi.f[i] << " ";
    for (auto &psi : wf.valence_orbitals)
      of << psi.f[i] << " ";
    of << "\n";
  }
  of << "\n# g block\n";
  for (std::size_t i = 0; i < wf.rgrid.num_points; i++) {
    of << wf.rgrid.r[i] << " ";
    for (auto &psi : wf.core_orbitals)
      of << psi.g[i] << " ";
    for (auto &psi : wf.valence_orbitals)
      of << psi.g[i] << " ";
    of << "\n";
  }
  of << "\n# density block\n";
  auto rho = wf.coreDensity();
  for (std::size_t i = 0; i < wf.rgrid.num_points; i++) {
    of << wf.rgrid.r[i] << " " << rho[i] << "\n";
  }
  of.close();
  std::cout << "Orbitals written to file: " << oname << "\n";
}

//******************************************************************************
void SecondOrder(const UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"lmax", "kmax"});
  auto kmax = input.get("kmax", 10);
  auto lmax = input.get("lmax", 10);

  for (const auto &v : wf.valence_orbitals) {
    if (v.l() > lmax)
      continue;

    double delta = 0.0;
    for (int k = 0; k <= kmax; ++k) {
      auto f = (2 * k + 1) * v.twojp1();
      double sigma_k = 0.0;
#pragma omp parallel for
      for (auto ib = 0ul; ib < wf.core_orbitals.size(); ib++) {
        double sigma_b = 0.0;
        const auto &b = wf.core_orbitals[ib];
        for (const auto &n : wf.basis) {
          if (wf.isInCore(n))
            continue;

          for (const auto &a : wf.core_orbitals) {
            auto zx = Coulomb::Zk_abcd(v, n, a, b, k) *
                      Coulomb::Xk_abcd(v, n, a, b, k);
            auto dele = v.en + n.en - a.en - b.en;
            sigma_b += zx / dele;
          }
          for (const auto &m : wf.basis) {
            if (wf.isInCore(m))
              continue;
            auto zx = Coulomb::Zk_abcd(m, n, v, b, k) *
                      Coulomb::Xk_abcd(m, n, v, b, k);
            auto dele = m.en + n.en - v.en - b.en;
            sigma_b -= zx / dele;
          }
        }

#pragma omp critical(sum)
        { sigma_k += sigma_b; }
      }
      delta += sigma_k / f;
    } // k
    std::cout << v.en << " + " << delta << " = " << (v.en + delta) << " = "
              << (v.en + delta) * PhysConst::Hartree_invcm << "\n";
  }
}

} // namespace Module
