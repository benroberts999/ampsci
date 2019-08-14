#include "Module_runModules.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "DiracOperator.hpp"
#include "HartreeFockClass.hpp"
#include "Module_fitParametric.hpp"
#include "Module_matrixElements.hpp"
#include "Operators.hpp"
#include "UserInput.hpp"
#include "Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <string>
//
#include "Physics/PhysConst_constants.hpp"
#include <cmath>

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
  auto module_name = module_input.name();
  if (module_name.substr(0, 14) == "MatrixElements") {
    matrixElements(module_input, wf);
  } else if (module_name == "Module::Tests") {
    Module_tests(module_input, wf);
  } else if (module_name == "Module::WriteOrbitals") {
    Module_WriteOrbitals(module_input, wf);
  } else if (module_name == "Module::AtomicKernal") {
    atomicKernal(module_input, wf);
  } else if (module_name == "Module::FitParametric") {
    fitParametric(module_input, wf);
  } else if (module_name == "Module::BohrWeisskopf") {
    Module_BohrWeisskopf(module_input, wf);
  } else {
    std::cerr << "\nWARNING: Module `" << module_name
              << "' not known. Spelling mistake?\n";
  }
}

//******************************************************************************
// Below: some basic modules:

//******************************************************************************
void Module_BohrWeisskopf(const UserInputBlock &input, const Wavefunction &wf)
//
{
  UserInputBlock point_in("MatrixElements::hfs", input);
  UserInputBlock ball_in("MatrixElements::hfs", input);
  UserInputBlock BW_in("MatrixElements::hfs", input);
  point_in.add("F(r)=pointlike");
  ball_in.add("F(r)=ball");
  if (wf.Anuc() % 2 == 0)
    BW_in.add("F(r)=doublyOddBW");
  else
    BW_in.add("F(r)=VolotkaBW");

  auto hp = generateOperator("hfs", point_in, wf);
  auto hb = generateOperator("hfs", ball_in, wf);
  auto hw = generateOperator("hfs", BW_in, wf);

  std::cout << "\nTabulate A (Mhz), Bohr-Weisskopf effect: " << wf.atom()
            << "\n"
            << "state :     point       ball         BW |   F_ball      F_BW\n";
  for (const auto &phi : wf.valence_orbitals) {
    auto Ap = HyperfineOperator::hfsA2(hp.get(), phi);
    auto Ab = HyperfineOperator::hfsA2(hb.get(), phi);
    auto Aw = HyperfineOperator::hfsA2(hw.get(), phi);
    auto Fball = ((Ab / Ap) - 1.0) * M_PI * PhysConst::c;
    auto Fbw = ((Aw / Ap) - 1.0) * M_PI * PhysConst::c;
    printf("%6s: %9.3f  %9.3f  %9.3f | %8.4f  %8.4f\n", phi.symbol().c_str(),
           Ap, Ab, Aw, Fball, Fbw);
  }
}

//******************************************************************************
void Module_tests(const UserInputBlock &input, const Wavefunction &wf) {
  std::string ThisModule = "Module::Tests";
  if (input.get("orthonormal", true))
    Module_Tests_orthonormality(wf);
  if (input.get("Hamiltonian", false))
    Module_Tests_Hamiltonian(wf);
}

//------------------------------------------------------------------------------
void Module_Tests_orthonormality(const Wavefunction &wf) {
  std::cout << "\nTest orthonormality: ";
  std::cout << "log10(|1 - <a|a>|) or log10(|<a|b>|)\n";
  std::cout << "(should all read zero).\n";

  for (int i = 0; i < 3; i++) {
    const auto &tmp_b = (i == 2) ? wf.valence_orbitals : wf.core_orbitals;
    const auto &tmp_a = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;

    if (tmp_b.empty() || tmp_a.empty())
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
          xo -= 1.0;
        if (xo == 0)
          printf("   0");
        else
          printf(" %+3.0f", log10(std::fabs(xo)));
      }
      std::cout << "\n";
    }
  }
}

//------------------------------------------------------------------------------
void Module_Tests_Hamiltonian(const Wavefunction &wf) {
  std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";

  DirectHamiltonian Hd(wf.vnuc, wf.vdir, wf.get_alpha());
  for (const auto tmp_orbs : {&wf.core_orbitals, &wf.valence_orbitals}) {
    for (const auto &psi : *tmp_orbs) {
      double Haa = Hd.matrixEl(psi, psi) + psi * wf.get_VexPsi(psi);
      double ens = psi.en;
      double fracdiff = (Haa - ens) / ens;
      printf("<%i% i|H|%i% i> = %17.11f, E = %17.11f; % .0e\n", psi.n, psi.k,
             psi.n, psi.k, Haa, ens, fracdiff);
    }
  }
}

//******************************************************************************
void Module_WriteOrbitals(const UserInputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::WriteOrbitals";

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
  for (std::size_t i = 0; i < wf.rgrid.ngp; i++) {
    of << wf.rgrid.r[i] << " ";
    for (auto &psi : wf.core_orbitals)
      of << psi.f[i] << " ";
    for (auto &psi : wf.valence_orbitals)
      of << psi.f[i] << " ";
    of << "\n";
  }
  of << "\n# g block\n";
  for (std::size_t i = 0; i < wf.rgrid.ngp; i++) {
    of << wf.rgrid.r[i] << " ";
    for (auto &psi : wf.core_orbitals)
      of << psi.g[i] << " ";
    for (auto &psi : wf.valence_orbitals)
      of << psi.g[i] << " ";
    of << "\n";
  }
  of << "\n# density block\n";
  auto rho = wf.coreDensity();
  for (std::size_t i = 0; i < wf.rgrid.ngp; i++) {
    of << wf.rgrid.r[i] << " " << rho[i] << "\n";
  }
  of.close();
  std::cout << "Orbitals written to file: " << oname << "\n";
}

} // namespace Module
