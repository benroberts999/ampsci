#include "Modules/Module_runModules.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/Hamiltonian.hpp"
#include "Dirac/Operators.hpp"
#include "Dirac/Wavefunction.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/UserInput.hpp"
#include "Modules/Module_fitParametric.hpp"
#include "Modules/Module_matrixElements.hpp"
#include "Physics/PhysConst_constants.hpp"
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
    Module_WriteOrbitals(module_input, wf);
  } else if (module_name == "Module::AtomicKernal") {
    atomicKernal(module_input, wf);
  } else if (module_name == "Module::FitParametric") {
    fitParametric(module_input, wf);
  } else if (module_name == "Module::BohrWeisskopf") {
    Module_BohrWeisskopf(module_input, wf);
  } else if (module_name == "Module::pnc") {
    Module_testPNC(module_input, wf);
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
            << "\nstate :         point         ball           SP |   F_ball   "
               "  F_BW   eps(BW)\n";
  for (const auto &phi : wf.valence_orbitals) {
    auto Ap = HyperfineOperator::hfsA(hp.get(), phi);
    auto Ab = HyperfineOperator::hfsA(hb.get(), phi);
    auto Aw = HyperfineOperator::hfsA(hw.get(), phi);
    auto Fball = ((Ab / Ap) - 1.0) * M_PI * PhysConst::c;
    auto Fbw = ((Aw / Ap) - 1.0) * M_PI * PhysConst::c;
    // printf("%6s: %9.1f  %9.1f  %9.1f | %8.4f  %8.4f   %9.6f\n",
    printf("%7s: %12.5e %12.5e %12.5e | %8.4f %8.4f %9.6f\n",
           phi.symbol().c_str(), Ap, Ab, Aw, Fball, Fbw,
           -Fbw / (M_PI * PhysConst::c));
  }

  if (!wf.basis.empty())
    std::cout << "\nTest hfs using basis (pointlike):\n";
  for (const auto &phi : wf.basis) {
    auto Abasis = HyperfineOperator::hfsA(hp.get(), phi);
    // auto Abasis = hp.get()->radialIntegral(phi, phi);
    printf("%7s: %12.5e ", phi.symbol().c_str(), Abasis);
    const auto *hf_phi = wf.getState(phi.n, phi.k);
    if (hf_phi != nullptr) {
      auto Ahf = HyperfineOperator::hfsA(hp.get(), *hf_phi);
      auto delta = 2.0 * (Abasis - Ahf) / (Abasis + Ahf);
      printf(" %12.5e  %8.1e", Ahf, delta);
    }
    std::cout << "\n";
  }
}

//******************************************************************************
void Module_tests(const UserInputBlock &input, const Wavefunction &wf) {
  std::string ThisModule = "Module::Tests";
  input.checkBlock(
      {"orthonormal", "orthonormal_all", "Hamiltonian", "boundaries"});
  auto othon = input.get("orthonormal", true);
  auto othon_all = input.get("orthonormal_all", false);
  if (othon || othon_all)
    Module_Tests_orthonormality(wf, othon_all);
  if (input.get("Hamiltonian", false))
    Module_Tests_Hamiltonian(wf);
  if (input.get("boundaries", false))
    Module_test_r0pinf(wf);
}

//------------------------------------------------------------------------------
void Module_test_r0pinf(const Wavefunction &wf) {
  std::cout << "\nTesting boundaries r0 and pinf: f(r)/f_max\n";
  std::cout << " State    f(r0)   f(pinf)   pinf/Rinf\n";
  // for (const auto &phi : wf.core_orbitals)
  for (const auto tmp_orbs : {&wf.core_orbitals, &wf.valence_orbitals}) {
    for (const auto &phi : *tmp_orbs) {
      auto ratios = phi.r0pinfratio();
      printf("%7s:  %.0e   %.0e   %5i/%6.2f\n", phi.symbol().c_str(),
             std::abs(ratios.first), std::abs(ratios.second), (int)phi.pinf,
             wf.rgrid.r[phi.pinf - 1]);
      // std::cout << ratios.first << " " << ratios.second << "\n";
    }
    std::cout << "--------------\n";
  }
}

//------------------------------------------------------------------------------
void Module_Tests_orthonormality(const Wavefunction &wf, const bool print_all) {
  std::cout << "\nTest orthonormality: ";
  if (print_all) {
    std::cout << "log10(|1 - <a|a>|) or log10(|<a|b>|)\n";
    std::cout << "(should all read zero).";
  }
  std::cout << "\n";

  std::stringstream buffer;
  for (int i = 0; i < 6; i++) {
    // const auto &tmp_b = (i == 2) ? wf.valence_orbitals : wf.core_orbitals;
    // const auto &tmp_a = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;

    const auto &tmp_b = (i == 2 || i == 4)
                            ? wf.valence_orbitals
                            : (i != 5) ? wf.core_orbitals : wf.basis;
    // core, core, valence, core, valence, basis

    const auto &tmp_a =
        (i == 0) ? wf.core_orbitals : (i < 3) ? wf.valence_orbitals : wf.basis;
    // core, valence, valence, basis, basis

    if (tmp_b.empty() || tmp_a.empty())
      continue;

    // Core-Core:
    if (print_all) {
      if (i == 0)
        std::cout << "\nCore-Core\n    ";
      else if (i == 1)
        std::cout << "\nValence-Core\n    ";
      else if (i == 2)
        std::cout << "\nValence-Valence\n    ";
      else if (i == 3)
        std::cout << "\nBasis-core\n    ";
      else if (i == 4)
        std::cout << "\nBasis-Valence\n    ";
      else if (i == 5)
        std::cout << "\nBasis-Basis\n    ";
    } else {
      if (i == 0)
        buffer << "cc ";
      else if (i == 1)
        buffer << "vc ";
      else if (i == 2)
        buffer << "vv ";
      else if (i == 3)
        buffer << "bc ";
      else if (i == 4)
        buffer << "bv ";
      else if (i == 5)
        buffer << "bb ";
    }

    auto worst_xo = 0.0;
    std::string worst_braket = "";
    if (print_all) {
      for (auto &psi_b : tmp_b)
        printf("%2i%2i", psi_b.n, psi_b.k);
      std::cout << "\n";
    }
    for (auto &psi_a : tmp_a) {
      if (print_all)
        printf("%2i%2i", psi_a.n, psi_a.k);
      for (auto &psi_b : tmp_b) {
        if (psi_b > psi_a && (&tmp_b == &tmp_a)) {
          if (print_all)
            std::cout << "    ";
          continue;
        }
        if (psi_a.k != psi_b.k) {
          if (print_all)
            std::cout << "    ";
          continue;
        }
        double xo = (psi_a * psi_b);
        if (psi_a.n == psi_b.n) {
          xo -= 1.0;
        } else {
          if (std::abs(xo) > std::abs(worst_xo)) {
            worst_xo = xo;
            worst_braket = "<" + psi_a.symbol() + "|" + psi_b.symbol() + ">";
          }
        }
        if (print_all) {
          if (xo == 0)
            printf("   0");
          else
            printf(" %+3.0f", std::log10(std::fabs(xo)));
        }
      } // psi_b
      if (print_all)
        std::cout << "\n";
    } // Psi_a
    if (worst_braket != "") {
      std::string eq = worst_xo > 0 ? " =  " : " = ";
      buffer << worst_braket << eq << std::setprecision(2) << std::scientific
             << worst_xo;
    }
    buffer << "\n";
  } // cc, cv, vv
  if (print_all)
    std::cout << "\n";
  std::cout << buffer.str();
}

//------------------------------------------------------------------------------
void Module_Tests_Hamiltonian(const Wavefunction &wf) {
  std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";

  // DirectHamiltonian Hd(wf.vnuc, wf.vdir, wf.get_alpha());
  auto Hd = RadialHamiltonian(wf.rgrid, wf.get_alpha());
  Hd.set_v(-1, wf.vnuc, wf.vdir); // same each kappa

  for (const auto tmp_orbs :
       {&wf.core_orbitals, &wf.valence_orbitals, &wf.basis}) {
    if (tmp_orbs->empty())
      continue;
    double worst_eps = 0.0;
    const DiracSpinor *worst_psi = nullptr;
    for (const auto &psi : *tmp_orbs) {
      double Haa_d = Hd.matrixEl(psi, psi);
      double Haa_x =
          (tmp_orbs != &wf.basis)
              ? psi * wf.get_VexPsi(psi)
              : psi * HartreeFock::vex_psia_any(psi, wf.core_orbitals);
      auto Haa = Haa_d + Haa_x;
      double ens = psi.en;
      double fracdiff = (Haa - ens) / ens;
      printf("<%2i% i|H|%2i% i> = %17.11f, E = %17.11f; % .0e\n", psi.n, psi.k,
             psi.n, psi.k, Haa, ens, fracdiff);
      if (std::abs(fracdiff) >= std::abs(worst_eps)) {
        worst_eps = fracdiff;
        worst_psi = &psi;
      }
    }
    std::cout << worst_psi->symbol() << ": eps=" << worst_eps << "\n";
  }
}

//******************************************************************************
void Module_WriteOrbitals(const UserInputBlock &input, const Wavefunction &wf) {
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

} // namespace Module
