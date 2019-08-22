#include "Module_runModules.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "Dirac/DiracOperator.hpp"
#include "HF/HartreeFockClass.hpp"
#include "Module_fitParametric.hpp"
#include "Module_matrixElements.hpp"
#include "Dirac/Operators.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "IO/UserInput.hpp"
#include "Dirac/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
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

  std::cout
      << "\nTabulate A (Mhz), Bohr-Weisskopf effect: " << wf.atom() << "\n"
      << "state :         point          ball            BW |   F_ball    "
         "  F_BW   "
         "  eps(BW)\n";
  for (const auto &phi : wf.valence_orbitals) {
    auto Ap = HyperfineOperator::hfsA(hp.get(), phi);
    auto Ab = HyperfineOperator::hfsA(hb.get(), phi);
    auto Aw = HyperfineOperator::hfsA(hw.get(), phi);
    auto Fball = ((Ab / Ap) - 1.0) * M_PI * PhysConst::c;
    auto Fbw = ((Aw / Ap) - 1.0) * M_PI * PhysConst::c;
    // printf("%6s: %9.1f  %9.1f  %9.1f | %8.4f  %8.4f   %9.6f\n",
    printf("%7s: %12.5e  %12.5e  %12.5e | %8.4f  %8.4f   %9.6f\n",
           phi.symbol().c_str(), Ap, Ab, Aw, Fball, Fbw,
           -Fbw / (M_PI * PhysConst::c));
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
          printf(" %+3.0f", std::log10(std::fabs(xo)));
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
void Module_testPNC(const UserInputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::PNC";

  // ChronoTimer("pnc");
  auto t_dflt = Nuclear::default_t;
  auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  auto c_dflt = Nuclear::c_hdr_formula_rrms_t(r_rms);
  auto t = input.get("t", t_dflt);
  auto c = input.get("c", c_dflt);

  auto transition_str = input.get<std::string>("transition");
  replace(transition_str.begin(), transition_str.end(), ',', ' ');
  auto ss = std::stringstream(transition_str);
  int na, ka, nb, kb;
  ss >> na >> ka >> nb >> kb;

  auto ncore = wf.maxCore_n();
  auto main_n = input.get("nmain", ncore + 4);

  PNCnsiOperator hpnc(c, t, wf.rgrid, -wf.Nnuc());
  E1Operator he1(wf.rgrid);

  const auto &a6s = wf.getState(na, ka);
  const auto &a7s = wf.getState(nb, kb);
  std::cout << "\nE_pnc: " << wf.atom() << ": " << a6s.symbol() << " -> "
            << a7s.symbol() << "\n";

  auto tja = a6s.twoj();
  auto tjb = a7s.twoj();
  auto twom = std::min(tja, tjb);
  auto c10 = Wigner::threej_2(tjb, 2, tja, -twom, 0, twom) *
             Wigner::threej_2(tja, 0, tja, -twom, 0, twom);
  auto c01 = Wigner::threej_2(tjb, 0, tjb, -twom, 0, twom) *
             Wigner::threej_2(tjb, 2, tja, -twom, 0, twom);

  double pnc = 0, core = 0, main = 0;
  for (int i = 0; i < 2; i++) {
    auto &tmp_orbs = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
    for (auto &np : tmp_orbs) {
      // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
      if (np == a7s || np == a6s)
        continue;
      if (hpnc.isZero(np, a6s) && hpnc.isZero(np, a7s))
        continue;
      double pnc1 = c10 * he1.reducedME(a7s, np) * hpnc.reducedME(np, a6s) /
                    (a6s.en - np.en);
      double pnc2 = c01 * hpnc.reducedME(a7s, np) * he1.reducedME(np, a6s) /
                    (a7s.en - np.en);
      std::cout << np.symbol() << ", pnc= ";
      printf("%12.5e + %12.5e = %12.5e\n", pnc1, pnc2, pnc1 + pnc2);
      pnc += pnc1 + pnc2;
      if (np.n == main_n)
        main = pnc - core;
    }
    if (i == 0)
      core = pnc;
  }
  std::cout << "Core= " << core << "\n";
  std::cout << "Main= " << main << "\n";
  std::cout << "Tail= " << pnc - main - core << "\n";
  std::cout << "Total= " << pnc << "\n";
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
