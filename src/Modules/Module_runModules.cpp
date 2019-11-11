#include "Modules/Module_runModules.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/Operators.hpp"
#include "Dirac/Wavefunction.hpp"
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
            << "\n"
            << "state :         point         ball           BW |   F_ball   "
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
             wf.rgrid.r[phi.pinf]);
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
  for (int i = 0; i < 3; i++) {
    const auto &tmp_b = (i == 2) ? wf.valence_orbitals : wf.core_orbitals;
    const auto &tmp_a = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;

    if (tmp_b.empty() || tmp_a.empty())
      continue;

    // Core-Core:
    if (print_all) {
      if (i == 0)
        std::cout << "\nCore-Core\n    ";
      else if (i == 1)
        std::cout << "\nValence-Core\n    ";
      else
        std::cout << "\nValence-Valence\n    ";
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
        if (psi_b > psi_a) {
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
             << worst_xo << "\n";
    }
  } // cc, cv, vv
  if (print_all)
    std::cout << "\n";
  std::cout << buffer.str();
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

  input.checkBlock({"t", "c", "transition", "nmain"});

  // std::cout << "\nPNC test:\n";
  // ChronoTimer("pnc");
  auto t_dflt = Nuclear::default_t;
  auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  auto c_dflt = Nuclear::c_hdr_formula_rrms_t(r_rms);
  auto t = input.get("t", t_dflt);
  auto c = input.get("c", c_dflt);

  auto transition_str = input.get<std::string>("transition");
  std::replace(transition_str.begin(), transition_str.end(), ',', ' ');
  auto ss = std::stringstream(transition_str);
  int na, ka, nb, kb;
  ss >> na >> ka >> nb >> kb;

  auto ncore = wf.maxCore_n();
  auto main_n = input.get("nmain", ncore + 4);

  PNCnsiOperator hpnc(c, t, wf.rgrid, -wf.Nnuc());
  E1Operator he1(wf.rgrid);
  auto alpha = wf.get_alpha();

  const auto &aA = wf.getState(na, ka);
  const auto &aB = wf.getState(nb, kb);
  std::cout << "\n********************************************** \n";
  std::cout << "E_pnc: " << wf.atom() << ":   A = " << aA.symbol()
            << "   -->   B = " << aB.symbol() << "\n\n";

  auto tja = aA.twoj();
  auto tjb = aB.twoj();
  auto twom = std::min(tja, tjb);
  auto c10 = he1.rme3js(tja, tjb, twom) * hpnc.rme3js(tjb, tjb, twom);
  auto c01 = hpnc.rme3js(tja, tja, twom) * he1.rme3js(tja, tjb, twom);

  {
    std::cout << "Sum-over-states method (HF valence):\n";
    std::cout << " <A|d|n><n|hw|B>/dEB + <A|hw|n><n|d|B>/dEA\n";
    double pnc = 0, core = 0, main = 0;
    for (int i = 0; i < 2; i++) {
      auto &tmp_orbs = (i == 0) ? wf.core_orbitals : wf.valence_orbitals;
      for (auto &np : tmp_orbs) {
        // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
        if (np == aB || np == aA)
          continue;
        if (hpnc.isZero(np.k, aA.k) && hpnc.isZero(np.k, aB.k))
          continue;
        double pnc1 = c10 * he1.reducedME(aA, np) * hpnc.reducedME(np, aB) /
                      (aB.en - np.en);
        double pnc2 = c01 * hpnc.reducedME(aA, np) * he1.reducedME(np, aB) /
                      (aA.en - np.en);
        printf("%7s, pnc= %12.5e + %12.5e = %12.5e\n", np.symbol().c_str(),
               pnc1, pnc2, pnc1 + pnc2);
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

  {
    auto v = NumCalc::add_vectors(wf.vnuc, wf.vdir);
    auto v1 = 0.0, v2 = 0.0;

    auto hA_dag = hpnc.reduced_lhs(-aA.k, aA);
    auto hB = hpnc.reduced_rhs(-aB.k, aB);

    auto e1A = he1.reduced_lhs(-aB.k, aA);
    auto e1B = he1.reduced_rhs(-aA.k, aB);

    auto del_A_dag = HartreeFock::solveMixedState(aA, hA_dag.k, 0, v, alpha,
                                                  wf.core_orbitals, hA_dag);
    auto del_B = HartreeFock::solveMixedState(aB, hB.k, 0, v, alpha,
                                              wf.core_orbitals, hB);

    auto omega = aB.en - aA.en;
    auto xA = HartreeFock::solveMixedState(aA, e1A.k, -omega, v, alpha,
                                           wf.core_orbitals, e1A);
    auto yB = HartreeFock::solveMixedState(aB, e1B.k, omega, v, alpha,
                                           wf.core_orbitals, e1B);

    auto pnc1_w = c01 * he1.reducedME(del_A_dag, aB);
    auto pnc2_w = c10 * he1.reducedME(aA, del_B);

    std::cout << "\nMixed states method: \n";
    std::cout << "<dA |d| B> + <A |d| dB> = ";
    std::cout << pnc1_w << " + " << pnc2_w << " = " << pnc1_w + pnc2_w << "\n";

    auto pnc1_d = c01 * hpnc.reducedME(aA, yB);
    auto pnc2_d = c10 * hpnc.reducedME(xA, aB);
    std::cout << "<A |h|Y_B> + <X_A|h| B> = ";
    std::cout << pnc1_d << " + " << pnc2_d << " = " << pnc1_d + pnc2_d << "\n";

    v1 = pnc1_w + pnc2_w;
    v2 = pnc1_d + pnc2_d;
    std::cout << "eps = " << (v1 - v2) / (0.5 * (v1 + v2)) << "\n";

    // orthog wrt core:
    wf.orthogonaliseWrtCore(del_A_dag);
    wf.orthogonaliseWrtCore(del_B);
    // Core contribution:
    auto pnc1_c = pnc1_w - c01 * he1.reducedME(del_A_dag, aB);
    auto pnc2_c = pnc2_w - c10 * he1.reducedME(aA, del_B);
    for (const auto &phiv : wf.valence_orbitals) {
      if (phiv.n > main_n)
        continue;
      if (phiv.k == del_A_dag.k) {
        del_A_dag -= (del_A_dag * phiv) * phiv;
      }
      if (phiv.k == del_B.k) {
        del_B -= (del_B * phiv) * phiv;
      }
    }

    // Main contribution:
    auto pnc1_m = pnc1_w - pnc1_c - c10 * he1.reducedME(del_A_dag, aB);
    auto pnc2_m = pnc2_w - pnc2_c - c01 * he1.reducedME(aA, del_B);

    std::cout << "\n<dA'|d| B>  +  <A |d|dB'>  (by force orthog):\n";
    std::cout << "core : " << pnc1_c + pnc2_c << "\n";
    std::cout << "main : " << pnc1_m + pnc2_m << "\n";
    std::cout << "tail : "
              << pnc1_w - pnc1_m - pnc1_c + pnc2_w - pnc2_m - pnc2_c << "\n";
    std::cout << "Total: " << v1 << "\n";
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
