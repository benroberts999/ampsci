#include "Modules/Module_tests.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/UserInput.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Hamiltonian.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace Module {

//******************************************************************************
void Module_tests(const UserInputBlock &input, const Wavefunction &wf) {
  using namespace Tests;
  std::string ThisModule = "Module::Tests";
  input.checkBlock({"orthonormal", "orthonormal_all", "Hamiltonian",
                    "boundaries", "sumRules"});
  auto othon = input.get("orthonormal", true);
  auto othon_all = input.get("orthonormal_all", false);
  if (othon || othon_all)
    Module_Tests_orthonormality(wf, othon_all);
  if (input.get("Hamiltonian", false))
    Module_Tests_Hamiltonian(wf);
  if (input.get("boundaries", false))
    Module_test_r0pinf(wf);
  if (input.get("sumRules", false))
    basisSumRules(wf);
}

namespace Tests {

//------------------------------------------------------------------------------
void basisSumRules(const Wavefunction &wf) {

  if (wf.basis.empty())
    return;

  std::cout << "\nTesting basis (sum rules):\n";
  std::cout << "(must include +ve energy states. Works best for pure Coloumb "
               "functions)\n";

  auto rhat = DiracOperator::E1(wf.rgrid);          // vector E1
  auto r2hat = DiracOperator::RadialF(wf.rgrid, 2); // scalar r^2

  auto comp_l = [](const auto &Fa, const auto &Fb) { return Fa.l() < Fb.l(); };
  auto max_l = std::max_element(wf.basis.begin(), wf.basis.end(), comp_l)->l();

  std::cout << "TKR sum rule\n";
  {
    auto Fa = wf.basis.front();
    for (int l = 0; l <= max_l; l++) {
      auto sum_el = 0.0;
      auto sum_p = 0.0;
      for (const auto &Fn : wf.basis) {
        if (Fn == Fa)
          continue;
        auto f = (Fn.k == l) ? l : (Fn.k == -l - 1) ? l + 1 : 0;
        if (f == 0)
          continue;
        auto w = Fn.en - Fa.en;
        // auto Ran = rhat.radialIntegral(Fa, Fn);
        auto Ran = Fa * (wf.rgrid.r * Fn);
        auto term = f * w * Ran * Ran / (2 * l + 1);
        if (Fn.n > 0)
          sum_el += term;
        else
          sum_p += term;
      }
      printf("l=%1i, sum = %10.6f%+10.6f = %8.1e\n", l, sum_el, sum_p,
             sum_el + sum_p);
    }
  }

  auto comp_ki = [](const auto &Fm, const auto &Fn) {
    return Fm.k_index() < Fn.k_index();
  };
  auto max_ki =
      std::max_element(wf.basis.begin(), wf.basis.end(), comp_ki)->k_index();

  std::cout << "Drake-Goldman sum rules: w^n |<a|r|b>|^2  (n=0,1,2)\n";
  for (int ki = 0; ki <= max_ki; ki++) {
    auto kappa = Angular::kappaFromIndex(ki);
    auto comp_k = [=](const auto &Fn) { return Fn.k == kappa; };
    auto Fa = *std::find_if(wf.basis.begin(), wf.basis.end(), comp_k);
    // need to have l_n = la+1 terms, or sum doesn't work:
    if (Fa.l() == max_l)
      continue;
    std::cout << "kappa: " << kappa << " (" << Fa.symbol() << ")\n";
    for (int i = 0; i < 3; i++) {
      auto sum = 0.0;
      for (const auto &Fn : wf.basis) {
        auto w = Fn.en - Fa.en;
        auto Ran = rhat.reducedME(Fa, Fn);
        double c = 1.0 / (2 * std::abs(Fa.k));
        auto term = std::pow(w, i) * Ran * Ran * c;
        sum += term;
      }
      if (i == 2)
        sum *= wf.get_alpha() * wf.get_alpha() / 3;
      auto s0 = (i == 0) ? r2hat.radialIntegral(Fa, Fa) : (i == 1) ? 0.0 : 1.0;
      printf("%i: sum=%9.6f, exact=%+9.6f, diff=%8.1e\n", i, sum, s0, sum - s0);
    }
  }
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
  Hd.set_v_mag(wf.Hse_mag);

  for (const auto tmp_orbs :
       {&wf.core_orbitals, &wf.valence_orbitals, &wf.basis}) {
    if (tmp_orbs->empty())
      continue;
    double worst_eps = 0.0;
    const DiracSpinor *worst_psi = nullptr;
    for (const auto &psi : *tmp_orbs) {
      double Haa_d = Hd.matrixEl(psi, psi);
      // double Haa_x =
      //     (tmp_orbs != &wf.basis)
      //         ? psi * wf.get_VexPsi(psi)
      //         : psi * vex_psia_any(psi, wf.core_orbitals);
      double Haa_x = psi * HF::vex_psia_any(psi, wf.core_orbitals);
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
    if (worst_psi != nullptr)
      std::cout << worst_psi->symbol() << ": eps=" << worst_eps << "\n";
  }
}

} // namespace Tests

} // namespace Module
