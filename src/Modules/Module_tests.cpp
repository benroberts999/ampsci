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
void Module_tests(const IO::UserInputBlock &input, const Wavefunction &wf) {
  using namespace Tests;
  std::string ThisModule = "Module::Tests";
  input.checkBlock({"orthonormal", "orthonormal_all", "Hamiltonian",
                    "boundaries", "basisTests"});
  auto othon = input.get("orthonormal", true);
  auto othon_all = input.get("orthonormal_all", false);
  if (othon || othon_all)
    Module_Tests_orthonormality(wf, othon_all);
  if (input.get("Hamiltonian", false))
    Module_Tests_Hamiltonian(wf);
  if (input.get("boundaries", false))
    Module_test_r0pinf(wf);
  if (input.get("basisTests", false))
    basisTests(wf);
}

namespace Tests {

namespace Helper {
int countNodes(const DiracSpinor &Fn)
// Just counts the number of times orbital (f) changes sign
{
  double sp = Fn.f[Fn.p0 + 3];
  int counted_nodes = 0;
  for (auto i = Fn.p0 + 4; i < Fn.pinf - 3; ++i) {
    if (sp * Fn.f[i] < 0) {
      ++counted_nodes;
      sp = Fn.f[i];
    }
  }
  return counted_nodes;
}
} // namespace Helper

//------------------------------------------------------------------------------
void basisTests(const Wavefunction &wf) {

  std::cout << "\nTesting basis/spectrum (sum rules):\n";
  std::cout << "(Must include +ve energy states.)\n";

  const auto &basis = wf.spectrum.empty() ? wf.basis : wf.spectrum;
  if (basis.empty())
    return;

  if (&basis == &(wf.spectrum))
    std::cout << "Using Sprectrum\n";
  else
    std::cout << "Using Basis\n";

  auto rhat = DiracOperator::E1(wf.rgrid);          // vector E1
  auto r2hat = DiracOperator::RadialF(wf.rgrid, 2); // scalar r^2

  auto comp_l = [](const auto &Fa, const auto &Fb) { return Fa.l() < Fb.l(); };
  auto max_l = std::max_element(basis.begin(), basis.end(), comp_l)->l();

  // check to see if there are any negative-energy states:
  const bool negative_statesQ = std::any_of(
      basis.cbegin(), basis.cend(), [](const auto &Fa) { return Fa.n < 0; });

  //----------------------------------------------------------------------------
  if (negative_statesQ) {
    std::cout << "\nTKR sum rule (should =0)\n";
    {
      const auto &Fa = basis.front();
      for (int l = 0; l <= max_l; l++) {
        auto sum_el = 0.0;
        auto sum_p = 0.0;
        for (const auto &Fn : basis) {
          if (Fn == Fa)
            continue;
          const auto f = (Fn.k == l) ? l : (Fn.k == -l - 1) ? l + 1 : 0;
          if (f == 0)
            continue;
          const auto Ran = Fa * (wf.rgrid.r * Fn);
          const auto term = f * (Fn.en - Fa.en) * Ran * Ran / (2 * l + 1);
          if (Fn.n > 0)
            sum_el += term;
          else
            sum_p += term;
        }
        printf("l=%1i, sum = %10.6f%+10.6f = %8.1e\n", l, sum_el, sum_p,
               sum_el + sum_p);
      }
    }

    //----------------------------------------------------------------------------
    std::cout << "\nDrake-Goldman sum rules: w^n |<a|r|b>|^2  (n=0,1,2)\n";
    std::cout << "(Only up to lmax-1, since need to have states with l'=l+1)\n";

    auto comp_ki = [](const auto &Fm, const auto &Fn) {
      return Fm.k_index() < Fn.k_index();
    };
    auto max_ki =
        std::max_element(basis.begin(), basis.end(), comp_ki)->k_index();
    int n_max_DG = 3; // wf.core.empty() ? 3 : 1;
    for (int ki = 0; ki <= max_ki; ki++) {
      auto kappa = Angular::kappaFromIndex(ki);
      auto comp_k = [=](const auto &Fn) { return Fn.k == kappa; };
      auto Fa = *std::find_if(basis.begin(), basis.end(), comp_k);
      // need to have l_n = la+1 terms, or sum doesn't work:
      if (Fa.l() == max_l)
        continue;
      std::cout << "kappa: " << kappa << " (" << Fa.symbol() << ")\n";
      for (int i = 0; i < n_max_DG; i++) {
        auto sum = 0.0;
        for (const auto &Fn : basis) {
          auto w = Fn.en - Fa.en;
          auto Ran = rhat.reducedME(Fa, Fn);
          double c = 1.0 / (2 * std::abs(Fa.k));
          auto term = std::pow(w, i) * Ran * Ran * c;
          sum += term;
        }
        if (i == 2)
          sum *= wf.alpha * wf.alpha / 3;
        auto s0 =
            (i == 0) ? r2hat.radialIntegral(Fa, Fa) : (i == 1) ? 0.0 : 1.0;
        printf("%i: sum=%11.6f, exact=%+11.6f, diff = %8.1e\n", i, sum, s0,
               sum - s0);
      }
    }
  }

  //----------
  const auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto mu = isotope.mu;
  const auto I_nuc = isotope.I_N;
  const auto hfs = DiracOperator::Hyperfine(
      mu, I_nuc, 0.0, wf.rgrid, DiracOperator::Hyperfine::pointlike_F());

  std::cout << "\nHFS and Energies: Basis cf HF:\n";
  std::cout << "    | A(HF)      Basis      eps   | En(HF)      "
               "Basis       eps   | nodes\n";
  int count = 0;
  for (const auto &Fn : basis) {
    if (Fn.n < 0)
      continue;
    const auto *hf_phi = wf.getState(Fn.n, Fn.k);
    const bool hfQ = hf_phi != nullptr;
    const auto Ahf = hfQ ? DiracOperator::Hyperfine::hfsA(&hfs, *hf_phi) : 0.0;
    const auto Ab = DiracOperator::Hyperfine::hfsA(&hfs, Fn);
    const auto Eb = Fn.en;
    const auto Ehf = hfQ ? hf_phi->en : 0.0;

    const auto nodes = Helper::countNodes(Fn);
    const int expected_nodes = Fn.n - Fn.l() - 1;

    if (hfQ) {
      count = 0;
      printf("%4s| %9.3e  %9.3e  %5.0e | ", Fn.shortSymbol().c_str(), Ahf, Ab,
             std::abs((Ahf - Ab) / Ab));
      printf("%10.3e  %10.3e  %5.0e | ", Ehf, Eb, std::abs((Ehf - Eb) / Eb));
      std::cout << nodes << "/" << expected_nodes << "\n";
    } else {
      count++;
      if (count >= 3)
        continue;
      printf("%4s|    ---     %9.3e   ---  | ", Fn.shortSymbol().c_str(), Ab);
      printf("    ---     %10.3e   ---  | ", Eb);
      std::cout << nodes << "/" << expected_nodes << "\n";
    }
  }
}

//******************************************************************************
void Module_test_r0pinf(const Wavefunction &wf) {
  std::cout << "\nTesting boundaries r0 and pinf: f(r)/f_max\n";
  std::cout << " State    f(r0)   f(pinf)   pinf/Rinf\n";
  // for (const auto &phi : wf.core)
  for (const auto tmp_orbs : {&wf.core, &wf.valence}) {
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
  for (int i = 0; i < 9; i++) {
    // const auto &tmp_b = (i == 2) ? wf.valence : wf.core;
    // const auto &tmp_a = (i == 0) ? wf.core : wf.valence;

    const auto &tmp_basis = i < 6 ? wf.basis : wf.spectrum;

    const auto &tmp_b =
        (i == 2 || i == 4 || i == 7)
            ? wf.valence
            : (i == 0 || i == 1 || i == 3 || i == 6) ? wf.core : tmp_basis;
    // core, core, valence, core, valence, basis

    const auto &tmp_a = (i == 0) ? wf.core : (i < 3) ? wf.valence : tmp_basis;
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
      else if (i == 6)
        std::cout << "\nSpectrum-core\n    ";
      else if (i == 7)
        std::cout << "\nSpectrum-Valence\n    ";
      else if (i == 8)
        std::cout << "\nSpectrum-Basis\n    ";
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
      else if (i == 6)
        buffer << "sc ";
      else if (i == 7)
        buffer << "sv ";
      else if (i == 8)
        buffer << "ss ";
    }

    auto worst_xo = 0.0;
    std::string worst_braket = "";
    if (print_all) {
      for (auto &psi_b : tmp_b)
        printf("%2i%2i", psi_b.n, psi_b.k);
      std::cout << "\n";
    }
    for (auto &psi_a : tmp_a) {
      if (psi_a.n < 0)
        continue;
      if (print_all)
        printf("%2i%2i", psi_a.n, psi_a.k);
      for (auto &psi_b : tmp_b) {
        if (psi_b.n < 0)
          continue;
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

  // DirectHamiltonian Hd(wf.vnuc, wf.vdir, wf.alpha);
  auto Hd = RadialHamiltonian(wf.rgrid, wf.alpha);
  Hd.set_v(-1, wf.get_Vlocal(0)); // same each kappa //?? XXX
  Hd.set_v_mag(wf.get_Hmag(0));

  const auto &basis = wf.spectrum.empty() ? wf.basis : wf.spectrum;

  for (const auto tmp_orbs : {&wf.core, &wf.valence, &basis}) {
    if (tmp_orbs->empty())
      continue;
    double worst_eps = 0.0;
    const DiracSpinor *worst_psi = nullptr;
    for (const auto &psi : *tmp_orbs) {
      double Haa_d = Hd.matrixEl(psi, psi);
      double Haa_x = psi * HF::vex_psia_any(psi, wf.core);
      auto Haa = Haa_d + Haa_x;
      // if (!wf.isInCore(psi.n, psi.k) && wf.getSigma() != nullptr) {
      //   Haa += psi * (*wf.getSigma())(psi);
      // }
      if (tmp_orbs != &wf.core && wf.getSigma() != nullptr) {
        Haa += psi * (*wf.getSigma())(psi);
      }
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
