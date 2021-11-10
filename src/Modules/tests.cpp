#include "Modules/tests.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace Module {

//******************************************************************************
void Module_tests(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace Tests;
  std::string ThisModule = "Module::Tests";
  input.checkBlock_old({"orthonormal", "orthonormal_all", "Hamiltonian",
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
  double sp = Fn.f(Fn.min_pt() + 3);
  int counted_nodes = 0;
  for (auto i = Fn.min_pt() + 4; i < Fn.max_pt() - 3; ++i) {
    if (sp * Fn.f(i) < 0) {
      ++counted_nodes;
      sp = Fn.f(i);
    }
  }
  return counted_nodes;
}
} // namespace Helper

//------------------------------------------------------------------------------
void basisTests(const Wavefunction &wf) {

  std::cout << "\nTesting basis/spectrum:\n";

  const auto &basis = wf.spectrum.empty() ? wf.basis : wf.spectrum;
  if (basis.empty())
    return;

  if (&basis == &(wf.spectrum))
    std::cout << "Using Sprectrum\n";
  else
    std::cout << "Using Basis\n";

  //----------
  const auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto mu = isotope.mu;
  const auto I_nuc = isotope.I_N;
  const auto hfs = DiracOperator::HyperfineA(
      mu, I_nuc, 0.0, *(wf.rgrid), DiracOperator::Hyperfine::pointlike_F());

  std::cout << "\nHFS and Energies: Basis cf HF:\n";
  std::cout << "    | A(HF)      Basis      eps   | En(HF)      "
               "Basis       eps   | norm\n"; // nodes\n";
  int count = 0;
  for (const auto &Fn : basis) {
    if (Fn.n < 0)
      continue;
    const auto *hf_phi = wf.getState(Fn.n, Fn.k);
    const bool hfQ = hf_phi != nullptr;
    const auto Ahf = hfQ ? DiracOperator::HyperfineA::hfsA(&hfs, *hf_phi) : 0.0;
    const auto Ab = DiracOperator::HyperfineA::hfsA(&hfs, Fn);
    const auto Eb = Fn.en();
    const auto Ehf = hfQ ? hf_phi->en() : 0.0;

    // const auto nodes = Helper::countNodes(Fn);
    // const int expected_nodes = Fn.n - Fn.l() - 1;

    if (hfQ) {
      count = 0;
      printf("%4s| %9.3e  %9.3e  %5.0e | ", Fn.shortSymbol().c_str(), Ahf, Ab,
             std::abs((Ahf - Ab) / Ab));
      printf("%10.3e  %10.3e  %5.0e | ", Ehf, Eb, std::abs((Ehf - Eb) / Eb));
      printf("%.0e", *hf_phi * Fn - 1.0);
    } else {
      count++;
      if (count >= 3)
        continue;
      printf("%4s|    ---     %9.3e   ---  | ", Fn.shortSymbol().c_str(), Ab);
      printf("    ---     %10.3e   ---  | ", Eb);
    }
    std::cout /*<< nodes << "/" << expected_nodes*/ << "\n";
  }

  std::cout << "\nCompleteness test:\n";
  std::cout << "Sum_n <a|r|n><n|1/r|a>  <a|r|n><n|r|a>\n";
  std::cout << "vs:   <a|a>             <a|r^2|a>\n";
  for (const auto orbs : {/*&wf.core,*/ &wf.valence}) {
    for (const auto &Fa : *orbs) {
      auto [e1, er2] = SplineBasis::r_completeness(Fa, basis, *wf.rgrid);
      printf("%4s   %10.2e         %10.2e\n", Fa.shortSymbol().c_str(), e1,
             er2);
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
             std::abs(ratios.first), std::abs(ratios.second), (int)phi.max_pt(),
             wf.rgrid->r()[phi.max_pt() - 1]);
      // std::cout << ratios.first << " " << ratios.second << "\n";
    }
    std::cout << "--------------\n";
  }
}

//------------------------------------------------------------------------------
void Module_Tests_orthonormality(const Wavefunction &wf, const bool) {
  std::cout << "\nTest orthonormality:\n";

  const std::vector orbs = {&wf.core, &wf.valence, &wf.basis, &wf.spectrum};
  const std::vector names = {'c', 'v', 'b', 's'};

  for (auto i = 0ul; i < orbs.size(); ++i) {
    if (orbs[i]->empty())
      continue;
    for (auto j = i; j < orbs.size(); ++j) {
      if (orbs[j]->empty())
        continue;
      const auto [eps, str] = DiracSpinor::check_ortho(*orbs[i], *orbs[j]);
      std::cout << names[i] << names[j] << " ";
      printf("%11s = %.1e\n", str.c_str(), eps);
      // std::cout << std::left << std::setw(11) << str << " = ";
      // std::cout << std::setprecision(1) << std::scientific << eps << "\n";
    }
  }
  // std::cout.flags(f);
}

//------------------------------------------------------------------------------
void Module_Tests_Hamiltonian(const Wavefunction &wf) {
  std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";

  // XXX Add Breit!

  // auto Hd = RadialHamiltonian(wf.rgrid, wf.alpha);
  // Hd.set_v(-1, wf.get_Vlocal(0)); // same each kappa //?? XXX
  // Hd.set_v_mag(wf.get_Hmag(0));

  const auto &basis = wf.spectrum.empty() ? wf.basis : wf.spectrum;

  for (const auto tmp_orbs : {&wf.core, &wf.valence, &basis}) {
    if (tmp_orbs->empty())
      continue;
    double worst_eps = 0.0;
    const DiracSpinor *worst_Fn = nullptr;
    for (const auto &Fn : *tmp_orbs) {
      double Haa_d = wf.Hab(Fn, Fn);
      double Haa_x = Fn * HF::vexFa(Fn, wf.core);
      auto Haa = Haa_d + Haa_x;
      if (tmp_orbs != &wf.core && wf.getSigma() != nullptr) {
        Haa += Fn * (*wf.getSigma())(Fn);
      }
      double ens = Fn.en();
      double fracdiff = (Haa - ens) / ens;
      printf("<%2i% i|H|%2i% i> = %17.11f, E = %17.11f; % .0e\n", Fn.n, Fn.k,
             Fn.n, Fn.k, Haa, ens, fracdiff);
      if (std::abs(fracdiff) >= std::abs(worst_eps)) {
        worst_eps = fracdiff;
        worst_Fn = &Fn;
      }
    }
    if (worst_Fn != nullptr)
      std::cout << worst_Fn->symbol() << ": eps=" << worst_eps << "\n";
  }
}

} // namespace Tests

} // namespace Module
