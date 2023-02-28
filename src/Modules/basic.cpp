#include "basic.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/modules_list.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Module {

//==============================================================================
static void write_orbitals(const std::string &fname,
                           const std::vector<DiracSpinor> &orbs, int l) {
  if (orbs.empty())
    return;
  const auto &gr = orbs.front().grid();

  std::ofstream of(fname);
  of << "r ";
  for (auto &psi : orbs) {
    if (psi.l() != l && l >= 0)
      continue;
    of << "\"" << psi.symbol(true) << "\" ";
  }
  of << "\n";

  of << "# f block\n";
  for (std::size_t i = 0; i < gr.num_points(); i++) {
    of << gr.r(i) << " ";
    for (auto &psi : orbs) {
      if (psi.l() != l && l >= 0)
        continue;
      of << psi.f(i) << " ";
    }
    of << "\n";
  }

  of << "\n# g block\n";
  for (std::size_t i = 0; i < gr.num_points(); i++) {
    of << gr.r(i) << " ";
    for (auto &psi : orbs) {
      if (psi.l() != l && l >= 0)
        continue;
      of << psi.g(i) << " ";
    }
    of << "\n";
  }

  of.close();
  std::cout << "Orbitals written to file: " << fname << "\n";
}

//==============================================================================
void writeOrbitals(const IO::InputBlock &input, const Wavefunction &wf) {
  const std::string ThisModule = "Module::WriteOrbitals";
  input.check({{"label", ""}, {"l", ""}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  std::cout << "\n Running: " << ThisModule << "\n";
  const auto label = input.get<std::string>("label", "");
  // to write only for specific l. l<0 means all
  auto l = input.get("l", -1);

  std::string oname = wf.atomicSymbol();
  if (label != "")
    oname += "_" + label;

  write_orbitals(oname + "_core.txt", wf.core(), l);
  write_orbitals(oname + "_valence.txt", wf.valence(), l);
  write_orbitals(oname + "_basis.txt", wf.basis(), l);
  write_orbitals(oname + "_spectrum.txt", wf.spectrum(), l);
}

//==============================================================================
void continuum(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"energy", "List. energy for cntm states (>0) [0.5]"},
       {"max_l", "maximum l. Will calculate orbital for each energy and each l "
                 "(and j) [0]"},
       {"filename",
        "filename to output continuum orbitals. If blank, will not write"},
       {"operator",
        "Operator to calculate matrix elements (e.g., E1) [blank by default]"},
       {"options", "options specific to operator [blank by default]"},
       {"rpa", "Include RPA (TDHF for now)? [false]"},
       {"omega", "Frequency for RPA [0.0]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto en_list = input.get("energy", std::vector{0.5});
  const auto lmax = input.get("max_l", 0);
  const auto fname = input.get("filename", std::string{""});

  // Method options for solveContinuumHF
  // auto force_rescale = input.get<bool>("force_rescale", false);
  // auto subtract_self = input.get<bool>("subtract_self", false);
  // auto force_orthog = input.get<bool>("force_orthog", false);

  // n=0 invalid -> nullptr
  auto n = input.get<int>("n", 0);
  auto k = input.get<int>("kappa", 0);
  // auto p_psi = wf.getState(n, k);
  auto p_psi = (n == 0) ? nullptr : wf.getState(n, k);

  if (p_psi != nullptr) {
    std::cout << "\n"
              << "State for n = " << n << " and kappa = " << k << ": "
              << p_psi->symbol() << "\n";
    // std::cout << "Bound state: " << p_psi->symbol() << "\n";
  }

  std::cout << "\nContinuum Orbitals:\n";
  auto cntm = ContinuumOrbitals(wf);

  std::cout << "For l up to l=" << lmax << "\n";
  std::cout << "At energy: ";
  for (const auto en_c : en_list) {
    std::cout << en_c << ", ";
    cntm.solveContinuumHF(en_c, 0, lmax, nullptr, false, false, false);
  }
  std::cout << "\n";

  //-----------------------------------------------
  std::cout << "\nCheck orthogonanilty:\n";
  cntm.check_orthog(true);

  // Matrix elements:
  //-----------------------------------------------
  const auto oper = input.get<std::string>("operator", "");
  if (oper != "") {

    // Get optional 'options' for operator
    const auto h_options = input.getBlock("options");
    const auto h = DiracOperator::generate(
        oper, h_options ? *h_options : IO::InputBlock{}, wf);

    std::cout << "\nMatrix elements of " << h->name() << "\n";

    bool eachFreqQ = false;

    const auto rpaQ = input.get("rpa", false);
    const auto omega = input.get("omega", 0.0);

    auto rpa = ExternalField::TDHF(h.get(), wf.vHF());
    const auto p_rpa = rpaQ ? &rpa : nullptr;

    const auto mes_c = ExternalField::calcMatrixElements(
        wf.core(), cntm.orbitals, h.get(), p_rpa, omega, eachFreqQ);

    const auto mes_v = ExternalField::calcMatrixElements(
        wf.valence(), cntm.orbitals, h.get(), p_rpa, omega, eachFreqQ);

    std::cout << (rpaQ ? ExternalField::MEdata::title() :
                         ExternalField::MEdata::title_noRPA())
              << "\n";
    if (!wf.core().empty()) {
      std::cout << "Core:\n";
    }
    for (const auto &me : mes_c) {
      std::cout << me << "\n";
    }
    if (!wf.valence().empty()) {
      std::cout << "Valence:\n";
    }
    for (const auto &me : mes_v) {
      std::cout << me << "\n";
    }
  }

  //-----------------------------------------------
  // Write continuum wavefunctions to file
  if (fname != "") {
    std::cout << "\nWriting orbitals to " << fname << "\n";
    std::ofstream of(fname);
    const auto &gr = wf.grid();
    of << "r ";
    for (const auto &Fc : cntm.orbitals) {
      of << "\"f" << Fc.symbol(true) << "(" << Fc.en() << ")\" ";
    }
    for (const auto &Fc : cntm.orbitals) {
      of << "\"g" << Fc.symbol(true) << "(" << Fc.en() << ")\" ";
    }
    of << "\n";
    for (std::size_t i = 0; i < gr.num_points(); i++) {
      of << gr.r(i) << " ";
      for (const auto &Fc : cntm.orbitals) {
        of << Fc.f(i) << " ";
      }
      for (const auto &Fc : cntm.orbitals) {
        of << Fc.g(i) << " ";
      }
      of << "\n";
    }
  }
}

//==============================================================================
void tests(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace Tests;
  std::string ThisModule = "Module::Tests";
  input.check({{"orthonormal", ""},
               {"orthonormal_all", ""},
               {"Hamiltonian", ""},
               {"boundaries", ""},
               {"basis", ""}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  auto othon = input.get("orthonormal", true);
  auto othon_all = input.get("orthonormal_all", false);
  if (othon || othon_all)
    orthonormality(wf, othon_all);
  if (input.get("Hamiltonian", false))
    Hamiltonian(wf);
  if (input.get("boundaries", false))
    r0pinf(wf);
  if (input.get("basis", false))
    basis(wf);
}

//==============================================================================
//==============================================================================
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
void basis(const Wavefunction &wf) {

  std::cout << "\nTesting basis/spectrum:\n";

  const auto &basis = wf.spectrum().empty() ? wf.basis() : wf.spectrum();
  if (basis.empty())
    return;

  if (&basis == &(wf.spectrum()))
    std::cout << "Using Sprectrum\n";
  else
    std::cout << "Using Basis\n";

  //----------
  const auto isotope = Nuclear::findIsotopeData(wf.Znuc(), wf.Anuc());
  const auto mu = isotope.mu;
  const auto I_nuc = isotope.I_N;
  const auto hfs = DiracOperator::hfs(1, mu / I_nuc, 0.0, wf.grid(),
                                      DiracOperator::Hyperfine::pointlike_F());

  std::cout << "\nHFS and Energies: Basis cf HF:\n";
  std::cout << "    | A(HF)      Basis      eps   | En(HF)      "
               "Basis       eps   | norm\n"; // nodes\n";
  int count = 0;
  for (const auto &Fn : basis) {
    if (Fn.n() < 0)
      continue;
    const auto *hf_phi = wf.getState(Fn.n(), Fn.kappa());
    const bool hfQ = hf_phi != nullptr;
    const auto Ahf = hfQ ? DiracOperator::Hyperfine::hfsA(&hfs, *hf_phi) : 0.0;
    const auto Ab = DiracOperator::Hyperfine::hfsA(&hfs, Fn);
    const auto Eb = Fn.en();
    const auto Ehf = hfQ ? hf_phi->en() : 0.0;

    // const auto nodes = Helper::countNodes(Fn);
    // const int expected_nodes = Fn.n() - Fn.l() - 1;

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
  for (const auto orbs : {/*&wf.core(),*/ &wf.valence()}) {
    for (const auto &Fa : *orbs) {
      auto [e1, er2] = SplineBasis::r_completeness(Fa, basis, wf.grid());
      printf("%4s   %10.2e         %10.2e\n", Fa.shortSymbol().c_str(), e1,
             er2);
    }
  }
}

//==============================================================================
void r0pinf(const Wavefunction &wf) {
  std::cout << "\nTesting boundaries r0 and pinf: f(r)/f_max\n";
  std::cout << " State    f(r0)   f(pinf)   pinf/Rinf\n";
  // for (const auto &phi : wf.core())
  for (const auto tmp_orbs : {&wf.core(), &wf.valence()}) {
    for (const auto &phi : *tmp_orbs) {
      auto ratios = phi.r0pinfratio();
      printf("%7s:  %.0e   %.0e   %5i/%6.2f\n", phi.symbol().c_str(),
             std::abs(ratios.first), std::abs(ratios.second), (int)phi.max_pt(),
             wf.grid().r()[phi.max_pt() - 1]);
      // std::cout << ratios.first << " " << ratios.second << "\n";
    }
    std::cout << "--------------\n";
  }
}

//------------------------------------------------------------------------------
void orthonormality(const Wavefunction &wf, const bool) {
  std::cout << "\nTest orthonormality:\n";

  const std::vector orbs = {&wf.core(), &wf.valence(), &wf.basis(),
                            &wf.spectrum()};
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
void Hamiltonian(const Wavefunction &wf) {
  std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";

  // XXX Add Breit!

  // auto Hd = RadialHamiltonian(wf.grid_sptr(), wf.alpha());
  // Hd.set_v(-1, wf.vlocal(0)); // same each kappa //?? XXX
  // Hd.set_v_mag(wf.Hmag(0));

  const auto &basis = wf.spectrum().empty() ? wf.basis() : wf.spectrum();

  for (const auto tmp_orbs : {&wf.core(), &wf.valence(), &basis}) {
    if (tmp_orbs->empty())
      continue;
    double worst_eps = 0.0;
    const DiracSpinor *worst_Fn = nullptr;
    for (const auto &Fn : *tmp_orbs) {
      double Haa_d = wf.Hab(Fn, Fn);
      double Haa_x = Fn * HF::vexFa(Fn, wf.core());
      auto Haa = Haa_d + Haa_x;
      if (tmp_orbs != &wf.core() && wf.Sigma() != nullptr) {
        Haa += Fn * (*wf.Sigma())(Fn);
      }
      double ens = Fn.en();
      double fracdiff = (Haa - ens) / ens;
      printf("<%2i% i|H|%2i% i> = %17.11f, E = %17.11f; % .0e\n", Fn.n(),
             Fn.kappa(), Fn.n(), Fn.kappa(), Haa, ens, fracdiff);
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