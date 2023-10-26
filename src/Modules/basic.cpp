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
#include "fmt/format.hpp"
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
void testBasis(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"E1", "Calculate E1 matrix elements. Test of large-r, but large "
              "output. [false]"},
       {"scale_rN", "Scale-factor for nuclear radius in hyperfine operator. "
                    "Uses Ball model, with rN given by charge radius. Set to "
                    "zero for pointlike. [1.0]"},
       {"sum_rules", "Do TKR and DG sum rules - only accurate if -ve energy "
                     "states included [false]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto do_E1 = input.get("E1", false);
  const auto scale_rN = input.get("scale_rN", 1.0);
  const auto sum_rules = input.get("sum_rules", false);

  // 0. Compare orthonormality
  auto comp_ortho = [&](const std::vector<DiracSpinor> &orbitals,
                        const std::vector<DiracSpinor> &basis) {
    fmt::print("a    <a|a>    b    <a|b>\n");
    for (auto &c : orbitals) {
      const auto b = DiracSpinor::find(c.n(), c.kappa(), basis);
      if (!b)
        continue;
      const auto norm = std::abs(*b * c - 1.0);
      double eps = 0.0;
      std::string worst;
      for (const auto &n : basis) {
        if (n.kappa() != c.kappa() || n.n() == c.n())
          continue;
        auto teps = std::abs(c * n);
        if (teps > eps) {
          eps = teps;
          worst = n.shortSymbol();
        }
      }
      fmt::print("{:4s} {:.1e}  {:4s} {:.1e}\n", c.shortSymbol(), norm, worst,
                 eps);
    }
  };

  // 1. Compare energies (and <r>) to HF orbitals
  auto comp_energies = [&](const std::vector<DiracSpinor> &orbitals,
                           const std::vector<DiracSpinor> &basis) {
    fmt::print("      <r>_HF       <r>_b         eps     E_HF          E_b   "
               "        eps\n");
    for (auto &c : orbitals) {
      const auto b = DiracSpinor::find(c.n(), c.kappa(), basis);
      if (b) {
        const auto r_c = c * (wf.grid().r() * c);
        const auto r_b = *b * (wf.grid().r() * *b);
        const auto eps_r = r_b / r_c - 1.0;
        // const auto norm = std::abs(*b * c - 1.0);
        const auto eps = b->en() / c.en() - 1.0;
        fmt::print("{:4s}  {:5e} {:5e} {:+.0e}  {:5e} {:5e} {:+.0e}\n",
                   c.shortSymbol(), r_c, r_b, eps_r, c.en(), b->en(), eps);
      }
    }
  };

  // 2. Compare E1
  const DiracOperator::E1 d(wf.grid());
  auto comp_E1 = [&](const std::vector<DiracSpinor> &orbitals,
                     const std::vector<DiracSpinor> &basis) {
    for (auto &a : orbitals) {
      for (auto &b : orbitals) {
        if (b > a || d.isZero(a, b))
          continue;

        const auto n = DiracSpinor::find(a.n(), a.kappa(), basis);
        const auto m = DiracSpinor::find(b.n(), b.kappa(), basis);
        if (!n || !m)
          continue;

        const auto dab = d.reducedME(a, b);
        const auto dnm = d.reducedME(*n, *m);
        const auto eps = dnm / dab - 1.0;

        fmt::print("<{:3s}||d||{:3s}> {:+.4e} {:+.4e} {:+.0e}\n",
                   a.shortSymbol(), b.shortSymbol(), dab, dnm, eps);
      }
    }
  };

  // 3. Compare HFS

  const auto rN_fm = scale_rN * std::sqrt(5.0 / 3) * wf.get_rrms();
  const auto rN_au = rN_fm / PhysConst::aB_fm;

  DiracOperator::hfs hfs(1, 1.0, rN_au, wf.grid(),
                         DiracOperator::Hyperfine::sphericalBall_F());
  auto comp_HFS = [&](const std::vector<DiracSpinor> &orbitals,
                      const std::vector<DiracSpinor> &basis) {
    for (auto &a : orbitals) {

      const auto n = DiracSpinor::find(a.n(), a.kappa(), basis);
      if (!n)
        continue;
      const auto A_const =
          DiracOperator::Hyperfine::convert_RME_to_AB(1, a.kappa(), a.kappa());
      const auto Aa = hfs.reducedME(a, a) * A_const;
      const auto An = hfs.reducedME(*n, *n) * A_const;
      const auto eps = Aa / An - 1.0;

      fmt::print("{:4s} {:+.4e} {:+.4e} {:+.0e}\n", a.shortSymbol(), Aa, An,
                 eps);
    }
  };

  // Test with Basis
  if (!wf.basis().empty()) {
    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Basis:\n";

    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Orthonormality: <a_HF|a_basis>, and worst <a_HF|b_basis>:\n";
    std::cout << "\nOrhtonromality: Basis vs core HF orbitals:\n";
    comp_ortho(wf.core(), wf.basis());
    std::cout << "\nOrhtonromality: Basis vs valence HF orbitals:\n";
    if (wf.Sigma()) {
      std::cout
          << "Note: Bruckner valence orbitals, not expected to be equal:\n";
    }
    comp_ortho(wf.valence(), wf.basis());

    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Energies and <r> expectation values:\n";

    std::cout << "\nEnergies: Basis vs core HF orbitals:\n";
    comp_energies(wf.core(), wf.basis());
    std::cout << "\nEnergies: Basis vs valence HF orbitals:\n";
    if (wf.Sigma()) {
      std::cout
          << "Note: Bruckner valence orbitals, not expected to be equal:\n";
    }
    comp_energies(wf.valence(), wf.basis());

    if (do_E1) {
      std::cout
          << "\n----------------------------------------------------------\n";
      std::cout << "E1 reduced ME (large r test):\n";

      std::cout << "\nE1: Basis vs core HF orbitals:\n";
      comp_E1(wf.core(), wf.basis());
      std::cout << "\nE1: Basis vs valence HF orbitals:\n";
      if (wf.Sigma()) {
        std::cout
            << "Note: Bruckner valence orbitals, not expected to be equal:\n";
      }
      comp_E1(wf.valence(), wf.basis());
    }

    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Hyperfine A (ball model, g=1) (small r test):\n";
    std::cout << "rN = " << rN_fm << " fm\n";

    std::cout << "\nHFS: Basis vs core HF orbitals:\n";
    comp_HFS(wf.core(), wf.basis());
    std::cout << "\nHFS: Basis vs valence HF orbitals:\n";
    if (wf.Sigma()) {
      std::cout
          << "Note: Bruckner valence orbitals, not expected to be equal:\n";
    }
    comp_HFS(wf.valence(), wf.basis());

    if (sum_rules) {
      std::cout
          << "\n----------------------------------------------------------\n";
      std::cout << "Sum rules:\n";
      std::cout << "(Only work if -ve energy states are included!)\n";

      std::cout << "\nTKR sum rule:\n";
      SplineBasis::sumrule_TKR(wf.basis(), wf.grid().r(), true);

      std::cout << "\nDrake-Gordon sum rule:\n";
      for (int nDG = 0; nDG < 3; ++nDG) {
        SplineBasis::sumrule_DG(nDG, wf.basis(), wf.grid(), wf.alpha(), true);
      }
    }
  }

  //====================================================================

  if (!wf.spectrum().empty()) {
    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "----------------------------------------------------------\n";
    std::cout << "Spectrum:\n";

    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Orthonormality: <a_HF|a_spect>, and worst <a_HF|b_spect>:\n";
    std::cout << "\nOrhtonromality: Spectrum vs core HF orbitals:\n";
    if (wf.Sigma()) {
      std::cout
          << "Note: Correlations in spectrum, not expected to be equal:\n";
    }
    comp_ortho(wf.core(), wf.spectrum());
    std::cout << "\nOrhtonromality: Spectrum vs valence HF orbitals:\n";
    comp_ortho(wf.valence(), wf.spectrum());

    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Energies and <r> expectation values:\n";

    std::cout << "\nEnergies: Spectrum vs core HF orbitals:\n";
    if (wf.Sigma()) {
      std::cout
          << "Note: Correlations in spectrum, not expected to be equal:\n";
    }
    comp_energies(wf.core(), wf.spectrum());
    std::cout << "\nEnergies: Spectrum vs valence HF orbitals:\n";
    comp_energies(wf.valence(), wf.spectrum());

    if (do_E1) {
      std::cout
          << "\n----------------------------------------------------------\n";
      std::cout << "E1 reduced ME (large r test):\n";

      std::cout << "\nE1: Spectrum vs core HF orbitals:\n";
      if (wf.Sigma()) {
        std::cout
            << "Note: Correlations in spectrum, not expected to be equal:\n";
      }
      comp_E1(wf.core(), wf.spectrum());
      std::cout << "\nE1: Spectrum vs valence HF orbitals:\n";
      comp_E1(wf.valence(), wf.spectrum());
    }

    std::cout
        << "\n----------------------------------------------------------\n";
    std::cout << "Hyperfine A (ball model, g=1) (small r test):\n";
    std::cout << "rN = " << rN_fm << " fm\n";

    std::cout << "\nHFS: Spectrum vs core HF orbitals:\n";
    if (wf.Sigma()) {
      std::cout
          << "Note: Correlations in spectrum, not expected to be equal:\n";
    }
    comp_HFS(wf.core(), wf.spectrum());
    std::cout << "\nHFS: Spectrum vs valence HF orbitals:\n";
    comp_HFS(wf.valence(), wf.spectrum());

    if (sum_rules) {
      std::cout
          << "\n----------------------------------------------------------\n";
      std::cout << "Sum rules:\n";
      std::cout << "(Only work if -ve energy states are included!)\n";

      std::cout << "\nTKR sum rule:\n";
      SplineBasis::sumrule_TKR(wf.spectrum(), wf.grid().r(), true);

      std::cout << "\nDrake-Gordon sum rule:\n";
      for (int nDG = 0; nDG < 3; ++nDG) {
        SplineBasis::sumrule_DG(nDG, wf.spectrum(), wf.grid(), wf.alpha(),
                                true);
      }
    }
  }
}

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
  input.check({{"label", "Optional label for output text file"},
               {"l", "If set, will only write out for specified l. If not set, "
                     "will write out all"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  std::cout << "\n Running: " << ThisModule << "\n";
  const auto label = input.get<std::string>("label", "");
  // to write only for specific l. l<0 means all
  auto l = input.get("l", -1);

  std::string oname = wf.identity();
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

  input.check(
      {{"orthonormal", "Checks orthonormalityies [true]"},
       {"Hamiltonian", "Compare eigen energies to Hamiltonian matrix elements"},
       {"boundaries", "Check boundaries: f(r0)/f_max, f(rmax)/fmax"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  auto othon = input.get("orthonormal", true);
  if (othon)
    orthonormality(wf, true);
  if (input.get("Hamiltonian", true))
    Hamiltonian(wf);
  if (input.get("boundaries", true))
    r0pinf(wf);
}

//==============================================================================
//==============================================================================
namespace Tests {

//==============================================================================
void r0pinf(const Wavefunction &wf) {
  std::cout << "\nTesting boundaries r0 and pinf: R = f(r)/f_max\n";
  std::cout << " State    R(r0)   R(pinf)   pinf/Rinf\n";
  // for (const auto &phi : wf.core())
  for (const auto tmp_orbs : {&wf.core(), &wf.valence()}) {
    for (const auto &phi : *tmp_orbs) {
      auto ratios = phi.r0pinfratio();
      printf("%7s:  %.0e   %.0e   %5i/%6.2f\n", phi.symbol().c_str(),
             std::abs(ratios.first), std::abs(ratios.second), (int)phi.max_pt(),
             wf.grid().r()[phi.max_pt() - 1]);
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
    }
  }
}

//------------------------------------------------------------------------------
void Hamiltonian(const Wavefunction &wf) {
  std::cout << "\nTesting wavefunctions: <n|H|n>  (numerical error)\n";

  for (const auto tmp_orbs : {&wf.core(), &wf.valence()}) {
    if (tmp_orbs->empty())
      continue;
    double worst_eps = 0.0;
    const DiracSpinor *worst_Fn = nullptr;
    for (const auto &Fn : *tmp_orbs) {
      double Haa_d = wf.H0ab(Fn, Fn);
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
    std::cout << "--------------\n";
  }
}

} // namespace Tests

} // namespace Module