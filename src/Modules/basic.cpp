#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
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

/*!
  @brief Runs basic numerical tests on HF orbitals and basis/spectrum.
  @details Tests orthonormality, Hamiltonian expectation values, radial
  boundaries, and basis/spectrum comparisons (energies, E1, HFS, sum rules).
  Basis/spectrum tests are skipped if the corresponding sets are empty.
*/
void tests(const IO::InputBlock &input, const Wavefunction &wf);

//! Writes all orbital sets to disk for plotting.
void writeOrbitals(const IO::InputBlock &input, const Wavefunction &wf);

/*!
  @brief Computes continuum orbitals and tests orthogonality.
  @details Solves for continuum states at each specified energy for l = 0..max_l.
  Optionally computes matrix elements with a chosen operator, with or without RPA.
  Orbitals can be written to file for plotting.
*/
void continuum(const IO::InputBlock &input, const Wavefunction &wf);

// Register the modules
namespace {
const Register r_tests{"tests", "Basic wavefunction and basis numerical tests",
                       &tests};
const Register r_writeOrbitals{
  "writeOrbitals", "Write orbitals to disk for plotting", &writeOrbitals};
const Register r_continuum{
  "continuum", "Compute and use continuum wavefunctions", &continuum};
} // namespace

//==============================================================================
// Sub-functions for "tests" (defined below)
namespace Tests {

// Tests orthonormality of core and valence orbitals
void orthonormality(const Wavefunction &wf);
// Tests Hamiltonian expectation values against stored energies
void Hamiltonian(const Wavefunction &wf);
// Tests boundary conditions: psi(r0) and psi(r_inf) ~ 0
void r0pinf(const Wavefunction &wf);
// Tests orthonormality of basis against orbitals
void basisOrtho(const std::vector<DiracSpinor> &orbs,
                const std::vector<DiracSpinor> &basis);
// Compares basis energies against HF orbital energies
void basisEnergies(const std::vector<DiracSpinor> &orbs,
                   const std::vector<DiracSpinor> &basis,
                   const std::vector<double> &r);
// Compares E1 matrix elements between basis and HF orbitals
void basisE1(const std::vector<DiracSpinor> &orbs,
             const std::vector<DiracSpinor> &basis, const DiracOperator::E1 &d);
// Compares HFS matrix elements between basis and HF orbitals
void basisHFS(const std::vector<DiracSpinor> &orbs,
              const std::vector<DiracSpinor> &basis,
              const DiracOperator::hfs &hfs);
// Tests basis completeness via sum rules
void sumRules(const std::vector<DiracSpinor> &basis, const Grid &grid,
              double alpha);
// Runs full suite of basis tests (ortho, energies, E1, HFS, sum rules)
void compareBasisSet(const std::string &label, const Wavefunction &wf,
                     const std::vector<DiracSpinor> &basis, bool do_E1,
                     double scale_rN, bool sum_rules,
                     const std::string &sigma_note,
                     bool sigma_note_before_core);

} // namespace Tests

//==============================================================================
void tests(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace Tests;

  input.check(
    {{"orthonormal", "Checks orthonormality [true]"},
     {"Hamiltonian",
      "Compare eigen energies to Hamiltonian matrix elements [true]"},
     {"boundaries", "Check boundaries: f(r0)/f_max, f(rmax)/fmax [true]"},
     {"E1", "Calculate E1 matrix elements for basis/spectrum. "
            "Test of large-r, but large output. [false]"},
     {"scale_rN", "Scale-factor for nuclear radius in hyperfine operator. "
                  "Uses ball model with rN from charge radius. "
                  "Set to 0 for point-like. [1.0]"},
     {"sum_rules", "Do TKR and DG sum rules (only accurate if -ve energy "
                   "states included) [false]"}});
  if (input.has_option("help"))
    return;

  if (input.get("orthonormal", true))
    orthonormality(wf);
  if (input.get("Hamiltonian", true))
    Hamiltonian(wf);
  if (input.get("boundaries", true))
    r0pinf(wf);

  const auto do_E1 = input.get("E1", false);
  const auto scale_rN = input.get("scale_rN", 1.0);
  const auto sum_rules = input.get("sum_rules", false);

  if (!wf.basis().empty())
    compareBasisSet(
      "Basis", wf, wf.basis(), do_E1, scale_rN, sum_rules,
      "Note: Bruckner valence orbitals, not expected to be equal:", false);
  if (!wf.spectrum().empty())
    compareBasisSet(
      "Spectrum", wf, wf.spectrum(), do_E1, scale_rN, sum_rules,
      "Note: Correlations in spectrum, not expected to be equal:", true);
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

void writeOrbitals(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({{"label", "Optional label for output text file"},
               {"l", "If set, will only write out for specified l. If not set, "
                     "will write out all"}});
  if (input.has_option("help"))
    return;

  std::cout << "\n Running: Module::WriteOrbitals\n";
  const auto label = input.get<std::string>("label", "");
  auto l = input.get("l", -1); // l<0 means all

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
  if (input.has_option("help"))
    return;

  const auto en_list = input.get("energy", std::vector{0.5});
  const auto lmax = input.get("max_l", 0);
  const auto fname = input.get("filename", std::string{""});

  // n=0 invalid -> nullptr
  auto n = input.get<int>("n", 0);
  auto k = input.get<int>("kappa", 0);
  auto p_psi = (n == 0) ? nullptr : wf.getState(n, k);

  if (p_psi != nullptr) {
    std::cout << "\nState for n = " << n << " and kappa = " << k << ": "
              << p_psi->symbol() << "\n";
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

  std::cout << "\nCheck orthogonanilty:\n";
  cntm.check_orthog(true);

  // Matrix elements
  const auto oper = input.get<std::string>("operator", "");
  if (oper != "") {

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
    if (!wf.core().empty())
      std::cout << "Core:\n";
    for (const auto &me : mes_c)
      std::cout << me << "\n";
    if (!wf.valence().empty())
      std::cout << "Valence:\n";
    for (const auto &me : mes_v)
      std::cout << me << "\n";
  }

  // Write continuum wavefunctions to file
  if (fname != "") {
    std::cout << "\nWriting orbitals to " << fname << "\n";
    std::ofstream of(fname);
    const auto &gr = wf.grid();
    of << "r ";
    for (const auto &Fc : cntm.orbitals)
      of << "\"f" << Fc.symbol(true) << "(" << Fc.en() << ")\" ";
    for (const auto &Fc : cntm.orbitals)
      of << "\"g" << Fc.symbol(true) << "(" << Fc.en() << ")\" ";
    of << "\n";
    for (std::size_t i = 0; i < gr.num_points(); i++) {
      of << gr.r(i) << " ";
      for (const auto &Fc : cntm.orbitals)
        of << Fc.f(i) << " ";
      for (const auto &Fc : cntm.orbitals)
        of << Fc.g(i) << " ";
      of << "\n";
    }
  }
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
void orthonormality(const Wavefunction &wf) {
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

//------------------------------------------------------------------------------
void basisOrtho(const std::vector<DiracSpinor> &orbs,
                const std::vector<DiracSpinor> &basis) {
  fmt::print("a    <a|a>    b    <a|b>\n");
  for (const auto &c : orbs) {
    const auto b = DiracSpinor::find(c.n(), c.kappa(), basis);
    if (!b)
      continue;
    const auto norm = std::abs(*b * c - 1.0);
    double eps = 0.0;
    std::string worst;
    for (const auto &n : basis) {
      if (n.kappa() != c.kappa() || n.n() == c.n())
        continue;
      const auto teps = std::abs(c * n);
      if (teps > eps) {
        eps = teps;
        worst = n.shortSymbol();
      }
    }
    fmt::print("{:4s} {:.1e}  {:4s} {:.1e}\n", c.shortSymbol(), norm, worst,
               eps);
  }
}

//------------------------------------------------------------------------------
void basisEnergies(const std::vector<DiracSpinor> &orbs,
                   const std::vector<DiracSpinor> &basis,
                   const std::vector<double> &r) {
  fmt::print("      <r>_HF       <r>_b         eps     E_HF          E_b   "
             "        eps\n");
  for (const auto &c : orbs) {
    const auto b = DiracSpinor::find(c.n(), c.kappa(), basis);
    if (!b)
      continue;
    const auto r_c = c * (r * c);
    const auto r_b = *b * (r * *b);
    const auto eps_r = r_b / r_c - 1.0;
    const auto eps = b->en() / c.en() - 1.0;
    fmt::print("{:4s}  {:5e} {:5e} {:+.0e}  {:5e} {:5e} {:+.0e}\n",
               c.shortSymbol(), r_c, r_b, eps_r, c.en(), b->en(), eps);
  }
}

//------------------------------------------------------------------------------
void basisE1(const std::vector<DiracSpinor> &orbs,
             const std::vector<DiracSpinor> &basis,
             const DiracOperator::E1 &d) {
  for (const auto &a : orbs) {
    for (const auto &b : orbs) {
      if (b > a || d.isZero(a, b))
        continue;
      const auto n = DiracSpinor::find(a.n(), a.kappa(), basis);
      const auto m = DiracSpinor::find(b.n(), b.kappa(), basis);
      if (!n || !m)
        continue;
      const auto dab = d.reducedME(a, b);
      const auto dnm = d.reducedME(*n, *m);
      const auto eps = dnm / dab - 1.0;
      fmt::print("<{:3s}||d||{:3s}> {:+.4e} {:+.4e} {:+.0e}\n", a.shortSymbol(),
                 b.shortSymbol(), dab, dnm, eps);
    }
  }
}

//------------------------------------------------------------------------------
void basisHFS(const std::vector<DiracSpinor> &orbs,
              const std::vector<DiracSpinor> &basis,
              const DiracOperator::hfs &hfs) {
  for (const auto &a : orbs) {
    const auto n = DiracSpinor::find(a.n(), a.kappa(), basis);
    if (!n)
      continue;
    const auto A_const = DiracOperator::Hyperfine::convert_RME_to_HFSconstant(
      1, a.kappa(), a.kappa());
    const auto Aa = hfs.reducedME(a, a) * A_const;
    const auto An = hfs.reducedME(*n, *n) * A_const;
    const auto eps = Aa / An - 1.0;
    fmt::print("{:4s} {:+.4e} {:+.4e} {:+.0e}\n", a.shortSymbol(), Aa, An, eps);
  }
}

//------------------------------------------------------------------------------
void sumRules(const std::vector<DiracSpinor> &basis, const Grid &grid,
              double alpha) {
  std::cout << "\nTKR sum rule:\n";
  SplineBasis::sumrule_TKR(basis, grid.r(), true);
  std::cout << "\nDrake-Gordon sum rule:\n";
  for (int nDG = 0; nDG < 3; ++nDG)
    SplineBasis::sumrule_DG(nDG, basis, grid, alpha, true);
}

//------------------------------------------------------------------------------
void compareBasisSet(const std::string &label, const Wavefunction &wf,
                     const std::vector<DiracSpinor> &basis, bool do_E1,
                     double scale_rN, bool sum_rules,
                     const std::string &sigma_note,
                     bool sigma_note_before_core) {
  const auto has_sigma = (wf.Sigma() != nullptr);
  const auto rN_fm = scale_rN * std::sqrt(5.0 / 3) * wf.get_rrms();
  const auto rN_au = rN_fm / PhysConst::aB_fm;
  const DiracOperator::E1 d(wf.grid());
  DiracOperator::hfs hfs(1, 1.0, rN_au, wf.grid(),
                         DiracOperator::Hyperfine::sphericalBall_F(1));

  auto run = [&](auto fn) {
    if (sigma_note_before_core && has_sigma)
      std::cout << sigma_note << "\n";
    fn(wf.core(), "core");
    if (!sigma_note_before_core && has_sigma)
      std::cout << sigma_note << "\n";
    fn(wf.valence(), "valence");
  };

  // "basis" for Basis, "spect" for Spectrum -- matches original output
  const auto lbl = (label == "Basis") ? "basis" : "spect";
  const auto div =
    "\n----------------------------------------------------------\n";

  std::cout << div << label << ":\n";

  std::cout << div << "Orthonormality: <a_HF|a_" << lbl
            << ">, and worst <a_HF|b_" << lbl << ">:\n";
  run([&](const auto &orbs, const std::string &which) {
    std::cout << "\nOrthonormality: " << label << " vs " << which
              << " HF orbitals:\n";
    basisOrtho(orbs, basis);
  });

  std::cout << div << "Energies and <r> expectation values:\n";
  run([&](const auto &orbs, const std::string &which) {
    std::cout << "\n" << label << " vs " << which << " HF orbitals:\n";
    basisEnergies(orbs, basis, wf.grid().r());
  });

  if (do_E1) {
    std::cout << div << "E1 reduced ME (large r test):\n";
    run([&](const auto &orbs, const std::string &which) {
      std::cout << "\nE1: " << label << " vs " << which << " HF orbitals:\n";
      basisE1(orbs, basis, d);
    });
  }

  std::cout << div << "Hyperfine A (ball model, g=1) (small r test):\n";
  std::cout << "rN = " << rN_fm << " fm\n";
  run([&](const auto &orbs, const std::string &which) {
    std::cout << "\nHFS: " << label << " vs " << which << " HF orbitals:\n";
    basisHFS(orbs, basis, hfs);
  });

  if (sum_rules) {
    std::cout << div << "Sum rules:\n";
    std::cout << "(Only work if -ve energy states are included!)\n";
    sumRules(basis, wf.grid(), wf.alpha());
  }
}

} // namespace Tests

} // namespace Module
