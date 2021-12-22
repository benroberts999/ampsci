#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

namespace UnitTest {

//******************************************************************************
//! Unit tests for B-Spline basis (finite basis of relativistic orbitals)
bool BSplineBasis(std::ostream &obuff) {
  bool pass = true;

  { // Check vs. Hartree-Fock (energies/ortho)
    // Create wavefunction object, solve HF for core + valence
    Wavefunction wf({2500, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear"},
                    {"Cs", -1, "Fermi"});
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d4f");

    // Compare the energy of a Dirac spinor to a double:
    const auto comp_en = [](const auto &Fa, double en) {
      return (en - Fa.en()) / Fa.en();
    };

    // Test splines compared to HF states (for various spline sets)
    for (std::size_t k = 7; k <= 9; k += 2) {
      const auto nlst = k == 7 ? std::vector{75, 85} : std::vector{70, 90};
      for (const auto n : nlst) {

        // Form spline basis:
        const std::string states = "30spdf";
        const auto r0 = 1.0e-5;
        const auto r0_eps = 0.0;
        const auto rmax = 75.0;
        const auto positronQ = false;
        const auto basis = SplineBasis::form_basis(
            {states, std::size_t(n), k, r0, r0_eps, rmax, positronQ}, wf);

        // Check orthonormality <a|b>:
        const auto [eps, str] = DiracSpinor::check_ortho(basis, basis);
        const auto [eps1, str1] = DiracSpinor::check_ortho(basis, wf.core);
        const auto [eps2, str2] = DiracSpinor::check_ortho(basis, wf.valence);

        // Find basis states corresponding to core/valence to compare energies
        // Note: Need large cavity and large basis for this
        std::vector<double> core_en;
        for (const auto &Fc : wf.core) {
          const auto &pFb = std::find(cbegin(basis), cend(basis), Fc);
          core_en.push_back(pFb->en());
        }
        std::vector<double> val_en;
        for (const auto &Fv : wf.valence) {
          const auto &pFb = std::find(cbegin(basis), cend(basis), Fv);
          val_en.push_back(pFb->en());
        }
        // Compare core+valence HF energies to the corresponding splines
        const auto [ce, cs] = qip::compare(wf.core, core_en, comp_en);
        const auto [ve, vs] = qip::compare(wf.valence, val_en, comp_en);

        const std::string label = "HFspl [75] " + std::to_string(n) + "/" +
                                  std::to_string(k) + " (30) ";

        pass &= qip::check_value(&obuff, label + "orth",
                                 qip::max_abs(eps, eps1, eps2), 0.0, 4.0e-5);

        if (std::abs(ce) > std::abs(ve)) {
          pass &= qip::check_value(&obuff, label + "E " + cs->shortSymbol(), ce,
                                   0.0, 3.0e-7);
        } else {
          pass &= qip::check_value(&obuff, label + "E " + vs->shortSymbol(), ve,
                                   0.0, 3.0e-7);
        }
      }
    }
  }

  { // Check low-r behavour (splines vs HF for hyperfine)
    Wavefunction wf({6000, 1.0e-7, 150.0, 0.33 * 150.0, "loglinear"},
                    {"Cs", -1, "Fermi"}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");

    const std::string states = "7sp5d4f";
    wf.solve_valence(states);

    const auto r0 = 1.0e-6;
    const auto r0eps = 0.0;
    const auto rmax = 60.0; // need large rmax, to match to large nl HF
    const int nspl = 70;    // need large n due to large rmax
    const int kspl = 7;
    wf.formBasis({states, nspl, kspl, r0, r0eps, rmax, false});

    // Hyperfine operator: Pointlike, g=1
    const auto h = DiracOperator::HyperfineA(
        1.0, 1.0, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());

    // Calculate A with HF and spline states, compare for each l:
    std::vector<std::vector<double>> hfs(4); // s,p,d,f
    for (const auto orbs : {&wf.core, &wf.valence}) {
      for (const auto &Fv : *orbs) {
        const auto pFb = std::find(cbegin(wf.basis), cend(wf.basis), Fv);
        const auto Ahf = h.hfsA(Fv);
        const auto Aspl = h.hfsA(*pFb);
        const auto eps = (Ahf - Aspl) / Ahf;
        hfs.at(std::size_t(Fv.l())).push_back(eps);
      }
    }
    for (std::size_t l = 0; l < hfs.size(); l++) {
      const auto eps =
          *std::max_element(cbegin(hfs[l]), cend(hfs[l]), qip::comp_abs);
      const auto tol = l >= 3 ? 1.0e-4 : 1.0e-5;
      pass &= qip::check_value(&obuff, "low-r (hfs) l=" + std::to_string(l),
                               eps, 0.0, tol);
    }
  }

  // Test with a Coulomb potential Z=2: TKR and DG sum rules:
  {
    Wavefunction wf({5000, 1.0e-6, 50.0, 0.33 * 50.0, "loglinear", -1.0},
                    {"2", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("Hartree", 0.0, "[]");

    std::string states = "spdfghi";
    std::size_t nspl = 40;
    std::size_t kspl = 9;
    const auto basis = SplineBasis::form_basis(
        {states, nspl, kspl, 1.0e-3, 0.0, 50.0, true}, wf);

    const std::string label =
        "Z=2 [50] " + std::to_string(nspl) + "/" + std::to_string(kspl);

    const auto [eps, str] = DiracSpinor::check_ortho(basis, basis);
    pass &= qip::check_value(&obuff, label + " orth ", eps, 0.0, 3.0e-8);
    {
      const auto tkr = SplineBasis::sumrule_TKR(basis, wf.rgrid->r(), false);

      const auto worst =
          std::max_element(cbegin(tkr), cend(tkr), qip::comp_abs);
      const auto best = std::min_element(cbegin(tkr), cend(tkr), qip::comp_abs);

      const auto blabel =
          label + " TKR(b) l=" + std::to_string(int(best - begin(tkr)));
      const auto wlabel =
          label + " TKR(w) l=" + std::to_string(int(worst - begin(tkr)));

      pass &= qip::check_value(&obuff, blabel, *best, 0.0, 4.0e-8);
      pass &= qip::check_value(&obuff, wlabel, *worst, 0.0, 1.0e-3);
    }

    for (int nDG = 0; nDG <= 2; ++nDG) {
      const auto dg =
          SplineBasis::sumrule_DG(nDG, basis, *wf.rgrid, wf.alpha, false);

      const auto worst = std::max_element(cbegin(dg), cend(dg), qip::comp_abs);
      const auto best = std::min_element(cbegin(dg), cend(dg), qip::comp_abs);

      const auto kib = AtomData::kappaFromIndex(int(best - begin(dg)));
      const auto kiw = AtomData::kappaFromIndex(int(worst - begin(dg)));
      const auto blabel =
          label + " DG" + std::to_string(nDG) + "(b) k=" + std::to_string(kib);
      const auto wlabel =
          label + " DG" + std::to_string(nDG) + "(w) k=" + std::to_string(kiw);

      pass &= qip::check_value(&obuff, blabel, *best, 0.0, 1.0e-9);
      pass &= qip::check_value(&obuff, wlabel, *worst, 0.0, 3.0e-5);
    }
  }

  // Test with a local potential (Hartree): DG sum rules:
  {
    Wavefunction wf({5000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("Hartree", 0.0, "[Xe]");

    std::string states = "spdfghi";
    std::size_t nspl = 50;
    std::size_t kspl = 9;
    const auto basis = SplineBasis::form_basis(
        {states, nspl, kspl, 1.0e-4, 0.0, 50.0, true}, wf);

    const std::string label =
        "Local [50] " + std::to_string(nspl) + "/" + std::to_string(kspl);

    const auto [eps, str] = DiracSpinor::check_ortho(basis, basis);
    pass &= qip::check_value(&obuff, label + " orth ", eps, 0.0, 1.0e-6);

    for (int nDG = 0; nDG <= 2; nDG += 2) {
      const auto dg =
          SplineBasis::sumrule_DG(nDG, basis, *wf.rgrid, wf.alpha, false);

      const auto worst = std::max_element(cbegin(dg), cend(dg), qip::comp_abs);
      const auto best = std::min_element(cbegin(dg), cend(dg), qip::comp_abs);

      const auto kib = AtomData::kappaFromIndex(int(best - begin(dg)));
      const auto kiw = AtomData::kappaFromIndex(int(worst - begin(dg)));
      const auto blabel =
          label + " DG" + std::to_string(nDG) + "(b) k=" + std::to_string(kib);
      const auto wlabel =
          label + " DG" + std::to_string(nDG) + "(w) k=" + std::to_string(kiw);

      pass &= qip::check_value(&obuff, blabel, *best, 0.0, 1.0e-8);
      pass &= qip::check_value(&obuff, wlabel, *worst, 0.0, 5.0e-5);
    }
  }

  return pass;
}

} // namespace UnitTest
