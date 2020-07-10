#pragma once
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
namespace helper {

// Move to DiracSpinor!
std::pair<double, std::string> ortho(const std::vector<DiracSpinor> &a,
                                     const std::vector<DiracSpinor> &b) {
  double worst_del = 0.0;
  std::string worst_F = "";
  for (const auto &Fa : a) {
    for (const auto &Fb : b) {
      if (Fb.k != Fa.k)
        continue;
      const auto del =
          Fa == Fb ? std::abs(std::abs(Fa * Fb) - 1.0) : std::abs(Fa * Fb);
      // nb: sometimes sign of Fb is wrong. Perhaps this is an issue??
      if (del > worst_del) {
        worst_del = del;
        worst_F = "<" + Fa.shortSymbol() + "|" + Fb.shortSymbol() + ">";
      }
    }
  }
  return {worst_del, worst_F};
}

} // namespace helper

//******************************************************************************
//******************************************************************************
bool BSplineBasis(std::ostream &obuff) {
  bool pass = true;
  {
    // Create wavefunction object
    Wavefunction wf({2500, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
    wf.hartreeFockValence("7sp5d4f");

    auto comp_en = [](const auto &Fa, const auto &en) {
      return (en - Fa.en) / Fa.en;
    };

    // test cf HF states
    for (std::size_t k = 7; k <= 9; k += 2) {
      std::size_t n0 = k == 7 ? 75 : 70;
      for (auto n = n0; n <= 90; n += 10) {
        std::string states = "30spdf";
        const auto basis = SplineBasis::form_basis(
            {states, n, k, 0.0, 1.0e-9, 75.0, false}, wf);

        // Orthonormality:
        const auto [eps, str] = helper::ortho(basis, basis);
        const auto [eps1, str1] = helper::ortho(basis, wf.core);
        const auto [eps2, str2] = helper::ortho(basis, wf.valence);

        // Find basis states corresponding to core/valence to compare energies
        // Note: Need large cavity and large basis for this
        std::vector<double> core_en;
        for (const auto &Fc : wf.core) {
          const auto &pFb = std::find(cbegin(basis), cend(basis), Fc);
          core_en.push_back(pFb->en);
        }
        std::vector<double> val_en;
        for (const auto &Fv : wf.valence) {
          const auto &pFb = std::find(cbegin(basis), cend(basis), Fv);
          val_en.push_back(pFb->en);
        }
        const auto [ce, cs] = qip::compare(wf.core, core_en, comp_en);
        const auto [ve, vs] = qip::compare(wf.valence, val_en, comp_en);

        std::string label = "HFspl [75] " + std::to_string(n) + "/" +
                            std::to_string(k) + " (30) ";

        pass &= qip::check_value(&obuff, label + "orth",
                                 std::max({eps, eps1, eps2}), 0.0, 4.0e-5);
        if (std::abs(ce) > std::abs(ve)) {
          pass &= qip::check_value(&obuff, label + "E_c " + cs->shortSymbol(),
                                   ce, 0.0, 3.0e-7);
        } else {
          pass &= qip::check_value(&obuff, label + "E_v " + vs->shortSymbol(),
                                   ve, 0.0, 3.0e-7);
        }
      }
    }
  }

  // Test with a Coulomb potential Z=2
  {
    Wavefunction wf({5000, 1.0e-6, 50.0, 0.33 * 50.0, "loglinear", -1.0},
                    {"2", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.hartreeFockCore("Hartree", 0.0, "[]");

    std::string states = "spdfghi";
    std::size_t nspl = 40;
    std::size_t kspl = 9;
    const auto basis = SplineBasis::form_basis(
        {states, nspl, kspl, 1.0e-3, 0.0, 50.0, true}, wf);

    const std::string label =
        "Z=2 [50] " + std::to_string(nspl) + "/" + std::to_string(kspl);

    const auto [eps, str] = helper::ortho(basis, basis);
    pass &= qip::check_value(&obuff, label + " orth ", eps, 0.0, 1.0e-8);
    {
      const auto tkr = SplineBasis::sumrule_TKR(basis, wf.rgrid->r, false);

      const auto worst =
          std::max_element(cbegin(tkr), cend(tkr), qip::comp_abs);
      const auto best = std::min_element(cbegin(tkr), cend(tkr), qip::comp_abs);

      const auto blabel =
          label + " TKR(b) l=" + std::to_string(int(best - begin(tkr)));
      const auto wlabel =
          label + " TKR(w) l=" + std::to_string(int(worst - begin(tkr)));

      pass &= qip::check_value(&obuff, blabel, *best, 0.0, 4.0e-8);
      pass &= qip::check_value(&obuff, wlabel, *worst, 0.0, 1.0e-4);
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

      pass &= qip::check_value(&obuff, blabel, *best, 0.0, 1.0e-10);
      pass &= qip::check_value(&obuff, wlabel, *worst, 0.0, 3.0e-5);
    }
  }

  // Test with a local potential (Hartree)
  {
    Wavefunction wf({5000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.hartreeFockCore("Hartree", 0.0, "[Xe]");

    std::string states = "spdfghi";
    std::size_t nspl = 50;
    std::size_t kspl = 9;
    const auto basis = SplineBasis::form_basis(
        {states, nspl, kspl, 1.0e-4, 0.0, 50.0, true}, wf);

    const std::string label =
        "Local [50] " + std::to_string(nspl) + "/" + std::to_string(kspl);

    const auto [eps, str] = helper::ortho(basis, basis);
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
