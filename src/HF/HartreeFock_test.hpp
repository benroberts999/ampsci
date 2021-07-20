#pragma once
#include "DiracOperator/Operators.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <string>
#include <tuple>

namespace UnitTest {

//! Unit tests for Hartree Fock equations
bool HartreeFock(std::ostream &obuff) {
  bool pass = true;

  // lambdas for comparing results to test data:
  const auto cmpr = [](const auto &ds, double v) { return (ds.en() - v) / v; };
  const auto cmpr2 = [](const auto &ds, const auto &ds2) {
    return (ds.en() - ds2.en()) / ds2.en();
  };

  { // Test case: Cs (typical)

    // Grid parameters:
    const auto r0 = 1.0e-6;
    const auto rmax = 120.0;
    const auto b = 0.3 * rmax;
    const auto x_Breit = 0.0; // do not include Breit

    // Create wavefunction object + solve HF
    Wavefunction wf({2000, r0, rmax, b, "loglinear"}, {"Cs", 133, "Fermi"});
    wf.solve_core("HartreeFock", x_Breit, "[Xe]");
    wf.solve_valence("6sp5d");

    // do again, this time with more points in the grid:
    Wavefunction wf2({7000, r0, rmax, b, "loglinear"}, {"Cs", 133, "Fermi"});
    wf2.solve_core("HartreeFock", x_Breit, "[Xe]");
    wf2.solve_valence("6sp5d");

    // Test data: core energies (from Dzuba code)
    const std::vector expected_VD{
        -1330.118692, -212.5644584, -199.4294701, -186.4365217, -45.96974213,
        -40.44831295, -37.89428811, -28.30949791, -27.77513323, -9.512822,
        -7.44628727,  -6.92099641,  -3.48561693,  -3.39689621,  -1.48980554,
        -0.90789792,  -0.84033918};
    // Test data: valence energies (from Dzuba code)
    const std::vector expected_VD_v{-0.12736810, -0.08561590, -0.08378548,
                                    -0.06441963, -0.06452976};

    { // Check HF core vs test data:
      const auto [eps, at] = qip::compare(wf.core, expected_VD, cmpr);
      const std::string worst = at == wf.core.end() ? "" : at->symbol();
      pass &= qip::check_value(&obuff, "HF core Cs vs. VD " + worst, eps, 0.0,
                               5.0e-6);
    }

    { // compare 'dense grid' version to 'normal grid' (core)
      const auto [eps, at] = qip::compare(wf.core, wf2.core, cmpr2);
      const std::string worst = at == wf.core.end() ? "" : at->symbol();
      pass &= qip::check_value(&obuff, "HF core Cs grid " + worst, eps, 0.0,
                               5.0e-6);
    }

    { // Check Valence vs test data:
      const auto [eps, at] = qip::compare(wf.valence, expected_VD_v, cmpr);
      const std::string worst = at == wf.valence.end() ? "" : at->symbol();
      pass &= qip::check_value(&obuff, "HF val Cs vs. VD " + worst, eps, 0.0,
                               3.0e-6);
    }

    { // compare 'dense grid' version to 'normal grid' (valence)
      const auto [eps, at] = qip::compare(wf.valence, wf2.valence, cmpr2);
      const std::string worst = at == wf.valence.end() ? "" : at->symbol();
      pass &=
          qip::check_value(&obuff, "HF val Cs grid " + worst, eps, 0.0, 3.0e-6);
    }
  }

  //****************************************************************************

  { // Test case: Yb+ (one that fails more easily)
    Wavefunction wf3({2000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                     {"Yb", 173, "Fermi", -1.0, -1.0}, 1.0);
    wf3.solve_core("HartreeFock", 0.0, "[Xe],4f14");
    wf3.solve_valence("6sp5d5f");

    const std::vector dzuba{-0.41366446, -0.30111301, -0.28830673, -0.30307172,
                            -0.30088656, -0.12509665, -0.12507894};
    const auto [eps, at] = qip::compare(wf3.valence, dzuba, cmpr);
    pass &= qip::check_value(&obuff, "HF Yb+ val", eps, 0.0, 7.0e-6);
  }

  //****************************************************************************

  { // Test case: Hyperfine constants (HF; no RPA) for Rb
    // cf: Grunefeld, Roberts, Ginges, Phys. Rev. A 100, 042506 (2019).
    Wavefunction wf({4000, 1.0e-7, 425.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Rb", 87, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Kr]");
    wf.solve_valence("12sp");

    // Test data [Phys. Rev. A 100, 042506 (2019)]
    const std::vector s{2183.0, 583.1, 238.8, 120.6,
                        69.24,  43.37, 28.95, 20.27};
    const std::vector p{236.8, 83.21, 38.60, 20.97, 12.64, 8.192, 5.610, 4.008};

    // Generate operator: use exact same parameters as paper:
    const auto h = DiracOperator::HyperfineA(
        2.751818, 1.5, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());

    // Calculate HFS A constant for each valence state (store s,p seperately)
    std::vector<double> sme, pme;
    for (const auto &Fv : wf.valence) {
      if (Fv.k == -1)
        sme.push_back(h.hfsA(Fv));
      else if (Fv.k == 1)
        pme.push_back(h.hfsA(Fv));
    }
    const auto [es, ats] = qip::compare_eps(sme, s);
    const auto [ep, atp] = qip::compare_eps(pme, p);
    pass &= qip::check_value(&obuff, "HF hfs Rb", qip::min_abs(es, ep), 0.0,
                             3.0e-4);
  }

  { // Test case: Hyperfine constants (HF; no RPA) for Cs
    // Grunefeld, Roberts, Ginges, Phys. Rev. A 100, 042506 (2019).
    Wavefunction wf({4000, 1.0e-7, 400.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", 133, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("12sp");
    const std::vector s{1434.0, 393.9, 164.5, 84.11, 48.71, 30.70, 20.59};
    const std::vector p{161.0, 57.65, 27.09, 14.85, 9.002, 5.863, 4.030};

    // Use exact same parameters as paper:
    const auto h = DiracOperator::HyperfineA(
        2.582025, 3.5, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());

    std::vector<double> sme, pme;
    for (const auto &Fv : wf.valence) {
      if (Fv.k == -1)
        sme.push_back(h.hfsA(Fv));
      else if (Fv.k == 1)
        pme.push_back(h.hfsA(Fv));
    }
    const auto [es, ats] = qip::compare_eps(sme, s);
    const auto [ep, atp] = qip::compare_eps(pme, p);
    pass &= qip::check_value(&obuff, "HF hfs Cs", qip::min_abs(es, ep), 0.0,
                             3.0e-4);
  }

  { // Test case: Hyperfine constants (HF; no RPA) for Fr
    // Grunefeld, Roberts, Ginges, Phys. Rev. A 100, 042506 (2019).
    Wavefunction wf({4000, 1.0e-7, 325.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Fr", 211, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Rn]");
    wf.solve_valence("12sp");
    const std::vector s{5929.0, 1520.0, 624.0, 316.8, 182.7, 114.9};
    const std::vector p{628.2, 222.9, 104.4, 57.13, 34.61, 22.53};

    // Use exact same parameters as paper:
    const auto h = DiracOperator::HyperfineA(
        4.00, 4.5, 0.0, *wf.rgrid, DiracOperator::Hyperfine::pointlike_F());

    std::vector<double> sme, pme;
    for (const auto &Fv : wf.valence) {
      if (Fv.k == -1)
        sme.push_back(h.hfsA(Fv));
      else if (Fv.k == 1)
        pme.push_back(h.hfsA(Fv));
    }
    const auto [es, ats] = qip::compare_eps(sme, s);
    const auto [ep, atp] = qip::compare_eps(pme, p);
    pass &= qip::check_value(&obuff, "HF hfs Fr", qip::min_abs(es, ep), 0.0,
                             3.0e-4);
  }

  return pass;
}

} // namespace UnitTest
