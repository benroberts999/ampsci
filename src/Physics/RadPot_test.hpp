#pragma once
#include "DiracOperator/Operators.hpp"
#include "Physics/RadiativePotential.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace UnitTest {

//******************************************************************************

namespace helper {
struct TabIVdata {
  // IV: J. S. M. Ginges and J. C. Berengut, J. Phys. B 49, 095001 (2016).
  // VI: J. S. M. Ginges and J. C. Berengut, Phys. Rev. A 93, 052509 (2016).
  std::string atom;
  std::string core;
  std::string val;
  std::vector<double> Ude0;
  std::vector<double> Ude;
  std::vector<double> SEde0;
  std::vector<double> SEde;
};
} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests for Ginges/Flambaum Radiative potential method
bool RadPot(std::ostream &obuff) {
  bool pass = true;

  // Test data: s,p,d energies from:
  // IV: J. S. M. Ginges and J. C. Berengut, J. Phys. B 49, 095001 (2016).
  // VI: J. S. M. Ginges and J. C. Berengut, Phys. Rev. A 93, 052509 (2016).
  // Note: I use slightly different A,B (internal rad pot parameters) than the
  // above for d (and higher) states, so exact match is not assumed.
  const auto data = std::vector<helper::TabIVdata>{
      {"Na",
       "[Ne]",
       "3spd",
       // Uehling, first-order shift:
       {-5.559E-07, -2.006E-10, -4.101E-11, -6.538E-18, -2.018E-18},
       // Uehling, + Hartree-Fock (relaxation):
       {-5.910E-07, 3.067E-08, 3.077E-08, -1.323E-10, -1.347E-10},
       // Self-energy, 0th order:
       {1.0680E-05 /*, -5.6980E-08, 1.1120E-07, -2.6440E-09, 1.7700E-09*/},
       // Self-energy + relaxation:
       {1.1250E-05 /*, -7.0880E-07, -4.8820E-07, 4.0850E-11, 3.4070E-09*/}},
      {"K",
       "[Ar]",
       "4sp3d",
       {-1.224E-06, -1.879E-09, -3.547E-10, -1.046E-14, -3.112E-15},
       {-1.333E-06, 7.262E-08, 7.396E-08, 1.395E-08, 1.376E-08},
       {1.8450E-05 /*, -1.0230E-07, 3.0720E-07, -1.9150E-08, 1.3820E-08*/},
       {1.9740E-05 /*, -1.4000E-06, -8.4110E-07, -2.5680E-07,
       -2.6470E-07*/}},
      {"Rb",
       "[Kr]",
       "5sp4d",
       {-4.648E-06, -3.234E-08, -4.881E-09, -3.865E-12, -1.048E-12},
       {-5.114E-06, 1.776E-07, 2.108E-07, 1.290E-07, 1.218E-07},
       {4.8360E-05 /*, 1.3810E-08, 1.3160E-06, -1.0980E-07, 1.0050E-07*/},
       {5.1360E-05 /*, -2.8540E-06, -1.0830E-06, -1.7910E-06, -1.7900E-06*/}},
      {"Cs",
       "[Xe]",
       "6sp5d",
       {-1.054E-05, -1.942E-07, -2.178E-08, -1.938E-10, -4.612E-11},
       {-1.160E-05, 2.288E-07, 4.518E-07, 1.127E-06, 1.007E-06},
       {8.1280E-05 /*, 1.0770E-06, 3.1830E-06, -6.0660E-07, 7.1740E-07*/},
       {8.4310E-05 /*, -3.8310E-06, -9.2030E-07, -1.2120E-05, -1.1150E-05*/}},
      {"Fr",
       "[Rn]",
       "7sp6d",
       {-4.997E-05, -2.569E-06, -1.356E-07, -3.579E-09, -6.921E-10},
       {-5.540E-05, -1.610E-06, 1.834E-06, 3.822E-06, 2.606E-06},
       {2.2010E-04 /*, 1.0680E-05, 1.0340E-05, -8.0460E-07, 1.9680E-06*/},
       {2.1660E-04 /*, 1.2760E-09, -5.8880E-08, -2.6680E-05, -2.2470E-05*/}}
      ,{"E119",
       "[Og]",
       "8sp7d",
       {-3.427E-04, -3.921E-05, -5.420E-07, -2.517E-08, -4.420E-09},
       {-4.046E-04, -4.701E-05, 1.120E-05, 8.715E-06, -1.374E-06},
       {7.1960E-04 /*, 7.2040E-05, 2.6970E-05, 4.5230E-07, 4.7170E-06*/},
       {6.8320E-04 /*, 5.2170E-05, -2.0010E-06, -4.0690E-05, -3.0860E-05*/}}
  };

  // comparators for custom struct
  auto cmpr = [](const auto &ds, const auto &gb) {
    return (std::get<1>(ds) - gb) / gb;
  };
  auto cmpr_hf = [](const auto &ds, const auto &gb) {
    return (std::get<2>(ds) - gb) / gb;
  };
  const bool rw = false; // don't read/write QED rad pot to file

  // Run test for each atom, and for Uehling (vac pol) and Self-energy
  // seperately
  for (const auto &at : data) {
    for (const std::string pot : {"Ueh", "SE"}) {

      // Ensure we use same parameters as above paper:
      const int aa = at.atom == "Fr" ? 211 : -1; // Fr-211 not 213
      // Must use exact same rms value as above papers
      const double rrms = at.atom == "E119" ? 6.5 : -1.0;

      // Construct wavefunction, solve HF core+valence (without QED):
      Wavefunction wf({5000, 1.0e-6, 150.0, 0.3 * 150.0, "loglinear", -1.0},
                      {at.atom, aa, "Fermi", rrms, 2.3}, 1.0);
      wf.hartreeFockCore("HartreeFock", 0.0, at.core);
      wf.hartreeFockValence(at.val);

      // Solve new wavefunction, WITH QED corrections (into Hartree-Fock):
      Wavefunction wf2({5000, 1.0e-6, 150.0, 0.3 * 150.0, "loglinear", -1.0},
                       {at.atom, aa, "Fermi", rrms, 2.3}, 1.0);
      const auto rcut = 15.0;
      if (pot == "Ueh")
        wf2.radiativePotential(0.0, 1.0, 0.0, 0.0, 0.0, rcut, 1.0, {1.0}, rw);
      else
        wf2.radiativePotential(0.0, 0.0, 1.0, 1.0, 1.0, rcut, 1.0, {1.0}, rw);
      wf2.hartreeFockCore("HartreeFock", 0.0, at.core);
      wf2.hartreeFockValence(at.val);

      // First, test the zeroth-order energy shifts.
      // That is <a|h|a>, where a is wavefunction WITHOUT QED, and h is QED
      // operator
      // Then, test (HF+QED) [includes relaxtion]
      auto h = DiracOperator::Hrad_el(wf2.vrad.get_Hel(0));
      auto hm = DiracOperator::Hrad_mag(wf2.vrad.get_Hmag(0));
      // result: {name, 0th, HF/relaxation}
      std::vector<std::tuple<std::string, double, double>> result;
      for (auto i = 0ul; i < wf2.valence.size(); ++i) {
        const auto &Fv = wf.valence[i];
        const auto &Fv2 = wf2.valence[i];
        if (pot == "SE" && Fv2.l() != 0)
          continue;
        result.emplace_back(
            Fv.symbol(), h.radialIntegral(Fv, Fv) + hm.radialIntegral(Fv, Fv),
            Fv2.en - Fv.en);
      }

      // Compare results to test data:
      if (pot == "SE") {
        const auto [eps, pos] = qip::compare(result, at.SEde0, cmpr);
        pass &= qip::check_value(&obuff, "SE(s) de0 (TabVI) " + at.atom, eps,
                                 0.0, 5.0e-4);
        const auto [eps2, pos2] = qip::compare(result, at.SEde, cmpr_hf);
        pass &= qip::check_value(&obuff, "SE(s) deHF (TabVI) " + at.atom, eps2,
                                 0.0, 2.0e-3);
      } else {
        const auto [eps, pos] = qip::compare(result, at.Ude0, cmpr);
        pass &= qip::check_value(&obuff, "Ueh de0 (TabIV) " + at.atom, eps, 0.0,
                                 5.0e-4);
        const auto [eps2, pos2] = qip::compare(result, at.Ude, cmpr_hf);
        // this is the limit of number of decimal places (for Na)
        pass &= qip::check_value(&obuff, "Ueh deHF (TabIV) " + at.atom, eps2,
                                 0.0, 1.0e-3);
      }
    }
  }

  return pass;
}

} // namespace UnitTest
