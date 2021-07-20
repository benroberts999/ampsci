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
//! Unit tests for correlations (second-order MBPT correlation energy/potential)
bool CorrelationPotential(std::ostream &obuff) {
  bool pass = true;

  { // Compare with  K. Beloy and A. Derevianko,
    // Comput. Phys. Commun. 179, 310 (2008).
    Wavefunction wf({4000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("6s");
    const auto &Fv = wf.valence.front();

    {
      // K. Beloy and A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
      const auto partial_KBAD =
          std::vector{-0.0000130, -0.0020027, -0.0105623, -0.0039347,
                      -0.0007563, -0.0002737, -0.0001182};
      // const auto total_KBAD = -0.0176609;
      const auto error = 0.0000001; // allow ony difference of 1 in last digit

      double prev = 0.0;
      std::vector<double> vals;
      wf.formBasis({"spdfghi", 100, 11, 0.0, 1.0e-6, 50.0, false});
      wf.formSigma(1, false);
      const auto Sigma = wf.getSigma();
      std::cout << "cf Table 2 from Beloy, Derevianko, Comput.Phys.Commun. "
                   "179, 310 (2008):\n";
      for (int l = 0; l <= 6; ++l) {
        const auto de = Sigma->SOEnergyShift(Fv, Fv, l);
        vals.push_back(de - prev);
        printf("%i %10.7f %10.7f  [%10.7f]\n", l, de, de - prev,
               partial_KBAD[std::size_t(l)]);
        prev = de;
      }
      for (auto l = 0ul; l <= 6; ++l) {
        auto del = vals[l] - partial_KBAD[l];
        pass &= qip::check_value(
            &obuff, "MBPT(2) vs. KB,AD " + std::to_string(l), del, 0.0, error);
      }
    }

    { // "smaller" basis set (not exactly same as Derev)
      wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
      wf.formSigma(1, false);
      const auto Sigma = wf.getSigma();
      const auto de = Sigma->SOEnergyShift(Fv, Fv);
      auto ok = de >= -0.01767 && de <= -0.01748 ? 1 : 0;
      pass &= qip::check_value(&obuff, "MBPT(2) 'small' Cs 6s", ok, 1, 0);
    }
  }

  //****************************************************************************
  { // Compare Dzuba, only using up to l=4 for splines
    // Note: Works pretty well up to f states (not sure if difference is ok)
    auto dzuba_g = std::vector{
        -0.0195938,  -0.00399679, -0.00770113, -0.00682331, -0.00214125,
        -0.00193494, -0.01400596, -0.01324942, -0.00033882, -0.00033866};
    std::sort(begin(dzuba_g), end(dzuba_g)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 120.0, 0.33 * 120.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d4f"); //"7sp5d4f"
    wf.formBasis({"30spdfg", 40, 7, 0.0, 1.0e-6, 40.0, false});
    wf.formSigma(3, true, 1.0e-4, 30.0, 14 /*stride*/);

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();

    for (const auto &Fv : wf.valence) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order

    auto [eps, at] = qip::compare_eps(dzuba_g, de);
    pass &= qip::check_value(&obuff, "Sigma2 Cs (spdfg)", eps, 0.0, 0.05);
  }
  //

  { // Compare Dzuba, using up to l=6 for splines
    auto dzuba_i = std::vector{
        -0.02013813, -0.00410942, -0.00792483, -0.00702407, -0.00220878,
        -0.00199737, -0.01551449, -0.01466935, -0.00035253, -0.00035234};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d4f");
    wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
    wf.formSigma(3, true, 1.0e-4, 30.0, 14 /*stride*/);

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();

    for (const auto &Fv : wf.valence) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      std::cout << dzuba_i[i] << " " << de[i] << "\n";
    }

    auto [eps, at] = qip::compare_eps(dzuba_i, de);
    pass &= qip::check_value(&obuff, "Sigma2 Cs (spdfghi)", eps, 0.0, 0.01);
  }

  { // Compare Dzuba, All-order sigma
    auto dzuba_i =
        std::vector{-0.14332871, -0.05844404, -0.09244689, -0.08985968,
                    -0.04392404, -0.04309476, -0.07812666, -0.07759564};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({4000, 1.0e-6, 120.0, 0.33 * 120.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d");
    wf.formBasis({"35spdfghi", 40, 9, 0.0, 1.0e-6, 40.0, false});

    const auto n_min_core = 3;
    const auto rmin = 1.0e-4;
    const auto rmax = 30.0;
    const auto stride =
        int(wf.rgrid->getIndex(30.0) - wf.rgrid->getIndex(1.0e-4)) / 150;

    const auto omre = -std::abs(0.33 * wf.energy_gap());
    const double w0 = 0.01;
    const double wratio = 1.5;
    const auto lmax = 6;

    const std::vector fk{0.71, 0.589, 0.84, 0.885, 0.95, 0.976, 0.991};
    // const std::vector fk{0.72, 0.62, 0.83, 0.89, 0.94, 1.0};

    // wf.formSigma(3, true, 1.0e-4, 30.0, 14 /*stride*/);
    wf.formSigma(n_min_core, true, rmin, rmax, stride, false, false, {}, fk,
                 "false", "false", true, true, true, lmax, false, false, omre,
                 w0, wratio);

    wf.hartreeFockBrueckner();

    std::vector<double> br;
    for (const auto &Fv : wf.valence) {
      br.push_back(Fv.en());
    }
    std::sort(begin(br), end(br)); // sort: don't depend on order

    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      std::cout << dzuba_i[i] << " " << br[i] << "\n";
    }

    auto [eps, at] = qip::compare_eps(dzuba_i, br);
    pass &= qip::check_value(&obuff, "Sigma all-orders Cs", eps, 0.0, 5e-04);
  }

  return pass;
}

} // namespace UnitTest
