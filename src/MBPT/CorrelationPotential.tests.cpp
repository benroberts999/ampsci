#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
//! Unit tests for second-order MBPT energy correction
TEST_CASE("MBPT: 2nd Order de", "[MBPT][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: 2nd Order de\n";

  { // Compare with  K. Beloy and A. Derevianko,
    // Comput. Phys. Commun. 179, 310 (2008).
    Wavefunction wf({4000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("6s");
    const auto &Fv = wf.valence().front();

    {
      // K. Beloy and A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
      const auto partial_KBAD =
          std::vector{-0.0000130, -0.0020027, -0.0105623, -0.0039347,
                      -0.0007563, -0.0002737, -0.0001182};
      // const auto total_KBAD = -0.0176609;
      const auto error = 0.0000002; // allow ony difference of 1.5 in last digit

      double prev = 0.0;
      std::vector<double> vals;
      wf.formBasis({"spdfghi", 100, 11, 0.0, 1.0e-6, 50.0, false});
      wf.formSigma(1, false);
      const auto Sigma = wf.Sigma();
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
        // pass &= qip::check_value(
        //     &obuff, "MBPT(2) vs. KB,AD " + std::to_string(l), del, 0.0,
        //     error);
        REQUIRE(std::abs(del) < error);
      }
    }

    { // "smaller" basis set (not exactly same as Derev)
      wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
      wf.formSigma(1, false);
      const auto Sigma = wf.Sigma();
      const auto de = Sigma->SOEnergyShift(Fv, Fv);
      auto ok = de >= -0.01767 && de <= -0.01748 ? 1 : 0;
      // pass &= qip::check_value(&obuff, "MBPT(2) 'small' Cs 6s", ok, 1, 0);
      REQUIRE(ok);
    }
  }
}

//! Unit tests for second-order correlation potential
TEST_CASE("MBPT: Correlation Potential: Sigma2",
          "[MBPT][Sigma2][slow][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: Correlation Potential: Sigma2, [MBPT][Sigma2][slow]\n";

  std::cout << "\n";
  //----------------------------------------------------------------------------
  // Test Sigma:
  // Cs:
  std::vector<double> first_run;
  { // Compare Dzuba, using up to l=6 for splines
    std::cout << "Test Sigma(2) Brueckner, for Cs:\n";
    auto dzuba_i = std::vector{
        -0.02013813, -0.00410942, -0.00792483, -0.00702407, -0.00220878,
        -0.00199737, -0.01551449, -0.01466935, -0.00035253, -0.00035234};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d4f");
    wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
    wf.formSigma(3, true, 1.0e-4, 30.0, 14 /*stride*/, false, false, {}, {}, {},
                 "false", "tmp_sigma_deleteme");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    std::cout << "delta Sigma(2) Bruckner, cf Dzuba:\n";
    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      const auto eps = std::abs((de[i] - dzuba_i[i]) / dzuba_i[i]);
      std::cout << de[i] << " [" << dzuba_i[i] << "] " << eps << "\n";
    }
    first_run = de; // copy this data, test the next run against:

    const auto [eps, at] = qip::compare_eps(dzuba_i, de);
    // pass &= qip::check_value(&obuff, "Sigma2 Cs", eps, 0.0, 0.01);
    REQUIRE(std::abs(eps) < 0.01);
  }

  std::cout << "\n";

  //----------------------------------------------------------------------------
  // Test reading in Sigma:
  // Cs:
  {
    std::cout << "Test reading in Sigma(2) Brueckner file, for Cs:\n";
    std::sort(begin(first_run), end(first_run)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d4f");
    // wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
    // Don't calculate Sigma, read it in from above example:
    wf.formSigma(1, true, 0.0, 0.0, 1 /*stride*/, false, false, {}, {}, {},
                 "tmp_sigma_deleteme", "false");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    for (auto i = 0ul; i < first_run.size(); ++i) {
      const auto eps = std::abs((de[i] - first_run[i]) / first_run[i]);
      std::cout << de[i] << " [" << first_run[i] << "] " << eps << "\n";
    }

    const auto [eps, at] = qip::compare_eps(first_run, de);
    // pass &= qip::check_value(&obuff, "Sigma2 Cs (read)", eps, 0.0, 1.0e-16);
    REQUIRE(std::abs(eps) < 1.0e-16);
  }

  std::cout << "\n";

  //----------------------------------------------------------------------------
  // Fr:
  { // Compare Dzuba, using up to l=6 for splines
    std::cout << "Test Sigma(2) Brueckner, for Fr:\n";
    auto dzuba_i =
        std::vector{-0.0245075, -0.0098094, -0.0069442, -0.0153430, -0.0133382};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Fr", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Rn]");
    wf.solve_valence("7sp6d");
    wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
    wf.formSigma(4, true, 1.0e-4, 30.0, 12 /*stride*/, false, false, {}, {}, {},
                 "false", "tmp_sigma_deleteme");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      const auto eps = std::abs((de[i] - dzuba_i[i]) / dzuba_i[i]);
      std::cout << de[i] << " [" << dzuba_i[i] << "] " << eps << "\n";
    }

    const auto [eps, at] = qip::compare_eps(dzuba_i, de);
    // pass &= qip::check_value(&obuff, "Sigma2 Fr", eps, 0.0, 0.02);
    REQUIRE(std::abs(eps) < 0.02);
  }
}

//==============================================================================
//! Unit tests for all-orders correlation potential
TEST_CASE("MBPT: Correlation Potential: SigmaAO",
          "[MBPT][SigmaAO][slow][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: Correlation Potential: SigmaAO, [MBPT][SigmaAO][slow]\n";

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
        int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;

    const auto omre = -std::abs(0.33 * wf.energy_gap());
    const double w0 = 0.01;
    const double wratio = 1.5;
    const auto lmax = 6;

    const std::vector fk{0.71, 0.589, 0.84, 0.885, 0.95, 0.976, 0.991};
    // const std::vector fk{0.72, 0.62, 0.83, 0.89, 0.94, 1.0};
    // wf.formSigma(3, true, 1.0e-4, 30.0, 14 /*stride*/);
    wf.formSigma(n_min_core, true, rmin, rmax, stride, false, false, {}, fk, {},
                 "false", "false", "", true, true, true, lmax, false, false,
                 omre, w0, wratio);

    wf.hartreeFockBrueckner();

    std::vector<double> br;
    for (const auto &Fv : wf.valence()) {
      br.push_back(Fv.en());
    }
    std::sort(begin(br), end(br)); // sort: don't depend on order

    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      std::cout << dzuba_i[i] << " " << br[i] << "\n";
    }

    auto [eps, at] = qip::compare_eps(dzuba_i, br);
    // pass &= qip::check_value(&obuff, "Sigma all-orders Cs", eps, 0.0, 5e-04);
    REQUIRE(std::abs(eps) < 5e-04);
  }
}
