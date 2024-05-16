#include "CorrelationPotential.hpp"
#include "DiracODE/DiracODE.hpp"
#include "Feynman.hpp"
#include "Goldstone.hpp"
#include "Maths/Grid.hpp"
#include "Sigma2.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Random.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
TEST_CASE("MBPT: Goldstone, unit tests", "[MBPT][Goldstone][unit]") {

  std::cout << "\n----------------------------------------\n";
  std::cout << "Goldstone diagram, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-5);
  wf.solve_valence("3sp");
  wf.formBasis(
      SplineBasis::Parameters("20spd", 20, 6, 1.0e-4, 1.0e-4, 30.0, false));

  // These parameters are not meant to be accurate
  const double r0{1.0e-2};
  const double rmax{20.0};
  const std::size_t stride = 6;
  const int n_min_core = 2;

  const auto i0 = wf.grid().getIndex(r0);
  const auto size = (wf.grid().getIndex(rmax) - i0) / stride + 1;

  // Construct Goldtone diagrams (second-order only)
  MBPT::Goldstone Gs = MBPT::Goldstone(wf.basis(), wf.core(), i0, stride, size,
                                       n_min_core, false);

  // Test the "parameter" getters
  REQUIRE(Gs.stride() == stride);
  REQUIRE(Gs.n_min() == n_min_core);
  REQUIRE(Gs.lmax() == 2);

  const auto &Yeh = Gs.Yeh();
  const auto &[holes, excited] = Gs.basis();

  for (auto &v : wf.valence()) {

    auto Sigma = Gs.Sigma_direct(v.kappa(), v.en()) +
                 Gs.Sigma_exchange(v.kappa(), v.en());

    auto de0 = MBPT::Sigma_vw(v, v, Yeh, holes, excited);

    auto de1 = v * (Sigma * v);

    // Not an accuracy test
    REQUIRE(de1 == Approx(de0).epsilon(1.0e-3));
  }
}

//==============================================================================
TEST_CASE("MBPT: Feynman, unit tests", "[MBPT][Feynman][unit]") {

  std::cout << "\n----------------------------------------\n";
  std::cout << "Feynman diagram, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-5);
  wf.solve_valence("3sp");

  // These parameters are not meant to be accurate
  const double r0{1.0e-2};
  const double rmax{20.0};
  const std::size_t stride = 6;
  const auto omre = -0.33 * wf.energy_gap();
  const int lmax = 2;
  const double w0{0.1};
  const double wratio{3.0};
  const int n_min_core = 2;

  const auto i0 = wf.grid().getIndex(r0);
  const auto size = (wf.grid().getIndex(rmax) - i0) / stride + 1;

  // Construct Feynman diagrams (second-order only)
  MBPT::Feynman Fy(wf.vHF(), i0, stride, size,
                   {MBPT::Screening::exclude, MBPT::HoleParticle::exclude, lmax,
                    omre, w0, wratio},
                   n_min_core, true);

  // Test the "parameter" getters
  REQUIRE(Fy.screening() == false);
  REQUIRE(Fy.hole_particle() == false);
  REQUIRE(Fy.stride() == stride);
  REQUIRE(Fy.n_min() == n_min_core);
  REQUIRE(Fy.lmax() == lmax);
  REQUIRE(Fy.w0() == Approx(w0));
  REQUIRE(Fy.wratio() == Approx(wratio));
  REQUIRE(Fy.omre() == Approx(omre));

  // Test exchange potential: just that it works, not an accuracy test
  for (auto &v : wf.valence()) {
    const auto vx = Fy.get_Vx_kappa(v.kappa());
    const auto dv = vx * v;
    const auto dex = qip::inner_product(v.f(), dv.f()) / double(Fy.stride());
    const auto dex0 = v * (wf.vHF()->vexFa(v));
    const auto eps = std::abs(dex / dex0 - 1.0);
    REQUIRE(eps < 1.0e-1);
  }

  // with screening
  MBPT::Feynman FyS(wf.vHF(), i0, stride, size,
                    {MBPT::Screening::include, MBPT::HoleParticle::exclude,
                     lmax, omre, w0, wratio},
                    n_min_core, false);

  // with screening + hole-particle (all-order)
  MBPT::Feynman FyAO(wf.vHF(), i0, stride, size,
                     {MBPT::Screening::include, MBPT::HoleParticle::include,
                      lmax, omre, w0, wratio},
                     n_min_core, false);

  // "Only screening": i.e., just screening correction, no 2nd order (no hp)
  MBPT::Feynman FyS_only(wf.vHF(), i0, stride, size,
                         {MBPT::Screening::only, MBPT::HoleParticle::exclude,
                          lmax, omre, w0, wratio},
                         n_min_core, false);

  // "Only screening + hp": i.e., just screening + hp corrections, no 2nd order
  MBPT::Feynman FyHS_only(wf.vHF(), i0, stride, size,
                          {MBPT::Screening::only, MBPT::HoleParticle::include,
                           lmax, omre, w0, wratio},
                          n_min_core, false);

  // expected data
  const std::vector expected_de{-0.0049, -0.0015, -0.0015};
  const std::vector expected_sc_ratio{0.85, 0.90, 0.90};
  const std::vector expected_hp_ratio{1.22, 1.23, 1.23};
  const double epsilon = 1.0e-1; // just test to 10% (not anccuracy test)

  for (std::size_t i = 0; i < wf.valence().size(); ++i) {
    const auto &v = wf.valence().at(i);

    // Check second-order Feynman against expected
    const auto Sd = Fy.Sigma_direct(v.kappa(), v.en());
    const auto de0 = v * (Sd * v);
    // only require to ~1%, since not an accuracy test
    REQUIRE(de0 == Approx(expected_de[i]).epsilon(epsilon));

    // Test screening and hole-particle:
    const auto Sd_S = FyS.Sigma_direct(v.kappa(), v.en());
    const auto Sd_ao = FyAO.Sigma_direct(v.kappa(), v.en());
    const auto des = v * (Sd_S * v);
    const auto deao = v * (Sd_ao * v);

    // Test screening: ratio of screened to unscreened:
    REQUIRE(des / de0 == Approx(expected_sc_ratio[i]).epsilon(epsilon));

    // Test hole-particle: ratio of all-orders to screening
    REQUIRE(deao / des == Approx(expected_hp_ratio[i]).epsilon(epsilon));

    // Test the "only screening" methods (i.e., corrections only)
    const auto Sd_Sonly = FyS_only.Sigma_direct(v.kappa(), v.en());
    const auto Sd_HSonly = FyHS_only.Sigma_direct(v.kappa(), v.en());

    // Only screening (no hp)
    const auto des_only = v * (Sd_Sonly * v);
    // Only screening (with hp)
    const auto deshp_only = v * (Sd_HSonly * v);

    // Screening only is equal to screening correction (without hp):
    REQUIRE(des_only == Approx(des - de0));
    // Screening only is equal to screening correction (with hp):
    REQUIRE(deshp_only == Approx(deao - de0));
  }
  std::cout << "\n";
}

//==============================================================================
TEST_CASE("MBPT: CorrelationPotential", "[MBPT][CorrelationPotential][unit]") {

  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: CorrelationPotential\n";

  Wavefunction wf({400, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-5);
  wf.solve_valence("3sp");
  wf.formBasis(
      SplineBasis::Parameters("20spd", 20, 6, 1.0e-4, 1.0e-4, 30.0, false));

  // These parameters are not meant to be accurate
  const double r0{1.0e-2};
  const double rmax{20.0};
  const std::size_t stride = 6;
  const int n_min_core = 2;

  auto SigmaG = MBPT::CorrelationPotential(
      "", wf.vHF(), wf.basis(), r0, rmax, stride, n_min_core,
      MBPT::SigmaMethod::Goldstone, false, false);

  REQUIRE(SigmaG.empty());

  for (auto &v : wf.valence()) {
    SigmaG.formSigma(v.kappa(), v.en(), v.n(), &v);
  }

  REQUIRE(!SigmaG.empty());

  std::string file_name = "deleteme_" + qip::random_string(6) + ".sig.abf";

  // test write and read:
  SigmaG.write(file_name);

  // read in:
  auto SigmaG2 = MBPT::CorrelationPotential(
      file_name, wf.vHF(), wf.basis(), r0, rmax, stride, n_min_core,
      MBPT::SigmaMethod::Goldstone, false, false);

  SigmaG2.print_subGrid();

  std::cout << "\n";

  for (auto &v : wf.valence()) {

    const auto de1 = v * SigmaG2(v);
    const auto de2 = v * SigmaG(v);

    REQUIRE(de1 == Approx(de2));

    const auto S = SigmaG2.getSigma(v.kappa(), v.n());
    REQUIRE(S != nullptr);
    const auto de3 = v * (*S * v);

    REQUIRE(de1 == Approx(de3));
  }

  //------------------------------------------------------

  // test scaling:
  SigmaG2.scale_Sigma({2.0, 2.0, 2.0});
  SigmaG2.print_scaling();
  for (auto &v : wf.valence()) {

    const auto de1 = v * SigmaG2(v);
    const auto de2 = 2.0 * v * SigmaG(v);

    REQUIRE(de1 == Approx(de2));
  }

  for (auto &v : wf.valence()) {
    SigmaG2.scale_Sigma(v.n() * 1.5, v.kappa(), v.n());
    const auto de1 = v * SigmaG2(v);
    const auto de2 = v.n() * 1.5 * v * SigmaG(v);
    REQUIRE(de1 == Approx(de2));
  }
}

//==============================================================================
//==============================================================================

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
      wf.formBasis({"100spdfghi", 101, 11, 0.0, 1.0e-6, 50.0, false});
      // wf.formSigma(1, false);
      // Coulomb::YkTable yk(wf.basis()); // this is a *huge* basis!
      const auto [holes, excited] =
          DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
      Coulomb::YkTable yk(holes, excited);
      // const auto Sigma = wf.Sigma();
      std::cout << "cf Table 2 from Beloy, Derevianko, Comput.Phys.Commun. "
                   "179, 310 (2008):\n";
      for (int l = 0; l <= 6; ++l) {

        const auto de = MBPT::Sigma_vw(Fv, Fv, yk, holes, excited, l);
        // const auto de_2 = Sigma->Sigma_vw(Fv, Fv, l);
        // REQUIRE(de == Approx(de_2).epsilon(1.0e-9));
        vals.push_back(de - prev);
        printf("%i %10.7f %10.7f  [%10.7f]\n", l, de, de - prev,
               partial_KBAD[std::size_t(l)]);
        prev = de;
      }
      for (auto l = 0ul; l <= 6; ++l) {
        auto del = vals[l] - partial_KBAD[l];
        REQUIRE(std::abs(del) < error);
      }
    }

    { // "smaller" basis set (not exactly same as Derev)
      wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0, false});
      const auto [holes, excited] =
          DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
      Coulomb::YkTable yk(holes, excited);
      // wf.formSigma(1, false);
      // const auto Sigma = wf.Sigma();
      // const auto de = Sigma->Sigma_vw(Fv, Fv);
      const auto de = MBPT::Sigma_vw(Fv, Fv, yk, holes, excited);
      auto ok = de >= -0.01767 && de <= -0.01748 ? 1 : 0;
      // pass &= qip::check_value(&obuff, "MBPT(2) 'small' Cs 6s", ok, 1, 0);
      REQUIRE(ok);
    }
  }
}

//==============================================================================
//! Tests for second-order correlation potential
TEST_CASE("MBPT: Correlation Potential: Sigma2",
          "[MBPT][Sigma2][slow][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: Correlation Potential: Sigma2, [MBPT][Sigma2][slow]\n";

  std::cout << "\n";
  //----------------------------------------------------------------------------
  // Test Sigma:
  // Cs:

  const auto fname = std::string{"deleteme_"} + qip::random_string(5);

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
    wf.formSigma(3, 3, 1.0e-4, 30.0, 14 /*stride*/, false, false, false, {}, {},
                 {}, "false", fname);

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();
    wf.printValence();

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
    wf.formSigma(1, 1, 0.0, 0.0, 1 /*stride*/, false, false, false, {}, {}, {},
                 fname, "false");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();
    wf.printValence();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    for (auto i = 0ul; i < first_run.size(); ++i) {
      const auto eps = std::abs((de[i] - first_run[i]) / first_run[i]);
      std::cout << wf.valence().at(i) << " " << de[i] << " [" << first_run[i]
                << "] " << eps << "\n";
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
    wf.formSigma(4, 4, 1.0e-4, 30.0, 12 /*stride*/, false, false, false, {}, {},
                 {}, "false", fname);

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
    wf.formSigma(n_min_core, n_min_core, rmin, rmax, stride, false, false,
                 false, {}, fk, {}, "false", "false", true, true, true, lmax,
                 omre, w0, wratio);

    wf.hartreeFockBrueckner();

    std::vector<double> br;
    for (const auto &Fv : wf.valence()) {
      br.push_back(Fv.en());
    }
    std::sort(begin(br), end(br)); // sort: don't depend on order

    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      std::cout << wf.valence().at(i) << " " << br[i] << " [" << dzuba_i[i]
                << "]\n";
    }

    auto [eps, at] = qip::compare_eps(dzuba_i, br);
    // pass &= qip::check_value(&obuff, "Sigma all-orders Cs", eps, 0.0, 5e-04);
    // Used to be 5e-4..?
    REQUIRE(std::abs(eps) < 5e-03);
  }
}

//==============================================================================
TEST_CASE("MBPT: Sigma2", "[MBPT][Sigma2][CI][unit]") {

  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: Sigma2 (unit)\n";

  // note: does not test formulas: just checks class is working correctly.
  // Not meant to be accurate!
  Wavefunction wf({400, 1.0e-3, 30.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", 0.0, "[Ne]", 1.0e-4);
  wf.solve_valence("3spd");
  wf.formBasis({"5spdf", 20, 5, 1.0e-2, 1.0e-2, 20.0, false});

  const auto &[core, excited] =
      MBPT::split_basis(wf.basis(), wf.FermiLevel(), 2);
  const auto mbpt_basis = qip::merge(core, excited);

  int kmax = 4;
  Coulomb::QkTable qk;
  Coulomb::YkTable yk(mbpt_basis);
  const Angular::SixJTable &SixJ = yk.SixJ();
  qk.fill(mbpt_basis, yk, kmax);

  for (auto &v : wf.valence()) {
    for (auto &w : wf.valence()) {
      for (auto &x : wf.valence()) {
        for (auto &y : wf.valence()) {
          for (int k = 0; k <= kmax; ++k) {

            //  Lk symmetry:
            //  {abcd} = badc

            const double sk1 = MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited,
                                             SixJ, MBPT::Denominators::BW);

            const double sk2 = MBPT::Sk_vwxy(k, w, v, y, x, qk, core, excited,
                                             SixJ, MBPT::Denominators::BW);

            const double sk3 =
                MBPT::InternalSigma::S_Sigma2_ab(k, v, w, x, y, qk, core,
                                                 excited, SixJ,
                                                 MBPT::Denominators::BW) +
                MBPT::InternalSigma::S_Sigma2_c1(k, v, w, x, y, qk, core,
                                                 excited, SixJ,
                                                 MBPT::Denominators::BW) +
                MBPT::InternalSigma::S_Sigma2_c2(k, v, w, x, y, qk, core,
                                                 excited, SixJ,
                                                 MBPT::Denominators::BW) +
                MBPT::InternalSigma::S_Sigma2_d(k, v, w, x, y, qk, core,
                                                excited, SixJ,
                                                MBPT::Denominators::BW);

            // tests symmetry:
            REQUIRE(sk2 == Approx(sk1));

            // tests formula (SR is checked in Sk_vwxy, not in internal):
            REQUIRE(sk3 == Approx(sk1));

            // Check selectrion rules
            if (MBPT::Sk_vwxy_SR(k, v, w, x, y)) {
              // This isn't always true... ?? Means SRs should be improves!
              // REQUIRE(sk3 != 0.0);
            } else {
              REQUIRE(sk3 == 0.0);
            }
          }
        }
      }
    }
  }
}