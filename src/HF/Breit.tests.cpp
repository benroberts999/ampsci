#include "HF/Breit.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "IO/ChronoTimer.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <string>
#include <utility>
#include <vector>

//! Unit tests for Breit
TEST_CASE("Breit (local)", "[Breit][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Breit (local)\n";

  // solve HF (local) without breit. Calc Breit seperately:
  Wavefunction wf({750, 1.0e-4, 120.0, 40.0, "loglinear", -1.0},
                  {"Rb", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", 0.0, "[Kr]");
  wf.solve_valence("5sp");

  HF::Breit Vb(1.0);

  REQUIRE(std::abs(Vb.scale_factor() - 1.0) < 1.0e-6);

  // energies:
  // expected generated with 10,000 pts
  const auto expected =
      std::vector{1.446617548413e-04, 9.357310567545e-05, 6.905846245113e-05};
  std::size_t count = 0;
  for (const auto &Fv : wf.valence()) {
    auto de = Fv * Vb.VbrFa(Fv, wf.core());
    std::cout << Fv << " - " << de << " [" << expected.at(count) << "]\n";
    const auto eps = std::abs((de - expected.at(count)) / expected.at(count));
    // REQUIRE(de2 == Approx(de));
    printf("%3s %.6e [%.6e] %.0e\n", Fv.shortSymbol().c_str(), de,
           expected.at(count), eps);
    REQUIRE(eps < 1.0e-3);
    ++count;
  }

  {
    IO::ChronoTimer t1("old");
    for (int i = 0; i < 100; ++i) {
      double de1 = 0.0;
      for (const auto &Fv : wf.core()) {
        de1 += Fv * Vb.VbrFa(Fv, wf.core());
      }
      for (const auto &Fv : wf.valence()) {
        de1 += Fv * Vb.VbrFa(Fv, wf.core());
      }
      Vb.update_scale();
    }
  }
}

//==============================================================================
TEST_CASE("Breit (HF)", "[Breit][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Breit (HF) [lowest-order only]\n";

  // solve HF (local) without breit. Calc Breit seperately:
  Wavefunction wf({2000, 1.0e-6, 150.0, 40.0, "loglinear", -1.0},
                  {"Cu", 63, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ar],3d10");
  wf.solve_valence("4sp");

  const HF::Breit Vb(1.0);

  const auto expected_B1 = std::vector{
      std::pair{"4s", 1.880e-4}, {"4p-", 7.015e-5}, {"4p+", 5.140e-5}};

  std::cout << "\nBreit correction to Cu:\n"
               "cf Table 8.6 of Atomic Structure Theory, W. R. Johnson\n";
  for (const auto &[state, expected] : expected_B1) {
    const auto v = wf.getState(state);
    REQUIRE(v != nullptr);
    const auto de0 = *v * Vb.VbrFa(*v, wf.core());
    const auto eps = std::abs((de0 - expected) / expected);
    printf("B1: %3s %.3e [%.3e] %.0e\n", state, de0, expected, eps);
    REQUIRE(de0 == Approx(expected).epsilon(1.0e-3));
  }
}

//==============================================================================
//==============================================================================
//! integration tests for Breit
TEST_CASE("Breit", "[Breit][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Breit\n";

  // Solve Hartree-Fock, including Breit
  std::cout << "\nSolving WF, with Breit:\n";
  Wavefunction wf({3500, 1.0e-6, 125.0, 40.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  const double x_Breit = 1.0;
  wf.solve_core("HartreeFock", x_Breit, "[Xe]");
  wf.solve_valence("7sp5d");

  // Solve Hartree-Fock, without Breit
  std::cout << "\nSolving WF, without Breit:\n";
  Wavefunction wf0({3500, 1.0e-6, 125.0, 40.0, "loglinear", -1.0},
                   {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf0.solve_core("HartreeFock", 0.0, "[Xe]");
  wf0.solve_valence("7sp5d");

  // Lambda to compare against (From Vladimir's code):
  // Overall sign difference between E1, E2, and PNC ME definition (and <ab>
  // vs <ba>) taken into account Store data in list of pairs. Sort the data by
  // states, to make comparisons easier
  using datav = std::vector<std::pair<std::string, double>>;
  auto sort_by_first = [](auto x, auto y) { return x.first < y.first; };

  //============================================================================
  { // Test Breit energies:

    // Test data: From Dzuba calculation: energies (with Breit)
    auto en_VD = datav{{"1s_1/2", -1326.96780222}, {"2s_1/2", -212.26379763},
                       {"2p_1/2", -198.91648966},  {"2p_3/2", -186.09010791},
                       {"3s_1/2", -45.92577676},   {"3p_1/2", -40.36742765},
                       {"3p_3/2", -37.84531317},   {"3d_3/2", -28.28379474},
                       {"3d_5/2", -27.76338201},   {"4s_1/2", -9.50647336},
                       {"4p_1/2", -7.43335889},    {"4p_3/2", -6.91467392},
                       {"4d_3/2", -3.48506227},    {"4d_5/2", -3.39873098},
                       {"5s_1/2", -1.48931551},    {"5p_1/2", -0.90678723},
                       {"5p_3/2", -0.84005768},    {"6s_1/2", -0.12735352},
                       {"7s_1/2", -0.05518246},    {"6p_1/2", -0.08558175},
                       {"6p_3/2", -0.08377240},    {"7p_1/2", -0.04200916},
                       {"7p_3/2", -0.04136327},    {"5d_3/2", -0.06446590},
                       {"5d_5/2", -0.06458303}};
    // Test data: From Dzuba calculation: energies (without Breit)
    auto en0_VD = datav{{"1s_1/2", -1330.11855083}, {"2s_1/2", -212.56443998},
                        {"2p_1/2", -199.42946902},  {"2p_3/2", -186.43652026},
                        {"3s_1/2", -45.96973802},   {"3p_1/2", -40.44831351},
                        {"3p_3/2", -37.89428850},   {"3d_3/2", -28.30949745},
                        {"3d_5/2", -27.77513275},   {"4s_1/2", -9.51282121},
                        {"4p_1/2", -7.44628735},    {"4p_3/2", -6.92099642},
                        {"4d_3/2", -3.48561698},    {"4d_5/2", -3.39689628},
                        {"5s_1/2", -1.48980539},    {"5p_1/2", -0.90789806},
                        {"5p_3/2", -0.84033929},    {"6s_1/2", -0.12736808},
                        {"7s_1/2", -0.05518736},    {"6p_1/2", -0.08561590},
                        {"6p_3/2", -0.08378548},    {"7p_1/2", -0.04202139},
                        {"7p_3/2", -0.04136804},    {"5d_3/2", -0.06441963},
                        {"5d_5/2", -0.06452976}};
    // Sort, so don't worry about order:
    std::sort(begin(en_VD), end(en_VD), sort_by_first);
    std::sort(begin(en0_VD), end(en0_VD), sort_by_first);
    // Dzuba breit corr:
    std::vector<double> dB_VD;
    for (auto i = 0ul; i < en_VD.size(); ++i) {
      dB_VD.push_back(en_VD[i].second - en0_VD[i].second);
    }

    // Get my values, with Breit:
    datav en{};
    for (const auto &Fc : wf.core())
      en.emplace_back(Fc.symbol(), Fc.en());
    for (const auto &Fv : wf.valence())
      en.emplace_back(Fv.symbol(), Fv.en());

    // Get my values, without Breit:
    datav en0{};
    for (const auto &Fc : wf0.core())
      en0.emplace_back(Fc.symbol(), Fc.en());
    for (const auto &Fv : wf0.valence())
      en0.emplace_back(Fv.symbol(), Fv.en());
    // sort, so in same order as VD:
    std::sort(begin(en), end(en), sort_by_first);
    std::sort(begin(en0), end(en0), sort_by_first);

    // My Breit corr:
    std::cout << "\nHF Breit energy correction:\n";
    std::string worst;
    double weps = 0.0;
    for (auto i = 0ul; i < en_VD.size(); ++i) {
      std::cout << en[i].first << ": ";
      const auto dBn = en[i].second - en0[i].second;
      const auto dBn_VD = en_VD[i].second - en0_VD[i].second;
      const auto eps = std::abs(0.5 * (dBn - dBn_VD) / (dBn + dBn_VD));
      printf("%11.4e [%11.4e] ; %.1e\n", dBn, dBn_VD, eps);
      if (eps > weps) {
        weps = eps;
        worst = en[i].first;
      }
    }

    // pass &= qip::check_value(&obuff, "dE(Br) " + worst, weps, 0.0, 1.0e-3);
    REQUIRE(std::abs(weps) < 1.0e-3);
  }

  //============================================================================
  {
    // Breit to E1:
    // Test data: from Dzuba code
    // Breit corrections to absolute value of E1
    // i.e., = |<E1_Br>| - |<E1_HF>|
    auto e1_VD_HF = datav{{"6s+6p-", 0.000350633},  {"6s+6p+", 0.000777315},
                          {"6s+7p-", 0.001811049},  {"6s+7p+", 0.00060308},
                          {"7s+6p-", 0.004605428},  {"7s+6p+", 0.001841538},
                          {"7s+7p-", -0.0010422},   {"7s+7p+", 0.00075283},
                          {"6p-6s+", 0.000350633},  {"6p-7s+", 0.004605428},
                          {"6p-5d-", -0.004413645}, {"6p+6s+", 0.000777315},
                          {"6p+7s+", 0.001841538},  {"6p+5d-", -0.002804759},
                          {"6p+5d+", -0.01009528},  {"7p-6s+", 0.001811049},
                          {"7p-7s+", -0.0010422},   {"7p-5d-", -0.020635544},
                          {"7p+6s+", 0.00060308},   {"7p+7s+", 0.00075283},
                          {"7p+5d-", -0.007751954}, {"7p+5d+", -0.026316928},
                          {"5d-6p-", -0.004413645}, {"5d-6p+", -0.002804759},
                          {"5d-7p-", -0.020635544}, {"5d-7p+", -0.007751954},
                          {"5d+6p+", -0.01009528},  {"5d+7p+", -0.026316928}};

    auto e1_VD_RPA = datav{{"6s+6p-", -2.49E-07},    {"6s+6p+", 0.000135474},
                           {"6s+7p-", 0.001614732},  {"6s+7p+", 0.000309938},
                           {"7s+6p-", 0.004542807},  {"7s+6p+", 0.00184885},
                           {"7s+7p-", -0.00112601},  {"7s+7p+", 0.00057199},
                           {"6p-6s+", -2.49E-07},    {"6p-7s+", 0.004542807},
                           {"6p-5d-", -0.00488065},  {"6p+6s+", 0.000135474},
                           {"6p+7s+", 0.00184885},   {"6p+5d-", -0.003054769},
                           {"6p+5d+", -0.01087686},  {"7p-6s+", 0.001614732},
                           {"7p-7s+", -0.00112601},  {"7p-5d-", -0.020318177},
                           {"7p+6s+", 0.000309938},  {"7p+7s+", 0.00057199},
                           {"7p+5d-", -0.007615348}, {"7p+5d+", -0.02586555},
                           {"5d-6p-", -0.00488065},  {"5d-6p+", -0.003054769},
                           {"5d-7p-", -0.020318177}, {"5d-7p+", -0.007615348},
                           {"5d+6p+", -0.01087686},  {"5d+7p+", -0.02586555}};
    std::sort(begin(e1_VD_HF), end(e1_VD_HF), sort_by_first);
    std::sort(begin(e1_VD_RPA), end(e1_VD_RPA), sort_by_first);

    std::cout
        << "\nBreit corrections to E1 matrix elements, cf expected (Dzuba)\n";

    // Solve TDHF with Breit (for RPA)
    const auto h{DiracOperator::E1(wf.grid())};
    auto rpa = ExternalField::TDHF(&h, wf.vHF());
    auto rpa0 = ExternalField::TDHF(&h, wf0.vHF());
    rpa.solve_core(0.0, 20);  // w=0
    rpa0.solve_core(0.0, 20); // w=0

    // Get my values:
    datav e1_me_HF, e1_me_RPA;
    for (const auto &Fv : wf.valence()) {
      for (const auto &Fw : wf.valence()) {
        if (h.isZero(Fv.kappa(), Fw.kappa()))
          continue;
        const auto &Fv0 = *wf0.getState(Fv.n(), Fv.kappa());
        const auto &Fw0 = *wf0.getState(Fw.n(), Fw.kappa());
        // hf:
        const auto e10 = std::abs(h.reducedME(Fv0, Fw0));
        const auto e1 = std::abs(h.reducedME(Fv, Fw));
        // rpa:
        const auto e1r0 = std::abs(h.reducedME(Fv0, Fw0) + rpa0.dV(Fv0, Fw0));
        const auto e1r = std::abs(h.reducedME(Fv, Fw) + rpa.dV(Fv, Fw));
        e1_me_HF.push_back({Fv.shortSymbol() + Fw.shortSymbol(), e1 - e10});
        e1_me_RPA.push_back({Fv.shortSymbol() + Fw.shortSymbol(), e1r - e1r0});
      }
    }
    std::sort(begin(e1_me_HF), end(e1_me_HF), sort_by_first);
    std::sort(begin(e1_me_RPA), end(e1_me_RPA), sort_by_first);

    assert(e1_me_HF.size() == e1_me_RPA.size() &&
           e1_me_RPA.size() == e1_VD_RPA.size() &&
           e1_VD_RPA.size() == e1_VD_HF.size());

    // Print results, and find worst offender:
    std::string worst, worstr;
    double weps{0.0}, wepsr{0.0};
    for (auto i = 0u; i < e1_me_HF.size(); ++i) {
      std::cout << e1_me_HF[i].first << ": ";

      const auto eps =
          std::min(std::abs((e1_me_RPA[i].second - e1_VD_RPA[i].second) /
                            e1_VD_RPA[i].second),
                   std::abs(e1_me_RPA[i].second - e1_VD_RPA[i].second));
      const auto eps0 =
          std::min(std::abs((e1_me_HF[i].second - e1_VD_HF[i].second) /
                            e1_VD_HF[i].second),
                   std::abs(e1_me_HF[i].second - e1_VD_HF[i].second));

      printf("%9.6f [%9.6f] ; %9.6f [%9.6f] : %.1e\n", e1_me_HF[i].second,
             e1_VD_HF[i].second, e1_me_RPA[i].second, e1_VD_RPA[i].second,
             std::max(eps0, eps));

      if (eps0 > weps) {
        weps = eps0;
        worst = e1_me_HF[i].first;
      }
      if (eps > wepsr) {
        wepsr = eps;
        worstr = e1_me_HF[i].first;
      }
    }
    REQUIRE(std::abs(weps) < 1.0e-6);
    REQUIRE(std::abs(wepsr) < 1.0e-4);
  }

  //============================================================================
  // nb: Have to do this at the end, since Sigma will be included into
  // wavefunction
  {
    // Test Breit energies, @ Sigma(2)

    // Test data:
    // A. Derevianko, Phys. Rev. A 65, 012106 (2001).
    // Also has E1, and Hyperfine!
    // Breit corrections to energies (at HF and Sigma(2)), cf. Derevianko [in
    // cm]
    const auto de0 =
        datav{{"6s+", 3.2}, {"7s+", 1.1}, {"6p-", 7.5},   {"7p-", 2.7},
              {"6p+", 2.9}, {"7p+", 1.0}, {"5d-", -10.2}, {"5d+", -11.8}};
    const auto de2 = datav{{"6s+", -2.6}, /*{"7s+", -0.26},*/ {"6p-", 7.1},
                           {"7p-", 2.5},  {"6p+", 0.84},
                           {"7p+", 0.38}, {"5d-", -22.0},
                           {"5d+", -26.0}};

    // My values (as a regression test, and Derevianko not necisarily
    // better..) nb: if this one fails, not neccisarily an issue, BUT should
    // be checked!
    const auto de2_me = datav{{"6s+", -2.878612320},  {"7s+", 0.085538550},
                              {"6p-", 7.305198685},   {"7p-", 2.571275498},
                              {"6p+", 0.571894669},   {"7p+", 0.434694080},
                              {"5d-", -25.752413108}, {"5d+", -30.679223816}};

    // First, compare the HF energies
    std::cout << "\nBreit corrections to HF energies \ncf. "
                 "Derevianko [Phys. Rev. A 65, 012106 (2001)] (/cm)\n";
    std::string worst;
    double weps = 0.0;
    for (auto [state, dBr] : de0) {
      const auto &Fv0 = *wf0.getState(state);
      const auto &Fv = *wf.getState(state);
      const auto de = (Fv.en() - Fv0.en()) * PhysConst::Hartree_invcm;
      const auto eps = std::abs((de - dBr) / dBr);
      std::cout << state << " : ";
      printf("%6.2f [%5.1f] ; %.0e\n", de, dBr, eps);
      if (eps > weps) {
        weps = eps;
        worst = state;
      }
    }

    // Then, calculate Sigma(2), compare those
    wf.formBasis({"30spdfghi", 40, 7, 1.0e-4, 1.0e-4, 40.0, false});
    wf.formSigma(3, true, 3.0e-4, 30.0, 25 /*stride*/, false, false, {}, {}, {},
                 "false", "false");

    wf0.formBasis({"30spdfghi", 40, 7, 1.0e-4, 1.0e-4, 40.0, false});
    wf0.formSigma(3, true, 3.0e-4, 30.0, 25 /*stride*/, false, false, {}, {},
                  {}, "false", "false");

    wf.hartreeFockBrueckner();
    wf0.hartreeFockBrueckner();

    std::cout << "\nBreit corrections to Sigma(2) energies \ncf."
                 "Derevianko [Phys. Rev. A 65, 012106 (2001)] (/cm)\n";
    std::string worst2;
    double weps2 = 0.0;
    for (auto [state, dBr] : de2) {
      const auto &Fv0 = *wf0.getState(state);
      const auto &Fv = *wf.getState(state);
      const auto de = (Fv.en() - Fv0.en()) * PhysConst::Hartree_invcm;
      const auto eps = std::abs((de - dBr) / dBr);
      std::cout << state << " : ";
      printf("%6.2f [%5.1f] ; %.0e\n", de, dBr, eps);
      if (eps > weps2) {
        weps2 = eps;
        worst2 = state;
      }
    }

    std::cout << "\nBreit corrections to Sigma(2) energies cf. "
                 "me [regression test] (/cm)\n";
    std::string worstme;
    double wepsme = 0.0;
    for (auto [state, dBr] : de2_me) {
      const auto &Fv0 = *wf0.getState(state);
      const auto &Fv = *wf.getState(state);
      const auto de = (Fv.en() - Fv0.en()) * PhysConst::Hartree_invcm;
      const auto eps = std::abs((de - dBr) / dBr);
      std::cout << state << " : ";
      printf("%8.9f [%8.4f] ; %.0e\n", de, dBr, eps);
      if (eps > wepsme) {
        wepsme = eps;
        worstme = state;
      }
    }

    REQUIRE(std::abs(weps) < 0.1);
    REQUIRE(std::abs(weps2) < 0.4);
    REQUIRE(std::abs(wepsme) < 1.0e-2);
  }

  // return pass;
}
