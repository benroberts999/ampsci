#include "HF/Breit.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
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
    CHECK(eps < 1.0e-3);
    ++count;
  }
}

//==============================================================================
TEST_CASE("Breit: Bk formulas", "[Breit][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Breit: Bk formulas\n";

  const auto radial_grid = std::make_shared<const Grid>(
      GridParameters{500, 1.0e-4, 250.0, 50.0, GridType::loglinear});
  const double zeff = 1.0;
  const int lmax = 6;

  // build set of H-like orbitals, one n for each kappa up to l=lmax
  std::vector<DiracSpinor> orbs;
  for (int l = 0; l <= lmax; ++l) {
    int n_min = l + 1;
    if (l != 0) {
      orbs.push_back(DiracSpinor::exactHlike(n_min, l, radial_grid, zeff));
    }
    orbs.push_back(DiracSpinor::exactHlike(n_min, -l - 1, radial_grid, zeff));
  }

  auto Br = HF::Breit{};

  // With scaling factors:
  const auto lamba = 0.5;
  auto BrGau = HF::Breit{};
  BrGau.update_scale(lamba, 1.0, 1.0, 0.0, 0.0);
  auto BrRet = HF::Breit{};
  BrRet.update_scale(lamba, 0.0, 0.0, 1.0, 1.0);

  {
    IO::ChronoTimer t("fill");
    Br.fill_gb(orbs);
  }

  // 1. check Symmetry:
  for (const auto &a : orbs) {
    for (const auto &b : orbs) {
      if (b <= a)
        continue;
      for (const auto &c : orbs) {
        if (c <= a)
          continue;
        for (const auto &d : orbs) {
          if (d <= a || d <= c)
            continue;
          for (int k = 0; k <= 12; ++k) {

            // Test symmetry (+ numerical error)
            const auto b1 = Br.Bk_abcd(k, a, b, c, d);
            const auto b2 = Br.Bk_abcd(k, b, a, d, c);
            const auto b3 = Br.Bk_abcd(k, c, d, a, b);
            const auto b4 = a * Br.Bkv_bcd(k, a.kappa(), b, c, d);
            const auto b5 = Br.Bk_abcd_2(k, a, b, c, d);
            REQUIRE(b1 == Approx(b2).margin(1.0e-15));
            REQUIRE(b1 == Approx(b3).margin(1.0e-15));
            REQUIRE(b1 == Approx(b4).margin(1.0e-15));
            REQUIRE(b1 == Approx(b5).margin(1.0e-15));

            // test selection rule (it's not inside Bk_abcd explicitely)
            if (b1 == 0.0) {
              REQUIRE(!HF::Breit::Bk_SR(k, a, b, c, d));
            } else {
              REQUIRE(HF::Breit::Bk_SR(k, a, b, c, d));
            }
            if (HF::Breit::Bk_SR(k, a, b, c, d)) {
              REQUIRE(k != 0);
            }

            const auto [k1, k2] = HF::Breit::k_minmax(a, b, c, d);
            if (b1 != 0) {
              REQUIRE(k >= k1);
              REQUIRE(k <= k2);
            }

            const auto [k3, k4] =
                HF::Breit::k_minmax_tj(a.twoj(), b.twoj(), c.twoj(), d.twoj());
            REQUIRE((k1 == k3 && k2 == k4));

            const auto b6 = Br.Bk_abcd(k, c, b, a, d);
            const auto b7 = Br.Bk_abcd(k, a, d, c, b);
            REQUIRE(b6 == Approx(b7).margin(1.0e-15));
            if ((a.l() + c.l() + k) % 2 == 0) {
              // M, O, P terms are anti-symmetric
              REQUIRE(b6 == Approx(-b1).margin(1.0e-15));
            } else {
              // N terms are symmetric
              REQUIRE(b6 == Approx(b1).margin(1.0e-15));
            }
            const auto s = Angular::neg1pow(a.l() + c.l() + k + 1);
            REQUIRE(b6 == Approx(s * b1).margin(1.0e-15));

            // Test Wk (and Pk) formulas
            const auto w1 = Br.BWk_abcd(k, a, b, c, d);
            auto w2 = b1;
            for (int l = 0; l <= 16; ++l) {
              w2 += Coulomb::sixj(a, c, k, b, d, l) * (2 * k + 1) *
                    Br.Bk_abcd(l, a, b, d, c);
            }
            REQUIRE(w2 == Approx(w1).margin(1.0e-15));
            const auto w3 = a * Br.BWkv_bcd(k, a.kappa(), b, c, d);
            REQUIRE(w3 == Approx(w1).margin(1.0e-15));
            const auto w4 = Br.BWk_abcd_2(k, a, b, c, d);
            REQUIRE(w4 == Approx(w1).margin(1.0e-15));

            // Test scaling factors
            if (b1 == 0.0)
              continue;
            const auto bG = BrGau.Bk_abcd(k, a, b, c, d);
            const auto bR = BrRet.Bk_abcd(k, a, b, c, d);
            REQUIRE(bG != bR);
            REQUIRE(bG != b1);
            REQUIRE(bG + bR == Approx(lamba * b1).margin(1.0e-15));
          }
        }
      }
    }
  }

  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  double total1 = 0.0;
  double total2 = 0.0;
  double t1 = 0.0;
  double t2 = 0.0;
  auto Br2 = HF::Breit{};

  Br2.fill_gb(orbs);
  {
    IO::ChronoTimer t("Old");

    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        for (const auto &c : orbs) {
          for (const auto &d : orbs) {
            for (int k = 0; k <= 12; ++k) {
              // Test symmetry (+ numerical error)
              total1 += Br.Bk_abcd(k, a, b, c, d);
            }
          }
        }
      }
    }
    t1 = t.reading_ms();
  }

  {
    IO::ChronoTimer t("new");

    for (const auto &a : orbs) {
      for (const auto &b : orbs) {
        for (const auto &c : orbs) {
          for (const auto &d : orbs) {
            for (int k = 0; k <= 12; ++k) {
              // Test symmetry (+ numerical error)
              total2 += Br2.Bk_abcd_2(k, a, b, c, d);
            }
          }
        }
      }
    }
    t2 = t.reading_ms();
  }

  std::cout << total1 << " " << total2 << "\n";
  std::cout << "Speedup: " << t1 / t2 << "\n";
  REQUIRE(t1 / t2 > 1.5);
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
    wf.formSigma(3, 3, 3.0e-4, 30.0, 25 /*stride*/, false, false, {}, {}, {},
                 "false", "false");

    wf0.formBasis({"30spdfghi", 40, 7, 1.0e-4, 1.0e-4, 40.0, false});
    wf0.formSigma(3, 3, 3.0e-4, 30.0, 25 /*stride*/, false, false, {}, {}, {},
                  "false", "false");

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

//==============================================================================
// Tests for Diragram RPA (integration/regression)
TEST_CASE("Breit: RPA Corrections",
          "[ExternalField][DiagramRPA][RPA][TDHF][Breit][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Breit: RPA Corrections\n";

  // No Breit:
  Wavefunction wf0({1600, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                   {"Cs", -1, "Fermi", -1.0, -1.0});
  wf0.solve_core("HartreeFock", 0.0, "[Xe]");
  wf0.solve_valence("7sp5d");
  wf0.formBasis({"25spd15f", 30, 7, 1.0e-4, 1.0e-6, 40.0, false});

  // With Breit:
  auto wfB = wf0;
  wfB.solve_core("HartreeFock", 1.0, "[Xe]");
  wfB.solve_valence("7sp5d");
  wfB.formBasis({"25spd15f", 30, 7, 1.0e-4, 1.0e-6, 40.0, false});

  // Format: "a", "b", Breit correction to <a|h|b>, Breit corr. to <a|h+dV|b>,
  // Derev. PRA 65 012106 (2001)
  // Corrections are as %, relative to lowest-order Hartree-Fock ME

  const auto expected_pnc = std::vector{std::tuple{"6s+", "6p-", -0.32, -0.80},
                                        {"6s+", "7p-", -0.31, -0.79},
                                        {"7s+", "6p-", -0.32, -0.81},
                                        {"7s+", "7p-", -0.31, -0.80}};
  auto h_pnc = DiracOperator::generate_pnc({}, wf0);

  const auto expected_E1 = std::vector{std::tuple{"6s+", "6p-", 0.007, 0.002},
                                       {"6s+", "6p+", 0.011, 0.005},
                                       {"6s+", "7p-", 0.48, 0.45},
                                       {"6s+", "7p+", 0.085, 0.056},
                                       {"7s+", "6p-", 0.10, 0.10},
                                       {"7s+", "6p+", 0.028, 0.028},
                                       {"7s+", "7p-", -0.010, -0.010},
                                       {"7s+", "7p+", 0.005, 0.004},
                                       {"5d-", "6p-", -0.049, -0.053},
                                       {"5d-", "6p+", -0.069, -0.074}};
  auto h_E1 = DiracOperator::generate_E1({}, wf0);

  for (auto &[h, expected] : std::vector{std::tuple{h_pnc.get(), expected_pnc},
                                         {h_E1.get(), expected_E1}}) {

    std::cout << "\nSolve RPA; using TDHF/diagram/Basis method, with/without "
                 "Breit:\n";
    std::cout << "For " << h->name() << "\n";

    // RPA, using TDHF method:
    auto tdhf_0 = ExternalField::TDHF(h, wf0.vHF());
    auto tdhf_B = ExternalField::TDHF(h, wfB.vHF());
    if (h->name() != "hfs") {
      tdhf_0.solve_core(0.0);
      tdhf_B.solve_core(0.0);
    }

    // RPA, using diaagram method:
    auto rpad_0 = ExternalField::DiagramRPA(h, wf0.basis(), wf0.vHF(), "");
    auto rpad_B = ExternalField::DiagramRPA(h, wfB.basis(), wfB.vHF(), "");
    rpad_0.solve_core(0.0);
    rpad_B.solve_core(0.0);

    // RPA, using diaagram method:
    auto rpab_0 = ExternalField::TDHFbasis(h, wf0.vHF(), wf0.basis());
    auto rpab_B = ExternalField::TDHFbasis(h, wfB.vHF(), wfB.basis());
    rpab_0.solve_core(0.0);
    rpab_B.solve_core(0.0);

    // Note: RPA(D/B) for PNC matrix elements depends stronly on basis!

    std::cout << "\nBreit corrections to " << h->name()
              << " matrix elements (%, relative to hab(HF)\n";
    std::cout << "Compare to: [Derevianko, Phys. Rev. A 65, 012106 (2001)]\n";
    std::cout << "          HF                TDHF    Diragr. Basis         \n";

    for (const auto &[a, b, d0, dRPA] : expected) {
      auto a0 = wf0.getState(a);
      auto b0 = wf0.getState(b);
      auto aB = wfB.getState(a);
      auto bB = wfB.getState(b);
      REQUIRE(a0 != nullptr);
      REQUIRE(b0 != nullptr);
      REQUIRE(aB != nullptr);
      REQUIRE(bB != nullptr);

      const auto hab_0 = h->reducedME(*a0, *b0);
      const auto hab_B = h->reducedME(*aB, *bB);

      const auto TDHF_ab_0 = hab_0 + tdhf_0.dV(*a0, *b0);
      const auto TDHF_ab_B = hab_B + tdhf_B.dV(*aB, *bB);

      const auto RPAD_ab_0 = hab_0 + rpad_0.dV(*a0, *b0);
      const auto RPAD_ab_B = hab_B + rpad_B.dV(*aB, *bB);

      const auto RPAB_ab_0 = hab_0 + rpab_0.dV(*a0, *b0);
      const auto RPAB_ab_B = hab_B + rpab_B.dV(*aB, *bB);

      // Express Breit corrections as percentage, relative to HF0 value
      const auto Breit_hf = 100.0 * (hab_B - hab_0) / hab_0;
      const auto Breit_tfhf =
          (h->name() == "hfs1") ? 0.0 : 100.0 * (TDHF_ab_B - TDHF_ab_0) / hab_0;
      const auto Breit_rpad = 100.0 * (RPAD_ab_B - RPAD_ab_0) / hab_0;
      const auto Breit_rpab = 100.0 * (RPAB_ab_B - RPAB_ab_0) / hab_0;

      printf("%3s %3s  %+.3f  [%+.3f]  %+.3f  %+.3f  %+.3f   [%+.3f]\n", a, b,
             Breit_hf, d0, Breit_tfhf, Breit_rpad, Breit_rpab, dRPA);

      // Compare Breit corrections to Derevianko, HF level, to 0.01%
      REQUIRE(Breit_hf == Approx(d0).margin(0.015));

      // Compare to Derevianko, RPA level, using TDHF: 0.1%
      REQUIRE(Breit_tfhf == Approx(dRPA).margin(0.1));

      // Compare to Derevianko, RPA level, using diagram RPA: 3%
      // nb: large margins for these caused by strong basis dependence for PNC
      CHECK(Breit_rpad == Approx(dRPA).margin(3.0));

      // Compare Breit between TDHF and RPAD
      CHECK(Breit_rpad == Approx(Breit_tfhf).margin(3.0));

      // Compare Breit between TDHF and RPA_basis
      CHECK(Breit_rpad == Approx(Breit_rpab).margin(0.1));
    }
  }
}

//==============================================================================
// Tests for Diragram RPA (integration/regression)
TEST_CASE("Breit: RPA Corrections - for HFS",
          "[ExternalField][DiagramRPA][Breit][integration][!mayfail]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Breit: RPA Corrections - for HFS\n";

  // No Breit:
  Wavefunction wf0({1000, 1.0e-6, 120.0, 0.33 * 120.0, "loglinear", -1.0},
                   {"Cs", -1, "Fermi", -1.0, -1.0});
  wf0.solve_core("HartreeFock", 0.0, "[Xe]");
  wf0.solve_valence("7sp5d");
  wf0.formBasis({"25spd15f", 30, 7, 1.0e-5, 1.0e-6, 40.0, false});

  // With Breit:
  auto wfB = wf0;
  wfB.solve_core("HartreeFock", 1.0, "[Xe]");
  wfB.solve_valence("7sp5d");
  wfB.formBasis({"25spd15f", 30, 7, 1.0e-5, 1.0e-6, 40.0, false});

  // Format: "a", "b", Breit correction to <a|h|b>, Breit corr. to <a|h+dV|b>,
  // Derev. PRA 65 012106 (2001)
  // Corrections are as %, relative to lowest-order Hartree-Fock ME
  const auto expected = std::vector{std::tuple{"6s+", "6s+", 0.0008, 0.29},
                                    {"7s+", "7s+", -0.0074, 0.27},
                                    {"6p-", "6p-", -0.42, -0.18},
                                    {"7p-", "7p-", -0.40, -0.16},
                                    {"6p+", "6p+", -0.25, 0.17},
                                    {"7p+", "7p+", -0.25, 0.19},
                                    {"5d-", "5d-", 0.54, 1.09}};

  auto h = DiracOperator::generate_hfs({}, wf0);

  std::cout << "\nSolve RPA; using diagram/Basis method, with/without Breit:\n";
  std::cout << "For " << h->name() << "\n";

  // RPA, using diagram method:
  auto rpad_0 = ExternalField::DiagramRPA(h.get(), wf0.basis(), wf0.vHF(), "");
  auto rpad_B = ExternalField::DiagramRPA(h.get(), wfB.basis(), wfB.vHF(), "");
  rpad_0.solve_core(0.0);
  rpad_B.solve_core(0.0);

  // RPA, using basis method:
  auto rpab_0 = ExternalField::TDHFbasis(h.get(), wf0.vHF(), wf0.basis());
  auto rpab_B = ExternalField::TDHFbasis(h.get(), wfB.vHF(), wfB.basis());
  rpab_0.solve_core(0.0);
  rpab_B.solve_core(0.0);

  std::cout << "\nBreit corrections to " << h->name()
            << " matrix elements (%, relative to hab(HF)\n";
  std::cout << "Compare to: [Derevianko, Phys. Rev. A 65, 012106 (2001)]\n";
  std::cout << "          HF                Diragr. Basis         \n";

  // bool passed = true;
  for (const auto &[a, b, d0, dRPA] : expected) {
    auto a0 = wf0.getState(a);
    auto b0 = wf0.getState(b);
    auto aB = wfB.getState(a);
    auto bB = wfB.getState(b);
    REQUIRE(a0 != nullptr);
    REQUIRE(b0 != nullptr);
    REQUIRE(aB != nullptr);
    REQUIRE(bB != nullptr);

    const auto hab_0 = h->reducedME(*a0, *b0);
    const auto hab_B = h->reducedME(*aB, *bB);

    const auto RPAD_ab_0 = hab_0 + rpad_0.dV(*a0, *b0);
    const auto RPAD_ab_B = hab_B + rpad_B.dV(*aB, *bB);

    const auto RPAB_ab_0 = hab_0 + rpab_0.dV(*a0, *b0);
    const auto RPAB_ab_B = hab_B + rpab_B.dV(*aB, *bB);

    // Express Breit corrections as percentage, relative to HF0 value
    const auto Breit_hf = 100.0 * (hab_B - hab_0) / hab_0;
    const auto Breit_rpad = 100.0 * (RPAD_ab_B - RPAD_ab_0) / hab_0;
    const auto Breit_rpab = 100.0 * (RPAB_ab_B - RPAB_ab_0) / hab_0;

    printf("%3s %3s  %+.3f  [%+.3f]  %+.3f  %+.3f   [%+.3f]\n", a, b, Breit_hf,
           d0, Breit_rpad, Breit_rpab, dRPA);

    // Compare Breit corrections to Derevianko, HF level, to 0.01%
    REQUIRE(Breit_hf == Approx(d0).margin(0.05));
    // Compare to Derevianko, RPA level, using diagram RPA: 3%
    CHECK(Breit_rpad == Approx(dRPA).margin(0.2)); // FAILS??
    // Compare Breit between TDHF and RPAD
    REQUIRE(Breit_rpad == Approx(Breit_rpab).margin(0.1));
  }
  // REQUIRE(passed); // print table even if fails
}
