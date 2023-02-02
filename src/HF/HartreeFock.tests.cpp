#include "DiracOperator/DiracOperator.hpp"
#include "HF/HartreeFock_test_data.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <string>
#include <tuple>

//! Unit tests for Hartree Fock equations
TEST_CASE("HartreeFock", "[HF][HartreeFock][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "HartreeFock\n";

  //============================================================================

  // Grid parameters (etc.):
  const auto grid_type = "loglinear";
  const auto points = 4000;
  const auto r0 = 1.0e-6;
  const auto rmax = 150.0;
  // const auto points = 10000;
  // const auto r0 = 1.0e-7;
  // const auto rmax = 170.0;
  const auto b = 0.3 * rmax;
  const auto x_Breit = 0.0; // do not include Breit
  const int A = -1;         // use default A
  const auto nucleus_type = "Fermi";

  // Regression test: Compare Energies, E1 and HFS to my own calcs (test data
  // run using very dense radial grids)
  //--------------------------------------------------------------------------
  {
    double worst_eps{0.0};
    std::string worst_case{""};

    double wE1_eps{0.0};
    std::string wE1_case{""};

    double wHFSsp_eps{0.0};
    std::string wHFSsp_case{""};

    double wHFSdf_eps{0.0};
    std::string wHFSdf_case{""};

    // Loop through each test case, run HF, compare to test data
    for (auto &[Atom, Core, Valence, EnergyData, E1Data, HFSData] :
         UnitTest::HF_test_data::regression_test_data) {

      Wavefunction wf({points, r0, rmax, b, grid_type},
                      {Atom, A, nucleus_type});

      const auto h = DiracOperator::HyperfineA(
          1.0, 0.5, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());
      const auto d = DiracOperator::E1(wf.grid());

      std::cout << "\n" << wf.atom() << "\n";
      wf.solve_core("HartreeFock", x_Breit, Core);
      wf.solve_valence(Valence);

      // // For generating test data:
      // for (auto &Fc : wf.core()) {
      //   printf("{\"%s\", %.12f},\n", Fc.shortSymbol().c_str(), Fc.en());
      // }

      // Test the energies:
      for (const auto &[state, energy] : EnergyData) {
        const auto &Fv = *wf.getState(state);
        const auto del = std::abs(Fv.en() - energy);
        const auto eps = std::abs(del / energy);
        const auto err = std::min(del, eps);
        printf("%3s %9.6f [%9.6f] %.1e\n", state, Fv.en(), energy, err);
        if (err > worst_eps) {
          worst_eps = err;
          worst_case = wf.atomicSymbol() + ":" + state;
        }
      }

      // Test the E1 matrix elements (no RPA):
      for (const auto &[fa, fb, e1] : E1Data) {
        const auto &Fa = *wf.getState(fa);
        const auto &Fb = *wf.getState(fb);
        const auto e1_me = d.reducedME(Fa, Fb);
        const auto del = std::abs(e1_me - e1);
        const auto eps = std::abs((e1_me - e1) / e1);
        const auto err = std::min(del, eps);
        printf("E1:%3s,%3s %9.5f [%9.5f] %.1e\n", fa, fb, e1_me, e1, err);
        // nb: don't test very small MEs - large eps despite high accuracy
        if (err > wE1_eps) {
          wE1_eps = err;
          wE1_case = wf.atomicSymbol() + ":" + fa + "," + fb;
        }
      }

      // Test the HFS constants (no RPA):
      for (const auto &[fv, Ahfs] : HFSData) {
        const auto &Fv = *wf.getState(fv);
        const auto Ahfs_me = h.hfsA(Fv);
        const auto eps = std::abs((Ahfs_me - Ahfs) / Ahfs);
        printf("HFS:%3s %11.5e [%11.5e] %.1e\n", fv, Ahfs_me, Ahfs, eps);
        if (Fv.l() <= 1 && eps > wHFSsp_eps) {
          wHFSsp_eps = eps;
          wHFSsp_case = wf.atomicSymbol() + ":" + fv;
        }
        if (Fv.l() > 1 && eps > wHFSdf_eps) {
          wHFSdf_eps = eps;
          wHFSdf_case = wf.atomicSymbol() + ":" + fv;
        }
      }
    }

    std::cout << "\nWorst cases:\n";
    std::cout << "En : " << worst_case << " " << worst_eps << "\n";
    std::cout << "E1 : " << wE1_case << " " << wE1_eps << "\n";
    std::cout << "HFSsp: " << wHFSsp_case << " " << wHFSsp_eps << "\n";
    std::cout << "HFSdf: " << wHFSdf_case << " " << wHFSdf_eps << "\n";

    REQUIRE(std::abs(worst_eps) < 3.0e-6);
    REQUIRE(std::abs(wE1_eps) < 1.0e-5);
    REQUIRE(std::abs(wHFSsp_eps) < 1.0e-5);
    REQUIRE(std::abs(wHFSdf_eps) < 1.0e-4);
  }

  // Accuracy test: test HF energies against Dzuba code
  //--------------------------------------------------------------------------
  {
    double worst_eps{0.0};
    std::string worst_case{""};

    // Loop through each test case, run HF, compare to test data
    for (auto &[Atom, Core, Valence, Data] :
         UnitTest::HF_test_data::compare_VD) {

      Wavefunction wf({points, r0, rmax, b, grid_type},
                      {Atom, A, nucleus_type});

      std::cout << "\n" << wf.atom() << "\n";
      wf.solve_core("HartreeFock", x_Breit, Core);
      wf.solve_valence(Valence);

      // test energies:
      for (auto &[state, energy] : Data) {
        const auto &Fv = *wf.getState(state);
        const auto eps = std::abs((Fv.en() - energy) / energy);
        if (eps > worst_eps) {
          worst_eps = eps;
          worst_case = wf.atomicSymbol() + ":" + state;
        }
        printf("%3s %9.6f [%9.6f] %.1e\n", state, Fv.en(), energy, eps);
      }
    }
    std::cout << worst_case << " " << worst_eps << "\n";

    REQUIRE(std::abs(worst_eps) < 1.0e-5);
  }
}

// NB: Though the following are not really "unit" tests, HF
// is so fundamental, these can be considered "unit" tests.
// These are pretty quick too.

// Same as above, but only for Cs, and with fewer points - quick for "unit"
// tests part
//============================================================================
TEST_CASE("HartreeFock - just Cs", "[HF][HartreeFock][Breit][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "HartreeFock - just Cs\n";

  //============================================================================

  for (auto method :
       {HF::Method::HartreeFock, HF::Method::ApproxHF, HF::Method::Hartree,
        HF::Method::KohnSham, HF::Method::Local}) {

    auto m_str = HF::parseMethod(method);
    REQUIRE(HF::parseMethod(m_str) == method);
  }

  // Grid parameters (etc.):
  const auto grid_type = "loglinear";
  const auto points = 2000;
  const auto r0 = 1.0e-6;
  const auto rmax = 150.0;
  const auto b = 0.3 * rmax;
  const auto x_Breit = 0.0; // do not include Breit
  const int A = -1;         // use default A
  const auto nucleus_type = "Fermi";

  double worst_eps{0.0};
  double wE1_eps{0.0};
  double wHFSsp_eps{0.0};
  double wHFSdf_eps{0.0};

  auto Cs_data = *std::find_if(
      UnitTest::HF_test_data::regression_test_data.begin(),
      UnitTest::HF_test_data::regression_test_data.end(),
      [](auto &d) { return std::get<0>(d) == std::string("Cs"); });

  const auto &[Atom, Core, Valence, EnergyData, E1Data, HFSData] = Cs_data;

  Wavefunction wf({points, r0, rmax, b, grid_type}, {Atom, A, nucleus_type});

  const auto h = DiracOperator::HyperfineA(
      1.0, 0.5, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());
  const auto d = DiracOperator::E1(wf.grid());

  std::cout << "\n" << wf.atom() << "\n";
  wf.solve_core("HartreeFock", x_Breit, Core);
  wf.solve_valence(Valence);

  // Test the energies:
  for (const auto &[state, energy] : EnergyData) {
    const auto &Fv = *wf.getState(state);
    const auto del = std::abs(Fv.en() - energy);
    const auto eps = std::abs(del / energy);
    const auto err = std::min(del, eps);
    printf("%3s %9.6f [%9.6f] %.1e\n", state, Fv.en(), energy, err);
    if (err > worst_eps) {
      worst_eps = err;
    }
  }

  // Test the E1 matrix elements (no RPA):
  for (const auto &[fa, fb, e1] : E1Data) {
    const auto &Fa = *wf.getState(fa);
    const auto &Fb = *wf.getState(fb);
    const auto e1_me = d.reducedME(Fa, Fb);
    const auto del = std::abs(e1_me - e1);
    const auto eps = std::abs((e1_me - e1) / e1);
    const auto err = std::min(del, eps);
    printf("E1:%3s,%3s %9.5f [%9.5f] %.1e\n", fa, fb, e1_me, e1, err);
    // nb: don't test very small MEs - large eps despite high accuracy
    if (err > wE1_eps) {
      wE1_eps = err;
    }
  }

  // Test the HFS constants (no RPA):
  for (const auto &[fv, Ahfs] : HFSData) {
    const auto &Fv = *wf.getState(fv);
    const auto Ahfs_me = h.hfsA(Fv);
    const auto eps = std::abs((Ahfs_me - Ahfs) / Ahfs);
    printf("HFS:%3s %11.5e [%11.5e] %.1e\n", fv, Ahfs_me, Ahfs, eps);
    if (Fv.l() <= 1 && eps > wHFSsp_eps) {
      wHFSsp_eps = eps;
    }
    if (Fv.l() > 1 && eps > wHFSdf_eps) {
      wHFSdf_eps = eps;
    }
  }

  REQUIRE(std::abs(worst_eps) < 1.0e-5);
  REQUIRE(std::abs(wE1_eps) < 1.0e-4);
  REQUIRE(std::abs(wHFSsp_eps) < 1.0e-4);
  REQUIRE(std::abs(wHFSdf_eps) < 1.0e-3);

  // Breit (regression test)
  std::cout << "Breit regression test: de\n";
  HF::Breit Vbr{1.0};
  // generated with large # points.
  const auto breit_data = std::vector{
      std::tuple{"6s+", 1.333977079113e-04}, {"7s+", 3.660284833750e-05},
      {"6p-", 6.839614663249e-05},           {"7p-", 2.450254202039e-05},
      {"6p+", 4.938700048099e-05},           {"7p+", 1.785296077599e-05},
      {"5d-", 5.683556145298e-05},           {"5d+", 4.292583833884e-05},
      {"4f-", 1.023721519428e-08},           {"4f+", 6.291692462117e-09}};
  REQUIRE(Vbr.scale_factor() == 1.0);
  for (auto &[state, de_t] : breit_data) {
    const auto &Fv = *wf.getState(state);
    const auto de = Fv * Vbr.VbrFa(Fv, wf.core());
    const auto eps = std::abs((de - de_t) / de_t);
    // printf("{\"%s\", %.12e}\n,", state, de);
    printf("%3s : %.5e [%.5e] %.1e\n", Fv.shortSymbol().c_str(), de, de_t, eps);
    REQUIRE(eps < 1.0e-3);
  }
}

//============================================================================
TEST_CASE("HartreeFock - KS Core-Hartree and ApproxHF",
          "[HF][HartreeFock][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "HartreeFock - KS Core-Hartree and ApproxHF\n";

  // Grid parameters (etc.):
  const auto grid_type = "loglinear";
  const auto points = 4000;
  const auto r0 = 1.0e-6;
  const auto rmax = 150.0;
  // const auto points = 10000;
  // const auto r0 = 1.0e-7;
  // const auto rmax = 170.0;
  const auto b = 0.3 * rmax;
  const auto nucleus_type = "Fermi";

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // For KS, Core-Hartree and ApproxHF
  //--------------------------------------------------------------------------
  {
    Wavefunction wf({points, r0, rmax, b, grid_type},
                    {"Cs", 133, nucleus_type});

    std::cout << "\n" << wf.atom() << ": Core-Hartree:\n";
    wf.solve_core("Hartree", 0.0, "[Xe]");
    wf.solve_valence("6sp");

    const auto CH_data = std::vector{std::tuple{"1s+", -1294.301178790881},
                                     {"2s+", -200.035308851470},
                                     {"2p-", -187.406456454656},
                                     {"2p+", -174.902524561300},
                                     {"3s+", -40.521041198710},
                                     {"3p-", -35.394491225319},
                                     {"3p+", -33.047465684328},
                                     {"3d-", -24.111351801170},
                                     {"3d+", -23.623370439065},
                                     {"4s+", -7.247914593744},
                                     {"4p-", -5.462936127458},
                                     {"4p+", -4.986496836797},
                                     {"4d-", -2.147711860432},
                                     {"4d+", -2.073607194631},
                                     {"5s+", -0.785050228744},
                                     {"5p-", -0.389407938461},
                                     {"5p+", -0.342378387138},
                                     {"6s+", -0.120056417166},
                                     {"6p-", -0.081300533231},
                                     {"6p+", -0.078769175539}};

    double worst_eps{0.0};
    std::string worst_case{""};
    // Test the energies:
    for (const auto &[state, energy] : CH_data) {
      const auto &Fv = *wf.getState(state);
      const auto eps = std::abs((Fv.en() - energy) / energy);
      printf("%3s %9.6f [%9.6f] %.1e\n", state, Fv.en(), energy, eps);
      if (eps > worst_eps) {
        worst_eps = eps;
        worst_case = wf.atomicSymbol() + ":" + state;
      }
    }

    // pass &= qip::check_value(&obuff, "CH regression, En: " + worst_case,
    //                          worst_eps, 0.0, 3.0e-5);
    REQUIRE(std::abs(worst_eps) < 3.0e-5);
  }
  //--------------------------------------------------------------------------
  {
    // nb: Kohn-Sham seems quite grid dependent.. set b=-1 => log grid
    Wavefunction wf({points, r0, rmax, 0, "logarithmic"},
                    {"Cs", 133, nucleus_type});

    std::cout << "\n" << wf.atom() << ": Kohn-Sham:\n";
    wf.solve_core("KohnSham", 0.0, "[Xe],6s1", 1.0e-13);
    wf.solve_valence("6sp5d");

    const auto KS_data = std::vector{std::tuple{"1s+", -1312.974817713126},
                                     {"2s+", -205.876011196460},
                                     {"2p-", -193.780516952252},
                                     {"2p+", -180.853127701472},
                                     {"3s+", -42.750869442020},
                                     {"3p-", -37.686361216099},
                                     {"3p+", -35.210672201287},
                                     {"3d-", -26.273459503137},
                                     {"3d+", -25.747297666300},
                                     {"4s+", -8.030659588002},
                                     {"4p-", -6.212801652799},
                                     {"4p+", -5.703625002649},
                                     {"4d-", -2.743417372249},
                                     {"4d+", -2.659126008550},
                                     {"5s+", -0.962265397955},
                                     {"5p-", -0.519598486174},
                                     {"5p+", -0.459140168906},
                                     {"6s+", -0.124015002195},
                                     {"6s+", -0.124015002195},
                                     {"6p-", -0.084886590220},
                                     {"6p+", -0.082779593449},
                                     {"5d-", -0.060460563100},
                                     {"5d+", -0.060149046241}};

    double worst_eps{0.0};
    std::string worst_case{""};
    // Test the energies:
    for (const auto &[state, energy] : KS_data) {
      const auto &Fv = *wf.getState(state);
      const auto eps = std::abs((Fv.en() - energy) / energy);
      printf("%3s %9.6f [%9.6f] %.1e\n", state, Fv.en(), energy, eps);
      if (eps > worst_eps) {
        worst_eps = eps;
        worst_case = wf.atomicSymbol() + ":" + state;
      }
    }

    // pass &= qip::check_value(&obuff, "KS regression, En: " + worst_case,
    //                          worst_eps, 0.0, 1.0e-5);
    REQUIRE(std::abs(worst_eps) < 1.0e-5);
  }
  //--------------------------------------------------------------------------
  {
    Wavefunction wf({points, r0, rmax, b, grid_type},
                    {"Cs", 133, nucleus_type});

    std::cout << "\n" << wf.atom() << ": ApproxHF:\n";
    wf.solve_core("ApproxHF", 0.0, "[Xe]");
    wf.solve_valence("6sp");

    const auto aHF_data = std::vector{std::tuple{"1s+", -1330.095428346559},
                                      {"2s+", -212.549839273448},
                                      {"2p-", -199.412267819184},
                                      {"2p+", -186.421440251763},
                                      {"3s+", -45.963541804192},
                                      {"3p-", -40.441520330164},
                                      {"3p+", -37.888138071025},
                                      {"3d-", -28.302116154339},
                                      {"3d+", -27.768071917649},
                                      {"4s+", -9.511433719799},
                                      {"4p-", -7.444887195937},
                                      {"4p+", -6.919796669697},
                                      {"4d-", -3.484491541740},
                                      {"4d+", -3.395859554961},
                                      {"5s+", -1.489655369673},
                                      {"5p-", -0.907778905846},
                                      {"5p+", -0.840244318069},
                                      {"6s+", -0.127364872198},
                                      {"6p-", -0.085613564779},
                                      {"6p+", -0.083783585251}};

    double worst_eps{0.0};
    std::string worst_case{""};
    // Test the energies:
    for (const auto &[state, energy] : aHF_data) {
      const auto &Fv = *wf.getState(state);
      const auto eps = std::abs((Fv.en() - energy) / energy);
      printf("%3s %9.6f [%9.6f] %.1e\n", state, Fv.en(), energy, eps);
      if (eps > worst_eps) {
        worst_eps = eps;
        worst_case = wf.atomicSymbol() + ":" + state;
      }
    }

    // pass &= qip::check_value(&obuff, "aHF regression, En: " + worst_case,
    //                          worst_eps, 0.0, 1.0e-5);
    REQUIRE(std::abs(worst_eps) < 1.0e-5);
  }
}

//============================================================================
TEST_CASE("HartreeFock - Hyperfine", "[HF][HartreeFock]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "HartreeFock - Hyperfine\n";

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Test hyperfine constants (HF; no RPA) against:
  // Grunefeld, Roberts, Ginges, Phys. Rev. A 100, 042506 (2019).
  // nb: only few digits given there
  //----------------------------------------------------------------------------
  std::cout << "\nHyperfine, cf Grunefeld Phys. Rev. A 100, 042506 (2019):\n";
  {
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
        2.751818, 1.5, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());

    // Calculate HFS A constant for each valence state (store s,p seperately)
    std::vector<double> sme, pme;
    for (const auto &Fv : wf.valence()) {
      if (Fv.kappa() == -1) {
        sme.push_back(h.hfsA(Fv));
      } else if (Fv.kappa() == 1) {
        pme.push_back(h.hfsA(Fv));
      }
    }

    std::cout << "Rb\n";
    assert(p.size() == s.size());
    for (auto i = 0ul; i < s.size(); ++i) {
      const auto epss = std::abs((sme.at(i) - s.at(i)) / s.at(i));
      const auto epsp = std::abs((pme.at(i) - p.at(i)) / p.at(i));
      printf("s %.3e [%.3e] %.0e\n", sme.at(i), s.at(i), epss);
      printf("p %.3e [%.3e] %.0e\n", pme.at(i), p.at(i), epsp);
    }
    const auto [es, ats] = qip::compare_eps(sme, s);
    const auto [ep, atp] = qip::compare_eps(pme, p);
    // pass &= qip::check_value(&obuff, "HF hfs Rb", qip::min_abs(es, ep), 0.0,
    //                          3.0e-4);
    REQUIRE(qip::min_abs(es, ep) < 3.0e-4);
  }

  //----------------------------------------------------------------------------
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
        2.582025, 3.5, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());

    std::vector<double> sme, pme;
    for (const auto &Fv : wf.valence()) {
      if (Fv.kappa() == -1)
        sme.push_back(h.hfsA(Fv));
      else if (Fv.kappa() == 1)
        pme.push_back(h.hfsA(Fv));
    }
    std::cout << "Cs\n";
    assert(p.size() == s.size());
    for (auto i = 0ul; i < s.size(); ++i) {
      const auto epss = std::abs((sme.at(i) - s.at(i)) / s.at(i));
      const auto epsp = std::abs((pme.at(i) - p.at(i)) / p.at(i));
      printf("s %.3e [%.3e] %.0e\n", sme.at(i), s.at(i), epss);
      printf("p %.3e [%.3e] %.0e\n", pme.at(i), p.at(i), epsp);
    }
    const auto [es, ats] = qip::compare_eps(sme, s);
    const auto [ep, atp] = qip::compare_eps(pme, p);
    // pass &= qip::check_value(&obuff, "HF hfs Cs", qip::min_abs(es, ep), 0.0,
    //                          3.0e-4);
    REQUIRE(qip::min_abs(es, ep) < 3.0e-4);
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
        4.00, 4.5, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());

    std::vector<double> sme, pme;
    for (const auto &Fv : wf.valence()) {
      if (Fv.kappa() == -1)
        sme.push_back(h.hfsA(Fv));
      else if (Fv.kappa() == 1)
        pme.push_back(h.hfsA(Fv));
    }

    std::cout << "Fr\n";
    assert(p.size() == s.size());
    for (auto i = 0ul; i < s.size(); ++i) {
      const auto epss = std::abs((sme.at(i) - s.at(i)) / s.at(i));
      const auto epsp = std::abs((pme.at(i) - p.at(i)) / p.at(i));
      printf("s %.3e [%.3e] %.0e\n", sme.at(i), s.at(i), epss);
      printf("p %.3e [%.3e] %.0e\n", pme.at(i), p.at(i), epsp);
    }
    const auto [es, ats] = qip::compare_eps(sme, s);
    const auto [ep, atp] = qip::compare_eps(pme, p);
    // pass &= qip::check_value(&obuff, "HF hfs Fr", qip::min_abs(es, ep), 0.0,
    //                          3.0e-4);
    REQUIRE(qip::min_abs(es, ep) < 3.0e-4);
  }
}
