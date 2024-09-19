#include "Wavefunction.hpp"
#include "ContinuumOrbitals.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "catch2/catch.hpp"

TEST_CASE("Wavefunction", "[wf][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Wavefunction\n";

  const unsigned long num_points = 500; // need no numerical
  const auto r0 = 1.0e-5;
  const auto rmax = 50.0;
  const auto b = 10.0;
  const auto basis_str = "5spd";
  const auto val_str = "3sp";

  Wavefunction wf({num_points, r0, rmax, b, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0, "0");
  wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-4);
  wf.formBasis({basis_str, 20, 7, 1.0e-3, 1.0e-3, 20.0, false});
  wf.solve_valence(val_str);

  wf.printCore();
  wf.printValence();
  wf.printBasis(wf.basis());

  // test copy:
  const auto wf2 = wf;
  REQUIRE(wf.core().size() == wf2.core().size());

  const auto core_en = wf.coreEnergyHF();
  REQUIRE(core_en < 0.0);

  // ensure grid set correctly, and grid shared-pointer working
  const auto grid_sp = wf.grid_sptr();
  const auto &grid = wf.grid();
  REQUIRE(&*grid_sp == &grid);
  REQUIRE(grid.num_points() == num_points);
  REQUIRE(std::abs(grid.r0() - r0) < 1.0e-10);
  REQUIRE(std::abs(grid.rmax() - rmax) < 1.0e-10);
  REQUIRE(grid.type() == GridType::loglinear);

  // ensure alpha set correctly (default)
  REQUIRE(std::abs(wf.alpha() - PhysConst::alpha) < 1.0e-10);

  // Should be set with default nuclear parameters:
  // For Na: Na-23 (Z=11, A=23)
  REQUIRE(wf.Znuc() == 11);
  REQUIRE(wf.Anuc() == 23);
  REQUIRE(std::abs(wf.get_rrms() - 2.9936) < 0.0001);
  const auto nuc = wf.nucleus();
  REQUIRE(wf.Znuc() == nuc.z());

  REQUIRE(DiracSpinor::state_config(wf.core()) == "2sp");
  REQUIRE(DiracSpinor::state_config(wf.valence()) == val_str);
  REQUIRE(DiracSpinor::state_config(wf.basis()) == basis_str);
  REQUIRE(wf.spectrum().size() == 0);

  REQUIRE(wf.vnuc().size() == num_points);
  REQUIRE(wf.vHF() != nullptr);
  REQUIRE(wf.vlocal().size() == num_points);
  // Do not have a vRad or Sigma (yet)
  REQUIRE(wf.Hmag().size() == 0);
  REQUIRE(wf.vrad() == nullptr);
  REQUIRE(wf.Sigma() == nullptr);

  // Na in V(n-1) shuold have Z-1 electrons
  REQUIRE(wf.Ncore() == wf.Znuc() - 1);

  const auto Fn = wf.getState("1s");
  REQUIRE(Fn != nullptr);
  REQUIRE(Fn->n() == 1);
  REQUIRE(Fn->kappa() == -1);
  REQUIRE(wf.isInCore(Fn->n(), Fn->kappa()));
  REQUIRE(!wf.isInValence(Fn->n(), Fn->kappa()));
  const auto Fm = wf.getState("3p+");
  REQUIRE(Fm != nullptr);
  REQUIRE(Fm->n() == 3);
  REQUIRE(Fm->kappa() == -2);
  REQUIRE(!wf.isInCore(Fm->n(), Fm->kappa()));
  REQUIRE(wf.isInValence(Fm->n(), Fm->kappa()));

  REQUIRE(wf.FermiLevel() < 0.0);
  REQUIRE(wf.energy_gap() > 0.0);

  REQUIRE(wf.coreConfiguration() == "1s2,2s2,2p6");
  REQUIRE(wf.coreConfiguration_nice() == "[Ne]");

  REQUIRE(wf.atom() == "Na, Z=11 A=23");
  REQUIRE(wf.atomicSymbol() == "Na");
  REQUIRE(wf.identity() == "Na1_0");

  REQUIRE(wf.Zion() == 1);

  const auto rho = wf.coreDensity();
  REQUIRE(rho.size() == num_points);
  const auto number_in_core_from_rho =
      NumCalc::integrate(wf.grid().du(), 0, 0, rho, wf.grid().drdu());
  REQUIRE(std::abs(number_in_core_from_rho - wf.Ncore()) < 1.0e-6);

  //============================================================================
  for (auto &v : wf.valence()) {
    REQUIRE(v.en() == Approx(wf.Hab(v, v)).epsilon(1.0e-7));
  }

  //============================================================================

  wf.formSigma(2, 2, 1.0e-3, 15.0, 8, false, false, false, 0, {}, {}, {},
               "false", "false");
  wf.hartreeFockBrueckner(true);
  wf.formSpectrum({basis_str, 20, 7, 1.0e-3, 1.0e-3, 20.0, false});
  wf.printValence();
  //
  // wf.SOEnergyShift();
  std::vector<double> br_en_list;
  for (const auto &Fv : wf.valence()) {
    br_en_list.push_back(Fv.en());
  }
  wf.fitSigma_hfBrueckner("", br_en_list);

  REQUIRE(wf.Sigma() != nullptr);

  auto Fs = wf.spectrum().front();
  const auto spectrum = wf.spectrum();
  auto [eps1, x1] = DiracSpinor::check_ortho({Fs}, wf.core());
  Fs = DiracSpinor::orthonormaliseWrt(Fs, wf.core());

  auto [eps2, x2] = DiracSpinor::check_ortho({Fs}, wf.core());
  std::cout << Fs << " " << eps1 << " " << eps2 << "\n";

  auto core = wf.core();
  auto [eps3, x3] = DiracSpinor::check_ortho(core, core);
  DiracSpinor::orthonormaliseOrbitals(core);
  auto [eps4, x4] = DiracSpinor::check_ortho(core, core);
  std::cout << eps3 << " " << eps4 << "\n";
  REQUIRE(eps4 < 1.0e-12);

  //============================================================================
  // With Brueckner
  wf.printValence();
  for (auto &v : wf.valence()) {
    REQUIRE(v.en() == Approx(wf.Hab(v, v)).epsilon(1.0e-7));
  }

  //============================================================================
  std::cout << "\nRadiative Potential inclusion (just wf unit test)\n";

  Wavefunction wfQED({num_points, r0, rmax, b, "loglinear", -1.0},
                     {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wfQED.set_HF("HartreeFock", 0.0, "[Ne]", 1.0e-4);
  wfQED.radiativePotential({1.0, 0.0, 1.0, 1.0, 1.0}, 1.0, 1.0, {1.0, 0.0},
                           false, true);
  wfQED.solve_core(true);
  wfQED.solve_valence(val_str);

  REQUIRE(wfQED.vrad() != nullptr);
  REQUIRE(wfQED.Hmag().size() == num_points);

  for (auto &v : wfQED.valence()) {
    REQUIRE(v.en() != Approx(wf.Hab(v, v)).epsilon(1.0e-7));
    REQUIRE(v.en() == Approx(wfQED.Hab(v, v)).epsilon(1.0e-7));
  }

  wfQED.radiativePotential(IO::InputBlock{}, false, true);

  REQUIRE(wfQED.vrad() != nullptr);
  REQUIRE(wfQED.Hmag().size() == num_points);

  //============================================================================
  std::cout << "\nBreit inclusion (just wf unit test)\n";

  Wavefunction wfBr({num_points, r0, rmax, b, "loglinear", -1.0},
                    {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wfBr.solve_core("HartreeFock", 1.0, "[Ne]", 1.0e-4);
  wfBr.solve_valence(val_str);

  REQUIRE(wfBr.vHF()->vBreit() != nullptr);

  for (auto &v : wfBr.valence()) {
    REQUIRE(v.en() != Approx(wf.Hab(v, v)).epsilon(1.0e-7));
    REQUIRE(v.en() == Approx(wfBr.Hab(v, v)).epsilon(1.0e-7));
  }

  //============================================================================
  std::cout << "\nH-like test:\n";
  Wavefunction hlike({2000, r0, 150, b, "loglinear", -1.0},
                     {"H", -1, "pointlike", -1.0, -1.0}, 1.0);
  hlike.solve_core("HartreeFock", 0.0, "[]"); // should do nothing!
  hlike.solve_valence("4spdf");
  hlike.printValence();

  for (auto &v : hlike.valence()) {
    REQUIRE(v.en() ==
            Approx(AtomData::diracen(1.0, v.n(), v.kappa(), hlike.alpha()))
                .epsilon(1.0e-9));
  }

  std::cout << "\nH-like test, non-relativistic limit:\n";
  Wavefunction hlike_nr({2000, r0, 150, b, "loglinear", -1.0},
                        {"H", -1, "pointlike", -1.0, -1.0}, 1.0e-9);
  // should be allowed to skip for H-like!
  // hlike_nr.solve_core("HartreeFock", 0.0, "[]");
  hlike_nr.solve_valence("4spdf");
  hlike_nr.printValence();
  for (auto &v : hlike_nr.valence()) {
    std::cout << v.en() << " "
              << AtomData::diracen(1.0, v.n(), v.kappa(), hlike.alpha()) << " "
              << AtomData::diracen(1.0, v.n(), v.kappa(), hlike_nr.alpha())
              << "\n";

    REQUIRE(v.en() !=
            Approx(AtomData::diracen(1.0, v.n(), v.kappa(), hlike.alpha()))
                .epsilon(1.0e-9));
    REQUIRE(v.en() ==
            Approx(AtomData::diracen(1.0, v.n(), v.kappa(), hlike_nr.alpha()))
                .epsilon(1.0e-9));
  }
}

//==============================================================================
TEST_CASE("Wavefunction - continuum", "[wf][cntm][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "Wavefunction - continuum\n";

  {
    std::cout << "\nNa, V^N-1\n";
    const unsigned long num_points = 500; // need no numerical
    const auto r0 = 1.0e-5;
    const auto rmax = 50.0;
    const auto b = 10.0;

    Wavefunction wf({num_points, r0, rmax, b, "loglinear", -1.0},
                    {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-8);

    std::cout << "ContinuumOrbitals\n";
    ContinuumOrbitals cntm(wf);

    cntm.solveContinuumHF(0.1, 0, 1, nullptr, false, false, false);
    const auto eps = cntm.check_orthog(true);
    REQUIRE(std::abs(eps) < 1.0e-4);

    cntm.clear();
    cntm.solveContinuumHF(0.1, 0, 1, nullptr, true, false, false);
    const auto epsX = cntm.check_orthog(true);
    REQUIRE(std::abs(epsX) < 1.0e-4);

    REQUIRE(cntm.orbitals.size() != 0);
    cntm.clear();
    REQUIRE(cntm.orbitals.size() == 0);

    cntm.solveContinuumHF(0.1, 0, 0, &wf.core().front(), false, true, true);
    cntm.check_orthog(true);
    auto ortho = std::abs(wf.core().front() * cntm.orbitals.front());
    std::cout << ortho << "\n";
    REQUIRE(std::abs(ortho) < 1.0e-15);
  }

  {
    std::cout << "\nNe, V^N\n";
    const unsigned long num_points = 500;
    const auto r0 = 1.0e-5;
    const auto rmax = 50.0;
    const auto b = 10.0;

    Wavefunction wf({num_points, r0, rmax, b, "loglinear", -1.0},
                    {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-8);

    std::cout << "ContinuumOrbitals\n";
    ContinuumOrbitals cntm(wf);

    cntm.solveContinuumHF(0.1, 0, 1, nullptr, false, false, false);
    const auto eps = cntm.check_orthog(true);
    REQUIRE(std::abs(eps) < 1.0e-4);

    cntm.clear();
    cntm.solveContinuumHF(0.1, 0, 1, nullptr, true, false, true);
    const auto epsX = cntm.check_orthog(true);
    REQUIRE(std::abs(epsX) < 1.0e-4);

    // REQUIRE(cntm.orbitals.size() != 0);
    cntm.clear();
    REQUIRE(cntm.orbitals.size() == 0);

    cntm.solveContinuumHF(0.1, 0, 0, &wf.core().front(), false, true, true);
    cntm.check_orthog(true);
    auto ortho = std::abs(wf.core().front() * cntm.orbitals.front());
    std::cout << ortho << "\n";
    REQUIRE(std::abs(ortho) < 1.0e-15);
  }

  {
    std::cout << "\nNe, V^N - ApproxHF\n";
    const unsigned long num_points = 500;
    const auto r0 = 1.0e-5;
    const auto rmax = 50.0;
    const auto b = 10.0;

    Wavefunction wf({num_points, r0, rmax, b, "loglinear", -1.0},
                    {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("ApproxHF", 0.0, "[Ne]", 1.0e-8);

    std::cout << "ContinuumOrbitals\n";
    ContinuumOrbitals cntm(wf);

    cntm.solveContinuumHF(0.1, 0, 1, nullptr, false, false, true);
    const auto eps = cntm.check_orthog(true);
    REQUIRE(std::abs(eps) < 1.0e-6);
  }

  {
    std::cout << "\nNe, V^N - KohnSham\n";
    const unsigned long num_points = 500;
    const auto r0 = 1.0e-5;
    const auto rmax = 50.0;
    const auto b = 10.0;

    Wavefunction wf({num_points, r0, rmax, b, "loglinear", -1.0},
                    {"Ne", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("KohnSham", 0.0, "[Ne]", 1.0e-8);

    std::cout << "ContinuumOrbitals\n";
    ContinuumOrbitals cntm(wf);

    cntm.solveContinuumHF(0.1, 0, 1, nullptr, false, false, false);
    const auto eps = cntm.check_orthog(true);
    REQUIRE(std::abs(eps) < 1.0e-6);
  }
}
