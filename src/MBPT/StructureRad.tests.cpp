#include "StructureRad.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"
#include <algorithm>
#include <map>
#include <string>

//==============================================================================
TEST_CASE("MBPT: Structure Rad + Norm, basic", "[StrucRad][MBPT][unit]") {

  // note: does not test formulas: just checks class is working correctly.
  // Other (integration) tests below check correctness/accuracy of formulas

  Wavefunction wf({300, 1.0e-3, 30.0, 10.0, "loglinear", -1.0},
                  {"B", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", 0.0, "[He],2s2", 1.0e-12, false);
  SplineBasis::Parameters params{"3spd4f", 60, 5, 1.0e-2, 1.0e-2, 20.0};
  params.verbose = false;
  wf.formBasis(params);

  MBPT::StructureRad srn(wf.basis(), wf.FermiLevel(), {0, 999}, "", {}, {},
                         true);

  REQUIRE(srn.core().size() == wf.core().size());
  REQUIRE(srn.excited().size() == wf.basis().size() - wf.core().size());
  REQUIRE(srn.basis().size() == wf.basis().size());
  REQUIRE(srn.me_Table().empty());

  // Do for two operators, different rank and parity
  const auto h1 = DiracOperator::E1(wf.grid());
  const auto h2 = DiracOperator::Ek(wf.grid(), 2);

  const auto rand_str = qip::random_string(4);
  const auto QK_fname = "deleteme_" + rand_str + ".qk.abf";

  for (const auto h :
       std::vector<const DiracOperator::TensorOperator *>{&h1, &h2}) {

    srn.solve_core(h);
    REQUIRE(!srn.me_Table().empty());

    // test srn_table and me_table
    std::cout << "Test SRN, SRN_table, and meTable: " << h->name() << "\n";
    const auto tab = srn.srn_table(h, srn.core(), srn.excited());
    const auto &me_tab = srn.me_Table();
    for (const auto &a : srn.core()) {
      for (const auto &b : srn.excited()) {
        if (h->isZero(a, b))
          continue;
        // test srn_table():
        const auto tab_ab = tab.get(a, b);
        const auto tab_ba = tab.get(b, a);
        REQUIRE(tab_ab != nullptr);
        REQUIRE(tab_ba != nullptr);
        const auto srab = *tab_ab;
        const auto srba = *tab_ba;

        const auto relative_sign = h->symm_sign(a, b);
        REQUIRE(srab == Approx(srba * relative_sign));

        // test srn():
        const auto norm = srn.norm(a, b, h);
        const auto srn0 = srn.SR(a, b) + norm;
        REQUIRE(srab == Approx(srn0));

        // test underlying meTable was filled correctly
        const auto me_ab_tab = me_tab.getv(a, b);
        const auto me_ab_cal = h->reducedME(a, b);
        REQUIRE(me_ab_tab == Approx(me_ab_cal));

        // Test norm direct vs. fia factors
        const auto f_norm_a = srn.f_norm(a);
        const auto f_norm_b = srn.f_norm(b);
        const auto norm2 = (f_norm_a + f_norm_b) * (me_ab_tab);
        REQUIRE(norm == Approx(norm2));
      }
    }

    std::cout << "Test SRN: Qk vs Yk table, including read/write and copy: "
              << h->name() << "\n";
    // TEST SRN read-write using Qk_file

    // Create Qk file:
    auto print = h == &h1; // only once
    MBPT::StructureRad srn2(wf.basis(), wf.FermiLevel(), {0, 999}, QK_fname, {},
                            {}, print);
    // Read in QK file:
    MBPT::StructureRad srn3(wf.basis(), wf.FermiLevel(), {0, 999}, QK_fname,
                            {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, print);
    srn2.solve_core(h);
    srn3.solve_core(h);
    // Copy existing
    const auto srn4 = srn3;
    for (const auto &a : wf.basis()) {
      for (const auto &b : srn.excited()) {
        if (a < b && !h->isZero(a, b)) {
          // test all three versions (direct Yk, direct Qk, read-in Qk)
          const auto srC0 = srn.SR(a, b);
          const auto srC2 = srn2.SR(a, b);
          const auto srC3 = srn3.SR(a, b);
          REQUIRE(srC0 == Approx(srC2));
          REQUIRE(srC0 == Approx(srC3));

          const auto srN2 = srn2.norm(a, b, h);
          const auto srN3 = srn3.norm(a, b, h);
          REQUIRE(srN2 == Approx(srN3));

          const auto srn_tot = srn2.srn(a, b, h);
          REQUIRE(srn_tot == Approx(srC2 + srN2));

          // test copy-constructor version mathces
          const auto srn_tot4 = srn4.srn(a, b, h);
          REQUIRE(srn_tot == Approx(srn_tot4));
        }
      }
    }
  }
}

//==============================================================================
TEST_CASE("MBPT: Structure Rad + Norm",
          "[StrucRad][MBPT][ExternalField][integration][slow]") {

  // Test for Na: (use splines for legs)
  Wavefunction wf({800, 1.0e-6, 90.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-12, false);
  wf.solve_valence("4s3p", false);
  // wf.printValence();
  auto basis_params =
    SplineBasis::Parameters{"25spdfgh", 30, 7, 1.0e-4, 0, 40.0};
  // SplineBasis::Parameters{"40spdfghi", 42, 7, 1.0e-4, 0, 40.0};
  basis_params.verbose = false;
  wf.formBasis(basis_params);

  // Find core/valence energy: allows distingush core/valence states
  const auto en_core = wf.FermiLevel();

  const auto h = DiracOperator::E1(wf.grid());

  MBPT::StructureRad srn(wf.basis(), en_core);
  srn.solve_core(&h);

  // Expected data, from:
  // Johnson et al, At.Dat.Nuc.Dat.Tables 64, 279 (1996), 10.1006/adnd.1996.0024
  // Table E (Na)
  // Z(1), RPA, BO, SR, Norm
  const std::map<std::string, std::array<double, 5>> JohnsonData = {
    {"3p-3s+", {3.6906, -0.0385 - 0.0047, -0.1020, +0.0030, -0.0050}},
    {"3p+3s+", {5.2188, -0.0544 - 0.0066, -0.1443, +0.0043, -0.0070}},
    {"3p-4s+", {3.6004, +0.0068 + 0.0007, -0.0185, -0.0010, -0.0021}},
    {"3p+4s+", {5.1012, +0.0095 + 0.0010, -0.0254, -0.0014, -0.0029}},
  };

  std::cout << "\nStructure Radiation + Norm of states.\n"
               "E1 for Na, cf Johnson et al. (1996)\n";

  std::cout << "(nb: I use quite a smaller basis for speed; "
               "get closer agreement with larger basis. "
               "BO is particularly sensitive.)\n";
  std::cout << "Basis here   : " << DiracSpinor::state_config(wf.basis())
            << "\n";
  std::cout << "Basis Johnson: 40spdfghik\n\n";

  std::cout << "Table E [At. Data Nucl. Data Tables 64, 279 (1996)]\n";
  fmt::print("{:3} {:3} :  {:7}  {:7}  {:7}  {:7}\n", "w", "v", " Z(1)", " BO",
             " SR", " Norm");
  for (const auto &v : wf.valence()) {
    for (const auto &w : wf.valence()) {
      // v -> w <w||v>
      const auto key = w.shortSymbol() + v.shortSymbol();
      const auto data_ptr = JohnsonData.find(key);
      if (data_ptr == JohnsonData.end())
        continue;
      const auto [Z1, RPA, BO, SR, NORM] = data_ptr->second;

      const auto z1_tmp = h.reducedME(w, v);
      // Johnson data: z1 is always positive
      const auto s = z1_tmp / std::abs(z1_tmp);

      // include omega, though it makes no difference!
      const double omega = 0.0; //w.en() - v.en();

      const auto z1 = s * z1_tmp;
      const auto bo = s * srn.BO(w, v);
      const auto sr = s * srn.SR(w, v, omega);
      const auto norm = s * srn.norm(w, v, &h, nullptr);

      fmt::print("{:3} {:3} :  {:+.4f}  {:+.4f}  {:+.4f}  {:+.4f}\n",
                 w.shortSymbol(), v.shortSymbol(), z1, bo, sr, norm);
      fmt::print("        [  {:+.4f}  {:+.4f}  {:+.4f}  {:+.4f} ]\n", Z1, BO,
                 SR, NORM);

      CHECK(z1 == Approx(Z1).margin(0.0001));
      CHECK(bo == Approx(BO).margin(0.0030)); // sensitive to basis
      CHECK(sr == Approx(SR).margin(0.0002));
      CHECK(norm == Approx(NORM).margin(0.0002));
    }
  }
}

//==============================================================================
TEST_CASE("MBPT: Structure Rad + Norm HFS",
          "[StrucRad][MBPT][ExternalField][integration][slow]") {

  // Test for Na: (use splines for legs)
  Wavefunction wf({800, 1.0e-6, 80.0, 0.33 * 100.0, "loglinear"},
                  {"Cs", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]", 1.0e-12, false);
  wf.solve_valence("6p", false);
  // wf.printValence();
  auto basis_params =
    SplineBasis::Parameters{"20spdfgh", 30, 7, 1.0e-4, 0, 40.0};
  // SplineBasis::Parameters{"40spdfghi", 42, 7, 1.0e-4, 0, 40.0};
  basis_params.verbose = false;
  wf.formBasis(basis_params);

  std::cout << "\nStructure Radiation + Norm of states.\n"
               "HFS for Cs, cf Johnson et al. (2004)\n";

  std::cout << "(nb: I use quite a small basis here, so not exact)\n\n";

  // Expected data, from:
  // Johnson et al, PhysRevA 70, 014501 (2004), 10.1103/PhysRevA.70.014501
  // Table I (Cs)
  // HF, RPA, BO, SR, Norm
  const std::map<std::string, std::array<double, 5>> JohnsonData = {
    {"6p-", {160.88, 40.66, 84.40, 5.43, -1.20}},
    {"6p+", {23.92, 18.84, 16.08, -7.51, -0.23}},
  };

  // Nuclear radius not stated in paper; neither was model for nuclear distro
  // Though it seems Ball was used; difference from pointlike very small
  const auto rN = std::sqrt(5.0 / 3) * 4.8041 / PhysConst::aB_fm;
  const auto F = DiracOperator::Hyperfine::sphericalBall_F(1);
  const auto gI = 0.7376; // not explicitely presented in paper

  const auto hfs = DiracOperator::hfs(1, gI, rN, wf.grid(), F);
  auto dV =
    ExternalField::DiagramRPA(&hfs, wf.basis(), wf.vHF(), "false", false);
  dV.solve_core(0.0, 100, false);

  // I assume they didn't include 'fitting' for BO, though they may have:
  MBPT::StructureRad srn(wf.basis(), wf.FermiLevel(), {3, 99}, "", {}, {},
                         false);
  srn.solve_core(&hfs, &dV);

  std::cout << "Table I [Phys. Rev. A 70, 014501 (2004)]\n";
  fmt::print("{:3} :  {:7}  {:6}  {:6}  {:5}  {:5}\n", "v", " HF", " RPA",
             " BO", " SR", " Norm");
  for (const auto &v : wf.valence()) {

    const auto key = v.shortSymbol();
    const auto data_ptr = JohnsonData.find(key);
    if (data_ptr == JohnsonData.end())
      continue;
    const auto [HF, RPA, BO, SR, NORM] = data_ptr->second;

    const auto rme_to_A = DiracOperator::Hyperfine::convert_RME_to_HFSconstant(
      1, v.kappa(), v.kappa());
    const auto hf = hfs.reducedME(v, v) * rme_to_A;
    const auto rpa = dV.dV(v, v) * rme_to_A;

    const auto bo = srn.BO(v, v) * rme_to_A;
    const auto sr = srn.SR(v, v) * rme_to_A;
    const auto norm = srn.norm(v, v, &hfs, &dV) * rme_to_A;

    fmt::print("{:3} :  {:+7.2f}  {:+6.2f}  {:+6.2f}  {:+5.2f}  {:+5.2f}\n",
               v.shortSymbol(), hf, rpa, bo, sr, norm);
    fmt::print("    [  {:+7.2f}  {:+6.2f}  {:+6.2f}  {:+5.2f}  {:+5.2f} ]\n",
               HF, RPA, BO, SR, NORM);

    CHECK(hf == Approx(HF).margin(0.05));
    CHECK(rpa == Approx(RPA).margin(0.5));
    CHECK(bo == Approx(BO).margin(1.5)); // sensitive to basis
    CHECK(sr == Approx(SR).margin(0.75));
    CHECK(norm == Approx(NORM).margin(0.01));
  }
}