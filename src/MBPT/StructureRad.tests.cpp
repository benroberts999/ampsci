#include "StructureRad.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"
#include <algorithm>
#include <string>

//==============================================================================
TEST_CASE("MBPT: Structure Rad + Norm, basic", "[StrucRad][MBPT][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: Structure Rad + Norm: basic unit\n";

  // note: does not test formulas: just checks class is working correctly.
  //  Other (integration) tests below check correctness/accuracy of formulas

  Wavefunction wf({200, 1.0e-2, 30.0, 10.0, "loglinear", -1.0},
                  {"Li", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", 0.0, "[He]");
  wf.formBasis({"5spdf", 20, 5, 1.0e-2, 1.0e-2, 20.0});

  const auto h = DiracOperator::E1(wf.grid());

  MBPT::StructureRad srn(wf.basis(), wf.FermiLevel());

  // test srn_table
  const auto tab = srn.srn_table(&h, srn.core(), srn.excited());
  for (const auto &a : srn.core()) {
    for (const auto &b : srn.excited()) {
      if (h.isZero(a, b))
        continue;
      const auto tab_ab = tab.get(a, b);
      const auto tab_ba = tab.get(b, a);
      REQUIRE(tab_ab != nullptr);
      REQUIRE(tab_ba != nullptr);
      const auto [srab, xab] = *tab_ab;
      const auto [srba, xba] = *tab_ba;
      REQUIRE(std::abs(srab - h.symm_sign(a, b) * srba) < 1.0e-10);

      const auto srn0 = srn.srTB(&h, a, b).first + srn.srC(&h, a, b).first +
                        srn.norm(&h, a, b).first;
      REQUIRE(std::abs(srab - srn0) < 1.0e-10);
    }
  }

  // TEST SRN read-write using Qk_file
  const auto rand_str = qip::random_string(6);
  const auto fname = "deleteme_" + rand_str + ".qk.abf";
  MBPT::StructureRad srn2(wf.basis(), wf.FermiLevel(), {0, 999}, fname);
  MBPT::StructureRad srn3(wf.basis(), wf.FermiLevel(), {0, 999}, fname);
  for (const auto &a : wf.basis()) {
    for (const auto &b : wf.basis()) {
      if (a < b && !h.isZero(a, b)) {
        const auto srC0 = srn.srC(&h, a, b).first;
        const auto srC2 = srn2.srC(&h, a, b).first;
        const auto srC3 = srn3.srC(&h, a, b).first;
        REQUIRE(std::abs(srC0 - srC2) < 1.0e-10);
        REQUIRE(std::abs(srC0 - srC3) < 1.0e-10);

        const auto srTB2 = srn2.srTB(&h, a, b).first;
        const auto srTB3 = srn3.srTB(&h, a, b).first;
        REQUIRE(std::abs(srTB2 - srTB3) < 1.0e-14);

        const auto srN2 = srn2.norm(&h, a, b).first;
        const auto srN3 = srn3.norm(&h, a, b).first;
        REQUIRE(std::abs(srN2 - srN3) < 1.0e-14);
      }
    }
  }
}

//==============================================================================
TEST_CASE("MBPT: Structure Rad + Norm",
          "[StrucRad][MBPT][ExternalField][integration][slow]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "MBPT: Structure Rad + Norm\n";

  { // Test for Na: (use splines for legs)
    Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                    {"Na", -1, "Fermi"}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Ne]");
    wf.solve_valence("4s3p");
    wf.formBasis({"20spdfgh", 30, 9, 1.0e-4, 1.0e-6, 60.0});

    // Find core/valence energy: allows distingush core/valence states
    const auto en_core = wf.FermiLevel();

    const auto h = DiracOperator::E1(wf.grid());

    int nmin = 1;
    int nmax = 99;
    MBPT::StructureRad sr(wf.basis(), en_core, {nmin, nmax});

    // Expected data, from: Johnson et al, At.Dat.Nuc.Dat.Tables 64, 279 (1996),
    // Table E (Na)
    // Data in form: {a, b, SR/t0, NS/t0}
    using sp = std::tuple<std::string, std::string, double, double>;
    const auto expected = std::vector<sp>{
      {"3p-", "3s+", 0.0030 / 3.6906, -0.0050 / 3.6906},
      {"3p+", "3s+", 0.0043 / 5.2188, -0.0070 / 5.2188},
      //{"3p-", "4s+", -0.0010 / 3.6004, -0.0021 / 3.6004},
      {"3p+", "4s+", -0.0014 / 5.1012, -0.0029 / 5.1012}
      // Skip this just in interest of time
    };

    std::cout << "Structure Radiation + Norm of states;\ncf Johnson et al, "
                 "At.Dat.Nuc.Dat.Tables 64, 279 (1996), Table E (Na)\n";
    double worst_sr = 0.0;
    double worst_ns = 0.0;
    std::string at_sr = "";
    std::string at_ns = "";
    for (const auto &[w_str, v_str, sr_exp, n_exp] : expected) {

      // find the right basis states for SR/N "legs"
      const auto ws = std::find_if(
        cbegin(wf.basis()), cend(wf.basis()),
        [&sym = w_str](auto &a) { return a.shortSymbol() == sym; });
      const auto vs = std::find_if(
        cbegin(wf.basis()), cend(wf.basis()),
        [&sym = v_str](auto &a) { return a.shortSymbol() == sym; });
      assert(ws != cend(wf.basis()) && vs != cend(wf.basis()));

      // My calculations:
      const auto t0 = h.reducedME(*ws, *vs); // splines here?
      const auto [tb, dv1] = sr.srTB(&h, *ws, *vs, 0.0);
      const auto [c, dv2] = sr.srC(&h, *ws, *vs);
      const auto [n, dv3] = sr.norm(&h, *ws, *vs);

      // output table as we go:
      std::cout << "<" << w_str << "||" << v_str << ">: ";
      printf("%8.1e  [%8.1e] | ", (tb + c) / t0, sr_exp);
      printf("%8.1e  [%8.1e]\n", n / t0, n_exp);

      // Compare relative shifts to Johnson, store worst
      const auto eps_sr = std::abs(((tb + c) / t0 - sr_exp) / sr_exp);
      const auto eps_ns = std::abs((n / t0 - n_exp) / n_exp);
      if (eps_sr > worst_sr) {
        worst_sr = eps_sr;
        at_sr = w_str + "-" + v_str;
      }
      if (eps_ns > worst_ns) {
        worst_ns = eps_ns;
        at_ns = w_str + "-" + v_str;
      }
    }

    REQUIRE(std::abs(worst_sr) < 0.1);
    REQUIRE(std::abs(worst_ns) < 0.05);
  }

  //============================================================================
  { // Test for Cs: use valence states for legs
    Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                    {"Cs", -1, "Fermi"}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7s6p");
    wf.formBasis({"20spdfgh", 30, 9, 1.0e-4, 1.0e-6, 60.0});
    // Note: We use a very small basis, so the test can run in reasonable time
    // However, we get pretty good comparison to Johnson, so this is fine!

    // Find core/valence energy: allows distingush core/valence states
    const auto en_core = wf.FermiLevel();

    const auto h = DiracOperator::E1(wf.grid());

    // Only include core states above+including n=3
    MBPT::StructureRad sr(wf.basis(), en_core, {3, 99});

    // Expected data, from: Johnson et al, At.Dat.Nuc.Dat.Tables 64, 279 (1996),
    // Table J (Cs)
    // Data in form: {a, b, SR/t0, NS/t0}
    using sp = std::tuple<std::string, std::string, double, double>;
    const auto expected = std::vector<sp>{
      {"6p-", "6s+", 0.0445 / 5.2777, -0.0508 / 5.2777},
      {"6p+", "6s+", 0.0593 / 7.4265, -0.0694 / 7.4265},
      //{"6p-", "7s+", -0.0120 / 4.4135, -0.0198 / 4.4135},
      {"6p+", "7s+", -0.0153 / 6.6716, -0.0281 / 6.6716}
      // Skip this one just in interest of time
    };

    std::cout << "Structure Radiation + Norm of states;\ncf Johnson et al, "
                 "At.Dat.Nuc.Dat.Tables 64, 279 (1996), Table J (Cs)\n";
    double worst_sr = 0.0;
    double worst_ns = 0.0;
    std::string at_sr = "";
    std::string at_ns = "";
    for (const auto &[w_str, v_str, sr_exp, n_exp] : expected) {

      // this time, use valence states:
      const auto ws = std::find_if(
        cbegin(wf.valence()), cend(wf.valence()),
        [&sym = w_str](auto &a) { return a.shortSymbol() == sym; });
      const auto vs = std::find_if(
        cbegin(wf.valence()), cend(wf.valence()),
        [&sym = v_str](auto &a) { return a.shortSymbol() == sym; });
      assert(ws != cend(wf.valence()) && vs != cend(wf.valence()));

      // My calculations:
      const auto t0 = h.reducedME(*ws, *vs); // splines here?
      const auto [tb, dv1] = sr.srTB(&h, *ws, *vs, 0.0);
      const auto [c, dv2] = sr.srC(&h, *ws, *vs);
      const auto [n, dv3] = sr.norm(&h, *ws, *vs);

      // output table as we go:
      std::cout << "<" << w_str << "||" << v_str << ">: ";
      printf("%9.2e  [%9.2e] | ", (tb + c) / t0, sr_exp);
      printf("%9.2e  [%9.2e]\n", n / t0, n_exp);

      // Compare relative shifts to Johnson, store worst
      const auto eps_sr = std::abs(((tb + c) / t0 - sr_exp) / sr_exp);
      const auto eps_ns = std::abs((n / t0 - n_exp) / n_exp);
      if (eps_sr > worst_sr) {
        worst_sr = eps_sr;
        at_sr = w_str + "-" + v_str;
      }
      if (eps_ns > worst_ns) {
        worst_ns = eps_ns;
        at_ns = w_str + "-" + v_str;
      }
    }

    // Aim for better than 5% for SR, and 1% for Norm
    // Note: Quite possible ours are more accurate, despite very small basis
    REQUIRE(std::abs(worst_sr) < 0.05);
    REQUIRE(std::abs(worst_ns) < 0.01);
  }
}
