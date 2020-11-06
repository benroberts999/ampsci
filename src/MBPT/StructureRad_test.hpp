#pragma once
#include "DiracOperator/Operators.hpp"
#include "ExternalField/TDHF.hpp"
#include "MBPT/StructureRad.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include <algorithm>
#include <string>

namespace UnitTest {

//******************************************************************************
//! Unit tests for StructureRad + Normalisation of states
bool StructureRad(std::ostream &obuff) {
  bool pass = true;

  { // Test for Na: (use splines for legs)
    Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                    {"Na", -1, "Fermi"}, 1.0);
    wf.hartreeFockCore("HartreeFock", 0.0, "[Ne]");
    wf.hartreeFockValence("4s3p");
    wf.formBasis({"20spdfgh", 30, 9, 1.0e-4, 1.0e-6, 60.0, false});

    // Find core/valence energy: allows distingush core/valence states
    const auto en_core = wf.en_coreval();

    const auto h = DiracOperator::E1(*wf.rgrid);

    MBPT::StructureRad sr(wf.basis, en_core);

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
          cbegin(wf.basis), cend(wf.basis),
          [&sym = w_str](auto &a) { return a.shortSymbol() == sym; });
      const auto vs = std::find_if(
          cbegin(wf.basis), cend(wf.basis),
          [&sym = v_str](auto &a) { return a.shortSymbol() == sym; });
      assert(ws != cend(wf.basis) && vs != cend(wf.basis));

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

    // Data only known to 2 digits (with leading 1); so best is ~10%
    pass &= qip::check_value(&obuff, "StructRad(E1,Na) " + at_sr, worst_sr, 0.0,
                             0.1);
    pass &= qip::check_value(&obuff, "NormStates(E1,Na) " + at_ns, worst_ns,
                             0.0, 0.05);
  }

  //****************************************************************************
  { // Test for Cs: use valence states for legs
    Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                    {"Cs", -1, "Fermi"}, 1.0);
    wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
    wf.hartreeFockValence("7s6p");
    wf.formBasis({"20spdfgh", 30, 9, 1.0e-4, 1.0e-6, 60.0, false});
    // Note: We use a very small basis, so the test can run in reasonable time
    // However, we get pretty good comparison to Johnson, so this is fine!

    // Find core/valence energy: allows distingush core/valence states
    const auto en_core = wf.en_coreval();

    const auto h = DiracOperator::E1(*wf.rgrid);

    // Only include core states above+including n=3
    MBPT::StructureRad sr(wf.basis, en_core, {3, 20});

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
          cbegin(wf.valence), cend(wf.valence),
          [&sym = w_str](auto &a) { return a.shortSymbol() == sym; });
      const auto vs = std::find_if(
          cbegin(wf.valence), cend(wf.valence),
          [&sym = v_str](auto &a) { return a.shortSymbol() == sym; });
      assert(ws != cend(wf.valence) && vs != cend(wf.valence));

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
    pass &= qip::check_value(&obuff, "StructRad(E1,Cs) " + at_sr, worst_sr, 0.0,
                             0.05);
    pass &= qip::check_value(&obuff, "NormStates(E1,Cs) " + at_ns, worst_ns,
                             0.0, 0.01);
  }

  return pass;
}
} // namespace UnitTest
