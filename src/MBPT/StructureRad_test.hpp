#pragma once
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
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
  const auto expected =
      std::vector<sp>{{"3p-", "3s+", 0.0030 / 3.6906, -0.0050 / 3.6906},
                      {"3p+", "3s+", 0.0043 / 5.2188, -0.0070 / 5.2188},
                      {"3p-", "4s+", -0.0010 / 3.6004, -0.0021 / 3.6004},
                      {"3p+", "4s+", -0.0014 / 5.1012, -0.0029 / 5.1012}};

  std::cout << "Structure Radiation + Norm of states;\ncf Johnson et al, "
               "At.Dat.Nuc.Dat.Tables 64, 279 (1996), Table E (Na)\n";
  double worst_sr = 0.0;
  double worst_ns = 0.0;
  std::string at_sr = "";
  std::string at_ns = "";
  for (const auto [w_str, v_str, sr_exp, n_exp] : expected) {

    // find the right basis states for SR/N "legs"
    const auto ws =
        std::find_if(cbegin(wf.basis), cend(wf.basis),
                     [w_str](auto &a) { return a.shortSymbol() == w_str; });
    const auto vs =
        std::find_if(cbegin(wf.basis), cend(wf.basis),
                     [v_str](auto &a) { return a.shortSymbol() == v_str; });
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
                           0.12);
  pass &= qip::check_value(&obuff, "NormStates(E1,Na) " + at_ns, worst_ns, 0.0,
                           0.11);

  return pass;
}
} // namespace UnitTest
