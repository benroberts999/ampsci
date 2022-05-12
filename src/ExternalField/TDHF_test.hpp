#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include "TDHF.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

namespace UnitTest {

//==============================================================================
namespace helper {

// Calculates core-polarisation correction to matrix elements between all
// valence states, returns a vector of {states(string), value}
inline std::vector<std::pair<std::string, double>>
dV_result(const Wavefunction &wf, const DiracOperator::TensorOperator &h,
          double ww);

} // namespace helper

//==============================================================================
//==============================================================================
//! Unit tests External Field (RPA equations using TDHF method)
bool TDHF(std::ostream &obuff) {
  bool pass = true;

  // Create wavefunction object, solve HF for core+valence
  Wavefunction wf({6000, 1.0e-6, 175.0, 20.0, "loglinear", -1.0},
                  {"Cs", 133, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.solve_valence("7sp5d4f");

  // Lambda to compare mine to test data
  auto cmpr = [](const auto &ds, const auto &v) { return (ds.second - v) / v; };

  { // E1, w = 0.0
    const auto ww = 0.00;
    auto h = DiracOperator::E1(*wf.rgrid);
    // Expected <a||dV||b> values from V. Dzuba code:
    std::vector<double> expected_VD = {
        -0.303269160, 0.413337460,  -0.133010110, 0.185988360,  -0.036226022,
        0.041204439,  -0.087804876, 0.117345560,  -0.303269160, -0.036226022,
        -0.338428260, -0.413337460, -0.041204439, 0.145969420,  -0.433068210,
        -0.133010110, -0.087804876, -0.107910900, -0.185988360, -0.117345560,
        0.048037122,  -0.141564780, 0.338428260,  0.145969420,  0.107910900,
        0.048037122,  -0.138879730, 0.433068210,  0.141564780,  0.036916812,
        -0.165144850, 0.138879730,  0.036916812,  0.165144850};
    // -ve: |e|r vs er = -|e|r
    qip::scale(&expected_VD, -1.0);
    // sort to allow easy comparison:
    std::sort(begin(expected_VD), end(expected_VD));

    const auto result = helper::dV_result(wf, h, ww);
    const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
    const std::string worst = at == result.end() ? "" : at->first;
    pass &= qip::check_value(&obuff, "RPA E1 w=0 " + worst, eps, 0.0, 5.0e-5);
  }

  { // E1, w = 0.05
    const auto ww = 0.05;
    auto h = DiracOperator::E1(*wf.rgrid);
    std::vector<double> expected_VD = {
        -0.303213300, 0.412942330,  -0.132817610, 0.185545850,  -0.037186418,
        0.042653995,  -0.087825599, 0.117263410,  -0.303213300, -0.037186418,
        -0.342497480, -0.412942330, -0.042653995, 0.147791470,  -0.439517730,
        -0.132817610, -0.087825599, -0.107218100, -0.185545850, -0.117263410,
        0.047727162,  -0.139992440, 0.342497480,  0.147791470,  0.107218100,
        0.047727162,  -0.139387570, 0.439517730,  0.139992440,  0.037028834,
        -0.165663560, 0.139387570,  0.037028834,  0.165663560};
    // -ve: |e|r vs er = -er
    qip::scale(&expected_VD, -1.0);
    // sort to allow easy comparison:
    std::sort(begin(expected_VD), end(expected_VD));

    const auto result = helper::dV_result(wf, h, ww);
    const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
    const std::string worst = at == result.end() ? "" : at->first;
    pass &=
        qip::check_value(&obuff, "RPA E1 w=0.05 " + worst, eps, 0.0, 5.0e-5);
  }

  { // PNC, w = 0.0
    // Note: even zero-order PNC disagrees at ~5th digit - possibly due to c,t?
    const auto ww = 0.0;
    const auto c = Nuclear::c_hdr_formula_rrms_t(wf.get_rrms());
    auto h = DiracOperator::PNCnsi(c, Nuclear::default_t, *wf.rgrid);
    std::vector<double> expected_VD = {
        1.5428e-04,  9.0137e-05,  7.9206e-05,  4.6296e-05,  -1.5428e-04,
        -7.9206e-05, 4.9201e-05,  -9.0137e-05, -4.6296e-05, 3.0130e-05,
        -4.9201e-05, -3.0130e-05, 3.3958e-06,  -3.3958e-06};
    std::sort(begin(expected_VD), end(expected_VD));

    const auto result = helper::dV_result(wf, h, ww);
    const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
    const std::string worst = at == result.end() ? "" : at->first;
    pass &= qip::check_value(&obuff, "RPA PNC w=0 " + worst, eps, 0.0, 5.0e-4);
  }

  return pass;
}

} // namespace UnitTest

//==============================================================================
inline std::vector<std::pair<std::string, double>>
UnitTest::helper::dV_result(const Wavefunction &wf,
                            const DiracOperator::TensorOperator &h, double ww) {

  // Form TDHF (RPA) object for this operator
  auto dV = ExternalField::TDHF(&h, wf.getHF());
  // Solve set of TDHF equations for core, with frequency ww
  const auto max_iterations = 150;
  const auto print_details = true;
  dV.solve_core(ww, max_iterations, print_details);

  std::vector<std::pair<std::string, double>> result;
  for (const auto &Fv : wf.valence) {
    for (const auto &Fm : wf.valence) {
      if (h.isZero(Fv.kappa(), Fm.kappa()))
        continue;

      result.emplace_back(Fv.shortSymbol() + Fm.shortSymbol(), dV.dV(Fv, Fm));
    }
  }
  std::sort(begin(result), end(result),
            [](auto a, auto b) { return a.second < b.second; });
  return result;
}
