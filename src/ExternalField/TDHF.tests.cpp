#include "DiracOperator/DiracOperator.hpp"
#include "IO/ChronoTimer.hpp"
#include "TDHF.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
// Tests for TDHF (basic unit)
TEST_CASE("External Field: TDHF - basic unit tests",
          "[ExternalField][TDHF][RPA][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: TDHF - basic unit tests\n";

  Wavefunction wf({800, 1.0e-4, 100.0, 20.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.solve_valence("6sp");

  auto dE1 = DiracOperator::E1(wf.grid());

  auto F6s = wf.getState("6s");
  auto F6p = wf.getState("6p-");
  REQUIRE(F6s != nullptr);
  REQUIRE(F6p != nullptr);

  auto rpa = ExternalField::TDHF(&dE1, wf.vHF());
  rpa.solve_core(0.0, 30);

  const auto dv1_0 = rpa.dV1(*F6s, *F6p);
  const auto dv_0 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv1_0 << " " << dv_0 << "\n";
  REQUIRE(dv1_0 != 0.0);
  REQUIRE(dv_0 != 0.0);
  REQUIRE(std::abs(dv_0) < std::abs(dv1_0));

  // Doing 1 iteration shuold be equiv to first-order dV
  auto rpa1 = ExternalField::TDHF(&dE1, wf.vHF());
  rpa1.solve_core(0.0, 1);
  const auto dv1_1 = rpa1.dV(*F6s, *F6p);
  REQUIRE(dv1_1 == Approx(dv1_0));

  const auto &Fc = wf.core().back();
  auto dPsis1 = rpa.solve_dPsis(Fc, 0.0, ExternalField::dPsiType::X);
  const auto &dPsis2 = rpa.get_dPsis(Fc, ExternalField::dPsiType::X);
  // should both be first-order:
  auto dPsis01 = rpa.get_dPsis_0(Fc, ExternalField::dPsiType::X);
  const auto &dPsis02 = rpa1.get_dPsis(Fc, ExternalField::dPsiType::X);
  REQUIRE(dPsis1.size() == dPsis2.size());
  REQUIRE(dPsis1.size() != 0);
  REQUIRE(dPsis01.size() == dPsis2.size());
  REQUIRE(dPsis01.size() != 0);
  for (std::size_t i = 0; i < dPsis1.size(); ++i) {
    // only true if TDHF converged!
    REQUIRE(dPsis1.at(i).norm2() ==
            Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
    const auto kappa = dPsis1.at(i).kappa();
    auto dPsi1 = rpa.solve_dPsi(Fc, 0.0, ExternalField::dPsiType::X, kappa);
    const auto &dPsi2 = rpa.get_dPsi_x(Fc, ExternalField::dPsiType::X, kappa);
    REQUIRE(dPsi1.norm2() == Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
    REQUIRE(dPsi2.norm2() == Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
    // only true if TDHF converged!
    REQUIRE(dPsis01.at(i).norm2() ==
            Approx(dPsis02.at(i).norm2()).epsilon(1.0e-2));
  }

  // should be zero after clearing:
  rpa.clear();
  const auto dv1_00 = rpa.dV1(*F6s, *F6p);
  const auto dv_00 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv1_00 << " " << dv_00 << "\n";
  REQUIRE(dv1_00 == 0.0);
  REQUIRE(dv_00 == 0.0);
}

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
TEST_CASE("External Field: TDHF (RPA)",
          "[ExternalField][TDHF][RPA][integration]") {
  {
    IO::ChronoTimer t{"TDHF"};
    std::cout << "-------------------------------------\n";
    std::cout << "External Field: TDHF (RPA)\n";

    // Create wavefunction object, solve HF for core+valence
    Wavefunction wf({6000, 1.0e-6, 175.0, 20.0, "loglinear", -1.0},
                    {"Cs", 133, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
    wf.solve_valence("7sp5d4f");

    // Lambda to compare mine to test data
    auto cmpr = [](const auto &ds, const auto &v) {
      return (ds.second - v) / v;
    };

    { // E1, w = 0.0
      const auto ww = 0.00;
      auto h = DiracOperator::E1(wf.grid());
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
      // pass &= qip::check_value(&obuff, "RPA E1 w=0 " + worst, eps,
      // 0.0, 5.0e-5);
      std::cout << "TDHF: RPA E1 w=0 " << worst << " " << eps << "\n";
      REQUIRE(std::abs(eps) < 5.0e-5);
    }

    { // E1, w = 0.05
      const auto ww = 0.05;
      auto h = DiracOperator::E1(wf.grid());
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
      // pass &=
      // qip::check_value(&obuff, "RPA E1 w=0.05 " + worst, eps, 0.0, 5.0e-5);
      std::cout << "TDHF: RPA E1 w=0.05 " << worst << " " << eps << "\n";
      REQUIRE(std::abs(eps) < 5.0e-5);
    }

    { // PNC, w = 0.0
      // Note: even zero-order PNC disagrees at ~5th digit - possibly due to c,t?
      const auto ww = 0.0;
      const auto c = Nuclear::c_hdr_formula_rrms_t(wf.get_rrms());
      auto h = DiracOperator::PNCnsi(c, Nuclear::default_t, wf.grid());
      std::vector<double> expected_VD = {
          1.5428e-04,  9.0137e-05,  7.9206e-05,  4.6296e-05,  -1.5428e-04,
          -7.9206e-05, 4.9201e-05,  -9.0137e-05, -4.6296e-05, 3.0130e-05,
          -4.9201e-05, -3.0130e-05, 3.3958e-06,  -3.3958e-06};
      std::sort(begin(expected_VD), end(expected_VD));

      const auto result = helper::dV_result(wf, h, ww);
      const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
      const std::string worst = at == result.end() ? "" : at->first;
      // pass &= qip::check_value(&obuff, "RPA PNC w=0 " + worst, eps,
      // 0.0, 5.0e-4);
      std::cout << "TDHF: RPA PNC w=0 " << worst << " " << eps << "\n";
      REQUIRE(std::abs(eps) < 5.0e-4);
    }
  }
}

//==============================================================================
inline std::vector<std::pair<std::string, double>>
helper::dV_result(const Wavefunction &wf,
                  const DiracOperator::TensorOperator &h, double ww) {

  // Form TDHF (RPA) object for this operator
  auto dV = ExternalField::TDHF(&h, wf.vHF());
  // Solve set of TDHF equations for core, with frequency ww
  const auto max_iterations = 150;
  const auto print_details = true;
  dV.solve_core(ww, max_iterations, print_details);

  std::vector<std::pair<std::string, double>> result;
  for (const auto &Fv : wf.valence()) {
    for (const auto &Fm : wf.valence()) {
      if (h.isZero(Fv.kappa(), Fm.kappa()))
        continue;

      result.emplace_back(Fv.shortSymbol() + Fm.shortSymbol(), dV.dV(Fv, Fm));
    }
  }
  std::sort(begin(result), end(result),
            [](auto a, auto b) { return a.second < b.second; });
  return result;
}
