
#include "TDHFbasis.hpp"
#include "DiagramRPA.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
// Tests for TDHFbasis (basic unit)
TEST_CASE("External Field: TDHFbasis - basic unit tests",
          "[ExternalField][TDHFbasis][RPA][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: TDHFbasis - basic unit tests\n";

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]");
  wf.solve_valence("3sp");
  // nb: use very small basis. Don't care about numerical results, just that eveything is working correctly.
  wf.formBasis({"10sp8d", 20, 7, 1.0e-3, 1.0e-3, 30.0, false});

  auto dE1 = DiracOperator::E1(wf.grid());

  auto F6s = wf.getState("3s");
  auto F6p = wf.getState("3p-");
  REQUIRE(F6s != nullptr);
  REQUIRE(F6p != nullptr);

  auto rpa = ExternalField::TDHFbasis(&dE1, wf.vHF(), wf.basis());
  rpa.solve_core(0.0, 15);

  const auto dv_0 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv_0 << "\n";
  REQUIRE(dv_0 != 0.0);

  const auto &Fc = wf.core().back();

  auto dPsis1 = rpa.form_dPsis(Fc, 0.0, ExternalField::dPsiType::X, wf.basis());
  auto dPsis2 = rpa.get_dPsis(Fc, ExternalField::dPsiType::X);
  REQUIRE(dPsis1.size() == dPsis2.size());
  REQUIRE(dPsis1.size() != 0);
  for (std::size_t i = 0; i < dPsis1.size(); ++i) {
    REQUIRE(dPsis1.at(i).norm2() ==
            Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
  }

  rpa.clear();
  const auto dv_00 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv_00 << "\n";
  REQUIRE(dv_00 == 0.0);
}

//==============================================================================
namespace helper {

// Calculates core-polarisation correction to matrix elements between all
// valence states, returns a vector of {states(string), value}
inline std::vector<std::pair<std::string, double>>
dV_result_basis(const Wavefunction &wf, const DiracOperator::TensorOperator &h,
                double ww);

} // namespace helper

//==============================================================================
//==============================================================================
//! Unit tests External Field (RPA equations using TDHF method)
TEST_CASE("External Field: TDHFbasis (RPA)",
          "[ExternalField][TDHFbasis][RPA][slow][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: TDHFbasis (RPA)\n";

  // Create wavefunction object, solve HF for core+valence
  Wavefunction wf({3500, 1.0e-6, 150.0, 50.0, "loglinear", -1.0},
                  {"Cs", 133, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.solve_valence("7sp5d4f");

  SplineBasis::Parameters bspl_param;
  {
    bspl_param.states = "50spd40f20g";
    bspl_param.n = 60;
    bspl_param.k = 7;
    bspl_param.r0 = 1.0e-5;
    bspl_param.reps = 0.0;
    bspl_param.rmax = 50.0;
  }
  wf.formBasis(bspl_param);

  // Lambda to compare mine to test data
  auto cmpr = [](const auto &ds, const auto &v) { return (ds.second - v) / v; };

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

    const auto result = helper::dV_result_basis(wf, h, ww);
    const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
    const std::string worst = at == result.end() ? "" : at->first;
    // pass &= qip::check_value(&obuff, "TDHF(basis) E1 w=0 " + worst, eps, 0.0,
    //  5.0e-5);
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

    const auto result = helper::dV_result_basis(wf, h, ww);
    const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
    const std::string worst = at == result.end() ? "" : at->first;
    // pass &= qip::check_value(&obuff, "TDHF(basis) E1 w=0.05 " + worst, eps,
    // 0.0, 5.0e-5);
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

    const auto result = helper::dV_result_basis(wf, h, ww);
    const auto [eps, at] = qip::compare(result, expected_VD, cmpr);
    const std::string worst = at == result.end() ? "" : at->first;
    // pass &= qip::check_value(&obuff, "TDHF(basis) PNC w=0 " + worst, eps,
    // 0.0, 5.0e-3);
    REQUIRE(std::abs(eps) < 5.0e-3);
  }
}

//==============================================================================
inline std::vector<std::pair<std::string, double>>
helper::dV_result_basis(const Wavefunction &wf,
                        const DiracOperator::TensorOperator &h, double ww) {

  // Form TDHF (RPA) object for this operator
  auto dV = ExternalField::TDHFbasis(&h, wf.vHF(), wf.basis());
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
