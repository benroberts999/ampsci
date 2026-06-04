#include "TDHF.hpp"
#include "DiracOperator/include.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "calcMatrixElements.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
TEST_CASE("External Field: TDHF basic unit tests",
          "[ExternalField][TDHF][RPA][unit]") {

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);
  wf.solve_valence("3sp", false);
  // nb: use very small basis.
  // Don't care about numerical results, just that eveything is working correctly.
  SplineBasis::Parameters params{"10sp8d", 20, 7, 1.0e-3, 1.0e-3, 30.0};
  params.verbose = false;
  wf.formBasis(params);

  auto dE1 = DiracOperator::E1(wf.grid());

  auto rpa = ExternalField::TDHF(&dE1, wf.vHF());
  auto rpa2 = ExternalField::make_rpa("TDHF", &dE1, wf.vHF(), true, {}, "");

  rpa.solve_core(0.0, 25);
  rpa2->solve_core(0.0, 25);

  for (const auto &v : wf.valence()) {
    for (const auto &w : wf.valence()) {

      const auto dv_0 = rpa.dV(v, w);
      const auto dv_1 = rpa.dV(w, v) * dE1.symm_sign(v, w);

      const auto dv_2 = rpa2->dV(v, w);

      if (!dE1.isZero(v, w))
        REQUIRE(dv_0 != 0.0);
      REQUIRE(dv_0 == Approx(dv_1));
      REQUIRE(dv_0 == Approx(dv_2));
    }
  }

  for (const auto &Fc : wf.core()) {
    auto dPsis1 = rpa.solve_dPsis(Fc, 0.0, ExternalField::dPsiType::X);
    auto dPsis2 = rpa.get_dPsis(Fc, ExternalField::dPsiType::X);
    REQUIRE(dPsis1.size() == dPsis2.size());
    REQUIRE(dPsis1.size() != 0);
    for (std::size_t i = 0; i < dPsis1.size(); ++i) {
      // only perfectly true if converged perfectly
      REQUIRE(dPsis1.at(i).norm2() ==
              Approx(dPsis2.at(i).norm2()).epsilon(1.0e-4));

      const auto kappa = dPsis1.at(i).kappa();
      auto dPsi1 = rpa.solve_dPsi(Fc, 0.0, ExternalField::dPsiType::X, kappa);
      const auto &dPsi2 = rpa.get_dPsi_x(Fc, ExternalField::dPsiType::X, kappa);
      REQUIRE(dPsi1.norm2() == Approx(dPsis2.at(i).norm2()).epsilon(1.0e-4));
      REQUIRE(dPsi2.norm2() == Approx(dPsis2.at(i).norm2()).epsilon(1.0e-4));
    }
  }

  auto F6s = wf.getState("3s");
  auto F6p = wf.getState("3p-");
  REQUIRE(F6s != nullptr);
  REQUIRE(F6p != nullptr);

  // dV_rhs(Fa) should return a DiracSpinor with the requested kappa
  for (int kappa : {-1, 1, -2, 2, -3}) {
    REQUIRE(rpa.dV_rhs(kappa, *F6s).kappa() == kappa);
    REQUIRE(rpa.dV_rhs(kappa, *F6p).kappa() == kappa);
  }

  // Fb * dV_rhs(kappa_b, Fa) should equal dV(Fb, Fa)
  {
    const auto dv_ps = rpa.dV(*F6p, *F6s);
    const auto dv_sp = rpa.dV(*F6s, *F6p);
    const auto rhs_ps = *F6p * rpa.dV_rhs(F6p->kappa(), *F6s);
    const auto rhs_sp = *F6s * rpa.dV_rhs(F6s->kappa(), *F6p);
    REQUIRE(rhs_ps == Approx(dv_ps));
    REQUIRE(rhs_sp == Approx(dv_sp));
  }

  rpa.clear();
  const auto dv_00 = rpa.dV(*F6s, *F6p);
  REQUIRE(dv_00 == 0.0);
}

//==============================================================================
struct TestData {
  std::string a, b;
  double h, rpa;
};

// All data with: Omega = 0.0 or 0.1

const std::vector<TestData> dzuba_e1_0{
  {"6p-", "6s+", 5.277685, 4.974416},   {"6p-", "7s+", -4.413137, -4.449363},
  {"6p+", "6s+", 7.426433, 7.013095},   {"6p+", "7s+", -6.671013, -6.712218},
  {"7p-", "6s+", 0.3717475, 0.2387373}, {"7p-", "7s+", 11.00887, 10.92106},
  {"7p+", "6s+", 0.6947538, 0.5087653}, {"7p+", "7s+", 15.34478, 15.22744},
  {"5d-", "6p-", -8.978332, -8.639903}, {"5d-", "6p+", -4.062459, -3.916489},
  {"5d-", "7p-", 4.039455, 4.147366},   {"5d-", "7p+", 1.688036, 1.736073},
  {"5d+", "6p+", -12.18643, -11.75336}, {"5d+", "7p+", 5.024624, 5.166189},
  {"4f-", "5d-", -10.65973, -10.52085}, {"4f-", "5d+", -2.840239, -2.803322},
  {"4f+", "5d+", -12.70344, -12.5383}};

const std::vector<TestData> dzuba_e1_w{
  {"6p-", "6s+", 5.277685, 4.971342},   {"6p-", "7s+", -4.413137, -4.451618},
  {"6p+", "6s+", 7.426433, 7.009477},   {"6p+", "7s+", -6.671013, -6.715528},
  {"7p-", "6s+", 0.3717475, 0.2377626}, {"7p-", "7s+", 11.00887, 10.92009},
  {"7p+", "6s+", 0.6947538, 0.5077115}, {"7p+", "7s+", 15.34478, 15.22634},
  {"5d-", "6p-", -8.978332, -8.627918}, {"5d-", "6p+", -4.062459, -3.911142},
  {"5d-", "7p-", 4.039455, 4.147084},   {"5d-", "7p+", 1.688036, 1.735976},
  {"5d+", "6p+", -12.18643, -11.73544}, {"5d+", "7p+", 5.024624, 5.164451},
  {"4f-", "5d-", -10.65973, -10.51826}, {"4f-", "5d+", -2.840239, -2.802684},
  {"4f+", "5d+", -12.70344, -12.53541}};

const std::vector<TestData> dzuba_e2_0{
  {"6p+", "6p-", 68.48282, 68.23295},   {"6p+", "6p+", -70.2262, -69.97999},
  {"7p-", "6p+", 49.3854, 49.48805},    {"7p+", "6p-", -42.52614, -42.63387},
  {"7p+", "6p+", 46.81519, 46.9199},    {"7p+", "7p-", 299.8769, 299.8039},
  {"7p+", "7p+", -305.114, -305.0414},  {"5d-", "6s+", -43.84649, -43.63158},
  {"5d-", "7s+", 80.74682, 80.78614},   {"5d-", "5d-", -67.06309, -66.83697},
  {"5d+", "6s+", -53.71202, -53.47216}, {"5d+", "7s+", 98.51161, 98.55041},
  {"5d+", "5d-", 43.84507, 43.70795},   {"5d+", "5d+", -87.57531, -87.30408}};

const std::vector<TestData> dzuba_e2_w{
  {"6p+", "6p-", 68.48282, 68.23044},   {"6p+", "6p+", -70.2262, -69.97613},
  {"7p-", "6p+", 49.3854, 49.49049},    {"7p+", "6p-", -42.52614, -42.63489},
  {"7p+", "6p+", 46.81519, 46.92167},   {"7p+", "7p-", 299.8769, 299.8032},
  {"7p+", "7p+", -305.114, -305.0403},  {"5d-", "6s+", -43.84649, -43.62903},
  {"5d-", "7s+", 80.74682, 80.7866},    {"5d-", "5d-", -67.06309, -66.83467},
  {"5d+", "6s+", -53.71202, -53.46243}, {"5d+", "7s+", 98.51161, 98.54769},
  {"5d+", "5d-", 43.84507, 43.70911},   {"5d+", "5d+", -87.57531, -87.30183}};

const std::vector<TestData> dzuba_m1_0{{"6p+", "6p-", 1.153521, 1.153505},
                                       {"7p-", "6p+", 0.03014502, 0.03016005},
                                       {"7p+", "6p-", 0.02855965, 0.02855698},
                                       {"7p+", "7p-", 1.153342, 1.153337},
                                       {"5d+", "5d-", 1.549161, 1.54951},
                                       {"4f+", "4f-", 1.851639, 1.85164}};

const std::vector<TestData> dzuba_m1_w{
  {"6p+", "6p-", 1.153511, 1.153558},    {"7p-", "6p+", 0.03013823, 0.0301901},
  {"7p+", "6p-", 0.0285655, 0.02859815}, {"7p+", "7p-", 1.153301, 1.153316},
  {"5d+", "5d-", 1.549149, 1.549712},    {"4f+", "4f-", 1.851568, 1.85157}};

const std::vector<TestData> dzuba_pnc{
  {"6p-", "6s+", -5.72823E-4, -7.270836e-4},
  {"6p-", "7s+", -3.002691e-4, -3.794671e-4},
  {"7p-", "6s+", -3.429023e-4, -4.330309e-4},
  {"7p-", "7s+", -1.797466e-4, -2.260384e-4}};

// with g = 1.0;
const std::vector<TestData> dzuba_hfs{{"6s+", "6s+", 1943.394, 2342.449},
                                      {"7s+", "6s+", 1018.711, 1227.127},
                                      {"7s+", "7s+", 533.9993, 642.53},
                                      {"6p-", "6p-", 218.2665, 273.257},
                                      {"6p+", "6p+", 32.41904, 58.06549},
                                      {"7p-", "6p-", 130.5002, 162.7951},
                                      {"7p-", "7p-", 78.14994, 97.1154},
                                      {"7p+", "6p+", 19.46554, 34.75759},
                                      {"7p+", "7p+", 11.71131, 20.81699},
                                      {"5d-", "5d-", 24.71128, 21.90427},
                                      {"5d+", "5d+", 10.1202, -33.08472},
                                      {"4f-", "4f-", 0.05106033, 0.04225674},
                                      {"4f+", "4f+", 0.02838632, -0.01949691}};

//-/-/----------------------------------------------------------------/////-----
void test_RPA(const Wavefunction &wf, DiracOperator::TensorOperator &h,
              double omega, const std::vector<TestData> &test_data,
              double target_1, double target_2) {
  std::cout << "\n";
  ExternalField::TDHF dV(&h, wf.vHF());
  if (h.freqDependantQ()) {
    h.updateFrequency(omega);
  }
  dV.solve_core(omega);

  fmt::print("{:9s} {:>12s} [{:>12s}] {:5s}  {:>10s} [{:>10s}] "
             "{:5s}\n",
             h.name(), "<h+dV>", "expected", "eps", "<dV>", "expected", "eps");
  for (const auto &[a, b, h0, rpa0] : test_data) {
    const auto &Fa = *wf.getState(a);
    const auto &Fb = *wf.getState(b);

    const auto f = h.name() == "hfs1" ?
                     DiracOperator::Hyperfine::convert_RME_to_HFSconstant(
                       1, Fa.kappa(), Fb.kappa()) :
                     1.0;
    const auto hab = f * h.reducedME(Fa, Fb);
    const auto dv = f * dV.dV(Fa, Fb);

    // don't worry about sign (tested elsewhere)
    const auto s1 = hab / std::abs(hab);
    const auto s2 = h0 / std::abs(h0);

    const auto delta1 = (s1 * (hab + dv) - s2 * (rpa0));
    const auto eps1 = std::abs(delta1 / rpa0);

    const auto delta2 = (s1 * dv - s2 * (rpa0 - h0));
    const auto eps2 = std::abs(delta2 / (rpa0 - h0));

    fmt::print("{:4s} {:4s} {:12.5e} [{:12.5e}] {:.0e}  {:10.3e} [{:10.3e}] "
               "{:.0e}\n",
               a, b, s1 * (hab + dv), s2 * rpa0, eps1, s1 * dv,
               s2 * (rpa0 - h0), eps2);

    REQUIRE(eps1 < target_1);
    REQUIRE(eps2 < target_2);
  }
}

void test_RPA2(const Wavefunction &wf, DiracOperator::TensorOperator &h,
               double omega, const std::vector<TestData> &test_data,
               double target_1, double target_2) {
  ExternalField::TDHF dV(&h, wf.vHF());
  if (h.freqDependantQ()) {
    h.updateFrequency(omega);
  }
  dV.set_eta(0.75);
  dV.solve_core(omega);

  fmt::print("{:4s} {:>10s} [{:>10s}] {:5s}  {:>10s} [{:>10s}] "
             "{:5s}\n",
             wf.atomicSymbol(), "A0 (MHz)", "expected", "eps", "RPA (MHz)",
             "expected", "eps");
  for (const auto &[a, b, h0, rpa0] : test_data) {
    const auto &Fa = *wf.getState(a);
    const auto &Fb = *wf.getState(b);

    const auto f = DiracOperator::Hyperfine::convert_RME_to_HFSconstant(
      1, Fa.kappa(), Fb.kappa());
    const auto hab = f * h.reducedME(Fa, Fb);
    const auto dv = f * dV.dV(Fa, Fb);

    const auto delta1 = (hab - h0);
    const auto eps1 = std::abs(delta1 / h0);

    const auto delta2 = (hab + dv - rpa0);
    const auto eps2 = std::abs(delta2 / rpa0);

    fmt::print("{:4s} {:10.3e} [{:10.3e}] {:.0e}  {:10.3e} [{:10.3e}] "
               "{:.0e}\n",
               a, hab, h0, eps1, hab + dv, rpa0, eps2);

    REQUIRE(eps1 < target_1);
    REQUIRE(eps2 < target_2);
  }
}

//==============================================================================
//==============================================================================

//! Unit tests External Field (RPA equations using TDHF method)
TEST_CASE("External Field: TDHF (RPA)",
          "[ExternalField][TDHF][RPA][HFS][integration]") {
  {
    IO::ChronoTimer t{"TDHF"};

    // Create wavefunction object, solve HF for core+valence
    Wavefunction wf({6000, 1.0e-6, 175.0, 20.0, "loglinear", -1.0},
                    {"Cs", 133, "Fermi", 4.8041, 2.3}, 1.0);
    wf.solve_core("HartreeFock", "[Xe]");
    wf.solve_valence("7sp5d4f");

    DiracOperator::E1 dE1(wf.grid());
    DiracOperator::Ek dE2(wf.grid(), 2);
    DiracOperator::M1 m1(wf.grid(), wf.alpha());
    DiracOperator::PNCnsi dpnc(5.67073, 2.3, wf.grid());
    DiracOperator::hfs hfs(1, 1.0, 0.0, wf.grid(),
                           DiracOperator::Hyperfine::pointlike_F());

    test_RPA(wf, dE1, 0.0, dzuba_e1_0, 1.0e-4, 1.0e-4);

    test_RPA(wf, dE1, 0.1, dzuba_e1_w, 1.0e-4, 1.0e-4);

    // RPA is very small for E2, so test less stringent (simply, number of digits!)
    test_RPA(wf, dE2, 0.0, dzuba_e2_0, 1.0e-5, 1.0e-3);
    test_RPA(wf, dE2, 0.1, dzuba_e2_w, 1.0e-5, 1.0e-3);

    // RPA is extremely small for M1, so test less stringent (simply, number of digits!)
    test_RPA(wf, m1, 0.0, dzuba_m1_0, 1.0e-4, 3.0e-1);
    test_RPA(wf, m1, 0.1, dzuba_m1_w, 1.0e-4, 3.0e-1);

    test_RPA(wf, hfs, 0.0, dzuba_hfs, 1.0e-3, 3.0e-3);

    test_RPA(wf, dpnc, 0.0, dzuba_pnc, 1.0e-3, 1.0e-3);
  }
}

//==============================================================================
//==============================================================================

// Test data, from: 10.1103/PhysRevA.100.042506

const std::vector<TestData> pra_rb{
  {"5s+", "5s+", 2183, 2644},     {"6s+", "6s+", 583.1, 704.7},
  {"7s+", "7s+", 238.8, 288.5},   {"8s+", "8s+", 120.6, 145.6},
  {"9s+", "9s+", 69.24, 83.59},   {"10s+", "10s+", 43.37, 52.35},
  {"11s+", "11s+", 28.95, 34.94}, {"12s+", "12s+", 20.27, 24.47},
  {"5p-", "5p-", 236.8, 299.6},   {"6p-", "6p-", 83.21, 104.6},
  {"7p-", "7p-", 38.60, 48.41},   {"8p-", "8p-", 20.97, 26.27},
  {"9p-", "9p-", 12.64, 15.82},   {"10p-", "10p-", 8.192, 10.25},
  {"11p-", "11p-", 5.610, 7.018}, {"12p-", "12p-", 4.008, 5.013}};

const std::vector<TestData> pra_cs{
  {"6s+", "6s+", 1434, 1728},     {"7s+", "7s+", 393.9, 474.0},
  {"8s+", "8s+", 164.5, 197.8},   {"9s+", "9s+", 84.11, 101.1},
  {"10s+", "10s+", 48.71, 58.53}, {"11s+", "11s+", 30.70, 36.89},
  {"12s+", "12s+", 20.59, 24.74}, {"6p-", "6p-", 161.0, 201.6},
  {"7p-", "7p-", 57.65, 71.64},   {"8p-", "8p-", 27.09, 33.57},
  {"9p-", "9p-", 14.85, 18.38},   {"10p-", "10p-", 9.002, 11.13},
  {"11p-", "11p-", 5.863, 7.248}, {"12p-", "12p-", 4.030, 4.980}};

const std::vector<TestData> pra_fr{
  {"7s+", "7s+", 5929, 7040},     {"8s+", "8s+", 1520, 1802},
  {"9s+", "9s+", 624.0, 739.5},   {"10s+", "10s+", 316.8, 375.3},
  {"11s+", "11s+", 182.7, 216.5}, {"12s+", "12s+", 114.9, 136.1},
  {"7p-", "7p-", 628.2, 777.2},   {"8p-", "8p-", 222.9, 273.8},
  {"9p-", "9p-", 104.4, 127.9},   {"10p-", "10p-", 57.13, 69.90},
  {"11p-", "11p-", 34.61, 42.31}, {"12p-", "12p-", 22.53, 27.53}};

//==============================================================================
//! Unit tests External Field (RPA equations using TDHF method)
TEST_CASE("External Field: TDHF (RPA) for hyperfine",
          "[ExternalField][TDHF][RPA][HFS2][integration]") {
  {
    IO::ChronoTimer t{"TDHF"};

    // XXX These are not perfect, and should be improved!

    // Match data from: 10.1103/PhysRevA.100.042506
    const auto atom = std::vector<std::string>{"Rb", "Cs", "Fr"};
    const auto core = std::vector<std::string>{"[Kr]", "[Xe]", "[Rn]"};
    const auto A = std::vector{87, 133, 211};
    const auto rrms = std::vector{4.1989, 4.8041, 5.5876};
    const auto I = std::vector{1.5, 3.5, 4.5};
    const auto mu = std::vector{2.751818, 2.582025, 4.00};

    const auto test_data = std::vector{pra_rb, pra_cs, pra_fr};

    for (std::size_t i = 0; i < atom.size(); ++i) {

      // Create wavefunction object, solve HF for core+valence
      Wavefunction wf({6000, 1.0e-6, 400.0, 20.0, "loglinear", -1.0},
                      {atom[i], A[i], "Fermi", rrms[i], 2.3}, 1.0);
      wf.solve_core("HartreeFock", core[i]);
      wf.solve_valence("12sp");

      DiracOperator::hfs hfs(1, mu[i] / I[i], 0.0, wf.grid(),
                             DiracOperator::Hyperfine::pointlike_F());

      std::cout << "\nTable I from PhysRevA.100.042506, for " << atom[i] << "-"
                << A[i] << "\n";
      test_RPA2(wf, hfs, 0.0, test_data[i], 5.0e-4, 1.0e-2);
    }
  }
}

//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================

/*
 * Tests for TDHF, with frequency-dependent operator:
 * Calculate PNC few different ways:
 * <a|d|\delta b> vs. * <Y_a|h|b> (+cc)
 * Test a->b vs b->a
 * Do for E1 (L) and E1 (v)
 * Good test, because also imaginary operator!

*/

// z-component PNC amplitude f <- i using TDHF mixed states.
// h_source perturbs the states; h_me is the matrix element operator.
// Returns {term1, term2} where:
//   term1 = <f|h_me|X_i>
//   term2 = <Yf|h_me|i>
std::pair<double, double>
calc_mixed_states(DiracOperator::TensorOperator &h_source,
                  DiracOperator::TensorOperator &h_me,
                  ExternalField::TDHF &tdhf_source, const DiracSpinor &Ff,
                  const DiracSpinor &Fi, double omega_source,
                  const ExternalField::TDHF *dV_me = nullptr) {
  using namespace ExternalField;
  const bool emission = Ff.en() < Fi.en();
  const auto conj = emission;
  const int tja = Ff.twoj(), tjb = Fi.twoj(), twom = std::min(tja, tjb);
  const bool incl_dV = (dV_me != nullptr);

  // For freq-dependent ME operator (E1v): update with signed omega = Ef - Ei
  if (h_me.freqDependantQ())
    h_me.updateFrequency(Ff.en() - Fi.en());

  // For non-static source with emission: solve Y-type (not X) for ket,
  // and X-type (not Y) for bra. Equivalent to using conjugate perturbation.
  const bool swap_xy = emission && (omega_source != 0.0);
  const auto XY_i = swap_xy ? dPsiType::Y : dPsiType::X;
  const auto XY_f = swap_xy ? dPsiType::X : dPsiType::Y;

  // |Xi_n>: ket-type perturbed states of i under h_source at |omega_source|
  const auto Xi = tdhf_source.solve_dPsis(Fi, std::abs(omega_source), XY_i,
                                          nullptr, StateType::ket, incl_dV);
  // <Yf_n|: bra-type perturbed states of f under h_source
  const auto Yf = tdhf_source.solve_dPsis(Ff, std::abs(omega_source), XY_f,
                                          nullptr, StateType::bra, incl_dV);

  double t1 = 0.0, t2 = 0.0;

  // t1 = sum_n C(jf,jn,jb) * <Yf_n|h_me|i>
  //    ~ sum_n <f|h_src|n><n|h_me|i> / (Ef - En)   [h_src "first"]
  for (const auto &yf : Yf) {
    // C(jf,jn,jb): angular factor projecting RMEs to the z-component amplitude
    const auto c =
      h_source.rme3js(tja, yf.twoj(), twom) * h_me.rme3js(yf.twoj(), tjb, twom);
    t1 +=
      c * (h_me.reducedME(yf, Fi) + (dV_me ? dV_me->dV(yf, Fi, conj) : 0.0));
  }

  // t2 = sum_n C(jf,jn,jb) * <f|h_me|Xi_n>
  //    ~ sum_n <f|h_me|n><n|h_src|i> / (Ei - En)   [h_src "second"]
  for (const auto &xi : Xi) {
    const auto c =
      h_me.rme3js(tja, xi.twoj(), twom) * h_source.rme3js(xi.twoj(), tjb, twom);
    t2 +=
      c * (h_me.reducedME(Ff, xi) + (dV_me ? dV_me->dV(Ff, xi, conj) : 0.0));
  }
  return {t1, t2};
}

// z-component PNC amplitude f <- i via sum-over-states using a basis.
// Returns {t1, t2} where:
//   t1 = sum_n c10 <f||d||n><n||h||i>/(Ei-En)   (d first)
//   t2 = sum_n c01 <f||h||n><n||d||i>/(Ef-En)   (d second)
std::pair<double, double>
calc_sos(DiracOperator::TensorOperator &hd, DiracOperator::TensorOperator &hpnc,
         const DiracSpinor &Ff, const DiracSpinor &Fi,
         const std::vector<DiracSpinor> &basis,
         const ExternalField::TDHF *dV_d = nullptr,
         const ExternalField::TDHF *dV_pnc = nullptr) {
  // For freq-dependent operators (E1v): set signed omega to match method A ME
  if (hd.freqDependantQ())
    hd.updateFrequency(Ff.en() - Fi.en());
  const int tja = Ff.twoj(), tjb = Fi.twoj(), twom = std::min(tja, tjb);
  double t1 = 0.0, t2 = 0.0;
  for (const auto &Fn : basis) {
    const auto tjn = Fn.twoj();
    const auto c10 = hd.rme3js(tja, tjn, twom) * hpnc.rme3js(tjn, tjb, twom);
    const auto me_d_fn = hd.reducedME(Ff, Fn) + (dV_d ? dV_d->dV(Ff, Fn) : 0.0);
    const auto me_pnc_ni =
      hpnc.reducedME(Fn, Fi) + (dV_pnc ? dV_pnc->dV(Fn, Fi) : 0.0);
    t1 += c10 * me_d_fn * me_pnc_ni / (Fi.en() - Fn.en());
    const auto c01 = hpnc.rme3js(tja, tjn, twom) * hd.rme3js(tjn, tjb, twom);
    const auto me_pnc_fn =
      hpnc.reducedME(Ff, Fn) + (dV_pnc ? dV_pnc->dV(Ff, Fn) : 0.0);
    const auto me_d_ni = hd.reducedME(Fn, Fi) + (dV_d ? dV_d->dV(Fn, Fi) : 0.0);
    t2 += c01 * me_pnc_fn * me_d_ni / (Ff.en() - Fn.en());
  }
  return {t1, t2};
}

//==============================================================================
//! PNC amplitude: s-s and s-d transitions, E1 and E1v, forward and reverse.
/*!
  Prints three equivalent formulas for the z-component PNC amplitude, with
  and without dV (TDHF core polarisation), for Cs 6s<->7s and 7s<->5d-:
    A: delta states from PNC, d as ME operator
    B: X/Y states from d, PNC as ME operator (Y first = d "first")
    SOS: direct sum over basis states
  Tests forward and reverse directions, which must differ by sign only.
  WIP: A/B equivalence for freq-dep (E1v) s-d transition under investigation.
*/
TEST_CASE("PNC: amplitude methods (A, B, SOS)",
          "[ExternalField][TDHF][TDHFf][PNC][MixedStates][gauge]") {

  Wavefunction wf({2200, 1.0e-6, 120.0, 20.0, "loglinear", -1.0},
                  {"Cs", 133, "Fermi", 4.8041, 2.3}, 1.0);
  wf.solve_core("HartreeFock", "[Xe]");
  wf.solve_valence("7sp5d");
  wf.formBasis({"95p", 100, 9, 1.0e-4, 0.0, 75.0});

  const auto *p6s = wf.getState("6s+");
  const auto *p7s = wf.getState("7s+");
  const auto *p5dm = wf.getState("5d-");
  REQUIRE(p6s != nullptr);
  REQUIRE(p7s != nullptr);
  REQUIRE(p5dm != nullptr);

  DiracOperator::E1 hE1(wf.grid());
  DiracOperator::E1v hE1v(wf.alpha());
  DiracOperator::E1v hE1v_m(wf.alpha());
  DiracOperator::PNCnsi hPNC(5.67073, 2.3, wf.grid());

  ExternalField::TDHF tdhf_pnc(&hPNC, wf.vHF());
  ExternalField::TDHF tdhf_E1(&hE1, wf.vHF());
  ExternalField::TDHF tdhf_E1v(&hE1v, wf.vHF(), &hE1v_m);

  fmt::print("\nPNC amplitude, Cs-133, no dV:\n");

  using FIPair = std::pair<const DiracSpinor *, const DiracSpinor *>;
  const std::array<FIPair, 2> trans = {FIPair{p7s, p6s},   // s-s
                                       FIPair{p5dm, p6s}}; // s-d

  using OpTriple =
    std::tuple<DiracOperator::TensorOperator *, DiracOperator::TensorOperator *,
               ExternalField::TDHF *>;
  for (const auto &[hd, hd_m, tdhf_d] : std::array<OpTriple, 2>{
         OpTriple{&hE1, &hE1, &tdhf_E1}, OpTriple{&hE1v, &hE1v_m, &tdhf_E1v}}) {

    for (const auto &[pFf, pFi] : trans) {
      const auto *pRf = pFi; // reverse
      const auto *pRi = pFf;
      const auto w = pFf->en() - pFi->en(); // signed transition frequency
      if (hd->freqDependantQ()) {
        hd->updateFrequency(std::abs(w));
        hd_m->updateFrequency(-std::abs(w));
      }

      const bool emission_f = pFf->en() < pFi->en();
      const bool emission_r = pRf->en() < pRi->en();
      const auto label_f = fmt::format(
        "{} -> {} {} (w={:.4f}) {}", pFi->shortSymbol(), pFf->shortSymbol(),
        hd->name(), w, emission_f ? " [emission]" : "");
      const auto label_r = fmt::format(
        "{} -> {} {} (w={:.4f}) {}", pRi->shortSymbol(), pRf->shortSymbol(),
        hd->name(), -w, emission_r ? " [emission]" : "");

      // Method A: PNC source, d as ME; method B: d source, PNC as ME
      const auto A_f = calc_mixed_states(hPNC, *hd, tdhf_pnc, *pFf, *pFi, 0.0);
      if (hd->freqDependantQ())
        hd->updateFrequency(std::abs(w));
      const auto B_f = calc_mixed_states(*hd, hPNC, *tdhf_d, *pFf, *pFi, w);
      const auto A_r = calc_mixed_states(hPNC, *hd, tdhf_pnc, *pRf, *pRi, 0.0);
      if (hd->freqDependantQ())
        hd->updateFrequency(std::abs(w));
      const auto B_r = calc_mixed_states(*hd, hPNC, *tdhf_d, *pRf, *pRi, -w);
      const auto sos_f = calc_sos(*hd, hPNC, *pFf, *pFi, wf.basis());
      const auto sos_r = calc_sos(*hd, hPNC, *pRf, *pRi, wf.basis());

      const auto [a1_f, a2_f] = A_f; // a1=<df|d|i>, a2=<f|d|di>
      const auto [b2_f, b1_f] = B_f;
      const auto [t1_f, t2_f] = sos_f;
      const auto Af = a1_f + a2_f, Bf = b2_f + b1_f, Sf = t1_f + t2_f;
      // PNC source vs d source
      const auto eps_AB_f = 2.0 * std::abs((Af - Bf) / (Af + Bf));
      // mixed states vs SOS
      const auto eps_SOS_f = 2.0 * std::abs((Af - Sf) / (Af + Sf));
      fmt::print("\n  {}\n", label_f);
      fmt::print(
        "  A: <f|d|di>  = {:+11.4e}  <df|d|i> = {:+11.4e}  sum = {:+11.4e}\n",
        a2_f, a1_f, Af);
      fmt::print(
        "  B: {}  = {:+11.4e}  {} = {:+11.4e}  sum = {:+11.4e}  eps={:.1e}\n",
        emission_f ? "<Xf|h|i>" : "<Yf|h|i>", b2_f,
        emission_f ? "<f|h|Yi>" : "<f|h|Xi>", b1_f, Bf, eps_AB_f);
      fmt::print("  S: t1        = {:+11.4e}  t2       = {:+11.4e}  sum = "
                 "{:+11.4e}  eps={:.1e}\n",
                 t1_f, t2_f, Sf, eps_SOS_f);

      const auto [a1_r, a2_r] = A_r;
      const auto [b2_r, b1_r] = B_r;
      const auto [t1_r, t2_r] = sos_r;
      const auto Ar = a1_r + a2_r, Br = b2_r + b1_r, Sr = t1_r + t2_r;
      const auto eps_AB_r =
        2.0 * std::abs((Ar - Br) / (Ar + Br)); // PNC source vs d source
      const auto eps_SOS_r =
        2.0 * std::abs((Ar - Sr) / (Ar + Sr)); // mixed states vs SOS
      fmt::print("\n  {}\n", label_r);
      fmt::print(
        "  A: <f|d|di>  = {:+11.4e}  <df|d|i> = {:+11.4e}  sum = {:+11.4e}\n",
        a2_r, a1_r, Ar);
      fmt::print(
        "  B: {}  = {:+11.4e}  {} = {:+11.4e}  sum = {:+11.4e}  eps={:.1e}\n",
        emission_r ? "<Xf|h|i>" : "<Yf|h|i>", b2_r,
        emission_r ? "<f|h|Yi>" : "<f|h|Xi>", b1_r, Br, eps_AB_r);
      fmt::print("  S: t1        = {:+11.4e}  t2       = {:+11.4e}  sum = "
                 "{:+11.4e}  eps={:.1e}\n",
                 t1_r, t2_r, Sr, eps_SOS_r);

      const auto eps_sign =
        2.0 * std::abs((Af + Ar) / (Af - Ar)); // A(fwd) + A(rev) ~ 0
      fmt::print("  fwd/rev eps = {:.1e}\n", eps_sign);
    }
  }

  // === With dV (TDHF core polarisation) ===
  fmt::print("\nPNC amplitude, Cs-133, with dV:\n");

  tdhf_pnc.solve_core(0.0);

  for (const auto &[hd, hd_m, tdhf_d] : std::array<OpTriple, 2>{
         OpTriple{&hE1, &hE1, &tdhf_E1}, OpTriple{&hE1v, &hE1v_m, &tdhf_E1v}}) {

    for (const auto &[pFf, pFi] : trans) {
      const auto *pRf = pFi; // reverse
      const auto *pRi = pFf;
      const auto w = pFf->en() - pFi->en();
      const bool emission_f = pFf->en() < pFi->en();
      const bool emission_r = pRf->en() < pRi->en();
      const auto label_f = fmt::format(
        "{} -> {} {} (w={:.4f}) {}", pFi->shortSymbol(), pFf->shortSymbol(),
        hd->name(), w, emission_f ? " [emission]" : "");
      const auto label_r = fmt::format(
        "{} -> {} {} (w={:.4f}) {}", pRi->shortSymbol(), pRf->shortSymbol(),
        hd->name(), -w, emission_r ? " [emission]" : "");

      if (hd->freqDependantQ()) {
        hd->updateFrequency(std::abs(w));
        hd_m->updateFrequency(-std::abs(w));
      }

      std::cout << "\n";
      tdhf_d->solve_core(w);

      // A uses hd as ME operator (may update hd's freq to signed w)
      const auto A_f =
        calc_mixed_states(hPNC, *hd, tdhf_pnc, *pFf, *pFi, 0.0, tdhf_d);
      if (hd->freqDependantQ())
        hd->updateFrequency(std::abs(w));
      const auto B_f =
        calc_mixed_states(*hd, hPNC, *tdhf_d, *pFf, *pFi, w, &tdhf_pnc);
      const auto A_r =
        calc_mixed_states(hPNC, *hd, tdhf_pnc, *pRf, *pRi, 0.0, tdhf_d);
      if (hd->freqDependantQ())
        hd->updateFrequency(std::abs(w));
      const auto B_r =
        calc_mixed_states(*hd, hPNC, *tdhf_d, *pRf, *pRi, -w, &tdhf_pnc);
      const auto sos_f =
        calc_sos(*hd, hPNC, *pFf, *pFi, wf.basis(), tdhf_d, &tdhf_pnc);
      const auto sos_r =
        calc_sos(*hd, hPNC, *pRf, *pRi, wf.basis(), tdhf_d, &tdhf_pnc);

      const auto [a1_f, a2_f] = A_f; // a1=<df|d|i>, a2=<f|d|di>
      const auto [b2_f, b1_f] = B_f;
      const auto [t1_f, t2_f] = sos_f;
      const auto Af = a1_f + a2_f, Bf = b2_f + b1_f, Sf = t1_f + t2_f;
      const auto eps_AB_f = 2.0 * std ::abs((Af - Bf) / (Af + Bf));
      const auto eps_SOS_f = 2.0 * std::abs((Af - Sf) / (Af + Sf));
      fmt::print("\n  {}\n", label_f);
      fmt::print(
        "  A: <f|d|di>  = {:+11.4e}  <df|d|i> = {:+11.4e}  sum = {:+11.4e}\n",
        a2_f, a1_f, Af);
      fmt::print(
        "  B: {}  = {:+11.4e}  {} = {:+11.4e}  sum = {:+11.4e}  eps={:.1e}\n",
        emission_f ? "<Xf|h|i>" : "<Yf|h|i>", b2_f,
        emission_f ? "<f|h|Yi>" : "<f|h|Xi>", b1_f, Bf, eps_AB_f);
      fmt::print("  S: t1        = {:+11.4e}  t2       = {:+11.4e}  sum = "
                 "{:+11.4e}  eps={:.1e}\n",
                 t1_f, t2_f, Sf, eps_SOS_f);

      const auto [a1_r, a2_r] = A_r;
      const auto [b2_r, b1_r] = B_r;
      const auto [t1_r, t2_r] = sos_r;
      const auto Ar = a1_r + a2_r, Br = b2_r + b1_r, Sr = t1_r + t2_r;
      const auto eps_AB_r =
        2.0 * std::abs((Ar - Br) / (Ar + Br)); // PNC source vs d source
      const auto eps_SOS_r =
        2.0 * std::abs((Ar - Sr) / (Ar + Sr)); // mixed states vs SOS
      fmt::print("\n  {}\n", label_r);
      fmt::print(
        "  A: <f|d|di>  = {:+11.4e}  <df|d|i> = {:+11.4e}  sum = {:+11.4e}\n",
        a2_r, a1_r, Ar);
      fmt::print(
        "  B: {}  = {:+11.4e}  {} = {:+11.4e}  sum = {:+11.4e}  eps={:.1e}\n",
        emission_r ? "<Xf|h|i>" : "<Yf|h|i>", b2_r,
        emission_r ? "<f|h|Yi>" : "<f|h|Xi>", b1_r, Br, eps_AB_r);
      fmt::print("  S: t1        = {:+11.4e}  t2       = {:+11.4e}  sum = "
                 "{:+11.4e}  eps={:.1e}\n",
                 t1_r, t2_r, Sr, eps_SOS_r);

      const auto eps_sign =
        2.0 * std::abs((Af + Ar) / (Af - Ar)); // A(fwd) + A(rev) ~ 0
      fmt::print("  fwd/rev eps = {:.1e}\n", eps_sign);

      CHECK(eps_AB_f < 5.0e-3);
      CHECK(eps_AB_r < 5.0e-3);
      CHECK(eps_sign < 5.0e-3);
    }
  }
}

//==============================================================================
//! Tests E1 gauge equivalence (length = velocity form) and hermiticity,
//! over all valence pairs, with TDHF core polarisation.
/*!
  Loops over all ordered valence pairs (a,b) with non-zero E1 matrix element.
  For each pair, omega = E_a - E_b (signed: positive = absorption, negative =
  emission). Computes E1 (length) and E1v (velocity) matrix elements with TDHF
  and checks gauge equivalence inline. Stores results; at end, checks that each
  pair (a,b) and its reverse (b,a) satisfy hermiticity:
    <a||h||b> = symm_sign(a,b) * <b||h||a>
  For E1v, hermiticity holds between (a,b) at +omega and (b,a) at -omega,
  since E1v scales as 1/omega and the TDHF correction is symmetric in omega.
*/
TEST_CASE("External Field: TDHF E1 length-velocity gauge and hermiticity",
          "[ExternalField][TDHF][TDHFf][gauge]") {

  Wavefunction wf({2000, 1.0e-5, 75.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-10, false);
  wf.solve_valence("3spd", false);

  DiracOperator::E1 h_L(wf.grid());

  DiracOperator::E1v h_V_p(wf.alpha(), 0.0);
  DiracOperator::E1v h_V_m(wf.alpha(), 0.0);

  ExternalField::TDHF rpa_L(&h_L, wf.vHF());
  ExternalField::TDHF rpa_V(&h_V_p, wf.vHF(), &h_V_m);

  fmt::print("\nE1 length-velocity gauge equivalence (Na, TDHF)\n");

  for (const auto &Fv : wf.valence()) {
    for (const auto &Fw : wf.valence()) {
      if (h_L.isZero(Fv, Fw))
        continue;
      if (Fv.en() < Fw.en())
        continue;

      double absME = 0.0;
      double emmME = 0.0;

      for (const auto abs : {true, false}) {
        const auto &Fa = abs ? Fv : Fw;
        const auto &Fb = abs ? Fw : Fv;

        const double omega = Fa.en() - Fb.en();
        h_V_p.updateFrequency(std::abs(omega));
        h_V_m.updateFrequency(-std::abs(omega));

        std::cout << "\n";
        fmt::print("{:3s} -> {:3s}  omega = {:+.5f}  ({})\n", Fa.shortSymbol(),
                   Fb.shortSymbol(), omega,
                   omega > 0.0 ? "absorption" : "emission");
        rpa_L.solve_core(omega);
        rpa_V.solve_core(omega);

        const auto me_L0 = h_L.reducedME(Fa, Fb);
        const auto me_V0 =
          abs ? h_V_p.reducedME(Fa, Fb) : h_V_m.reducedME(Fa, Fb);
        const auto dV_L = rpa_L.dV(Fa, Fb);
        const auto dV_V = rpa_V.dV(Fa, Fb);
        const auto me_L = me_L0 + dV_L;
        const auto me_V = me_V0 + dV_V;
        const auto eps = std::abs(2.0 * (me_L - me_V) / (me_L + me_V));

        fmt::print("E1(L)   : {:10.6f} + {:10.6f} = {:10.6f}\n", me_L0, dV_L,
                   me_L);
        fmt::print("E1(V)   : {:10.6f} + {:10.6f} = {:10.6f}\n", me_V0, dV_V,
                   me_V);
        fmt::print("eps(L/V): {:.1e}\n", eps);

        REQUIRE(me_L == Approx(me_V).epsilon(1.0e-4));
        if (abs) {
          absME = me_V;
        } else {
          emmME = h_V_p.symm_sign(Fa, Fb) * me_V;
        }
      }
      const auto eps = std::abs(2.0 * (absME - emmME) / (absME + emmME));
      fmt::print("\nHermicity:\neps(+/-): {:.1e}\n", eps);
      // Goes to zero if we re-start RPA from scratch!
      REQUIRE(eps < 1.0e-6);
    }
  }
}