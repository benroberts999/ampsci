#include "TDHF.hpp"
#include "DiracOperator/include.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
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

  const auto dv_0 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv_0 << "\n";
  REQUIRE(dv_0 != 0.0);

  // // Doing 1 iteration shuold be equiv to first-order dV
  // auto rpa1 = ExternalField::TDHF(&dE1, wf.vHF());
  // rpa1.solve_core(0.0, 1);

  const auto &Fc = wf.core().back();
  auto dPsis1 = rpa.solve_dPsis(Fc, 0.0, ExternalField::dPsiType::X);
  const auto &dPsis2 = rpa.get_dPsis(Fc, ExternalField::dPsiType::X);
  // should both be first-order:
  // const auto &dPsis02 = rpa1.get_dPsis(Fc, ExternalField::dPsiType::X);
  REQUIRE(dPsis1.size() == dPsis2.size());
  REQUIRE(dPsis1.size() != 0);
  for (std::size_t i = 0; i < dPsis1.size(); ++i) {
    // only true if TDHF converged!
    REQUIRE(dPsis1.at(i).norm2() ==
            Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
    const auto kappa = dPsis1.at(i).kappa();
    auto dPsi1 = rpa.solve_dPsi(Fc, 0.0, ExternalField::dPsiType::X, kappa);
    const auto &dPsi2 = rpa.get_dPsi_x(Fc, ExternalField::dPsiType::X, kappa);
    REQUIRE(dPsi1.norm2() == Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
    REQUIRE(dPsi2.norm2() == Approx(dPsis2.at(i).norm2()).epsilon(1.0e-2));
  }

  // should be zero after clearing:
  rpa.clear();
  const auto dv_00 = rpa.dV(*F6s, *F6p);
  std::cout << *F6s << "-" << *F6p << ": " << dv_00 << "\n";
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

    const auto f =
      h.name() == "hfs1" ?
        DiracOperator::Hyperfine::convert_RME_to_AB(1, Fa.kappa(), Fb.kappa()) :
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

    const auto f =
      DiracOperator::Hyperfine::convert_RME_to_AB(1, Fa.kappa(), Fb.kappa());
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
    std::cout << "-------------------------------------\n";
    std::cout << "External Field: TDHF (RPA)\n";

    // Create wavefunction object, solve HF for core+valence
    Wavefunction wf({6000, 1.0e-6, 175.0, 20.0, "loglinear", -1.0},
                    {"Cs", 133, "Fermi", 4.8041, 2.3}, 1.0);
    wf.solve_core("HartreeFock", 0.0, "[Xe]");
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
    std::cout << "-------------------------------------\n";
    std::cout << "External Field: TDHF (RPA) for hyperfine\n";

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
      wf.solve_core("HartreeFock", 0.0, core[i]);
      wf.solve_valence("12sp");

      DiracOperator::hfs hfs(1, mu[i] / I[i], 0.0, wf.grid(),
                             DiracOperator::Hyperfine::pointlike_F());

      std::cout << "\nTable I from PhysRevA.100.042506, for " << atom[i] << "-"
                << A[i] << "\n";
      test_RPA2(wf, hfs, 0.0, test_data[i], 5.0e-4, 1.0e-2);
    }
  }
}