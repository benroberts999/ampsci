#pragma once
#include "DiracOperator/Operators.hpp"
#include "Physics/RadPot.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace UnitTest {

//******************************************************************************
namespace helper {

struct QEDData {
  struct Atom {
    std::string name;
    int A;
    std::string core;
    std::string val;
  };
  struct Data {
    std::string state;
    double de0;
    double de;
  };
  Atom at;
  std::vector<Data> d;
};

std::pair<std::string, double> compare_QED(const std::vector<QEDData> &QEDdata,
                                           QED::RadPot::Scale scale);

bool FGRadPot(std::ostream &obuff);

} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests for Ginges/Flambaum Radiative potential method
bool RadPot(std::ostream &obuff) {
  bool pass = true;

  // Run unit tests (tests integrals vs Mathematica)
  pass &= helper::FGRadPot(obuff);

  // Uehling data, from: Ginges, Berengut, J. Phys. B 49, 095001 (2016).
  const std::vector<helper::QEDData> Ueh_data = //
      {{{"Na", -1, "[Ne]", "3spd"},
        {{"3s+", -5.5590E-07, -5.9100E-07},
         {"3p-", -2.0060E-10, 3.0670E-08},
         {"3p+", -4.1010E-11, 3.0770E-08},
         {"3d-", -6.5380E-18, -1.3230E-10},
         {"3d+", -2.0180E-18, -1.3470E-10}}},
       {{"K", -1, "[Ar]", "4sp3d"},
        {{"4s+", -1.2240E-06, -1.3330E-06},
         {"4p-", -1.8790E-09, 7.2620E-08},
         {"4p+", -3.5470E-10, 7.3960E-08},
         {"3d-", -1.0460E-14, 1.3950E-08},
         {"3d+", -3.1120E-15, 1.3760E-08}}},
       {{"Rb", -1, "[Kr]", "5sp4d"},
        {{"5s+", -4.6480E-06, -5.1140E-06},
         {"5p-", -3.2340E-08, 1.7760E-07},
         {"5p+", -4.8810E-09, 2.1080E-07},
         {"4d-", -3.8650E-12, 1.2900E-07},
         {"4d+", -1.0480E-12, 1.2180E-07}}},
       {{"Cs", -1, "[Xe]", "6sp5d"},
        {{"6s+", -1.0540E-05, -1.1600E-05},
         {"6p-", -1.9420E-07, 2.2880E-07},
         {"6p+", -2.1780E-08, 4.5180E-07},
         {"5d-", -1.9380E-10, 1.1270E-06},
         {"5d+", -4.6120E-11, 1.0070E-06}}},
       {{"Fr", 211, "[Rn]", "7sp6d"},
        {{"7s+", -4.9970E-05, -5.5400E-05},
         {"7p-", -2.5690E-06, -1.6100E-06},
         {"7p+", -1.3560E-07, 1.8340E-06},
         {"6d-", -3.5790E-09, 3.8220E-06},
         {"6d+", -6.9210E-10, 2.6060E-06}}},
       {{"E119", -1, "[Og]", "8sp7d"},
        {{"8s+", -3.4270E-04, -4.0460E-04},
         {"8p-", -3.9210E-05, -4.7010E-05},
         {"8p+", -5.4200E-07, 1.1200E-05},
         {"7d-", -2.5170E-08, 8.7150E-06},
         {"7d+", -4.4200E-09, -1.3740E-06}}}};

  // Self-Energy data, from: Ginges, Berengut, Phys. Rev. A 93, 052509 (2016).
  const std::vector<helper::QEDData> SE_data = //
      {{{"Na", -1, "[Ne]", "3spd"},
        {{"3s+", 1.068E-05, 1.125E-05},
         // Exclude these: so small, can never reproduce exactly
         // do reproduce pretty well, however
         // ,{"3d-", -2.644E-09, 4.085E-11}
         // ,{"3d+", 1.770E-09, 3.407E-09}
         {"3p-", -5.698E-08, -7.088E-07},
         {"3p+", 1.112E-07, -4.882E-07}}},
       {{"K", -1, "[Ar]", "4sp3d"},
        {{"4s+", 1.845E-05, 1.974E-05},
         {"4p-", -1.023E-07, -1.400E-06},
         {"4p+", 3.072E-07, -8.411E-07},
         {"3d-", -1.915E-08, -2.568E-07},
         {"3d+", 1.382E-08, -2.647E-07}}},
       {{"Rb", -1, "[Kr]", "5sp4d"},
        {{"5s+", 4.836E-05, 5.136E-05},
         {"5p-", 1.381E-08, -2.854E-06},
         {"5p+", 1.316E-06, -1.083E-06},
         {"4d-", -1.098E-07, -1.791E-06},
         {"4d+", 1.005E-07, -1.790E-06}}},
       {{"Cs", -1, "[Xe]", "6sp5d"},
        {{"6s+", 8.128E-05, 8.431E-05},
         {"6p-", 1.077E-06, -3.831E-06},
         {"6p+", 3.183E-06, -9.203E-07},
         {"5d-", -6.066E-07, -1.212E-05},
         {"5d+", 7.174E-07, -1.115E-05}}},
       {{"Fr", 211, "[Rn]", "7sp6d"},
        {{"7s+", 2.201E-04, 2.166E-04},
         // Exclude these: so small, can never reproduce exactly
         // do reproduce pretty well, however
         // {"7p-", 1.068E-05, 1.276E-09},
         // {"7p+", 1.034E-05, -5.888E-08},
         {"6d-", -8.046E-07, -2.668E-05},
         {"6d+", 1.968E-06, -2.247E-05}}},
       {{"E119", -1, "[Og]", "8sp7d"},
        {{"8s+", 7.196E-04, 6.832E-04},
         {"8p-", 7.204E-05, 5.217E-05},
         {"8p+", 2.697E-05, -2.001E-06},
         {"7d-", 4.523E-07, -4.069E-05},
         {"7d+", 4.717E-06, -3.086E-05}}}};

  // "Regresseion" tests: Uehling:
  {
    std::cout << "Test Uehling: \nReproduce Table IV from Ginges, Berengut, J. "
                 "Phys. B 49, 095001 (2016)\n";
    std::cout
        << "        de(0)     [ expected ] |  de(HF)    [ expected ]  eps\n";
    const auto [at_worst, eps_worst] =
        helper::compare_QED(Ueh_data, {1.0, 0.0, 0.0, 0.0, 0.0});
    pass &= qip::check_value(&obuff, "Uehling (TabIV) " + at_worst, eps_worst,
                             0.0, 6.0e-4);
  }

  // "Regresseion" tests: Self-energy:
  {
    std::cout << "\nTest Self-energy: \nReproduce Table VI from Ginges, "
                 "Berengut, Phys. Rev. A 93, 052509 (2016)\n";
    std::cout
        << "        de(0)     [ expected ] |  de(HF)    [ expected ]  eps\n";
    const auto [at_worst, eps_worst] =
        helper::compare_QED(SE_data, {0.0, 1.0, 1.0, 1.0, 0.0});
    pass &= qip::check_value(&obuff, "Self-Energy (TabVI) " + at_worst,
                             eps_worst, 0.0, 6.0e-4);
  }

  // "Regresseion" tests: Self-energy (parts):
  {
    std::cout << "\nSelf-energy components for Cs (Core-Hartree)\n";
    std::cout << "c.set_f(). Ginges, Berengut, Phys. Rev. A 93, 052509 (2016), "
                 "TabII\n";
    // Construct wavefunction, solve HF core+valence (without QED):
    Wavefunction wf({2000, 1.0e-6, 120.0, 0.3 * 120.0, "loglinear", -1.0},
                    {"Cs", 133, "Fermi", -1, 2.3}, 1.0);
    wf.solve_core("Hartree", 0.0, "[Xe]", 0.0, false);
    wf.solve_valence("6s", false);
    std::vector pt = {1.391, 6.583, 0.762};
    std::vector st = {1.391, 6.577, 0.762};

    const auto r_N_au = std::sqrt(5.0 / 3.0) * 1.0 *
                        wf.get_nuclearParameters().r_rms / PhysConst::aB_fm;
    // Rad pot, pointlike and step nucleus
    const auto rp0 = QED::RadPot(wf.rgrid->r(), wf.Znuc(), 0.0, 15.0,
                                 QED::RadPot::Scale{0.0, 1.0, 1.0, 1.0, 0.0},
                                 std::vector{1.0}, false, false);
    const auto rp = QED::RadPot(wf.rgrid->r(), wf.Znuc(), r_N_au, 15.0,
                                QED::RadPot::Scale{0.0, 1.0, 1.0, 1.0, 0.0},
                                std::vector{1.0}, false, false);

    std::cout
        << "                 mag   [expct], high  [expct], low   [expct]\n";
    const auto lam = [&](auto &exp, auto &RP, std::string name) {
      const auto hm = DiracOperator::Hrad({}, RP.Hmag(0));
      const auto hh = DiracOperator::Hrad(RP.Vh(0), {});
      const auto hl = DiracOperator::Hrad(RP.Vl(0), {});
      const auto &Fv = wf.valence.at(0);

      const auto m = hm.radialIntegral(Fv, Fv) * 1.0e5;
      const auto h = hh.radialIntegral(Fv, Fv) * 1.0e5;
      const auto l = hl.radialIntegral(Fv, Fv) * 1.0e5;

      const auto eps = (h + m + l) / (exp[0] + exp[1] + exp[2]) - 1.0;

      printf("%15s: %.3f [%.3f], %.3f [%.3f], %.3f [%.3f]\n", name.c_str(), m,
             exp[0], h, exp[1], l, exp[2]);

      pass &= qip::check_value(&obuff, "Self-Energy (II) " + name, eps, 0.0,
                               1.0e-4);
    };

    lam(pt, rp0, "point");
    lam(st, rp, "step");
  }

  return pass;
}

//******************************************************************************
std::pair<std::string, double>
helper::compare_QED(const std::vector<QEDData> &QEDdata,
                    QED::RadPot::Scale scale) {
  // Store worst comparison for unit test
  double eps_worst = 0.0;
  std::string at_worst;

  // Loop through atoms:
  for (const auto &[atom, data] : QEDdata) {
    std::cout << atom.name << ":\n";

    // Must use exact same rms value as above papers
    const double rrms = atom.name == "E119" ? 6.5 : -1.0;

    // Construct wavefunction, solve HF core+valence (without QED):
    Wavefunction wf0({5000, 1.0e-6, 120.0, 0.3 * 120.0, "loglinear", -1.0},
                     {atom.name, atom.A, "Fermi", rrms, 2.3}, 1.0);
    wf0.solve_core("HartreeFock", 0.0, atom.core, 0.0, false);
    wf0.solve_valence(atom.val, false);

    // Construct wavefunction, solve HF core+valence (with QED):
    Wavefunction wf({5000, 1.0e-6, 120.0, 0.3 * 120.0, "loglinear", -1.0},
                    {atom.name, atom.A, "Fermi", rrms, 2.3}, 1.0);
    const auto rcut = 15.0;
    wf.radiativePotential(scale, rcut, 1.0, {1.0}, false, false);
    wf.solve_core("HartreeFock", 0.0, atom.core, 0.0, false);
    wf.solve_valence(atom.val, false);

    // Loop through each state (line in table) and compare data:
    for (auto [state, t_de0, t_de] : data) {

      const auto Fv0 = wf0.getState(state);
      const auto Fv = wf.getState(state);
      if (!Fv0 || !Fv)
        continue;

      // Operator (for first-order shift) - l-dependent
      auto h = DiracOperator::Hrad(wf.qed->Vel(Fv->l()), wf.qed->Hmag(Fv->l()));

      // Calculate first-order Uehling shift (matrix element):
      const auto de0 = h.radialIntegral(*Fv, *Fv);
      // Find Hartree-Fock Uehling shift:
      const auto de = Fv->en() - Fv0->en();
      // print results:
      printf("%4s | %10.3e [%10.3e] | %10.3e [%10.3e]", state.c_str(), de0,
             t_de0, de, t_de);
      const auto eps0 = (de0 - t_de0) / (t_de0);
      const auto eps = (de - t_de) / (t_de);
      printf("  %6.0e %6.0e\n", eps0, eps);

      if (std::abs(eps) > eps_worst) {
        eps_worst = std::abs(eps);
        at_worst = atom.name + " " + state;
      }
    }
  }
  return {at_worst, eps_worst};
}

//******************************************************************************
//******************************************************************************
// Test t integrals compared to Mathematica:
bool helper::FGRadPot(std::ostream &obuff) {
  bool pass = true;

  {
    // Compare Uehling integral (J function) to Mathematica
    const double rN = 1.0e-5;
    const auto r_list =
        std::vector{1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2};
    const auto expected = std::vector{0.1705223366, 1.699101548, 10.65108690,
                                      5.825117247,  1.748103780, 0.02520380190};
    std::vector<double> res;
    for (auto r : r_list) {
      res.push_back(FGRP::t_integral(FGRP::Uehling::J_Ueh_gsl, {r, rN}));
    }
    const auto [worst, at] = qip::compare_eps(res, expected);
    pass &=
        qip::check_value(&obuff, "FGRP::Ueh (cf Mathem)", worst, 0.0, 1.0e-6);
  }

  {
    // Compare Mag integral (J function) to Mathematica
    const double rN = 1.0e-5;
    const auto r_list =
        std::vector{1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 0.1};
    const auto expected =
        std::vector{-2.441802648e-8, -0.00001909634519, -0.001578897130,
                    -0.07268613382,  -0.8494339412,     -1.000000000};
    std::vector<double> res;
    for (auto r : r_list) {
      res.push_back(FGRP::t_integral(FGRP::Magnetic::J_mag_gsl, {r, rN}));
    }
    const auto [worst, at] = qip::compare(res, expected);
    pass &=
        qip::check_value(&obuff, "FGRP::Mag (cf Mathem)", worst, 0.0, 1.0e-7);
  }

  {
    // Compare Self-energy (low-f) integral (F function) to Mathematica
    const double z = 55.0;
    const double rN = 1.0e-4;
    const auto r_list =
        std::vector{1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 0.1};
    const auto expected =
        std::vector{0.995883787662840, 0.995856739485738, 0.993424136749548,
                    0.946383897482915, 0.576945209193922, 0.00408677930550273};
    std::vector<double> res;
    for (auto r : r_list) {
      res.push_back(FGRP::SE::F_SEl(z, r, rN));
    }
    const auto [worst, at] = qip::compare_eps(res, expected);
    pass &=
        qip::check_value(&obuff, "FGRP::SEl (cf Mathem)", worst, 0.0, 1.0e-11);
  }

  {
    // Compare Self-energy (high-f) integral (J function) to Mathematica
    // nb: This is not great agreement; not sure what the problem is..
    const double z = 55.0;
    const double rN = 1.0e-5;
    const auto r_list =
        std::vector{1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2};
    const auto expected =
        std::vector{0.0552394653230691, 0.551391255370585, 4.63893363508357,
                    10.3617772959789,   3.70455183897125,  0.0350009486189897};
    std::vector<double> res;
    for (auto r : r_list) {
      double x = FGRP::t_integral(FGRP::SE::J_SE_gsl, {r, rN, z});
      res.push_back(x);
    }
    const auto [worst, at] = qip::compare_eps(res, expected);
    pass &=
        qip::check_value(&obuff, "FGRP::SEh (cf Mathem)", worst, 0.0, 1.0e-3);
  }

  {
    // Test of integration method
    const auto expected = 0.0285214358159077;
    const auto expected2 = 0.0231175704164569;
    const auto f1 = [](auto x) {
      return std::sin(35.0 * x) * std::exp(-1.7 * x);
    };
    const auto f2 = [](auto x) {
      return std::sin(35.0 * x) * std::exp(-17.0 * x);
    };
    const auto a1 = FGRP::r_integral(f1, 0.0, 4.2, 1000);
    const auto a2 = FGRP::r_integral<FGRP::IntType::Log>(f2, 1.0e-7, 1.0, 500);

    const auto worst = std::max(std::abs(a1 - expected) / expected,
                                std::abs(a2 - expected2) / expected2);

    pass &=
        qip::check_value(&obuff, "FGRP::Int (cf Mathem)", worst, 0.0, 1.0e-6);
  }

  return pass;
}

} // namespace UnitTest
