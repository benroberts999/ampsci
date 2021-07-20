#pragma once
#include "DiracOperator/Operators.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <string>
#include <utility>
#include <vector>

namespace UnitTest {

//! Unit tests for Breit
bool Breit(std::ostream &obuff) {
  bool pass = true;

  // Solve Hartree-Fock, including Breit
  Wavefunction wf({3500, 1.0e-6, 125.0, 40.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  const double x_Breit = 1.0;
  wf.solve_core("HartreeFock", x_Breit, "[Xe]");
  wf.solve_valence("7sp");

  // Lambda to compare against (From Vladimir's code):
  // Overall sign difference between E1, E2, and PNC ME definition (and <ab> vs
  // <ba>) taken into account
  // Store data in list of pairs.
  // Sort the data by states, to make comparisons easier
  using datav = std::vector<std::pair<std::string, double>>;
  auto sort_by_first = [](auto x, auto y) { return x.first < y.first; };
  auto eps_second = [](const auto &x, const auto &y) {
    return (x.second - y.second) / y.second;
  };

  //****************************************************************************
  { // Test Breit energies:

    // Test data: From Dzuba calculation: energies (with Breit)
    auto en_VD = datav{{"6s_1/2", -0.12735352},    {"7s_1/2", -0.05518246},
                       {"6p_1/2", -0.08558175},    {"6p_3/2", -0.08377240},
                       {"7p_1/2", -0.04200916},    {"7p_3/2", -0.04136327},
                       {"1s_1/2", -1326.96780222}, {"2s_1/2", -212.26379763},
                       {"2p_1/2", -198.91648966},  {"2p_3/2", -186.09010791},
                       {"3s_1/2", -45.92577676},   {"3p_1/2", -40.36742765},
                       {"3p_3/2", -37.84531317},   {"3d_3/2", -28.28379474},
                       {"3d_5/2", -27.76338201},   {"4s_1/2", -9.50647336},
                       {"4p_1/2", -7.43335889},    {"4p_3/2", -6.91467392},
                       {"4d_3/2", -3.48506227},    {"4d_5/2", -3.39873098},
                       {"5s_1/2", -1.48931551},    {"5p_1/2", -0.90678723},
                       {"5p_3/2", -0.84005768}};
    // Sort, so don't worry about order:
    std::sort(begin(en_VD), end(en_VD), sort_by_first);

    // Get my values:
    datav en{};
    for (const auto &Fc : wf.core)
      en.emplace_back(Fc.symbol(), Fc.en());
    for (const auto &Fv : wf.valence)
      en.emplace_back(Fv.symbol(), Fv.en());
    std::sort(begin(en), end(en), sort_by_first);

    const auto [eps, at] = qip::compare(en, en_VD, eps_second);
    pass &= qip::check_value(&obuff, "Energies " + at->first, eps, 0.0, 3.0e-6);
  }

  //****************************************************************************
  { // Test E1

    // Test data (Dzuba), w/ Breit, no RPA
    auto e1_VD =
        datav{{"6p-6s+", -5.27804}, {"6p+6s+", 7.42721},  {"7p-6s+", -0.37356},
              {"7p+6s+", 0.695359}, {"6p-7s+", 4.41774},  {"6p+7s+", -6.67286},
              {"7p-7s+", -11.0078}, {"7p+7s+", 15.3455},  {"6s+6p-", -5.27804},
              {"7s+6p-", 4.41774},  {"6s+6p+", -7.42721}, {"7s+6p+", 6.67286},
              {"6s+7p-", -0.37356}, {"7s+7p-", -11.0078}, {"6s+7p+", -0.695359},
              {"7s+7p+", -15.3455}};
    std::sort(begin(e1_VD), end(e1_VD), sort_by_first);

    // Test data (Dzuba): E1 + RPA
    auto e1_VD_RPA = datav{
        {"6p-6s+", -4.97441},  {"6p+6s+", 7.01323},  {"7p-6s+", -0.240353},
        {"7p+6s+", 0.509077},  {"6p-7s+", 4.45391},  {"6p+7s+", -6.71407},
        {"7p-7s+", -10.9199},  {"7p+7s+", 15.228},   {"6s+6p-", -4.97441},
        {"7s+6p-", 4.45391},   {"6s+6p+", -7.01323}, {"7s+6p+", 6.71407},
        {"6s+7p-", -0.240353}, {"7s+7p-", -10.9199}, {"6s+7p+", -0.509077},
        {"7s+7p+", -15.228}};
    std::sort(begin(e1_VD_RPA), end(e1_VD_RPA), sort_by_first);

    // Solve TDHF with Breit (for RPA)
    const auto h{DiracOperator::E1(*wf.rgrid)};
    auto rpa = ExternalField::TDHF(&h, wf.getHF());
    rpa.solve_core(0.0, 20); // w=0

    // Get my values:
    datav me{}, me_RPA{};
    for (const auto &Fv : wf.valence) {
      for (const auto &Fw : wf.valence) {
        if (h.isZero(Fv.k, Fw.k))
          continue;
        auto h_vw = h.reducedME(Fv, Fw);
        auto rpa_vw = h_vw + rpa.dV(Fv, Fw);
        me.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(), h_vw);
        me_RPA.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(), rpa_vw);
      }
    }
    std::sort(begin(me), end(me), sort_by_first);
    std::sort(begin(me_RPA), end(me_RPA), sort_by_first);

    const auto [eps1, at1] = qip::compare(me, e1_VD, eps_second);
    const auto [eps2, at2] = qip::compare(me_RPA, e1_VD_RPA, eps_second);
    pass &=
        qip::check_value(&obuff, "E1       " + at1->first, eps1, 0.0, 5.0e-5);
    pass &=
        qip::check_value(&obuff, "E1 (RPA) " + at2->first, eps2, 0.0, 1.0e-4);
  }

  //****************************************************************************
  {
    // E2 (w/ Breit, no RPA)
    auto e2_VD =
        datav{{"6p+6p-", 68.5272},  {"7p+6p-", -42.5844}, {"6p-6p+", -68.5272},
              {"6p+6p+", 70.2513},  {"7p-6p+", 49.371},   {"7p+6p+", -46.8283},
              {"6p+7p-", -49.371},  {"7p+7p-", 300.011},  {"6p-7p+", 42.5844},
              {"6p+7p+", -46.8283}, {"7p-7p+", -300.011}, {"7p+7p+", 305.189}};
    std::sort(begin(e2_VD), end(e2_VD), sort_by_first);

    // E2 + RPA (w/ Breit and RPA)
    auto e2_VD_RPA =
        datav{{"6p+6p-", 68.2769},  {"7p+6p-", -42.6923}, {"6p-6p+", -68.2769},
              {"6p+6p+", 70.0047},  {"7p-6p+", 49.4739},  {"7p+6p+", -46.9332},
              {"6p+7p-", -49.4739}, {"7p+7p-", 299.938},  {"6p-7p+", 42.6923},
              {"6p+7p+", -46.9332}, {"7p-7p+", -299.938}, {"7p+7p+", 305.117}};
    std::sort(begin(e2_VD_RPA), end(e2_VD_RPA), sort_by_first);

    auto h = DiracOperator::Ek(*wf.rgrid, 2);
    auto rpa = ExternalField::TDHF(&h, wf.getHF());
    rpa.solve_core(0.0, 20); // w=0

    datav me{}, me_RPA{};
    for (const auto &Fv : wf.valence) {
      for (const auto &Fw : wf.valence) {
        if (h.isZero(Fv.k, Fw.k))
          continue;
        auto h_vw = h.reducedME(Fv, Fw);
        auto rpa_vw = h_vw + rpa.dV(Fv, Fw);
        me.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(), h_vw);
        me_RPA.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(), rpa_vw);
      }
    }
    std::sort(begin(me), end(me), sort_by_first);
    std::sort(begin(me_RPA), end(me_RPA), sort_by_first);

    const auto [eps1, at1] = qip::compare(me, e2_VD, eps_second);
    const auto [eps2, at2] = qip::compare(me_RPA, e2_VD_RPA, eps_second);
    pass &=
        qip::check_value(&obuff, "E2       " + at1->first, eps1, 0.0, 1.0e-5);
    pass &=
        qip::check_value(&obuff, "E2 (RPA) " + at2->first, eps2, 0.0, 2.0e-4);
  }

  //****************************************************************************
  {
    // PNC (w/ Breit, no RPA)
    auto pnc_VD = datav{{"6p-6s+", -5.70956e-4}, {"7p-6s+", -3.41811e-4},
                        {"6p-7s+", -2.99278e-4}, {"7p-7s+", -1.79167e-4},
                        {"6s+6p-", 5.70956e-4},  {"7s+6p-", 2.99278e-4},
                        {"6s+7p-", 3.41811e-4},  {"7s+7p-", 1.79167e-4}};
    std::sort(begin(pnc_VD), end(pnc_VD), sort_by_first);
    // PNC + RPA (w/ Breit and RPA)
    auto pnc_VD_RPA = datav{{"6p-6s+", -7.22577e-4}, {"7p-6s+", -4.30375e-4},
                            {"6p-7s+", -3.77096e-4}, {"7p-7s+", -2.24641e-4},
                            {"6s+6p-", 7.22577e-4},  {"7s+6p-", 3.77096e-4},
                            {"6s+7p-", 4.30375e-4},  {"7s+7p-", 2.24641e-4}};
    std::sort(begin(pnc_VD_RPA), end(pnc_VD_RPA), sort_by_first);

    // nb: use exact same 'c' and 't' params used for Dzuba test data:
    auto h = DiracOperator::PNCnsi(5.674800, 2.3, *wf.rgrid);
    auto rpa = ExternalField::TDHF(&h, wf.getHF());
    rpa.solve_core(0.0, 20); // w=0

    datav me{}, me_RPA{};
    for (const auto &Fv : wf.valence) {
      for (const auto &Fw : wf.valence) {
        if (h.isZero(Fv.k, Fw.k))
          continue;
        auto h_vw = h.reducedME(Fv, Fw);
        auto rpa_vw = h_vw + rpa.dV(Fv, Fw);
        me.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(), h_vw);
        me_RPA.emplace_back(Fv.shortSymbol() + Fw.shortSymbol(), rpa_vw);
      }
    }
    std::sort(begin(me), end(me), sort_by_first);
    std::sort(begin(me_RPA), end(me_RPA), sort_by_first);

    const auto [eps1, at1] = qip::compare(me, pnc_VD, eps_second);
    const auto [eps2, at2] = qip::compare(me_RPA, pnc_VD_RPA, eps_second);
    pass &=
        qip::check_value(&obuff, "PNC      " + at1->first, eps1, 0.0, 8.0e-5);
    pass &=
        qip::check_value(&obuff, "PNC(RPA) " + at2->first, eps2, 0.0, 2.0e-4);
  }

  return pass;
}
} // namespace UnitTest
