#pragma once
#include "DiracOperator/Operators.hpp"
#include "Physics/FGRadPot.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>
#include <utility>
#include <vector>

bool FGRadPot(std::ostream &obuff) {
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
