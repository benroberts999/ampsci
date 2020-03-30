#include "MixedStates.hpp"
#include "Angular/Angular_369j.hpp"
#include "DiracODE/Adams_Greens.hpp"
#include "DiracODE/DiracODE.hpp"
#include "Coulomb/Coulomb.hpp"
#include "HF/HartreeFockClass.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace HF {

//******************************************************************************
DiracSpinor solveMixedState(const int k, const DiracSpinor &Fa,
                            const double omega, const std::vector<double> &vl,
                            const double alpha,
                            const std::vector<DiracSpinor> &core,
                            const DiracSpinor &hFa, const double eps_target) {
  auto dF = DiracSpinor(0, k, *(Fa.p_rgrid));
  solveMixedState(dF, Fa, omega, vl, alpha, core, hFa, eps_target);
  return dF;
}
//------------------------------------------------------------------------------
DiracSpinor solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa,
                            const double omega, const std::vector<double> &vl,
                            const double alpha,
                            const std::vector<DiracSpinor> &core,
                            const DiracSpinor &hFa, const double eps_target)
// Solves:  (H - e - w)X = -h*Fa for X
{
  auto sp = IO::Profile::safeProfiler(__func__);
  auto damper = rampedDamp(0.8, 0.33, 3, 15);
  const int max_its = 100;

  auto dF20 = std::abs(dF * dF); // monitor convergance
  auto dF0 = dF;

  const std::vector<double> H_mag = {}; // XXX Add Magnetic FF (QED)?
  if (dF20 == 0) {
    DiracODE::solve_inhomog(dF, Fa.en + omega, vl, H_mag, alpha, -1 * hFa);
  } else {
    const auto vx0 = form_approx_vex_any(dF, core);
    const auto v0 = NumCalc::add_vectors(vl, vx0);
    const auto rhs0 = (vx0 * dF) - vex_psia_any(dF, core) - hFa;
    DiracODE::solve_inhomog(dF, Fa.en + omega, v0, H_mag, alpha, rhs0);
    // const auto a = 0.0;
    // const auto l = (1.0 - a);
    // dF = l * dF + a * dF0;
  }

  dF20 = std::abs(dF * dF);
  dF0 = dF;

  for (int its = 0; true; its++) {
    const auto vx = form_approx_vex_any(dF, core);
    const auto v = NumCalc::add_vectors(vl, vx);

    const auto rhs = (vx * dF) - vex_psia_any(dF, core) - hFa;
    DiracODE::solve_inhomog(dF, Fa.en + omega, v, H_mag, alpha, rhs);

    const auto a = damper(its);
    const auto l = (1.0 - a);
    dF = l * dF + a * dF0;
    dF0 = dF;

    auto dF2 = std::abs(dF * dF);
    auto eps = std::abs((dF2 - dF20) / dF2);
    if constexpr (print_each_eps) {
      std::cout << __LINE__ << "| " << Fa.symbol() << " " << its << " " << eps
                << "\n";
    }
    if (eps < eps_target || its == max_its) {
      if constexpr (print_final_eps) {
        std::cout << __LINE__ << "| " << Fa.symbol() << " " << its << " " << eps
                  << "   (<dF|dF>=" << dF2 << ")\n";
        if (its == max_its)
          std::cout << "************\n";
      }
      break;
    }
    dF20 = dF2;
  }
  return dF;
}

} // namespace HF
