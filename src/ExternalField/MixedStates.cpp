#include "MixedStates.hpp"
#include "Angular/Angular_369j.hpp"
#include "Coulomb/Coulomb.hpp"
#include "DiracODE/Adams_Greens.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace ExternalField {
using namespace HF;

//******************************************************************************
DiracSpinor solveMixedState(const int k, const DiracSpinor &Fa,
                            const double omega, const std::vector<double> &vl,
                            const double alpha,
                            const std::vector<DiracSpinor> &core,
                            const DiracSpinor &hFa, const double eps_target,
                            const MBPT::CorrelationPotential *const Sigma,
                            const HF::Breit *const VBr,
                            const std::vector<double> &H_mag) {
  DiracSpinor dF{0, k, Fa.rgrid};
  solveMixedState(dF, Fa, omega, vl, alpha, core, hFa, eps_target, Sigma, VBr,
                  H_mag);
  return dF;
}
//------------------------------------------------------------------------------
void solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa, const double omega,
                     const std::vector<double> &vl, const double alpha,
                     const std::vector<DiracSpinor> &core,
                     const DiracSpinor &hFa, const double eps_target,
                     const MBPT::CorrelationPotential *const Sigma,
                     const HF::Breit *const VBr,
                     const std::vector<double> &H_mag)
// Solves:  (H - e - w)X = -h*Fa for X
{
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  auto damper = rampedDamp(0.8, 0.33, 3, 15);
  const int max_its = eps_target < 1.0e-8 ? 100 : 30;

  if (std::abs(dF * dF) == 0) {
    // If dF is not yet a solution, solve from scratch:
    DiracODE::solve_inhomog(dF, Fa.en() + omega, vl, H_mag, alpha, -1.0 * hFa);
  }

  // monitor convergance:
  auto dF20 = std::abs(dF * dF);
  auto dF0 = dF;

  for (int its = 0; true; its++) {
    const auto vx = vex_approx(dF, core);
    const auto v = qip::add(vl, vx);
    auto rhs = (vx * dF) - vexFa(dF, core) - hFa;
    if (Sigma)
      rhs -= (*Sigma)(dF);
    if (VBr)
      rhs -= (*VBr)(dF);
    DiracODE::solve_inhomog(dF, Fa.en() + omega, v, H_mag, alpha, rhs);

    const auto a = its == 0 ? 0.0 : damper(its);
    dF = (1.0 - a) * dF + a * dF0;
    dF0 = dF;

    auto dF2 = std::abs(dF * dF);
    auto eps = std::abs((dF2 - dF20) / dF2);

    if constexpr (print_each_eps) {
      std::cout << __LINE__ << "| " << Fa.symbol() << " " << its << " " << eps
                << "\n";
    }
    if constexpr (print_final_eps) {
      if (eps < eps_target || its == max_its) {
        std::cout << __LINE__ << "| " << Fa.symbol() << " " << its << " " << eps
                  << "   (<dF|dF>=" << dF2 << ")\n";
      }
    }

    if (eps < eps_target || its == max_its) {
      break;
    }
    dF20 = dF2;
  }
}

} // namespace ExternalField
