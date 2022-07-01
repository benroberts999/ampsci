#include "MixedStates.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/Adams_Greens.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace ExternalField {
using namespace HF;

inline auto rampedDamp(double a_beg, double a_end, int beg, int end) {
  return [=](int i) {
    if (i >= end)
      return a_end;
    if (i <= beg)
      return a_beg;
    return (a_end * (i - beg) + a_beg * (end - i)) / (end - beg);
  };
}

//==============================================================================
DiracSpinor solveMixedState(const int k, const DiracSpinor &Fa,
                            const double omega, const std::vector<double> &vl,
                            const double alpha,
                            const std::vector<DiracSpinor> &core,
                            const DiracSpinor &hFa, const double eps_target,
                            const MBPT::CorrelationPotential *const Sigma,
                            const HF::Breit *const VBr,
                            const std::vector<double> &H_mag) {
  DiracSpinor dF{0, k, Fa.grid_sptr()};
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
      rhs -= VBr->VbrFa(dF, core);
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

//------------------------------------------------------------------------------
void solveMixedState_v2(DiracSpinor &dF, const DiracSpinor &Fa,
                        const double omega, const std::vector<double> &vl,
                        const double alpha,
                        const std::vector<DiracSpinor> &core,
                        const DiracSpinor &hFa, const double eps_target,
                        const MBPT::CorrelationPotential *const Sigma,
                        const HF::Breit *const VBr,
                        const std::vector<double> &H_mag)
// Solves:  (H - e - w)X = -h*Fa for X
{

  /*
  void solve_inhomog(DiracSpinor &Fa, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_mag, const double alpha,
                     const DiracSpinor &source, const DiracSpinor *const VxFa,
                     const DiracSpinor *const Fa0, double zion)
  */

  using namespace qip::overloads;
  DiracODE::solve_inhomog(dF, Fa.en() + omega, vl, H_mag, alpha, -1.0 * hFa);
  for (int its = 0; its < 100; its++) {
    const auto vx = 1.0 * vex_approx(dF, core);
    auto VxdF = (vexFa(dF, core) - vx * dF);
    if (Sigma)
      VxdF += (*Sigma)(dF);
    if (VBr)
      VxdF += VBr->VbrFa(dF, core);
    const auto dF0 = dF;

    DiracODE::solve_inhomog(dF, Fa.en() + omega, vl + vx, H_mag, alpha,
                            -1.0 * hFa, &VxdF, &dF0, 1);

    const auto a = its == 0 ? 0.0 : 0.35;
    if (its != 0)
      dF = (1.0 - a) * dF + a * dF0;
    auto dF2 = std::abs(dF * dF);
    auto dF20 = std::abs(dF0 * dF0);
    auto eps = std::abs((dF2 - dF20) / dF2);
    // std::cout << __LINE__ << "| " << Fa.symbol() << " - " << dF.symbol()
    // << "
    // "
    //           << its << " " << eps << "\n";
    // converges to !1.0e-14!
    if (eps < eps_target || its == 99) {
      break;
    }
  }
  return;
}

} // namespace ExternalField
