#include "MixedStates.hpp"
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

//==============================================================================
DiracSpinor solveMixedState(const DiracSpinor &Fa, double omega,
                            const std::vector<double> &vl, double alpha,
                            const std::vector<DiracSpinor> &core,
                            const DiracSpinor &hFa, double eps_target,
                            const MBPT::CorrelationPotential *const Sigma,
                            const HF::Breit *const VBr,
                            const std::vector<double> &H_mag) {
  DiracSpinor dF{0, hFa.kappa(), Fa.grid_sptr()};
  solveMixedState(dF, Fa, omega, vl, alpha, core, hFa, eps_target, Sigma, VBr,
                  H_mag);
  return dF;
}

//==============================================================================
void solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa, const double omega,
                     const std::vector<double> &vl, const double alpha,
                     const std::vector<DiracSpinor> &core,
                     const DiracSpinor &hFa, const double eps_target,
                     const MBPT::CorrelationPotential *const Sigma,
                     const HF::Breit *const VBr,
                     const std::vector<double> &H_mag) {
  using namespace qip::overloads;
  assert(dF.kappa() == hFa.kappa());

  const auto eta_damp = 0.5;
  const int max_its = (eps_target < 1.0e-7) ? 100 : 30;

  if (std::abs(dF * dF) == 0.0) {
    // If dF is not yet a solution, solve from scratch:
    DiracODE::solve_inhomog(dF, Fa.en() + omega, vl, H_mag, alpha, -1.0 * hFa);
    dF.orthog(Fa);
  }

  // Old/previous value of dF, used for damping
  auto dF0 = dF;

  int its{0};
  double eps{};
  for (;; its++) {
    // Include approximate exchange on LHS of equation (therefore also on RHS)
    // This makes v_local have correct long-range behaviour
    const auto vx = HF::vex_approx(dF, core);
    const auto v = vl + vx;
    auto rhs = (vx * dF) - HF::vexFa(dF, core) - hFa;
    if (Sigma)
      rhs -= (*Sigma)(dF);
    if (VBr)
      rhs -= VBr->VbrFa(dF, core);

    DiracODE::solve_inhomog(dF, Fa.en() + omega, v, H_mag, alpha, rhs);

    // damp the solution
    if (its != 0)
      dF = (1.0 - eta_damp) * dF + eta_damp * dF0;

    // Force orthogonality
    dF.orthog(Fa);

    // Check convergence:
    const auto dF2 = dF * dF;
    eps = (dF - dF0).norm2() / dF2;

    if (eps < eps_target || its == max_its) {
      break;
    }

    // store old values (for damping and convergence)
    dF0 = dF;
  }

  dF.its() = its;
  dF.eps() = eps;

  return;
}

//==============================================================================
DiracSpinor solveMixedState(const DiracSpinor &Fa, double omega,
                            const DiracSpinor &hFa,
                            const HF::HartreeFock *const hf, double eps_target,
                            const MBPT::CorrelationPotential *const Sigma) {
  return solveMixedState(Fa, omega, hf->vlocal(Fa.l()), hf->alpha(), hf->core(),
                         hFa, eps_target, Sigma, hf->vBreit(),
                         hf->Hmag(Fa.l()));
}

void solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa, double omega,
                     const DiracSpinor &hFa, const HF::HartreeFock *const hf,
                     double eps_target,
                     const MBPT::CorrelationPotential *const Sigma) {
  return solveMixedState(dF, Fa, omega, hf->vlocal(Fa.l()), hf->alpha(),
                         hf->core(), hFa, eps_target, Sigma, hf->vBreit(),
                         hf->Hmag(Fa.l()));
}

//==============================================================================
DiracSpinor solveMixedState_basis(const DiracSpinor &Fa, const DiracSpinor &hFa,
                                  double omega,
                                  const std::vector<DiracSpinor> &basis) {
  DiracSpinor dFa = 0.0 * hFa;

  for (const auto &n : basis) {
    if (n == Fa || n.kappa() != hFa.kappa())
      continue;
    dFa += ((n * hFa) / (Fa.en() - n.en() + omega)) * n;
  }
  return dFa;
}

} // namespace ExternalField
