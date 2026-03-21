#include "MixedStates.hpp"
#include "DiracODE/include.hpp"
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

  const auto eta_damp = 0.45;
  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;

  if (std::abs(dF * dF) == 0.0) {
    // If dF is not yet a solution, solve from scratch:
    DiracODE::solve_inhomog(dF, Fa.en() + omega, vl, H_mag, alpha, -1.0 * hFa);
    dF.orthog(Fa);
  }

  // Old/previous value of dF, used for damping
  auto dF0 = dF;

  int its{0};
  double eps{};
  for (; its < max_its; its++) {
    // Include approximate exchange on LHS of equation (therefore also on RHS)
    // This makes v_local have correct long-range behaviour
    const auto vx = HF::vex_approx(dF, core);
    const auto v = vl + vx;
    auto rhs = (vx * dF) - HF::vexFa(dF, core) - hFa;
    if (VBr)
      rhs -= VBr->VbrFa(dF, core);
    if (Sigma)
      rhs -= (*Sigma)(dF);

    DiracODE::solve_inhomog(dF, Fa.en() + omega, v, H_mag, alpha, rhs);

    // damp the solution
    if (its != 0)
      dF = (1.0 - eta_damp) * dF + eta_damp * dF0;

    // Force orthogonality
    dF.orthog(Fa);

    // Check convergence:
    eps = std::sqrt((dF - dF0).norm2() / dF0.norm2());

    if (eps < eps_target) {
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
DiracSpinor solveMixedState_continuum(const DiracSpinor &Fa, double omega,
                                      const std::vector<double> &vl,
                                      double alpha,
                                      const std::vector<DiracSpinor> &core,
                                      const DiracSpinor &hFa,
                                      double eps_target,
                                      const MBPT::CorrelationPotential *const Sigma,
                                      const HF::Breit *const VBr,
                                      const std::vector<double> &H_mag) {
  DiracSpinor dF{0, hFa.kappa(), Fa.grid_sptr()};
  solveMixedState_continuum(dF, Fa, omega, vl, alpha, core, hFa, eps_target,
                            Sigma, VBr, H_mag);
  return dF;
}

//==============================================================================
void solveMixedState_continuum(DiracSpinor &dF, const DiracSpinor &Fa,
                               const double omega,
                               const std::vector<double> &vl,
                               const double alpha,
                               const std::vector<DiracSpinor> &core,
                               const DiracSpinor &hFa,
                               const double eps_target,
                               const MBPT::CorrelationPotential *const Sigma,
                               const HF::Breit *const VBr,
                               const std::vector<double> &H_mag) {
  using namespace qip::overloads;
  assert(dF.kappa() == hFa.kappa());

  const auto eta_damp = 0.45;
  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;

  if (std::abs(dF * dF) == 0.0) {
    // Start from scratch: zeroth-order (no exchange correction)
    DiracODE::solve_inhomog_continuum(dF, Fa.en() + omega, vl, H_mag, alpha,
                                      -1.0 * hFa);
  }

  auto dF0 = dF;

  int its{0};
  double eps{};
  for (; its < max_its; its++) {
    // Include approximate exchange on LHS (therefore also on RHS)
    const auto vx = HF::vex_approx(dF, core);
    const auto v = vl + vx;
    auto rhs = (vx * dF) - HF::vexFa(dF, core) - hFa;
    if (VBr)
      rhs -= VBr->VbrFa(dF, core);
    if (Sigma)
      rhs -= (*Sigma)(dF);

    DiracODE::solve_inhomog_continuum(dF, Fa.en() + omega, v, H_mag, alpha,
                                      rhs);

    // Damp the solution
    if (its != 0)
      dF = (1.0 - eta_damp) * dF + eta_damp * dF0;

    // No orthogonalisation: chi_a is in the continuum, naturally orthogonal
    // to the bound core orbital Fa

    // Convergence: use overlap with source to avoid large-r oscillation
    const auto me_new = dF * hFa;
    const auto me_old = dF0 * hFa;
    eps = (std::abs(me_old) > 0.0) ? std::abs(me_new - me_old) / std::abs(me_old)
                                    : 1.0;

    if (eps < eps_target)
      break;

    dF0 = dF;
  }

  dF.its() = its;
  dF.eps() = eps;

  return;
}

//==============================================================================
DiracSpinor
solveMixedState_continuum(const DiracSpinor &Fa, double omega,
                          const DiracSpinor &hFa,
                          const HF::HartreeFock *const hf, double eps_target,
                          const MBPT::CorrelationPotential *const Sigma) {
  return solveMixedState_continuum(Fa, omega, hf->vlocal(Fa.l()), hf->alpha(),
                                   hf->core(), hFa, eps_target, Sigma,
                                   hf->vBreit(), hf->Hmag(Fa.l()));
}

void solveMixedState_continuum(DiracSpinor &dF, const DiracSpinor &Fa,
                               double omega, const DiracSpinor &hFa,
                               const HF::HartreeFock *const hf,
                               double eps_target,
                               const MBPT::CorrelationPotential *const Sigma) {
  return solveMixedState_continuum(dF, Fa, omega, hf->vlocal(Fa.l()),
                                   hf->alpha(), hf->core(), hFa, eps_target,
                                   Sigma, hf->vBreit(), hf->Hmag(Fa.l()));
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
