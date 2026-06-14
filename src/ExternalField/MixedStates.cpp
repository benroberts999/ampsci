#include "MixedStates.hpp"
#include "DiracODE/include.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace ExternalField {

//==============================================================================
// Bound states of the solve channel that make (h_l - e0) near-singular.
std::vector<const DiracSpinor *>
conditioning_states(const std::vector<DiracSpinor> &core, const DiracSpinor &Fa,
                    int kappa, double e0) {
  std::vector<const DiracSpinor *> states;
  for (const auto &Fm : core) {
    if (Fm.kappa() != kappa || Fm == Fa)
      continue;
    // relative nearness: e_m within ~20% of e0. Catches fine-structure partners
    // (which must be conditioned) while leaving well-separated states to be
    // found naturally by the solve.
    if (std::abs(e0 - Fm.en()) < 0.2 * std::abs(e0 + Fm.en())) {
      states.push_back(&Fm);
    }
  }
  // Fa only when its kappa matches the channel (else the projection/orthog is
  // a cross-kappa no-op, but the inner product is not guarded -- so skip it).
  if (Fa.kappa() == kappa) {
    states.push_back(&Fa);
  }
  return states;
}

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

  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;
  const auto e0 = Fa.en() + omega;

  // (h_l - e0) is (near-)singular for components along same-kappa bound states
  // with e_m ~ e0 (fine-structure partners, and the diagonal Fa). We project
  // the source orthogonal to these, solve for the well-conditioned remainder
  // while forcing it orthogonal to them, then restore the off-diagonal
  // components analytically via first-order PT. See conditioning_states().
  const auto resonant = conditioning_states(core, Fa, dF.kappa(), e0);

  const auto orthogonalise = [&resonant](DiracSpinor &dF_v) {
    for (const auto *Fm : resonant) {
      dF_v.orthog(*Fm);
    }
  };

  // Project the source orthogonal to the conditioning set; keep the amplitudes
  // <m|src> for the analytic restoration below.
  auto src = hFa;
  std::vector<double> amps;
  amps.reserve(resonant.size());
  for (const auto *Fm : resonant) {
    const auto cm = (*Fm) * src;
    amps.push_back(cm);
    src -= cm * (*Fm);
  }

  if (std::abs(dF * dF) == 0.0) {
    // If dF is not yet a solution, solve from scratch:
    DiracODE::solve_inhomog(dF, e0, vl, H_mag, alpha, -1.0 * src);
  }
  orthogonalise(dF);

  const auto eta_damp = 0.85;

  // Local exchange on the LHS of the equation (also added to the RHS, so it
  // cancels at the fixed point -- it only conditions the iteration). We use the
  // density-based Kohn-Sham/Slater local exchange: being orbital-independent it
  // conditions all channels uniformly and is computed once. (The alternative,
  // vex_approx, divides by the perturbation and is poorly conditioned for
  // awkward channels.)
  const auto Ux = HF::vex_KS(core);
  const auto v = vl + Ux;

  // Old/previous value of dF, used for damping
  auto dF0 = dF;

  int its{0};
  double eps{};
  for (; its < max_its; its++) {
    auto rhs = (Ux * dF) - HF::vexFa(dF, core) - src;
    if (VBr) {
      rhs -= VBr->VbrFa(dF, core);
    }
    if (Sigma) {
      rhs -= (*Sigma)(dF);
    }

    DiracODE::solve_inhomog(dF, e0, v, H_mag, alpha, rhs);

    // damp the solution
    if (its != 0) {
      dF = (1.0 - eta_damp) * dF + eta_damp * dF0;
    }

    // Force orthogonality
    orthogonalise(dF);

    // Check convergence:
    eps = std::sqrt((dF - dF0).norm2() / dF0.norm2());

    if (eps < eps_target) {
      break;
    }

    // store old values (for damping and convergence)
    dF0 = dF;
  }

  // Restore the projected components analytically: <m|dF> = <m|src>/(e0 - e_m).
  // Skip Fa (the left-orthogonality constraint -- left orthogonal) and any
  // genuinely resonant denominator (Pauli-blocked -- left orthogonal).
  for (std::size_t i = 0; i < resonant.size(); ++i) {
    const auto *Fm = resonant[i];
    if (*Fm == Fa)
      continue;
    const auto denom = e0 - Fm->en();
    // This might need to be checked for some systems.. probably fine
    if (std::abs(denom) < 1.0e-5)
      continue;
    dF += (amps[i] / denom) * (*Fm);
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
