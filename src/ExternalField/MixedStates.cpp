#include "MixedStates.hpp"
#include "DiracODE/include.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "LinAlg/Solvers.hpp"
#include "LinAlg/Vector.hpp"
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
void solveMixedState_cntm(DiracSpinor &dF, const DiracSpinor &Fa,
                          const double omega, const std::vector<double> &vl,
                          const double alpha,
                          const std::vector<DiracSpinor> &core,
                          const DiracSpinor &hFa, const double eps_target,
                          const HF::Breit *const VBr,
                          const std::vector<double> &H_mag,
                          const DiracSpinor *const Fhole) {
  // Bound mixed-state solve for the continuum (photoionisation) TDHF path. Same
  // physics as solveMixedState(), but solves the LINEAR equation
  // (h_HF - e0)dF = -src with Anderson/Pulay acceleration instead of damped
  // fixed-point iteration. At the high frequencies needed for photoionisation a
  // bound channel's e0 can sit near a SPURIOUS eigenvalue of the local
  // preconditioner (h_local + Ux - e0), where the damped iteration diverges even
  // though (h_HF - e0) is non-singular and the response is finite. Anderson
  // mixing (= preconditioned GMRES for this linear system) converges on the
  // conditioning of the TRUE operator, independent of the preconditioner's
  // spurious modes. Kept separate so the standard bound solveMixedState() (and
  // the bound RPA path) is untouched. Also supports the V^{N-1} hole term via
  // @p Fhole (see solveMixedState() note).
  using namespace qip::overloads;
  assert(dF.kappa() == hFa.kappa());

  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;
  const auto e0 = Fa.en() + omega;

  // Near-singular same-kappa bound states: project the source orthogonal to
  // them, solve the well-conditioned remainder, restore analytically. (Same as
  // solveMixedState(); see conditioning_states().)
  const auto resonant = conditioning_states(core, Fa, dF.kappa(), e0);
  const auto orthogonalise = [&resonant](DiracSpinor &dF_v) {
    for (const auto *Fm : resonant)
      dF_v.orthog(*Fm);
  };
  auto src = hFa;
  std::vector<double> amps;
  amps.reserve(resonant.size());
  for (const auto *Fm : resonant) {
    const auto cm = (*Fm) * src;
    amps.push_back(cm);
    src -= cm * (*Fm);
  }

  // Seed (or re-seed if dF holds nan from a failed solve at a nearby omega).
  if (std::abs(dF * dF) == 0.0 || !std::isfinite(dF * dF)) {
    DiracODE::solve_inhomog(dF, e0, vl, H_mag, alpha, -1.0 * src);
  }
  orthogonalise(dF);

  // Local exchange on the LHS for conditioning (cancels on the RHS at the fixed
  // point); orbital-independent Kohn-Sham/Slater exchange, as solveMixedState().
  const auto Ux = HF::vex_KS(core);
  const auto v = vl + Ux;

  // Preconditioned fixed-point map G(dF) (fixed point = the solution).
  const auto G = [&](const DiracSpinor &x) {
    auto vnl = HF::vexFa(x, core);
    // V^{N-1}: remove one electron's exchange with the hole orbital
    if (Fhole)
      vnl -= HF::vexFa_1el(x, *Fhole);
    auto rhs = (Ux * x) - vnl - src;
    if (VBr)
      rhs -= VBr->VbrFa(x, core);
    auto gx = x; // solve_inhomog uses x as the starting guess
    DiracODE::solve_inhomog(gx, e0, v, H_mag, alpha, rhs);
    orthogonalise(gx);
    return gx;
  };

  constexpr std::size_t and_dim = 8;     // Anderson history depth
  constexpr double cond_floor = 1.0e-14; // drop-oldest threshold on B
  std::vector<DiracSpinor> g_hist;       // map outputs g_k = G(x_k)
  std::vector<DiracSpinor> r_hist;       // residuals r_k = g_k - x_k

  int its{0};
  double eps{};
  for (; its < max_its; its++) {
    const auto g = G(dF);
    const auto r = g - dF;
    eps = std::sqrt(r.norm2() / (dF.norm2() + 1.0e-300));
    if (eps < eps_target) {
      dF = g;
      break;
    }

    g_hist.push_back(g);
    r_hist.push_back(r);
    if (g_hist.size() > and_dim) {
      g_hist.erase(g_hist.begin());
      r_hist.erase(r_hist.begin());
    }

    // First step: nothing to mix yet -- plain fixed-point step.
    auto m = g_hist.size();
    if (m == 1) {
      dF = g;
      continue;
    }

    // Anderson coefficients: minimise || sum_i c_i r_i ||^2 s.t. sum_i c_i = 1,
    // via the bordered system [B 1; 1' 0][c; lam] = [0; 1], B_ij = <r_i|r_j>.
    // If B is ill-conditioned (saturated history), drop the oldest and re-build.
    LinAlg::Vector<double> c;
    while (true) {
      LinAlg::Matrix<double> B(m + 1, m + 1);
      LinAlg::Vector<double> bvec(m + 1);
      double dmax = 0.0, dmin = 1.0e300;
      for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j)
          B(i, j) = r_hist[i] * r_hist[j];
        B(i, m) = 1.0;
        B(m, i) = 1.0;
        bvec(i) = 0.0;
        dmax = std::max(dmax, B(i, i));
        dmin = std::min(dmin, B(i, i));
      }
      B(m, m) = 0.0;
      bvec(m) = 1.0;
      if (m > 2 && dmin < cond_floor * dmax) {
        g_hist.erase(g_hist.begin());
        r_hist.erase(r_hist.begin());
        --m;
        continue;
      }
      c = LinAlg::solve_Axeqb<double>(B, bvec);
      // No exceptions in this codebase: guard a singular/near-singular LU by
      // checking the coefficients are finite; if not, drop the oldest entry and
      // re-build. (The ill-conditioning check above normally prevents this.)
      bool finite = true;
      for (std::size_t i = 0; i < m; ++i) {
        if (!std::isfinite(c(i)))
          finite = false;
      }
      if (finite)
        break;
      g_hist.erase(g_hist.begin());
      r_hist.erase(r_hist.begin());
      --m;
      if (m == 0)
        break;
    }
    if (m == 0)
      continue;

    // New iterate: dF = sum_i c_i g_i (DIIS extrapolation).
    dF = c(0) * g_hist[0];
    for (std::size_t i = 1; i < m; ++i)
      dF += c(i) * g_hist[i];
    orthogonalise(dF);
  }

  // Restore the projected components analytically: <m|dF> = <m|src>/(e0 - e_m).
  for (std::size_t i = 0; i < resonant.size(); ++i) {
    const auto *Fm = resonant[i];
    if (*Fm == Fa)
      continue;
    const auto denom = e0 - Fm->en();
    if (std::abs(denom) < 1.0e-5)
      continue;
    dF += (amps[i] / denom) * (*Fm);
  }

  dF.its() = its;
  dF.eps() = eps;
}

//==============================================================================
void solveContinuumMixedState(DiracSpinor &phi, DiracSpinor &Freg,
                              DiracSpinor &Firr, const DiracSpinor &Fa,
                              double omega, const std::vector<double> &vl,
                              double alpha, const DiracSpinor &Fs) {
  using namespace qip::overloads;
  assert(phi.kappa() == Fs.kappa());

  const auto en_plus = Fa.en() + omega;
  assert(en_plus > 0.0 && "solveContinuumMixedState requires en_+ > 0");

  const auto kappa_alpha = phi.kappa();

  // Regular (energy-normalised) homogeneous solution at en_+, and its
  // irregular partner, in the local potential vl.
  Freg = DiracSpinor(0, kappa_alpha, Fa.grid_sptr());
  DiracODE::solveContinuum(Freg, en_plus, vl, alpha);
  Firr = DiracSpinor(0, kappa_alpha, Fa.grid_sptr());
  DiracODE::solveContinuumIrregular(Firr, Freg, en_plus, vl, alpha);

  // Particular solution of (h_r - en_+) phi = -Fs with the standing-wave BC,
  // by forward (outward) integration + F_reg subtraction. (Not the global
  // Green's function: F_irr diverges at the origin like r^{-l-1}, which the
  // outward-integrated solution avoids -- F_irr is only sampled in the outer
  // window for the c,K extraction.)
  DiracODE::solveContinuumForward(phi, Freg, Firr, en_plus, vl, alpha,
                                  -1.0 * Fs);
}

//==============================================================================
void solveContinuumMixedState(DiracSpinor &phi, DiracSpinor &Freg,
                              DiracSpinor &Firr, const DiracSpinor &Fa,
                              const double omega, const std::vector<double> &vl,
                              const double alpha,
                              const std::vector<DiracSpinor> &core,
                              const DiracSpinor &Fs, const double eps_target,
                              const DiracSpinor *const Fhole) {
  using namespace qip::overloads;
  assert(phi.kappa() == Fs.kappa());

  const auto en_plus = Fa.en() + omega;
  assert(en_plus > 0.0 && "solveContinuumMixedState requires en_+ > 0");

  const auto kappa = phi.kappa();
  const auto eta_damp = 0.85;
  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;

  // Local exchange on the LHS for conditioning (added to v, cancels on the RHS
  // at the fixed point). Orbital-independent density-based Kohn-Sham/Slater
  // exchange, exactly as in the bound solveMixedState(): the alternative
  // vex_approx divides by the (oscillatory) continuum orbital and spikes at its
  // nodes, which destabilises the iteration. Being orbital-independent, Ux (and
  // hence v) is fixed across iterations, so F_reg/F_irr are built only once
  // (which also removes the per-rebuild F_reg renormalisation noise).
  std::vector<double> Ux(vl.size(), 0.0);
  if (!core.empty())
    Ux = HF::vex_KS(core);
  const auto v = vl + Ux;

  // Energy-normalised regular continuum orbital F_reg and its irregular partner
  // F_irr at en_+ in the (fixed) local potential v.
  Freg = DiracSpinor(0, kappa, Fa.grid_sptr());
  DiracODE::solveContinuum(Freg, en_plus, v, alpha);
  Firr = DiracSpinor(0, kappa, Fa.grid_sptr());
  DiracODE::solveContinuumIrregular(Firr, Freg, en_plus, v, alpha);

  // Standing-wave particular solution of (h_r^{(v)} - en_+) F = S, by forward
  // (outward) integration + F_reg subtraction (regular at the origin; F_irr is
  // only sampled in the outer window, avoiding its r^{-l-1} divergence).
  const auto inhom_solve = [&](const DiracSpinor &S) {
    DiracSpinor out{0, kappa, Fa.grid_sptr()};
    DiracODE::solveContinuumForward(out, Freg, Firr, en_plus, v, alpha, S);
    return out;
  };

  // Project out occupied core orbitals of the same kappa (continuum
  // normalisation unaffected: the overlap with the decaying bound orbital is
  // finite). Replaces the bound solver's single dF.orthog(Fa).
  const auto orthog_to_core = [&](DiracSpinor &F) {
    for (const auto &Fb : core) {
      if (Fb.kappa() == kappa)
        F -= (F * Fb) * Fb;
    }
  };

  // First pass (no nonlocal exchange in the source): seeds the iteration, and
  // is the complete answer when there is no exchange.
  // (also re-seed if phi holds nan from a failed previous solve)
  if (std::abs(phi * phi) == 0.0 || !std::isfinite(phi * phi)) {
    phi = inhom_solve(-1.0 * Fs);
    orthog_to_core(phi);
  }
  if (core.empty()) {
    phi.its() = 0;
    phi.eps() = 0.0;
    return;
  }

  // Iterate the nonlocal exchange, mirroring the bound solveMixedState():
  //   (h_r^{(v)} - en_+) phi = Ux*phi - V^exch*phi - Fs,  v = vl + Ux.
  auto phi0 = phi;
  int its = 0;
  double eps = 0.0;
  for (; its < max_its; ++its) {
    auto vnl = HF::vexFa(phi, core);
    // V^{N-1}: remove one electron's exchange with the hole orbital (the caller
    // pairs this with vl -> vl - y^0_hole,hole for the direct part)
    if (Fhole)
      vnl -= HF::vexFa_1el(phi, *Fhole);
    const auto src = (Ux * phi) - vnl - Fs;
    phi = inhom_solve(src);

    if (its != 0)
      phi = (1.0 - eta_damp) * phi + eta_damp * phi0;
    orthog_to_core(phi);

    eps = std::sqrt((phi - phi0).norm2() / phi0.norm2());
    if (eps < eps_target)
      break;
    phi0 = phi;
  }
  phi.its() = its;
  phi.eps() = eps;
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
