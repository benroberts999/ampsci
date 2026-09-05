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
std::vector<double>
anderson_coefficients(const std::vector<std::vector<double>> &gram) {
  // Bordered system [B 1; 1^T 0][c; lambda] = [0; 1] over the kept entries.
  // Drops the oldest entries while B is ill-conditioned (saturated history)
  // or the solve fails; there are no exceptions in this codebase, so a
  // failed LU is detected from non-finite coefficients.
  constexpr double cond_floor = 1.0e-14;
  const auto n = gram.size();
  for (auto m = n; m > 0; --m) {
    const auto i0 = n - m; // oldest kept entry
    LinAlg::Matrix<double> B(m + 1, m + 1);
    LinAlg::Vector<double> rhs(m + 1);
    double dmax = 0.0, dmin = 1.0e300;
    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < m; ++j) {
        B(i, j) = gram[i0 + i][i0 + j];
      }
      B(i, m) = 1.0;
      B(m, i) = 1.0;
      rhs(i) = 0.0;
      dmax = std::max(dmax, B(i, i));
      dmin = std::min(dmin, B(i, i));
    }
    B(m, m) = 0.0;
    rhs(m) = 1.0;
    if (m > 2 && dmin < cond_floor * dmax)
      continue;
    const auto c = LinAlg::solve_Axeqb<double>(B, rhs);
    std::vector<double> coefs(m);
    bool finite = true;
    for (std::size_t i = 0; i < m; ++i) {
      coefs[i] = c(i);
      finite = finite && std::isfinite(coefs[i]);
    }
    if (finite)
      return coefs;
  }
  return {};
}

//==============================================================================
std::vector<double>
anderson_coefficients(const std::vector<DiracSpinor> &residuals) {
  const auto n = residuals.size();
  std::vector<std::vector<double>> gram(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      gram[i][j] = residuals[i] * residuals[j];
      gram[j][i] = gram[i][j];
    }
  }
  return anderson_coefficients(gram);
}

//==============================================================================
void solveMixedState_cntm(DiracSpinor &dF, const DiracSpinor &Fa,
                          const double omega, const std::vector<double> &vl,
                          const double alpha,
                          const std::vector<DiracSpinor> &core,
                          const DiracSpinor &hFa, const double eps_target,
                          const HF::Breit *const VBr,
                          const std::vector<double> &H_mag) {
  // As solveMixedState(), with Anderson mixing of the preconditioned
  // fixed-point map in place of damped iteration (see header).
  using namespace qip::overloads;
  assert(dF.kappa() == hFa.kappa());

  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;
  const auto e0 = Fa.en() + omega;

  // Near-singular same-kappa bound states: project the source orthogonal to
  // them, solve the well-conditioned remainder, restore analytically. (Same
  // as solveMixedState(); see conditioning_states().)
  const auto resonant = conditioning_states(core, Fa, dF.kappa(), e0);
  const auto orthogonalise = [&resonant](DiracSpinor &dF_v) {
    for (const auto *Fm : resonant) {
      dF_v.orthog(*Fm);
    }
  };
  auto src = hFa;
  std::vector<double> amps;
  amps.reserve(resonant.size());
  for (const auto *Fm : resonant) {
    const auto cm = (*Fm) * src;
    amps.push_back(cm);
    src -= cm * (*Fm);
  }

  // Seed (or re-seed if dF holds nan from a failed solve at a nearby omega)
  if (std::abs(dF * dF) == 0.0 || !std::isfinite(dF * dF)) {
    DiracODE::solve_inhomog(dF, e0, vl, H_mag, alpha, -1.0 * src);
  }
  orthogonalise(dF);

  // Local exchange on the LHS for conditioning (cancels on the RHS at the
  // fixed point); orbital-independent Kohn-Sham exchange, as solveMixedState()
  const auto Ux = HF::vex_KS(core);
  const auto v = vl + Ux;

  // Preconditioned fixed-point map G(dF) (fixed point = the solution)
  const auto G = [&](const DiracSpinor &x) {
    auto rhs = (Ux * x) - HF::vexFa(x, core) - src;
    if (VBr) {
      rhs -= VBr->VbrFa(x, core);
    }
    auto gx = x; // solve_inhomog uses x as the starting guess
    DiracODE::solve_inhomog(gx, e0, v, H_mag, alpha, rhs);
    orthogonalise(gx);
    return gx;
  };

  constexpr std::size_t history_depth = 8;
  std::vector<DiracSpinor> g_hist; // map outputs g_k = G(x_k)
  std::vector<DiracSpinor> r_hist; // residuals r_k = g_k - x_k

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
    if (g_hist.size() > history_depth) {
      g_hist.erase(g_hist.begin());
      r_hist.erase(r_hist.begin());
    }

    // Anderson extrapolation over the kept history: dF = sum_i c_i g_i
    const auto c = anderson_coefficients(r_hist);
    const auto n_drop = long(r_hist.size() - c.size());
    g_hist.erase(g_hist.begin(), g_hist.begin() + n_drop);
    r_hist.erase(r_hist.begin(), r_hist.begin() + n_drop);
    if (c.empty())
      continue;
    dF = c[0] * g_hist[0];
    for (std::size_t i = 1; i < c.size(); ++i) {
      dF += c[i] * g_hist[i];
    }
    orthogonalise(dF);
  }

  // Restore the projected components analytically: <m|dF> = <m|src>/(e0 - e_m)
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
void solveContinuumMixedState(DiracSpinor *phi, DiracSpinor *Freg,
                              DiracSpinor *Firr, double *K,
                              const DiracSpinor &Fa, const double omega,
                              const std::vector<double> &vl, const double alpha,
                              const std::vector<DiracSpinor> &core,
                              const DiracSpinor &Fs, const double eps_target,
                              const DiracSpinor *const Fhole) {
  using namespace qip::overloads;
  assert(phi != nullptr && Freg != nullptr && Firr != nullptr);
  assert(phi->kappa() == Fs.kappa());

  const auto en_plus = Fa.en() + omega;
  assert(en_plus > 0.0 && "solveContinuumMixedState requires en_+ > 0");

  const auto kappa = phi->kappa();
  const int max_its = (eps_target < 1.0e-8) ? 256 : 128;

  // Orbital-independent Kohn-Sham exchange on the LHS for conditioning
  // (cancels in the source at the fixed point), exactly as the bound
  // solveMixedState(). Being orbital-independent, v is fixed across the
  // iterations, so the homogeneous pair is built only once.
  std::vector<double> Ux(vl.size(), 0.0);
  if (!core.empty()) {
    Ux = HF::vex_KS(core);
  }
  const auto v = vl + Ux;

  // Energy-normalised regular continuum orbital F_reg and its irregular
  // partner F_irr at en_+ in v; reused if the caller already holds the pair
  // at this energy (cached across the outer TDHF iterations at fixed omega)
  const bool reuse = Freg->en() == en_plus && Freg->kappa() == kappa &&
                     Freg->norm2() != 0.0 && Firr->norm2() != 0.0;
  if (!reuse) {
    *Freg = DiracSpinor(0, kappa, Fa.grid_sptr());
    DiracODE::solveContinuum(*Freg, en_plus, v, alpha);
    *Firr = DiracSpinor(0, kappa, Fa.grid_sptr());
    DiracODE::solveContinuumIrregular(*Firr, *Freg, en_plus, v, alpha);
  }

  // Standing-wave particular solution of (h_r^{(v)} - en_+) F = S, by
  // outward integration plus F_reg subtraction. Records the K amplitude of
  // the latest solve.
  double K_last = 0.0;
  const auto inhom_solve = [&](const DiracSpinor &S) {
    DiracSpinor out{0, kappa, Fa.grid_sptr()};
    K_last =
      DiracODE::solveContinuumForward(out, *Freg, *Firr, en_plus, v, alpha, S);
    return out;
  };

  // Norm conservation: project out the diagonal Fa component only, and only
  // when the channel shares Fa's kappa (as the bound solver). All other
  // occupied components are kept: their dV contributions cancel pairwise
  // across channels (<b|X_a> against <a|X_b>), and removing them here while
  // the bound channels keep them would break that cancellation. No occupied
  // state is near-resonant at en_+ > 0.
  const auto orthog_diag = [&](DiracSpinor &F) {
    if (Fa.kappa() == kappa) {
      F -= (F * Fa) * Fa;
    }
  };

  // First pass (no non-local exchange in the source) seeds the iteration,
  // and is the complete answer when there is no exchange. Also re-seeds if
  // phi holds nan from a failed previous solve.
  if (std::abs(*phi * *phi) == 0.0 || !std::isfinite(*phi * *phi)) {
    *phi = inhom_solve(-1.0 * Fs);
    orthog_diag(*phi);
  }
  if (core.empty()) {
    phi->its() = 0;
    phi->eps() = 0.0;
    if (K) {
      *K = K_last;
    }
    return;
  }

  // Iterate the non-local exchange,
  //   (h_r^{(v)} - en_+) phi = Ux*phi - V^exch*phi - Fs,  v = vl + Ux,
  // with Anderson mixing, as solveMixedState_cntm(): the equation is
  // linear, and just above a threshold (especially behind a centrifugal
  // barrier) the damped iteration's multiplier can exceed 1. K is linear in
  // phi, so it mixes with the same coefficients.
  const auto G = [&](const DiracSpinor &x) {
    auto vnl = HF::vexFa(x, core);
    // V^{N-1}: remove one electron's exchange with the hole orbital (the
    // caller pairs this with vl -> vl - y^0_hole,hole for the direct part)
    if (Fhole) {
      vnl -= HF::vexFa_1el(x, *Fhole);
    }
    const auto src = (Ux * x) - vnl - Fs;
    auto gx = inhom_solve(src); // sets K_last
    orthog_diag(gx);
    return gx;
  };

  constexpr std::size_t history_depth = 8;
  std::vector<DiracSpinor> g_hist; // map outputs g_k = G(x_k)
  std::vector<double> k_hist;      // their standing-wave amplitudes
  std::vector<DiracSpinor> r_hist; // residuals r_k = g_k - x_k

  double K_phi = K_last;
  int its = 0;
  double eps = 0.0;
  for (; its < max_its; ++its) {
    const auto g = G(*phi);
    const auto Kg = K_last;
    const auto r = g - *phi;
    eps = std::sqrt(r.norm2() / (phi->norm2() + 1.0e-300));
    if (eps < eps_target) {
      *phi = g;
      K_phi = Kg;
      break;
    }

    g_hist.push_back(g);
    k_hist.push_back(Kg);
    r_hist.push_back(r);
    if (g_hist.size() > history_depth) {
      g_hist.erase(g_hist.begin());
      k_hist.erase(k_hist.begin());
      r_hist.erase(r_hist.begin());
    }

    // Anderson extrapolation over the kept history: phi = sum_i c_i g_i, and K
    // with the same coefficients
    const auto c = anderson_coefficients(r_hist);
    const auto n_drop = long(r_hist.size() - c.size());
    g_hist.erase(g_hist.begin(), g_hist.begin() + n_drop);
    k_hist.erase(k_hist.begin(), k_hist.begin() + n_drop);
    r_hist.erase(r_hist.begin(), r_hist.begin() + n_drop);
    if (c.empty())
      continue;
    *phi = c[0] * g_hist[0];
    K_phi = c[0] * k_hist[0];
    for (std::size_t i = 1; i < c.size(); ++i) {
      *phi += c[i] * g_hist[i];
      K_phi += c[i] * k_hist[i];
    }
    orthog_diag(*phi);
  }
  phi->its() = its;
  phi->eps() = eps;
  if (K) {
    *K = K_phi;
  }
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
