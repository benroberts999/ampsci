#include "MixedStatesContinuum.hpp"
#include "DiracODE/include.hpp"
#include "HF/HartreeFock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

namespace ExternalField {

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
  // with Anderson mixing, as solveMixedState_Anderson(): the equation is
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

  // Anderson mixing over the recent iterates, K carried with each entry.
  // The iteration also ends when the residual stops improving: the
  // round-off floor of the solve can lie above a tight target, and the best
  // iterate is then the answer
  constexpr std::size_t history_depth = 8;
  AndersonHistory history(history_depth);
  double K_phi = K_last;
  auto best = *phi;
  double best_K = K_phi;
  double best_eps = 1.0e30;
  int count_worse = 0;

  int its = 0;
  double eps = 0.0;
  for (; its < max_its; ++its) {
    auto g = G(*phi);
    const auto Kg = K_last;
    auto r = g - *phi;
    const auto r2 = r.norm2();
    eps = std::sqrt(r2 / (phi->norm2() + 1.0e-300));
    if (eps < eps_target) {
      *phi = std::move(g);
      K_phi = Kg;
      break;
    }
    if (eps < 0.99 * best_eps) {
      best_eps = eps;
      best = g;
      best_K = Kg;
      count_worse = 0;
    } else if (++count_worse > 5) {
      *phi = std::move(best);
      K_phi = best_K;
      eps = best_eps;
      break;
    }

    // Anderson extrapolation over the kept history: phi = sum_i c_i g_i, and
    // K with the same coefficients
    history.push(std::move(g), std::move(r), r2, Kg);
    history.mix(phi, &K_phi);
    orthog_diag(*phi);
  }
  phi->its() = its;
  phi->eps() = eps;
  if (K) {
    *K = K_phi;
  }
}

} // namespace ExternalField
