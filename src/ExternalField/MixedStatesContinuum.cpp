#include "MixedStatesContinuum.hpp"
#include "DiracODE/include.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "LinAlg/Solvers.hpp"
#include "LinAlg/Vector.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

namespace ExternalField {

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
void AndersonHistory::push(DiracSpinor g, DiracSpinor r, double r2, double k) {
  if (m_g.size() >= m_depth) {
    drop_oldest(m_g.size() + 1 - m_depth);
  }
  // New row and column of the Gram matrix: one inner product per kept entry
  const auto n = m_r.size();
  std::vector<double> row(n + 1);
  for (std::size_t i = 0; i < n; ++i) {
    row[i] = m_r[i] * r;
    m_gram[i].push_back(row[i]);
  }
  row[n] = r2;
  m_gram.push_back(std::move(row));
  m_g.push_back(std::move(g));
  m_r.push_back(std::move(r));
  m_k.push_back(k);
}

//------------------------------------------------------------------------------
void AndersonHistory::drop_oldest(std::size_t n) {
  n = std::min(n, m_g.size());
  if (n == 0)
    return;
  const auto ln = long(n);
  m_g.erase(m_g.begin(), m_g.begin() + ln);
  m_r.erase(m_r.begin(), m_r.begin() + ln);
  m_k.erase(m_k.begin(), m_k.begin() + ln);
  m_gram.erase(m_gram.begin(), m_gram.begin() + ln);
  for (auto &row : m_gram) {
    row.erase(row.begin(), row.begin() + ln);
  }
}

//------------------------------------------------------------------------------
void AndersonHistory::mix(DiracSpinor *x, double *k) {
  assert(x != nullptr && !m_g.empty());
  const auto c = anderson_coefficients(m_gram);
  if (c.empty()) {
    // No usable entry: the plain step, and the history starts afresh
    *x = std::move(m_g.back());
    if (k) {
      *k = m_k.back();
    }
    drop_oldest(m_g.size());
    return;
  }
  drop_oldest(m_g.size() - c.size());
  *x = m_g[0];
  *x *= c[0];
  double k_mixed = c[0] * m_k[0];
  for (std::size_t i = 1; i < c.size(); ++i) {
    x->add_scaled(c[i], m_g[i]);
    k_mixed += c[i] * m_k[i];
  }
  if (k) {
    *k = k_mixed;
  }
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

  // Anderson mixing over the recent iterates. The iteration also ends when
  // the residual stops improving: the round-off floor of the solve can lie
  // above a tight target, and the best iterate is then the answer
  constexpr std::size_t history_depth = 8;
  AndersonHistory history(history_depth);
  auto best = dF;
  double best_eps = 1.0e30;
  int count_worse = 0;

  int its{0};
  double eps{};
  for (; its < max_its; its++) {
    auto g = G(dF);
    auto r = g - dF;
    const auto r2 = r.norm2();
    eps = std::sqrt(r2 / (dF.norm2() + 1.0e-300));
    if (eps < eps_target) {
      dF = std::move(g);
      break;
    }
    if (eps < 0.99 * best_eps) {
      best_eps = eps;
      best = g;
      count_worse = 0;
    } else if (++count_worse > 5) {
      dF = std::move(best);
      eps = best_eps;
      break;
    }

    // Anderson extrapolation over the kept history: dF = sum_i c_i g_i
    history.push(std::move(g), std::move(r), r2);
    history.mix(&dF);
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
