#include "DiracODE/BoundState.hpp"
#include "DiracODE/AsymptoticSpinor.hpp"
#include "LinAlg/include.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace DiracODE {

using namespace Internal;

//==============================================================================
void boundState(DiracSpinor &Fn, const double en0, const std::vector<double> &v,
                const std::vector<double> &H_mag, const double alpha,
                double eps_goal, const DiracSpinor *const VxFa,
                const DiracSpinor *const Fa0, double zion, double mass) {

  assert(Fn.l() < Fn.n() && "Must have valid kappa given n");

  // orbital should have (n-l-1) nodes:
  const int required_nodes = Fn.n() - Fn.l() - 1;

  // Tracks energy guesses + solution state, passed between the steps below.
  TrackEnGuess sofar;
  // start from initial energy guess
  sofar.en = en0;
  // index after last nonzero point
  sofar.pinf = Fn.f().size();
  // g_out - g_in near ctp (for PT)
  sofar.dg.resize(2 * Param::d_ctp + 1);

  // Scratch for the inward solution, allocated once and reused across the
  // energy iterations below (avoids re-allocating two grid-sized vectors each
  // iteration).
  std::vector<double> f_in(Fn.f().size()), g_in(Fn.g().size());

  // Iterate: build a trial solution at the current energy, then adjust the
  // energy (PT step if # nodes correct, large step otherwise) until converged.
  for (; sofar.its < Param::max_its; ++sofar.its) {
    sofar = trialSolution(&Fn, v, H_mag, alpha, VxFa, Fa0, zion, mass, f_in,
                          g_in, sofar);
    sofar = adjustEnergy(&Fn, required_nodes, alpha, sofar);
    if ((sofar.eps < eps_goal && sofar.correct_nodes) || std::isnan(sofar.en))
      break;
  }

  // Extend the decaying tail past pinf and zero the boundaries.
  adjustTail(&Fn, v, H_mag, alpha, zion, mass, sofar);

  // Record results and normalise.
  Fn.en() = sofar.en;
  Fn.eps() = sofar.eps;
  Fn.its() = sofar.its;
  Fn.normalise();

  // Solve failed (nan): zero the orbital and flag as unconverged.
  if (std::isnan(sofar.en)) {
    Fn.en() = 0.0;
    Fn.eps() = 1.0 / 0.0;
    Fn.max_pt() = Fn.grid().num_points();
    Fn.its() = 0;
    Fn *= 0.0;
  }
}

//==============================================================================
void regularAtOrigin(DiracSpinor &Fa, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_mag, const double alpha,
                     const DiracSpinor *const VxFa,
                     const DiracSpinor *const Fa0, double zion, double mass) {

  const auto &gr = Fa.grid();
  if (en != 0.0)
    Fa.en() = en;
  const auto pinf =
    Internal::findPracticalInfinity(Fa.en(), v, gr.r(), Param::cALR, mass);
  Internal::DiracDerivative Hd(gr, v, Fa.kappa(), Fa.en(), alpha, H_mag, VxFa,
                               Fa0, zion, mass);
  Internal::solve_Dirac_outwards(Fa.f(), Fa.g(), Hd, pinf);
  Fa.min_pt() = 0;
  Fa.max_pt() = pinf;
  // for safety: make sure zerod! (I may re-use existing orbitals!)
  Fa.zero_boundaries();
}

//==============================================================================
void regularAtInfinity(DiracSpinor &Fa, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_mag, const double alpha,
                       const DiracSpinor *const VxFa,
                       const DiracSpinor *const Fa0, double zion, double mass) {

  const auto &gr = Fa.grid();
  if (en < 0)
    Fa.en() = en;
  const auto pinf =
    Internal::findPracticalInfinity(Fa.en(), v, gr.r(), Param::cALR, mass);
  Internal::DiracDerivative Hd(gr, v, Fa.kappa(), Fa.en(), alpha, H_mag, VxFa,
                               Fa0, zion, mass);
  Internal::solve_Dirac_inwards(Fa.f(), Fa.g(), Hd, 0, pinf, mass);
  Fa.min_pt() = 0;
  Fa.max_pt() = pinf;
  // for safety: make sure zerod! (I may re-use existing orbitals!)
  Fa.zero_boundaries();
}

//==============================================================================
//==============================================================================
namespace Internal {

//==============================================================================
// Build the trial solution at the current energy (sofar.en): find the practical
// infinity and the matching point near the classical turning point, then
// integrate outwards and inwards and join the two solutions. The asymptotic
// boundary conditions are set inside the integrators (using a previous solution
// Fa0 if supplied, otherwise the standard expansions).
// Returns sofar with pinf, ctp, and dg updated.
TrackEnGuess trialSolution(DiracSpinor *Fn_ptr, const std::vector<double> &v,
                           const std::vector<double> &H_mag, double alpha,
                           const DiracSpinor *VxFa, const DiracSpinor *Fa0,
                           double zion, double mass, std::vector<double> &f_in,
                           std::vector<double> &g_in, TrackEnGuess sofar) {
  auto &Fn = *Fn_ptr;
  const auto &rgrid = Fn.grid();

  sofar.pinf = findPracticalInfinity(sofar.en, v, rgrid.r(), Param::cALR, mass);
  const auto t_ctp =
    findClassicalTurningPoint(sofar.en, v, sofar.pinf, Param::d_ctp);
  sofar.ctp = (1 * sofar.pinf + 4 * t_ctp) / 5;

  trialDiracSolution(Fn.f(), Fn.g(), sofar.dg, f_in, g_in, sofar.en, Fn.kappa(),
                     v, H_mag, rgrid, sofar.ctp, Param::d_ctp, sofar.pinf,
                     alpha, VxFa, Fa0, zion, mass);

  // Keep min/max_pt in sync with the solved region, so the norm (Fn*Fn) used
  // by the PT energy step integrates over exactly [0, pinf). Matters when
  // re-using an orbital whose previous max_pt is smaller than the new pinf.
  Fn.min_pt() = 0;
  Fn.max_pt() = sofar.pinf;
  return sofar;
}

//==============================================================================
// Count the nodes of the current trial solution and update the energy guess:
// correct node count => small perturbation-theory step (mark correct_nodes);
// wrong node count => large (bracketing) step.
// Returns sofar with en, eps, and correct_nodes updated.
TrackEnGuess adjustEnergy(DiracSpinor *Fn_ptr, int required_nodes, double alpha,
                          TrackEnGuess sofar) {
  auto &Fn = *Fn_ptr;

  const int counted_nodes = countNodes(Fn.f(), sofar.pinf);
  const double en_old = sofar.en;
  if (counted_nodes == required_nodes) {
    sofar.correct_nodes = true;
    const double anorm = Fn * Fn;
    sofar.en = smallEnergyChangePT(en_old, anorm, Fn.f(), sofar.dg, sofar.ctp,
                                   Param::d_ctp, alpha, &sofar);
  } else {
    sofar.correct_nodes = false;
    const bool toomany_nodes = (counted_nodes > required_nodes);
    largeEnergyChange(&sofar.en, &sofar, Param::lfrac_de, toomany_nodes);
  }
  sofar.eps = std::abs((sofar.en - en_old) / en_old);
  return sofar;
}

//==============================================================================
// Extend the exponentially-decaying tail past the practical infinity (on the
// converged, not-yet-normalised solution) so the orbital ends at the cutoff,
// and set max_pt. extendTail zeros everything beyond the new pinf. On the rare
// failure where the node count never converged there is no valid tail to
// extend, so we just set pinf and zero the (possibly stale) tail.
void adjustTail(DiracSpinor *Fn_ptr, const std::vector<double> &v,
                const std::vector<double> &H_mag, double alpha, double zion,
                double mass, const TrackEnGuess &sofar) {
  auto &Fn = *Fn_ptr;

  // No valid solution: zero the (possibly stale) tail, since we may be
  // re-using the orbital.
  if (!sofar.correct_nodes || std::isnan(sofar.en)) {
    Fn.max_pt() = sofar.pinf;
    Fn.zero_boundaries();
    return;
  }

  // Extend/trim the decaying tail; extendTail zeros everything beyond new pinf.
  const DiracDerivative Hd_tail(Fn.grid(), v, Fn.kappa(), sofar.en, alpha,
                                H_mag, nullptr, nullptr, zion, mass);
  Fn.max_pt() =
    extendTail(Fn.f(), Fn.g(), Hd_tail, sofar.pinf, Param::tail_cut);
}

//==============================================================================
void largeEnergyChange(double *en, TrackEnGuess *sofar_ptr, double frac_de,
                       bool toomany_nodes) {

  // wf did not have correct number of nodes. Make a large energy adjustment
  // toomany_nodes=true means there were too many nodes

  auto &sofar = *sofar_ptr; // for ease of typing only
  auto etemp = *en;
  if (toomany_nodes) {
    ++sofar.count_toomany;
    if ((sofar.count_toomany == 1) || (*en < sofar.high_en))
      sofar.high_en = *en;
    etemp *= (1.0 + frac_de);
    if ((sofar.count_toofew != 0) && (etemp < sofar.low_en))
      etemp = 0.5 * (sofar.high_en + sofar.low_en);
  } else {
    ++sofar.count_toofew;
    if ((sofar.count_toofew == 1) || (*en > sofar.low_en))
      sofar.low_en = *en;
    etemp *= (1.0 - frac_de);
    if ((sofar.count_toomany != 0) && (etemp > sofar.high_en))
      etemp = 0.5 * (sofar.high_en + sofar.low_en);
  }
  *en = etemp;
}

//==============================================================================
double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, std::size_t ctp,
                           std::size_t d_ctp, const double alpha,
                           TrackEnGuess *sofar_ptr) {
  auto &sofar = *sofar_ptr;

  // delta E = c*f(r)*[g_out(r)-g_in(r)] - evaluate at ctp
  // nb: wf not yet normalised (anorm is input param)!

  assert(ctp + d_ctp < f.size());
  assert(ctp >= d_ctp);
  double p_del_q = f[ctp] * dg[d_ctp];
  double denom = 1.0;
  // weighted average around ctp:
  for (std::size_t i = 1; i <= d_ctp; i++) {
    const auto w = Param::weight(i);
    p_del_q +=
      0.5 * (f[ctp + i] * dg[d_ctp + i] + f[ctp - i] * dg[d_ctp - i]) * w;
    denom += w;
  }

  const double de = p_del_q / (alpha * anorm * denom);

  // Degenerate trial solution (e.g., anorm = 0): PT step is meaningless.
  // Bisect the bracket if we have one; otherwise just nudge the energy
  // (must differ from en, so eps != 0 and we don't falsely "converge").
  if (!std::isfinite(de)) {
    return (sofar.count_toofew != 0 && sofar.count_toomany != 0) ?
             0.5 * (sofar.low_en + sofar.high_en) :
             0.9 * en;
  }

  // Node count is correct, so the sign of de says which side of the
  // eigenvalue en lies on: de > 0 means en is below it (a new lower bound),
  // de < 0 above it (a new upper bound). Tighten the brackets accordingly.
  if (de > 0.0) {
    if (sofar.count_toofew == 0 || en > sofar.low_en)
      sofar.low_en = en;
    ++sofar.count_toofew;
  } else if (de < 0.0) {
    if (sofar.count_toomany == 0 || en < sofar.high_en)
      sofar.high_en = en;
    ++sofar.count_toomany;
  }

  double new_en = en + de;

  if ((sofar.count_toofew != 0) && (new_en < sofar.low_en)) {
    new_en = 0.5 * (en + sofar.low_en);
  } else if ((sofar.count_toomany != 0) && (new_en > sofar.high_en)) {
    new_en = 0.5 * (en + sofar.high_en);
  } else if (new_en > 0) {
    // This only happens v. rarely. nodes correct, but P.T. gives silly result!
    new_en = (de > 0) ? 0.9 * en : 1.1 * en;
  }

  return new_en;
}

//==============================================================================
// Find the practical infinity 'pinf'
// Step backwards from the last point (num_points-1) until
// mass*(V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
// The mass factor accounts for the decay constant ~ sqrt(2*mass*(V-E)):
// without it, pinf for exotic (heavy) particles sits ~sqrt(mass) times too
// many decay lengths out, wasting points and risking underflow in the tail.
std::size_t findPracticalInfinity(const double en, const std::vector<double> &v,
                                  const std::vector<double> &r,
                                  const double alr, const double mass) {
  auto pinf = r.size();
  while (pinf > 10 &&
         mass * (v[pinf - 1] - en) * r[pinf - 1] * r[pinf - 1] > alr) {
    assert(pinf > 1);
    --pinf;
  }
  return pinf;
}

//==============================================================================
std::size_t findClassicalTurningPoint(const double en,
                                      const std::vector<double> &v,
                                      std::size_t pinf, std::size_t d_ctp) {

  // Finds classical turning point 'ctp'
  // Enforced to be between (0+ctp) and (pinf-ctp)
  //  V(r) > E        [nb: both V and E are <0]

  const auto low = std::lower_bound(v.begin() + long(d_ctp + 1),
                                    v.begin() + long(pinf - d_ctp - 1), en);

  const auto distance = std::distance(v.begin(), low);

  return distance > 20 ? std::size_t(distance - 1) : ((pinf - d_ctp + 20) / 2);
}

//==============================================================================
int countNodes(const std::vector<double> &f, std::size_t pinf) {

  // Just counts the number of times orbital (f) changes sign
  int counted_nodes = 0;
  for (std::size_t i = 2; i < pinf; ++i) {
    if (f[i - 1] * f[i] < 0.0) {
      ++counted_nodes;
    }
  }
  return counted_nodes;
}

//==============================================================================
void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, std::vector<double> &f_in,
                        std::vector<double> &g_in, const double en,
                        const int ka, const std::vector<double> &v,
                        const std::vector<double> &H_mag, const Grid &gr,
                        std::size_t ctp, std::size_t d_ctp, std::size_t pinf,
                        const double alpha, const DiracSpinor *const VxFa,
                        const DiracSpinor *const Fa0, double zion,
                        double mass) {

  // Performs inward (from pinf) and outward (from r0) integrations for given
  // energy. Intergated in/out towards ctp +/- d_ctp [class. turn. point]
  // Then, joins solutions, including weighted meshing around ctp +/ d_ctp
  // Also: stores dg [the difference: (gout-gin)], which is used for PT

  DiracDerivative Hd(gr, v, ka, en, alpha, H_mag, VxFa, Fa0, zion, mass);
  solve_Dirac_outwards(f, g, Hd, ctp + d_ctp + 1);
  // f_in/g_in are caller-owned scratch (allocated once per orbital, reused
  // across energy iterations). solve_Dirac_inwards fully (re)writes
  // [ctp-d_ctp, pinf) and zeros beyond, covering everything joinInOutSolutions
  // reads, so no pre-clear is required.
  solve_Dirac_inwards(f_in, g_in, Hd, ctp - d_ctp, pinf, mass);
  joinInOutSolutions(f, g, dg, f_in, g_in, ctp, d_ctp, pinf);
}

//==============================================================================
void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, std::size_t ctp,
                        std::size_t d_ctp, std::size_t pinf) {

  // Find the re-scaling factor (for inward soln)
  // Weighted average of the pointwise ratios f_out/f_in over ctp+/-d_ctp.
  // Skip any point where f_in is zero to avoid division by zero; each side
  // contributes half the weight (matching the symmetric average).
  double rescale = 0.0;
  double denom = 0.0;
  if (f_in[ctp] != 0.0) {
    rescale += f[ctp] / f_in[ctp];
    denom += 1.0;
  }
  for (std::size_t i = 1; i <= d_ctp; i++) {
    const auto w = Param::weight(i);
    if (f_in[ctp + i] != 0.0) {
      rescale += 0.5 * (f[ctp + i] / f_in[ctp + i]) * w;
      denom += 0.5 * w;
    }
    if (f_in[ctp - i] != 0.0) {
      rescale += 0.5 * (f[ctp - i] / f_in[ctp - i]) * w;
      denom += 0.5 * w;
    }
  }
  // Fallback (degenerate: f_in ~ 0 across the whole matching region).
  rescale = (denom > 0.0) ? rescale / denom : 1.0;

  // store difference between in/out solutions (for g) - after re-scaling
  // Used later for P.T.
  for (auto i = 0ul; i < dg.size(); i++) {
    dg[i] = g[ctp - d_ctp + i] - g_in[ctp - d_ctp + i] * rescale;
  }

  // Join the in and outward solutions, meshed smoothly across ctp +/- d_ctp.
  // b is the weight of the inward solution; it runs 0 -> 1 across the mesh
  // region via a smoothstep (3x^2 - 2x^3), which has zero slope at both ends.
  // This makes the join C1-continuous with the pure outward solution below
  // ctp-d_ctp and the pure inward solution above ctp+d_ctp.
  for (std::size_t i = ctp - d_ctp; i <= ctp + d_ctp; i++) {
    const double x =
      double(long(i) - long(ctp) + long(d_ctp)) / double(2 * d_ctp);
    const double b = x * x * (3.0 - 2.0 * x);
    const double a = 1.0 - b;
    f[i] = a * f[i] + b * f_in[i] * rescale;
    g[i] = a * g[i] + b * g_in[i] * rescale;
  }
  for (std::size_t i = ctp + d_ctp + 1; i < pinf; i++) {
    f[i] = f_in[i] * rescale;
    g[i] = g_in[i] * rescale;
  }
}

//==============================================================================
void solve_Dirac_outwards(std::vector<double> &f, std::vector<double> &g,
                          const DiracDerivative &Hd, std::size_t t_pinf) {

  const auto &r = Hd.pgr->r();
  const auto du = Hd.pgr->du();
  const auto &v = *(Hd.v);
  const auto alpha = Hd.alpha;
  const auto ka = Hd.k;
  const auto Z_eff = (-1.0 * v[0] * r[0]);
  const double az0 = Z_eff * alpha;
  const auto ka2 = (double)(ka * ka);
  const double ga0 = std::sqrt(ka2 - az0 * az0);
  const auto pinf = t_pinf == 0 ? f.size() : t_pinf;

  // initial wf values

  // Set initial value:
  double f0{0.0}, g0{0.0};
  if (Hd.Fa0 && Hd.Fa0->f(0) != 0.0) {
    // If we have a previous solution, use that as initial point
    f0 = Hd.Fa0->f(0);
    g0 = Hd.Fa0->g(0);
  } else {
    // Otherwise, use H-like form for f and ratio of f/g
    // f0 is arbitrary, but it's nice to be the correct order-of-magnitude

    // Energy-dependent ratio (exact H-like series); only valid for bound
    // states (en<0), and may fail numerically - fall back to leading-order
    // (r->0) ratio in that case
    const auto gf0 = (ka > 0) ? (ka + ga0) / az0 : az0 / (ka - ga0);
    const auto gf1 =
      Hd.en < 0.0 ?
        DiracHydrogen::gfratio(r[0], ka, Z_eff, alpha, Hd.en, Hd.mass) :
        gf0;
    const auto g_f_ratio = std::isfinite(gf1) ? gf1 : gf0;

    f0 = std::pow(r[0], ga0);
    g0 = f0 * g_f_ratio;
  }

  AdamsMoulton::ODESolver2D<Param::K_Adams, std::size_t, double> ode{du, &Hd};

  ode.solve_initial_K(0, f0, g0);
  for (std::size_t i = 0; i < ode.K_steps(); ++i) {
    f[i] = ode.f[i];
    g[i] = ode.g[i];
  }
  for (std::size_t i = ode.K_steps(); i < pinf; ++i) {
    ode.drive(i);
    f[i] = ode.last_f();
    g[i] = ode.last_g();
  }
}

//==================================================================
void solve_Dirac_inwards(std::vector<double> &f, std::vector<double> &g,
                         const DiracDerivative &Hd, std::size_t nf,
                         std::size_t pinf, double mass) {

  // Program to start the INWARD integration.
  // Starts from Pinf, and uses an expansion [WKB approx] to go to (pinf-K_Adams)
  // i.e., gets last K_Adams points
  // Then, it then call ADAMS-MOULTON, to finish (from num_loops*K_Adams+1
  //   to nf = ctp-d_ctp)

  // short-cuts
  const auto alpha = Hd.alpha;
  const auto ka = Hd.k;
  const auto en = Hd.en;
  const auto &r = Hd.pgr->r();
  const auto &v = *(Hd.v);
  const auto du = Hd.pgr->du();

  const auto Zeff = -v[pinf - 1] * r[pinf - 1];

  static constexpr std::size_t K_AMO = Param::K_Adams;
  AdamsMoulton::ODESolver2D<K_AMO, std::size_t, double> ode{-du, &Hd};

  ode.S_scale = 0.0;

  const auto Rasym = AsymptoticSpinor{ka, Zeff, en, alpha, Param::nx_eps, mass};

  // nb: can use AsymptoticWavefunction for more r values?
  for (std::size_t i0 = pinf - 1, i = 0; i < ode.K_steps(); ++i) {
    const auto [f0, g0] = Rasym.fg(r[i0]);
    ode.f[i] = f0;
    ode.g[i] = g0;
    ode.df[i] = ode.dfdt(f0, g0, i0);
    ode.dg[i] = ode.dgdt(f0, g0, i0);
    ode.t[i] = i0;
    --i0;
  }

  // If we have an existing solution (and an inhomogenous term)
  // we could use existing solution for boundary condition.
  // However, this is numerically unstable, since f(r) is very small.
  // Instead, we re-scale the inhomogenous term so that the normalisation is
  // consistant between the two
  if (Hd.VxFa) {
    const auto i0 = ode.last_t();
    const auto f0 = ode.last_f();
    const auto Xscl =
      (Hd.VxFa && Hd.Fa0->f(i0) != 0.0) ? f0 / Hd.Fa0->f(i0) : 0.0;
    ode.S_scale = Xscl;
    // nb: what if pinf is very different between current and previous soln?
  }

  for (std::size_t i = 0; i < ode.K_steps(); ++i) {
    f[ode.t[i]] = ode.f[i];
    g[ode.t[i]] = ode.g[i];
  }
  const auto i_start = ode.last_t();
  for (std::size_t i = i_start - 1; i >= nf; --i) {
    ode.drive(i);
    f[i] = ode.last_f();
    g[i] = ode.last_g();
    if (i == 0)
      break;
  }
  for (std::size_t i = i_start + ode.K_steps(); i < f.size(); ++i) {
    f[i] = 0.0;
    g[i] = 0.0;
  }
}

//==============================================================================
std::size_t extendTail(std::vector<double> &f, std::vector<double> &g,
                       const DiracDerivative &Hd, std::size_t pinf,
                       double tail_cut) {

  // Adjusts pinf so that |f(pinf-1)| ~ tail_cut*max|f| (see header).
  // If |f| at pinf is still above the cutoff, continues the decaying tail outward
  // (preferring outward Adams-Moulton integration, falling back to the large-r
  // asymptotic exp-decay expansion if that becomes unstable). If |f| at pinf is
  // already below the cutoff, instead reduces pinf back to where |f| ~ cutoff.

  const auto num_points = f.size();
  if (pinf < 2 || pinf > num_points)
    return pinf;

  // Cutoff is relative to the largest |f| over the solved region.
  double fmax = 0.0;
  for (std::size_t i = 0; i < pinf; ++i)
    fmax = std::max(fmax, std::abs(f[i]));
  if (fmax == 0.0)
    return pinf;
  const double f_cut = tail_cut * fmax;
  // only reduce if ~10x smaller than required.
  const double f_cut_reduce = 0.1 * tail_cut * fmax;

  // --- Reduce: pinf already sits below the cutoff; trim the near-zero tail ---
  // Walk inward to the first point (from the outside) with |f| >= f_cut.
  if (std::abs(f[pinf - 1]) < f_cut_reduce) {
    std::size_t new_pinf = pinf;
    while (new_pinf > 1 && std::abs(f[new_pinf - 2]) < f_cut_reduce) {
      --new_pinf;
    }
    for (std::size_t i = new_pinf; i < num_points; ++i) {
      f[i] = 0.0;
      g[i] = 0.0;
    }
    return new_pinf;
  }

  // Extend: continue the decaying tail outward from pinf

  // Nothing to do if already at the grid edge, or too few points to seed the
  // multi-step method.
  if (pinf >= num_points) {
    return pinf;
  }

  // large-r asymptotic (exp-decay) expansion
  // Scaled to match the solved orbital at pinf-1 for a continuous join.
  const auto &r = Hd.pgr->r();
  const auto Zeff = -(*Hd.v)[pinf - 1] * r[pinf - 1];
  const auto Rasym =
    AsymptoticSpinor{Hd.k, Zeff, Hd.en, Hd.alpha, Param::nx_eps, Hd.mass};
  const auto f_a = Rasym.fg(r[pinf - 1]).first;
  const auto scale = (f_a != 0.0) ? f[pinf - 1] / f_a : 0.0;
  std::size_t t_pinf = num_points;
  for (std::size_t i = pinf; i < num_points; ++i) {
    const auto [fi_a, gi_a] = Rasym.fg(r[i]);
    f[i] = scale * fi_a;
    g[i] = scale * gi_a;
    if (std::abs(f[i]) < f_cut) {
      t_pinf = i + 1;
      break;
    }
  }

  // Zero everything beyond the new practical infinity.
  for (std::size_t i = t_pinf; i < num_points; ++i) {
    f[i] = 0.0;
    g[i] = 0.0;
  }

  return t_pinf;
}

//==============================================================================
DiracDerivative::DiracDerivative(
  const Grid &in_grid, const std::vector<double> &in_v, const int in_k,
  const double in_en, const double in_alpha,
  const std::vector<double> &V_off_diag, const DiracSpinor *const iVxFa,
  const DiracSpinor *const iFa0, double izion, double in_mass)
  : pgr(&in_grid),
    v(&in_v),
    Hmag(V_off_diag.empty() ? nullptr : &V_off_diag),
    VxFa(iVxFa),
    Fa0(iFa0),
    zion(izion),
    k(in_k),
    en(in_en),
    alpha(in_alpha),
    cc(1.0 / in_alpha),
    mass(in_mass) {}

double DiracDerivative::a(std::size_t i) const {
  const auto h_mag = (Hmag == nullptr) ? 0.0 : (*Hmag)[i];
  return (double(-k)) * pgr->drduor(i) + alpha * h_mag * pgr->drdu(i);
}
double DiracDerivative::b(std::size_t i) const {
  return (alpha * en + 2.0 * mass * cc - alpha * (*v)[i]) * pgr->drdu(i);
}
double DiracDerivative::c(std::size_t i) const {
  return alpha * ((*v)[i] - en) * pgr->drdu(i);
}
double DiracDerivative::d(std::size_t i) const { return -a(i); }

double DiracDerivative::Sf(std::size_t i) const {
  return VxFa ? -alpha * VxFa->g(i) * pgr->drdu(i) : 0.0;
}
double DiracDerivative::Sg(std::size_t i) const {
  return VxFa ? alpha * VxFa->f(i) * pgr->drdu(i) : 0.0;
}

} // namespace Internal

} // namespace DiracODE
