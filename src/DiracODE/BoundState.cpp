#include "DiracODE/BoundState.hpp"
#include "DiracODE/AsymptoticSpinor.hpp"
#include "LinAlg/include.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
//
#include "Physics/DiracHydrogen.hpp"

namespace DiracODE {

using namespace Internal;

//==============================================================================
void boundState(DiracSpinor &Fn, const double en0, const std::vector<double> &v,
                const std::vector<double> &H_mag, const double alpha,
                double eps_goal, const DiracSpinor *const VxFa,
                const DiracSpinor *const Fa0, double zion, double mass)
/*
Solves local, spherical bound state dirac equation using Adams-Moulton
method. Based on method presented in book by W. R. Johnson:
  W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007)
I have added a few extensions to this method. In particular, I integrate past
the classical turning point. (See below)

Rough description of method:
1. Start with initial 'guess' of energy
2. Find "practical infinity" (Fn~0), and Classical turning point (e=v)
3. Performs 'inward' integration (Adams Moulton).
   Integrates inwards from the practical infinity inwards to d_ctp points
   past the classical turning point (ctp).
4. Performs 'outward' integration (Adams Moulton). Integrates from 0
   outwards to d_ctp points past the ctp.
5. Matches the two functions around ctp, by re-scaling the 'inward' solution
   (for f). Uses a weighted mean, with weights given by distance from ctp.
6. Checks the number of nodes the wf has. If too many or too few nodes, makes
   a large change to the energy and tries again (from step 2).
   If the correct number of nodes, uses perturbation theory to make minor
   corrections to the energy to 'zoom in' (matching the in/out solution for
   g), then re-starts from step 2.
Continues until this energy adjustment falls below a prescribed threshold.

Orbitals defined:
  Fn := (1/r) {f O_k, ig O_(-k)}
*/
{

  assert(Fn.l() < Fn.n() && "Must have valid kappa given n");

  const auto &rgrid = Fn.grid();

  // orbital should have (n-l-1) nodes:
  const int required_nodes = Fn.n() - Fn.l() - 1;
  bool correct_nodes = false;
  TrackEnGuess sofar; // track higest/lowest energy guesses etc.

  // Energy values; starts from initial guess
  double t_en = en0;
  // value for practical infinity
  // nb: pinf is index *after* final nonz-zero; must not eval array at pinf
  std::size_t t_pinf = Fn.f().size();
  // epsilon: convergance parameter
  double t_eps = 1.0;
  // Normalisation constant; kept for eficiancy
  double anorm = 1.0;

  // Iterations
  int t_its = 1;
  for (; t_its < Param::max_its; ++t_its) {

    t_pinf = Internal::findPracticalInfinity(t_en, v, rgrid.r(), Param::cALR);

    const std::size_t t_ctp =
        Internal::findClassicalTurningPoint(t_en, v, t_pinf, Param::d_ctp);
    const auto ctp = (1 * t_pinf + 4 * t_ctp) / 5;

    // Find solution (f,g) to DE for given energy:
    // Also stores dg (gout-gin) for PT [used for PT to find better e]
    std::vector<double> dg(2 * Param::d_ctp + 1);

    Internal::trialDiracSolution(Fn.f(), Fn.g(), dg, t_en, Fn.kappa(), v, H_mag,
                                 rgrid, ctp, Param::d_ctp, t_pinf, alpha, VxFa,
                                 Fa0, zion, mass);

    const int counted_nodes = Internal::countNodes(Fn.f(), t_pinf);

    // If correct number of nodes, use PT to make minor energy adjustment.
    // Otherwise, make large adjustmunt until correct # of nodes
    const double en_old = t_en;
    if (counted_nodes == required_nodes) {
      correct_nodes = true;
      anorm = Fn * Fn;
      t_en = Internal::smallEnergyChangePT(en_old, anorm, Fn.f(), dg, ctp,
                                           Param::d_ctp, alpha, sofar);

    } else {
      correct_nodes = false;
      const bool toomany_nodes = (counted_nodes > required_nodes);
      Internal::largeEnergyChange(&t_en, &sofar, Param::lfrac_de,
                                  toomany_nodes);
    }
    t_eps = std::abs((t_en - en_old) / en_old);

    if ((t_eps < eps_goal && correct_nodes) || std::isnan(t_en))
      break;

  } // END itterations

  // If we never got correct # of nodes, never calc'd norm constant.
  // This is rare - means a failure. Occurs when energy guess is too wrong.
  // It sometimes occurs on the first few HF iterations. Usually, failure will
  // go away after a few more HF iterations.
  // But if we don't normalise the wf, HF will fail.
  if (!correct_nodes) {
    anorm = Fn * Fn;
  }

  // store energy etc.
  Fn.en() = t_en;
  Fn.eps() = t_eps;
  Fn.max_pt() = t_pinf;
  Fn.its() = t_its;

  // Explicitely set 'tail' to zero (we may be re-using orbital)
  Fn.zero_boundaries();
  // normalises the orbital (alrady cal'd anorm)
  Fn *= 1.0 / std::sqrt(anorm);

  if (std::isnan(t_en)) {
    Fn.en() = 0.0;
    Fn.eps() = 1.0 / 0.0;
    Fn.max_pt() = rgrid.num_points();
    Fn.its() = 0;
    Fn *= 0.0;
  }

  return;
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
      Internal::findPracticalInfinity(Fa.en(), v, gr.r(), Param::cALR);
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
      Internal::findPracticalInfinity(Fa.en(), v, gr.r(), Param::cALR);
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
void largeEnergyChange(double *en, TrackEnGuess *sofar_ptr, double frac_de,
                       bool toomany_nodes)
// wf did not have correct number of nodes. Make a large energy adjustment
// toomany_nodes=true means there were too many nodes
{
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
                           const TrackEnGuess &sofar)
// delta E = c*f(r)*[g_out(r)-g_in(r)] - evaluate at ctp
// nb: wf not yet normalised (anorm is input param)!
{
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
std::size_t findPracticalInfinity(const double en, const std::vector<double> &v,
                                  const std::vector<double> &r,
                                  const double alr)
// Find the practical infinity 'pinf'
// Step backwards from the last point (num_points-1) until
// (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
{
  auto pinf = r.size();
  while (pinf > 100 && (v[pinf - 1] - en) * r[pinf - 1] * r[pinf - 1] > alr) {
    assert(pinf > 1);
    --pinf;
  }
  return pinf;
}

//==============================================================================
std::size_t findClassicalTurningPoint(const double en,
                                      const std::vector<double> &v,
                                      std::size_t pinf, std::size_t d_ctp)
// Finds classical turning point 'ctp'
// Enforced to be between (0+ctp) and (pinf-ctp)
//  V(r) > E        [nb: both V and E are <0]
{
  const auto low = std::lower_bound(v.begin() + long(d_ctp + 1),
                                    v.begin() + long(pinf - d_ctp - 1), en);

  const auto distance = std::distance(v.begin(), low);

  return distance > 20 ? std::size_t(distance - 1) : ((pinf - d_ctp + 20) / 2);
}

//==============================================================================
int countNodes(const std::vector<double> &f, std::size_t pinf)
// Just counts the number of times orbital (f) changes sign
{
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
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v,
                        const std::vector<double> &H_mag, const Grid &gr,
                        std::size_t ctp, std::size_t d_ctp, std::size_t pinf,
                        const double alpha, const DiracSpinor *const VxFa,
                        const DiracSpinor *const Fa0, double zion, double mass)
// Performs inward (from pinf) and outward (from r0) integrations for given
// energy. Intergated in/out towards ctp +/- d_ctp [class. turn. point]
// Then, joins solutions, including weighted meshing around ctp +/ d_ctp
// Also: stores dg [the difference: (gout-gin)], which is used for PT
{

  DiracDerivative Hd(gr, v, ka, en, alpha, H_mag, VxFa, Fa0, zion, mass);
  solve_Dirac_outwards(f, g, Hd, ctp + d_ctp + 1);
  std::vector<double> f_in(gr.num_points()), g_in(gr.num_points());
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
  // Average of points ctp+/-d_ctp [using a weighted average]
  double rescale = f[ctp] / f_in[ctp];
  double denom = 1;
  for (std::size_t i = 1; i <= d_ctp; i++) {
    auto w = Param::weight(i);
    rescale +=
        0.5 * (f[ctp + i] / f_in[ctp + i] + f[ctp - i] / f_in[ctp - i]) * w;
    denom += w;
  }
  rescale /= denom;

  // store difference between in/out solutions (for g) - after re-scaling
  // Used later for P.T.
  for (auto i = 0ul; i < dg.size(); i++) {
    dg[i] = g[ctp - d_ctp + i] - g_in[ctp - d_ctp + i] * rescale;
  }

  // Join the in and outward solutions. "Meshed" around ctp +/- d_ctp
  for (std::size_t i = ctp - d_ctp; i <= ctp + d_ctp; i++) {
    //"Mesh" in the intermediate region, using weighted av.
    const auto b =
        (i - ctp > 0) ? 1.0 - Param::weight(i - ctp) : Param::weight(i - ctp);
    const auto a = 1.0 - b;
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

    const auto g_f_ratio = (ka > 0) ? (ka + ga0) / az0 : az0 / (ka - ga0);

    // const auto g_f_ratio =
    //     DiracHydrogen::gfratio(r[0], ka, Z_eff, alpha, Hd.en, Hd.mass);

    f0 = 2.0 * std::pow(r[0], ga0);
    g0 = f0 * g_f_ratio;
  }

  AdamsMoulton::ODESolver2D<Param::K_Adams, std::size_t, double> ode{du, &Hd};

  ode.solve_initial_K(0, f0, g0);
  for (std::size_t i = 0; i < ode.K_steps(); ++i) {
    f.at(i) = ode.f.at(i);
    g.at(i) = ode.g.at(i);
  }
  for (std::size_t i = ode.K_steps(); i < pinf; ++i) {
    ode.drive(i);
    f.at(i) = ode.last_f();
    g.at(i) = ode.last_g();
  }
}

//==================================================================
void solve_Dirac_inwards(std::vector<double> &f, std::vector<double> &g,
                         const DiracDerivative &Hd, std::size_t nf,
                         std::size_t pinf, double mass)
// Program to start the INWARD integration.
// Starts from Pinf, and uses an expansion [WKB approx] to go to (pinf-K_Adams)
// i.e., gets last K_Adams points
// Then, it then call ADAMS-MOULTON, to finish (from num_loops*K_Adams+1
//   to nf = ctp-d_ctp)
{

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
    f.at(ode.t.at(i)) = ode.f.at(i);
    g.at(ode.t.at(i)) = ode.g.at(i);
  }
  const auto i_start = ode.last_t();
  for (std::size_t i = i_start - 1; i >= nf; --i) {
    ode.drive(i);
    f.at(i) = ode.last_f();
    g.at(i) = ode.last_g();
    if (i == 0)
      break;
  }
  for (std::size_t i = i_start + ode.K_steps(); i < f.size(); ++i) {
    f.at(i) = 0.0;
    g.at(i) = 0.0;
  }
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
