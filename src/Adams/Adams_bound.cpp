#include "Adams_bound.hpp"
#include "DiracSpinor.hpp"
#include "Grid.hpp"
#include "Matrix_linalg.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

/*
Program to solve single-electron bound-state Dirac problem for a (given)
local, central potential.
Based on method presented in book by W. Johnson.
Employs the Adams-Moulton method.
solveDBS is the main routine that is called from elsewhere.
All other functions called by solveDBS.
*/

namespace Adams {

#define DO_DEBUG false
#if DO_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif // DEBUG

static const AdamsCoefs<AMO> AMcoef; // Adamns-Moulton coeficients

// weighting function for meshing in/out solutions
// nb: must be positive, but i may be negative (?) [ctp - d_ctp]
static auto weight = [](int i) { return 1. / double(i * i + 1); };
static const double cALR = 600; // 'assymptotically large r [kinda..]' (=800)

//******************************************************************************
void solveDBS(DiracSpinor &psi, const std::vector<double> &v, const Grid &rgrid,
              const double alpha, int log_dele)
// Solves local, spherical bound state dirac equation using Adams-Moulton
// method. Based on method presented in book by W. R. Johnson:
//   W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007)
// I have added a few extensions to this method. In particular, I integrate past
// the classical turning point. (See below)
//
// See also:
//  * https://en.wikipedia.org/wiki/Linear_multistep_method
//  * Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (1993),
//    Solving ordinary differential equations I: Nonstiff problems (2nd ed.),
//    Berlin: Springer Verlag, ISBN 978-3-540-56670-0.
//  * Quarteroni, Alfio; Sacco, Riccardo; Saleri, Fausto (2000),
//    Matematica Numerica, Springer Verlag, ISBN 978-88-470-0077-3.
//  * http://mathworld.wolfram.com/AdamsMethod.html
//
// Rough description of method:
// 1. Start with initial 'guess' of energy
// 2. Find "practical infinity" (psi~0), and Classical turning point (e=v)
// 3. Performs 'inward' integration (Adams Moulton).
//    Integrates inwards from the practical infinity inwards to d_ctp points
//    past the classical turning point (ctp).
// 4. Performs 'outward' integration (Adams Moulton). Integrates from 0
//    outwards to d_ctp points past the ctp.
// 5. Matches the two functions around ctp, by re-scaling the 'inward' solution
//    (for f). Uses a weighted mean, with weights given by distance from ctp.
// 6. Checks the number of nodes the wf has. If too many or too few nodes, makes
//    a large change to the energy and tries again (from step 2).
//    If the correct number of nodes, uses perturbation theory to make minor
//    corrections to the energy to 'zoom in' (matching the in/out solution for
//    g), then re-starts from step 2.
// Continues until this energy adjustment falls below a prescribed threshold.
//
// Orbitals defined:
//   psi := (1/r) {f O_k, ig O_(-k)}
//
{
  // Parameters.
  static const int max_its = 99; // Max # attempts at converging [sove bs] (30)
  static const double lfrac_de = 0.12; // 'large' energy variations (0.1 => 10%)
  static const int d_ctp = 4;          // Num points past ctp +/- d_ctp.

  // Convergance goal. Default: 1e-14
  const double eps_goal = (log_dele > 0) ? std::pow(10, -log_dele) : 1.e-14;

  DEBUG(if (!(std::abs(psi.k) <= psi.n && psi.k != psi.n)) {
    std::cerr << "\nFail96 in Adams: bad state " << psi.symbol() << "\n";
    return;
  })

  DEBUG(std::cerr << "Start: " << psi.symbol() << ", en=" << psi.en << "\n";)

  // orbital should have (n-l-1) nodes:
  const int required_nodes = psi.n - psi.l() - 1;
  bool correct_nodes = false;
  // Number of times there was too many/too few nodes:
  int count_toomany = 0, count_toofew = 0;
  // Upper and lower energy window before correct # nodes
  double high_en = 0, low_en = 0;

  // Start eigenvalue iterations:
  double t_en = psi.en;
  int t_pinf = 0;
  double t_eps = 1;
  double t_eps_prev = 1;
  double anorm = 0;
  int t_its = 1;
  for (; t_its < max_its; ++t_its) {
    t_pinf = findPracticalInfinity(t_en, v, rgrid.r, cALR);
    const int ctp = findClassicalTurningPoint(t_en, v, t_pinf, d_ctp);

    // Find solution (f,g) to DE for given energy:
    // Also stores dg (gout-gin) for PT [used for PT to find better e]
    std::vector<double> dg(2 * d_ctp + 1);
    trialDiracSolution(psi.f, psi.g, dg, t_en, psi.k, v, rgrid, ctp, d_ctp,
                       t_pinf, alpha);

    const int counted_nodes = countNodes(psi.f, t_pinf);

    // If correct number of nodes, use PT to make minor energy adjustment.
    // Otherwise, make large adjustmunt until correct # of nodes
    double en_old = t_en;
    if (counted_nodes == required_nodes) {
      correct_nodes = true;
      anorm = calcNorm(psi.f, psi.g, rgrid.drdu, rgrid.du, t_pinf);
      t_en = smallEnergyChangePT(en_old, anorm, psi.f, dg, ctp, d_ctp, alpha,
                                 count_toofew, count_toomany, low_en, high_en);
    } else {
      correct_nodes = false;
      bool toomany_nodes = (counted_nodes > required_nodes) ? true : false;
      largeEnergyChange(t_en, count_toomany, count_toofew, high_en, low_en,
                        lfrac_de, toomany_nodes);
    }
    t_eps = fabs((t_en - en_old) / en_old);

    DEBUG(std::cerr << " :: it=" << t_its << " nodes:" << counted_nodes << "/"
                    << required_nodes << " new_en = " << t_en
                    << " delta=" << t_eps * t_en << " eps=" << t_eps << "\n";
          std::cin.get();)

    auto getting_worse =
        (t_its > 10 && t_eps >= 1.2 * t_eps_prev && correct_nodes);
    auto converged = (t_eps < eps_goal && correct_nodes);
    if (converged || getting_worse)
      break;
    t_eps_prev = t_eps;
  } // END itterations

  // If never got correct nodes, never calc'd norm constant.
  // This is rare - means a failure. But hopefully, failure will go away after
  // a few count_toomany HF iterations..If we don't norm wf, HF will fail.
  if (!correct_nodes) {
    anorm = calcNorm(psi.f, psi.g, rgrid.drdu, rgrid.du, t_pinf);
    DEBUG(std::cerr << "\nFAIL-148: wrong nodes:" << countNodes(psi.f, t_pinf)
                    << "/" << required_nodes << " for " << psi.symbol()
                    << "\n";)
  }

  // store energy etc.
  psi.en = t_en;
  psi.eps = t_eps;
  psi.pinf = (std::size_t)t_pinf;
  psi.its = t_its;

  // normalises the orbital (and zero's after pinf)
  const double an = 1. / std::sqrt(anorm);
  for (auto i = 0ul; i < psi.pinf; i++) {
    psi.f[i] = an * psi.f[i];
    psi.g[i] = an * psi.g[i];
  }
  for (auto i = psi.pinf; i < rgrid.ngp; i++) {
    psi.f[i] = 0;
    psi.g[i] = 0;
  }

  return;
}

//******************************************************************************
void diracODE_regularAtOrigin(DiracSpinor &phi, const double en,
                              const std::vector<double> &v,
                              const double alpha) {
  const auto &gr = phi.p_rgrid;
  if (en != 0)
    phi.en = en;
  auto pinf = findPracticalInfinity(phi.en, v, gr->r, cALR);
  outwardAM(phi.f, phi.g, phi.en, v, phi.k, gr->r, gr->drdu, gr->du, pinf - 1,
            alpha);
  // for safety: make sure zerod! (I may re-use existing orbitals!)
  for (std::size_t i = pinf; i < gr->ngp; i++) {
    phi.f[i] = 0;
    phi.g[i] = 0;
  }
  phi.pinf = pinf;
}
//******************************************************************************
void diracODE_regularAtInfinity(DiracSpinor &phi, const double en,
                                const std::vector<double> &v,
                                const double alpha) {
  const auto &gr = phi.p_rgrid;
  if (en < 0)
    phi.en = en;
  auto pinf = findPracticalInfinity(phi.en, v, gr->r, cALR);
  inwardAM(phi.f, phi.g, phi.en, v, phi.k, gr->r, gr->drdu, gr->du, 0, pinf - 1,
           alpha);
  for (std::size_t i = pinf; i < gr->ngp; i++) {
    phi.f[i] = 0;
    phi.g[i] = 0;
  }
  phi.pinf = pinf;
}

//******************************************************************************
void largeEnergyChange(double &en, int &count_toomany, int &count_toofew,
                       double &high_en, double &low_en, double lfrac_de,
                       bool toomany_nodes)
// wf did not have correct number of nodes. Make a large energy adjustment
// toomany_nodes=true means there were too many nodes
{
  if (toomany_nodes) {
    ++count_toomany;
    if ((count_toomany == 1) || (en < high_en))
      high_en = en;
    double etemp = (1. + lfrac_de) * en;
    if ((count_toofew != 0) && (etemp < low_en))
      etemp = 0.5 * (high_en + low_en);
    en = etemp;
  } else {
    ++count_toofew;
    if ((count_toofew == 1) || (en > low_en))
      low_en = en;
    double etemp = (1. - lfrac_de) * en;
    if ((count_toomany != 0) && (etemp > high_en))
      etemp = 0.5 * (high_en + low_en);
    en = etemp;
  }
}

//******************************************************************************
double calcNorm(const std::vector<double> &f, const std::vector<double> &g,
                const std::vector<double> &drdu, const double du,
                const int pinf) {
  double anormF = NumCalc::integrate({&f, &f, &drdu}, 1.0, 0, pinf);
  double anormG = NumCalc::integrate({&g, &g, &drdu}, 1.0, 0, pinf);
  return (anormF + anormG) * du;
}

//******************************************************************************
double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, int ctp, int d_ctp,
                           double alpha, int count_toofew, int count_toomany,
                           double low_en, double high_en)
// delta E = c*f(r)*[g_out(r)-g_in(r)] - evaluate at ctp
// nb: wf not yet normalised (anorm is input param)!
{
  double p_del_q = f[ctp] * dg[d_ctp];
  double denom = 1.;

  // weighted average around ctp:
  for (int i = 1; i <= d_ctp; i++) {
    auto w = weight(i);
    p_del_q +=
        0.5 * (f[ctp + i] * dg[d_ctp + i] + f[ctp - i] * dg[d_ctp - i]) * w;
    denom += w;
  }

  double de = p_del_q / (alpha * anorm * denom);
  double new_en = en + de;

  if ((count_toofew != 0) && (new_en < low_en)) {
    new_en = 0.5 * (en + low_en);
  } else if ((count_toomany != 0) && (new_en > high_en)) {
    new_en = 0.5 * (en + high_en);
  } else if (new_en > 0) {
    // This only happens v. rarely. nodes correct, but P.T. gives silly result!
    if (de > 0)
      new_en = 0.9 * en;
    else
      new_en = 1.1 * en;
  }

  return new_en;
}

//******************************************************************************
int findPracticalInfinity(const double en, const std::vector<double> &v,
                          const std::vector<double> &r, const double alr)
// Find the practical infinity 'pinf'
// Step backwards from the last point (ngp-1) until
// (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
// XXX Note: unsafe, and a little slow. Would be better to use
// std::lower_bound.. but would need lambda or something for r^2 part?
{
  auto pinf = r.size() - 1;
  while ((en - v[pinf]) * r[pinf] * r[pinf] + alr < 0) {
    --pinf;
    DEBUG(if (pinf == 0) {
      std::cerr << "Fail290 in WF: pinf underflowed?\n";
      std::cin.get();
    })
  }

  return (int)pinf;
}

//******************************************************************************
int findClassicalTurningPoint(const double en, const std::vector<double> &v,
                              const int pinf, const int d_ctp)
// Finds classical turning point 'ctp'
// Enforced to be between (0+ctp) and (pinf-ctp)
//  V(r) > E        [nb: both V and E are <0]
{
  auto low =
      std::lower_bound(v.begin() + d_ctp + 1, v.begin() + pinf - d_ctp, en);
  return (int)(low - v.begin()) - 1;
}

//******************************************************************************
int countNodes(const std::vector<double> &f, const int pinf)
// Just counts the number of times orbital (f) changes sign
{
  double sp = f[1];
  int counted_nodes = 0;
  for (int i = 2; i < pinf; ++i) {
    if (sp * f[i] < 0) { // ok if f[i]==0 ?? (-0<0?) XXX
      ++counted_nodes;
      sp = f[i];
    }
  }
  return counted_nodes;
}

//******************************************************************************
void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const double en, const int ka,
                        const std::vector<double> &v, const Grid &gr,
                        const int ctp, const int d_ctp, const int pinf,
                        const double alpha)
// Performs inward (from pinf) and outward (from r0) integrations for given
// energy. Intergated in/out towards ctp +/- d_ctp [class. turn. point]
// Then, joins solutions, including weighted meshing around ctp +/ d_ctp
// Also: stores dg [the difference: (gout-gin)], which is used for PT
{
  outwardAM(f, g, en, v, ka, gr.r, gr.drdu, gr.du, ctp + d_ctp, alpha);
  std::vector<double> f_in(gr.ngp), g_in(gr.ngp);
  inwardAM(f_in, g_in, en, v, ka, gr.r, gr.drdu, gr.du, ctp - d_ctp, pinf - 1,
           alpha);
  joinInOutSolutions(f, g, dg, f_in, g_in, ctp, d_ctp, pinf);
}

//******************************************************************************
void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg,
                        const std::vector<double> &f_in,
                        const std::vector<double> &g_in, const int ctp,
                        const int d_ctp, const int pinf) {

  // Find the re-scaling factor (for inward soln)
  // Average of points ctp+/-d_ctp [using a weighted average]
  double rescale = f[ctp] / f_in[ctp];
  double denom = 1;
  for (int i = 1; i <= d_ctp; i++) {
    auto w = weight(i);
    rescale +=
        0.5 * (f[ctp + i] / f_in[ctp + i] + f[ctp - i] / f_in[ctp - i]) * w;
    denom += w;
  }
  rescale /= denom;

  // store difference between in/out solutions (for g) - after re-scaling
  // Used later for P.T.
  for (auto i = 0; i < (int)dg.size(); i++) {
    dg[i] = g[ctp - d_ctp + i] - g_in[ctp - d_ctp + i] * rescale;
  }

  // Join the in and outward solutions. "Meshed" around ctp +/- d_ctp
  for (int i = ctp - d_ctp; i <= ctp + d_ctp; i++) {
    //"Mesh" in the intermediate region, using weighted av.
    auto b = (i - ctp > 0) ? 1. - weight(i - ctp) : weight(i - ctp);
    auto a = 1. - b;
    f[i] = a * f[i] + b * f_in[i] * rescale;
    g[i] = a * g[i] + b * g_in[i] * rescale;
  }
  for (int i = ctp + d_ctp + 1; i < pinf; i++) {
    f[i] = f_in[i] * rescale;
    g[i] = g_in[i] * rescale;
  }
}

// Adams putward Coeficients
static const int NOL = 1; // # of outdir runs [finds first NOL*AMO+1 points (3)]
//******************************************************************************
void outwardAM(std::vector<double> &f, std::vector<double> &g, const double en,
               const std::vector<double> &v, const int ka,
               const std::vector<double> &r, const std::vector<double> &drdu,
               const double du, const int nf, const double alpha)
// Program to start the OUTWARD integration.
// Starts from 0, and uses an expansion(?) to go to (NOL*AMO).
// Then, it then call ADAMS-MOULTON, to finish
// (from NOL*AMO+1 to nf = ctp+d_ctp)
{
  const double c2 = 1. / (alpha * alpha);

  double az = -1 * v[1] * r[1] * alpha; //  Z = -1 * v[AMO] * r[AMO]
  // take average? Or maybe update each point?? Or lowest?
  double ga = std::sqrt(ka * ka - az * az);

  // initial wf values
  // P(r) = r^gamma u(r) // f(r) = P(r)
  // Q(r) = r^gamma v(r) // g(r) = -Q(r)
  double u0 = 1;
  double v0 = (ka > 0) ? -(ga + ka) / az : az / (ga - ka);
  f[0] = std::pow(r[0], ga) * u0;
  g[0] = -std::pow(r[0], ga) * v0;

  // loop through and find first NOL*AMO points of wf
  for (int ln = 0; ln < NOL; ln++) {
    int i0 = ln * AMO + 1;

    // defines/populates em expansion coeficients (then inverts)
    std::array<double, AMO> coefa, coefb, coefc, coefd;
    Matrix::SqMatrix em(AMO);
    const auto oid_du = AMcoef.OId * du;
    for (int i = 0; i < AMO; i++) {
      double dror = drdu[i + i0] / r[i + i0];
      coefa[i] = (-oid_du * (ga + ka) * dror);
      coefb[i] = (-oid_du * (en + 2 * c2 - v[i + i0]) * drdu[i + i0] * alpha);
      coefc[i] = (oid_du * (en - v[i + i0]) * drdu[i + i0] * alpha);
      coefd[i] = (-oid_du * (ga - ka) * dror);
      for (int j = 0; j < AMO; j++) {
        em[i][j] = AMcoef.OIe[i][j];
      }
      em[i][i] -= coefd[i];
    }
    // from here on, em is the inverted matrix
    em.invert();

    // defines/populates fm, s coefs
    std::array<double, AMO> s;
    Matrix::SqMatrix fm(AMO);
    for (int i = 0; i < AMO; i++) {
      s[i] = -AMcoef.OIa[i] * u0;
      for (int j = 0; j < AMO; j++) {
        fm[i][j] = AMcoef.OIe[i][j] - coefb[i] * em[i][j] * coefc[j];
        s[i] -= coefb[i] * em[i][j] * AMcoef.OIa[j] * v0;
      }
      fm[i][i] -= coefa[i];
    }
    // from here on, fm is the inverted matrix
    fm.invert();

    // writes u(r) in terms of coefs and the inverse of fm
    // P(r) = r^gamma u(r)
    std::array<double, AMO> us;
    for (int i = 0; i < AMO; i++) {
      us[i] = 0;
      for (int j = 0; j < AMO; j++) {
        us[i] += fm[i][j] * s[j];
      }
    }

    // writes v(r) in terms of coefs + u(r)
    // Q(r) = r^gamma v(r)
    std::array<double, AMO> vs;
    for (int i = 0; i < AMO; i++) {
      vs[i] = 0;
      for (int j = 0; j < AMO; j++) {
        vs[i] += em[i][j] * (coefc[j] * us[j] - AMcoef.OIa[j] * v0);
      }
    }

    // writes wavefunction: P= r^gamma u(r) etc..
    for (int i = 0; i < AMO; i++) {
      double r_ga = std::pow(r[i + i0], ga);
      f[i + i0] = r_ga * us[i];
      g[i + i0] = -r_ga * vs[i];
    }

    // re-sets 'starting point' for next ln
    u0 = us.back();
    v0 = vs.back();

  } // END for (int ln=0; ln<NOL; ln++)  [loop through outint `NOL' times]

  // Call adamsmoulton to finish integration from (NOL*AMO) to nf = ctp+d_ctp
  int na = NOL * AMO + 1;
  if (nf > na)
    adamsMoulton(f, g, en, v, ka, r, drdu, du, na, nf, alpha);

  return;
}

// order of the expansion coeficients in 'inwardAM'  (15 orig.)
static const int NX = 15;
// PRIMARY convergance for expansion in `inwardAM'' (10^-8)
static const double NXEPSP = 1.e-10;
//******************************************************************
void inwardAM(std::vector<double> &f, std::vector<double> &g, const double en,
              const std::vector<double> &v, const int ka,
              const std::vector<double> &r, const std::vector<double> &drdu,
              const double du, const int nf, const int pinf, const double alpha)
// Program to start the INWARD integration.
// Starts from Pinf, and uses an expansion(?) to go to (pinf-AMO)
// Then, it then call ADAMS-MOULTON, to finish
// (from NOL*AMO+1 to nf = ctp-d_ctp)
{

  const double alpha2 = alpha * alpha;
  const double cc = 1. / alpha;
  const double c2 = 1. / alpha2;

  const double lambda = std::sqrt(-en * (2. + en * alpha2));
  const double zeta = -v[pinf] * r[pinf];
  const double sigma = (1. + en * alpha2) * (zeta / lambda);
  const double Ren = en + c2; // total relativistic energy

  // Generates the expansion coeficients for asymptotic wf
  // up to order NX (NX is 'param')
  std::array<double, NX> bx;
  std::array<double, NX> ax;
  bx[0] = (ka + (zeta / lambda)) * (0.5 * alpha);
  double ka2 = (double)(ka * ka);
  double zeta2 = zeta * zeta;
  for (int i = 0; i < NX; i++) {
    ax[i] = (ka + (i + 1 - sigma) * Ren * alpha2 - zeta * lambda * alpha2) *
            bx[i] * cc / ((i + 1) * lambda);
    if (i < (NX - 1))
      bx[i + 1] =
          (ka2 - std::pow((double(i + 1) - sigma), 2) - zeta2 * alpha2) *
          bx[i] / (2 * (i + 1) * lambda);
  }

  // Generates last `AMO' points for P and Q [actually AMO+1?]
  double f1 = std::sqrt(1. + en * alpha2 * 0.5);
  double f2 = std::sqrt(-en * 0.5) * alpha;
  for (int i = pinf; i >= (pinf - AMO); i--) {
    double rfac = std::pow(r[i], sigma) * exp(-lambda * r[i]);
    double ps = 1.;
    double qs = 0.;
    double rk = 1.;
    double xe = 1.;
    for (int k = 0; k < NX; k++) {
      // this will loop until a) converge, or b) k=NX
      rk *= r[i];
      ps += (ax[k] / rk);
      qs += (bx[k] / rk);
      xe = fmax(std::fabs(ax[k] / ps), std::fabs(bx[k] / qs)) / rk;
      if (xe < NXEPSP)
        break;
    }
    DEBUG({
      const double nxepss = 1.e-3;
      if (xe > nxepss)
        std::cerr << "WARNING: Asymp. expansion in inwardAM didn't converge: "
                  << i << " " << xe << "\n";
    })
    f[i] = rfac * (f1 * ps + f2 * qs);
    g[i] = rfac * (f1 * qs - f2 * ps);
  }

  if ((pinf - AMO - 1) >= nf)
    adamsMoulton(f, g, en, v, ka, r, drdu, du, pinf - AMO - 1, nf, alpha);

  return;
}

//******************************************************************************
void adamsMoulton(std::vector<double> &f, std::vector<double> &g,
                  const double en, const std::vector<double> &v, const int ka,
                  const std::vector<double> &r, const std::vector<double> &drdu,
                  const double du, const int ni, const int nf,
                  const double alpha)
// program finishes the INWARD/OUTWARD integrations (ADAMS-MOULTON)
//   * ni is starting (initial) point for integration
//   * nf is end (final) point for integration (nf=ctp+/-d_ctp)
{
  double c2 = 1. / (alpha * alpha); // c^2 - just to shorten code

  int inc;     //'increment' for integration (+1 for forward, -1 for backward)
  int nosteps; // number of steps integration should make
  if (nf > ni) {
    inc = 1;
    nosteps = nf - ni + 1;
  } else if (nf < ni) {
    inc = -1;
    nosteps = ni - nf + 1;
  } else {
    DEBUG(std::cerr
              << "FAIL 611: ni=nf in adamsmoulton.. no further integration\n";
          std::cin.get();)
    return;
  }

  // create arrays for wf derivatives
  auto ngp = r.size();
  std::vector<double> df(ngp), dg(ngp);
  std::array<double, AMO> amcoef;
  int k1 = ni - inc * AMO; // nb: k1 is iterated
  for (auto i = 0; i < AMO; i++) {
    double dror = drdu[k1] / r[k1];
    df[i] = inc * ((-ka * dror * f[k1]) +
                   (alpha * ((en + 2 * c2) - v[k1]) * drdu[k1] * g[k1]));
    dg[i] = inc * (ka * dror * g[k1] - alpha * (en - v[k1]) * drdu[k1] * f[k1]);
    amcoef[i] = du * AMcoef.AMd * AMcoef.AMa[i];
    k1 += inc;
  }

  // integrates the function from ni to the c.t.p
  double a0 = du * AMcoef.AMd * AMcoef.AMaa;
  int k2 = ni;
  for (int i = 0; i < nosteps; i++) {
    double dror = drdu[k2] / r[k2];
    double dai = -inc * (ka * dror);
    double dbi = -inc * alpha * (en + 2 * c2 - v[k2]) * drdu[k2];
    double dci = inc * alpha * (en - v[k2]) * drdu[k2];
    double ddi = -dai;
    double det_inv = 1. / (1. - a0 * a0 * (dbi * dci - dai * ddi));
    double sp = f[k2 - inc];
    double sq = -g[k2 - inc];
    for (int l = 0; l < AMO; l++) {
      sp += amcoef[l] * df[l];
      sq -= amcoef[l] * dg[l];
    }
    f[k2] = (sp + a0 * (dbi * sq - ddi * sp)) * det_inv;
    g[k2] = -(sq + a0 * (dci * sp - dai * sq)) * det_inv;
    // loads next 'first' k values:
    for (int l = 0; l < (AMO - 1); l++) {
      df[l] = df[l + 1];
      dg[l] = dg[l + 1];
    }
    // loads next 'first' deriv's:
    df[AMO - 1] = dai * f[k2] - dbi * g[k2];
    dg[AMO - 1] = ddi * g[k2] - dci * f[k2];
    k2 += inc;
  }

  return;
} // END adamsmoulton

} // namespace Adams
