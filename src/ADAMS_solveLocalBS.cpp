#include "ADAMS_solveLocalBS.h"
#include "DiracSpinor.h"
#include "Grid.h"
#include "Matrix_linalg.h"
#include "NumCalc_quadIntegrate.h"
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

namespace ADAMS {

#define DO_DEBUG false
#if DO_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif // DEBUG

static const AdamsCoefs<AMO> AMcoef; // coeficients

//******************************************************************************
void solveDBS(DiracSpinor &psi, const std::vector<double> &v, const Grid &rgrid,
              double alpha, int log_dele)
/*
Solves local, spherical bound state dirac equation using Adams-Moulton method.
Based on method presented in book by W. R. Johnson:
  W. R. Johnson, Atomic Structure Theory (Springer, New York, 2007)
I have added a few extensions to this method. In particular, I integrate past
the classical turning point. (See below)

See also:
 * https://en.wikipedia.org/wiki/Linear_multistep_method
 * Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (1993),
   Solving ordinary differential equations I: Nonstiff problems (2nd ed.),
   Berlin: Springer Verlag, ISBN 978-3-540-56670-0.
 * Quarteroni, Alfio; Sacco, Riccardo; Saleri, Fausto (2000),
   Matematica Numerica, Springer Verlag, ISBN 978-88-470-0077-3.
 * http://mathworld.wolfram.com/AdamsMethod.html

Rough description of method:
1. Start with initial 'guess' of energy
2. Find the "practical infinity" (psi~0), and the Classical turning point [e=v]
3. Performs 'inward' integration (Adams Moulton). Integrates from the practical
   infinity inwards to d_ctp points past the classical turning point (ctp).
4. Performs 'outward' integration (Adams Moulton). Integrates from 0
   outwards to d_ctp points past the ctp.
5. Matches the two functions around ctp, by re-scaling the 'inward' solution
   (for f). Uses a weighted mean, with weights given by distance from ctp.
6. Checks the number of nodes the wf has. If too many or too few nodes, makes
   a large change to the energy and tries again (from step 2).
   If the correct number of nodes, uses perturbation theory to make minor
   corrections to the energy to 'zoom in' (matching the in/out solution for g),
   then re-starts from step 2.
Continues until this energy adjustment falls below a prescribed threshold.

Orbitals defined:
  psi := (1/r) {f O_k, ig O_(-k)}

Was originaly in terms of (p,q). Now in terms of (f,g).
Other "sub"-programs still use (p,q).
Defn: f = p, g = -q. (My g includes alpha)

*/
{
  // Parameters.
  const int max_its = 16;       // Max # attempts at converging [sove bs] (30)
  const double alr = 800;       // ''assymptotically large r [kinda..]''  (=800)
  const double lfrac_de = 0.15; // 'large' energy variations (0.1 => 10%)
  const int d_ctp = 4;          // Num points past ctp +/- d_ctp.

  // int d_ctp = d_ctp_in; // from tests..

  // Temporary refs to make transition easier. Should remove these
  // + propogate changes through proplerly XXX XXX XXX
  auto ngp = rgrid.ngp;
  auto &r = rgrid.r;
  auto &drdu = rgrid.drdu;
  auto du = rgrid.du;

  // Temporary refs to make transition easier. Should remove these
  // + propogate changes through proplerly XXX XXX XXX
  auto &f = psi.f;
  auto &g = psi.g;
  auto &en_inout = psi.en;
  auto n = psi.n;
  auto ka = psi.k;
  auto &pinf_out = psi.pinf;
  auto &its_out = psi.its;
  auto &eps_out = psi.eps;

  // Convergance goal. Default: 1e-15
  double eps_goal = (log_dele > 0) ? 1. / pow(10, log_dele) : 1e-14;

  DEBUG( // Checks to see if legal n is requested.
      if (!((abs(ka) <= n) && (ka != n))) {
        std::cerr << "\nFail96 in Adams: bad state n,k=" << n << "," << ka
                  << "\n";
        return;
      })

  // Find 'l' from 'kappa' (ang. momentum Q number) for # of nodes
  int l = (ka > 0) ? ka : -ka - 1;
  int required_nodes = n - l - 1;
  bool correct_nodes = false;

  double en = en_inout;
  int pinf = 0;
  double anorm = 0; // normalisation constant
  double tmp_eps_en = 1;

  DEBUG(std::cerr << "Start: n=" << n << ", k=" << ka << ", en=" << en << "\n";)

  // Some parameters used by the Adams Moulton method:
  // Number of times there was too many/too few nodes:
  int more = 0, less = 0;
  // Highest/lowest energies tried before correct # of nodes:
  double eupper = 0, elower = 0;
  int its = 0;
  for (its = 1; its < max_its; its++) {
    // Find the practical infinity 'pinf' [(V(r) - E)*r^2 >  alr]
    pinf = findPracticalInfinity(en, v, r, alr);
    // Find classical turning point 'ctp' [V(r) > E ]
    int ctp = findClassicalTurningPoint(en, v, pinf, d_ctp);

    // Find solution (f,g) to DE for given energy:
    // (Inward + outward solutions joined at ctp, merged over ctp+/-d_ctp)
    // Also stores dg (gout-gin) for PT [used for PT to find better e]
    std::vector<double> dg(size_t(2 * d_ctp + 1));
    trialDiracSolution(f, g, dg, en, ka, v, r, drdu, du, ctp, d_ctp, pinf,
                       alpha);

    // Count the number of nodes (zeros) the wf has.
    int counted_nodes = countNodes(f, pinf);

    // If correct number of nodes, use PT to make minor energy adjustment.
    // Otherwise, make large adjustmunt until correct # of nodes
    double en_old = en;
    if (counted_nodes == required_nodes) {
      correct_nodes = true;
      anorm = calcNorm(f, g, drdu, du, pinf);
      en = smallEnergyChangePT(en_old, anorm, f, dg, ctp, d_ctp, alpha, less,
                               more, elower, eupper);
    } else {
      correct_nodes = false; // can happen?
      bool more_nodes = (counted_nodes > required_nodes) ? true : false;
      largeEnergyChange(en, more, less, eupper, elower, lfrac_de, more_nodes);
    }
    tmp_eps_en = fabs((en - en_old) / en_old);

    DEBUG(std::cerr << " :: it=" << its << " nodes:" << counted_nodes << "/"
                    << required_nodes << " new_en = " << en << " delta="
                    << tmp_eps_en * en << " eps=" << tmp_eps_en << "\n";
          std::cin.get();)

    if (tmp_eps_en < eps_goal && correct_nodes)
      break;
  } // end itterations

  // If never got correct nodes, never calc'd constant.
  // This is rare - means a failure. But hopefully, failure will go away after
  // a few more HF iterations..If we don't norm wf, HF will fail.
  if (!correct_nodes) {
    anorm = calcNorm(f, g, drdu, du, pinf);
    DEBUG(std::cerr << "\nFAIL-148: wrong nodes:" << countNodes(f, pinf) << "/"
                    << required_nodes << " for n,k=" << n << "," << ka << "\n";)
  }

  eps_out = tmp_eps_en;
  en_inout = en;
  pinf_out = (size_t)pinf;
  its_out = its;

  // normalises the wavefunction
  double an = 1. / sqrt(anorm);
  for (size_t i = 0; i < pinf_out; i++) {
    f[i] = an * f[i];
    g[i] = an * g[i];
  }
  for (auto i = pinf_out; i < ngp; i++) {
    f[i] = 0;
    g[i] = 0;
  }

  return;
}

//******************************************************************************
void largeEnergyChange(double &en, int &more, int &less, double &eupper,
                       double &elower, double lfrac_de, bool more_nodes)
/*
wf did not have correct number of nodes. Make a large energy adjustment
more_nodes=true means there were too many nodes
*/
{
  if (more_nodes) {
    // Too many nodes:
    more++;
    if ((more == 1) || (en < eupper))
      eupper = en;
    double etemp = (1. + lfrac_de) * en;
    if ((less != 0) && (etemp < elower))
      etemp = 0.5 * (eupper + elower);
    en = etemp;
  } else {
    // too few nodes:
    less++;
    if ((less == 1) || (en > elower))
      elower = en;
    double etemp = (1. - lfrac_de) * en;
    if ((more != 0) && (etemp > eupper))
      etemp = 0.5 * (eupper + elower);
    en = etemp;
  }
}

//******************************************************************************
double calcNorm(const std::vector<double> &f, const std::vector<double> &g,
                const std::vector<double> &drdu, const double du,
                const int pinf) {
  double anormF = NumCalc::integrate(f, f, drdu, 1., 0, pinf);
  double anormG = NumCalc::integrate(g, g, drdu, 1., 0, pinf);
  return (anormF + anormG) * du;
}

//******************************************************************************
double smallEnergyChangePT(const double en, const double anorm,
                           const std::vector<double> &f,
                           const std::vector<double> &dg, int ctp, int d_ctp,
                           double alpha, int less, int more, double elower,
                           double eupper)
/*
Uses PT to calculate small change in energy.
Also calculates (+outputs) norm constant (but doesn't normalise orbital!)
*/
{

  // Use perturbation theory to work out delta En
  // delta E = c*P(r)*[Qin(r)-Qout(r)] - evaluate at ctp
  // nb: wf not yet normalised!
  double p_del_q = f[ctp] * dg[d_ctp + 0];
  double denom = 1;

  // weighted average around ctp:
  for (int i = 1; i <= d_ctp; i++) {
    double w = 1. / (i * i + 1.0);
    p_del_q +=
        0.5 * (f[ctp + i] * dg[d_ctp + i] + f[ctp - i] * dg[d_ctp - i]) * w;
    denom += w;
  }

  p_del_q /= denom;
  double de = (1. / alpha) * p_del_q / anorm;
  double new_en = en + de;

  if ((less != 0) && (new_en < elower)) {
    new_en = 0.5 * (en + elower);
  } else if ((more != 0) && (new_en > eupper)) {
    new_en = 0.5 * (en + eupper);
  } else if (new_en > 0) {
    // This only happens v. rarely. nodes correct, but P.T. gives silly result!
    // Is this OK? Seems to work
    if (de > 0)
      new_en = 0.9 * en;
    else
      new_en = 1.1 * en;
  }

  return new_en;
}

//******************************************************************************
int findPracticalInfinity(double en, const std::vector<double> &v,
                          const std::vector<double> &r, double alr)
// Find the practical infinity 'pinf'
// Step backwards from the last point (ngp-1) until
// (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
// XXX Note: unsafe, and a little slow. Would be better to use
// std::lower_bound.. but would need lambda or something for r^2 part?
{
  auto pinf = r.size() - 1;
  while ((en - v[pinf]) * r[pinf] * r[pinf] + alr < 0)
    --pinf;

  return (int)pinf;
}

//******************************************************************************
int findClassicalTurningPoint(double en, const std::vector<double> &v, int pinf,
                              int d_ctp)
// Finds classical turning point 'ctp'
// Enforced to be between (0+ctp) and (pinf-ctp)
//  V(r) > E        [nb: both V and E are <0]
{
  auto low =
      std::lower_bound(v.begin() + d_ctp + 1, v.begin() + pinf - d_ctp, en);
  return (int)(low - v.begin()) - 1;
}

//******************************************************************************
int countNodes(const std::vector<double> &f, int maxi)
// Just counts the number of times wf changes sign
{
  int sizeof_f = (int)f.size();
  if (maxi == 0 || maxi > sizeof_f)
    maxi = sizeof_f;

  double sp = f[1];
  double spn;
  int counted_nodes = 0;
  for (int i = 2; i < maxi; i++) {
    spn = f[i];
    if (sp * spn < 0)
      ++counted_nodes;
    if (spn != 0)
      sp = spn;
  }
  return counted_nodes;
}

//******************************************************************************
void trialDiracSolution(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, double en, int ka,
                        const std::vector<double> &v,
                        const std::vector<double> &r,
                        const std::vector<double> &drdu, double du, int ctp,
                        int d_ctp, int pinf, double alpha) {
  int ngp = (int)f.size();
  // Temporary vectors for in/out integrations:
  std::vector<double> pin(ngp), qin(ngp), pout(ngp), qout(ngp);
  // Perform the "inwards integration":
  inwardAM(pin, qin, en, v, ka, r, drdu, du, ctp - d_ctp, pinf, alpha);
  // Perform the "outwards integration"
  outwardAM(pout, qout, en, v, ka, r, drdu, du, ctp + d_ctp, alpha);
  // Join in/out solutions into f,g + store dg (gout-gin) for PT
  joinInOutSolutions(f, g, dg, pin, qin, pout, qout, ctp, d_ctp, pinf);
}

//******************************************************************************
void joinInOutSolutions(std::vector<double> &f, std::vector<double> &g,
                        std::vector<double> &dg, const std::vector<double> &pin,
                        const std::vector<double> &qin,
                        const std::vector<double> &pout,
                        const std::vector<double> &qout, int ctp, int d_ctp,
                        int pinf) {
  // Find the re-scaling factor (for inward soln) using a weighted average:
  double rescale = pout[ctp] / pin[ctp];
  double denom = 1;
  for (int i = 1; i <= d_ctp; i++) {
    // re-scale from weigted average.
    double w = 1. / (i * i + 1.);
    rescale +=
        0.5 * (pout[ctp + i] / pin[ctp + i] + pout[ctp - i] / pin[ctp - i]) * w;
    denom += w;
  }
  rescale /= denom;

  // Join the in and outward solutions. "Meshed" around ctp +/- d_ctp
  for (int i = 0; i < ctp - d_ctp; i++) {
    f[i] = pout[i];
    g[i] = -1 * qout[i];
  }
  for (int i = ctp + d_ctp + 1; i <= pinf; i++) {
    f[i] = pin[i] * rescale;
    g[i] = -1 * qin[i] * rescale;
  }
  for (int i = ctp - d_ctp; i <= ctp + d_ctp; i++) {
    //"Mesh" in the intermediate region, using weighted av.
    int n = i - ctp;
    double B = (1. / (n * n + 2));
    if (n > 0)
      B = 1 - B;
    double A = 1 - B;
    f[i] = A * pout[i] + B * pin[i] * rescale;
    g[i] = -1 * (A * qout[i] + B * qin[i] * rescale);
  }

  // store difference between in/out solutions (for q) - after re-scaling
  // Used later for P.T.
  for (size_t i = 0; i < dg.size(); i++)
    dg[i] = qin[ctp - d_ctp + i] * rescale - qout[ctp - d_ctp + i];
}

// Adams putward Coeficients
const int NOL = 1; // # of outdir runs [finds first NOL*AMO+1 points (3)]
const static auto &OIE = AMcoef.OIe;
const static auto &OIA = AMcoef.OIa;
const static auto OID = AMcoef.OId;
//******************************************************************************
void outwardAM(std::vector<double> &p, std::vector<double> &q, double en,
               const std::vector<double> &v, int ka,
               const std::vector<double> &r, const std::vector<double> &drdu,
               double du, int nf, double alpha)
/*
Program to start the OUTWARD integration.
Starts from 0, and uses an expansion(?) to go to (NOL*AMO).
Then, it then call ADAMS-MOULTON, to finish (from NOL*AMO+1 to nf = ctp+d_ctp)
*/
{
  double az = -1 * v[AMO] * r[AMO] * alpha; //  Z = -1 * v[AMO] * r[AMO]
  double c2 = 1. / (alpha * alpha);
  double ga = sqrt(ka * ka - az * az);

  // initial wf values
  // P(r) = r^gamma u(r)
  // Q(r) = r^gamma v(r)
  double u0 = 1;
  double v0 = (ka > 0) ? -(ga + ka) / az : az / (ga - ka);
  p[0] = 0;
  q[0] = 0;

  // loop through and find first NOL*AMO points of wf
  for (int ln = 0; ln < NOL; ln++) {
    // re-work out ga (from az) in here? Poss. slightly diff. z_eff (?XX)
    int i0 = ln * AMO + 1;

    // defines/populates em coefs
    std::array<double, AMO> coefa, coefb, coefc, coefd;
    // std::array<std::array<double, AMO>, AMO> em;
    Matrix::SqMatrix em(AMO);
    for (int i = 0; i < AMO; i++) {
      double dror = drdu[i + i0] / r[i + i0];
      coefa[i] = (-OID * du * (ga + ka) * dror);
      coefb[i] = (-OID * du * (en + 2 * c2 - v[i + i0]) * drdu[i + i0] * alpha);
      coefc[i] = (OID * du * (en - v[i + i0]) * drdu[i + i0] * alpha);
      coefd[i] = (-OID * du * (ga - ka) * dror);
      for (int j = 0; j < AMO; j++)
        em[i][j] = OIE[i][j];
      em[i][i] = em[i][i] - coefd[i];
    }
    // //inverts the em matrix
    // em = Matrix::invert(em); // from here on, em is the inverted matrix
    em.invert();

    // defines/populates fm, s coefs
    std::array<double, AMO> s;
    // std::array<std::array<double, AMO>, AMO> fm;
    Matrix::SqMatrix fm(AMO);
    for (int i = 0; i < AMO; i++) {
      s[i] = -OIA[i] * u0;
      for (int j = 0; j < AMO; j++) {
        fm[i][j] = OIE[i][j] - coefb[i] * em[i][j] * coefc[j];
        s[i] = s[i] - coefb[i] * em[i][j] * OIA[j] * v0;
      }
      fm[i][i] = fm[i][i] - coefa[i];
    }
    // inverts the matrix!  fm =-> Inv(fm)
    // fm = Matrix::invert(fm); // from here on, fm is the inverted matrix
    fm.invert();

    // writes u(r) in terms of coefs and the inverse of fm
    // P(r) = r^gamma u(r)
    std::array<double, AMO> us;
    for (int i = 0; i < AMO; i++) {
      us[i] = 0;
      for (int j = 0; j < AMO; j++)
        us[i] = us[i] + fm[i][j] * s[j];
    }

    // writes v(r) in terms of coefs + u(r)
    // Q(r) = r^gamma v(r)
    std::array<double, AMO> vs;
    for (int i = 0; i < AMO; i++) {
      vs[i] = 0;
      for (int j = 0; j < AMO; j++)
        vs[i] = vs[i] + em[i][j] * (coefc[j] * us[j] - OIA[j] * v0);
      //??? here: is this the large cancellation?
    }

    // writes wavefunction: P= r^gamma u(r) etc..
    for (int i = 0; i < AMO; i++) {
      p[i + i0] = pow(r[i + i0], ga) * us[i];
      q[i + i0] = pow(r[i + i0], ga) * vs[i];
    }

    // re-sets 'starting point' for next ln
    u0 = us[AMO - 1];
    v0 = vs[AMO - 1];

  } // END for (int ln=0; ln<NOL; ln++)  [loop through outint `NOL' times]

  // Call adamsmoulton to finish integration from (NOL*AMO+1) to nf = ctp+d_ctp
  int na = NOL * AMO + 1;
  if (nf > na)
    adamsMoulton(p, q, en, v, ka, r, drdu, du, na, nf, alpha);

  return;
}

// order of the expansion coeficients in 'inwardAM'  (15 orig.)
const static int NX = 15;
// PRIMARY convergance for expansion in `inwardAM'' (10^-8)
const static double NXEPSP = 1.e-10;
//******************************************************************
void inwardAM(std::vector<double> &p, std::vector<double> &q, double en,
              const std::vector<double> &v, int ka,
              const std::vector<double> &r, const std::vector<double> &drdu,
              double du, int nf, int pinf, double alpha)
/*
Program to start the INWARD integration.
Starts from Pinf, and uses an expansion(?) to go to (pinf-AMO)
Then, it then call ADAMS-MOULTON, to finish (from NOL*AMO+1 to nf = ctp-d_ctp)
*/
{

  double alpha2 = alpha * alpha; //(alpha, 2);
  double cc = 1. / alpha;
  double c2 = 1. / alpha2;

  double lambda = sqrt(-en * (2. + en * alpha2));
  double zeta = -v[pinf] * r[pinf];
  double sigma = (1. + en * alpha2) * (zeta / lambda);
  double Ren = en + c2; // total relativistic energy

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
      bx[i + 1] = (ka2 - pow((double(i + 1) - sigma), 2) - zeta2 * alpha2) *
                  bx[i] / (2 * (i + 1) * lambda);
  }

  // Generates last `AMO' points for P and Q [actually AMO+1?]
  double f1 = sqrt(1. + en * alpha2 / 2.);
  double f2 = sqrt(-en / 2.) * alpha;
  for (int i = pinf; i >= (pinf - AMO); i--) {
    double rfac = pow(r[i], sigma) * exp(-lambda * r[i]);
    double ps = 1.;
    double qs = 0.;
    double rk = 1.;
    double xe = 1.;
    for (int k = 0; k < NX; k++) { // this will loop until a) converge, b) k=NX
      rk = rk * r[i];
      ps = ps + (ax[k] / rk);
      qs = qs + (bx[k] / rk);
      xe = fmax(fabs((ax[k] / rk) / ps), fabs((bx[k] / rk) / qs));
      if (xe < NXEPSP)
        break; // reached convergance
    }
    DEBUG(
        // SECONDARY convergance for expansion in `inint'' (10e-3):
        const double nxepss = 1.e-3;
        if (xe > nxepss) std::cerr
        << "WARNING: Asymp. expansion in ININT didn't converge: " << i << " "
        << xe << "\n";)
    p[i] = rfac * (f1 * ps + f2 * qs);
    q[i] = rfac * (f2 * ps - f1 * qs); //??? here? the 'small cancellation'(?)
  }

  // calls adams-moulton
  if ((pinf - AMO - 1) >= nf)
    adamsMoulton(p, q, en, v, ka, r, drdu, du, pinf - AMO - 1, nf, alpha);

  return;
}

// AM coefs.
static const auto &AMA = AMcoef.AMa;
static const auto AMD = AMcoef.AMd;
static const auto AMAA = AMcoef.AMaa;
//******************************************************************************
void adamsMoulton(std::vector<double> &p, std::vector<double> &q, double en,
                  const std::vector<double> &v, int ka,
                  const std::vector<double> &r, const std::vector<double> &drdu,
                  double du, int ni, int nf, double alpha)
/*
program finishes the INWARD/OUTWARD integrations (ADAMS-MOULTON)
  //- ni is starting (initial) point for integration
  //- nf is end (final) point for integration (nf=ctp+/-d_ctp)
*/
{
  double c2 = 1. / (alpha * alpha); // c^2 - just to shorten code

  int ngp = (int)r.size();

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
  std::vector<double> dp(ngp), dq(ngp);
  std::array<double, AMO> amcoef;
  int k1 = ni - inc * AMO;
  for (int i = 0; i < AMO; i++) { // nb: k1 is iterated
    double dror = drdu[k1] / r[k1];
    dp[i] = inc * (-ka * dror * p[k1] -
                   alpha * ((en + 2 * c2) - v[k1]) * drdu[k1] * q[k1]);
    dq[i] = inc * (ka * dror * q[k1] + alpha * (en - v[k1]) * drdu[k1] * p[k1]);
    amcoef[i] = du * AMD * AMA[i];
    k1 += inc;
  }

  // integrates the function from ni to the c.t.p
  double a0 = du * AMD * AMAA;
  int k2 = ni;
  for (int i = 0; i < nosteps; i++) {
    double dror = drdu[k2] / r[k2];
    double dai = -inc * (ka * dror);
    double dbi = -inc * alpha * (en + 2 * c2 - v[k2]) * drdu[k2];
    double dci = inc * alpha * (en - v[k2]) * drdu[k2];
    double ddi = -dai;
    double det_inv = 1. / (1. - a0 * a0 * (dbi * dci - dai * ddi));
    double sp = p[k2 - inc];
    double sq = q[k2 - inc];
    for (int l = 0; l < AMO; l++) {
      sp = sp + amcoef[l] * dp[l];
      sq = sq + amcoef[l] * dq[l];
    }
    p[k2] = (sp + a0 * (dbi * sq - ddi * sp)) * det_inv;
    q[k2] = (sq + a0 * (dci * sp - dai * sq)) * det_inv;
    for (int l = 0; l < (AMO - 1); l++) { // loads next 'first' k values (?)
      dp[l] = dp[l + 1];
      dq[l] = dq[l + 1];
    }
    dp[AMO - 1] = dai * p[k2] + dbi * q[k2]; // loads next 'first' deriv's (?)
    dq[AMO - 1] = dci * p[k2] + ddi * q[k2];
    k2 += inc;
  }

  return;
} // END adamsmoulton

} // namespace ADAMS
