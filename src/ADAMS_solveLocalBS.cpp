#include "ADAMS_solveLocalBS.h"
#include "Matrix_linalg.h"
#include "NumCalc_quadIntegrate.h"
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

/*
Thu 06 Nov 2014 23:29:57 AEDT. Updated 2018 to go past ctp.
Program to solve single-electron bound-state Dirac problem for a (given)
local, central potential.
Based on method presented in book by W. Johnson.
Employs the Adams-Moulton method.
solveDBS is the main routine that is called from elsewhere.
All other functions called by solveDBS.

###==  To Do / ISSUES  ###==
  * If not correct number of nodes, will finish without saying anything
   ... maybe, if not correct nodes, allow program to keep trying!
  * Gets the adams and out coeficients many times!! Not efficient!!
    (may actually be significant! Allocating memory!)

*/

namespace ADAMS {

#define DO_DEBUG false
#if DO_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif // DEBUG

//******************************************************************************
int solveDBS(std::vector<double> &f, std::vector<double> &g, double &en_inout,
             const std::vector<double> &v, int n, int ka,
             const std::vector<double> &r, const std::vector<double> &drdt,
             double h, int &pinf_out, int &its_out, double &eps_out,
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
  const int max_its = 32;      // Max # attempts at converging [sove bs] (30)
  const double alr = 800;      // ''assymptotically large r [kinda..]''  (=800)
  const double lfrac_de = 0.2; // 'large' energy variations (0.1 => 10%)
  const int d_ctp_in = 4;      // Num points past ctp +/- d_ctp.

  int d_ctp = d_ctp_in; // from tests..

  int ngp = (int)r.size();

  // Convergance goal. Default: 1e-15
  double dele_goal = (log_dele > 0) ? 1. / pow(10, log_dele) : 1e-15;

  DEBUG( // Checks to see if legal n is requested.
      if (!((abs(ka) <= n) && (ka != n))) {
        std::cerr << "\nFail96 in Adams: bad state n,k=" << n << "," << ka
                  << "\n";
        return 1;
      })

  // Find 'l' from 'kappa' (ang. momentum Q number) for # of nodes
  int l = (ka > 0) ? ka : -ka - 1;
  int required_nodes = n - l - 1;
  bool correct_nodes = false;

  double en = en_inout;
  int pinf = -1;
  double anorm = 0; // normalisation constant
  double delta_en = 1;

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
    int ctp = findClassicalTurningPoint(en, v, pinf - d_ctp);
    // Check we won't go past any bounds because of d_ctp:
    // int d_ctp = d_ctp_in;
    // if (d_ctp >= pinf - ctp) {
    //   d_ctp = pinf - ctp - 1;
    //   if (d_ctp < 0)
    //     d_ctp = 0;
    // }

    // Find solution (f,g) to DE for given energy:
    // (Inward + outward solutions joined at ctp, merged over ctp+/-d_ctp)
    // Also stores dg (gout-gin) for PT
    std::vector<double> dg(2 * d_ctp + 1); // used for PT to find better e
    trialDiracSolution(f, g, dg, en, ka, v, r, drdt, h, ctp, d_ctp, pinf,
                       alpha);

    // Count the number of nodes (zeros) the wf has.
    int counted_nodes = countNodes(f, pinf);

    // If correct number of nodes, use PT to make minor energy adjustment.
    // Otherwise, make large adjustmunt until correct # of nodes
    if (counted_nodes == required_nodes) {
      correct_nodes = true;
      // do PT, AND find norm const:
      delta_en =
          smallEnergyChangePT(anorm, en, f, g, dg, pinf, ctp, d_ctp, alpha,
                              drdt, h, less, more, elower, eupper);
      // XXX - Update routine! Make const. Output de! ? or something..
      // XXX get 'anorm' should be sepperate!
      // XXX Npte: if never correct nodes, anorm not gotten!
    } else {
      correct_nodes = false; // can happen?
      bool more_nodes = (counted_nodes > required_nodes) ? true : false;
      delta_en = largeEnergyChange(en, more, less, eupper, elower, lfrac_de,
                                   more_nodes);
    }

    DEBUG(std::cerr << " :: it=" << its << " nodes:" << counted_nodes << "/"
                    << required_nodes << " new_en = " << en << " delta="
                    << delta_en * en << " eps=" << delta_en << "\n";
          std::cin.get();)

    if (delta_en < dele_goal && correct_nodes)
      break;
  } // end itterations

  DEBUG(if (!correct_nodes) {
    // Does this ever happen?? Yes, but rarely. Also: in fitParametric
    int counted_nodes = countNodes(f, pinf);
    std::cerr << "\nFAIL-148: wrong nodes:" << counted_nodes << "/"
              << required_nodes << " for n,k=" << n << "," << ka << "\n";
  })

  eps_out = delta_en;
  en_inout = en;
  pinf_out = pinf;
  its_out = its;

  // normalises the wavefunction
  double an = 1. / sqrt(anorm);
  for (int i = 0; i < pinf; i++) {
    f[i] = an * f[i];
    g[i] = an * g[i];
  }
  for (int i = pinf; i < ngp; i++) {
    // kills remainders (just for safety)
    f[i] = 0;
    g[i] = 0;
  }

  int ret_code = -1;
  if (correct_nodes && eps_out < dele_goal)
    ret_code = 0;
  else if (correct_nodes)
    ret_code = 1;
  else
    ret_code = 2;

  return ret_code;
}

//******************************************************************************
double largeEnergyChange(double &en, int &more, int &less, double &eupper,
                         double &elower, double lfrac_de, bool more_nodes)
/*
wf did not have correct number of nodes. Make a large energy adjustment
more_nodes=true means there were too many nodes
*/
{
  double en_old = en;

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

  return fabs((en_old - en) / en_old);
}

//******************************************************************************
double
smallEnergyChangePT(double &anorm, double &en, const std::vector<double> &f,
                    const std::vector<double> &g, const std::vector<double> &dg,
                    int pinf, int ctp, int d_ctp, double alpha,
                    const std::vector<double> &drdt, double h, int less,
                    int more, double elower, double eupper)
/*
Uses PT to calculate small change in energy.
Also calculates (+outputs) norm constant (but doesn't normalise orbital!)
*/
{
  double anormF = NumCalc::integrate(f, f, drdt, 1., 0, pinf);
  double anormG = NumCalc::integrate(g, g, drdt, 1., 0, pinf);
  anorm = (anormF + anormG) * h;

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
  double delta_en = fabs(de / en);
  double etemp = en + de;

  if ((less != 0) && (etemp < elower)) {
    delta_en = fabs((en - 0.5 * (en + elower)) / en);
    en = 0.5 * (en + elower);
  } else if ((more != 0) && (etemp > eupper)) {
    delta_en = fabs((en - 0.5 * (en + eupper)) / en);
    en = 0.5 * (en + eupper);
  } else if (etemp > 0) {
    // This only happens v. rarely. nodes correct, but P.T. gives silly result!
    // Is this OK? Seems to work
    if (de > 0)
      etemp = 0.9 * en;
    else
      etemp = 1.1 * en;
    delta_en = fabs((en - etemp) / en);
    en = etemp;
  } else {
    en = etemp;
  }

  return delta_en;
}

//******************************************************************************
int findPracticalInfinity(double en, const std::vector<double> &v,
                          const std::vector<double> &r, double alr)
/*
//Find the practical infinity 'pinf'
//Step backwards from the last point (ngp-1) until
// (V(r) - E)*r^2 >  alr    (alr = "asymptotically large r")
*/
{
  int ngp = (int)r.size();
  int pinf = ngp - 1;
  while ((en - v[pinf]) * pow(r[pinf], 2) + alr < 0)
    pinf--;

  DEBUG(if (pinf == ngp - 1) std::cerr
            << "WARNING 281: pract. inf. exceeds grid for en=" << en << "\n";)

  return pinf;
}

//******************************************************************************
int findClassicalTurningPoint(double en, const std::vector<double> &v,
                              int pinf) {
  // Finds classical turning point 'ctp'
  // Step backwards from the "practical infinity" until
  //  V(r) > E        [nb: both V and E are <0]
  int ctp = pinf;
  while ((en - v[ctp]) < 0) {
    --ctp;
    /*if (ctp<=0){
      //fails if ctp<0, (or ctp>pinf?)
      printf("FAILURE 96 in solveDBS: No classical region?\n");
    return 1;
    }*/
  }
  if (ctp >= pinf) {
    // Didn't find ctp! Does this ever happen? Yes, if energy guess too wrong
    DEBUG(printf("FAILURE: Turning point at or after pract. inf. \n");
          printf("ctp=%i, pinf=%i, ngp=%lu\n", ctp, pinf, v.size());)
    ctp = pinf - 5; //??? ok?
  }
  return ctp;
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
                        const std::vector<double> &drdt, double h, int ctp,
                        int d_ctp, int pinf, double alpha) {
  int ngp = (int)f.size();
  // Temporary vectors for in/out integrations:
  std::vector<double> pin(ngp), qin(ngp), pout(ngp), qout(ngp);
  // Perform the "inwards integration":
  inwardAM(pin, qin, en, v, ka, r, drdt, h, ctp - d_ctp, pinf, alpha);
  // Perform the "outwards integration"
  outwardAM(pout, qout, en, v, ka, r, drdt, h, ctp + d_ctp, alpha);
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

//******************************************************************************
int outwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
              const std::vector<double> &v, int ka,
              const std::vector<double> &r, const std::vector<double> &drdt,
              double h, int nf, double alpha)
/*
Program to start the OUTWARD integration.
Starts from 0, and uses an expansion(?) to go to (nol*AMO).
Then, it then call ADAMS-MOULTON, to finish (from nol*AMO+1 to nf = ctp+d_ctp)
*/
{
  const int nol = 1; // # of outdir runs [finds first nol*AMO+1 points (3)]

  double Z = -1 * v[AMO] * r[AMO]; //??? is this OK?

  double az = Z * alpha;
  double c2 = 1. / pow(alpha, 2);
  double ga = sqrt(pow(ka, 2) - pow(az, 2));

  // initial wf values
  // P(r) = r^gamma u(r)
  // Q(r) = r^gamma v(r)
  double u0 = 1;
  double v0;
  if (ka > 0)
    v0 = -(ga + ka) / az;
  else
    v0 = az / (ga - ka);
  p[0] = 0;
  q[0] = 0;

  // Coeficients used by the method
  std::vector<std::vector<double>> ie(AMO, std::vector<double>(AMO));
  std::vector<double> ia(AMO);
  double id;
  getOutwardCoefs(ie, ia, id);

  // loop through and find first nol*AMO points of wf
  for (int ln = 0; ln < nol; ln++) {
    int i0 = ln * AMO + 1;

    // defines/populates em coefs
    std::array<double, AMO> coefa, coefb, coefc, coefd;
    std::array<std::array<double, AMO>, AMO> em;
    for (int i = 0; i < AMO; i++) {
      double dror = drdt[i + i0] / r[i + i0];
      coefa[i] = (-id * h * (ga + ka) * dror);
      coefb[i] = (-id * h * (en + 2 * c2 - v[i + i0]) * drdt[i + i0] * alpha);
      coefc[i] = (id * h * (en - v[i + i0]) * drdt[i + i0] * alpha);
      coefd[i] = (-id * h * (ga - ka) * dror);
      for (int j = 0; j < AMO; j++)
        em[i][j] = ie[i][j];
      em[i][i] = em[i][i] - coefd[i];
    }
    // //inverts the em matrix
    em = Matrix::invert(em); // from here on, em is the inverted matrix

    // defines/populates fm, s coefs
    std::array<double, AMO> s;
    std::array<std::array<double, AMO>, AMO> fm;
    for (int i = 0; i < AMO; i++) {
      s[i] = -ia[i] * u0;
      for (int j = 0; j < AMO; j++) {
        fm[i][j] = ie[i][j] - coefb[i] * em[i][j] * coefc[j];
        s[i] = s[i] - coefb[i] * em[i][j] * ia[j] * v0;
      }
      fm[i][i] = fm[i][i] - coefa[i];
    }
    // inverts the matrix!  fm =-> Inv(fm)
    fm = Matrix::invert(fm); // from here on, fm is the inverted matrix

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
        vs[i] = vs[i] + em[i][j] * (coefc[j] * us[j] - ia[j] * v0);
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

  } // END for (int ln=0; ln<nol; ln++)  [loop through outint `nol' times]

  // Call adamsmoulton to finish integration from (nol*AMO+1) to nf = ctp+d_ctp
  int na = nol * AMO + 1;
  if (nf > na)
    adamsMoulton(p, q, en, v, ka, r, drdt, h, na, nf, alpha);

  return 0;
}

//******************************************************************
int inwardAM(std::vector<double> &p, std::vector<double> &q, double &en,
             const std::vector<double> &v, int ka, const std::vector<double> &r,
             const std::vector<double> &drdt, double h, int nf, int pinf,
             double alpha)
/*
Program to start the INWARD integration.
Starts from Pinf, and uses an expansion(?) to go to (pinf-AMO)
Then, it then call ADAMS-MOULTON, to finish (from nol*AMO+1 to nf = ctp-d_ctp)
*/
{
  // order of the expansion coeficients in 'inint'  (15 orig.)
  const int nx = 15;
  // PRIMARY convergance for expansion in `inint'' (10^-8)
  const double nxepsp = 1e-10;

  double alpha2 = pow(alpha, 2);
  double cc = 1. / alpha;
  double c2 = 1. / alpha2;

  double lambda = sqrt(-en * (2. + en * alpha2));
  double zeta = -v[pinf] * r[pinf];
  double sigma = (1. + en * alpha2) * (zeta / lambda);
  double Ren = en + c2; // total relativistic energy

  // Generates the expansion coeficients for asymptotic wf
  // up to order NX (nx is 'param')
  std::array<double, nx> bx;
  std::array<double, nx> ax;
  bx[0] = (ka + (zeta / lambda)) * (alpha / 2);
  for (int i = 0; i < nx; i++) {
    ax[i] = (ka + (i + 1 - sigma) * Ren * alpha2 - zeta * lambda * alpha2) *
            bx[i] * cc / ((i + 1) * lambda);
    if (i < (nx - 1))
      bx[i + 1] =
          (pow(ka, 2) - pow((i + 1 - sigma), 2) - pow(zeta, 2) * alpha2) *
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
    for (int k = 0; k < nx; k++) { // this will loop until a) converge, b) k=nx
      rk = rk * r[i];
      ps = ps + (ax[k] / rk);
      qs = qs + (bx[k] / rk);
      xe = fmax(fabs((ax[k] / rk) / ps), fabs((bx[k] / rk) / qs));
      if (xe < nxepsp)
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
    adamsMoulton(p, q, en, v, ka, r, drdt, h, pinf - AMO - 1, nf, alpha);

  return 0;
}

//******************************************************************************
int adamsMoulton(std::vector<double> &p, std::vector<double> &q, double &en,
                 const std::vector<double> &v, int ka,
                 const std::vector<double> &r, const std::vector<double> &drdt,
                 double h, int ni, int nf, double alpha)
/*
program finishes the INWARD/OUTWARD integrations (ADAMS-MOULTON)
  //- ni is starting (initial) point for integration
  //- nf is end (final) point for integration (nf=ctp+/-d_ctp)
*/
{
  double c2 = 1. / pow(alpha, 2); // c^2 - just to shorten code
  // AM coefs.
  std::vector<double> ama(AMO);
  double amd, amaa;
  getAdamsCoefs(ama, amd, amaa);

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
    return 1;
  }

  // create arrays for wf derivatives
  std::vector<double> dp(ngp), dq(ngp);
  std::array<double, AMO> amcoef;
  int k1 = ni - inc * AMO;
  for (int i = 0; i < AMO; i++) { // nb: k1 is iterated
    double dror = drdt[k1] / r[k1];
    dp[i] = inc * (-ka * dror * p[k1] -
                   alpha * ((en + 2 * c2) - v[k1]) * drdt[k1] * q[k1]);
    dq[i] = inc * (ka * dror * q[k1] + alpha * (en - v[k1]) * drdt[k1] * p[k1]);
    amcoef[i] = (h / amd) * ama[i];
    k1 += inc;
  }

  // integrates the function from ni to the c.t.p
  double a0 = (h / amd) * amaa;
  int k2 = ni;
  for (int i = 0; i < nosteps; i++) {
    double dror = drdt[k2] / r[k2];
    double dai = -inc * (ka * dror);
    double dbi = -inc * alpha * (en + 2 * c2 - v[k2]) * drdt[k2];
    double dci = inc * alpha * (en - v[k2]) * drdt[k2];
    double ddi = -dai;
    double det = 1 - a0 * a0 * (dbi * dci - dai * ddi);
    double sp = p[k2 - inc];
    double sq = q[k2 - inc];
    for (int l = 0; l < AMO; l++) {
      sp = sp + amcoef[l] * dp[l];
      sq = sq + amcoef[l] * dq[l];
    }
    p[k2] = (sp + a0 * (dbi * sq - ddi * sp)) / det;
    q[k2] = (sq + a0 * (dci * sp - dai * sq)) / det;
    for (int l = 0; l < (AMO - 1); l++) { // loads next 'first' k values (?)
      dp[l] = dp[l + 1];
      dq[l] = dq[l + 1];
    }
    dp[AMO - 1] = dai * p[k2] + dbi * q[k2]; // loads next 'first' deriv's (?)
    dq[AMO - 1] = dci * p[k2] + ddi * q[k2];
    k2 += inc;
  }

  return 0;
} // END adamsmoulton

//**************************************************************************
//**************************************************************************
//**************************************************************************
//**************************************************************************

//******************************************************************
int getAdamsCoefs(std::vector<double> &mia, double &mid, double &miaa)
/*
coeficients for the ADAMS-MOULTON routine
*/
{
  // XXX XXX XXX Make array of MAX domension! then, just 'input' the correct #

  if (AMO == 8) {
    double tia[8] = {-33953,   312874,  -1291214, 3146338,
                     -5033120, 5595358, -4604594, 4467094};
    mid = 3628800;
    miaa = 1070017;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 7) {
    double tia[7] = {1375, -11351, 41499, -88547, 123133, -121797, 139849};
    mid = 120960;
    miaa = 36799;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 6) {
    double tia[6] = {-863, 6312, -20211, 37504, -46461, 65112};
    mid = 60480;
    miaa = 19087;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 5) {
    double tia[5] = {27, -173, 482, -798, 1427};
    mid = 1440;
    miaa = 475;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 4) {
    double tia[4] = {-19, 106, -264, 646};
    mid = 720;
    miaa = 251;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 3) {
    double tia[3] = {1, -5, 19};
    mid = 24;
    miaa = 9;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 2) {
    double tia[2] = {-1, 8};
    mid = 12;
    miaa = 5;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else if (AMO == 1) {
    double tia[1] = {1};
    mid = 2;
    miaa = 1;
    for (int i = 0; i < AMO; i++) {
      mia[i] = tia[i];
    }
  } else {
    printf("FAILURE: No Adams-Moulton coeficients. Check AMO\n");
    return 1;
  }

  return 0;
} // END AMcoefs

//******************************************************************
int getOutwardCoefs(std::vector<std::vector<double>> &oie,
                    std::vector<double> &oia, double &oid)
/*
coeficients for the OUTINT routine
*/
{
  // XXX XXX XXX Make array of MAX domension! then, just 'input' the correct #
  // ?? little harder here..
  // XX ? only down to 5?

  if (AMO == 8) {
    double tie[8][8] = {{-1338, 2940, -2940, 2450, -1470, 588, -140, 15},
                        {-240, -798, 1680, -1050, 560, -210, 48, -5},
                        {60, -420, -378, 1050, -420, 140, -30, 3},
                        {-32, 168, -672, 0, 672, -168, 32, -3},
                        {30, -140, 420, -1050, 378, 420, -60, 5},
                        {-48, 210, -560, 1050, -1680, 798, 240, -15},
                        {140, -588, 1470, -2450, 2940, -2940, 1338, 105},
                        {-960, 3920, -9408, 14700, -15680, 11760, -6720, 2283}};
    double tia[8] = {-105, 15, -5, 3, -3, 5, -15, 105};
    oid = 840;
    for (int i = 0; i < AMO; i++) {
      oia[i] = tia[i];
      // oia.push_back(tia[i]);
      for (int j = 0; j < AMO; j++) {
        oie[i][j] = tie[i][j];
      }
    }
  } else if (AMO == 7) {
    double tie[7][7] = {{-609, 1260, -1050, 700, -315, 84, -10},
                        {-140, -329, 700, -350, 140, -35, 4},
                        {42, -252, -105, 420, -126, 28, -3},
                        {-28, 126, -420, 105, 252, -42, 4},
                        {35, -140, 350, -700, 329, 140, -10},
                        {-84, 315, -700, 1050, -1260, 609, 60},
                        {490, -1764, 3675, -4900, 4410, -2940, 1089}};
    double tia[7] = {-60, 10, -4, 3, -4, 10, -60};
    oid = 420;
    for (int i = 0; i < AMO; i++) {
      oia[i] = tia[i];
      for (int j = 0; j < AMO; j++) {
        oie[i][j] = tie[i][j];
      }
    }
  } else if (AMO == 6) {
    double tie[6][6] = {
        {-77, 150, -100, 50, -15, 2}, {-24, -35, 80, -30, 8, -1},
        {9, -45, 0, 45, -9, 1},       {-8, 30, -80, 35, 24, -2},
        {15, -50, 100, -150, 77, 10}, {-72, 225, -400, 450, -360, 147}};
    double tia[6] = {-10, 2, -1, 1, -2, 10};
    oid = 60;
    for (int i = 0; i < AMO; i++) {
      oia[i] = tia[i];
      for (int j = 0; j < AMO; j++) {
        oie[i][j] = tie[i][j];
      }
    }
  } else if (AMO == 5) {
    double tie[5][5] = {{-65, 120, -60, 20, -3},
                        {-30, -20, 60, -15, 2},
                        {15, -60, 20, 30, -3},
                        {-20, 60, -120, 65, 12},
                        {75, -200, 300, -300, 137}};
    double tia[5] = {-12, 3, -2, 3, -12};
    oid = 60;
    for (int i = 0; i < AMO; i++) {
      oia[i] = tia[i];
      for (int j = 0; j < AMO; j++) {
        oie[i][j] = tie[i][j];
      }
    }
  } else {
    printf("FAILURE: No Adams-Moulton (OUTINT) coeficients. Check AMO\n");
    return 1;
  }

  return 0;
} //   END OIcoefs

} // namespace ADAMS
