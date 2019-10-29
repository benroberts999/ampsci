#include "Adams_continuum.hpp"
#include "Adams_bound.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "DiracODE.hpp"
#include "Maths/Grid.hpp"
#include <cmath>

namespace DiracODE {
using namespace Adams;

// Program to solve single-electron continuum-state Dirac problem for a
// (given) local, central potential. Uses "outwardAM" (from
// adamsSolveLocalBS.cpp) to solve dirac equation
//
// Normalises wavefunction by solving all the way out to asymptotic region,
// where solution should be sinosoidal. Then, matches amplitude with low-r
// expansion. Therefore, need grid to go very far out (v. large r). Also, need
// at least ~10 points per half-period. More points for higher energy!
//
// ###== To do ###==
//   * Find asymptotic region + normalise...better?

//******************************************************************************
int solveContinuum(DiracSpinor &phi, const std::vector<double> &v,
                   const Grid &ext_grid, std::size_t i_asym, double alpha)
// Solves Dirac equation for continuum state, for given energy, ec
// by integrating outwards from 0
// ec > 0
// Normalises wf. by comparison w/ H-like solution at v. large r
// num_pointsb id the regular (bound-state) grid.
// num_pointsc is grid for continuum (only for solving). num_pointsc >>
// num_pointsb
{
  auto num_pointsb = phi.p_rgrid->num_points;
  auto num_pointsc = ext_grid.num_points;

  // Perform the "outwards integration"
  DiracSpinor psic(phi.n, phi.k, ext_grid);
  psic.en = phi.en;

  DiracMatrix Hd(ext_grid, v, psic.k, psic.en, alpha);
  outwardAM(psic.f, psic.g, Hd, (int)num_pointsc - 1);

  // Find a better (lower) asymptotic region:
  i_asym = findAsymptoticRegion(psic.f, ext_grid.r, num_pointsb, num_pointsc,
                                i_asym);

  // Find amplitude of large-r (asymptotic region) sine-like wf
  double amp = findSineAmplitude(psic.f, ext_grid.r, num_pointsc, i_asym);

  // Calculate normalisation coeficient, D, and re-scaling factor:
  // D = Sqrt[alpha/(pi*eps)] <-- Amplitude of large-r p(r)
  // eps = Sqrt[en/(en+2mc^2)]
  double al2 = std::pow(alpha, 2);
  double ceps = std::sqrt(psic.en / (psic.en * al2 + 2.)); // c*eps = eps/alpha
  double D = 1. / std::sqrt(M_PI * ceps);
  double sf = D / amp; // re-scale factor

  // Normalise the wfs, and transfer back to shorter arrays:
  for (std::size_t i = 0; i < num_pointsb; i++) {
    phi.f[i] = sf * psic.f[i];
    phi.g[i] = sf * psic.g[i];
  }

  return 0;
}

namespace Adams {
//******************************************************************************
double findSineAmplitude(std::vector<double> &pc, const std::vector<double> &rc,
                         std::size_t num_pointsc, std::size_t i_asym)
//  Find "maximum" amplitude, by using a quadratic fit to 2 nearest points
//  Scale by ratio of this maximum to max of analytic soln
{
  int ntry = 0, maxtry = 5;
  double amp = 0;
  while (ntry < maxtry) {
    // find first zero after r_asym
    for (std::size_t i = i_asym; i < num_pointsc; i++) {
      if (pc[i] * pc[i - 1] < 0) {
        i_asym = i;
        break;
      }
    }
    // find max:
    double y0 = 0, y1 = 0, y2 = 0, y3 = 0, y4 = 0;
    double x0 = 0, x1 = 0, x2 = 0, x3 = 0, x4 = 0;
    for (std::size_t i = i_asym + 1; i < num_pointsc - 1; i++) {
      if (std::fabs(pc[i]) < std::fabs(pc[i - 1])) {
        y0 = std::fabs(pc[i - 3]);
        y1 = std::fabs(pc[i - 2]);
        y2 = std::fabs(pc[i - 1]);
        y3 = std::fabs(pc[i]);
        y4 = std::fabs(pc[i + 1]);
        x0 = rc[i - 3];
        x1 = rc[i - 2];
        x2 = rc[i - 1];
        x3 = rc[i];
        x4 = rc[i + 1];
        break;
      }
    }
    i_asym++;
    ntry++;
    double out1 = fitQuadratic(x1, x2, x3, y1, y2, y3);
    double out2 = fitQuadratic(x0, x2, x4, y0, y2, y4);
    amp += 0.5 * (out1 + out2);
  }

  return amp /= maxtry;
}

//******************************************************************************
std::size_t findAsymptoticRegion(std::vector<double> &pc,
                                 const std::vector<double> &rc,
                                 std::size_t num_pointsb,
                                 std::size_t num_pointsc, std::size_t i_asym)
// Finds a 'better' guess for the asymptotic region, by looking for where
// the period of oscilations becomes constant
// (variation in periof drops below certain value)
//
// Note: this method works well, but not perfectly.
//  a) am I not going out far enough?
//  b) Or, it's just not that accurate? Seems acurate to 1 - 0.1% ?
{
  // Find the r's for psi=zero, two consec => period
  // Once period is converged enough, can normalise by comparison with
  // exact (asymptotic) solution (??)
  // nb: looks for convergence between r(num_points) and r_asym
  // If doesn't 'converge' in this region, uses r_asym
  double xa = 1, xb = pc[num_pointsb]; // num_pointsb is #pts in regular grid
  double wk1 = -1, wk2 = 0;
  for (std::size_t i = num_pointsb; i < i_asym; i++) {
    xa = xb;
    xb = pc[i];
    if (xb * xa < 0) {
      // Use linear extrapolation to find exact r of the zero:
      double r1 = (rc[i] * pc[i - 1] - rc[i - 1] * pc[i]) / (pc[i - 1] - pc[i]);
      // double ya = xb,
      auto yb = xb;
      for (std::size_t j = i + 1; j < num_pointsc; j++) {
        auto ya = yb;
        yb = pc[j];
        if (ya * yb < 0) {
          // Use linear extrapolation to find exact r of the next zero:
          double r2 =
              (rc[j] * pc[j - 1] - rc[j - 1] * pc[j]) / (pc[j - 1] - pc[j]);
          wk1 = wk2;
          wk2 = r2 - r1;
          break;
        }
      }
      if (std::fabs(wk1 - wk2) < 1.e-4) {
        // check for "convergence"
        i_asym = i;
        break;
      }
    }
  }
  return i_asym;
}

//******************************************************************************
double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3)
// Takes in three points, and fits them to a quadratic function.
// Returns y-value for vertex of quadratic.
// Used for finding the amplitude of a sine/cosine function, given thee
// points. i.e., will return amplitude of since function. Note: the given 3
// points _MUST_ be close to maximum, otherwise, fit wont work
{
  if (y1 < 0)
    y1 = std::fabs(y1);
  if (y2 < 0)
    y2 = std::fabs(y2);
  if (y3 < 0)
    y3 = std::fabs(y3);

  double d = (x1 - x2) * (x1 - x3) * (x2 - x3);
  double Ad = x3 * (x2 * (x2 - x3) * y1 + x1 * (-x1 + x3) * y2) +
              x1 * (x1 - x2) * x2 * y3;
  double Bd = x3 * x3 * (y1 - y2) + x1 * x1 * (y2 - y3) + x2 * x2 * (-y1 + y3);
  double Cd = x3 * (-y1 + y2) + x2 * (y1 - y3) + x1 * (-y2 + y3);
  double y0 = (Ad / d) - Bd * Bd / (4. * Cd * d);

  // Find largest input y:
  double ymax = y2;
  if (y1 > ymax)
    ymax = y1;
  if (y3 > ymax)
    ymax = y3;

  if (ymax > y0)
    y0 = ymax; // y0 can't be less than (y1,y2,y3)

  return y0;
}

} // namespace Adams
} // namespace DiracODE
