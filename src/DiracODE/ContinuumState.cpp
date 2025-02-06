#include "DiracODE/ContinuumState.hpp"
#include "DiracODE/BoundState.hpp"
#include "DiracODE/include.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include <algorithm>
#include <cmath>

namespace DiracODE {
using namespace DiracODE::Internal;

//==============================================================================
void solveContinuum(DiracSpinor &Fa, double en, const std::vector<double> &v,
                    double alpha, const DiracSpinor *const VxFa,
                    const DiracSpinor *const Fa0) {

  Fa.en() = en;
  const auto &gr = Fa.grid();
  const auto num_points = gr.num_points();
  Fa.max_pt() = num_points;

  // Rough expression for wavelenth at large r
  // nb: sin(kr + \eta * log(kr)), so not exactly constant
  const double approx_wavelength = 2.0 * M_PI / std::sqrt(2.0 * Fa.en());

  // The solution on the radial grid must be reasonable at large r
  // We ensure there is at least N (N=10) points per wavelength in this region
  // If not, solution unstable; write zeros and return.
  const int N_ppw = 10;
  const auto dr0 = gr.drdu().back() * gr.du();
  const double dr0_target = approx_wavelength / N_ppw;
  if (dr0 > dr0_target) {
    fmt2::styled_print(fg(fmt::color::red), "\nERROR 104: ");
    fmt::print(
        "Grid not dense enough for continuum state with e={:.2f} (kappa={}); \n"
        "Try increasing points (to ~ du < {:.3f})\n"
        "Writing zeros to Spinor for this state.\n",
        en, Fa.kappa(), dr0_target);
    Fa *= 0.0;
    return;
  }

  // Solve on regular grid - not yet normalised:
  DiracDerivative Hd(Fa.grid(), v, Fa.kappa(), Fa.en(), alpha, {}, VxFa, Fa0);
  solve_Dirac_outwards(Fa.f(), Fa.g(), Hd);

  // Now, normalise the solution.
  // Keep solving ODE outwards, on linearly-spaced grid
  // Use "H-like" derivative: assume exchange etc. negligable here
  // Calculate running:
  //    - wavelength
  //    - amplitude
  // Until they stabilise; find amplitude in this region
  // Re-scale wavefunction so tha large-r amplitude matches analytic expression

  // Step-size for large-r solution: Uses linear grid
  // Require at least 40 points per wavelength (20 per half-wave)
  // but limit to ~100 pts per half-wave (? no benefit after this)
  auto dr = std::min(dr0, approx_wavelength / 40.0);
  dr = std::max(dr, approx_wavelength / 200.0);

  // Parameters needed to solve f(r) to large r:
  const auto ztmp = -1.0 * v.back() * gr.r().back();
  const auto Zeff = std::max(1.0, ztmp);
  const auto r_final = gr.r().back();
  const auto f_final = Fa.f().back();
  const auto g_final = Fa.g().back();

  // Find large-r amplitude:
  const auto [amp, eps_amp] = numerical_f_amplitude(
      en, Fa.kappa(), alpha, Zeff, f_final, g_final, r_final, dr);

  Fa.max_pt() = num_points;
  Fa.eps() = eps_amp;

  // Calculate normalisation coeficient, D, and re-scaling factor:
  const auto D = analytic_f_amplitude(en, alpha);
  Fa *= (D / amp);
}

//==============================================================================
std::pair<double, double> numerical_f_amplitude(double en, int kappa,
                                                double alpha, double Zeff,
                                                double f_final, double g_final,
                                                double r_final, double dr) {

  // Now, normalise the solution.
  // Keep solving ODE outwards, on linearly-spaced grid
  // Use "H-like" derivative: assume exchange etc. negligable here
  // Calculate running:
  //    - wavelength
  //    - amplitude
  // Until they stabilise; find amplitude in this region

  // Set-up H-like solver for linearly-spaced grid
  DiracContinuumDerivative Heff(Zeff, kappa, en, alpha);
  AdamsMoulton::ODESolver2D<Param::K_Adams, double, double> ode{dr, &Heff};
  // start from existing final position:
  ode.solve_initial_K(r_final, f_final, g_final);

  // convergence targets, for amplitude and wavelength
  const double eps_target_amp = 1.0e-7;
  const double eps_target_lambda = 2.0e-5;

  // amplitude:
  double amp = 0.0;
  double eps_amp = 1.0;
  double amp_maxf = 0.0;

  // wavelength
  double lambda = 0.0;
  double eps_lambda = 1.0;
  double r_node_0 = 0.0;

  // previous points:
  double rm1 = 0.0;
  double fm1 = 0.0;
  double r0 = ode.last_t();
  double f0 = ode.last_f();

  // Reasonable guess at fully asymptotic region:
  // (Only reaches max_r if amp/wavelength don't converge)
  const double max_r = 5000.0 / std::sqrt(2 * en);

  while (ode.last_t() < max_r) {
    ode.drive();
    const auto r1 = ode.last_t();
    const auto f1 = ode.last_f();

    // Update wavelength:
    if (f1 * f0 < 0.0) {
      // linear interp: find rx (position of node)
      const double r_node_1 = (r1 * f0 - r0 * f1) / (f0 - f1);
      const double lambda1 = 2.0 * (r_node_1 - r_node_0);
      eps_lambda = std::abs(2.0 * (lambda1 - lambda) / (lambda1 + lambda));
      lambda = lambda1;
      r_node_0 = r_node_1;
    }

    // find amplitude via maximum point:
    if (std::abs(f1) > amp_maxf) {
      amp_maxf = std::abs(f1);
    }

    // find amplitude by fitting around maximum:
    if (std::abs(f1) < std::abs(f0) && std::abs(f0) > std::abs(fm1)) {

      const auto amp_fit = fitQuadratic(rm1, r0, r1, fm1, f0, f1);
      const auto eps_amp_fit =
          std::abs(2.0 * (amp_fit - amp) / (amp_fit + amp));
      amp = amp_fit;

      const auto eps_amp_max =
          std::abs(2.0 * (amp_fit - amp_maxf) / (amp_fit + amp_maxf));

      const auto eps_amp_best = std::min({eps_amp_fit, eps_amp_max});
      const auto eps_amp_worst = std::max({eps_amp_fit, eps_amp_max});
      eps_amp = eps_amp_worst;

      if (eps_amp_best < eps_target_amp &&
          eps_amp_worst < 100.0 * eps_target_amp &&
          eps_lambda < eps_target_lambda) {
        amp = std::max(amp_fit, amp_maxf);
        break;
      }
    }
    // update "previous" values
    rm1 = r0;
    fm1 = f0;
    r0 = r1;
    f0 = f1;
  }
  assert(amp != 0.0);
  return {amp, eps_amp};
}

//==============================================================================
double analytic_f_amplitude(double en, double alpha) {
  // D = Sqrt[alpha/(pi*eps)] <-- Amplitude of large-r f(r)
  // eps = Sqrt[en/(en+2mc^2)]
  // ceps = c*eps = eps/alpha
  const double al2 = std::pow(alpha, 2);
  const double ceps = std::sqrt(en / (en * al2 + 2.0));
  return 1.0 / std::sqrt(M_PI * ceps);
}

//==============================================================================
double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3)
// Takes in three points, and fits them to a quadratic function.
// Returns y-value for vertex of quadratic.
// Used for finding the amplitude of a sine/cosine function, given thee
// points. i.e., will return amplitude of since function. Note: the given 3
// points _MUST_ be close to maximum, otherwise, fit wont work
{
  if (y1 < 0)
    y1 = std::abs(y1);
  if (y2 < 0)
    y2 = std::abs(y2);
  if (y3 < 0)
    y3 = std::abs(y3);

  const auto d = (x1 - x2) * (x1 - x3) * (x2 - x3);
  const auto Ad = x3 * (x2 * (x2 - x3) * y1 + x1 * (-x1 + x3) * y2) +
                  x1 * (x1 - x2) * x2 * y3;
  const auto Bd =
      x3 * x3 * (y1 - y2) + x1 * x1 * (y2 - y3) + x2 * x2 * (-y1 + y3);
  const auto Cd = x3 * (-y1 + y2) + x2 * (y1 - y3) + x1 * (-y2 + y3);
  auto y0 = (Ad / d) - Bd * Bd / (4.0 * Cd * d);

  // Find largest input y:
  const auto ymax = std::max({y1, y2, y3});
  // if (y1 > ymax)
  //   ymax = y1;
  // if (y3 > ymax)
  //   ymax = y3;

  if (ymax > y0)
    y0 = ymax; // y0 can't be less than (y1,y2,y3)

  return y0;
}

} // namespace DiracODE
