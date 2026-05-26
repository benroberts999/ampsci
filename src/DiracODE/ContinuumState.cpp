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
  // Require at least 50 points per wavelength (25 per half-wave)
  auto dr = std::min(dr0, approx_wavelength / 50.0);

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
  Fa *= (amp != 0.0 ? (D / amp) : 0.0);
}

//==============================================================================
std::pair<double, double> numerical_f_amplitude(double en, int kappa,
                                                double alpha, double Zeff,
                                                double f_final, double g_final,
                                                double r_final, double dr) {

  // In the asymptotic region (V negligible) the Dirac ODE reduces to:
  //   f'' = -k^2 f,  k^2 = en * (alpha^2 * en + 2)
  // giving f = A sin(kr+phi), g = A*rho*cos(kr+phi),
  // where rho = sqrt(alpha^2*en / (alpha^2*en + 2)).
  // Combining: A^2 = f^2 + c_g^2 * g^2,  c_g^2 = 1 + 2/(alpha^2 * en).
  // This gives A at every ODE step without requiring peak or node detection.
  const double c_g2 = 1.0 + 2.0 / (alpha * alpha * en);

  DiracContinuumDerivative Heff(Zeff, kappa, en, alpha);
  AdamsMoulton::ODESolver2D<Param::K_Adams, double, double> ode{dr, &Heff};
  ode.solve_initial_K(r_final, f_final, g_final);

  const double eps_target = 1.0e-6;
  // Relative to r_final: for high energy, 5000/sqrt(2*en) << r_final

  const double approx_wavelength = 2.0 * M_PI / std::sqrt(2.0 * en);
  const double max_r = r_final + 250 * approx_wavelength;

  // Average A over windows of ~one wavelength (dr ~ lambda/40, so N=40 steps).
  // Compare consecutive window averages for convergence.
  const int N_window = 40;
  double amp = 0.0;
  double eps_amp = 1.0;
  double window_sum = 0.0;
  int window_count = 0;
  double amp_prev = 0.0;
  bool have_prev = false;
  int count_num_windows = 0;

  while (ode.last_t() < max_r) {
    ode.drive();
    const double f = ode.last_f();
    const double g = ode.last_g();
    window_sum += std::sqrt(f * f + c_g2 * g * g);
    ++window_count;

    if (window_count == N_window) {
      const double amp_curr = window_sum / N_window;
      if (have_prev) {
        eps_amp = std::abs(amp_curr - amp_prev) / (0.5 * (amp_curr + amp_prev));
        amp = amp_curr;
        if (eps_amp < eps_target && count_num_windows > 5)
          break;
      }
      amp = amp_curr;
      amp_prev = amp_curr;
      have_prev = true;
      window_sum = 0.0;
      window_count = 0;
      ++count_num_windows;
    }
  }

  if (eps_amp > 1.0e-2)
    amp = 0.0;

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
