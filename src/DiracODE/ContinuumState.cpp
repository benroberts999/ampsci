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

  // Rough expression for wavelength at large r
  // nb: sin(kr + \eta * log(kr)), so not exactly constant
  const double approx_wavelength = 2.0 * M_PI / std::sqrt(2.0 * Fa.en());

  // The solution on the radial grid must be reasonable at large r
  // We ensure there is at least N (N=10) points per wavelength in this region
  // If not, solution unstable; write zeros and return.
  const int N_ppw = 15;
  const auto dr0 = gr.drdu().back() * gr.du();
  const double dr0_target = approx_wavelength / N_ppw;
  if (dr0 > dr0_target) {
    // fmt2::styled_print(fg(fmt::color::red), "\nERROR 104: ");
    // fmt::print(
    //   "Grid not dense enough for continuum state with e={:.2f} (kappa={}); \n"
    //   "Try increasing points (to ~ du < {:.3f})\n"
    //   "Writing zeros to Spinor for this state.\n",
    //   en, Fa.kappa(), dr0_target);
    Fa *= 0.0;
    return;
  }

  // Solve on regular grid - not yet normalised:
  DiracDerivative Hd(Fa.grid(), v, Fa.kappa(), Fa.en(), alpha, {}, VxFa, Fa0);
  solve_Dirac_outwards(Fa.f(), Fa.g(), Hd);

  // Now, normalise the solution.
  // Keep solving ODE outwards, on linearly-spaced grid
  // Use "H-like" derivative: assume exchange etc. negligible here
  // Calculate running:
  //    - wavelength
  //    - amplitude
  // Until they stabilise; find amplitude in this region
  // Re-scale wavefunction so that large-r amplitude matches analytic expression

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

  Fa.eps() = eps_amp;

  // Calculate normalisation coefficient, D, and re-scaling factor:
  const auto D = analytic_f_amplitude(en, alpha);
  Fa *= (amp != 0.0 ? (D / amp) : 0.0);
}

//==============================================================================
std::pair<double, double> numerical_f_amplitude(double en, int kappa,
                                                double alpha, double Zeff,
                                                double f_final, double g_final,
                                                double r_final, double dr) {

  // Continue solving on a linearly-spaced grid using an H-like potential
  // (exchange negligible at large r). Track wavelength and amplitude
  // half-cycle by half-cycle until both converge.

  // Set-up H-like solver for linearly-spaced grid
  DiracContinuumDerivative Heff(Zeff, kappa, en, alpha);
  AdamsMoulton::ODESolver2D<Param::K_Adams, double, double> ode{dr, &Heff};
  // start from existing final position:
  ode.solve_initial_K(r_final, f_final, g_final);

  // Convergence targets for amplitude and wavelength
  const double eps_target_amp = 1.0e-7;
  const double eps_target_lambda = 2.0e-5;

  // Asymptotic form of the large component:
  //   f(r) → D * sin(k*r + η*log(2kr) + δ_κ)
  // where k = sqrt(en*(en+2/α²)) and η = -Z*α/(k*α²) is the Coulomb parameter.
  // At large r, both the wavelength λ = 2π/k and envelope amplitude D converge
  // to constants.
  //
  // Strategy: track λ and max|f| half-cycle by half-cycle.
  //   Node positions r_i satisfy k*r_i + ... = i*π  ⟹  r_{i+1} - r_i ≈ λ/2
  //   ⟹  λ_i = 2*(r_{i+1} - r_i)
  //   A_i = max|f| between node i and node i+1 (one half-period)
  // Converged when both λ and A change by < eps between consecutive half-cycles.

  // Wavelength tracking.
  // We need two consecutive node positions to form λ = 2*(r_{i+1} - r_i),
  // and a third to compute eps_lambda = |λ_{i+1} - λ_i| / <λ>.
  // n_nodes counts crossings seen so far; λ and eps_lambda are valid only
  // once n_nodes >= 2 (third crossing).
  int n_nodes = 0;
  double r_node_prv = 0.0; // position of most-recently seen node
  double lambda = 0.0;     // most recent half-cycle estimate of λ
  double eps_lambda = 1.0; // |Δλ/λ| between successive half-cycles

  // Per-half-cycle amplitude tracking.
  // A half-cycle is the interval between two consecutive nodes.
  // amp_this_half = running max|f| in the current half-cycle (resets at each node)
  // amp_prev_half = max|f| in the completed half-cycle before the latest node
  // eps_amp       = |A_this - A_prev| / ½(A_this + A_prev)
  //
  // The half-cycle max is accurate to O(dr/λ)² ≈ 0.06% at 40 pts/wavelength,
  // which is more than sufficient for normalisation purposes.
  double amp_this_half = std::abs(ode.last_f()); // include starting point
  double amp_prev_half = 0.0;
  double eps_amp = 1.0;

  double r0 = ode.last_t();
  double f0 = ode.last_f();

  // Fallback limit: a fully asymptotic region is expected within O(1000)
  // de Broglie wavelengths of the grid boundary.
  const double max_r = 5000.0 / std::sqrt(2.0 * en);

  while (ode.last_t() < max_r) {
    ode.drive();
    const double r1 = ode.last_t();
    const double f1 = ode.last_f();

    // Track running half-cycle amplitude
    amp_this_half = std::max(amp_this_half, std::abs(f1));

    // Node crossing: f changes sign between r0 and r1
    if (f1 * f0 < 0.0) {
      // Linear interpolation for the node position r_node where f(r_node) = 0:
      //   f(r) ≈ f0 + (f1-f0)/(r1-r0) * (r - r0) = 0
      //   ⟹  r_node = (r1*f0 - r0*f1) / (f0 - f1)
      const double r_node = (r1 * f0 - r0 * f1) / (f0 - f1);

      if (n_nodes >= 1) {
        // λ_i = 2*(r_node_i - r_node_{i-1})  [distance between adjacent nodes = λ/2]
        const double lambda1 = 2.0 * (r_node - r_node_prv);

        if (n_nodes >= 2) {
          // Third and later nodes: proper relative-change convergence check.
          // eps_λ = |λ_new - λ_old| / ½(λ_new + λ_old)
          eps_lambda = std::abs(2.0 * (lambda1 - lambda) / (lambda1 + lambda));
        }
        // (On the second node, lambda was not yet set; just initialise it.)
        lambda = lambda1;

        // Half-cycle amplitude convergence:
        //   eps_A = |A_this - A_prev| / ½(A_this + A_prev)
        if (amp_prev_half > 0.0) {
          eps_amp = std::abs(2.0 * (amp_this_half - amp_prev_half) /
                             (amp_this_half + amp_prev_half));
          if (eps_amp < eps_target_amp && eps_lambda < eps_target_lambda)
            break;
        }
      }

      n_nodes++;
      r_node_prv = r_node;
      amp_prev_half = amp_this_half;
      amp_this_half = 0.0; // reset for the next half-cycle
    }

    r0 = r1;
    f0 = f1;
  }

  // At break, amp_this_half is the just-completed half-cycle; amp_prev_half
  // is the one before it. Both have converged — use the most recent.
  // If the loop ran to max_r without converging, take whichever is larger.
  const double amp = std::max(amp_prev_half, amp_this_half);
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
                    double y3) {

  // Fit p(x) = A*x² + B*x + C through three points and return the vertex
  // amplitude y* = C - B²/(4A).
  //
  // Previously used by numerical_f_amplitude to refine the per-half-cycle peak
  // estimate: given the grid point with the largest |f| and its two neighbours,
  // the quadratic vertex corrects for the O(dr/λ)² discretisation bias.
  // Replaced by the plain half-cycle max, which is already accurate to
  // O(dr/λ)² ≈ 0.06% at 40 pts/wavelength and sufficient for normalisation;
  // the quadratic refinement added complexity without measurable benefit.
  //
  // The y inputs are taken as |y| — sign doesn't matter for the envelope.
  // Coefficients via Lagrange, writing d = (x1-x2)(x1-x3)(x2-x3):
  //   A·d ≡ Cd = x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)
  //   B·d ≡ Bd = x3²*(y1-y2) + x1²*(y2-y3) + x2²*(y3-y1)
  //   C·d ≡ Ad = x3*x2*(x2-x3)*y1 + x1*x3*(x3-x1)*y2 + x1*x2*(x1-x2)*y3
  //   y*  = C - B²/(4A) = Ad/d - Bd²/(4·Cd·d)
  //
  // A < 0 (downward parabola, genuine maximum) iff Cd·d < 0.
  // If A ≥ 0 the points don't bracket a peak; return the largest input value.

  y1 = std::abs(y1);
  y2 = std::abs(y2);
  y3 = std::abs(y3);

  const double d = (x1 - x2) * (x1 - x3) * (x2 - x3);
  const double Ad =
    x3 * (x2 * (x2 - x3) * y1 + x1 * (x3 - x1) * y2) + x1 * (x1 - x2) * x2 * y3;
  const double Bd =
    x3 * x3 * (y1 - y2) + x1 * x1 * (y2 - y3) + x2 * x2 * (y3 - y1);
  const double Cd = x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2); // = A·d

  // Guard: only trust the vertex when the parabola opens downward (A < 0)
  const double ymax = std::max({y1, y2, y3});
  if (Cd * d >= 0.0)
    return ymax;

  return Ad / d - Bd * Bd / (4.0 * Cd * d);
}

} // namespace DiracODE
