#include "DiracODE/ContinuumState.hpp"
#include "DiracODE/AsymptoticSpinor.hpp"
#include "DiracODE/BoundState.hpp"
#include "DiracODE/include.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
#include <algorithm>
#include <cmath>

namespace DiracODE {
using namespace DiracODE::Internal;

namespace {
// Dirac derivative matrix for a local potential v(r), linearly interpolated
// off the radial grid. Used to solve through the region where the radial
// grid under-resolves the continuum oscillations (fine fixed-step grid).
struct InterpPotentialDerivative
  : AdamsMoulton::DerivativeMatrix<double, double> {

  InterpPotentialDerivative(const Grid &in_gr, const std::vector<double> &in_v,
                            int in_kappa, double in_en, double in_alpha)
    : gr(in_gr),
      pv(in_v),
      kappa(in_kappa),
      en(in_en),
      alpha(in_alpha),
      cc(1.0 / alpha) {}

  const Grid &gr;
  const std::vector<double> &pv;
  int kappa;
  double en, alpha, cc;

  double v(double r) const {
    const auto i0 = gr.getIndex(r);
    const auto i1 = std::min(i0 + 1, gr.num_points() - 1);
    if (i1 == i0) {
      return pv[i0];
    }
    const auto t = (r - gr.r(i0)) / (gr.r(i1) - gr.r(i0));
    return (1.0 - t) * pv[i0] + t * pv[i1];
  }
  double a(double r) const final { return double(-kappa) / r; }
  double b(double r) const final { return alpha * (en - v(r)) + 2.0 * cc; }
  double c(double r) const final { return -alpha * (en - v(r)); }
  double d(double r) const final { return -a(r); }
};
} // namespace

//==============================================================================
void solveContinuum(DiracSpinor &Fa, double en, const std::vector<double> &v,
                    double alpha, const DiracSpinor *const VxFa,
                    const DiracSpinor *const Fa0, bool average_tail) {

  Fa.en() = en;
  const auto &gr = Fa.grid();
  const auto num_points = gr.num_points();

  // Rough expression for wavelenth at large r:
  // nb: sin(kr + nu * log(kr)), so not exactly constant
  const double k_wave = std::sqrt(en * (2.0 + alpha * alpha * en));
  const double approx_wavelength = 2.0 * M_PI / k_wave;

  // The solution is only usable where the radial grid resolves the
  // oscillations. The grid spacing dr grows with r, so (for high energy)
  // this bounds the usable radius: solve and store the solution only up to
  // there (max_pt), and zero the remaining tail. For low energies, this is
  // the entire grid. Require:
  //  - at least N_ppw points per wavelength for the stored solution
  //  - a denser N_ppw_norm at the point where the normalisation
  //    integration starts (its accuracy is limited by the accumulated ODE
  //    error at that point)
  const int N_ppw = 10;
  const int N_ppw_norm = 40;
  const auto resolved_until = [&](double dr_max) {
    auto ip = num_points;
    while (ip > 0 && gr.drdu(ip - 1) * gr.du() > dr_max) {
      --ip;
    }
    return ip;
  };
  const auto i_max = resolved_until(approx_wavelength / N_ppw);
  const auto i_norm = resolved_until(approx_wavelength / N_ppw_norm);

  // If even the low-r region is unresolvable, write zeros and return:
  if (i_norm < 2 * Param::K_Adams) {
    Fa *= 0.0;
    Fa.eps() = 1.0;
    return;
  }

  // Solve on regular grid, where it is dense enough (up to i_norm):
  DiracDerivative Hd(Fa.grid(), v, Fa.kappa(), Fa.en(), alpha, {}, VxFa, Fa0);
  solve_Dirac_outwards(Fa.f(), Fa.g(), Hd, i_norm);

  // Values at the start of the normalisation integration (moved to the end
  // of the marginal band below, if there is one):
  auto r_norm = gr.r(i_norm - 1);
  auto f_norm = Fa.f(i_norm - 1);
  auto g_norm = Fa.g(i_norm - 1);

  // Marginal band [i_norm, i_max): the grid resolves the oscillations well
  // enough to store/use the solution (>= N_ppw), but not well enough to
  // integrate the ODE accurately on it (< N_ppw_norm; errors grow quickly,
  // worse for high energy and high kappa). Instead, integrate through this
  // region on a fine fixed-step grid, with the (linearly interpolated)
  // potential, and sample the solution back onto the grid points using
  // cubic Hermite (values + ODE derivatives). The Hermite error,
  // ~(2*pi/N_band)^4/384, sets the pointwise accuracy: ~1e-10 for N=400.
  // nb: exchange (VxFa) and inhomogeneous (Fa0) terms are neglected in the
  // band; it only exists at high energy, where they are negligible.
  if (i_max > i_norm) {
    const auto dr_band = approx_wavelength / 400.0;
    // Start a few grid points early, so the Adams starting block (first
    // K_Adams fine points) lies below the first sampled grid point:
    const auto i_begin = i_norm - 1 - Param::K_Adams;
    InterpPotentialDerivative Hband(gr, v, Fa.kappa(), en, alpha);
    AdamsMoulton::ODESolver2D<Param::K_Adams, double, double> ode{dr_band,
                                                                  &Hband};
    ode.solve_initial_K(gr.r(i_begin), Fa.f(i_begin), Fa.g(i_begin));

    const auto df = [&](double t, double ff, double gg) {
      return Hband.a(t) * ff + Hband.b(t) * gg;
    };
    const auto dg = [&](double t, double ff, double gg) {
      return Hband.c(t) * ff + Hband.d(t) * gg;
    };

    auto t_prev = ode.last_t();
    auto f_prev = ode.last_f();
    auto g_prev = ode.last_g();
    std::size_t i_next = i_norm;
    while (i_next < i_max) {
      ode.drive();
      const auto t = ode.last_t();
      const auto ff = ode.last_f();
      const auto gg = ode.last_g();
      // fill any grid points inside (t_prev, t], cubic Hermite:
      while (i_next < i_max && gr.r(i_next) <= t) {
        const auto h = t - t_prev;
        const auto x = (gr.r(i_next) - t_prev) / h;
        const auto x2 = x * x, x3 = x2 * x;
        const auto h00 = 2.0 * x3 - 3.0 * x2 + 1.0;
        const auto h10 = x3 - 2.0 * x2 + x;
        const auto h01 = -2.0 * x3 + 3.0 * x2;
        const auto h11 = x3 - x2;
        Fa.f()[i_next] = h00 * f_prev + h10 * h * df(t_prev, f_prev, g_prev) +
                         h01 * ff + h11 * h * df(t, ff, gg);
        Fa.g()[i_next] = h00 * g_prev + h10 * h * dg(t_prev, f_prev, g_prev) +
                         h01 * gg + h11 * h * dg(t, ff, gg);
        ++i_next;
      }
      t_prev = t;
      f_prev = ff;
      g_prev = gg;
    }
    // normalisation starts from the band end:
    r_norm = t_prev;
    f_norm = f_prev;
    g_norm = g_prev;
  }

  Fa.max_pt() = i_max;
  std::fill(Fa.f().begin() + long(i_max), Fa.f().end(), 0.0);
  std::fill(Fa.g().begin() + long(i_max), Fa.g().end(), 0.0);

  // Now, normalise the solution.
  // Keep solving ODE outwards from r_norm, on linearly-spaced grid
  // Use "H-like" derivative (local effective charge): assumes the potential
  // is Coulomb-like, and exchange etc. negligable, beyond r_norm.
  // Average the amplitude estimate over full oscillation cycles, extrapolate
  // the cycle means to r -> infinity, and re-scale the wavefunction so the
  // large-r amplitude matches the analytic (energy-normalised) expression

  // Step-size for large-r solution: uses linear grid, with fixed number of
  // points per wavelength (cycle means are sensitive to node location, which
  // is resolved to ~dr^2 by interpolation)
  const auto dr = approx_wavelength / 800.0;

  const auto ztmp = -1.0 * v.at(i_max - 1) * gr.r(i_max - 1);
  const auto Zeff = std::max(1.0, ztmp);

  // Find large-r amplitude:
  const auto [amp, eps_amp] = numerical_f_amplitude(en, Fa.kappa(), alpha, Zeff,
                                                    f_norm, g_norm, r_norm, dr);

  Fa.eps() = eps_amp;

  // Calculate normalisation coeficient, D, and re-scaling factor:
  const auto D = analytic_f_amplitude(en, alpha);
  Fa *= (amp != 0.0 ? (D / amp) : 0.0);

  // Optionally, fill the unresolved (zeroed) tail with local averages:
  if (average_tail) {
    averageTail(Fa, v, alpha);
  }
}

//==============================================================================
GridRequirements RequiredContinuumGrid(double en, const Grid &gr, double N_ppw,
                                       double alpha) {

  // Storing the wavefunction pointwise requires at least N_ppw points per
  // wavelength at the largest grid spacing (the last grid point).
  // The constraint is always at large r, where the grid is coarsest: at
  // small r the spacing shrinks linearly with r, faster than the local
  // wavelength does (k(r) grows only as sqrt(Z/r)).
  // relativistic wavenumber:
  const double k = std::sqrt(en * (2.0 + alpha * alpha * en));
  const double lambda = 2.0 * M_PI / k;
  const double dr_target = lambda / N_ppw;

  const auto n = gr.num_points();
  const double r0 = gr.r(0);
  const double rmax = gr.r(n - 1);

  // Number of points required, geometry (r0, rmax, type, b) unchanged.
  // Convert the target r-spacing at rmax to the uniform u-step:
  // dr = drdu * du, and drdu(rmax) is fixed by the geometry:
  const double du_req = dr_target / gr.drdu(n - 1);
  const auto npts_req =
    Grid::calc_num_points_from_du(r0, rmax, du_req, gr.type(), gr.loglin_b());

  // For a loglinear grid, u = r + b*ln(r):
  //   dr(rmax) = du * rmax/(rmax + b),
  //   du = [(rmax - r0) + b*ln(rmax/r0)] / (n - 1)
  // Solving dr(rmax) = dr_target for b (with n fixed) is exact; dr(rmax)
  // grows with b, so this is the largest b that suffices. If below b_min
  // (including negative: no b > 0 is small enough), clamp to b_min:
  const double b_min = 0.01;
  const double L = std::log(rmax / r0);
  const double tn = dr_target * double(n - 1);
  const double b_exact = rmax * (tn - (rmax - r0)) / (L * rmax - tn);
  const double b_req = std::max(b_exact, b_min);

  // Number of points required at the (possibly clamped) suggested b.
  // If not clamped, current num_points is (exactly) sufficient at b_req:
  const auto npts_b =
    b_exact >= b_min ?
      n :
      Grid::calc_num_points_from_du(r0, rmax, dr_target * (rmax + b_req) / rmax,
                                    GridType::loglinear, b_req);

  return {npts_req, b_req, npts_b};
}

//==============================================================================
std::size_t averageTail(DiracSpinor &Fa, const std::vector<double> &v,
                        double alpha) {

  const auto &gr = Fa.grid();
  const auto num_points = gr.num_points();
  const auto en = Fa.en();

  // Points-per-wavelength below which the pointwise solution cannot be
  // stored on the grid, and local averaging takes over. Must sit just
  // above the ~10 ppw zeroing threshold of solveContinuum: where the grid
  // still stores the solution pointwise (>= 12 ppw), averaging would only
  // suppress a perfectly usable oscillation (biasing integrals), so it
  // must not engage there:
  const int N_ppw_avg = 12;
  // Kernel width: sigma(r) = eta * (dr(r) - wavelength / N_ppw_avg), which
  // grows smoothly from zero at the seam. Larger eta suppresses the
  // unstorable oscillation more strongly (wider average):
  const double eta = 10.0;

  // relativistic wavelength (as in solveContinuum):
  const double approx_wavelength =
    2.0 * M_PI / std::sqrt(en * (2.0 + alpha * alpha * en));

  // Seam: first grid point where spacing exceeds wavelength / N_ppw_avg:
  std::size_t i_avg = num_points;
  for (std::size_t i = 0; i < num_points; ++i) {
    if (gr.drdu(i) * gr.du() > approx_wavelength / N_ppw_avg) {
      i_avg = i;
      break;
    }
  }

  // Nothing to do (grid resolves the entire solution), or no valid
  // pointwise solution to continue from:
  if (i_avg == num_points || i_avg < Param::K_Adams + 2 ||
      i_avg > Fa.max_pt()) {
    return num_points;
  }
  const std::size_t i_start = i_avg - 1;
  if (Fa.f(i_start) == 0.0 && Fa.g(i_start) == 0.0) {
    return num_points;
  }

  // Fine fixed-step grid; same step as the band solve in solveContinuum:
  const double h = approx_wavelength / 400.0;

  const auto sigma_at = [&](std::size_t i) {
    const auto s = eta * (gr.drdu(i) * gr.du() - approx_wavelength / N_ppw_avg);
    // at least one fine step wide (harmless near the seam, avoids an
    // empty kernel window):
    return std::max(s, h);
  };

  // Solve the Dirac ODE on the fine grid, from the last pointwise grid
  // point through to the end of the radial grid. Kernels of the last few
  // points reach slightly beyond r_max; the interpolated potential clamps
  // to v.back() there. Exchange is neglected (tail only exists at high
  // energy, where it is negligible).
  InterpPotentialDerivative Hband(gr, v, Fa.kappa(), en, alpha);
  AdamsMoulton::ODESolver2D<Param::K_Adams, double, double> ode{h, &Hband};
  ode.solve_initial_K(gr.r(i_start), Fa.f(i_start), Fa.g(i_start));

  const double r_end = gr.r(num_points - 1) + 5.0 * sigma_at(num_points - 1);
  std::vector<double> ts, fs, gs;
  const auto n_est = std::size_t((r_end - gr.r(i_start)) / h) + 16;
  ts.reserve(n_est);
  fs.reserve(n_est);
  gs.reserve(n_est);
  while (ode.last_t() < r_end) {
    ode.drive();
    ts.push_back(ode.last_t());
    fs.push_back(ode.last_f());
    gs.push_back(ode.last_g());
  }

  // Store the Gaussian local average at each remaining grid point,
  // kernel truncated at +/- 5 sigma:
  const double t0 = ts.front();
  for (std::size_t i = i_avg; i < num_points; ++i) {
    const double ri = gr.r(i);
    const double sigma = sigma_at(i);
    const auto j0 = std::size_t(std::max(0.0, (ri - 5.0 * sigma - t0) / h));
    const auto j1 =
      std::min(ts.size(), std::size_t((ri + 5.0 * sigma - t0) / h) + 1);
    double sw = 0.0, swf = 0.0, swg = 0.0;
    for (std::size_t j = j0; j < j1; ++j) {
      const double x = (ts[j] - ri) / sigma;
      const double w = std::exp(-0.5 * x * x);
      sw += w;
      swf += w * fs[j];
      swg += w * gs[j];
    }
    if (sw > 0.0) {
      Fa.f()[i] = swf / sw;
      Fa.g()[i] = swg / sw;
    }
  }

  Fa.max_pt() = num_points;
  return i_avg;
}

//==============================================================================
void solveContinuumIrregular(DiracSpinor &Firr, const DiracSpinor &Freg,
                             double en, const std::vector<double> &v,
                             double alpha) {
  // Build the irregular standing-wave partner of F_reg by seeding the inward
  // integration at the outermost grid points (see header).
  // Method B (Coulomb-series projection): project F_reg onto the asymptotic
  // Dirac-Coulomb standing pair {F^C, G^C} (1/r series), F_reg = a F^C + b
  // G^C, and seed F_irr = b F^C - a G^C -- the exact 90-degree partner to
  // series order, w[F_reg, F_irr] = (a^2+b^2) alpha/pi. This removes the
  // O(1/(p r_box)) phase-dependent F_reg admixture of the leading-order
  // component swap (method A), which otherwise contaminates the standing-wave
  // solution with an omega-oscillating on-shell admixture. Falls back to the
  // component swap if the projection is poor (series not converged: very low
  // p r_box, or non-Coulomb tail).

  Firr.en() = en;
  const auto &gr = Freg.grid();
  const auto kappa = Freg.kappa();
  const auto pinf = Freg.max_pt();
  Firr.max_pt() = pinf;

  // small/large amplitude ratio beta = sqrt(en/(en+2c^2))
  const auto c2 = 1.0 / (alpha * alpha);
  const double beta = std::sqrt(en / (en + 2.0 * c2));

  // Residual-ion charge seen at the box edge (Coulomb tail of v)
  const auto Zion = std::max(0.0, -v[pinf - 1] * gr.r(pinf - 1));
  const AsymptoticSpinorContinuum<15> asy{kappa, Zion, en, alpha};

  // Project F_reg = a F^C + b G^C pointwise over an outer window (2x2 in the
  // two components; W[F^C,G^C] = fC gG - fG gC = -alpha/pi), and average.
  const std::size_t n_w = std::min<std::size_t>(24, pinf / 2);
  double a_sum = 0.0, b_sum = 0.0;
  for (std::size_t i = pinf - n_w; i < pinf; ++i) {
    const auto [fC, gC, fG, gG] = asy.fg(gr.r(i));
    const auto Wcg = fC * gG - fG * gC;
    a_sum += (Freg.f(i) * gG - Freg.g(i) * fG) / Wcg;
    b_sum += (fC * Freg.g(i) - gC * Freg.f(i)) / Wcg;
  }
  const auto a = a_sum / double(n_w);
  const auto b = b_sum / double(n_w);
  // Both F_reg and {F^C,G^C} are energy-normalised, so a^2+b^2 = 1 when the
  // series represents F_reg well; use as the quality check for the fallback
  const auto q = a * a + b * b;
  const bool good_projection = std::abs(q - 1.0) < 0.1;

  // Homogeneous radial Dirac operator at en (same local potential as F_reg)
  DiracDerivative Hd(gr, v, kappa, en, alpha);

  const auto du = gr.du();
  AdamsMoulton::ODESolver2D<Param::K_Adams, std::size_t, double> ode{-du, &Hd};
  ode.S_scale = 0.0;

  auto &f = Firr.f();
  auto &g = Firr.g();

  // Seed the outermost K_steps points: method B projection, or method A
  // (component swap on F_reg) as the fallback
  for (std::size_t i0 = pinf - 1, i = 0; i < ode.K_steps(); ++i) {
    double f0{}, g0{};
    if (good_projection) {
      const auto [fC, gC, fG, gG] = asy.fg(gr.r(i0));
      f0 = b * fC - a * fG;
      g0 = b * gC - a * gG;
    } else {
      f0 = -Freg.g(i0) / beta;
      g0 = beta * Freg.f(i0);
    }
    ode.f[i] = f0;
    ode.g[i] = g0;
    ode.df[i] = ode.dfdt(f0, g0, i0);
    ode.dg[i] = ode.dgdt(f0, g0, i0);
    ode.t[i] = i0;
    --i0;
  }
  for (std::size_t i = 0; i < ode.K_steps(); ++i) {
    f.at(ode.t.at(i)) = ode.f.at(i);
    g.at(ode.t.at(i)) = ode.g.at(i);
  }

  // Integrate inwards to the origin
  const auto i_start = ode.last_t();
  for (std::size_t i = i_start - 1;; --i) {
    ode.drive(i);
    f.at(i) = ode.last_f();
    g.at(i) = ode.last_g();
    if (i == 0)
      break;
  }

  // Zero anything beyond the seeded region
  for (std::size_t i = pinf; i < f.size(); ++i) {
    f.at(i) = 0.0;
    g.at(i) = 0.0;
  }
}

//==============================================================================
double solveContinuumForward(DiracSpinor &phi, const DiracSpinor &Freg,
                             const DiracSpinor &Firr, double en,
                             const std::vector<double> &v, double alpha,
                             const DiracSpinor &Sr) {
  // Forward integration plus F_reg subtraction (see header). Outward
  // integration suppresses the irregular admixture (r^{-l} vs r^{l+1} near
  // the origin); the F_reg admixture from the zero seed is removed by the
  // c-subtraction, so the start is not critical.
  phi.en() = en;
  const auto &gr = Freg.grid();
  const auto kappa = Freg.kappa();
  const auto pinf = Freg.max_pt();
  phi.max_pt() = pinf;

  // (h_r - en)phi = Sr; DiracDerivative solves (h_r - en)phi = -VxFa.
  // solve_Dirac_outwards seeds the usual H-like regular form at r0 and
  // carries the source; the (arbitrary-scale) F_reg admixture this start
  // introduces is exactly what the c-subtraction below removes.
  const auto mSr = -1.0 * Sr;
  DiracDerivative Hd(gr, v, kappa, en, alpha, {}, &mSr);
  solve_Dirac_outwards(phi.f(), phi.g(), Hd, pinf);

  auto &f = phi.f();
  auto &g = phi.g();
  for (std::size_t i = pinf; i < f.size(); ++i) {
    f.at(i) = 0.0;
    g.at(i) = 0.0;
  }

  // Beyond the source: phi~ = c*F_reg + K*F_irr. Extract (c, K) pointwise
  // from the 2x2 component system over the outer grid, and average.
  const auto iw0 = std::size_t(0.70 * double(pinf));
  const auto iw1 = std::size_t(0.95 * double(pinf));
  double c_sum = 0.0, K_sum = 0.0;
  int n_w = 0;
  for (std::size_t i = iw0; i < iw1; ++i) {
    const auto w_i = Freg.f(i) * Firr.g(i) - Firr.f(i) * Freg.g(i);
    if (w_i == 0.0)
      continue;
    c_sum += (f[i] * Firr.g(i) - g[i] * Firr.f(i)) / w_i;
    K_sum += (Freg.f(i) * g[i] - Freg.g(i) * f[i]) / w_i;
    ++n_w;
  }
  const auto c = (n_w > 0) ? c_sum / n_w : 0.0;
  const auto K = (n_w > 0) ? K_sum / n_w : 0.0;

  // Impose the oscillating boundary condition: phi -> K*F_irr
  phi -= c * Freg;

  return K;
}

//==============================================================================
std::pair<double, double> numerical_f_amplitude(double en, int kappa,
                                                double alpha, double Zeff,
                                                double f_final, double g_final,
                                                double r_final, double dr) {

  // In the asymptotic region (V ~ -Zeff/r) the Dirac ODE reduces to:
  //   f'' = -k^2 f,  k^2 = en * (alpha^2 * en + 2)
  // giving f = A sin(theta), g = A*rho*cos(theta),
  // where rho = sqrt(alpha^2*en / (alpha^2*en + 2)), and
  // theta = kr + nu*ln(2kr) + const is the Coulomb phase.
  // Combining: A^2 = f^2 + c_g^2 * g^2,  c_g^2 = 1 + 2/(alpha^2 * en).
  // This is exact only asymptotically: at finite r the estimator oscillates
  // within each cycle (relative amplitude ~ nu/kr), and its mean approaches
  // the true amplitude algebraically. So:
  //  - average over full cycles, delimited by zero-crossings of f (this
  //    stays phase-exact despite the log-phase, unlike fixed-length windows);
  //    node positions are located to sub-step accuracy by linear
  //    interpolation (cycle means are otherwise noise-limited by the
  //    O(dr) boundary jitter)
  //  - the cycle means converge as A(r) = A_inf + c2/r^2 + c3/r^3 + ...
  //    (cycle-averaging kills the oscillating 1/r term); fit four well
  //    separated cycle means to this form (through c4/r^4) to obtain A_inf
  // Converged when the extrapolated estimates using the full and half
  // integration ranges agree (consecutive estimates share abscissae, so
  // their agreement badly underestimates the remaining common bias).
  const double c_g2 = 1.0 + 2.0 / (alpha * alpha * en);

  DiracContinuumDerivative Heff(Zeff, kappa, en, alpha);
  AdamsMoulton::ODESolver2D<Param::K_Adams, double, double> ode{dr, &Heff};
  ode.solve_initial_K(r_final, f_final, g_final);

  const double eps_target = 1.0e-7;
  const std::size_t min_cycles = 8;
  const std::size_t max_cycles = 256;

  // Hard cap on total steps (guards against node detection failing):
  const double approx_wavelength = 2.0 * M_PI / std::sqrt(2.0 * en);
  const auto steps_per_cycle = std::size_t(approx_wavelength / dr) + 1;
  const auto max_steps = 2 * max_cycles * steps_per_cycle;

  // Cycle means A_j, at cycle centres r_j:
  std::vector<double> Aj, rj;
  Aj.reserve(max_cycles);
  rj.reserve(max_cycles);

  // Extrapolate r -> infinity: fit A(r) = A_inf + c2/r^2 + c3/r^3 + c4/r^4
  // through cycles (2j/5, 3j/5, 4j/5, j). Gaussian elimination with partial
  // pivoting; only A_inf (x[0]) is needed:
  const auto extrapolate = [&](std::size_t j) {
    const std::size_t idx[4] = {(2 * j) / 5, (3 * j) / 5, (4 * j) / 5, j};
    double M[4][5];
    for (int a = 0; a < 4; ++a) {
      const auto r = rj[idx[a]];
      M[a][0] = 1.0;
      M[a][1] = 1.0 / (r * r);
      M[a][2] = M[a][1] / r;
      M[a][3] = M[a][2] / r;
      M[a][4] = Aj[idx[a]];
    }
    for (int c = 0; c < 3; ++c) {
      int pv = c;
      for (int a = c + 1; a < 4; ++a) {
        if (std::abs(M[a][c]) > std::abs(M[pv][c])) {
          pv = a;
        }
      }
      if (pv != c) {
        std::swap_ranges(std::begin(M[c]), std::end(M[c]), std::begin(M[pv]));
      }
      for (int a = c + 1; a < 4; ++a) {
        const auto fac = M[a][c] / M[c][c];
        for (int cc = c; cc < 5; ++cc) {
          M[a][cc] -= fac * M[c][cc];
        }
      }
    }
    double x[4];
    for (int a = 3; a >= 0; --a) {
      double s = M[a][4];
      for (int cc = a + 1; cc < 4; ++cc) {
        s -= M[a][cc] * x[cc];
      }
      x[a] = s / M[a][a];
    }
    return x[0];
  };

  // Extrapolated estimate after each completed cycle (0 until min_cycles):
  std::vector<double> Ej;
  Ej.reserve(max_cycles);

  double amp = 0.0;
  double eps_amp = 1.0;

  // First zero-crossing of f starts the first cycle (discard partial cycle);
  // every second crossing after that completes a cycle. Trapezoid-integrate
  // A over each cycle, locating the cycle boundaries (nodes of f) to
  // sub-step accuracy by linear interpolation:
  bool started = false;
  int crossings = 0;
  double I_cycle = 0.0;
  double len = 0.0;
  double r_start = 0.0;
  double f_prev = ode.last_f();
  double A_prev = 0.0;

  for (std::size_t step = 0; step < max_steps; ++step) {
    ode.drive();
    const double f = ode.last_f();
    const double g = ode.last_g();
    const double A = std::sqrt(f * f + c_g2 * g * g);
    if (started) {
      I_cycle += 0.5 * (A_prev + A) * dr;
      len += dr;
    }
    if (f_prev * f < 0.0) {
      // node at r_node, fraction t through this step:
      const double t = f_prev / (f_prev - f);
      const double A_node = A_prev + t * (A - A_prev);
      const double r_node = ode.last_t() - (1.0 - t) * dr;
      // contribution of this step's segment beyond the node:
      const double I_post = 0.5 * (A_node + A) * (1.0 - t) * dr;
      const double len_post = (1.0 - t) * dr;
      if (!started) {
        started = true;
        I_cycle = I_post;
        len = len_post;
        r_start = r_node;
      } else if (++crossings == 2) {
        Aj.push_back((I_cycle - I_post) / (len - len_post));
        rj.push_back(r_start + 0.5 * (len - len_post));
        I_cycle = I_post;
        len = len_post;
        r_start = r_node;
        crossings = 0;

        const auto j = Aj.size() - 1;
        Ej.push_back(j + 1 >= min_cycles ? extrapolate(j) : 0.0);
        if (Ej[j] != 0.0) {
          amp = Ej[j];
          // compare against the estimate from half the range:
          if (Ej[j / 2] != 0.0) {
            eps_amp =
              std::abs(Ej[j] - Ej[j / 2]) / (0.5 * std::abs(Ej[j] + Ej[j / 2]));
            if (eps_amp < eps_target)
              break;
          }
        }
        if (Aj.size() == max_cycles)
          break;
      }
    }
    f_prev = f;
    A_prev = A;
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
