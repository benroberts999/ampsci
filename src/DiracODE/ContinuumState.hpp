#pragma once
#include "AdamsMoulton.hpp"
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

/*!
  @brief Solves Dirac equation for a continuum state (en > 0) with energy normalisation.
  @details
  Normalisation is achieved by continuing the ODE integration to very large r and
  comparing the asymptotic amplitude to that of the analytic solution.
  Only the solution on the regular grid is kept; the extended part is discarded.
  @param Fa     Output spinor (result stored here).
  @param en     Continuum energy (must be > 0).
  @param v      Local potential v(r).
  @param alpha  Fine-structure constant.
  @param VxFa   Optional exchange potential. If nullptr, ignored.
  @param Fa0    Optional inhomogeneous source spinor. If nullptr, ignored.
*/
void solveContinuum(DiracSpinor &Fa, double en, const std::vector<double> &v,
                    double alpha, const DiracSpinor *const VxFa = nullptr,
                    const DiracSpinor *const Fa0 = nullptr);

/*!
  @brief Builds the irregular continuum partner F_irr of a regular continuum
  orbital F_reg, by projection onto the asymptotic Dirac-Coulomb pair
  (method B), with the component swap (method A) as fallback.
  @details
  Given the regular (outward-integrated, energy-normalised) continuum orbital
  F_reg at energy en>0 in the local potential v, constructs the linearly
  independent irregular oscilating partner F_irr that obeys the outer
  (90-degree phase-shifted) boundary condition. It is needed as the
  "regular-at-infinity" homogeneous solution F^inf in the continuum Green's
  function (replacing the bound decaying solution).

  Method B (Coulomb-series projection): over the outer grid points, F_reg is
  projected onto the energy-normalised asymptotic Dirac-Coulomb standing pair
  {F^C, G^C} (1/r series, @ref AsymptoticSpinorContinuum), using both spinor
  components pointwise:

  \f[ F_{\rm reg} = a\,F^C + b\,G^C
      \quad\Longrightarrow\quad F_{\rm irr} = b\,F^C - a\,G^C , 
  \f]

  the exact 90-degree-shifted partner to series order (no phase fit; the
  short-range phase is contained in a, b). These values seed an inward
  integration of the homogeneous radial Dirac equation
  \f$ (h_r - \en)F_{\rm irr} = 0 \f$ (same local potential v as F_reg).

  Method A (component swap) fallback: asymptotically
  F_reg = (A_L cos X, A_S sin X) with A_S = beta A_L,
  beta = sqrt(en/(en+2c^2)), so the 90-degree partner follows by eliminating
  the phase X: f_irr = -g_reg/beta, g_irr = beta f_reg. Used when the series
  projection is poor (a^2+b^2 far from 1: very small p*r_box, or a
  non-Coulomb tail). Method A is exact only to leading asymptotic order: the
  seed carries an O(1/(p*r_box)) phase-dependent admixture of F_reg, which
  contaminates the oscilating Green's solution with an omega-oscillating
  on-shell term (worst near threshold) -- the reason method B is preferred.

  The result has the conserved Wronskian
  \f$ w[F_{\rm reg},F_{\rm irr}] = f_{\rm reg}g_{\rm irr} - f_{\rm irr}g_{\rm reg}
  = \alpha/\pi \f$ (times a^2+b^2 = 1 for method B). It need not be
  normalised: the overall scale cancels in the Green's-function prefactor
  alpha/w.

  @param Firr   Output: the irregular partner (overwritten; kappa taken from Freg).
  @param Freg   Input: the regular, energy-normalised continuum orbital.
  @param en     Continuum energy (must be > 0; should equal Freg.en()).
  @param v      Local potential v(r) (the same used to solve Freg; its tail
                sets the residual-ion charge Z_ion for the Coulomb series).
  @param alpha  Fine-structure constant.

  @warning Requires en>0 and a grid dense enough at large r that F_reg is
           well-resolved there (same condition as solveContinuum).
*/
void solveContinuumIrregular(DiracSpinor &Firr, const DiracSpinor &Freg,
                             double en, const std::vector<double> &v,
                             double alpha);

/*!
  @brief Solves the inhomogeneous continuum equation (h_r - en)phi = S by
  forward integration plus F_reg subtraction; alternative to the Green's
  function method.
  @details
  Integrates the inhomogeneous radial Dirac equation

  \f[ (h_r - \en)\,\tilde\varphi = S, \qquad \en > 0, \f]

  outwards from the origin with regular (zero) initial conditions, giving a
  particular solution \f$ \tilde\varphi \f$ that is determined only up to an
  arbitrary admixture of the regular homogeneous solution F_reg (itself
  regular at the origin). The oscilating boundary condition is then imposed
  a posteriori by subtraction,

  \f[ \varphi(r) = \tilde\varphi(r) - c\,F_{\rm reg}(r), \f]

  with c chosen such that \f$ \varphi \propto F_{\rm irr} \f$ as r -> infinity.
  Beyond the (short-ranged) source, \f$ \tilde\varphi = c F_{\rm reg}
  + K F_{\rm irr} \f$ exactly; c and K are extracted pointwise over the outer
  grid by solving the 2x2 component system with the local Wronskian
  \f$ w_i = f_{\rm reg} g_{\rm irr} - f_{\rm irr} g_{\rm reg} \f$
  (no phase fit), and averaged.

  This is mathematically identical to the Green's-function construction
  [@ref Internal::GreenSolution with (F_reg, F_irr)]: the subtracted
  \f$ c F_{\rm reg} \f$ is exactly the admixture the Green's function never
  builds. Useful as an independent cross-check of the continuum mixed-state
  solve.

  @param phi    Output: the oscilating particular solution (kappa from Freg).
  @param Freg   Regular homogeneous solution at en (energy-normalised).
  @param Firr   Irregular homogeneous partner (see solveContinuumIrregular).
  @param en     Continuum energy (must be > 0).
  @param v      Local potential v(r) (same as used for Freg/Firr).
  @param alpha  Fine-structure constant.
  @param Sr     Source spinor S.

  @return The K-matrix amplitude: phi -> K * F_irr at large r.

  @note The c,K extraction window is the outer part of the grid; the source
        must have decayed there (it is built from bound orbitals, so this
        holds whenever the box is large enough for the Green's method too).
  @warning The outward integration is numerically delicate when the
           oscillations are dense (very near threshold); prefer the Green's
           route in production (see Methods, "TDHF for continuum states").
*/
double solveContinuumForward(DiracSpinor &phi, const DiracSpinor &Freg,
                             const DiracSpinor &Firr, double en,
                             const std::vector<double> &v, double alpha,
                             const DiracSpinor &Sr);

/*!
  @brief Analytic amplitude of f(r) at very large r for an H-like Dirac continuum state.
  @param en     Continuum energy.
  @param alpha  Fine-structure constant.
  @return Analytic asymptotic amplitude.
*/
double analytic_f_amplitude(double en, double alpha);

/*!
  @brief Finds the numerical amplitude and phase of f(r) for a continuum Dirac solution at large r.
  @details
  Continues ODE integration beyond the regular grid until both the wavelength
  and amplitude become constant. Assumes an H-like potential (-Zeff/r) and a
  linearly-spaced extension grid with step dr.
  @param en       Continuum energy.
  @param kappa    Orbital kappa quantum number.
  @param alpha    Fine-structure constant.
  @param Zeff     Effective nuclear charge.
  @param f_final  Value of f at the end of the regular grid.
  @param g_final  Value of g at the end of the regular grid.
  @param r_final  Radial position at the end of the regular grid.
  @param dr       Step size for the extended linear grid.
  @return {amplitude, phase} of the asymptotic f(r) oscillation.
*/
std::pair<double, double> numerical_f_amplitude(double en, int kappa,
                                                double alpha, double Zeff,
                                                double f_final, double g_final,
                                                double r_final, double dr);

/*!
  @brief H-like Dirac derivative matrix for continuum states at large r.
  @details
  Implements AdamsMoulton::DerivativeMatrix<double, double>, using r directly
  as the argument type. Valid for H-like potential (-Zeff/r); used to extend
  continuum integration beyond the regular grid for normalisation.
  @note Non-copyable.
*/
struct DiracContinuumDerivative
  : AdamsMoulton::DerivativeMatrix<double, double> {

  /*!
    @brief Constructs the H-like continuum derivative matrix.
    @param in_Zeff   Effective nuclear charge.
    @param in_kappa  Orbital kappa quantum number.
    @param in_en     Continuum energy.
    @param in_alpha  Fine-structure constant.
  */
  DiracContinuumDerivative(double in_Zeff, const int in_kappa,
                           const double in_en, const double in_alpha)
    : Zeff(in_Zeff),
      kappa(in_kappa),
      en(in_en),
      alpha(in_alpha),
      cc(1.0 / alpha) {}

  double Zeff = 1.0;
  int kappa;
  double en, alpha, cc;

  double a(double r) const final { return double(-kappa) / r; }
  double b(double r) const final {
    return (alpha * en + 2.0 * cc + Zeff * alpha / r);
  }
  double c(double r) const final { return -alpha * (Zeff / r + en); }
  double d(double r) const final { return -a(r); }
};

/*!
  @brief Fits a quadratic to three points and returns the interpolated maximum.
  @details
  Assumes |y2| = max(|y1|, |y2|, |y3|); used to find the amplitude of a
  sinusoidal oscillation. The three points must be close to the maximum.
  @param x1, x2, x3  x-coordinates of the three points.
  @param y1, y2, y3  y-coordinates of the three points.
  @return Interpolated maximum value.
*/
double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

} // namespace DiracODE
