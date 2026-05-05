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
