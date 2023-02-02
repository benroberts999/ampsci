#pragma once
#include "AdamsMoulton.hpp"
#include <utility>
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

//! For given energy en (en > 0), solves Dirac eq. for continuum state (with
//! energy normalisation).
/*! @details
Normalisation is achieved by continuuing solving ODE to very large r, and 
comparing asymptotic amplitude to that of analytic solution.
We only keep solution on regular grid; extended part is not kept.
*/
void solveContinuum(DiracSpinor &Fa, double en, const std::vector<double> &v,
                     double alpha, const DiracSpinor *const VxFa = nullptr,
                     const DiracSpinor *const Fa0 = nullptr);

//! Analytic amplitude of f(r) at very large r, for H-like Dirac continuum
double analytic_f_amplitude(double en, double alpha);

//! Finds the (numerical) amplitude of f(r) continuum Dirac solution at large r
/*! @details
It does this by continuing ODE integration to large r until: \n
  (a) wavelength and, \n
  (b) amplitude \n
become constant.
It starts from a solution hat has been already integrated out to rmax of regular
radial grid. It uses linearly-spaced grid (dr), and assumes H-like potential
(-Z/r).
*/
std::pair<double, double> numerical_f_amplitude(double en, int kappa,
                                                double alpha, double Zeff,
                                                double f_final, double g_final,
                                                double r_final, double dr);

//! Derivative function for H-like; valid for continuum states at large r
struct DiracContinuumDerivative
    : AdamsMoulton::DerivativeMatrix<double, double> {

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

//! Fits a quadratic to three points {x,y}, assuming |y2| = max(|y1|,|y2|,|y3|).
//! @details Used to find amplitude of sin(x); points must be close to max
double fitQuadratic(double x1, double x2, double x3, double y1, double y2,
                    double y3);

} // namespace DiracODE
