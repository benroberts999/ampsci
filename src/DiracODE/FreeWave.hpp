#pragma once
#include <cstddef>
#include <memory>
#include <vector>
class DiracSpinor;
class Grid;

namespace DiracODE {

/*!
  @brief Momentum, normalisation, and resolved extent of the free (V = 0)
  Dirac waves at one energy on a grid; see freeWaves().
*/
struct FreeWaveParameters {
  //! Momentum, k = sqrt(en (2 + alpha^2 en)), in au
  double k{};
  //! Amplitude N = k D of f = N r j_l(kr): D the energy-normalised large-r
  //! amplitude (analytic_f_amplitude)
  double norm{};
  //! Small-to-large ratio of a free particle, k alpha / (2 + alpha^2 en)
  double g_ratio{};
  //! Extent of the stored solution: the grid resolves the oscillations (at
  //! least ~10 points per wavelength, as solveContinuum) below this index
  std::size_t max_pt{};
};

//! The parameters of the free waves at energy en (> 0, au) on the grid
FreeWaveParameters freeWaveParameters(double en, const Grid &grid,
                                      double alpha);

/*!
  @brief Free (V = 0) Dirac spherical waves at energy en, for every kappa
  with l in [min_l, max_l]; energy normalised.

  @details
  Exact solutions of the free Dirac equation, regular at the origin,
  \f[
  \begin{align}
    f(r) &= N\, r\, j_l(kr), \\
    g(r) &= \frac{\kappa}{|\kappa|}\, N\, \frac{k\alpha}{2 + \alpha^2\en}\,
            r\, j_{\tilde l}(kr),
  \end{align}
  \f]
  with
  \f[
  \begin{align}
    k &= \sqrt{\en(2 + \alpha^2\en)}, \\
    \tilde l &= l(-\kappa), \\
    N &= k D,
  \end{align}
  \f]
  (\f$ \tilde l \f$ is l+1 for kappa < 0, l-1 for kappa > 0), where D is the
  energy-normalised large-r amplitude (analytic_f_amplitude()),
  \f[
    f \to D\sin(kr - l\pi/2),
  \f]
  so that
  \f[
    \int (f_\en f_{\en'} + g_\en g_{\en'})\,dr = \delta(\en - \en').
  \f]
  Same (f, g) convention as solveContinuum(); the non-relativistic limit of
  f is DiracContinuum::P_el with zero charge. Standing waves (not outgoing):
  only |amplitude|^2 summed over channels is meaningful.

  As solveContinuum(), the solution is stored only where the grid resolves
  the oscillations (at least ~10 points per wavelength): max_pt() is set
  accordingly and the tail zeroed, so that free and distorted waves at the
  same energy are truncated alike (as required when completing a truncated
  multipole sum with plane waves).

  All l are evaluated together, from one Bessel recurrence per grid point;
  for a single kappa use freeWave().

  @param en     Continuum (kinetic) energy, > 0, in au.
  @param min_l  Minimum orbital l.
  @param max_l  Maximum orbital l.
  @param grid   Radial grid.
  @param alpha  Fine-structure constant.
  @return The waves, in kappa-index order (-1, 1, -2, 2, ...), with n = 0.

  @note Do not obtain free waves from solveContinuum() with a zero
        potential: its normalisation assumes a Coulomb tail.
*/
std::vector<DiracSpinor> freeWaves(double en, int min_l, int max_l,
                                   std::shared_ptr<const Grid> grid,
                                   double alpha);

//! Free (V = 0) Dirac spherical wave of a single kappa; see freeWaves()
DiracSpinor freeWave(double en, int kappa, std::shared_ptr<const Grid> grid,
                     double alpha);

} // namespace DiracODE
