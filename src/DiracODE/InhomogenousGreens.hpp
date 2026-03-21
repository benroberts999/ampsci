#pragma once
#include <vector>
class DiracSpinor;

namespace DiracODE {

//==============================================================================

//! @brief Solves inhomogeneous Dirac equation
/*! @details
\f[ (H_0 + v -\epsilon_a)F_a = S \f]
with `source' term, S. Solves for \f$\psi_\kappa\f$ with angular momentum kappa.
en = \f$\epsilon\f$ is given. Note sign of S.
Uses Green's method (see Method documentation).
*/
DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v,
                          const std::vector<double> &H_mag, const double alpha,
                          const DiracSpinor &source,
                          const DiracSpinor *const VxFa = nullptr,
                          const DiracSpinor *const Fa0 = nullptr,
                          double zion = 1, double mass = 1.0);

//! @brief Solves inhomogeneous Dirac equation
/*! @details
As above. Overload to accept/overwrite solution to Fa. kappa is taken from Fa.
*/
void solve_inhomog(DiracSpinor &Fa, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source,
                   const DiracSpinor *const VxFa = nullptr,
                   const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                   double mass = 1.0);

//! @brief Solves inhomogeneous Dirac equation
/*! @details
As above. Overload to accept/overwrite solution to Fa.
All these routines solve also for Fzero, Finf, which are solutions to
homogeneous equation  (H-en)Fa = 0 [reg @ origin, and infinity, respectively].
  - The first two throw these solutions away, the third keeps them (in some
cases they can be re-used)
  - These Spinors are solved internally and over-written, they don't need to be
solved first (i.e., they are out parameters, not in/out parameters)
*/
void solve_inhomog(DiracSpinor &Fa, DiracSpinor &Fzero, DiracSpinor &Finf,
                   const double en, const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source,
                   const DiracSpinor *const VxFa = nullptr,
                   const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                   double mass = 1.0);

//==============================================================================

//! @brief Solves inhomogeneous Dirac equation for continuum energies (en > 0)
/*! @details
As solve_inhomog, but for en > 0. Uses the real principal-value Green's
function with regularAtOrigin (\f$f_{\rm reg}\f$) and irregularAtOrigin
(\f$f_{\rm irr}\f$) as the two homogeneous solutions:
\f[
  \chi(r) = \frac{\alpha}{W}\left[
    f_{\rm irr}(r)\int_0^r f_{\rm reg}\cdot S\,dr'
    + f_{\rm reg}(r)\int_r^\infty f_{\rm irr}\cdot S\,dr'
  \right]
\f]
where \f$W\f$ is the Wronskian. Since \f$S\f$ is a localised core orbital,
both integrals converge and \f$\chi\f$ is real (Johnson's method).
Integration is limited to source.max_pt(), avoiding the oscillatory large-r region.
*/
DiracSpinor solve_inhomog_continuum(const int kappa, const double en,
                                    const std::vector<double> &v,
                                    const std::vector<double> &H_mag,
                                    const double alpha,
                                    const DiracSpinor &source,
                                    const DiracSpinor *const VxFa = nullptr,
                                    const DiracSpinor *const Fa0 = nullptr,
                                    double zion = 1, double mass = 1.0);

//! @brief As above. Overload: accepts/overwrites Fa; kappa taken from Fa.
void solve_inhomog_continuum(DiracSpinor &Fa, const double en,
                             const std::vector<double> &v,
                             const std::vector<double> &H_mag,
                             const double alpha, const DiracSpinor &source,
                             const DiracSpinor *const VxFa = nullptr,
                             const DiracSpinor *const Fa0 = nullptr,
                             double zion = 1, double mass = 1.0);

//==============================================================================

namespace Internal {

// Takes solution regular at infinity (Finf), and that regular at zero (Fzero),
// and the inhomogenous source term, Sr, to find particular solution, Fa.
// If w2_override != 0, uses it as the Wronskian instead of computing internally.
void GreenSolution(DiracSpinor &Fa, const DiracSpinor &Finf,
                   const DiracSpinor &Fzero, const double alpha,
                   const DiracSpinor &Sr, double w2_override = 0.0);

} // namespace Internal
} // namespace DiracODE
