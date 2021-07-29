#pragma once
#include <vector>
class DiracSpinor;
class Grid;

/*!
@brief
Solve Dirac equation. Functions defined in "<DiracODE/DiracODE.hpp>"
@details
  - v is local potential (e.g., v = v_dir + v_nuc)
  - H_mag is off-diagonal magnetic form factor (QED radiative potential); enter
as '{}' to not use.
  - alpha: \f$\alpha = \lambda\alpha_0\f$ is the effective value of
fine-structure constant
*/

namespace DiracODE {

//! @brief Solves bound-state problem for local potential (en < 0)
/*! @details
\f[ (H_0 + v - \epsilon_a)F_a = 0\f]
en0 is initial energy guess (must be reasonably good).
log_eps: log10(eps); eps is convergence target for energy.
*/
void boundState(DiracSpinor &Fa, const double en0, const std::vector<double> &v,
                const std::vector<double> &H_mag, const double alpha,
                int log_eps = 14, const DiracSpinor *const VxFa = nullptr,
                const DiracSpinor *const Fa0 = nullptr, double zion = 1);

//! @brief For given energy en, solves (local) DE with correct boundary
//! conditions at the origin
void regularAtOrigin(DiracSpinor &Fa, const double en,
                     const std::vector<double> &v,
                     const std::vector<double> &H_mag, const double alpha);

//! @brief For given energy en, solves (local) DE with correct boundary
//! conditions at infinity
void regularAtInfinity(DiracSpinor &Fa, const double en,
                       const std::vector<double> &v,
                       const std::vector<double> &H_mag, const double alpha);

//! @brief For given energy en (en > 0), solves (local) DE for continuum state
//! (with energy normalisation).
/*! @details
ext_grid is an 'extended grid' (see Grid class); needed since we need to solve
to large r to enforce boundary conditions (especially for large en). r_asym0 is
initial guess for 'asymptotic region'. ext_grid must extend past r_asym0.
*/
void solveContinuum(DiracSpinor &Fa, const double en,
                    const std::vector<double> &v, const Grid &ext_grid,
                    const double r_asym0, const double alpha,
                    const DiracSpinor *const VxFa = nullptr,
                    const DiracSpinor *const Fa0 = nullptr);

//******************************************************************************

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
                          const DiracSpinor &source);

//! @brief Solves inhomogeneous Dirac equation
/*! @details
As above. Overload to accept/overwrite solution to Fa. kappa is taken from Fa.
*/
void solve_inhomog(DiracSpinor &Fa, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source);

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
                   const DiracSpinor &source);

} // namespace DiracODE
