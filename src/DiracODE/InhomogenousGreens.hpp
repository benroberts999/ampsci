#pragma once
#include <vector>
class DiracSpinor;

namespace DiracODE {

//==============================================================================

/*!
  @brief Solves the inhomogeneous Dirac equation, returning the solution spinor.
  @details
  Solves \f$ (H_0 + v - \epsilon_a) F_a = S \f$ for \f$ \psi_\kappa \f$ using
  Green's method (see Methods documentation). Note the sign convention for S.
  @param kappa   Angular momentum kappa of the solution.
  @param en      Orbital energy \f$ \epsilon \f$.
  @param v       Local potential v(r).
  @param H_mag   Off-diagonal (magnetic) potential.
  @param alpha   Fine-structure constant.
  @param source  Inhomogeneous source term S.
  @param VxFa    Optional exchange potential. If nullptr, ignored.
  @param Fa0     Optional inhomogeneous source spinor. If nullptr, ignored.
  @param zion    Effective ionic charge (default 1).
  @param mass    Effective particle mass in atomic units (default 1 = m_e).
  @return Solution spinor Fa.
*/
DiracSpinor solve_inhomog(const int kappa, const double en,
                          const std::vector<double> &v,
                          const std::vector<double> &H_mag, const double alpha,
                          const DiracSpinor &source,
                          const DiracSpinor *const VxFa = nullptr,
                          const DiracSpinor *const Fa0 = nullptr,
                          double zion = 1, double mass = 1.0);

/*!
  @brief Solves the inhomogeneous Dirac equation, overwriting Fa.
  @details
  As above; kappa is taken from Fa.
*/
void solve_inhomog(DiracSpinor &Fa, const double en,
                   const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source,
                   const DiracSpinor *const VxFa = nullptr,
                   const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                   double mass = 1.0);

/*!
  @brief Solves the inhomogeneous Dirac equation, overwriting Fa and exposing the homogeneous solutions.
  @details
  As above, but also returns the homogeneous solutions Fzero (regular at origin)
  and Finf (regular at infinity), which satisfy \f$ (H_0 + v - \epsilon)F = 0 \f$.
  The first two overloads discard these; this one keeps them for potential reuse.
  Fzero and Finf are out parameters -- they are overwritten internally and do not
  need to be initialised before calling.
*/
void solve_inhomog(DiracSpinor &Fa, DiracSpinor &Fzero, DiracSpinor &Finf,
                   const double en, const std::vector<double> &v,
                   const std::vector<double> &H_mag, const double alpha,
                   const DiracSpinor &source,
                   const DiracSpinor *const VxFa = nullptr,
                   const DiracSpinor *const Fa0 = nullptr, double zion = 1,
                   double mass = 1.0);

//==============================================================================

namespace Internal {

//! Constructs the particular solution Fa from homogeneous solutions Finf, Fzero and source Sr.
void GreenSolution(DiracSpinor &Fa, const DiracSpinor &Finf,
                   const DiracSpinor &Fzero, const double alpha,
                   const DiracSpinor &Sr);

} // namespace Internal
} // namespace DiracODE
