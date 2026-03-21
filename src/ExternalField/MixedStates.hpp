#pragma once
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;
namespace MBPT {
class CorrelationPotential;
}
namespace HF {
class HartreeFock;
class Breit;
} // namespace HF

//! External field: Mixed-states + Core Polarisation
namespace ExternalField {

constexpr bool print_final_eps = false;
constexpr bool print_each_eps = false;

//! Solves Mixed States (TDHF) equation, inhomogenous equation, with
//! Hartree-Fock Hamiltonian, including exchange
/*! @details
Solves
\f[ 
  (H_{\rm HF} - \epsilon - \omega)\delta\phi + F_S = 0
\f]
for \f$\delta\phi\f$ (dF). Typically
\f[ 
  F_S = (\hat h  + \delta V_h - \delta \varepsilon) \phi.
\f]
Requires unperturbed orbital, Fa, a local potential (vl = vnuc + vdir), set of core electrons (for exchange). 
kappa (Dirac Q. number) of solution will have that of Fs.
Note sign on hFa 
  (this is \f$\hat h \phi\f$, not \f$-\hat h \phi\f$). 
eps_target is convergance goal for solving the inhomogenous dif.
equation.
Can optionally include Sigma (correlations), Breit, and Magnetic part of QED radiative potential (electric part should be included in vl).
*/
DiracSpinor solveMixedState(
  const DiracSpinor &Fa, double omega, const std::vector<double> &vl,
  double alpha, const std::vector<DiracSpinor> &core, const DiracSpinor &Fs,
  double eps_target = 1.0e-9,
  const MBPT::CorrelationPotential *const Sigma = nullptr,
  const HF::Breit *const VBr = nullptr, const std::vector<double> &H_mag = {});

//! Solves Mixed States (TDHF) equation, inhomogenous equation, staring from existing approximate solition, dF.
/*! @details
As above [solveMixedState], but starts with existing solution dF (may be 'zero'). If the existing solution is already approximate solution, this allows equation to be solved
much quicker.
*/
void solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa, double omega,
                     const std::vector<double> &vl, double alpha,
                     const std::vector<DiracSpinor> &core,
                     const DiracSpinor &Fs, double eps_target = 1.0e-9,
                     const MBPT::CorrelationPotential *const Sigma = nullptr,
                     const HF::Breit *const VBr = nullptr,
                     const std::vector<double> &H_mag = {});

//! Solves Mixed States (TDHF) equation. Overload; takes hf object
DiracSpinor
solveMixedState(const DiracSpinor &Fa, double omega, const DiracSpinor &Fs,
                const HF::HartreeFock *const hf, double eps_target = 1.0e-9,
                const MBPT::CorrelationPotential *const Sigma = nullptr);

//! Solves Mixed States (TDHF) equation. Overload; takes hf object
void solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa, double omega,
                     const DiracSpinor &Fs, const HF::HartreeFock *const hf,
                     double eps_target = 1.0e-9,
                     const MBPT::CorrelationPotential *const Sigma = nullptr);

//! @brief Solves the continuum mixed-states (TDHF) equation for en > 0
/*! @details
Solves the inhomogeneous Dirac equation for the RPA response of a core
orbital \f$\phi_a\f$ to an external perturbation at photoionisation
frequency \f$\omega\f$:
\f[
  (H_{\rm HF} - \epsilon_a - \omega)\,\chi_a = -F_S,
  \qquad \epsilon_a + \omega > 0
\f]
Uses the real principal-value Green's function (Johnson's method), so that
\f$\chi_a\f$ is real. Exchange is included iteratively as in solveMixedState.
No orthogonalisation is applied (\f$\chi_a\f$ is in the continuum and is
naturally orthogonal to bound core orbitals). Convergence is measured via
\f$\langle F_S|\chi_a\rangle\f$ to avoid the oscillatory large-r behaviour.
*/
DiracSpinor solveMixedState_continuum(
  const DiracSpinor &Fa, double omega, const std::vector<double> &vl,
  double alpha, const std::vector<DiracSpinor> &core, const DiracSpinor &Fs,
  double eps_target = 1.0e-9,
  const MBPT::CorrelationPotential *const Sigma = nullptr,
  const HF::Breit *const VBr = nullptr, const std::vector<double> &H_mag = {});

//! @brief As above, updating existing dF in-place (faster if dF is already approximate)
void solveMixedState_continuum(
  DiracSpinor &dF, const DiracSpinor &Fa, double omega,
  const std::vector<double> &vl, double alpha,
  const std::vector<DiracSpinor> &core, const DiracSpinor &Fs,
  double eps_target = 1.0e-9,
  const MBPT::CorrelationPotential *const Sigma = nullptr,
  const HF::Breit *const VBr = nullptr, const std::vector<double> &H_mag = {});

//! @brief Solves continuum mixed-states equation. Overload; takes hf object.
DiracSpinor
solveMixedState_continuum(const DiracSpinor &Fa, double omega,
                          const DiracSpinor &Fs,
                          const HF::HartreeFock *const hf,
                          double eps_target = 1.0e-9,
                          const MBPT::CorrelationPotential *const Sigma = nullptr);

//! @brief Solves continuum mixed-states equation. Overload; takes hf object.
void solveMixedState_continuum(DiracSpinor &dF, const DiracSpinor &Fa,
                               double omega, const DiracSpinor &Fs,
                               const HF::HartreeFock *const hf,
                               double eps_target = 1.0e-9,
                               const MBPT::CorrelationPotential *const Sigma = nullptr);

//! Directly defines dF via explicit sum over basis: mainly for tests.
//! \f$ \delta F = \sum_n |n\rangle\langle n|h|Fa\rangle / (e_a - e_n + \omega) \f$
DiracSpinor solveMixedState_basis(const DiracSpinor &Fa, const DiracSpinor &hFa,
                                  double omega,
                                  const std::vector<DiracSpinor> &basis);

} // namespace ExternalField
