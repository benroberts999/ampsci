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

/*!
  @brief Solves the inhomogeneous TDHF (mixed-states) equation for perturbed orbital dF.
  @details
  Solves
  \f[
    (h_{\rm HF} - \en_a \mp \omega)\delta F + F_S = 0
  \f]
  for \f$ \delta F \f$, where \f$ F_S \f$ is the source term. Typically
  \f[
    F_S = (t_\pm + \delta V_\pm - \delta\en^a_\pm)\phi_a.
  \f]
  
  - The angular momentum \f$ \kappa \f$ of the solution is that of @p Fs.
  - \f$ t \f$: Extenral field operator
  - \f$ \delta V_\pm \f$: core polarisation correction (see @ref CorePolarisation)
  - Solved iteratively using the Green's function method.

  @param Fa          Unperturbed orbital \f$ \phi_a \f$.
  @param omega       External-field frequency \f$ \omega \f$.
  @param vl          Local potential (nuclear + direct).
  @param alpha       Fine-structure constant.
  @param core        Core electrons (for exchange).
  @param Fs          Source term \f$ F_S \f$ (note sign: this is \f$ h\phi_a \f$,
                     not \f$ -h\phi_a \f$).
  @param eps_target  Convergence goal for the inhomogeneous ODE solver.
  @param Sigma       Optional correlation potential.
  @param VBr         Optional Breit interaction.
  @param H_mag       Magnetic part of QED radiative potential (electric part
                     should be included in @p vl).

  @return Perturbed orbital \f$ \delta F \f$.
*/
DiracSpinor solveMixedState(
  const DiracSpinor &Fa, double omega, const std::vector<double> &vl,
  double alpha, const std::vector<DiracSpinor> &core, const DiracSpinor &Fs,
  double eps_target = 1.0e-9,
  const MBPT::CorrelationPotential *const Sigma = nullptr,
  const HF::Breit *const VBr = nullptr, const std::vector<double> &H_mag = {});

/*!
  @brief As solveMixedState(), but updates an existing solution @p dF in place.
  @details
  Starts from @p dF as an initial guess rather than zero; converges faster if
  @p dF is already an approximate solution (e.g., from a nearby frequency).
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

/*!
  @brief Solves for dF via explicit sum over basis; mainly for tests.
  @details
  \f[
    \delta F = \sum_n \frac{\ket{n}\matel{n}{F_S}{a}}{\en_a - \en_n \pm \omega}
  \f]
  where @p hFa is the already-evaluated source spinor \f$ F_S \f$.
*/
DiracSpinor solveMixedState_basis(const DiracSpinor &Fa, const DiracSpinor &hFa,
                                  double omega,
                                  const std::vector<DiracSpinor> &basis);

} // namespace ExternalField
