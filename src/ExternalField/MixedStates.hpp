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

  @note Near-resonant channels are handled automatically: \f$ (h_{\rm HF} -
  \en_a \mp \omega) \f$ is (near-)singular for components along any same-kappa
  bound state with \f$ \en_m \approx \en_a \pm \omega \f$ (the diagonal @p Fa,
  fine-structure partners, etc.). These are projected out of the source,
  the solution is forced orthogonal to them, and the off-diagonal components are
  restored analytically. The caller need not pre-condition @p Fs.
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

/*!
  @brief Find bound states of the solve channel that make (h_l - e0) near-singular.
  @details
  For the channel of kappa @p kappa and energy \f$ e_0 = \en_a \pm \omega \f$,
  the radial operator \f$ (h_l - e_0) \f$ is (near-)singular for components
  along any bound core state \f$ \phi_m \f$ of the same kappa with
  \f$ \en_m \approx e_0 \f$: exactly singular for the diagonal (\f$ \phi_a \f$,
  \f$ \omega = 0 \f$) case, near-singular for e.g. fine-structure partners.
  Those components cannot be resolved reliably by the Green's-function solve, so
  they are projected out of the source (forces the solution orthogonal
  to them), they should be restore the off-diagonal ones analytically afterwards.

  The set is the same-kappa core states satisfying a relative nearness criterion
  \f$ |e_0 - \en_m| < \eta\,|e_0 + \en_m| \f$ (\f$ \eta = 0.2 \f$), plus @p Fa
  itself when it shares the channel kappa (the \f$ \matel{a}{\delta F}{} = 0 \f$
  / left-orthogonality constraint).
*/
std::vector<const DiracSpinor *>
conditioning_states(const std::vector<DiracSpinor> &core, const DiracSpinor &Fa,
                    int kappa, double e0);

} // namespace ExternalField
