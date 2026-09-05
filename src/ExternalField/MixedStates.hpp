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

/*!
  @brief Anderson mixing coefficients for a fixed-point
  iteration x -> G(x).
  @details
  Given the residuals \f$ r_k = G(x_k) - x_k \f$ of the stored iterates
  (oldest first) through their Gram matrix \f$ B_{kl} = \braket{r_k}{r_l} \f$,
  returns the \f$ c_k \f$ minimising
  \f[ \Big\|\sum_k c_k r_k\Big\|^2 \quad\text{subject to}\quad \sum_k c_k = 1, \f]
  via the bordered system [B 1; 1^T 0][c; lambda] = [0; 1]. The next iterate
  is \f$ x = \sum_k c_k G(x_k) \f$. For a linear map this is equivalent to
  (truncated) GMRES: it converges wherever the linear problem is
  non-singular, where damped iteration may not.

  A saturated history makes B ill-conditioned: the oldest entries are then
  dropped, and the returned vector holds one coefficient per KEPT entry, the
  last c.size() entries of the history (the caller drops the same entries
  from its own history). Empty if no entry is usable (take a plain step).
*/
std::vector<double>
anderson_coefficients(const std::vector<std::vector<double>> &gram);

//! As above, from the residual spinors directly (Gram matrix of their inner
//! products).
std::vector<double>
anderson_coefficients(const std::vector<DiracSpinor> &residuals);

/*!
  @brief Bound mixed-state solve with Anderson acceleration; the bound
  channels of @ref TDHFcntm.
  @details
  Same physics and conditioning as @ref solveMixedState (in-place overload),
  but the linear equation \f$ (h_{\rm HF} - \en_0)\,\delta F = -F_S \f$ is
  solved by Anderson mixing (@ref anderson_coefficients) of the preconditioned
  fixed-point map instead of damped iteration. At the high frequencies of
  ionisation a bound channel's \f$ \en_0 = \en_a + \omega \f$ can sit near a
  spurious eigenvalue of the local preconditioner
  \f$ (h_{\rm local} + U_x - \en_0) \f$, where damped iteration diverges
  although \f$ (h_{\rm HF} - \en_0) \f$ is non-singular; Anderson mixing
  converges on the conditioning of the true operator. Parameters as
  @ref solveMixedState.
*/
void solveMixedState_cntm(DiracSpinor &dF, const DiracSpinor &Fa, double omega,
                          const std::vector<double> &vl, double alpha,
                          const std::vector<DiracSpinor> &core,
                          const DiracSpinor &Fs, double eps_target = 1.0e-9,
                          const HF::Breit *const VBr = nullptr,
                          const std::vector<double> &H_mag = {});

/*!
  @brief Continuum (en_+ > 0) mixed-state solve with the standing-wave
  boundary condition and the non-local exchange iterated in the source: the
  open channels of @ref TDHFcntm.
  @details
  Solves
  \f[ 
    (h_r^{(\kappa)} + V^{\rm nl} - \en_+)\,\varphi = -F_S ,
      \qquad 
    \en_+ = \en_a + \omega > 0 , 
  \f]
  with \f$ \varphi \to K\,F_{\rm irr} \f$ at large r, by outward integration
  plus F_reg subtraction (@ref DiracODE::solveContinuumForward) on the two
  real homogeneous solutions at en_+: the regular, energy-normalised
  continuum orbital F_reg and its irregular partner F_irr
  (@ref DiracODE::solveContinuumIrregular), in the local potential
  \f$ v = v_l + U_x \f$, with \f$ U_x \f$ the orbital-independent Kohn-Sham
  exchange (@ref HF::vex_KS) for conditioning (cancelled in the source at the
  fixed point). The non-local exchange remainder is iterated in the source
  with Anderson mixing (@ref anderson_coefficients); K is linear in phi and
  is mixed with the same coefficients. After each solve, phi is
  orthogonalised to Fa when the channel shares its kappa (norm conservation,
  as the bound solver); other occupied components are kept, since their dV
  contributions cancel pairwise across channels.

  The hole-particle treatment is the caller's: pass @p vl already adjusted
  (\f$ v_l - y^0_{aa} \f$) and the hole orbital as @p Fhole, so that the
  iterated exchange is \f$ V^{\rm exch} - X_a \f$ (@ref HF::vexFa_1el);
  together these put the ejected electron in the V^{N-1} Hamiltonian (do
  both or neither).

  If @p Freg already holds the pair at en_+ (Freg->en() == en_+, nonzero),
  the homogeneous solutions are reused (they depend only on the channel and
  omega); otherwise they are built and written.

  @param phi    In/out: the correction orbital (kappa = channel kappa, as
                @p Fs). Used as the starting guess if nonzero.
  @param Freg   In/out: regular energy-normalised continuum orbital at en_+.
  @param Firr   In/out: irregular partner.
  @param K      Output (if non-null): standing-wave amplitude, phi -> K F_irr.
  @param Fa     Core orbital phi_a (energy and grid).
  @param omega  Frequency (en_+ = Fa.en() + omega must be > 0).
  @param vl     Local potential (including any hole-particle term).
  @param alpha  Fine-structure constant.
  @param core   Core orbitals (exchange).
  @param Fs     Source spinor F_S (sign as solveMixedState: +h phi_a).
  @param eps_target  Convergence goal of the exchange iteration.
  @param Fhole  Optional hole orbital (V^{N-1} exchange part).

  @warning Requires en_+ > 0 and a grid dense enough at large r; a reused
           pair must have been built in the same local potential (not
           checked).
*/
void solveContinuumMixedState(DiracSpinor *phi, DiracSpinor *Freg,
                              DiracSpinor *Firr, double *K,
                              const DiracSpinor &Fa, double omega,
                              const std::vector<double> &vl, double alpha,
                              const std::vector<DiracSpinor> &core,
                              const DiracSpinor &Fs, double eps_target = 1.0e-9,
                              const DiracSpinor *const Fhole = nullptr);

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
