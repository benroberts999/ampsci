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
  @brief Bound mixed-state solve for the continuum (photoionisation) TDHF path:
  Anderson/Pulay-accelerated, with optional V^{N-1} hole term.
  @details
  Same physics and conditioning as @ref solveMixedState (in-place overload), but
  solves the LINEAR equation \f$ (h_{\rm HF} - \en_0)\,\delta F = -F_S \f$ by
  Anderson/Pulay mixing (equivalent to preconditioned GMRES) instead of damped
  fixed-point iteration. At the high frequencies of photoionisation a bound
  channel's \f$ \en_0 = \en_a + \omega \f$ can sit near a SPURIOUS eigenvalue of
  the local preconditioner \f$ (h_{\rm local} + U_x - \en_0) \f$, where damped
  iteration diverges even though \f$ (h_{\rm HF} - \en_0) \f$ is non-singular and
  the response is finite; Anderson mixing converges on the true operator's
  conditioning, independent of those spurious modes.

  Kept separate from @ref solveMixedState so the standard bound RPA path is
  untouched. Used by @ref TDHF::solve_core_cntm for the closed channels and the
  bound (Y/-) partner of an ionised orbital.

  @note If @p Fhole is given, the iterated non-local exchange is
  \f$ V^{\rm exch} - X_{\rm hole} \f$ with \f$ X_{\rm hole} \f$ the one-electron
  self-exchange (@ref HF::vexFa_1el) -- the exchange part of the V^{N-1}
  (residual ion) Hamiltonian. The caller must pair this with the direct part,
  \f$ v_l \to v_l - y^0_{\rm hole,hole} \f$, in @p vl (do both or neither).
*/
void solveMixedState_cntm(DiracSpinor &dF, const DiracSpinor &Fa, double omega,
                          const std::vector<double> &vl, double alpha,
                          const std::vector<DiracSpinor> &core,
                          const DiracSpinor &Fs, double eps_target = 1.0e-9,
                          const HF::Breit *const VBr = nullptr,
                          const std::vector<double> &H_mag = {},
                          const DiracSpinor *const Fhole = nullptr);

/*!
  @brief Solves the continuum (en>0) mixed-state equation by forward
  integration (bare local solve: no exchange, no dV).
  @details
  For photoionisation the X/+ response energy en_+ = en_a + omega can be
  positive (continuum). Then the bound inhomogeneous solver is singular and is
  replaced by a continuum solve: solves
  \f[ (h_r^{(\kappa_\alpha)} - \en_+)\,\varphi_+ = -F_S \f]
  with the standing-wave (K-matrix) boundary condition, by forward (outward)
  integration plus F_reg subtraction, using the two real homogeneous solutions
  at en_+ -- the regular continuum orbital F_reg (energy-normalised) and its
  irregular partner F_irr -- in place of the bound (regular-at-origin /
  decaying-at-infinity) pair. Internally calls @ref DiracODE::solveContinuum,
  @ref DiracODE::solveContinuumIrregular and
  @ref DiracODE::solveContinuumForward. (The global continuum Green's function
  is avoided: F_irr diverges at the origin like r^{-l-1}; the outward-integrated
  solution is regular there, and F_irr is only sampled in the outer window for
  the c,K extraction.)

  This overload exposes the homogeneous solutions F_reg, F_irr (out-params);
  @p phi.kappa() selects the continuum channel kappa_alpha and must equal
  @p Fs.kappa().

  @param phi    Output: continuum correction orbital (kappa = continuum channel).
  @param Freg   Output: regular energy-normalised continuum orbital at en_+.
  @param Firr   Output: irregular partner.
  @param Fa     Bound orbital phi_a (provides en_a and the grid).
  @param omega  External-field frequency (en_+ = Fa.en() + omega must be > 0).
  @param vl     Local potential v(r) (homogeneous solutions are eigenstates of it).
  @param alpha  Fine-structure constant.
  @param Fs     Source spinor F_S (note sign: +h*phi_a, as in solveMixedState).

  @warning Requires en_+ > 0 and a grid dense enough at large r (see
           @ref DiracODE::solveContinuum); F_reg is zeroed if too sparse.
*/
void solveContinuumMixedState(DiracSpinor &phi, DiracSpinor &Freg,
                              DiracSpinor &Firr, const DiracSpinor &Fa,
                              double omega, const std::vector<double> &vl,
                              double alpha, const DiracSpinor &Fs);

/*!
  @brief Continuum (en>0) mixed-state solve with nonlocal exchange iterated in
  the source: the continuum X/+ channel for the TDHF iteration.
  @details
  As solveContinuumMixedState() above, but with Hartree-Fock exchange iterated
  to self-consistency, exactly as in the bound @ref solveMixedState. Solves
  \f[ (h_r^{(\kappa_\alpha)} + V^{\rm nl} - \en_+)\,\varphi_+ = -F_S \f]
  with the standing-wave (K-matrix) boundary condition. The orbital-independent
  Kohn-Sham local exchange @ref HF::vex_KS is folded into the local potential
  \f$ v = v_l + U_x \f$ for conditioning (cancelled on the right at the fixed
  point); the nonlocal exchange remainder is carried in the source and iterated.
  Because \f$ U_x \f$ is orbital-independent, v is fixed across iterations, so
  F_reg and F_irr are built only once. (vex_approx is NOT used here: it divides
  by the oscillatory continuum orbital and spikes at its nodes.) After each
  solve \f$ \varphi_+ \f$ is orthogonalised to the same-kappa occupied core
  orbitals.

  The hole-particle (self-interaction) treatment is the caller's: pass @p vl
  already adjusted (e.g. \f$ v_l - y^0_{aa} \f$), and pass the hole orbital as
  @p Fhole so the iterated exchange becomes \f$ V^{\rm exch} - X_a \f$ (with
  \f$ X_a \f$ the one-electron self-exchange, @ref HF::vexFa_1el). Together
  these put the photoelectron in the V^{N-1} Hamiltonian (do both or neither).

  @param phi    Output: continuum correction orbital. Used as the starting guess
                if nonzero (faster at a nearby omega).
  @param Freg   Output: regular energy-normalised continuum orbital at en_+.
  @param Firr   Output: irregular partner, final iteration.
  @param Fa     Bound orbital phi_a (provides en_a and the grid).
  @param omega  External-field frequency (en_+ = Fa.en() + omega must be > 0).
  @param vl     Local potential v_l(r) (caller includes any hole-particle term).
  @param alpha  Fine-structure constant.
  @param core   Core orbitals (nonlocal exchange and orthogonalisation).
  @param Fs     Source spinor F_S (note sign: +h*phi_a).
  @param eps_target  Convergence goal for the exchange self-consistency.
  @param Fhole  Optional hole orbital: subtract its one-electron self-exchange
                from the iterated non-local exchange (V^{N-1} exchange part).

  @note If @p core is empty, reduces to a single bare local continuum (forward)
        solve.
  @warning Requires en_+ > 0 and a grid dense enough at large r.
*/
void solveContinuumMixedState(DiracSpinor &phi, DiracSpinor &Freg,
                              DiracSpinor &Firr, const DiracSpinor &Fa,
                              double omega, const std::vector<double> &vl,
                              double alpha,
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
