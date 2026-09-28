#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <cstddef>
#include <vector>
namespace HF {
class Breit;
}

//! External field: the channel solvers of TDHFcntm (mixed states above the
//! ionisation threshold, and the bound channels solved alongside them), and
//! the Anderson mixing they use
namespace ExternalField {

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
  @brief Anderson mixing history of a fixed-point iteration x -> G(x) on
  spinors, as used by the mixed-state solvers.
  @details
  Holds the map outputs \f$ g_k = G(x_k) \f$, the residuals
  \f$ r_k = g_k - x_k \f$ and their Gram matrix, and a scalar carried with
  each entry (an amplitude linear in the state, mixed with the same
  coefficients; zero if unused). The Gram matrix is kept incrementally: a
  new entry costs one inner product per kept entry, not a rebuild of the
  whole matrix. At most @p depth entries are kept: the oldest is dropped
  when the history is full, and (see @ref anderson_coefficients) while the
  Gram matrix is ill-conditioned.
*/
class AndersonHistory {
public:
  explicit AndersonHistory(std::size_t depth) : m_depth(depth) {
    assert(depth > 0 && "AndersonHistory needs a depth of at least one");
  }

  //! Adds the entry (g, r), with |r|^2 = r2 (which the caller has from its
  //! convergence test) and the carried scalar k; drops the oldest entry if
  //! the history is full
  void push(DiracSpinor g, DiracSpinor r, double r2, double k = 0.0);

  /*!
    @brief The Anderson extrapolation over the kept entries,
    \f$ x = \sum_k c_k g_k \f$ and (if @p k is non-null)
    \f$ k = \sum_k c_k k_k \f$, with the coefficients of
    @ref anderson_coefficients; entries that had to be dropped for
    conditioning are dropped from the history too.
    @details If no entry is usable (the coefficients could not be solved
    for), x is the plain step: the latest map output, and the history is
    cleared.
  */
  void mix(DiracSpinor *x, double *k = nullptr);

  //! Number of kept entries
  std::size_t size() const { return m_g.size(); }

private:
  std::size_t m_depth;
  std::vector<DiracSpinor> m_g{};
  std::vector<DiracSpinor> m_r{};
  std::vector<double> m_k{};
  // Gram matrix of the residuals, <r_i|r_j>, over the kept entries
  std::vector<std::vector<double>> m_gram{};

  // Drops the n oldest entries (and their Gram rows and columns)
  void drop_oldest(std::size_t n);
};

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
  converges on the conditioning of the true operator. The iteration also
  ends when the residual stops improving (the round-off floor of the solve,
  which can lie above a tight @p eps_target): the best iterate is then
  returned, with its residual in dF.eps(). Parameters as
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
  is mixed with the same coefficients. The iteration also ends when the
  residual stops improving (the round-off floor of the solve, which can lie
  above a tight @p eps_target): the best iterate is then returned, with its
  residual in phi->eps(). After each solve, phi is
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

} // namespace ExternalField
