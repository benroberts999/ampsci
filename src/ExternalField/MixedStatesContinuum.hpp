#pragma once
#include "ExternalField/Anderson.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

namespace ExternalField {

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
