#pragma once
#include "Coulomb/meTable.hpp"
#include <iostream>
#include <utility>
#include <vector>
class DiracSpinor;
namespace DiracOperator {
class TensorOperator;
}
namespace ExternalField {
class CorePolarisation;
class TDHF;
} // namespace ExternalField
namespace MBPT {
class CorrelationPotential;
}

namespace Amplitudes {

/*!
 Second-order (in the external field) amplitudes for a single-valence atom.

 For a transition \f$ a \to b \f$ due to two one-body operators, \f$ t \f$
 (rank \f$ k_t \f$, frequency \f$ \omega \f$) and \f$ s \f$ (rank
 \f$ k_s \f$, frequency \f$ \omega_s \f$), the amplitude with definite
 projections \f$ q_1, q_2 \f$ of the two operators is

 \f[
   A^{k_tk_s}_{q_1q_2} = \sum_n \left[
     \frac{\matel{b}{t_{q_1}}{n}\matel{n}{s_{q_2}}{a}}{\en_a + \omega_s - \en_n}
   + \frac{\matel{b}{s_{q_2}}{n}\matel{n}{t_{q_1}}{a}}{\en_a + \omega - \en_n}
   \right],
 \f]

 the sum over \f$ n \f$ running over the magnetic quantum numbers too, and
 \f$ \en_b = \en_a + \omega + \omega_s \f$. The operators are coupled to rank
 \f$ K \f$, with \f$ Q = q_1+q_2 = m_b-m_a \f$ and \f$ [K] \equiv 2K+1 \f$,

 \f[
   A^K_Q = \sum_{q_1q_2}\braket{k_tq_1\,k_sq_2}{KQ}\,A^{k_tk_s}_{q_1q_2}
         = (-1)^{k_t-k_s+Q}\sqrt{[K]}\sum_{q_1q_2}
           \threej{k_t}{k_s}{K}{q_1}{q_2}{-Q}\,A^{k_tk_s}_{q_1q_2},
 \f]

 and \f$ A^K \f$ follows from the Wigner-Eckart theorem,

 \f[
   A^K_Q = (-1)^{j_b-m_b}\threej{j_b}{K}{j_a}{-m_b}{Q}{m_a}\,A^K
   \qquad {\rm with} \qquad
   A^K = \sum_n \left[
     c_1(j_n)\,
     \frac{\redmatel{b}{t}{n}\redmatel{n}{s}{a}}{\en_a + \omega_s - \en_n}
   + c_2(j_n)\,
     \frac{\redmatel{b}{s}{n}\redmatel{n}{t}{a}}{\en_a + \omega - \en_n}
   \right],
 \f]

 the coefficients being those of CI::A_K_coefs (evaluated with the
 single-particle j in place of J), which also gives the sign convention of the
 coupling and the specific cases; CI::z_component converts \f$ A^K \f$ to the
 z-component. With \f$ t = s = d \f$ (E1) and \f$ [j] \equiv 2j+1 \f$:

 \f[
   \alpha_0 = \frac{A^0}{\sqrt{3[j_a]}},
   \qquad
   \alpha_2 = -\sqrt{\frac{2j(2j-1)}{3(j+1)(2j+1)(2j+3)}}\;A^2,
   \qquad
   \beta = \frac{A^1}{\sqrt{2}\,\redmatel{b}{\bm\sigma}{a}},
 \f]

 for \f$ K = 0, 2, 1 \f$ (\f$ \alpha_2 \f$ requires
 \f$ j_b = j_a = j \ge 1 \f$), and with \f$ t = d \f$, \f$ s = h_W \f$ (PNC,
 \f$ k_s = 0 \f$, so \f$ K = 1 \f$),

 \f[
   E_{\rm PNC} = A^1_0
     = (-1)^{j_b-m}\threej{j_b}{1}{j_a}{-m}{0}{m}\,A^1 .
 \f]

 This is the single-valence analogue of CI::A_K; it covers static, dynamic,
 and transition polarisabilities (\f$ t = s = E1 \f$), and PNC amplitudes
 (\f$ s \f$ = PNC operator).

 Two methods are provided for the valence sum: sum-over-states (SOS) over a
 given spectrum, and mixed states (MS, solving the inhomogeneous equation via
 ExternalField::TDHF), which is complete (no truncation of the sum). The
 contribution of core excitations (closed core: K = 0, diagonal only) is
 likewise available both ways. SOS and MS must agree (to basis completeness);
 the comparison is a strong check on the numerics.
*/

//! Is rank K allowed: triangle rules for the operators (kt, ks) and states
[[nodiscard]] bool allowed_K(int K, int kt, int ks, int twoJb, int twoJa);

//! The smallest rank K allowed for the amplitude; negative if there is none
[[nodiscard]] int smallest_allowed_K(int kt, int ks, int twoJb, int twoJa);

/*!
  @brief Valence part of the second-order amplitude A^K, by sum-over-states.
  @details
  Evaluates the sum above directly over the states of @p spectrum.

  If the spectrum contains states below the Fermi level (e.g.
  Wavefunction::spectrum() does), the core-valence (Pauli blocking) part of
  the amplitude is included automatically through those terms, as in the
  polarisability module. The polarisation of the closed core is separate:
  see @ref sos_core.

  Matrix elements: taken from @p t_me / @p s_me if present in the table,
  otherwise calculated directly as \f$ \redmatel{}{h}{} + \delta V \f$. The
  tables (if given) should be formed at the operator's frequency, and may
  contain RPA and structure radiation. RPA enters on both vertices, as
  \f$ \delta V(\omega) \f$ of the operator's own frequency at every
  vertex (the emission component for a negative frequency; see
  CorePolarisation::dV), whichever side of the leg the intermediate state
  lies.

  @param K       Rank of the amplitude.
  @param Fb,Fa   Final and initial valence states.
  @param t,s     The two operators.
  @param omega   Frequency of t. For a real transition carried entirely by t
                 this is e_b - e_a, and s is static.
  @param omega_s Frequency of s. Energy conservation requires
                 omega + omega_s = e_b - e_a; it is -omega for a dynamic
                 polarisability (b = a).
  @param spectrum  Intermediate states summed over.
  @param dVt,dVs RPA for each operator, solved at the frequency of that
                 operator. May be nullptr. Ignored for table entries.
  @param t_me,s_me  Optional tables of single-particle reduced matrix
                 elements; empty (default) to calculate directly.
  @return \f$ A^K \f$ (valence part).

  @note Degenerate denominators arise for even-parity operator pairs (an
        intermediate state degenerate with \f$ \en_a + \omega \f$; e.g.
        n = a for a static diagonal amplitude). Such terms diverge, and
        nothing is skipped here: remove the offending states from the sum
        and treat them separately (cf. the project-out treatment in the
        pnc module).
*/
[[nodiscard]] double
sos_valence(int K, const DiracSpinor &Fb, const DiracSpinor &Fa,
            const DiracOperator::TensorOperator *t,
            const DiracOperator::TensorOperator *s, double omega,
            double omega_s, const std::vector<DiracSpinor> &spectrum,
            const ExternalField::CorePolarisation *dVt = nullptr,
            const ExternalField::CorePolarisation *dVs = nullptr,
            const Coulomb::meTable<double> &t_me = {},
            const Coulomb::meTable<double> &s_me = {});

/*!
  @brief Contribution to A^K from the polarisation of the closed core, by
  sum-over-states.
  @details
  Delegates to CI::A_K_core, which is the same quantity: a core electron
  excited by one operator and de-excited by the other,

  \f[
    A^0_{\rm core} = \sqrt{[J]}\sum_c \sqrt{[j_c]} \sum_m \left[
      c_1\,\frac{\redmatel{c}{t}{m}\redmatel{m}{s}{c}}{\en_c+\omega_s-\en_m}
    + c_2\,\frac{\redmatel{c}{s}{m}\redmatel{m}{t}{c}}{\en_c+\omega-\en_m}
    \right].
  \f]

  Non-zero only for K = 0 (closed core), which requires kt = ks, and only
  for a diagonal amplitude (b = a); the caller is responsible for the
  diagonal condition.

  @param K       Rank; zero returned unless K = 0.
  @param twoJ    2J of the valence state (the sqrt([J]) prefactor).
  @param t,s     The two operators.
  @param omega,omega_s  Frequency of each operator; see @ref sos_valence.
  @param core    Core states c.
  @param excited Particle states m: basis states above the Fermi level.
  @param dVt,dVs RPA for each operator. May be nullptr.
  @return \f$ A^0_{\rm core} \f$.

  @note RPA enters once: of the two matrix elements, only the one acting on
        the core orbital is dressed (see the note on CI::A_K_core).
*/
[[nodiscard]] double
sos_core(int K, int twoJ, const DiracOperator::TensorOperator *t,
         const DiracOperator::TensorOperator *s, double omega, double omega_s,
         const std::vector<DiracSpinor> &core,
         const std::vector<DiracSpinor> &excited,
         const ExternalField::CorePolarisation *dVt = nullptr,
         const ExternalField::CorePolarisation *dVs = nullptr);

/*!
  @brief Valence part of the second-order amplitude A^K, evaluated with
  mixed states (TDHF method).
  @details
  The sums over intermediate states are performed with mixed states
  (solutions of the inhomogeneous Dirac equation, via
  ExternalField::TDHF::solve_dPsi), so they are complete: no truncation of
  the spectrum. All intermediate states of a given kappa share the same
  angular coefficient, so one mixed state per kappa channel per term is
  required.

  Each sum is formed in two independent ways: with the mixed states of
  \f$ s \f$, and with those of \f$ t \f$, returned as the two elements of
  the pair. They must agree; the difference is a check on the numerics
  (cf. CI::A_K).

  RPA enters on both vertices: the mixed states are solved with
  \f$ t + \delta V \f$ (the inner vertex), and the outer matrix element is
  dressed with the \f$ \delta V \f$ of the outer operator; each is
  \f$ \delta V(\omega) \f$ at the frequency of its own operator (cf.
  @ref sos_valence). With RPA not solved, \f$ \delta V = 0 \f$ and the
  amplitude is at the HF level.

  @param K       Rank of the amplitude.
  @param Fb,Fa   Final and initial valence states.
  @param t,s     The two operators.
  The core-valence (Pauli blocking) part of the sum is included (the mixed
  states are complete). It may be separated by projecting the mixed states
  onto the span of the core states: pass the HF core as @p project_onto to
  obtain just that part (cf. the orthogonality treatment in the pnc module).

  @param omega,omega_s  Frequency of each operator; see @ref sos_valence.
  @param dVt,dVs TDHF object for each operator (required, not null): it
                 provides the mixed-state solver even when RPA is off
                 (unsolved TDHF gives dV = 0). Solve at the operator's
                 frequency for RPA. Must be TDHF (or TDHFbasis): the
                 diagram method cannot drive mixed states.
  @param Sigma   Optional correlation potential, included in the
                 mixed-state solutions (use with Brueckner valence states).
  @param project_onto  If non-empty, each mixed state is projected onto the
                 span of these states before the outer matrix element; pass
                 the HF core for the core-valence part of the amplitude.
                 Empty (default): the full sum.
  @param outstream  Stream for warnings.
  @return The two evaluations of \f$ A^K \f$: with the mixed states of s,
          and with those of t.

  @note If an intermediate state is degenerate with a denominator, the
        mixed-state equation is singular and the amplitude divergent. This
        cannot happen for odd-parity operators (the intermediate states
        have opposite parity to a and b); see the note on CI::A_K.
*/
[[nodiscard]] std::pair<double, double>
ms_valence(int K, const DiracSpinor &Fb, const DiracSpinor &Fa,
           const DiracOperator::TensorOperator *t,
           const DiracOperator::TensorOperator *s, double omega, double omega_s,
           const ExternalField::TDHF *dVt, const ExternalField::TDHF *dVs,
           const MBPT::CorrelationPotential *Sigma = nullptr,
           const std::vector<DiracSpinor> &project_onto = {},
           std::ostream &outstream = std::cout);

/*!
  @brief Contribution to A^K from the polarisation of the closed core,
  evaluated with mixed states.
  @details
  As @ref sos_core, with the sum over excited states m performed with mixed
  states of the core orbitals (generalises the TDHF core polarisability of
  the polarisability module to two operators). Non-zero only for K = 0 and
  a diagonal amplitude (caller's responsibility).

  RPA enters once, through the mixed-state solve (the vertex acting on the
  core orbital); the outer matrix element is bare - the same counting as
  @ref sos_core.

  The mixed states of a core orbital include (occupied) core intermediate
  states; for the diagonal amplitude (omega_s = -omega) their contributions
  cancel in the sum over the core, so this agrees with @ref sos_core, which
  sums over excited states only.

  @param K       Rank; zero returned unless K = 0.
  @param twoJ    2J of the valence state.
  @param t,s     The two operators.
  @param omega,omega_s  Frequency of each operator.
  @param core    Core states c.
  @param dVt,dVs TDHF object for each operator (required); see
                 @ref ms_valence.
  @param Sigma   Optional correlation potential for the mixed states. The
                 intermediate states here are core excitations, for which
                 the valence Sigma is not appropriate: normally nullptr
                 (matching @ref sos_core, which uses the HF basis).
  @return \f$ A^0_{\rm core} \f$.
*/
[[nodiscard]] double
ms_core(int K, int twoJ, const DiracOperator::TensorOperator *t,
        const DiracOperator::TensorOperator *s, double omega, double omega_s,
        const std::vector<DiracSpinor> &core, const ExternalField::TDHF *dVt,
        const ExternalField::TDHF *dVs,
        const MBPT::CorrelationPotential *Sigma = nullptr);

} // namespace Amplitudes
