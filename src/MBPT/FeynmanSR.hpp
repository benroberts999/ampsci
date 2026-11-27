#pragma once
#include "LinAlg/Matrix.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include <complex>
#include <optional>
#include <utility>
#include <vector>

namespace DiracOperator {
class TensorOperator;
}

namespace ExternalField {
class CorePolarisation;
}

namespace MBPT {

//------------------------------------------------------------------------------
/*!
  @brief Structure radiation diagrams via the Feynman technique.

  @details
  Uses Feynman diagrams to calculate the structure radiation (vertex) correction 
  for general one-body operators, using the existing @ref MBPT::Feynman class.

  The correction to the matrix element, due to externat field operator \f$ h \f$ 
  with frequency \f$ \omega \f$ is:
  \f[ 
    \bra{w}\delta\Sigma(\varepsilon_v)\ket{v} 
  \f],
  with
  \f$ \varepsilon_w = \varepsilon_v + \omega \f$. 
  To simplify, only this (absorption) form is used; omega may be negative.
  Other matrix element can be obtained by Hermicity,

  The perturbed Green's function is
  \f[
    \delta G_{12}(\varepsilon) =
      G_{1i}(\varepsilon + \omega)\, h_i\, G_{i2}(\varepsilon) .
  \f]
  Coordinate integration is implied over repeated Latin indices.

  - dSigma_a  : field on the outer Green's function line
  - dSigma_bc : field on either Green's function line of the polarisation
                loop (perturbed polarisation operator, see dPi_b)
  - dSigma_ax_cx : exchange diagrams, field on the internal line adjacent
                to the final state w (a_x), or to the initial state v (c_x)
  - dSigma_bx : exchange diagram, field on the middle internal line (b_x)

  The diagrams are formed as matrices on the radial sub-grid, by
  integration along the imaginary energy axis.
  After angular reduction, each gives the reduced matrix element as
  `Fw * (dSigma * Fv)`.

  Screening of the Coulomb lines and the hole-particle interaction are
  included as set in FeynmanOptions (as for the correlation potential):
  - dSigma_a  : both, through Q*Pi*Q
  - dSigma_bc : both Coulomb lines screened; hole-particle in the
                perturbed polarisation loop
  - dSigma_ax_cx : one Coulomb line screened; hole-particle in the
                exchange operator
  - dSigma_bx : one Coulomb line screened; hole-particle in the
                perturbed exchange operator
  See each function for what is approximated.

  Sigma_SR_direct, Sigma_SR_direct2, and the functions they use are for the
  E1 operator only (hard-coded).
*/
class FeynmanSR : public Feynman {

public:
  using Feynman::Feynman;

  //----------------------------------------------------------------------------
  /*!
    @brief Direct SR diagram (a): 
    external field on regular Greens function of the direct diagram.

    @details
    \f[
      \delta\Sigma^{(a)}_{12}(\varepsilon_v) 
        = \int\frac{d\varepsilon}{2\pi}\,
        \delta G_{12}(\varepsilon_v + \varepsilon)\,
        Q_{1i}\,\Pi_{ij}(\varepsilon)\,Q_{j2} .
    \f]

    @param h        External field operator. If frequency dependent, must
                    already be set to the frequency
    @param kappa_w  Kappa of the final state, w (on the left)
    @param kappa_v  Kappa of the initial state, v (on the right)
    @param en_v     Energy of the initial state
    @param omega    Frequency of the external field,
                    \f$ \varepsilon_w - \varepsilon_v \f$
    @param dV       Optional: core polarisation (already solved, at the
                    frequency of the caller's choice). If given, the vertex
                    is h + dV (see dV_matrices). For negative @p omega, the
                    conjugate component of dV is used

    @return Matrix on the sub-grid; the reduced matrix element is
    `Fw * (dSigma * Fv)`.

    @note Forms Q*Pi*Q if it has not yet been formed (see Feynman::form_qpiq).

    @warning The operator must be local in r (see perturbed_greens).
  
    @note With screening and/or hole-particle (FeynmanOptions), these are
          included through Q*Pi*Q: the polarisation loop with the
          hole-particle interaction, and one Coulomb line screened.
  */
  GMatrix dSigma_a(const DiracOperator::TensorOperator &h, int kappa_w,
                   int kappa_v, double en_v, double omega,
                   const ExternalField::CorePolarisation *dV = nullptr);

  //----------------------------------------------------------------------------
  /*!
    @brief Direct SR diagrams (b) and (c): 
    external field inside the polarisation loop of the direct diagram.

    @details
    \f[
      \delta\Sigma^{(b+c)}_{12}(\varepsilon_v) = \int\frac{d\varepsilon}{2\pi}\,
        G_{12}(\varepsilon_v + \varepsilon)\,
        Q_{1i}\,\delta\Pi_{ij}(\varepsilon)\,Q_{j2},
    \f]
    with the perturbed polarisation operator
    \f$ \delta\Pi = \delta\Pi^{(b)} + \delta\Pi^{(c)} \f$ (see dPi_b), where
    \f[
      \delta\Pi^{(c)}_{12}(\varepsilon) =
        \delta\Pi^{(b)}_{21}(\omega - \varepsilon) .
    \f]

    Parameters as for dSigma_a.

    @return {dSigma^(b), dSigma^(c)}, from the two terms of
    \f$ \delta\Pi \f$: for each, the reduced matrix element is
    `Fw * (dSigma * Fv)`. The two nearly cancel.

    @warning The operator must be local in r (see perturbed_greens).

    @note Hole-particle (FeynmanOptions): included in the perturbed
          polarisation loop (see dPi_b).

    @warning Screening (FeynmanOptions): both Coulomb lines are the screened
          lines, taken at the integration energy. The line next to w is at
          that energy less omega: the difference is neglected (none for
          omega = 0).

    @note Forms Q*Pi*Q if screening is on and it has not yet been formed
          (see Feynman::form_qpiq).
  */
  std::pair<GMatrix, GMatrix>
  dSigma_bc(const DiracOperator::TensorOperator &h, int kappa_w, int kappa_v,
            double en_v, double omega,
            const ExternalField::CorePolarisation *dV = nullptr);

  //----------------------------------------------------------------------------
  /*!
    @brief Exchange SR diagrams (a_x) and (c_x): 
    field acts on internal line next to final (a_x), or initial (c_x) state.

    @details
    \f[
      \delta\Sigma^{(a_{\rm x})}_{12}(\varepsilon_v)
        = -\int\frac{d\varepsilon}{2\pi}\,
        \delta G_{1i}(\varepsilon_v + \varepsilon)\,
        Q_{1j}\,\Gamma_{ij2}(\varepsilon)\,Q_{i2} ,
    \f]
    with the exchange operator
    \f[
      \Gamma_{123}(\varepsilon) = i\sum_a \left[
        P^a_{12}\, G^{\rm e}_{23}(\varepsilon_a - \varepsilon)
        + G^{\rm e}_{12}(\varepsilon_a + \varepsilon)\, P^a_{23} \right],
    \f]
    summed over the polarised core states a.

    (c_x) is the mirror image of (a_x): the same diagram with w and v
    swapped and \f$ \omega \to -\omega \f$, transposed, times the symmetry
    sign of the operator (TensorOperator::symm_sign). The two need the same
    Green's functions, so are formed together. For a static field
    (\f$ \omega = 0 \f$) between the same partial waves the mirror image is
    (a_x) itself: (c_x) is then the transpose of (a_x), and is not
    calculated separately.

    Screening (FeynmanOptions): the line \f$ Q_{1j} \f$ is the screened
    line (it is at the integration energy). \f$ Q_{i2} \f$ is at the
    energy of the integral inside \f$ \Gamma \f$, and stays bare.

    Hole-particle (FeynmanOptions): the excited Green's functions of
    \f$ \Gamma \f$ include the interaction with the hole in a.

    @param h        External field operator
    @param kappa_w  Kappa of the final state (on the left)
    @param kappa_v  Kappa of the initial state (on the right)
    @param en_v     Energy of the initial state
    @param omega    Frequency of the external field
    @param dV       Optional: core polarisation; vertex is then h + dV
                    (see dSigma_a)

    @return {dSigma^(ax), dSigma^(cx)}: for each, the reduced matrix
    element is `Fw * (dSigma * Fv)`.

    @note Local operators only (see perturbed_greens).
    @note The g components of the Green's functions are contracted only if
          FeynmanOptions::include_G_exchange (as for Feynman::Sigma_exchange).
    @note Forms Q*Pi*Q if screening is on and it has not yet been formed
          (see Feynman::form_qpiq).
  */
  std::pair<GMatrix, GMatrix>
  dSigma_ax_cx(const DiracOperator::TensorOperator &h, int kappa_w, int kappa_v,
               double en_v, double omega,
               const ExternalField::CorePolarisation *dV = nullptr);

  //----------------------------------------------------------------------------
  /*!
    @brief Perturbed Green's function (see perturbed_greens) for each pair
    of partial waves, with its angular factors, summed over the partial
    wave on the final-state side. For the line next to the final state in
    the exchange diagrams.

    @details
    For each partial wave \f$ \beta \f$ (next to the exchange operator) and
    multipolarity k of the Coulomb line at the final state:
    \f[
      \delta G^k_\beta = \sum_\alpha L\, C^k_{w\alpha}\,
        G_\alpha\, \langle\alpha||h||\beta\rangle\, G_\beta ,
    \f]
    with L from L_ang_dSR.

    @param h        External field operator
    @param kappa_w  Kappa of the final state (on the left)
    @param kappa_v  Kappa of the initial state (on the right)
    @param G_alpha  Green's functions to the left of h, for each kappa index
    @param G_beta   Green's functions to the right of h, for each kappa index
    @param dV_kappa Core polarisation matrices, from dV_matrices

    @return Indexed (kappa index of beta, k); empty where zero by selection
    rules.
  */
  LinAlg::Matrix<std::optional<ComplexGMatrix>> perturbed_greens_summed(
    const DiracOperator::TensorOperator &h, int kappa_w, int kappa_v,
    const std::vector<ComplexGMatrix> &G_alpha,
    const std::vector<ComplexGMatrix> &G_beta,
    const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dV_kappa) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Contribution of one core state a, at one integration energy, to
    the integrand of the exchange diagram (a_x) (see dSigma_ax_cx).

    @details
    \f[
      {\rm Re}\left[ i\,d\varepsilon\;
        \delta G_{1i}\, Q_{1j}\, \Gamma^a_{ij2}(\varepsilon)\, Q_{i2} \right],
    \f]
    with
    \f[
      \Gamma^a_{123}(\varepsilon) =
        P^a_{12}\, G^{\rm e}_{23}(\varepsilon_a - \varepsilon)
        + G^{\rm e}_{12}(\varepsilon_a + \varepsilon)\, P^a_{23} .
    \f]

    @param dG_beta_k    Perturbed Green's function (perturbed_greens_summed)
    @param kappa_v      Kappa of the initial state (on the right)
    @param Fa           Core state a
    @param G_delta      \f$ G^{\rm e}(\varepsilon_a - \varepsilon) \f$, for
                        each kappa index
    @param G_gamma      \f$ G^{\rm e}(\varepsilon_a + \varepsilon) \f$, for
                        each kappa index
    @param ie           Index of the integration energy (for the screened
                        Coulomb line, see dQ_screening)
    @param i_d_epsilon  i times the integration weight \f$ d\varepsilon \f$
  */
  GMatrix dSigma_ax_core_term(
    const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dG_beta_k, int kappa_v,
    const DiracSpinor &Fa, const std::vector<ComplexGMatrix> &G_delta,
    const std::vector<ComplexGMatrix> &G_gamma, std::size_t ie,
    std::complex<double> i_d_epsilon) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Exchange structure radiation diagram (b_x), for any external
    field operator: the field acts on the middle internal line.

    @details
    \f[
      \delta\Sigma^{(b_{\rm x})}_{12}(\varepsilon_v)
        = -\int\frac{d\varepsilon}{2\pi}\,
        G_{1i}(\varepsilon_v + \varepsilon)\,
        Q_{1j}\,\delta\Gamma_{ij2}(\varepsilon)\,Q_{i2} ,
    \f]
    with the perturbed exchange operator
    \f[
      \delta\Gamma_{123}(\varepsilon) =
        \int\frac{d\varepsilon'}{2\pi}\,
        \delta G_{12}(\varepsilon + \varepsilon')\,
        G_{23}(\varepsilon' + \omega) .
    \f]
    The inner integral is the sum of the residues at the core poles of
    each of its three Green's functions:
    \f[
      \delta\Gamma_{123}(\varepsilon) = i\sum_a \Big[
        [P^a h\, G(\varepsilon_a - \omega)]_{12}\,
          G_{23}(\varepsilon_a - \varepsilon)
        + [G(\varepsilon_a + \omega)\, h P^a]_{12}\,
          G_{23}(\varepsilon_a - \varepsilon + \omega)
        + [G(\varepsilon_a + \varepsilon)\, h\,
          G(\varepsilon_a + \varepsilon - \omega)]_{12}\, P^a_{23}
      \Big],
    \f]
    summed over the polarised core states a, with
    \f$ G = G^{\rm e} + G^{\rm c} \f$ (see green_loop).

    The parts with three core lines sum to zero. For operators with large
    core matrix elements (e.g., hyperfine) each is much larger than the
    diagram, so they are subtracted from every term: the result does not
    rely on their cancellation.

    If h joins a core state to itself, the first two terms each contain a
    part with \f$ P^a h P^a \f$ that has the factor \f$ 1/\omega \f$
    (as in dPi_b). The two are replaced by the single term
    \f[
      -[P^a h P^a]_{12}\,
      \left[G^{\rm e}(\varepsilon_a - \varepsilon)\,
        G^{\rm e}(\varepsilon_a - \varepsilon + \omega)\right]_{23} .
    \f]

    In the first two terms the perturbed line does not depend on the
    integration energy: it is the core state times a perturbed core state,
    \f$ X_a = G(\varepsilon_a + \omega)\, h F_a \f$ and
    \f$ Y_a = F_a h\, G(\varepsilon_a - \omega) \f$, formed once.

    Screening (FeynmanOptions): as for dSigma_ax_cx, the line
    \f$ Q_{1j} \f$ is the screened line, \f$ Q_{i2} \f$ stays bare. The
    screened line is taken at the integration energy
    \f$ \varepsilon \f$ (it is at \f$ \varepsilon - \omega \f$: the
    difference is neglected; none for \f$ \omega = 0 \f$).

    Hole-particle (FeynmanOptions): the excited part of every Green's
    function of \f$ \delta\Gamma \f$ includes the interaction with the
    hole in a.

    @param h        External field operator
    @param kappa_w  Kappa of the final state (on the left)
    @param kappa_v  Kappa of the initial state (on the right)
    @param en_v     Energy of the initial state
    @param omega    Frequency of the external field
    @param dV       Optional: core polarisation; vertex is then h + dV
                    (see dSigma_a)

    @return dSigma^(bx): the reduced matrix element is
    `Fw * (dSigma * Fv)`.

    @note Local operators only (see perturbed_greens).
    @note The g components of the Green's functions are contracted only if
          FeynmanOptions::include_G_exchange (as for Feynman::Sigma_exchange).
    @note Forms Q*Pi*Q if screening is on and it has not yet been formed
          (see Feynman::form_qpiq).
    @note Converges more slowly with the sub-grid than (a_x): the integrand
          falls only as the inverse square of the integration energy.
  */
  GMatrix dSigma_bx(const DiracOperator::TensorOperator &h, int kappa_w,
                    int kappa_v, double en_v, double omega,
                    const ExternalField::CorePolarisation *dV = nullptr);

  //----------------------------------------------------------------------------
  /*!
    @brief Term of the exchange diagram (b_x) (see dSigma_bx) in which the
    perturbed line is a product of a function of each of its coordinates;
    one integration energy.

    @details
    \f[
      {\rm Re}\left[ i\,d\varepsilon\;
        G_{1i}\, Q_{1j}\, [F_i(r_i)\, F_j(r_j)]\, G_{j2}\, Q_{i2} \right],
    \f]
    summed over the partial waves of \f$ G_{1i} \f$ and \f$ G_{j2} \f$
    and the multipolarities of the two Coulomb lines.

    @param K              Rank of the external field
    @param kappa_w        Kappa of the final state (on the left)
    @param kappa_v        Kappa of the initial state (on the right)
    @param kappa_gamma_i  Partial wave of the perturbed line at vertex i
    @param kappa_gamma_j  Partial wave of the perturbed line at vertex j
    @param F_i            The function of \f$ r_i \f$
    @param Fj_G_delta     The line from j to 2 with the function of
                          \f$ r_j \f$ attached (see spinor_times_greens),
                          for each kappa index; empty for those not
                          included
    @param G_beta         Green's function of the line next to w, for each
                          kappa index
    @param ie             Index of the integration energy (for the screened
                          Coulomb line, see dQ_screening)
    @param i_d_epsilon    i times the integration weight
                          \f$ d\varepsilon \f$
  */
  GMatrix dSigma_bx_separable_term(
    int K, int kappa_w, int kappa_v, int kappa_gamma_i, int kappa_gamma_j,
    const DiracSpinor &F_i,
    const std::vector<std::vector<ComplexRMatrix>> &Fj_G_delta,
    const std::vector<ComplexGMatrix> &G_beta, std::size_t ie,
    std::complex<double> i_d_epsilon) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Term of the exchange diagram (b_x) (see dSigma_bx) in which the
    line next to the initial state is the core state a; one integration
    energy.

    @details
    \f[
      {\rm Re}\left[ i\,d\varepsilon\;
        G_{1i}\, Q_{1j}\, \delta G_{ij}\, P^a_{j2}\, Q_{i2} \right],
    \f]
    summed over the partial waves of \f$ G_{1i} \f$ and
    \f$ \delta G_{ij} \f$ and the multipolarities of the two Coulomb
    lines.

    @param K            Rank of the external field
    @param kappa_w      Kappa of the final state (on the left)
    @param kappa_v      Kappa of the initial state (on the right)
    @param Fa           Core state a
    @param dG           Perturbed line (see perturbed_greens), indexed by
                        the kappa indices of its partial waves at vertices
                        i and j; empty where zero by selection rules
    @param G_beta       Green's function of the line next to w, for each
                        kappa index
    @param ie           Index of the integration energy (for the screened
                        Coulomb line, see dQ_screening)
    @param i_d_epsilon  i times the integration weight \f$ d\varepsilon \f$
  */
  GMatrix
  dSigma_bx_core_term(int K, int kappa_w, int kappa_v, const DiracSpinor &Fa,
                      const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dG,
                      const std::vector<ComplexGMatrix> &G_beta, std::size_t ie,
                      std::complex<double> i_d_epsilon) const;

  //----------------------------------------------------------------------------
  /*!
    @brief A Green's function with a radial function attached to its first
    coordinate.

    @details
    \f[
      [F\,G]^{\nu}_{12} = \sum_{\mu} F^{\mu}(r_1)\, G^{\mu\nu}_{12} ,
    \f]
    with no coordinate integration. Spinor components as in the exchange
    diagrams (f only, unless FeynmanOptions::include_G_exchange).

    @param F  Radial function (on the full grid)
    @param G  Green's function

    @return One radial matrix for each spinor component \f$ \nu \f$.
  */
  std::vector<ComplexRMatrix>
  spinor_times_greens(const DiracSpinor &F, const ComplexGMatrix &G) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Perturbed polarisation operator, diagram (b): the polarisation
    loop with the external field on the line running from vertex 1 to 2.

    @details
    \f[
      \delta\Pi^{(b)}_{12}(\varepsilon) =
        {\rm Tr}\int\frac{d\varepsilon'}{2\pi}\,
        \delta G_{21}(\varepsilon' + \varepsilon)\,
        G_{12}(\varepsilon' + \omega) .
    \f]
    The integral is the sum of the residues at the core poles of each of
    the three Green's functions:
    \f[
      \delta\Pi^{(b)}_{12}(\varepsilon) = i\,{\rm Tr}\sum_a \Big[
        P^a_{12}\, G_{2i}(\varepsilon_a + \varepsilon)\, h_i\,
          G_{i1}(\varepsilon_a + \varepsilon - \omega)
        + G_{12}(\varepsilon_a - \varepsilon)\, P^a_{2i}\, h_i\,
          G_{i1}(\varepsilon_a - \omega)
        + G_{12}(\varepsilon_a + \omega - \varepsilon)\,
          G_{2i}(\varepsilon_a + \omega)\, h_i\, P^a_{i1}
      \Big],
    \f]
    with \f$ P^a = \ket{a}\bra{a} \f$ for the polarised core states a, and
    \f$ G = G^{\rm e} + G^{\rm c} \f$ (see green_loop). Expanding each G
    gives the terms with one and two core lines; those with three core
    lines sum to zero. Each term has the form
    \f$ {\rm Tr}[A_{12}\,(B\,h\,C)_{21}] \f$ (see trace_product,
    perturbed_greens), with angular factor L_ang_dPi.

    If h joins a core state to itself, the second and third terms each
    contain a part with \f$ P^a h P^a \f$ (from the core part of the inner
    Green's function), which has the factor \f$ 1/\omega \f$. Only their
    sum is regular, since
    \f[
      G(\varepsilon_1) - G(\varepsilon_2) =
        (\varepsilon_2 - \varepsilon_1)\,G(\varepsilon_1)\,G(\varepsilon_2) .
    \f]
    The two are replaced by the single term
    \f[
      -\left[G(\varepsilon_a - \varepsilon)\,
        G(\varepsilon_a + \omega - \varepsilon)\right]_{12}
      \left[P^a\, h\, P^a\right]_{21} ,
    \f]
    which is the same for any \f$ \omega \f$, and finite for a static field
    (e.g., diagonal hyperfine matrix elements).

    @param h        External field operator
    @param epsilon  Complex energy
    @param omega    Frequency of the external field
    @param dV_kappa Core polarisation matrices, from dV_matrices (all
                    empty for none)

    @return \f$ \delta\pi^{(b)}_{kl}(r_1, r_2) \f$, indexed [k][l]: k, l are
    the multipolarities at vertices 1 and 2. No integration measures.
 
    @note Hole-particle (FeynmanOptions): the excited part of every loop
          Green's function includes the interaction with the hole in a,
          the core state of that term (see green_loop). In the terms where
          h joins two core states the hole changes along the loop, so this
          is approximate there.
  */
  LinAlg::Matrix<ComplexRMatrix>
  dPi_b(const DiracOperator::TensorOperator &h, std::complex<double> epsilon,
        double omega,
        const LinAlg::Matrix<std::optional<ComplexGMatrix>> &dV_kappa) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Product gA * <A||h||B> * gB, for the external field h: the
    perturbed Green's function, for gA and gB Green's functions.

    @details
    t is the radial part of the operator @p h between partial waves A and
    B, a 2x2 matrix in the spinor components at each r (see
    DiracOperator::TensorOperator::radial_rhs):
    \f[
      [g_A\, t\, g_B]^{\mu\nu}_{12} = \sum_i \sum_{\lambda\sigma}
        (g_A)^{\mu\lambda}_{1i}\, t^{\lambda\sigma}(r_i)\, dr_i\,
        (g_B)^{\sigma\nu}_{i2} .
    \f]
    The result is this times the angular factor of h
    (TensorOperator::angularF). With core polarisation, the non-local
    \f$ g_A\, \delta V_{AB}\, g_B \f$ is added (see dV_matrices).

    @param h        External field operator
    @param kappa_A  Kappa of the partial wave of @p gA
    @param kappa_B  Kappa of the partial wave of @p gB
    @param gA       Matrix to the left of the operator
    @param gB       Matrix to the right of the operator
    @param dV_AB    Core polarisation matrix between partial waves A and B
                    (from dV_matrices); empty for none

    @warning The operator must be local in r (multiplicative): not valid
             for operators that involve derivatives.

    @note Without g components, only the ff part of the operator is
          included: zero for operators that only mix f and g.
  */
  ComplexGMatrix
  perturbed_greens(const DiracOperator::TensorOperator &h, int kappa_A,
                   int kappa_B, const ComplexGMatrix &gA,
                   const ComplexGMatrix &gB,
                   const std::optional<ComplexGMatrix> &dV_AB = {}) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Core polarisation correction to the external field, dV, as a
    matrix on the sub-grid for each pair of partial waves.

    @details
    dV is non-local, so is expanded over a complete set of states n of
    partial wave B (as for the exchange potential):
    \f[
      \delta V_{AB}(r_1, r_2) = \sum_n [\delta V F_n]_A(r_1)\, F_n^\dagger(r_2) ,
    \f]
    including both integration measures. The reduced matrix element is
    \f$ \langle a || \delta V || b \rangle = F_a \cdot (\delta V_{AB} F_b) \f$.
    Direct and exchange parts of dV are treated together.

    @param h    External field operator (for the selection rules)
    @param dV   Core polarisation, already solved. nullptr: all entries
                empty
    @param conj false: dV for a field of positive frequency (absorption);
                true: its conjugate, for a field of negative frequency (see
                ExternalField::CorePolarisation::dV). Makes no difference
                if dV was solved at zero frequency

    @return Matrices indexed by kappa index, (A, B); empty where h is zero
            by selection rules.

    @note The complete set of states is formed on the first call, and kept.
  */
  LinAlg::Matrix<std::optional<ComplexGMatrix>>
  dV_matrices(const DiracOperator::TensorOperator &h,
              const ExternalField::CorePolarisation *dV, bool conj = false);

  //----------------------------------------------------------------------------
  /*!
    @brief Screening correction to the Coulomb line of multipolarity k, at
    integration-grid energy index ie; zero if screening is off.

    @details
    \f[
      \delta q^k(\varepsilon) = - i\,[q\,\pi\,\tilde q]^k(\varepsilon) ,
    \f]
    so the screened line is
    \f[
      \tilde q^k(\varepsilon) = q^k + \delta q^k(\varepsilon) .
    \f]
    From the stored Q*Pi*Q (must already be formed if screening is on).
    Integration measure as for Feynman::get_qk.
  */
  ComplexRMatrix dQ_screening(int k, std::size_t ie) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Green's function for the lines of the polarisation loop,
    G = G^e + G^c: excited states and polarised core states.

    @details
    \f[
      G(\varepsilon) = G^{\rm e}(\varepsilon)
        + \sum_a \frac{\ket{a}\bra{a}}{\varepsilon - \varepsilon_a},
    \f]
    where the sum is over the polarised core states (n >= n_min) only. If
    all core states are polarised, this is the full Hartree-Fock Green's
    function.

    @param kappa  Partial wave
    @param en     Complex energy
    @param Fc_hp  If given, the excited part includes the hole-particle
                  interaction for this core state (see
                  Feynman::green_excited)
  */
  ComplexGMatrix green_loop(int kappa, std::complex<double> en,
                            const DiracSpinor *Fc_hp = nullptr) const;

  //----------------------------------------------------------------------------
  /*!
    @brief Core part of the loop Green's function (see green_loop).

    @details
    \f[
      G^{\rm c}(\varepsilon) =
        \sum_a \frac{\ket{a}\bra{a}}{\varepsilon - \varepsilon_a},
    \f]
    where the sum is over the polarised core states (n >= n_min) of this
    partial wave.

    @param kappa    Partial wave
    @param en       Complex energy
    @param Fc_skip  If given, this core state is left out of the sum
  */
  ComplexGMatrix green_loop_core(int kappa, std::complex<double> en,
                                 const DiracSpinor *Fc_skip = nullptr) const;

  //----------------------------------------------------------------------------
  // E1 operator only:

  /*!
    @brief Direct structure radiation diagram with the external field on the
    outer Green's function line; E1 only.

    @details
    The outer Green's function of the direct diagram is replaced by the perturbed
    Green's function (see perturbed_greens),
    \f[
      \delta G 
      = G_A(\varepsilon + \omega + \Omega)\, h\, G_B(\varepsilon + \omega),
    \f]
    summed over the partial waves A, B of the two Green's functions, and
    the multipolarity k of Q*Pi*Q. The matrix element is
    `Fv * (Sigma * Fw)`: v is the state on the left, w on the right.

    @param kv     Kappa of the state on the left
    @param env    Energy of the state on the right
    @param kw     Kappa of the state on the right
    @param Omega  Frequency of the external field

    @note Forms Q*Pi*Q if it has not yet been formed (see Feynman::form_qpiq).
  */
  GMatrix Sigma_SR_direct(int kv, double env, int kw, double Omega);

  /*!
    @brief Direct structure radiation diagrams with the external field
    inside the polarisation loop; E1 only.

    @details
    The polarisation operator is replaced by the perturbed polarisation
    operator (see perturbed_polarisation_kl), between two bare Coulomb
    lines; the outer line is the unperturbed Green's function. Sum of the
    two diagrams (field on either line of the loop): with
    \f$ \delta\Pi \f$ that of perturbed_polarisation_kl, the second is
    \f[
      (-1)^{k+l+K}\, \delta\Pi^{lk}(r_2, r_1; \Omega - \omega) .
    \f]
    Same result as dSigma_bc for the E1 operator.

    Parameters as for Sigma_SR_direct.
  */
  GMatrix Sigma_SR_direct2(int kv, double env, int kw, double Omega);

  /*!
    @brief Perturbed Green's function, dG = gA * h * gB; E1 only.

    @details
    \f[
      \delta G_{ij} = \sum_k (g_A)_{ik}\, r_k\, dr_k\, (g_B)_{kj},
    \f]
    with spinor indices contracted. Radial part only (h = r): the angular
    factor of h is not included.
  */
  ComplexGMatrix perturbed_greens(ComplexGMatrix gA, ComplexGMatrix gB) const;

  /*!
    @brief Perturbed polarisation operator, for multipolarities k and l of
    the two Coulomb lines; E1 only.

    @details
    The polarisation loop with the external field acting on one of its
    Green's function lines. Summed over the polarised core states and the
    partial waves of the two internal lines; includes the angular factor
    (see L_ang_dSR_QdPiQ). The lines of the loop are the loop Green's
    functions (see green_loop). Same as dPi_b for the E1 operator.

    @param k               Multipolarity of the first Coulomb line
    @param l               Multipolarity of the second Coulomb line
    @param omega           Complex frequency
    @param Omega           Frequency of the external field
    @param hole_particle   Include the hole-particle interaction
  */
  ComplexRMatrix perturbed_polarisation_kl(int k, int l,
                                           std::complex<double> omega,
                                           double Omega,
                                           bool hole_particle) const;

private:
  // Complete set of states used to write dV as a matrix (see dV_matrices);
  // formed when first needed
  std::vector<DiracSpinor> m_coordinate_basis{};
};

//==============================================================================
//==============================================================================

//------------------------------------------------------------------------------
/*!
  @brief Angular factor for structure radiation diagram (a)
  (FeynmanSR::dSigma_a).

  @details
  \f[
    (-1)^{j_A + j_B + k + K}
    \begin{Bmatrix} j_w & K & j_v \\ j_B & k & j_A \end{Bmatrix},
  \f]
  with K the rank of the external field, w the final state (on the left),
  v the initial state (on the right), A and B the partial waves of the two
  Green's functions (A next to w), and k the multipolarity of the Coulomb
  lines.
*/
double L_ang_dSR(int K, int kw, int kv, int kB, int kA, int k);

//----------------------------------------------------------------------------
/*!
  @brief Angular factor for a term of the perturbed polarisation operator
  (FeynmanSR::dPi_b).

  @details
  For a term \f$ {\rm Tr}[A_{12}\,(B\,h\,C)_{21}] \f$, with multipolarity
  k at vertex 1 and l at vertex 2, and K the rank of h:
  \f[
    (-1)^{k+l+K+1}
    \begin{Bmatrix} k & l & K \\ j_B & j_C & j_A \end{Bmatrix}
    C^k_{CA}\, C^l_{AB} .
  \f]
  The angular factor of h, \f$ \langle B\|h\|C\rangle \f$, is not included.
*/
double L_ang_dPi(int K, int k, int l, int kappa_A, int kappa_B, int kappa_C);

/*!
  @brief Adds one term of the perturbed polarisation operator
  (FeynmanSR::dPi_b) to each multipolarity.

  @details
  \f[
    \delta\pi_{kl}(r_1, r_2) \mathrel{+}=
      L^{kl}_{ABC}\,{\rm Tr}[A_{12}\,(B\,h\,C)_{21}] ,
  \f]
  with the angular factor L_ang_dPi, and the trace from trace_product.
  Added in place: only the few k, l allowed by the angular factor are
  touched.

  @param dPi_kl   Perturbed polarisation operator, indexed [k][l]; the term
                  is added to it
  @param K        Rank of the external field
  @param kappa_A  Partial wave of A
  @param kappa_B  Partial wave of B
  @param kappa_C  Partial wave of C
  @param A        Line running from vertex 1 to 2
  @param BhC      Perturbed line, B*<B||h||C>*C (see
                  FeynmanSR::perturbed_greens)
*/
void add_dPi_term(LinAlg::Matrix<ComplexRMatrix> *dPi_kl, int K, int kappa_A,
                  int kappa_B, int kappa_C, const ComplexGMatrix &A,
                  const ComplexGMatrix &BhC);

/*!
  @brief Angular factor coupling the perturbed polarisation operator to the
  outer Green's function line, for structure radiation diagrams (b) + (c)
  (FeynmanSR::dSigma_bc).

  @details
  \f[
    (-1)^{j_w + j_v}
    \begin{Bmatrix} j_w & K & j_v \\ l & j_\alpha & k \end{Bmatrix},
  \f]
  with w the final state (on the left), v the initial state (on the right),
  \f$ \alpha \f$ the partial wave of the outer Green's function, and k, l the
  multipolarities of the left and right Coulomb lines. The factors
  \f$ C^k_{w\alpha}\, C^l_{\alpha v} \f$ are not included.
*/
double L_ang_dSigma_bc(int K, int kappa_w, int kappa_v, int kappa_alpha, int k,
                       int l);

//------------------------------------------------------------------------------
/*!
  @brief Angular factor for the exchange part of structure radiation
  diagram (a_x) (FeynmanSR::dSigma_ax_cx).

  @details
  \f[
    (-1)^{j_\gamma + j_\delta + k + l}
    \begin{Bmatrix} j_\beta & l & j_\gamma \\ j_\delta & k & j_v \end{Bmatrix}
    C^l_{\beta\gamma}\, C^k_{\gamma\delta}\, C^l_{\delta v} ,
  \f]
  with v the initial state (on the right), \f$ \beta \f$ the partial wave
  of the Green's function next to the exchange operator,
  \f$ \gamma, \delta \f$ those of the two lines of the exchange operator
  (\f$ \delta \f$ next to v), l the multipolarity of the Coulomb line
  joining the \f$ (\beta, \gamma) \f$ and \f$ (\delta, v) \f$ vertices, and
  k that of the line joining the \f$ (\gamma, \delta) \f$ vertex to the
  final state.
*/
double L_ang_dSigma_ax(int kappa_v, int kappa_beta, int kappa_gamma,
                       int kappa_delta, int k, int l);

//------------------------------------------------------------------------------
/*!
  @brief Angular factor for the exchange structure radiation diagram (b_x)
  (FeynmanSR::dSigma_bx).

  @details
  \f[
    (-1)^{j_w + j_v + k + l + K}
    \begin{Bmatrix}
      j_\beta & k & j_w \\
      j_{\gamma_i} & j_{\gamma_j} & K \\
      l & j_\delta & j_v
    \end{Bmatrix}
    C^k_{w\beta}\, C^l_{\beta\gamma_i}\, C^k_{\gamma_j\delta}\,
    C^l_{\delta v} ,
  \f]
  with K the rank of the external field, w the final state (on the left),
  v the initial state (on the right), \f$ \beta \f$ the partial wave of
  the Green's function next to w, \f$ \gamma_i, \gamma_j \f$ those of
  the perturbed line either side of the field vertex (\f$ \gamma_i \f$
  next to \f$ \beta \f$), \f$ \delta \f$ that of the line next to v, k
  the multipolarity of the Coulomb line joining the \f$ (w, \beta) \f$
  and \f$ (\gamma_j, \delta) \f$ vertices, and l that of the line
  joining the \f$ (\beta, \gamma_i) \f$ and \f$ (\delta, v) \f$
  vertices.

  The angular factor of h,
  \f$ \langle \gamma_i\|h\|\gamma_j\rangle \f$, is not included.
*/
double L_ang_dSigma_bx(int K, int kappa_w, int kappa_v, int kappa_beta,
                       int kappa_gamma_i, int kappa_gamma_j, int kappa_delta,
                       int k, int l);

//------------------------------------------------------------------------------
/*!
  @brief Spinor trace of the product of two spinor matrices with opposite
  coordinate order, Tr[A_12 X_21].

  @details
  \f[
    T_{12} = \sum_{\mu\nu} A^{\mu\nu}_{12}\, X^{\nu\mu}_{21} .
  \f]
  No coordinate integration. With \f$ A = P^a = \ket{a}\bra{a} \f$, this is
  \f$ F_a^\dagger(r_2)\, X(r_2, r_1)\, F_a(r_1) \f$.
*/
ComplexRMatrix trace_product(const ComplexGMatrix &A, const ComplexGMatrix &X);

//------------------------------------------------------------------------------
/*!
  @brief Angular factor for the perturbed polarisation operator
  (FeynmanSR::perturbed_polarisation_kl); E1 only.

  @details
  (-1)^(l + k + K + 1) C^k_{ca} C^l_{ab} {K l k; j_a j_c j_b}, with K the rank
  of the external field, k and l the multipolarities of the two Coulomb
  lines, and a, b, c the three lines of the loop (the field acts between b
  and c).
*/
double L_ang_dSR_QdPiQ(int K, int k, int l, int ka, int kb, int kc);

} // namespace MBPT
