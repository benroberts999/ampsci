#pragma once
#include "MBPT/Feynman.hpp"
#include "MBPT/RadialMatrix.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include <complex>

namespace MBPT {

//------------------------------------------------------------------------------
/*!
  @brief Structure radiation diagrams via the Feynman technique.

  @details
  Extends Feynman with the correction to the correlation potential that is
  first order in an external field, h, of frequency Omega: the field acts on
  an internal line of the direct Sigma diagram. The diagrams are formed as
  matrices on the radial sub-grid, by integration along Im(w), as for
  Feynman::Sigma_direct.

  - Sigma_SR_direct  : field on the valence (outer) Green's function line
  - Sigma_SR_direct2 : field on either Green's function line of the
                       polarisation loop

  @warning The external field is the E1 operator (rank 1, radial part r):
           hard-coded.
*/
class FeynmanSR : public Feynman {

public:
  using Feynman::Feynman;

  /*!
    @brief Direct structure radiation diagram with the external field on the
    valence line.

    @details
    The valence line of the direct diagram is replaced by the perturbed
    Green's function (see perturbed_greens),
    \f[
      \delta G = G_A(\varepsilon + \omega + \Omega)\, h\, G_B(\varepsilon + \omega),
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
    inside the polarisation loop.

    @details
    The polarisation operator is replaced by the perturbed polarisation
    operator (see perturbed_polarisation_kl); the valence line is the
    unperturbed Green's function. Sum of the two diagrams (field on either
    line of the loop).

    Parameters as for Sigma_SR_direct.
  */
  GMatrix Sigma_SR_direct2(int kv, double env, int kw, double Omega);

  /*!
    @brief Perturbed Green's function, dG = gA * h * gB.

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
    the two Coulomb lines.

    @details
    The polarisation loop with the external field acting on one of its
    Green's function lines. Summed over the polarised core states and the
    partial waves of the two internal lines; includes the angular factor
    (see L_ang_dSR_QdPiQ).

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
};

/*!
  @brief Angular factor for the direct structure radiation diagram with the
  external field on the valence line (FeynmanSR::Sigma_SR_direct).

  @details
  (-1)^(j_A + j_B + k + K) {j_w K j_v; j_B k j_A}, with K = 1 the rank of
  the external field. A, B are the partial waves of the two Green's
  functions, k the multipolarity of the Coulomb lines.
*/
double L_ang_dSR(int kw, int kv, int kB, int kA, int k);

/*!
  @brief Angular factor for the perturbed polarisation operator
  (FeynmanSR::perturbed_polarisation_kl).

  @details
  (-1)^(l + k + 1) C^k_{ca} C^l_{ab} {K l k; j_a j_c j_b}, with K the rank
  of the external field, k and l the multipolarities of the two Coulomb
  lines, and a, b, c the three lines of the loop (the field acts between b
  and c).
*/
double L_ang_dSR_QdPiQ(int K, int k, int l, int ka, int kb, int kc);

/*!
  @brief Element-wise product of two spinor matrices, block by block.

  @details
  \f[
    C^{\mu\nu}_{ij} = A^{\mu\nu}_{ij} B^{\mu\nu}_{ij} .
  \f]
  No sum over the spinor indices. If @p b has no g parts, the g blocks of
  the result are zero.
*/
ComplexGMatrix mult_elements_spinor(ComplexGMatrix a, const ComplexGMatrix &b);

} // namespace MBPT
