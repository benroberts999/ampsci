#pragma once
#include "Coulomb/meTable.hpp"
#include <array>
#include <vector>

namespace IO {
class InputBlock;
}
class Wavefunction;
class DiracSpinor;

namespace Module {

//! Module: third-order structure radiation for each valence transition,
//! split into the six Feynman diagrams
void SRBreakdown(const IO::InputBlock &input, const Wavefunction &wf);

/*!
  @brief Structure radiation diagrams a, c, ax, cx: reduced matrix elements
  <w||h_sr||v>.

  @details
  Feynman diagram labels. Direct (Coulomb line, polarisation loop, Coulomb
  line):
  - a: field on the outer (valence) line
  - b: field on the loop line running from the vertex of the w-side Coulomb
    line to that of the v-side line
  - c: field on the loop line running from the v-side to the w-side vertex

  Exchange (three internal lines between v and w):
  - ax: field on the internal line adjacent to the final state, w
  - bx: field on the middle internal line
  - cx: field on the internal line adjacent to the initial state, v

  After angular reduction, the Goldstone diagrams (time orderings) of all six
  combine into 12 terms. Six have two antisymmetrised Coulomb lines,
  W = Q + P, one attached to v and one to w; their four parts are diagrams
  a (Q_v Q_w), ax (P_v Q_w), cx (Q_v P_w), and c (P_v P_w). The other six
  give b and bx (see SR_b_bx).

  The energy of the final state appears as e_v + omega, so omega must be
  e_w - e_v for the physical transition.

  @param w, v     Final and initial states; must be basis states
  @param omega    e_w - e_v
  @param core     Core states (internal lines)
  @param excited  Excited states (internal lines)
  @param qk       Coulomb integrals for core + excited (YkTable or QkTable)
  @param h        Operator matrix elements for core + excited (may include
                  RPA)
  @param K        Operator rank
  @param fk       Effective Coulomb screening factors, for each k: each
                  Coulomb line of multipolarity k is scaled by sqrt(f_k)
                  (see Q_screened, P_screened). 1 for k beyond the list
  @param etak     Effective hole-particle factors, for each k: each term is
                  scaled by eta_mu, with mu the multipolarity of the two
                  Coulomb lines (both Q and P parts). 1 for k beyond the list

  @return {a, c, ax, cx}
*/
template <typename CoulombTable>
std::array<double, 4>
SR_a_c_ax_cx(const DiracSpinor &w, const DiracSpinor &v, double omega,
             const std::vector<DiracSpinor> &core,
             const std::vector<DiracSpinor> &excited, const CoulombTable &qk,
             const Coulomb::meTable<double> &h, int K,
             const std::vector<double> &fk, const std::vector<double> &etak);

/*!
  @brief Structure radiation diagrams b, bx: reduced matrix elements
  <w||h_sr||v>.

  @details
  Diagram labels as in SR_a_c_ax_cx. These come from the six terms (after
  angular reduction) with one direct (Q) and one antisymmetrised Coulomb
  line, W = Q + P; its two parts are diagrams b (Q) and bx (P).
  Parameters as in SR_a_c_ax_cx, except etak: each term is scaled by eta_mu,
  with mu the multipolarity of the direct (Q) line.

  @return {b, bx}
*/
template <typename CoulombTable>
std::array<double, 2>
SR_b_bx(const DiracSpinor &w, const DiracSpinor &v, double omega,
        const std::vector<DiracSpinor> &core,
        const std::vector<DiracSpinor> &excited, const CoulombTable &qk,
        const Coulomb::meTable<double> &h, int K, const std::vector<double> &fk,
        const std::vector<double> &etak);

//! f[k], or 1 if f has no entry for k (screening factors fk, etak)
double factor_k(const std::vector<double> &f, int k);

//! Coulomb integral Q^k_abcd, with effective screening: scaled by sqrt(f_k)
template <typename CoulombTable>
double Q_screened(const CoulombTable &qk, const std::vector<double> &fk, int k,
                  const DiracSpinor &a, const DiracSpinor &b,
                  const DiracSpinor &c, const DiracSpinor &d);

//! Exchange Coulomb integral P^k_abcd = [k] sum_l {a c k; b d l} Q^l_abdc,
//! with effective screening: each Q^l as in Q_screened
template <typename CoulombTable>
double P_screened(const CoulombTable &qk, const std::vector<double> &fk, int k,
                  const DiracSpinor &a, const DiracSpinor &b,
                  const DiracSpinor &c, const DiracSpinor &d);

} // namespace Module
