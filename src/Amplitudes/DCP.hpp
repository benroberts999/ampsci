#pragma once
#include "Coulomb/meTable.hpp"
#include <vector>
class DiracSpinor;

namespace Amplitudes {

/*!
  @brief Double core polarisation (DCP) correction to a second-order
  amplitude: centre diagrams with one hole, C1.

  @details
  For a transition \f$ v \to w \f$ due to two one-body operators,
  \f$ t \f$ (rank \f$ k_t \f$, frequency \f$ \omega_t \f$) and \f$ s \f$
  (\f$ k_s \f$, \f$ \omega_s \f$), coupled to rank \f$ K \f$ (see
  Amplitudes::sos_valence), the DCP correction \f$ \delta A^K \f$ is the
  lowest-order Coulomb correction in which neither \f$ t \f$ nor
  \f$ s \f$ acts on the valence states. It is the sum of twelve terms:
  dcp_C1() + dcp_C2() + dcp_R12() + dcp_L12().

  C/R/L: the Coulomb interaction is in the centre/right/left of the
  diagram; 1/2: one/two holes; \f$ ts \f$/\f$ st \f$: which operator acts
  first. With \f$ T_{ij} = \redmatel{i}{t}{j} \f$ and
  \f$ S_{ij} = \redmatel{i}{s}{j} \f$ (which may include core
  polarisation), \f$ a, b \f$ core states and \f$ n, m \f$ excited states,
  this function returns
  \f$ \delta A^K(C1_{ts}) + \delta A^K(C1_{st}) \f$:
  \f[
    \delta A^K(C1_{ts}) = \sum_{amn}\frac{(-1)^{j_w+j_m}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_m}{j_n}{j_a}
      \frac{T_{na}S_{am}W^{K}_{wmvn}}
           {(\en_a-\en_n+\omega_t)(\en_a-\en_m-\omega_s)} ,
  \f]
  \f[
    \delta A^K(C1_{st}) = \sum_{amn}
      \frac{(-1)^{K+k_s+k_t+j_w+j_m}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_n}{j_m}{j_a}
      \frac{T_{am}S_{na}W^{K}_{wmvn}}
           {(\en_a-\en_n+\omega_s)(\en_a-\en_m-\omega_t)} .
  \f]

  @param w        Final valence state.
  @param v        Initial valence state.
  @param omega_t  Frequency of \f$ t \f$.
  @param omega_s  Frequency of \f$ s \f$.
  @param K        Total rank.
  @param kt       Rank of \f$ t \f$.
  @param ks       Rank of \f$ s \f$.
  @param T        Table of \f$ \redmatel{i}{t}{j} \f$, for all pairs of
                  @p core and @p excited states, both orderings; see
                  Amplitudes::me_table.
  @param S        The same, for \f$ s \f$.
  @param core     Core states of the basis.
  @param excited  Excited states of the basis.

  @return \f$ \delta A^K(C1_{ts}) + \delta A^K(C1_{st}) \f$.
*/
double dcp_C1(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
              double omega_s, int K, int kt, int ks,
              const Coulomb::meTable<double> &T,
              const Coulomb::meTable<double> &S,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited);

/*!
  @brief DCP correction to a second-order amplitude: centre diagrams with
  two holes, C2.

  @details
  Returns \f$ \delta A^K(C2_{ts}) + \delta A^K(C2_{st}) \f$:
  \f[
    \delta A^K(C2_{ts}) = \sum_{abn}
      \frac{(-1)^{K+k_s+k_t+j_w-j_a}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_b}{j_a}{j_n}
      \frac{T_{na}S_{bn}W^{K}_{wavb}}
           {(\en_a-\en_n+\omega_t)(\en_b-\en_n-\omega_s)} ,
  \f]
  \f[
    \delta A^K(C2_{st}) = \sum_{abn}\frac{(-1)^{j_w-j_a}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_a}{j_b}{j_n}
      \frac{T_{bn}S_{na}W^{K}_{wavb}}
           {(\en_a-\en_n+\omega_s)(\en_b-\en_n-\omega_t)} .
  \f]
  Notation and arguments as for dcp_C1().
*/
double dcp_C2(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
              double omega_s, int K, int kt, int ks,
              const Coulomb::meTable<double> &T,
              const Coulomb::meTable<double> &S,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited);

/*!
  @brief DCP correction to a second-order amplitude: right diagrams, R1
  (one hole) and R2 (two holes).

  @details
  Returns the sum of the four terms (they share a Coulomb integral):
  \f[
    \delta A^K(R1_{ts}) = \sum_{amn}
      \frac{(-1)^{K+k_s+k_t+j_w+j_a}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_m}{j_a}{j_n}
      \frac{T_{na}S_{mn}W^{K}_{wavm}}
           {(\en_a-\en_n+\omega_t)(\en_a-\en_m+\omega_t+\omega_s)} ,
  \f]
  \f[
    \delta A^K(R1_{st}) = \sum_{amn}\frac{(-1)^{j_w+j_a}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_a}{j_m}{j_n}
      \frac{T_{mn}S_{na}W^{K}_{wavm}}
           {(\en_a-\en_n+\omega_s)(\en_a-\en_m+\omega_s+\omega_t)} ,
  \f]
  \f[
    \delta A^K(R2_{ts}) = \sum_{abn}\frac{(-1)^{j_w-j_b}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_b}{j_n}{j_a}
      \frac{T_{na}S_{ab}W^{K}_{wbvn}}
           {(\en_a-\en_n+\omega_t)(\en_b-\en_n+\omega_t+\omega_s)} ,
  \f]
  \f[
    \delta A^K(R2_{st}) = \sum_{abn}
      \frac{(-1)^{K+k_s+k_t+j_w-j_b}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_n}{j_b}{j_a}
      \frac{T_{ab}S_{na}W^{K}_{wbvn}}
           {(\en_a-\en_n+\omega_s)(\en_b-\en_n+\omega_s+\omega_t)} .
  \f]
  Notation and arguments as for dcp_C1().

  @note R2 contains matrix elements between two core states.
*/
double dcp_R12(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
               double omega_s, int K, int kt, int ks,
               const Coulomb::meTable<double> &T,
               const Coulomb::meTable<double> &S,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited);

/*!
  @brief DCP correction to a second-order amplitude: left diagrams, L1
  (one hole) and L2 (two holes).

  @details
  Returns the sum of the four terms (they share a Coulomb integral):
  \f[
    \delta A^K(L1_{ts}) = \sum_{amn}
      \frac{(-1)^{K+k_s+k_t+j_w+j_n}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_a}{j_n}{j_m}
      \frac{T_{mn}S_{am}W^{K}_{wnva}}
           {(\en_a-\en_n-\omega_t-\omega_s)(\en_a-\en_m-\omega_s)} ,
  \f]
  \f[
    \delta A^K(L1_{st}) = \sum_{amn}\frac{(-1)^{j_w+j_n}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_n}{j_a}{j_m}
      \frac{T_{am}S_{mn}W^{K}_{wnva}}
           {(\en_a-\en_n-\omega_s-\omega_t)(\en_a-\en_m-\omega_t)} ,
  \f]
  \f[
    \delta A^K(L2_{ts}) = \sum_{abn}\frac{(-1)^{j_w-j_n}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_n}{j_a}{j_b}
      \frac{T_{ab}S_{bn}W^{K}_{wnva}}
           {(\en_a-\en_n-\omega_t-\omega_s)(\en_b-\en_n-\omega_s)} ,
  \f]
  \f[
    \delta A^K(L2_{st}) = \sum_{abn}
      \frac{(-1)^{K+k_s+k_t+j_w-j_n}}{\sqrt{[K]}}
      \sixj{k_t}{k_s}{K}{j_a}{j_n}{j_b}
      \frac{T_{bn}S_{ab}W^{K}_{wnva}}
           {(\en_a-\en_n-\omega_s-\omega_t)(\en_b-\en_n-\omega_t)} .
  \f]
  Notation and arguments as for dcp_C1().

  @note L2 contains matrix elements between two core states.
*/
double dcp_L12(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
               double omega_s, int K, int kt, int ks,
               const Coulomb::meTable<double> &T,
               const Coulomb::meTable<double> &S,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited);

/*!
  @brief Total DCP correction to a second-order amplitude: sum of all
  twelve terms.

  @details
  \f[
    \delta A^K = \delta A^K(C1) + \delta A^K(C2)
               + \delta A^K(R1) + \delta A^K(R2)
               + \delta A^K(L1) + \delta A^K(L2) ,
  \f]
  each with its \f$ ts \f$ and \f$ st \f$ terms: dcp_C1() + dcp_C2() +
  dcp_R12() + dcp_L12(). Notation and arguments as for dcp_C1().

  @return \f$ \delta A^K \f$.
*/
double dcp(const DiracSpinor &w, const DiracSpinor &v, double omega_t,
           double omega_s, int K, int kt, int ks,
           const Coulomb::meTable<double> &T, const Coulomb::meTable<double> &S,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited);

} // namespace Amplitudes
