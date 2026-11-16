#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <optional>
#include <vector>

namespace MBPT {

/*!
  @brief Full ladder integral summed over all diagrams.

  @details
  Computes
  \f[
    L^k_{mnij} = L1^k_{mnij} + L2^k_{mnij} + L3^k_{mnij} [+ L4^k_{mnij}]
  \f]
  where \f$ L3^k_{mnij} = L2^k_{nmji} \f$ and \f$ L4 \f$ involves
  core--core intermediate states. @p Lk points to the ladder table from
  the previous iteration; pass nullptr on the first iteration.

  Lkmnij() and L1(), L2(), L3(), L4() evaluate the ladder diagrams one
  integral at a time: an independent evaluation of each ladder integral
  L^k_mnij as a direct sum over the intermediate states, with no recoupling.
  They are used to test the ladder iteration [calculate_Lk_mnib()], which
  evaluates all the integrals at once as matrix products. Far too slow for a
  production calculation.

  @param k          Multipole rank
  @param m,n        Excited (particle) orbitals
  @param i,j        Core (hole) or valence orbitals
  @param qk         Coulomb \f$ Q^k \f$ integral table
  @param core       Core orbitals
  @param excited    Excited orbitals
  @param include_L4 Include the core--core diagram L4
  @param SJ         6j symbol table
  @param Lk         Ladder table from previous iteration (nullptr on first)
  @param e_i        Optional: used in place of i.en() in energy denominators
  @param e_m        Optional: used in place of m.en() in energy denominators
  @return \f$ L^k_{mnij} \f$
  @note To evaluate at a fixed external energy (for a correlation potential),
        use @p e_i or @p e_m (for external line in the i or m slot): the energy
        only enters the denominators, never the integral lookups.
*/
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ,
              const Coulomb::LkTable *const Lk = nullptr,
              std::optional<double> e_i = {}, std::optional<double> e_m = {});

/*!
  @brief Particle--particle ladder diagram L1.

  @details
  \f[
    L1^k_{mnij} = \sum_{rs,ul} A^{kul}_{mnrsij}
      \frac{Q^u_{mnrs}\,(Q+L)^l_{rsij}}{\epsilon_{ij} - \epsilon_{rs}}
  \f]
  with the angular coefficient
  \f[
    A^{kul}_{mnrsij} = (-1)^{m+n+r+s+i+j+1}\,[k]
      \sixj{m}{i}{k}{l}{u}{r}\sixj{n}{j}{k}{l}{u}{s}
  \f]
  Intermediate states \f$ r,s \f$ run over excited orbitals.

  @param k       Multipole rank
  @param m,n     Excited (particle) orbitals
  @param i,j     Core (hole) or valence orbitals
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param excited Excited orbitals
  @param SJ      6j symbol table
  @param Lk      Ladder table from previous iteration (nullptr on first)
  @param e_i     Optional: used in place of i.en() in energy denominator
  @return \f$ L1^k_{mnij} \f$
*/
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          std::optional<double> e_i = {});

/*!
  @brief Particle--hole ladder diagram L2.

  @details
  \f[
    L2^k_{mnij} = \sum_{rc,ul} (-1)^{k+u+l+1} A^{klu}_{mjcrin}
      \frac{Q^u_{cnir}\,(Q+L)^l_{mrcj}}{\epsilon_{cj} - \epsilon_{mr}}
  \f]
  Intermediate states: \f$ r \f$ runs over excited, \f$ c \f$ over core.
  The diagram \f$ L3 \f$ is the exchange partner
  \f$ L3^k_{mnij} = L2^k_{nmji} \f$.

  @param k       Multipole rank
  @param m,n     Excited (particle) orbitals
  @param i,j     Core (hole) or valence orbitals
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @param SJ      6j symbol table
  @param Lk      Ladder table from previous iteration (nullptr on first)
  @param e_j     Optional: used in place of j.en() in energy denominator
  @param e_m     Optional: used in place of m.en() in energy denominator
  @return \f$ L2^k_{mnij} \f$
*/
double L2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          std::optional<double> e_j = {}, std::optional<double> e_m = {});

/*!
  @brief Exchange partner of L2; equals L2 with m,n and i,j swapped.

  @details
  \f[ L3^k_{mnij} = L2^k_{nmji} \f]
  (@p e_i optional: used in place of i.en() in energy denominator)
*/
inline double
L3(int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
   const DiracSpinor &j, const Coulomb::QkTable &qk,
   const std::vector<DiracSpinor> &core,
   const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
   const Coulomb::LkTable *const Lk = nullptr, std::optional<double> e_i = {}) {
  return L2(k, n, m, j, i, qk, core, excited, SJ, Lk, e_i);
}

/*!
  @brief Core--core (hole--hole) ladder diagram L4.

  @details
  \f[
    L4^k_{mnij} = \sum_{cd,ul} A^{klu}_{mncdij}
      \frac{Q^u_{cdij}\,(Q+L)^l_{mncd}}{\epsilon_{cd} - \epsilon_{mn}}
  \f]
  with the angular coefficient \f$ A \f$ as in L1().
  Intermediate states run over core orbitals only, making this the
  hole--hole counterpart of the particle--particle diagram L1.
  Enable via @p include_L4 in Lkmnij().

  @param k    Multipole rank
  @param m,n  Excited (particle) orbitals
  @param i,j  Core (hole) or valence orbitals
  @param qk   Coulomb \f$ Q^k \f$ integral table
  @param core Core orbitals
  @param SJ   6j symbol table
  @param Lk   Ladder table from previous iteration (nullptr on first)
  @param e_m  Optional: used in place of m.en() in energy denominator
  @return \f$ L4^k_{mnij} \f$
*/
double L4(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr,
          std::optional<double> e_m = {});

/*!
  @brief Screening and hole-particle interaction of the ladder diagrams,
         P^k_mnij: pieces (1)-(3) of the particle-hole term of the doubles
         equation, with their partners.

  @details
  In the m-scheme, with X = g + l the dressed integral (l from @p Lk),
  \f[
    P_{mnij} = \sum_{rc}\left(
      \frac{(g_{cnrj} - g_{cnjr})\,X_{mric}}{\epsilon_{ic} - \epsilon_{mr}}
      - \frac{g_{cnrj}\,X_{mrci}}{\epsilon_{ic} - \epsilon_{mr}}\right)
      + [m \leftrightarrow n,\, i \leftrightarrow j].
  \f]
  The three terms are pieces (1), (2) and (3) [P1(), P2(), P3()]. The fourth
  piece of the particle-hole term, \f$ +g_{cnjr}X_{mrci} \f$, is in the
  ladder integral as L2 + L3, written with i and j exchanged and the sign
  reversed, \f$ -g_{cnir}X_{mrcj} \f$. Each piece is diagonal in k, and the
  partner has the same k:
  \f[
    P^k_{mnij} = (P1 + P2 + P3)^k_{mnij} + (P1 + P2 + P3)^k_{nmji}.
  \f]
  Evaluated one integral at a time, as Lkmnij(), to test the matrix form
  [P_multipole()].

  @param k       Multipole rank
  @param m,n     Excited (particle) orbitals
  @param i,j     Core (hole) or valence orbitals
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @param SJ      6j symbol table
  @param Lk      Ladder table for the dressed integral (nullptr: X = g)
  @return \f$ P^k_{mnij} \f$
*/
double Pkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited,
              const Angular::SixJTable &SJ,
              const Coulomb::LkTable *const Lk = nullptr);

/*!
  @brief Piece (1) of P [see Pkmnij()], one side: the polarisation loop
         (screening).

  @details
  The m-scheme term \f$ +g_{cnrj}X_{mric}/(\epsilon_{ic}-\epsilon_{mr}) \f$
  summed over r and c:
  \f[
    P1^k_{mnij} = \frac{(-1)^k}{[k]} \sum_{rc}
      \frac{Q^k_{cnrj}\,(Q+L)^k_{mric}}{\epsilon_{ic} - \epsilon_{mr}}.
  \f]
  Intermediate states: \f$ r \f$ runs over excited, \f$ c \f$ over core.
  Exists only at the Coulomb-allowed k of the pairs (m,i) and (n,j). The
  partner is P1^k_nmji.

  @param k       Multipole rank
  @param m,n     Excited (particle) orbitals
  @param i,j     Core (hole) or valence orbitals
  @param qk      Coulomb \f$ Q^k \f$ integral table
  @param core    Core orbitals
  @param excited Excited orbitals
  @param Lk      Ladder table for the dressed integral (nullptr: X = g)
  @return \f$ P1^k_{mnij} \f$
*/
double P1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited,
          const Coulomb::LkTable *const Lk = nullptr);

/*!
  @brief Piece (2) of P [see Pkmnij()], one side.

  @details
  The m-scheme term \f$ -g_{cnjr}X_{mric}/(\epsilon_{ic}-\epsilon_{mr}) \f$
  summed over r and c:
  \f[
    P2^k_{mnij} = (-1)^k \sum_{rc,u} \sixj{c}{r}{k}{n}{j}{u}
      \frac{Q^u_{cnjr}\,(Q+L)^k_{mric}}{\epsilon_{ic} - \epsilon_{mr}}.
  \f]
  Intermediate states: \f$ r \f$ runs over excited, \f$ c \f$ over core.
  Both parities of k, the other one through \f$ L^k_{mric} \f$ only. The
  partner is P2^k_nmji. Parameters as P1().
*/
double P2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr);

/*!
  @brief Piece (3) of P [see Pkmnij()], one side.

  @details
  The m-scheme term \f$ -g_{cnrj}X_{mrci}/(\epsilon_{ic}-\epsilon_{mr}) \f$
  summed over r and c:
  \f[
    P3^k_{mnij} = (-1)^k \sum_{rc,l} \sixj{m}{i}{k}{r}{c}{l}
      \frac{Q^k_{cnrj}\,(Q+L)^l_{mrci}}{\epsilon_{ic} - \epsilon_{mr}},
  \f]
  with l over both parities (the L part). Intermediate states: \f$ r \f$ runs
  over excited, \f$ c \f$ over core. Exists only at the Coulomb-allowed k of
  the pairs (m,i) and (n,j). The partner is P3^k_nmji. Parameters as P1().
*/
double P3(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk = nullptr);

} // namespace MBPT
