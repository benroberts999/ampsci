#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <cstddef>
#include <vector>

namespace MBPT {

/*!
  @brief Quantum numbers and energy of an orbital: all that the ladder
         iteration needs from a DiracSpinor.
*/
struct LadderOrbital {
  //! Index of the orbital in the Coulomb integral table, see
  //! DiracSpinor::nk_index()
  Coulomb::nkIndex nk_index;
  //! Index of the orbital in the ladder integral table. Normally the same as
  //! nk_index. Differs for an orbital that stands in for an external line: the
  //! same radial orbital, but with the external energy in place of its own,
  //! so its ladder integrals are different quantities and need their own label
  Coulomb::nkIndex label;
  //! Dirac quantum number
  int kappa;
  //! 2*j
  int twoj;
  //! Orbital angular momentum
  int l;
  //! Orbital energy
  double en;
};

//! Extracts the quantum numbers and energies of a set of orbitals
std::vector<LadderOrbital>
ladder_orbitals(const std::vector<DiracSpinor> &orbitals);

//! Largest 2*j of a set of orbitals (zero if empty)
int max_twoj(const std::vector<LadderOrbital> &orbitals);

//! Parity of a pair of orbitals: (l_x + l_y) mod 2
inline int pair_parity(const LadderOrbital &x, const LadderOrbital &y) {
  return (x.l + y.l) % 2;
}

//! The sign (-1)^(rank + j_x + j_y), for integer rank
inline int pair_sign(int rank, const LadderOrbital &x, const LadderOrbital &y) {
  return Angular::neg1pow_2(2 * rank + x.twoj + y.twoj);
}

//==============================================================================
/*!
  @brief Dense table of 6j symbols with four half-integer and two integer
         arguments.

  @details
  Holds
  \f[
    \sixj{j_a}{j_b}{c}{j_d}{j_e}{f}
  \f]
  for all half-integer \f$j \le\f$ @p max_tj/2 and all integer
  \f$c, f \le\f$ @p max_rank, in a flat array: look-up is a direct array read,
  with no hashing. The symbols are calculated directly, so the table does not
  depend on the range of any other 6j table.
*/
class SixJArray {
public:
  /*!
    @param max_tj    Largest 2*j of any half-integer argument
    @param max_rank  Largest value of the integer arguments
  */
  SixJArray(int max_tj, int max_rank);

  /*!
    @brief Returns the 6j symbol {ja jb c; jd je f}.
    @param tja,tjb,tjd,tje  2*j of the half-integer arguments
    @param c,f              The integer arguments
    @warning No range checking: arguments must not exceed the limits given to
             the constructor.
  */
  double get(int tja, int tjb, int c, int tjd, int tje, int f) const {
    return m_data[index(std::size_t(tja) / 2, std::size_t(tjb) / 2,
                        std::size_t(c), std::size_t(tjd) / 2,
                        std::size_t(tje) / 2, std::size_t(f))];
  }

private:
  std::size_t m_num_j;
  std::size_t m_num_rank;
  std::vector<double> m_data;

  std::size_t index(std::size_t ia, std::size_t ib, std::size_t c,
                    std::size_t id, std::size_t ie, std::size_t f) const {
    return ((((ia * m_num_j + ib) * m_num_rank + c) * m_num_j + id) * m_num_j +
            ie) *
             m_num_rank +
           f;
  }
};

//==============================================================================
/*!
  @brief Index (numbering) of the ordered orbital pairs (x, y) that can couple
         to each angular momentum (rank) and parity.

  @details
  The pairs are sorted into blocks of fixed rank and parity. The pair (x, y),
  with x from the first set and y from the second, belongs to each block
  (rank, parity) for which the triangle rule
  \f$ |j_x - j_y| \le {\rm rank} \le j_x + j_y \f$ holds and
  parity \f$= (l_x + l_y) \bmod 2\f$. Within a block the pairs are numbered
  consecutively from zero; this number is the row (or column) of the pair in
  the matrices of that block.

  The rank is the angular momentum the two orbitals are coupled to: the total
  angular momentum J of a two-particle (or two-hole) pair, or the rank K of a
  particle-hole pair.
*/
class OrbitalPairIndex {
public:
  /*!
    @param first     Orbitals for the first member of each pair
    @param second    Orbitals for the second member of each pair
    @param max_rank  Blocks with rank above this are not formed
    @param allowed   Optional, size first.size()*second.size(): pair (ix, iy)
                     is included only if allowed[ix*second.size() + iy] is
                     non-zero. Empty: all pairs are included.
  */
  OrbitalPairIndex(const std::vector<LadderOrbital> &first,
                   const std::vector<LadderOrbital> &second, int max_rank,
                   const std::vector<char> &allowed = {});

  //! Largest rank for which blocks are formed
  int max_rank() const { return m_max_rank; }

  //! Number of blocks: 2*(max_rank + 1)
  std::size_t num_blocks() const { return m_count.size(); }

  //! Number of pairs in the block (rank, parity)
  std::size_t size(int rank, int parity) const {
    return m_count[block(rank, parity)];
  }

  /*!
    @brief Position of the pair (ix, iy) within its block of the given rank.
    @return -1 if the pair is not a member of a block of this rank.
  */
  int index(int rank, std::size_t ix, std::size_t iy) const {
    return m_index[std::size_t(rank)][ix * m_size_second + iy];
  }

  //! Single index for the block (rank, parity): 2*rank + parity
  static std::size_t block(int rank, int parity) {
    return std::size_t(2 * rank + parity);
  }

private:
  std::size_t m_size_second;
  int m_max_rank;
  // m_index[rank][ix*m_size_second + iy]
  std::vector<std::vector<int>> m_index;
  // m_count[block(rank, parity)]
  std::vector<std::size_t> m_count;
};

//==============================================================================
//! Size of the arrays used to hold all multipoles of one integral, indexed
//! directly by multipolarity k. Larger than any k of a realistic basis.
constexpr std::size_t ladder_max_multipoles = 32;

/*!
  @brief All multipoles of one integral: value[k] for k_min <= k <= k_max.
  @details Elements outside [k_min, k_max] are not set. Empty range if
           k_min > k_max.
*/
struct Multipoles {
  int k_min;
  int k_max;
  std::array<double, ladder_max_multipoles> value;
};

/*!
  @brief Coulomb integrals Q^k_abcd for all k, from the table.
  @details k runs over the Coulomb range [see Coulomb::k_minmax_Q()], in steps
           of two; the elements in between are set to zero.
*/
Multipoles coulomb_multipoles(const LadderOrbital &a, const LadderOrbital &b,
                              const LadderOrbital &c, const LadderOrbital &d,
                              const Coulomb::QkTable &qk);

/*!
  @brief Rung integrals T^k_abcd = Q^k_abcd + L^k_abcd for all k.
  @details k runs over the ladder range [see k_minmax_L()], both parities. Q
           contributes only at the Coulomb-allowed k. The ladder integral is
           looked up (by orbital label) in @p lk and in @p lk_other: an
           integral is held by at most one of them. Either may be nullptr; if
           both are, T = Q.
*/
Multipoles rung_multipoles(const LadderOrbital &a, const LadderOrbital &b,
                           const LadderOrbital &c, const LadderOrbital &d,
                           const Coulomb::QkTable &qk,
                           const Coulomb::LkTable *lk,
                           const Coulomb::LkTable *lk_other = nullptr);

//==============================================================================
/*!
  @brief Which ladder integrals L^k_mnij to calculate: the orbitals allowed in
         each of the four slots.

  @details
  The integrals are those with \f$m\f$ in @p m_orbitals, \f$n\f$ in the
  excited set, \f$i\f$ any of the @p ij_orbitals marked in @p in_i_slot,
  and \f$j\f$ any of those marked in @p in_j_slot.

  The usual case is @p m_orbitals = the excited set. Then
  @p symmetric_mn must be true: the partners \f$L^k_{nmji}\f$ (equal by
  symmetry) are part of the set, and only \f$n \ge m\f$ is calculated.
*/
struct LadderIntegralSet {
  //! Orbitals for the first index, m
  std::vector<LadderOrbital> m_orbitals;
  //! True if m_orbitals is the excited set itself
  bool symmetric_mn;
  //! Orbitals for the third and fourth indices, i and j
  std::vector<LadderOrbital> ij_orbitals;
  //! Non-zero: the orbital is allowed as the third index, i
  std::vector<char> in_i_slot;
  //! Non-zero: the orbital is allowed as the fourth index, j
  std::vector<char> in_j_slot;

  //! Is L^k_mnij, with i = ij_orbitals[ii] and j = ij_orbitals[ij], in the
  //! set? (Either directly, or as the partner of an integral in the set)
  bool required(std::size_t ii, std::size_t ij) const {
    return (in_i_slot[ii] && in_j_slot[ij]) ||
           (symmetric_mn && in_j_slot[ii] && in_i_slot[ij]);
  }
};

//! Which of the four ladder diagrams to include
struct LadderDiagrams {
  //! Particle-particle ladder, L1
  bool L1{true};
  //! Particle-hole ladder, L2
  bool L2{true};
  //! Particle-hole ladder, L3 [L3^k_mnij = L2^k_nmji]
  bool L3{true};
  //! Hole-hole ladder, L4
  bool L4{true};
};

//! Largest value, and largest change, of the integrals written in one pass
struct LadderChange {
  double max_value{0.0};
  double max_change{0.0};
};

/*!
  @brief Adds the entries for the given set of ladder integrals to the table,
         with value zero.
  @details Entries already present are left unchanged. Each integral is
           stored once (table symmetry \f$L^k_{mnij} = L^k_{nmji}\f$), under
           the orbital labels.
  @param lk         Ladder table (entries are added in place)
  @param excited    Excited orbitals (second index, n)
  @param integrals  The integrals required
  @param max_k      Maximum multipolarity k of the entries
*/
void add_Lk_entries(Coulomb::LkTable *lk,
                    const std::vector<LadderOrbital> &excited,
                    const LadderIntegralSet &integrals, int max_k);

/*!
  @brief Calculates a set of ladder integrals L^k_mnij at once, and writes
         them into the table.

  @details
  General form of calculate_Lk_mnib() (see there for the method): the first
  index need not be an excited orbital, and the third and fourth indices are
  any required pairs. This allows an external line in either slot: an orbital
  that carries the external energy in place of its own, and its own label in
  the ladder table.

  The rung \f$T = Q + L\f$ takes \f$L\f$ from @p lk_prev and @p lk_fixed
  (looked up by label; each integral is held by at most one of the two).

  @param lk          Ladder table to write to (only existing entries)
  @param qk          Coulomb integral table
  @param excited     Excited orbitals: second index n, and intermediate states
  @param core        Core orbitals: intermediate states
  @param integrals   The integrals to calculate
  @param diagrams    Which of the diagrams L1..L4 to include
  @param lk_prev     Ladder table for the rung (may be nullptr). Must not be
                     @p lk
  @param lk_fixed    Second ladder table for the rung (may be nullptr)
  @param max_k       Largest multipolarity k of the integrals to write
  @param only_new    If true, only entries whose current value is exactly
                     zero are written
  @param a_damp      Damping factor, as in calculate_Lk_mnib()
  @return Largest absolute value and change of the integrals written
*/
LadderChange calculate_Lk(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                          const std::vector<LadderOrbital> &excited,
                          const std::vector<LadderOrbital> &core,
                          const LadderIntegralSet &integrals,
                          const LadderDiagrams &diagrams,
                          const Coulomb::LkTable *lk_prev,
                          const Coulomb::LkTable *lk_fixed, int max_k,
                          bool only_new, double a_damp);

/*!
  @brief Particle-hole ladder in the cross-coupled scheme: the matrices
         Y^K_{pj,qi} for every rank K and parity.

  @details
  \f[
    Y^K_{pj,qi} = (-1)^K \sum_{rc} B^K_{pj,rc} A^K_{rc,qi},
  \f]
  with \f$A\f$, \f$B\f$ as given in calculate_Lk_mnib(); \f$r\f$ excited,
  \f$c\f$ core. The rows are the (particle, hole) pairs \f$(p,j)\f$ of
  @p row_pairs, the columns the pairs \f$(q,i)\f$ of @p column_pairs.
  Then \f$L2^k_{mnij}\f$ follows from \f$Y^K_{mj,ni}\f$, and
  \f$L3^k_{mnij}\f$ from \f$Y^K_{ni,mj}\f$.

  @param row_particles     Particles p of the row pairs
  @param row_holes         Orbitals j of the row pairs
  @param row_pairs         Index of the pairs (p, j)
  @param column_particles  Particles q of the column pairs
  @param column_holes      Orbitals i of the column pairs
  @param column_pairs      Index of the pairs (q, i)
  @param excited           Excited orbitals (intermediate r)
  @param core              Core orbitals (intermediate c)
  @param qk                Coulomb integral table
  @param lk_prev,lk_fixed  Ladder tables for the rung, see rung_multipoles()
  @param sixj              6j symbols
  @return One matrix for each block, indexed by OrbitalPairIndex::block();
          empty if the block has no rows, columns or intermediate pairs
*/
std::vector<LinAlg::Matrix<double>>
ladder_particle_hole(const std::vector<LadderOrbital> &row_particles,
                     const std::vector<LadderOrbital> &row_holes,
                     const OrbitalPairIndex &row_pairs,
                     const std::vector<LadderOrbital> &column_particles,
                     const std::vector<LadderOrbital> &column_holes,
                     const OrbitalPairIndex &column_pairs,
                     const std::vector<LadderOrbital> &excited,
                     const std::vector<LadderOrbital> &core,
                     const Coulomb::QkTable &qk,
                     const Coulomb::LkTable *lk_prev,
                     const Coulomb::LkTable *lk_fixed, const SixJArray &sixj);

/*!
  @brief The usual set of ladder integrals, L^k_mnib: m,n excited; i in
         i_orbs; b in core (and the partners L^k_nmbi).
  @details Orbitals of @p i_orbs that are neither in @p core nor in
           @p excited are ignored.
*/
LadderIntegralSet ladder_integrals_mnib(const std::vector<DiracSpinor> &excited,
                                        const std::vector<DiracSpinor> &core,
                                        const std::vector<DiracSpinor> &i_orbs);

//==============================================================================
/*!
  @brief Adds the entries for all required ladder integrals L^k_mnib to the
         table, with value zero.

  @details
  The required integrals are those with \f$m, n\f$ in @p excited,
  \f$i\f$ in @p i_orbs (and in @p excited or @p core), \f$b\f$ in @p core, and
  \f$k\f$ allowed by k_minmax_L(), up to @p max_k. Each is stored once (table
  symmetry \f$L^k_{mnib} = L^k_{nmbi}\f$). Entries already present are left
  unchanged.

  @param lk       Ladder table (entries are added in place)
  @param excited  Excited orbitals
  @param core     Core orbitals
  @param i_orbs   Orbitals for the \f$i\f$ index
  @param max_k    Maximum multipolarity k of the entries
*/
void add_Lk_mnib_entries(Coulomb::LkTable *lk,
                         const std::vector<DiracSpinor> &excited,
                         const std::vector<DiracSpinor> &core,
                         const std::vector<DiracSpinor> &i_orbs, int max_k);

/*!
  @brief Calculates the ladder integrals L^k_mnib for all excited m,n at once,
         and writes them into the table.

  @details
  Calculates
  \f[
    L^k_{mnib} = L1^k_{mnib} + L2^k_{mnib} + L3^k_{mnib} \, [+ L4^k_{mnib}]
  \f]
  [see L1(), L2(), L3(), L4()] for all \f$m, n\f$ in @p excited, \f$i\f$ in
  @p i_orbs, \f$b\f$ in @p core and \f$k \le\f$ @p max_k, using the rung
  \f$T = Q + L\f$ with \f$L\f$ from @p lk_prev. The result is the same as
  evaluating Lkmnij() for each integral, but each ladder diagram is evaluated
  as a matrix product, in the angular-momentum coupling scheme in which it is
  diagonal.

  Particle-particle (L1) and hole-hole (L4) diagrams: the orbital pairs
  \f$(m,n)\f$ and \f$(i,b)\f$ are each coupled to total angular momentum
  \f$J\f$. For any integral \f$X^k\f$ (\f$Q\f$, \f$L\f$ or \f$T\f$),
  \f[
    X^J_{ab,cd} = (-1)^{J+j_a+j_c}\sum_k (-1)^k
                  \sixj{j_a}{j_b}{J}{j_d}{j_c}{k} X^k_{abcd},
  \f]
  \f[
    X^k_{abcd} = (-1)^k [k]\sum_J [J] (-1)^{J+j_a+j_c}
                 \sixj{j_a}{j_b}{J}{j_d}{j_c}{k} X^J_{ab,cd},
  \f]
  and the diagrams are
  \f[
    L1^J_{mn,ib} = \sum_{rs} Q^J_{mn,rs}
                   \frac{T^J_{rs,ib}}{\epsilon_{ib} - \epsilon_{rs}},
  \f]
  \f[
    L4^J_{mn,ib} = \sum_{cd} \frac{T^J_{mn,cd}}{\epsilon_{cd} - \epsilon_{mn}}
                   Q^J_{cd,ib},
  \f]
  where \f$r,s\f$ run over all ordered pairs of excited orbitals, and
  \f$c,d\f$ over all ordered pairs of core orbitals.

  Particle-hole diagrams (L2, L3): the pairs \f$(m,b)\f$ and \f$(n,i)\f$ are
  each coupled to rank \f$K\f$. With
  \f[
    A^K_{rc,ni} = \sum_u \sixj{j_c}{j_r}{K}{j_n}{j_i}{u} Q^u_{cnir},
  \f]
  \f[
    B^K_{mb,rc} = [K]\sum_l \sixj{j_m}{j_b}{K}{j_r}{j_c}{l}
                  \frac{T^l_{mrcb}}{\epsilon_{cb} - \epsilon_{mr}},
  \f]
  the diagram is one matrix product for each \f$K\f$,
  \f[
    Y^K_{mb,ni} = (-1)^K \sum_{rc} B^K_{mb,rc} A^K_{rc,ni},
  \f]
  \f[
    L2^k_{mnib} = [k]\sum_K \sixj{j_m}{j_b}{K}{j_n}{j_i}{k} Y^K_{mb,ni},
  \f]
  with \f$r\f$ excited and \f$c\f$ core, and
  \f$L3^k_{mnib} = L2^k_{nmbi}\f$.

  All \f$J\f$ and \f$K\f$ are kept: a maximum multipolarity applies only to
  the integrals \f$Q^k\f$, \f$L^k\f$ in the tables.

  The integrals with four excited orbitals, \f$Q^J_{mn,rs}\f$, are never
  stored: they are formed from the Qk table for one \f$m\f$ and a few
  \f$n\f$ at a time.

  @param lk          Ladder table to write to. Only entries already in the
                     table are written [see add_Lk_mnib_entries()]
  @param qk          Coulomb \f$Q^k\f$ integral table
  @param excited     Excited orbitals
  @param core        Core orbitals
  @param i_orbs      Orbitals for the \f$i\f$ index of the integrals to
                     calculate (those not in @p core must be in @p excited)
  @param include_L4  Include the hole-hole diagram L4
  @param lk_prev     Ladder table for the rung \f$(Q+L)\f$; nullptr for the
                     lowest order, \f$L(Q,Q)\f$. Must not be @p lk
  @param max_k       Largest multipolarity \f$k\f$ of the integrals to write
  @param only_new    If true, only entries whose current value is exactly
                     zero are written
  @param a_damp      Damping factor [0,1), used if @p only_new is false:
                     \f$L \to a\,L + (1-a)\,L_{\rm new}\f$, with \f$L\f$ the
                     current entry of @p lk. No damping for entries whose
                     current value is exactly zero

  @note Integrals that vanish because no term is allowed may come out as
        rounding noise (relative to the terms that cancel) instead of exactly
        zero.
*/
void calculate_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                       const Coulomb::LkTable *lk_prev, int max_k,
                       bool only_new, double a_damp);

} // namespace MBPT
