#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <utility>
#include <vector>

namespace MBPT {

//! Largest multipolarity k, pair angular momentum J and particle-hole rank K
//! kept in any ladder sum
constexpr int ladder_max_k = 31;
//! Size of the arrays indexed directly by k, J or K
constexpr std::size_t ladder_max_multipoles = ladder_max_k + 1;

/*!
  @brief Range of k of the ladder integral L^k_abcd.

  @details
  Triangle rules {a,c,k} and {b,d,k}, and l_a + l_b + l_c + l_d even. Unlike
  Q^k there is no parity rule for each pair, so k takes both parities: not safe
  to step k by two. Cut off at ladder_max_k.

  @return {k_min, k_max}; an empty range (k_min > k_max) if no k is allowed
*/
std::pair<int, int> k_minmax_L(int kappa_a, int kappa_b, int kappa_c,
                               int kappa_d);

//! k_minmax_L() from the orbitals
std::pair<int, int> k_minmax_L(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d);

//! Parity of the orbital pair (a,b): (l_a + l_b) mod 2
inline int pair_parity(const DiracSpinor &a, const DiracSpinor &b) {
  return (a.l() + b.l()) % 2;
}

//! The sign (-1)^(rank + j_a + j_c), for integer rank
inline int pair_sign(int rank, const DiracSpinor &a, const DiracSpinor &c) {
  return Angular::neg1pow_2(2 * rank + a.twoj() + c.twoj());
}

//! Range {min, max} of the rank that couples both pairs (a,b) and (c,d):
//! triangle rules {j_a, j_b, rank} and {j_c, j_d, rank}, cut at max_rank.
//! Empty (min > max) if there is none
inline std::pair<int, int>
coupled_rank_range(const DiracSpinor &a, const DiracSpinor &b,
                   const DiracSpinor &c, const DiracSpinor &d, int max_rank) {
  return {
    std::max(std::abs(a.twoj() - b.twoj()), std::abs(c.twoj() - d.twoj())) / 2,
    std::min({a.twoj() + b.twoj(), c.twoj() + d.twoj(), 2 * max_rank}) / 2};
}

/*!
  @brief Which ladder diagrams to include.

  @details
  L2 or L3 alone does not have the table symmetry L_mnij = L_nmji [L3^k_mnij
  = L2^k_nmji]: for a table with that symmetry the two go together
  [calculate_Lk_mnib(), calculate_Lk_mnxb()]. With no table symmetry they
  need not [calculate_Lk_xnab()].

  P1, P2 and P3 are pieces (1)-(3) of P, the screening and hole-particle
  interaction of the ladder diagrams [see Pkmnij()]. L2 + L3 is the fourth
  piece, with i and j exchanged and the sign reversed (the direct slot). With
  all of them, the particle-hole term differs from that of the doubles
  equation of the SD method only by a part symmetric under i <-> j: it drops
  out of l_mnij - l_mnji, and so out of every energy, at all orders. Each
  piece of P, with its partner, has the table symmetry by itself. Off by
  default; not with the external line.

  hole_particle is the hole-particle interaction as the Feynman direct
  diagram includes it, through V^(a): -y^0_cc acts on the particle of each
  pair (n,c) whose hole c is a core orbital, with the hole unchanged (the
  monopole only) [see Vhp()]. With P1 it gives P_F, the part of P[X] the
  Feynman direct diagram contains (x_F = P_F[g + x_F]). It is part of piece
  (2): not with P2. Off by default; not with the external line.
*/
struct LadderDiagrams {
  //! Particle-particle ladder
  bool L1{true};
  //! Particle-hole ladder
  bool L2{true};
  //! Particle-hole ladder, L3^k_mnij = L2^k_nmji
  bool L3{true};
  //! Hole-hole ladder
  bool L4{true};
  //! Piece (1) of P, +g_cnrj X_mric: the polarisation loop (screening)
  bool P1{false};
  //! Piece (2) of P, -g_cnjr X_mric
  bool P2{false};
  //! Piece (3) of P, -g_cnrj X_mrci
  bool P3{false};
  //! Hole-particle interaction of V^(a): -y^0_jj on n, -y^0_ii on m (i core)
  bool hole_particle{false};
};

/*!
  @brief The external line of Sigma_L(e_v): the orbitals x, at the valence
         energy en_v.

  @details
  In the ladder integrals of Sigma_L(e_v), the external line x carries en_v in
  the energy denominators, not its own energy. The dressed integrals QL that
  contain it take L from the integrals with the same external line, not from
  the ladder table. Exactly one of the two tables is set, and it says where
  the line is in L^k_mnij:
  - L_mnxb: the line is the i or j orbital;
  - L_xnab: the line is the m orbital.
*/
struct ExternalLine {
  const std::vector<DiracSpinor> &x_orbitals;
  double en_v;
  const Coulomb::LkTable *L_mnxb;
  const Coulomb::NkTable *L_xnab;

  //! True if a is the line as the i or j orbital
  bool is_ij_line(const DiracSpinor &a) const;
  //! True if a is the line as the m orbital
  bool is_m_line(const DiracSpinor &a) const;
};

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
  with no hashing.
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

    @details
    tja, tjb, tjd, tje are 2*j of the half-integer arguments; c, f the integer
    arguments. No range checking: arguments must not exceed the limits given
    to the constructor.
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
  @brief Numbers the ordered orbital pairs (x,y) within each (rank, parity)
         block: the rows and columns of the coupled matrices.

  @details
  The pair (x,y), x from the first set and y from the second, belongs to each
  block (rank, parity) with \f$ |j_x - j_y| \le {\rm rank} \le j_x + j_y \f$
  and parity \f$= (l_x + l_y) \bmod 2\f$. Within a block the pairs are
  numbered consecutively from zero.

  The rank is the angular momentum the two orbitals are coupled to: the total
  angular momentum J of a particle-particle (or hole-hole) pair, or the rank K
  of a particle-hole pair.
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
  OrbitalPairIndex(const std::vector<DiracSpinor> &first,
                   const std::vector<DiracSpinor> &second, int max_rank,
                   const std::vector<char> &allowed = {});

  //! Single index for the block (rank, parity): 2*rank + parity
  static std::size_t block(int rank, int parity) {
    return std::size_t(2 * rank + parity);
  }

  //! Number of blocks: 2*(max_rank + 1)
  std::size_t num_blocks() const { return m_count.size(); }

  //! Number of pairs in the block (rank, parity)
  std::size_t size(int rank, int parity) const {
    return m_count[block(rank, parity)];
  }

  //! Row/column of (first[ix], second[iy]) in its block of this rank; -1 if
  //! not a member
  int index(int rank, std::size_t ix, std::size_t iy) const {
    return m_index[std::size_t(rank)][ix * m_size_second + iy];
  }

  //! Largest rank for which blocks are formed
  int max_rank() const { return m_max_rank; }

private:
  std::size_t m_size_second;
  int m_max_rank;
  // m_index[rank][ix*m_size_second + iy]
  std::vector<std::vector<int>> m_index;
  // m_count[block(rank, parity)]
  std::vector<std::size_t> m_count;
};

//! One zero matrix per (rank, parity) block in which both row_pairs and
//! column_pairs have members (rows row_pairs, columns column_pairs), indexed
//! by OrbitalPairIndex::block(); the other blocks are empty (0 x 0)
std::vector<LinAlg::Matrix<double>>
block_matrices(const OrbitalPairIndex &row_pairs,
               const OrbitalPairIndex &column_pairs);

//==============================================================================
/*!
  @brief All multipoles of one integral: value[k] for k = k_min, k_min +
         k_step, ..., k_max.

  @details
  Empty if k_min > k_max. Elements off the step are not used.
*/
struct Multipoles {
  int k_min = 0;
  int k_max = -1;
  int k_step = 1;
  std::array<double, ladder_max_multipoles> value{};
};

//! Q^k_abcd for all Coulomb-allowed k (k_step = 2), one table key, cut at
//! ladder_max_k
Multipoles coulomb_multipoles(const DiracSpinor &a, const DiracSpinor &b,
                              const DiracSpinor &c, const DiracSpinor &d,
                              const Coulomb::QkTable &qk);

//! (Q+L)^k_abcd for all k of k_minmax_L() (k_step = 1); Q only at the
//! Coulomb-allowed k; lk = nullptr: Q only. lk is a Coulomb::LkTable, or the
//! Coulomb::NkTable ExternalLine::L_xnab
template <typename LadderTable>
Multipoles QL_multipoles(const DiracSpinor &a, const DiracSpinor &b,
                         const DiracSpinor &c, const DiracSpinor &d,
                         const Coulomb::QkTable &qk, const LadderTable *lk);

/*!
  @brief W^k_abcd = Q^k_abcd + P^k_abcd for all k of k_minmax_L() (k_step =
         1): the antisymmetrised Coulomb integral [Coulomb::QkTable::W].

  @details
  The exchange part from the multipoles of Q_abdc:
  \f[
    P^k_{abcd} = [k] \sum_l \sixj{j_a}{j_c}{k}{j_b}{j_d}{l} Q^l_{abdc}.
  \f]
  @p sixj must hold every 6j symbol of the four orbitals.
*/
Multipoles W_multipoles(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d,
                        const Coulomb::QkTable &qk, const SixJArray &sixj);

//! Pair-coupled form L^J_{ab,cd} of one integral from its multipoles L^k_abcd
//! (Q, L or QL)
double pair_coupled(int J, const DiracSpinor &a, const DiracSpinor &b,
                    const DiracSpinor &c, const DiracSpinor &d,
                    const Multipoles &L_k, const SixJArray &sixj);

//! L^k_abcd from the pair-coupled L^J_{ab,cd}, J_min <= J <= J_max (inverse
//! of pair_coupled)
double multipole_from_pair_coupled(
  int k, const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
  const DiracSpinor &d, const std::array<double, ladder_max_multipoles> &L_J,
  int J_min, int J_max, const SixJArray &sixj);

//! Particle-hole coupled form L^K_{ad,bc} of one integral from L^k_abcd
double particle_hole_coupled(int K, const DiracSpinor &a, const DiracSpinor &b,
                             const DiracSpinor &c, const DiracSpinor &d,
                             const Multipoles &L_k, const SixJArray &sixj);

//! L^k_abcd from the particle-hole coupled L^K_{ad,bc}, K_min <= K <= K_max
//! (inverse of particle_hole_coupled)
double multipole_from_particle_hole_coupled(
  int k, const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
  const DiracSpinor &d, const std::array<double, ladder_max_multipoles> &L_K,
  int K_min, int K_max, const SixJArray &sixj);

//==============================================================================
/*!
  @brief Orbitals of the last two indices (i,j) of L^k_mnij: the core, then the
         i orbitals not in the core.

  @details
  For an i orbital not in the core, its copy in @p excited is used, so it has
  the same energy as an i line as it has as an intermediate line. i orbitals
  in neither core nor excited are ignored.
*/
std::vector<DiracSpinor>
ladder_ij_orbitals(const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const std::vector<DiracSpinor> &i_orbs);

//! Mask over (i,j) pairs of ij_orbitals [index_i*num_ij + index_j]: required
//! if one is in i_orbs and the other in core (both orders: L_mnib and its
//! partner L_mnbi)
std::vector<char> required_ij_pairs(const std::vector<DiracSpinor> &ij_orbitals,
                                    const std::vector<DiracSpinor> &core,
                                    const std::vector<DiracSpinor> &i_orbs);

//! Positions in `orbitals` grouped by kappa, groups in order of first
//! appearance
std::vector<std::vector<std::size_t>>
group_by_kappa(const std::vector<DiracSpinor> &orbitals);

//==============================================================================
//! L1 right factor QL^J_{rs,ij}/Delta e_{ij,rs}: one matrix per (J, parity)
//! block; rows rs_pairs, columns ij_pairs. external: the line is i or j
std::vector<LinAlg::Matrix<double>> QL_particle_particle_over_de(
  const std::vector<DiracSpinor> &excited,
  const std::vector<DiracSpinor> &ij_orbitals,
  const std::vector<char> &required_ij, const OrbitalPairIndex &rs_pairs,
  const OrbitalPairIndex &ij_pairs, const Coulomb::QkTable &qk,
  const Coulomb::LkTable *lk_prev, const SixJArray &sixj,
  const ExternalLine *external = nullptr);

//! L4 right factor Q^J_{cd,ij}: one matrix per (J, parity) block; rows
//! cd_pairs, columns ij_pairs
std::vector<LinAlg::Matrix<double>>
Q_hole_hole(const std::vector<DiracSpinor> &core,
            const std::vector<DiracSpinor> &ij_orbitals,
            const std::vector<char> &required_ij,
            const OrbitalPairIndex &cd_pairs, const OrbitalPairIndex &ij_pairs,
            const Coulomb::QkTable &qk, const SixJArray &sixj);

//! L1 left factor Q^J_{mn,rs} for one m and the n of one kappa group: one
//! matrix per J (rows n, columns rs_pairs); never stored for all (m,n)
std::vector<LinAlg::Matrix<double>> Q_particle_particle_rows(
  const DiracSpinor &m, const std::vector<std::size_t> &n_group,
  const std::vector<DiracSpinor> &excited, const OrbitalPairIndex &rs_pairs,
  const OrbitalPairIndex &ij_pairs, const Coulomb::QkTable &qk,
  const SixJArray &sixj);

//! L4 left factor QL^J_{mn,cd}/Delta e_{cd,mn} for one m and the n of one
//! kappa group: one matrix per J (rows n, columns cd_pairs). external: the
//! line is m
std::vector<LinAlg::Matrix<double>> QL_hole_hole_rows_over_de(
  const DiracSpinor &m, const std::vector<std::size_t> &n_group,
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const OrbitalPairIndex &cd_pairs, const OrbitalPairIndex &ij_pairs,
  const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
  const SixJArray &sixj, const ExternalLine *external = nullptr);

//! L2 left factor QL^K_{mj,rc}/Delta e_{cj,mr}: rows (m excited, j in
//! ij_orbitals) by ph_pairs; columns (r excited, c core) by rc_pairs.
//! external: the line is j (ExternalLine::L_mnxb) or m (ExternalLine::L_xnab)
std::vector<LinAlg::Matrix<double>> QL_particle_hole_over_de(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &ph_pairs,
  const OrbitalPairIndex &rc_pairs, const Coulomb::QkTable &qk,
  const Coulomb::LkTable *lk_prev, const SixJArray &sixj,
  const ExternalLine *external = nullptr);

//! L2 right factor Q^K_{cr,ni}: rows (r excited, c core) by rc_pairs; columns
//! (n in n_orbitals, i in i_orbitals) by ni_pairs
std::vector<LinAlg::Matrix<double>> Q_particle_hole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &n_orbitals,
  const std::vector<DiracSpinor> &i_orbitals, const OrbitalPairIndex &rc_pairs,
  const OrbitalPairIndex &ni_pairs, const Coulomb::QkTable &qk,
  const SixJArray &sixj);

//! P left factor QL^k_mric/Delta e_{ic,mr}, k of both parities, in the block
//! of rank k: rows (m excited, i in ij_orbitals) by ph_pairs; columns (r
//! excited, c core) by rc_pairs. No recoupling: the multipole k couples m
//! with i, and r with c
std::vector<LinAlg::Matrix<double>> QL_multipole_over_de(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &ph_pairs,
  const OrbitalPairIndex &rc_pairs, const Coulomb::QkTable &qk,
  const Coulomb::LkTable *lk_prev);

//! P right factor Q^k_cnrj, in the block of rank k: rows (r excited, c core)
//! by rc_pairs; columns (n excited, j in ij_orbitals) by ph_pairs. Coulomb k
//! only: the blocks of the other parity are empty. No recoupling
std::vector<LinAlg::Matrix<double>> Q_multipole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &rc_pairs,
  const OrbitalPairIndex &ph_pairs, const Coulomb::QkTable &qk);

/*!
  @brief Hole-particle factor of V^(a), Y_{rc,nj} = -delta_cj y^(c)_nr, in the
         blocks of every rank K allowed for (n,j): rows (r excited, c core) by
         rc_pairs; columns (n excited, j in ij_orbitals) by ph_pairs.

  @details
  \f[
    y^{(c)}_{nr} = \langle n|y^0_{cc}|r\rangle = R^0_{ncrc},
  \f]
  nonzero only for \f$\kappa_r = \kappa_n\f$: the monopole of the
  interaction with the hole c, which is left unchanged (orbital and
  projection). No recoupling: it does not change the rank of the pair, so it
  enters every block of both parities. Columns (n,j) with j not in the core
  are zero. ij_orbitals must start with the core, in the same order
  [ladder_ij_orbitals()].
*/
std::vector<LinAlg::Matrix<double>> Vhp_matrix(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &rc_pairs,
  const OrbitalPairIndex &ph_pairs, const Coulomb::QkTable &qk);

//! Left factor of calculate_Lk_vnab_F(): QL^k_vrac/Delta e_{ac,vr}, k of both
//! parities, in the block of rank k: rows (v in valence, a core) by va_pairs;
//! columns (r excited, c core) by rc_pairs. QL = Q + L, with L from
//! L_vnab_prev (nullptr: Q only). As QL_multipole_over_de(), with the rows
//! (v,a) of the closing amplitude
std::vector<LinAlg::Matrix<double>> QL_multipole_vrac_over_de(
  const std::vector<DiracSpinor> &valence,
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const OrbitalPairIndex &va_pairs, const OrbitalPairIndex &rc_pairs,
  const Coulomb::QkTable &qk, const Coulomb::NkTable *L_vnab_prev);

//==============================================================================
/*!
  @brief Particle-particle ladder L1^J_{mn,ij} for one m and the n of one
         kappa group: one matrix per J (rows n, columns ij_pairs).

  @details
  \f[
    L1^J_{mn,ij} = \sum_{rs} Q^J_{mn,rs}
                   \frac{QL^J_{rs,ij}}{\Delta\epsilon_{ij,rs}},
  \f]
  one matrix product per J: Q_particle_particle_rows() times
  @p QLJ_rs_ij_over_de [QL_particle_particle_over_de()].
*/
std::vector<LinAlg::Matrix<double>> L1_particle_particle(
  const DiracSpinor &m, const std::vector<std::size_t> &n_group,
  const std::vector<DiracSpinor> &excited, const OrbitalPairIndex &rs_pairs,
  const OrbitalPairIndex &ij_pairs,
  const std::vector<LinAlg::Matrix<double>> &QLJ_rs_ij_over_de,
  const Coulomb::QkTable &qk, const SixJArray &sixj);

/*!
  @brief Hole-hole ladder L4^J_{mn,ij} for one m and the n of one kappa group:
         one matrix per J (rows n, columns ij_pairs).

  @details
  \f[
    L4^J_{mn,ij} = \sum_{cd} \frac{QL^J_{mn,cd}}{\Delta\epsilon_{cd,mn}}
                   Q^J_{cd,ij},
  \f]
  one matrix product per J: QL_hole_hole_rows_over_de() times @p QJ_cd_ij
  [Q_hole_hole()]. @p external is passed to QL_hole_hole_rows_over_de().
*/
std::vector<LinAlg::Matrix<double>>
L4_hole_hole(const DiracSpinor &m, const std::vector<std::size_t> &n_group,
             const std::vector<DiracSpinor> &excited,
             const std::vector<DiracSpinor> &core,
             const OrbitalPairIndex &cd_pairs, const OrbitalPairIndex &ij_pairs,
             const std::vector<LinAlg::Matrix<double>> &QJ_cd_ij,
             const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
             const SixJArray &sixj, const ExternalLine *external = nullptr);

/*!
  @brief Particle-hole ladder L2^K_{mj,ni}: rows (m,j) in mj_pairs (m
         excited, j in j_orbitals), columns (n,i) in ni_pairs (n in
         n_orbitals, i in i_orbitals).

  @details
  One matrix product per (K, parity) block:
  \f[
    L2^K_{mj,ni} = (-1)^K \sum_{rc}
      \frac{QL^K_{mj,rc}}{\Delta\epsilon_{cj,mr}} Q^K_{cr,ni},
  \f]
  r excited, c core: QL_particle_hole_over_de() times Q_particle_hole(). The
  sign (-1)^K is the net sign of the minus sign of L2 in the m-scheme ladder
  equation and the phases of its three particle-hole recouplings: of
  Q^u_{cnir} [pairs (c,i),(n,r)] to Q^K_{cr,ni}, of QL^l_{mrcj} [pairs
  (m,c),(r,j)] to QL^K_{mj,rc}, and of L2^K_{mj,ni} back to L2^k_{mnij} [pairs
  (m,i),(n,j)].

  In the ladder equation the matrix is square (n_orbitals = excited,
  i_orbitals = j_orbitals, ni_pairs = mj_pairs), and L3 is read from the same
  matrix: L3^k_mnij from L2^K_{ni,mj}. n and i enter only the bare Q^K, so
  they may be any orbitals, indexed by position (one may have the key of a
  core orbital): the matrix is then rectangular, and gives L2 only
  [lowest_order_Lk_mnxb(), lowest_order_Lk_xnab()]. @p external is passed to
  QL_particle_hole_over_de().
*/
std::vector<LinAlg::Matrix<double>> L2_particle_hole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &j_orbitals, const OrbitalPairIndex &mj_pairs,
  const std::vector<DiracSpinor> &n_orbitals,
  const std::vector<DiracSpinor> &i_orbitals, const OrbitalPairIndex &ni_pairs,
  const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
  const SixJArray &sixj, const ExternalLine *external = nullptr);

//! Adds each matrix of term to the matrix of the same index in sum; an empty
//! matrix in sum takes the term's matrix. Both indexed the same way
void add_blocks(std::vector<LinAlg::Matrix<double>> *sum,
                std::vector<LinAlg::Matrix<double>> term);

//! (-1)^rank times the product left*right of each (rank, parity) block, all
//! indexed by OrbitalPairIndex::block(): the matrix products of the pieces of
//! P. Empty where either factor is empty
std::vector<LinAlg::Matrix<double>>
signed_block_products(const std::vector<LinAlg::Matrix<double>> &left,
                      const std::vector<LinAlg::Matrix<double>> &right);

/*!
  @brief Piece (1) of P, the polarisation loop, one side: P1^k_{mi,nj}, rows
         (m,i) and columns (n,j), in the block of rank k.

  @details
  \f[
    P1^k_{mi,nj} = \frac{(-1)^k}{[k]} \sum_{rc}
      \frac{QL^k_{mric}}{\Delta\epsilon_{ic,mr}} Q^k_{cnrj},
  \f]
  r excited, c core: @p QLk_mi_rc_over_de [QL_multipole_over_de()] times
  @p Qk_rc_nj [Q_multipole()], one matrix product per block. See
  calculate_Lk_mnib().
*/
std::vector<LinAlg::Matrix<double>>
P1_multipole(const std::vector<LinAlg::Matrix<double>> &QLk_mi_rc_over_de,
             const std::vector<LinAlg::Matrix<double>> &Qk_rc_nj);

/*!
  @brief Piece (2) of P, one side: P2^k_{mi,nj}, rows (m,i) and columns (n,j),
         in the block of rank k.

  @details
  \f[
    P2^k_{mi,nj} = (-1)^k \sum_{rc}
      \frac{QL^k_{mric}}{\Delta\epsilon_{ic,mr}} Q^K_{cr,nj},
    \qquad K = k,
  \f]
  r excited, c core: @p QLk_mi_rc_over_de [QL_multipole_over_de()] times
  @p QK_cr_nj [Q_particle_hole(), the L2 right factor], one matrix product
  per block. See calculate_Lk_mnib().
*/
std::vector<LinAlg::Matrix<double>>
P2_multipole(const std::vector<LinAlg::Matrix<double>> &QLk_mi_rc_over_de,
             const std::vector<LinAlg::Matrix<double>> &QK_cr_nj);

/*!
  @brief Piece (3) of P, one side: P3^k_{mi,nj}, rows (m,i) and columns (n,j),
         in the block of rank k.

  @details
  \f[
    P3^k_{mi,nj} = (-1)^k \sum_{rc}
      \frac{QL^K_{mi,rc}}{\Delta\epsilon_{ic,mr}} Q^k_{cnrj},
    \qquad K = k,
  \f]
  r excited, c core: @p QLK_mi_rc_over_de [QL_particle_hole_over_de(), the L2
  left factor] times @p Qk_rc_nj [Q_multipole()], one matrix product per
  block. See calculate_Lk_mnib().
*/
std::vector<LinAlg::Matrix<double>>
P3_multipole(const std::vector<LinAlg::Matrix<double>> &QLK_mi_rc_over_de,
             const std::vector<LinAlg::Matrix<double>> &Qk_rc_nj);

/*!
  @brief Hole-particle interaction of V^(a), one side: H^k_{mi,nj}, rows (m,i)
         and columns (n,j), in the block of rank k.

  @details
  \f[
    H^k_{mi,nj} = \sum_{rc} \frac{QL^k_{mric}}{\Delta\epsilon_{ic,mr}}
      Y_{rc,nj}
      = -\sum_{r} y^{(j)}_{nr} \frac{QL^k_{mrij}}{\Delta\epsilon_{ij,mr}},
  \f]
  r excited (\f$\kappa_r = \kappa_n\f$), j core: @p QLk_mi_rc_over_de
  [QL_multipole_over_de()] times @p Y_rc_nj [Vhp_matrix()], one matrix
  product per block, with no sign and no 1/[k]. See calculate_Lk_mnib().
*/
std::vector<LinAlg::Matrix<double>>
Vhp_multipole(const std::vector<LinAlg::Matrix<double>> &QLk_mi_rc_over_de,
              const std::vector<LinAlg::Matrix<double>> &Y_rc_nj);

/*!
  @brief The pieces of P switched on in @p diagrams, one side, summed:
         P^k_{mi,nj}, rows (m,i) and columns (n,j) by ph_pairs (m,n excited;
         i,j in ij_orbitals), in the block of rank k.

  @details
  P1_multipole() + P2_multipole() + P3_multipole() + Vhp_multipole(); each
  factor is formed once, and only if a piece switched on needs it. The table
  entry is P^k_mnij = P^k_{mi,nj} + P^k_{nj,mi} [store_Lk_mnij()]. See
  calculate_Lk_mnib().
*/
std::vector<LinAlg::Matrix<double>> P_multipole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &ph_pairs,
  const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
  const SixJArray &sixj, const LadderDiagrams &diagrams);

//==============================================================================
//! (L1+L4)^J_{mn,ij} for J_min <= J <= J_max (zero elsewhere): row `row` (the
//! n) of L1_L4_J_mn_ij[J], column of (ij_orbitals[index_i],
//! ij_orbitals[index_j]) in ij_pairs. Zero for J with an empty matrix
std::array<double, ladder_max_multipoles>
read_L1_L4_J(const std::vector<LinAlg::Matrix<double>> &L1_L4_J_mn_ij,
             std::size_t row, const OrbitalPairIndex &ij_pairs,
             std::size_t index_i, std::size_t index_j, int J_min, int J_max);

//! {L2^K_{mj,ni}, L2^K_{ni,mj}} for K_min <= K <= K_max (zero elsewhere): the
//! two elements of L2K_mj_ni [L2_particle_hole()] that give L2^k_mnij and
//! L3^k_mnij. m = excited[index_m] etc.; j = ij_orbitals[index_j] etc.
std::pair<std::array<double, ladder_max_multipoles>,
          std::array<double, ladder_max_multipoles>>
read_L2K(const std::vector<LinAlg::Matrix<double>> &L2K_mj_ni,
         const OrbitalPairIndex &ph_pairs, int parity_mj, std::size_t index_m,
         std::size_t index_j, std::size_t index_n, std::size_t index_i,
         int K_min, int K_max);

/*!
  @brief Element L2^K_{mj,ni} of a square or rectangular L2_particle_hole()
         matrix, for K_min <= K <= K_max (zero elsewhere).

  @details
  Row (m,j) in @p mj_pairs, column (n,i) in @p ni_pairs; index_m etc. are
  positions in the orbitals of the pair indices. Gives L2^k_mnij.
*/
std::array<double, ladder_max_multipoles>
read_L2K_element(const std::vector<LinAlg::Matrix<double>> &L2K_mj_ni,
                 const OrbitalPairIndex &mj_pairs,
                 const OrbitalPairIndex &ni_pairs, int parity_mj,
                 std::size_t index_m, std::size_t index_j, std::size_t index_n,
                 std::size_t index_i, int K_min, int K_max);

//! P^k_mnij = P^k_{mi,nj} + P^k_{nj,mi} for k_min <= k <= k_max (zero
//! elsewhere): element (m,i),(n,j) of Pk_mi_nj [P_multipole()], one side, plus
//! its partner [m <-> n, i <-> j], the transposed element. m = excited[index_m]
//! etc.; i = ij_orbitals[index_i] etc.
std::array<double, ladder_max_multipoles>
read_Pk(const std::vector<LinAlg::Matrix<double>> &Pk_mi_nj,
        const OrbitalPairIndex &ph_pairs, int parity_mi, std::size_t index_m,
        std::size_t index_i, std::size_t index_n, std::size_t index_j,
        int k_min, int k_max);

//! New value of a table entry, damped as CoulombTable::update:
//! a L_old + (1-a) L_new, with a = a_damp clamped to [0,1]; L_new if L_old is
//! exactly zero
double damped(double L_old, double L_new, double a_damp);

/*!
  @brief Core relaxation R^k_mnij, for all Coulomb-allowed k (k_step = 2):
         the source term of the core single excitations rho_ma in the doubles
         equation.

  @details
  \f[
    R^k_{mnij} = [i\in{\rm core}] \sum_{r} \rho_{ri}\,Q^k_{mnrj}
               + [j\in{\rm core}] \sum_{r} \rho_{rj}\,Q^k_{mnir}
               - \sum_{c} \rho_{mc}\,Q^k_{cnij}
               - \sum_{c} \rho_{nc}\,Q^k_{mcij},
  \f]
  r excited and c core, each with the kappa of the orbital it replaces. It is
  the first-order change of Q^k_mnij when each core hole relaxes,
  \f$ a \to a + \sum_r \rho_{ra}\,r \f$, and each particle is kept
  orthogonal to the core, \f$ m \to m - \sum_c \rho_{mc}\,c \f$. A valence
  hole is not changed: its single excitations are the Brueckner orbital, not
  part of the doubles equation. The single excitations are diagonal in kappa
  and the same for every projection, so there is no angular factor, and
  R^k_mnij = R^k_nmji.

  @param index_m,index_n  m = excited[index_m], n = excited[index_n]
  @param index_i,index_j  i = ij_orbitals[index_i], j = ij_orbitals[index_j]
  @param excited          Excited orbitals: the rows of rho_ma
  @param ij_orbitals      Orbitals of the (i,j) indices, starting with the
                          core: the first rho_ma.cols() of them are the core
                          orbitals, the columns of rho_ma
                          [ladder_ij_orbitals()]
  @param qk               Coulomb integral table
  @param rho_ma           Core single excitations rho_ma, rows excited and
                          columns core (m-scheme value; zero unless kappa_m =
                          kappa_a) [core_singles()]
*/
Multipoles relaxation_multipoles(std::size_t index_m, std::size_t index_n,
                                 std::size_t index_i, std::size_t index_j,
                                 const std::vector<DiracSpinor> &excited,
                                 const std::vector<DiracSpinor> &ij_orbitals,
                                 const Coulomb::QkTable &qk,
                                 const LinAlg::Matrix<double> &rho_ma);

/*!
  @brief Writes L^k_mnij = (L1+L4) + L2 + L3 + P [+ R] for one m and the n of
         one kappa group into the table.

  @details
  Only entries already in the table are written. For m = n only the
  normal-ordered (i,j): the two orders are the same integral. A diagram whose
  matrices are empty is not included. The core relaxation R
  [relaxation_multipoles()] is included if @p rho_ma is given.

  @param lk             Ladder table to write to
  @param index_m        m = excited[index_m]
  @param n_group        Positions in excited of the n (one kappa)
  @param excited        Excited orbitals
  @param ij_orbitals    Orbitals of the (i,j) indices [ladder_ij_orbitals()]
  @param required_ij    Required (i,j) pairs [required_ij_pairs()]
  @param ij_pairs       Column index of the (i,j) pairs
  @param ph_pairs       Row and column index of L2K_mj_ni
  @param L1_L4_J_mn_ij  Pair-coupled (L1+L4)^J_{mn,ij}: one matrix per J
                        (indexed by J alone), row r is the n of n_group[r],
                        columns ij_pairs
  @param L2K_mj_ni      Particle-hole coupled L2^K_{mj,ni} [L2_particle_hole()]:
                        indexed by OrbitalPairIndex::block(K, parity)
  @param Pk_mi_nj       P^k_{mi,nj} [P_multipole()], rows and columns
                        ph_pairs: indexed by OrbitalPairIndex::block(k,
                        parity); P^k_mnij = P^k_{mi,nj} + P^k_{nj,mi}
                        [read_Pk()]
  @param qk             Coulomb integral table (for R)
  @param rho_ma         Core single excitations for R [see
                        relaxation_multipoles()]; nullptr: no R
  @param sixj           6j symbols
  @param max_k          Largest multipolarity k written
  @param only_new       If true, only entries that are exactly zero are written
  @param a_damp         Damping factor [see damped()]; not used if only_new
*/
void store_Lk_mnij(Coulomb::LkTable *lk, std::size_t index_m,
                   const std::vector<std::size_t> &n_group,
                   const std::vector<DiracSpinor> &excited,
                   const std::vector<DiracSpinor> &ij_orbitals,
                   const std::vector<char> &required_ij,
                   const OrbitalPairIndex &ij_pairs,
                   const OrbitalPairIndex &ph_pairs,
                   const std::vector<LinAlg::Matrix<double>> &L1_L4_J_mn_ij,
                   const std::vector<LinAlg::Matrix<double>> &L2K_mj_ni,
                   const std::vector<LinAlg::Matrix<double>> &Pk_mi_nj,
                   const Coulomb::QkTable &qk,
                   const LinAlg::Matrix<double> *rho_ma, const SixJArray &sixj,
                   int max_k, bool only_new, double a_damp);

//==============================================================================
//! Adds zero entries for all L^k_mnib (m,n excited; i in i_orbs; b core;
//! k <= max_k), each stored once (L^k_mnib = L^k_nmbi); existing entries are
//! left unchanged
void add_Lk_mnib_entries(Coulomb::LkTable *lk,
                         const std::vector<DiracSpinor> &excited,
                         const std::vector<DiracSpinor> &core,
                         const std::vector<DiracSpinor> &i_orbs, int max_k);

/*!
  @brief One sweep of the ladder equation: all L^k_mnib (m,n excited, i in
         i_orbs, b core) and their partners L^k_nmbi at once.

  @details
  The ladder integral is the sum of four diagrams,
  L^k_mnij = L1 + L2 + L3 + L4, each with one bare Coulomb integral Q and one
  dressed integral QL = Q + L (L from @p lk_prev). Each diagram is a matrix
  product in the coupling scheme in which it is diagonal. Only the orbitals'
  quantum numbers and energies are used: integrals are looked up in the
  tables by DiracSpinor::nk_index().

  Pair coupling, pairs (a,b) and (c,d), for any integral (Q, L or QL):
  \f[
    L^J_{ab,cd} = (-1)^{J+j_a+j_c}\sum_k (-1)^k
                  \sixj{j_a}{j_b}{J}{j_d}{j_c}{k} L^k_{abcd},
  \f]
  \f[
    L^k_{abcd} = (-1)^k [k]\sum_J [J] (-1)^{J+j_a+j_c}
                 \sixj{j_a}{j_b}{J}{j_d}{j_c}{k} L^J_{ab,cd}.
  \f]
  Particle-hole coupling, pairs (a,d) and (b,c):
  \f[
    L^K_{ad,bc} = \sum_k \sixj{j_a}{j_d}{K}{j_b}{j_c}{k} L^k_{abcd},
    \qquad
    L^k_{abcd} = [k] \sum_K [K] \sixj{j_a}{j_d}{K}{j_b}{j_c}{k} L^K_{ad,bc}.
  \f]
  Superscript J always means pair coupling, superscript K particle-hole
  coupling, and superscript k the multipole form. With the energy difference
  \f$ \Delta\epsilon_{ab,cd} = \epsilon_a + \epsilon_b - \epsilon_c
  - \epsilon_d \f$, the diagrams are
  \f[
    L1^J_{mn,ij} = \sum_{rs} Q^J_{mn,rs}
                   \frac{QL^J_{rs,ij}}{\Delta\epsilon_{ij,rs}},
  \f]
  \f[
    L4^J_{mn,ij} = \sum_{cd} \frac{QL^J_{mn,cd}}{\Delta\epsilon_{cd,mn}}
                   Q^J_{cd,ij},
  \f]
  \f[
    L2^K_{mj,ni} = (-1)^K \sum_{rc}
                   \frac{QL^K_{mj,rc}}{\Delta\epsilon_{cj,mr}} Q^K_{cr,ni},
  \f]
  \f[
    L3^k_{mnij} = L2^k_{nmji},
  \f]
  where r,s run over all ordered excited pairs and c,d over ordered core
  pairs. Each diagram has its own function: L1_particle_particle(),
  L2_particle_hole() and L4_hole_hole(); L3 is read from L2^K_{ni,mj}, the same
  matrix, in store_Lk_mnij(). The integrals with four excited orbitals,
  Q^J_{mn,rs}, are never stored: they are formed for one m and one kappa group
  of n at a time. J and K are cut off at ladder_max_k.

  The screening and hole-particle interaction, pieces (1)-(3) of P [see
  Pkmnij(); off by default], are diagonal in the multipole form itself: k
  couples m with i, n with j and r with c. Each is one matrix product at rank
  k, rows (m,i) and columns (n,j), with no recoupling of the result:
  \f[
    P1^k_{mi,nj} = \frac{(-1)^k}{[k]} \sum_{rc}
                   \frac{QL^k_{mric}}{\Delta\epsilon_{ic,mr}} Q^k_{cnrj},
  \f]
  \f[
    P2^k_{mi,nj} = (-1)^k \sum_{rc}
                   \frac{QL^k_{mric}}{\Delta\epsilon_{ic,mr}} Q^K_{cr,nj},
  \f]
  \f[
    P3^k_{mi,nj} = (-1)^k \sum_{rc}
                   \frac{QL^K_{mi,rc}}{\Delta\epsilon_{ic,mr}} Q^k_{cnrj},
  \f]
  \f[
    P^k_{mnij} = (P1 + P2 + P3)^k_{mi,nj} + (P1 + P2 + P3)^k_{nj,mi},
  \f]
  with K = k. Pieces (2) and (3) contain the exchange multipoles of
  (c,n,r,j) and (m,r,i,c): these are [k] times the particle-hole coupled
  forms at K = k, the two factors of L2. The factor (-1)^k/[k] is from the sum
  over the projections of r and c; the minus signs of pieces (2) and (3)
  cancel those of the exchange multipoles. P^k_{nj,mi}, the partner [m <-> n,
  i <-> j], is the transposed element, at the same k. Each piece has its own
  function: P1_multipole(), P2_multipole() and P3_multipole()
  [P_multipole()].

  The hole-particle interaction of V^(a) [hole_particle; see Vhp(); off by
  default] has the same form, with the factor
  Y_{rc,nj} = -\delta_{cj} y^{(c)}_{nr} [Vhp_matrix()] in place of the Coulomb
  integral:
  \f[
    H^k_{mi,nj} = \sum_{rc}
                  \frac{QL^k_{mric}}{\Delta\epsilon_{ic,mr}} Y_{rc,nj}
                = -\sum_r y^{(j)}_{nr}
                  \frac{QL^k_{mrij}}{\Delta\epsilon_{ij,mr}},
  \f]
  [Vhp_multipole()]: -y^0_jj acts on n, with the hole j unchanged, and the
  rank is unchanged. Its partner H^k_{nj,mi} acts on m, and is zero unless i
  is a core orbital (Y has no columns with a valence hole).

  The core relaxation R [relaxation_multipoles(); with @p rho_ma, off by
  default] is a source term: the core single excitations act on the bare
  Coulomb integral, R^k_mnij = R[rho_ma], with no dressed integral. It is
  added to the new entry with the diagrams, before the damping; from the next
  sweep on it is part of the dressed integral QL, as every entry is. A hole
  is dressed only if it is a core orbital, at either position of the stored
  entry: in the valence entries the valence hole is stored at i or at j.

  Factors shared by all (m,n) are indexed by OrbitalPairIndex::block(J or K,
  parity); factors for one m and one kappa group of n by J alone (their parity
  is fixed).

  @param lk        Ladder table to write to. Only entries already in the table
                   are written [see add_Lk_mnib_entries()]
  @param qk        Coulomb integral table
  @param excited   Excited orbitals
  @param core      Core orbitals
  @param i_orbs    Orbitals for the i index (those not in core must be in
                   excited; others are ignored)
  @param diagrams  Which diagrams to include. L2 and L3 must go together:
                   L_mnij and L_nmji are one table entry, and L3^k_mnij =
                   L2^k_nmji. P and hole_particle only without @p external;
                   hole_particle not with P2 (it is part of piece (2))
  @param lk_prev   Ladder table for the dressed integral QL = Q + L; nullptr
                   for the lowest order. Must not be lk
  @param max_k     Largest multipolarity k of the integrals to write
  @param only_new  If true, only entries whose current value is exactly zero
                   are written
  @param a_damp    Damping factor a in [0,1), used if only_new is false:
                   L -> a L + (1-a) L_new, with L the current entry of lk. No
                   damping for entries whose current value is exactly zero
  @param external  Optional: the i orbitals are the external line of
                   Sigma_L(e_v) [ExternalLine::L_mnxb; see calculate_Lk_mnxb()]
  @param rho_ma    Optional: core single excitations rho_ma, rows excited
                   and columns core [core_singles()]; includes the core
                   relaxation R. Not with @p external

  @note Integrals that vanish because no term is allowed may come out as
        rounding noise (relative to the terms that cancel) instead of exactly
        zero.
*/
void calculate_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &i_orbs,
                       const LadderDiagrams &diagrams,
                       const Coulomb::LkTable *lk_prev, int max_k,
                       bool only_new, double a_damp,
                       const ExternalLine *external = nullptr,
                       const LinAlg::Matrix<double> *rho_ma = nullptr);

//==============================================================================
//! Adds zero entries for all L'^k_vnab (v in valence, n excited, a,b core,
//! k <= max_k), keyed (v,n,a,b), to the closing table of
//! calculate_Lk_vnab_F(); existing entries are left unchanged
void add_Lk_vnab_F_entries(Coulomb::NkTable *L_vnab,
                           const std::vector<DiracSpinor> &excited,
                           const std::vector<DiracSpinor> &core,
                           const std::vector<DiracSpinor> &valence, int max_k);

/*!
  @brief One sweep of the closing amplitude of x_F, for all valence states at
         once: x_F_vnab without the hole-particle term on the pair (v,a).

  @details
  In the second term of the direct (Q) closing of x_F [de_valence_F()],
  \f$ \sum_{nab} g_{abvn} x_{vnab}/(\epsilon_{vn} - \epsilon_{ab}) \f$, the
  pair (v,a) is part of the valence line, not a polarisation loop: V^(a)
  never acts on it, at any step after the pair is created. The core amplitude
  x_F_mnab (m = v included) has the hole-particle term on both pairs, since
  it also enters as an inner amplitude, where (v,a) is a loop. So the closing
  has its own amplitude x', for each valence state v:
  \f[
    x'_{vnab} = \sum_{rc}\frac{g_{cnrb}\,X'_{vrac}}{\epsilon_{ac}-\epsilon_{vr}}
      + \sum_{rc}\frac{g_{cvra}\,X^F_{nrbc}}{\epsilon_{bc}-\epsilon_{nr}}
      - \sum_r \frac{y^{(b)}_{nr}\,X'_{vrab}}{\epsilon_{ab}-\epsilon_{vr}},
  \f]
  with X' = g + x' and X^F = g + x_F, r excited and c core. The first term
  is the loop on (n,b), which keeps (v,a); the second the loop on (v,a), the
  step that creates (v,a) from a core amplitude; the third the hole-particle
  term on (n,b). In the multipole form, as calculate_Lk_mnib():
  \f[
    L'^k_{vnab} = \Big[A'\Big(\frac{(-1)^k}{[k]} D + Y\Big)\Big]_{va,nb}
      + \Big[A\,\frac{(-1)^k}{[k]} D\Big]_{nb,va},
  \f]
  in the block of rank k: A' = QL'^k_vrac/De from @p L_vnab_prev
  [QL_multipole_vrac_over_de()], A = QL^k_nrbc/De from @p lk_F
  [QL_multipole_over_de()], D = Q^k_cnrb [Q_multipole()] and Y
  [Vhp_matrix()]. Without Vhp, x' is the core amplitude x_F_vnab.

  The table is a Coulomb::NkTable: with n = w another valence state, x'_vwab
  (of v) and x'_wvba (of w) are different amplitudes.

  @param L_vnab       Closing table to write to, keyed (v,n,a,b). Only
                      entries already in the table are written [see
                      add_Lk_vnab_F_entries()]
  @param qk           Coulomb integral table
  @param excited      Excited orbitals
  @param core         Core orbitals
  @param valence      Valence orbitals: members of excited
  @param diagrams     P1 (polarisation loop) and hole_particle only
  @param lk_F         x_F [calculate_Lk_mnib() with P1 and hole_particle]:
                      its core entries; nullptr for the lowest order
  @param L_vnab_prev  Closing table for X' = g + x'; nullptr for the lowest
                      order. Must not be L_vnab
  @param a_damp       Damping factor [see damped()]
*/
void calculate_Lk_vnab_F(Coulomb::NkTable *L_vnab, const Coulomb::QkTable &qk,
                         const std::vector<DiracSpinor> &excited,
                         const std::vector<DiracSpinor> &core,
                         const std::vector<DiracSpinor> &valence,
                         const LadderDiagrams &diagrams,
                         const Coulomb::LkTable *lk_F,
                         const Coulomb::NkTable *L_vnab_prev, double a_damp);

//==============================================================================
/*!
  @brief Factor of a sum over the projections of a closed shell:
         f_xp = (-1)^(j_x + j_p + 1) sqrt([j_p]/[j_x]).

  @details
  For kappa_x = kappa_y and kappa_p = kappa_q, and one projection of x (the
  same for y):
  \f[
    \sum_{\mu_p} \widetilde g_{x p y q} = f_{xp}\,W^0_{xpyq},
  \f]
  with the projection of q that of p, and \f$ \widetilde g_{xpyq} = g_{xpyq}
  - g_{xpqy} \f$. Independent of the projection of x.
*/
double closed_shell_factor(const DiracSpinor &x, const DiracSpinor &p);

/*!
  @brief Source of the core single excitations rho_ma: the terms of the
         singles equation for the core that contain the doubles.

  @details
  The singles equation of the SD method for the core, m excited (valence
  states included) and a core, is
  \f[
    (\epsilon_a - \epsilon_m)\rho_{ma} = \sum_{bn} \widetilde g_{mban}
      \rho_{nb} + \sum_{bnr} g_{mbnr}\,\widetilde\rho_{nrab}
      - \sum_{bcn} g_{bcan}\,\widetilde\rho_{mnbc},
  \f]
  with \f$ \rho_{mnab} = X_{mnab}/\Delta\epsilon_{ab,mn} \f$, X = g + x the
  dressed integral, and \f$ \widetilde\rho_{mnab} = \rho_{mnab} -
  \rho_{mnba} \f$. This returns the last two terms; the first is added by
  core_singles(). In the multipole form, with X^k = Q^k + L^k (L from @p lk,
  both parities of k) and W = Q + P [W_multipoles()]:
  \f[
    pp_{ma} = \frac{1}{[j_a]} \sum_k \frac{1}{[k]} \sum_{bnr}
      \frac{W^k_{mbnr}\,X^k_{nrab}}{\Delta\epsilon_{ab,nr}},
  \f]
  \f[
    hh_{ma} = -\frac{1}{[j_a]} \sum_k \frac{1}{[k]} \sum_{bcn}
      \frac{W^k_{bcan}\,X^k_{mnbc}}{\Delta\epsilon_{bc,mn}},
  \f]
  b, c core and n, r excited. Single excitations are diagonal in kappa and
  the same for every projection: zero unless kappa_m = kappa_a. With X = g
  (lk = nullptr) the source is the second-order correlation potential
  Sigma_ma, and rho_ma = Sigma_ma/(e_a - e_m) is the lowest order.

  @param qk       Coulomb integral table
  @param excited  Excited orbitals: rows
  @param core     Core orbitals: columns
  @param lk       Ladder table x, core entries L^k_mnab (nullptr: X = g)
  @return pp + hh, rows excited and columns core
*/
LinAlg::Matrix<double> core_singles_source(
  const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
  const std::vector<DiracSpinor> &core, const Coulomb::LkTable *lk);

/*!
  @brief Core single excitations rho_ma from the source and the previous
         rho_ma: one step of the singles equation for the core.

  @details
  \f[
    \rho_{ma} = \frac{1}{\epsilon_a - \epsilon_m}\Big( s_{ma}
      + \sum_{bn} f_{ab}\,W^0_{mban}\,\rho^{\rm prev}_{nb} \Big),
    \qquad
    f_{ab} = (-1)^{j_a+j_b+1}\sqrt{[j_b]/[j_a]},
  \f]
  s the source [core_singles_source()], b core and n excited with kappa_n =
  kappa_b: the m-scheme term \f$ \sum_{bn} \widetilde g_{mban}\rho_{nb} \f$,
  summed over the projections of b. Cheap: the singles equation for a fixed
  x is solved by repeating this with one source.

  @param source    Source s_ma, rows excited and columns core
  @param qk        Coulomb integral table
  @param excited   Excited orbitals: rows
  @param core      Core orbitals: columns
  @param rho_prev  Previous rho_ma (nullptr: no first term)
  @return rho_ma, rows excited and columns core (m-scheme value; zero unless
          kappa_m = kappa_a)
*/
LinAlg::Matrix<double> core_singles(const LinAlg::Matrix<double> &source,
                                    const Coulomb::QkTable &qk,
                                    const std::vector<DiracSpinor> &excited,
                                    const std::vector<DiracSpinor> &core,
                                    const LinAlg::Matrix<double> *rho_prev);

//! Change of the core single excitations, relative to their size:
//! max|rho_ma - rho_prev| / max|rho_ma| (zero if rho_ma is zero). Both of one
//! size
double core_singles_change(const LinAlg::Matrix<double> &rho_ma,
                           const LinAlg::Matrix<double> &rho_prev);

} // namespace MBPT
