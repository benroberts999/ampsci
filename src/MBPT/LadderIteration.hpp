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

//! Which ladder diagrams to include. L2 and L3 go together: alone, neither
//! has the table symmetry L_mnij = L_nmji
struct LadderDiagrams {
  //! Particle-particle ladder
  bool L1{true};
  //! Particle-hole ladders [L3^k_mnij = L2^k_nmji]
  bool L2_L3{true};
  //! Hole-hole ladder
  bool L4{true};
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
//! (n excited, i in ij_orbitals) by ph_pairs
std::vector<LinAlg::Matrix<double>> Q_particle_hole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &rc_pairs,
  const OrbitalPairIndex &ph_pairs, const Coulomb::QkTable &qk,
  const SixJArray &sixj);

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
  @brief Particle-hole ladder L2^K_{mj,ni} for all (m,j), (n,i) in ph_pairs
         (m,n excited; i,j in ij_orbitals).

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

  L3 is read from the same matrix: L3^k_mnij from L2^K_{ni,mj}.
  @p external is passed to QL_particle_hole_over_de().
*/
std::vector<LinAlg::Matrix<double>> L2_particle_hole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &ph_pairs,
  const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
  const SixJArray &sixj, const ExternalLine *external = nullptr);

//! Adds each matrix of term to the matrix of the same index in sum; an empty
//! matrix in sum takes the term's matrix. Both indexed the same way
void add_blocks(std::vector<LinAlg::Matrix<double>> *sum,
                std::vector<LinAlg::Matrix<double>> term);

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

//! New value of a table entry, damped as CoulombTable::update:
//! a L_old + (1-a) L_new, with a = a_damp clamped to [0,1]; L_new if L_old is
//! exactly zero
double damped(double L_old, double L_new, double a_damp);

/*!
  @brief Writes L^k_mnij = (L1+L4) + L2 + L3 for one m and the n of one kappa
         group into the table.

  @details
  Only entries already in the table are written. For m = n only the
  normal-ordered (i,j): the two orders are the same integral.

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
                   const SixJArray &sixj, int max_k, bool only_new,
                   double a_damp);

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
  @param diagrams  Which diagrams to include
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
                       const ExternalLine *external = nullptr);

} // namespace MBPT
