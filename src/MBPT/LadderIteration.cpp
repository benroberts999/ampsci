#include "MBPT/LadderIteration.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/Ladder.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

namespace MBPT {

//==============================================================================
std::vector<LadderOrbital>
ladder_orbitals(const std::vector<DiracSpinor> &orbitals) {
  std::vector<LadderOrbital> out;
  out.reserve(orbitals.size());
  for (const auto &Fn : orbitals) {
    out.push_back(
      {Fn.nk_index(), Fn.nk_index(), Fn.kappa(), Fn.twoj(), Fn.l(), Fn.en()});
  }
  return out;
}

//==============================================================================
int max_twoj(const std::vector<LadderOrbital> &orbitals) {
  int max_tj = 0;
  for (const auto &x : orbitals) {
    max_tj = std::max(max_tj, x.twoj);
  }
  return max_tj;
}

//==============================================================================
SixJArray::SixJArray(int max_tj, int max_rank)
  : m_num_j(std::size_t(max_tj + 1) / 2),
    m_num_rank(std::size_t(max_rank) + 1),
    m_data(m_num_j * m_num_j * m_num_j * m_num_j * m_num_rank * m_num_rank) {
#pragma omp parallel for collapse(2)
  for (std::size_t ia = 0; ia < m_num_j; ++ia) {
    for (std::size_t ib = 0; ib < m_num_j; ++ib) {
      for (std::size_t c = 0; c < m_num_rank; ++c) {
        for (std::size_t id = 0; id < m_num_j; ++id) {
          for (std::size_t ie = 0; ie < m_num_j; ++ie) {
            for (std::size_t f = 0; f < m_num_rank; ++f) {
              m_data[index(ia, ib, c, id, ie, f)] =
                Angular::sixj_2(int(2 * ia + 1), int(2 * ib + 1), int(2 * c),
                                int(2 * id + 1), int(2 * ie + 1), int(2 * f));
            }
          }
        }
      }
    }
  }
}

//==============================================================================
OrbitalPairIndex::OrbitalPairIndex(const std::vector<LadderOrbital> &first,
                                   const std::vector<LadderOrbital> &second,
                                   int max_rank,
                                   const std::vector<char> &allowed)
  : m_size_second(second.size()),
    m_max_rank(max_rank),
    m_index(std::size_t(max_rank) + 1,
            std::vector<int>(first.size() * second.size(), -1)),
    m_count(2 * (std::size_t(max_rank) + 1), 0) {
  assert(allowed.empty() || allowed.size() == first.size() * second.size());
  for (std::size_t ix = 0; ix < first.size(); ++ix) {
    for (std::size_t iy = 0; iy < second.size(); ++iy) {
      if (!allowed.empty() && allowed[ix * m_size_second + iy] == 0)
        continue;
      const auto &x = first[ix];
      const auto &y = second[iy];
      const auto parity = pair_parity(x, y);
      const auto [rank_min, rank_max] = Coulomb::k_minmax_tj(x.twoj, y.twoj);
      for (int rank = rank_min; rank <= std::min(rank_max, max_rank); ++rank) {
        m_index[std::size_t(rank)][ix * m_size_second + iy] =
          int(m_count[block(rank, parity)]);
        ++m_count[block(rank, parity)];
      }
    }
  }
}

//==============================================================================
Multipoles coulomb_multipoles(const LadderOrbital &a, const LadderOrbital &b,
                              const LadderOrbital &c, const LadderOrbital &d,
                              const Coulomb::QkTable &qk) {
  Multipoles Q;
  std::tie(Q.k_min, Q.k_max) =
    Coulomb::k_minmax_Q(a.kappa, b.kappa, c.kappa, d.kappa);
  if (Q.k_max < Q.k_min)
    return Q;
  assert(Q.k_max < int(ladder_max_multipoles));
  const auto key =
    qk.NormalOrder(a.nk_index, b.nk_index, c.nk_index, d.nk_index);
  for (int k = Q.k_min; k <= Q.k_max; ++k) {
    Q.value[std::size_t(k)] = ((k - Q.k_min) % 2 == 0) ? qk.Q(k, key) : 0.0;
  }
  return Q;
}

//==============================================================================
Multipoles rung_multipoles(const LadderOrbital &a, const LadderOrbital &b,
                           const LadderOrbital &c, const LadderOrbital &d,
                           const Coulomb::QkTable &qk,
                           const Coulomb::LkTable *lk,
                           const Coulomb::LkTable *lk_other) {
  Multipoles T;
  std::tie(T.k_min, T.k_max) = k_minmax_L(a.kappa, b.kappa, c.kappa, d.kappa);
  if (T.k_max < T.k_min)
    return T;
  assert(T.k_max < int(ladder_max_multipoles));
  const auto [q_min, q_max] =
    Coulomb::k_minmax_Q(a.kappa, b.kappa, c.kappa, d.kappa);
  const auto Qkey =
    qk.NormalOrder(a.nk_index, b.nk_index, c.nk_index, d.nk_index);
  // same key for both ladder tables (the normal ordering is a property of the
  // table symmetry, not of the table)
  const Coulomb::LkTable *const any_lk = lk ? lk : lk_other;
  const auto Lkey =
    any_lk ? any_lk->NormalOrder(a.label, b.label, c.label, d.label) : 0ul;
  for (int k = T.k_min; k <= T.k_max; ++k) {
    const bool has_Q = k >= q_min && k <= q_max && (k - q_min) % 2 == 0;
    T.value[std::size_t(k)] = (has_Q ? qk.Q(k, Qkey) : 0.0) +
                              (lk ? lk->Q(k, Lkey) : 0.0) +
                              (lk_other ? lk_other->Q(k, Lkey) : 0.0);
  }
  return T;
}

//==============================================================================
LadderIntegralSet
ladder_integrals_mnib(const std::vector<DiracSpinor> &excited,
                      const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &i_orbs) {
  LadderIntegralSet integrals;
  integrals.m_orbitals = ladder_orbitals(excited);
  integrals.symmetric_mn = true;

  // "ij" orbitals: every orbital that appears in the (i, j) slots of the
  // required L^k_mnij: the core, and any i orbitals that are not core (valence
  // states; these are also members of the excited set). An integral is
  // required if one of (i, j) is an i orbital and the other is core [both
  // orders are needed, since L^k_mnij = L^k_nmji].
  const auto core_orbitals = ladder_orbitals(core);
  integrals.ij_orbitals = core_orbitals;
  std::vector<char> is_core(core.size(), 1);
  std::vector<char> is_i_orb(core.size(), 0);
  for (const auto &Fi : i_orbs) {
    const auto same_orbital = [&Fi](const LadderOrbital &x) {
      return x.nk_index == Fi.nk_index();
    };
    const auto in_core =
      std::find_if(core_orbitals.cbegin(), core_orbitals.cend(), same_orbital);
    if (in_core != core_orbitals.cend()) {
      is_i_orb[std::size_t(std::distance(core_orbitals.cbegin(), in_core))] = 1;
      continue;
    }
    const auto in_excited = std::find_if(
      integrals.m_orbitals.cbegin(), integrals.m_orbitals.cend(), same_orbital);
    if (in_excited != integrals.m_orbitals.cend() &&
        std::none_of(integrals.ij_orbitals.cbegin(),
                     integrals.ij_orbitals.cend(), same_orbital)) {
      integrals.ij_orbitals.push_back(*in_excited);
      is_core.push_back(0);
      is_i_orb.push_back(1);
    }
  }
  integrals.in_i_slot = is_i_orb;
  integrals.in_j_slot = is_core;
  return integrals;
}

//==============================================================================
void add_Lk_entries(Coulomb::LkTable *lk,
                    const std::vector<LadderOrbital> &excited,
                    const LadderIntegralSet &integrals, int max_k) {

  const auto &m_orbitals = integrals.m_orbitals;
  const auto &ij_orbitals = integrals.ij_orbitals;
  const auto num_ij = ij_orbitals.size();

  const auto num_k = std::size_t(max_k) + 1;
  if ((*lk)->size() < num_k) {
    (*lk)->resize(num_k);
  }

  // Each k has its own map: so can insert in parallel over k.
  // The table symmetry is L^k_mnij = L^k_nmji. If m runs over the excited set
  // (symmetric_mn), both (m,n,i,j) and (n,m,j,i) are in the loops below, and
  // are the same entry: take the normal-ordered one.
#pragma omp parallel for
  for (std::size_t k = 0; k < num_k; ++k) {
    auto &Lk_map = (*lk)->at(k);
    // First pass counts (to reserve space); second pass inserts
    for (const bool insert : {false, true}) {
      std::size_t count = 0;
      for (std::size_t ii = 0; ii < num_ij; ++ii) {
        const auto &i = ij_orbitals[ii];
        for (std::size_t ij = 0; ij < num_ij; ++ij) {
          const auto &j = ij_orbitals[ij];
          if (!integrals.required(ii, ij))
            continue;
          for (const auto &m : m_orbitals) {
            for (const auto &n : excited) {
              const auto [k0, kI] =
                k_minmax_L(m.kappa, n.kappa, i.kappa, j.kappa);
              if (int(k) < k0 || int(k) > kI)
                continue;
              if (integrals.symmetric_mn &&
                  !lk->is_NormalOrdered(m.label, n.label, i.label, j.label))
                continue;
              if (insert) {
                Lk_map.insert(
                  {lk->NormalOrder(m.label, n.label, i.label, j.label), 0.0});
              } else {
                ++count;
              }
            }
          }
        }
      }
      if (!insert) {
        Lk_map.reserve(count);
      }
    }
  }
}

//==============================================================================
void add_Lk_mnib_entries(Coulomb::LkTable *lk,
                         const std::vector<DiracSpinor> &excited,
                         const std::vector<DiracSpinor> &core,
                         const std::vector<DiracSpinor> &i_orbs, int max_k) {
  add_Lk_entries(lk, ladder_orbitals(excited),
                 ladder_integrals_mnib(excited, core, i_orbs), max_k);
}

//==============================================================================
void calculate_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                       const Coulomb::LkTable *lk_prev, int max_k,
                       bool only_new, double a_damp) {
  if (excited.empty() || core.empty() || i_orbs.empty())
    return;
  calculate_Lk(lk, qk, ladder_orbitals(excited), ladder_orbitals(core),
               ladder_integrals_mnib(excited, core, i_orbs),
               {true, true, true, include_L4}, lk_prev, nullptr, max_k,
               only_new, a_damp);
}

//==============================================================================
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
                     const Coulomb::LkTable *lk_fixed, const SixJArray &sixj) {
  using Matrix = LinAlg::Matrix<double>;
  const auto max_K = row_pairs.max_rank();
  assert(column_pairs.max_rank() == max_K);
  const OrbitalPairIndex rc_pairs(excited, core, max_K);

  std::vector<Matrix> A_rcqi(row_pairs.num_blocks());
  std::vector<Matrix> B_pjrc(row_pairs.num_blocks());
  std::vector<Matrix> Y_pjqi(row_pairs.num_blocks());
  for (int K = 0; K <= max_K; ++K) {
    for (int parity = 0; parity < 2; ++parity) {
      const auto rows = row_pairs.size(K, parity);
      const auto inner = rc_pairs.size(K, parity);
      const auto cols = column_pairs.size(K, parity);
      if (rows != 0 && inner != 0 && cols != 0) {
        const auto block = OrbitalPairIndex::block(K, parity);
        A_rcqi[block] = Matrix(inner, cols);
        B_pjrc[block] = Matrix(rows, inner);
        Y_pjqi[block] = Matrix(rows, cols);
      }
    }
  }

  // A^K_{rc,qi} = sum_u {jc jr K; jq ji u} Q^u_cqir
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ir = 0; ir < excited.size(); ++ir) {
    const auto &r = excited[ir];
    for (std::size_t ic = 0; ic < core.size(); ++ic) {
      const auto &c = core[ic];
      const auto parity = pair_parity(r, c);
      for (std::size_t iq = 0; iq < column_particles.size(); ++iq) {
        const auto &q = column_particles[iq];
        for (std::size_t ii = 0; ii < column_holes.size(); ++ii) {
          const auto &i = column_holes[ii];
          if (pair_parity(q, i) != parity)
            continue;
          const auto Q = coulomb_multipoles(c, q, i, r, qk);
          if (Q.k_max < Q.k_min)
            continue;
          const auto K_min =
            std::max(std::abs(r.twoj - c.twoj), std::abs(q.twoj - i.twoj)) / 2;
          const auto K_max =
            std::min({r.twoj + c.twoj, q.twoj + i.twoj, 2 * max_K}) / 2;
          for (int K = K_min; K <= K_max; ++K) {
            auto &A = A_rcqi[OrbitalPairIndex::block(K, parity)];
            if (A.rows() == 0)
              continue;
            double AK = 0.0;
            for (int u = Q.k_min; u <= Q.k_max; u += 2) {
              AK += sixj.get(c.twoj, r.twoj, K, q.twoj, i.twoj, u) *
                    Q.value[std::size_t(u)];
            }
            const auto row = std::size_t(rc_pairs.index(K, ir, ic));
            const auto col = std::size_t(column_pairs.index(K, iq, ii));
            A(row, col) = AK;
          }
        }
      }
    }
  }

  // B^K_{pj,rc} = [K] sum_l {jp jj K; jr jc l} T^l_prcj / (e_c+e_j-e_p-e_r)
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ip = 0; ip < row_particles.size(); ++ip) {
    const auto &p = row_particles[ip];
    for (std::size_t ij = 0; ij < row_holes.size(); ++ij) {
      const auto &j = row_holes[ij];
      const auto parity = pair_parity(p, j);
      for (std::size_t ir = 0; ir < excited.size(); ++ir) {
        const auto &r = excited[ir];
        for (std::size_t ic = 0; ic < core.size(); ++ic) {
          const auto &c = core[ic];
          if (pair_parity(r, c) != parity)
            continue;
          const auto T = rung_multipoles(p, r, c, j, qk, lk_prev, lk_fixed);
          if (T.k_max < T.k_min)
            continue;
          const auto inv_de = 1.0 / (c.en + j.en - p.en - r.en);
          const auto K_min =
            std::max(std::abs(p.twoj - j.twoj), std::abs(r.twoj - c.twoj)) / 2;
          const auto K_max =
            std::min({p.twoj + j.twoj, r.twoj + c.twoj, 2 * max_K}) / 2;
          for (int K = K_min; K <= K_max; ++K) {
            auto &B = B_pjrc[OrbitalPairIndex::block(K, parity)];
            if (B.rows() == 0)
              continue;
            double BK = 0.0;
            for (int l = T.k_min; l <= T.k_max; ++l) {
              BK += sixj.get(p.twoj, j.twoj, K, r.twoj, c.twoj, l) *
                    T.value[std::size_t(l)];
            }
            const auto row = std::size_t(row_pairs.index(K, ip, ij));
            const auto col = std::size_t(rc_pairs.index(K, ir, ic));
            B(row, col) = (2 * K + 1) * BK * inv_de;
          }
        }
      }
    }
  }

  // One matrix product for each block (BLAS is multi-threaded here)
  for (int K = 0; K <= max_K; ++K) {
    for (int parity = 0; parity < 2; ++parity) {
      const auto block = OrbitalPairIndex::block(K, parity);
      if (Y_pjqi[block].rows() == 0)
        continue;
      LinAlg::GEMM(B_pjrc[block], A_rcqi[block], &Y_pjqi[block]);
      Y_pjqi[block] *= double(Angular::neg1pow(K));
    }
  }
  return Y_pjqi;
}

//==============================================================================
LadderChange calculate_Lk(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                          const std::vector<LadderOrbital> &excited,
                          const std::vector<LadderOrbital> &core,
                          const LadderIntegralSet &integrals,
                          const LadderDiagrams &diagrams,
                          const Coulomb::LkTable *lk_prev,
                          const Coulomb::LkTable *lk_fixed, int max_k,
                          bool only_new, double a_damp) {
  assert(lk != lk_prev && lk != lk_fixed);
  const auto &m_orbitals = integrals.m_orbitals;
  const auto &ij_orbitals = integrals.ij_orbitals;
  const auto symmetric_mn = integrals.symmetric_mn;
  if (m_orbitals.empty() || excited.empty() || core.empty() ||
      ij_orbitals.empty())
    return {};

  const auto num_m = m_orbitals.size();
  const auto num_excited = excited.size();
  const auto num_core = core.size();
  const auto num_ij = ij_orbitals.size();
  assert(integrals.in_i_slot.size() == num_ij);
  assert(integrals.in_j_slot.size() == num_ij);
  assert(!symmetric_mn || num_m == num_excited);
  std::vector<char> required_ij(num_ij * num_ij, 0);
  for (std::size_t ii = 0; ii < num_ij; ++ii) {
    for (std::size_t ij = 0; ij < num_ij; ++ij) {
      required_ij[ii * num_ij + ij] = integrals.required(ii, ij);
    }
  }

  // Largest angular momenta: J of an (i,j) pair, K of a particle-hole pair
  const auto max_tj_particle =
    std::max(max_twoj(excited), max_twoj(m_orbitals));
  const auto max_tj_hole = std::max(max_twoj(ij_orbitals), max_twoj(core));
  const auto max_tj = std::max(max_tj_particle, max_tj_hole);
  const int max_J = max_twoj(ij_orbitals);
  const int max_K = (max_tj_particle + max_tj_hole) / 2;
  // every multipolarity and rank is at most (2*ja + 2*jb)/2 <= max_tj
  assert(max_tj < int(ladder_max_multipoles));
  const SixJArray sixj(max_tj, max_tj);

  // Damping: same convention as CoulombTable::update
  const auto damp_old = std::clamp(a_damp, 0.0, 1.0);
  const auto damp_new = 1.0 - a_damp;

  using Matrix = LinAlg::Matrix<double>;

  //----------------------------------------------------------------------------
  // Particle-particle rung: X^J_{rs,ij} = T^J_{rs,ij} / (e_i + e_j - e_r - e_s)
  // rows: ordered excited pairs (r,s); columns: required pairs (i,j)
  const OrbitalPairIndex rs_pairs(excited, excited, max_J);
  const OrbitalPairIndex ij_pairs(ij_orbitals, ij_orbitals, max_J, required_ij);
  std::vector<Matrix> X_rsij(rs_pairs.num_blocks());
  for (int J = 0; J <= max_J; ++J) {
    for (int parity = 0; parity < 2; ++parity) {
      const auto rows = rs_pairs.size(J, parity);
      const auto cols = ij_pairs.size(J, parity);
      if (diagrams.L1 && rows != 0 && cols != 0) {
        X_rsij[OrbitalPairIndex::block(J, parity)] = Matrix(rows, cols);
      }
    }
  }
  const auto num_r_for_rung = diagrams.L1 ? num_excited : 0;
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ir = 0; ir < num_r_for_rung; ++ir) {
    const auto &r = excited[ir];
    for (std::size_t is = 0; is < num_excited; ++is) {
      const auto &s = excited[is];
      const auto parity = pair_parity(r, s);
      for (std::size_t ii = 0; ii < num_ij; ++ii) {
        const auto &i = ij_orbitals[ii];
        for (std::size_t ij = 0; ij < num_ij; ++ij) {
          const auto &j = ij_orbitals[ij];
          if (!required_ij[ii * num_ij + ij] || pair_parity(i, j) != parity)
            continue;
          const auto T = rung_multipoles(r, s, i, j, qk, lk_prev, lk_fixed);
          if (T.k_max < T.k_min)
            continue;
          const auto inv_de = 1.0 / (i.en + j.en - r.en - s.en);
          const auto J_min =
            std::max(std::abs(r.twoj - s.twoj), std::abs(i.twoj - j.twoj)) / 2;
          const auto J_max =
            std::min({r.twoj + s.twoj, i.twoj + j.twoj, 2 * max_J}) / 2;
          for (int J = J_min; J <= J_max; ++J) {
            double TJ = 0.0;
            for (int l = T.k_min; l <= T.k_max; ++l) {
              TJ += Angular::neg1pow(l) *
                    sixj.get(r.twoj, s.twoj, J, j.twoj, i.twoj, l) *
                    T.value[std::size_t(l)];
            }
            const auto row = std::size_t(rs_pairs.index(J, ir, is));
            const auto col = std::size_t(ij_pairs.index(J, ii, ij));
            X_rsij[OrbitalPairIndex::block(J, parity)](row, col) =
              pair_sign(J, r, i) * TJ * inv_de;
          }
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // Hole-hole: Q^J_{cd,ij}; rows: ordered core pairs (c,d); columns: (i,j)
  const OrbitalPairIndex cd_pairs(core, core, max_J);
  std::vector<Matrix> Q_cdij(cd_pairs.num_blocks());
  if (diagrams.L4) {
    for (int J = 0; J <= max_J; ++J) {
      for (int parity = 0; parity < 2; ++parity) {
        const auto rows = cd_pairs.size(J, parity);
        const auto cols = ij_pairs.size(J, parity);
        if (rows != 0 && cols != 0) {
          Q_cdij[OrbitalPairIndex::block(J, parity)] = Matrix(rows, cols);
        }
      }
    }
#pragma omp parallel for schedule(dynamic)
    for (std::size_t ic = 0; ic < num_core; ++ic) {
      const auto &c = core[ic];
      for (std::size_t id = 0; id < num_core; ++id) {
        const auto &d = core[id];
        const auto parity = pair_parity(c, d);
        for (std::size_t ii = 0; ii < num_ij; ++ii) {
          const auto &i = ij_orbitals[ii];
          for (std::size_t ij = 0; ij < num_ij; ++ij) {
            const auto &j = ij_orbitals[ij];
            if (!required_ij[ii * num_ij + ij] || pair_parity(i, j) != parity)
              continue;
            const auto Q = coulomb_multipoles(c, d, i, j, qk);
            if (Q.k_max < Q.k_min)
              continue;
            const auto J_min =
              std::max(std::abs(c.twoj - d.twoj), std::abs(i.twoj - j.twoj)) /
              2;
            const auto J_max =
              std::min({c.twoj + d.twoj, i.twoj + j.twoj, 2 * max_J}) / 2;
            for (int J = J_min; J <= J_max; ++J) {
              double QJ = 0.0;
              for (int u = Q.k_min; u <= Q.k_max; u += 2) {
                QJ += Angular::neg1pow(u) *
                      sixj.get(c.twoj, d.twoj, J, j.twoj, i.twoj, u) *
                      Q.value[std::size_t(u)];
              }
              const auto row = std::size_t(cd_pairs.index(J, ic, id));
              const auto col = std::size_t(ij_pairs.index(J, ii, ij));
              Q_cdij[OrbitalPairIndex::block(J, parity)](row, col) =
                pair_sign(J, c, i) * QJ;
            }
          }
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // Particle-hole: L2^k_mnij needs Y^K_{mj,ni}; L3^k_mnij needs Y^K_{ni,mj}.
  // Only two rectangles of Y are needed, since i is an orbital of the i slot
  // and j one of the j slot:
  //   Y_mj_ni: rows (m, j in j-slot), columns (n, i in i-slot)
  //   Y_ni_mj: rows (n, i in i-slot), columns (m, j in j-slot)
  // [For the partner integrals (i in the j-slot set, j in the i-slot set; if
  // m runs over the excited set) the same two serve, with the roles swapped.]
  std::vector<LadderOrbital> i_slot_orbitals, j_slot_orbitals;
  // position of each ij orbital in those lists (-1: not a member)
  std::vector<int> i_slot_index(num_ij, -1), j_slot_index(num_ij, -1);
  for (std::size_t ih = 0; ih < num_ij; ++ih) {
    if (integrals.in_i_slot[ih]) {
      i_slot_index[ih] = int(i_slot_orbitals.size());
      i_slot_orbitals.push_back(ij_orbitals[ih]);
    }
    if (integrals.in_j_slot[ih]) {
      j_slot_index[ih] = int(j_slot_orbitals.size());
      j_slot_orbitals.push_back(ij_orbitals[ih]);
    }
  }
  const OrbitalPairIndex mj_pairs(m_orbitals, j_slot_orbitals, max_K);
  const OrbitalPairIndex ni_pairs(excited, i_slot_orbitals, max_K);
  // (empty matrices if the diagram is not included)
  const auto Y_mj_ni =
    diagrams.L2 ?
      ladder_particle_hole(m_orbitals, j_slot_orbitals, mj_pairs, excited,
                           i_slot_orbitals, ni_pairs, excited, core, qk,
                           lk_prev, lk_fixed, sixj) :
      std::vector<Matrix>(mj_pairs.num_blocks());
  const auto Y_ni_mj =
    diagrams.L3 ?
      ladder_particle_hole(excited, i_slot_orbitals, ni_pairs, m_orbitals,
                           j_slot_orbitals, mj_pairs, excited, core, qk,
                           lk_prev, lk_fixed, sixj) :
      std::vector<Matrix>(ni_pairs.num_blocks());

  //----------------------------------------------------------------------------
  // Main loop: for each m, and each group of excited n with the same kappa
  // (if m runs over the excited set: only n >= m; the other half follows from
  // L^k_mnij = L^k_nmji):
  // form Q^J_{mn,rs} for all (r,s) from the Qk table; then
  //   L1^J_{mn,ij} = sum_{rs} Q^J_{mn,rs} X^J_{rs,ij}
  // and add L4, L2, L3; transform back to L^k_mnij, and store.
  // The n of one group share all angular factors and the block structure, so
  // form one matrix (rows n) for each J.
  std::vector<std::vector<std::size_t>> kappa_groups;
  for (std::size_t in = 0; in < num_excited; ++in) {
    bool found_group = false;
    for (auto &group : kappa_groups) {
      if (excited[group.front()].kappa == excited[in].kappa) {
        group.push_back(in);
        found_group = true;
        break;
      }
    }
    if (!found_group) {
      kappa_groups.push_back({in});
    }
  }
  const auto num_groups = kappa_groups.size();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

  double max_value = 0.0;
  double max_change = 0.0;

#pragma omp parallel for collapse(2) schedule(dynamic)                         \
  reduction(max : max_value, max_change)
  for (std::size_t im = 0; im < num_m; ++im) {
    for (std::size_t ig = 0; ig < num_groups; ++ig) {
      const auto &m = m_orbitals[im];

      // n in this kappa group (with n >= m, if m runs over the excited set)
      std::vector<std::size_t> n_list;
      for (const auto in : kappa_groups[ig]) {
        if (!symmetric_mn || in >= im) {
          n_list.push_back(in);
        }
      }
      if (n_list.empty())
        continue;
      const auto num_n = n_list.size();
      // quantum numbers shared by all n in the group:
      const auto &n_any = excited[n_list.front()];
      const auto parity_mn = pair_parity(m, n_any);
      const auto J_min_mn = std::abs(m.twoj - n_any.twoj) / 2;
      const auto J_max_mn = std::min(m.twoj + n_any.twoj, 2 * max_J) / 2;

      // Pair-coupled part, L1^J + L4^J: for each J, one matrix with rows n
      // and columns (i,j)
      std::vector<Matrix> LJ_mnij(std::size_t(max_J) + 1);

      // L1: Q^J_{mn,rs}, rows n, columns (r,s)
      std::vector<Matrix> QJ_mnrs(std::size_t(max_J) + 1);
      bool any_J = false;
      for (int J = J_min_mn; J <= J_max_mn; ++J) {
        const auto num_rs = rs_pairs.size(J, parity_mn);
        if (diagrams.L1 && num_rs != 0 && ij_pairs.size(J, parity_mn) != 0) {
          QJ_mnrs[std::size_t(J)] = Matrix(num_n, num_rs);
          any_J = true;
        }
      }
      if (any_J) {
        for (std::size_t ir = 0; ir < num_excited; ++ir) {
          const auto &r = excited[ir];
          for (std::size_t is = 0; is < num_excited; ++is) {
            const auto &s = excited[is];
            if (pair_parity(r, s) != parity_mn)
              continue;
            const auto J_min =
              std::max(J_min_mn, std::abs(r.twoj - s.twoj) / 2);
            const auto J_max = std::min(J_max_mn, (r.twoj + s.twoj) / 2);
            if (J_max < J_min)
              continue;
            for (std::size_t row = 0; row < num_n; ++row) {
              const auto &n = excited[n_list[row]];
              const auto Q = coulomb_multipoles(m, n, r, s, qk);
              // k range depends on kappas only: same for all n in the group
              if (Q.k_max < Q.k_min)
                break;
              for (int J = J_min; J <= J_max; ++J) {
                auto &QJ = QJ_mnrs[std::size_t(J)];
                if (QJ.rows() == 0)
                  continue;
                double sum = 0.0;
                for (int u = Q.k_min; u <= Q.k_max; u += 2) {
                  sum += Angular::neg1pow(u) *
                         sixj.get(m.twoj, n.twoj, J, s.twoj, r.twoj, u) *
                         Q.value[std::size_t(u)];
                }
                const auto col = std::size_t(rs_pairs.index(J, ir, is));
                QJ(row, col) = pair_sign(J, m, r) * sum;
              }
            }
          }
        }
        for (int J = J_min_mn; J <= J_max_mn; ++J) {
          const auto &QJ = QJ_mnrs[std::size_t(J)];
          if (QJ.rows() == 0)
            continue;
          const auto &X = X_rsij[OrbitalPairIndex::block(J, parity_mn)];
          LJ_mnij[std::size_t(J)] = Matrix(num_n, X.cols());
          LinAlg::GEMM(QJ, X, &LJ_mnij[std::size_t(J)]);
        }
      }

      // L4: sum_{cd} [T^J_{mn,cd} / (e_c + e_d - e_m - e_n)] Q^J_{cd,ij}
      if (diagrams.L4) {
        for (int J = J_min_mn; J <= J_max_mn; ++J) {
          const auto &QJ_cdij = Q_cdij[OrbitalPairIndex::block(J, parity_mn)];
          if (QJ_cdij.rows() == 0)
            continue;
          Matrix X_mncd(num_n, QJ_cdij.rows());
          for (std::size_t row = 0; row < num_n; ++row) {
            const auto &n = excited[n_list[row]];
            for (std::size_t ic = 0; ic < num_core; ++ic) {
              const auto &c = core[ic];
              for (std::size_t id = 0; id < num_core; ++id) {
                const auto &d = core[id];
                const auto col = cd_pairs.index(J, ic, id);
                if (col < 0 || pair_parity(c, d) != parity_mn)
                  continue;
                const auto T =
                  rung_multipoles(m, n, c, d, qk, lk_prev, lk_fixed);
                double TJ = 0.0;
                for (int l = T.k_min; l <= T.k_max; ++l) {
                  TJ += Angular::neg1pow(l) *
                        sixj.get(m.twoj, n.twoj, J, d.twoj, c.twoj, l) *
                        T.value[std::size_t(l)];
                }
                X_mncd(row, std::size_t(col)) =
                  pair_sign(J, m, c) * TJ / (c.en + d.en - m.en - n.en);
              }
            }
          }
          Matrix L4J_mnij(num_n, QJ_cdij.cols());
          LinAlg::GEMM(X_mncd, QJ_cdij, &L4J_mnij);
          if (LJ_mnij[std::size_t(J)].rows() == 0) {
            LJ_mnij[std::size_t(J)] = std::move(L4J_mnij);
          } else {
            LJ_mnij[std::size_t(J)] += L4J_mnij;
          }
        }
      }

      // Transform back to the multipole form, add L2 + L3, and store
      std::array<double, ladder_max_multipoles> LJ{};
      std::array<double, ladder_max_multipoles> YK_mjni{};
      std::array<double, ladder_max_multipoles> YK_nimj{};
      for (std::size_t row = 0; row < num_n; ++row) {
        const auto in = n_list[row];
        const auto &n = excited[in];
        for (std::size_t ii = 0; ii < num_ij; ++ii) {
          const auto &i = ij_orbitals[ii];
          for (std::size_t ij = 0; ij < num_ij; ++ij) {
            const auto &j = ij_orbitals[ij];
            if (!required_ij[ii * num_ij + ij] ||
                pair_parity(i, j) != parity_mn)
              continue;
            // For m = n the two orders (i,j) and (j,i) are the same integral
            if (symmetric_mn && im == in &&
                !lk->is_NormalOrdered(m.label, n.label, i.label, j.label))
              continue;
            const auto [k_min, k_max_L] =
              k_minmax_L(m.kappa, n.kappa, i.kappa, j.kappa);
            const auto k_max = std::min(k_max_L, max_k);
            if (k_max < k_min)
              continue;

            // [J] (-1)^(J+jm+ji) (L1 + L4)^J_{mn,ij}
            const auto J_min =
              std::max(J_min_mn, std::abs(i.twoj - j.twoj) / 2);
            const auto J_max = std::min(J_max_mn, (i.twoj + j.twoj) / 2);
            for (int J = J_min; J <= J_max; ++J) {
              const auto &LJ_block = LJ_mnij[std::size_t(J)];
              LJ[std::size_t(J)] =
                LJ_block.rows() == 0 ?
                  0.0 :
                  (2 * J + 1) * pair_sign(J, m, i) *
                    LJ_block(row, std::size_t(ij_pairs.index(J, ii, ij)));
            }

            // Y^K_{mj,ni} (for L2) and Y^K_{ni,mj} (for L3)
            const auto parity_mj = pair_parity(m, j);
            const auto K_min =
              std::max(std::abs(m.twoj - j.twoj), std::abs(n.twoj - i.twoj)) /
              2;
            const auto K_max =
              std::min({m.twoj + j.twoj, n.twoj + i.twoj, 2 * max_K}) / 2;
            // i in the i-slot set and j in the j-slot set: as above.
            // Otherwise (partner integral; m and n both excited) the two
            // rectangles swap roles: (m,j) is a row of Y_ni_mj, etc.
            const bool direct =
              integrals.in_i_slot[ii] && integrals.in_j_slot[ij];
            const auto &Y_for_L2 = direct ? Y_mj_ni : Y_ni_mj;
            const auto &Y_for_L3 = direct ? Y_ni_mj : Y_mj_ni;
            const auto &mj_index = direct ? mj_pairs : ni_pairs;
            const auto &ni_index = direct ? ni_pairs : mj_pairs;
            const auto j_position =
              std::size_t(direct ? j_slot_index[ij] : i_slot_index[ij]);
            const auto i_position =
              std::size_t(direct ? i_slot_index[ii] : j_slot_index[ii]);
            for (int K = K_min; K <= K_max; ++K) {
              const auto block = OrbitalPairIndex::block(K, parity_mj);
              const auto mj = std::size_t(mj_index.index(K, im, j_position));
              const auto ni = std::size_t(ni_index.index(K, in, i_position));
              YK_mjni[std::size_t(K)] =
                Y_for_L2[block].rows() == 0 ? 0.0 : Y_for_L2[block](mj, ni);
              YK_nimj[std::size_t(K)] =
                Y_for_L3[block].rows() == 0 ? 0.0 : Y_for_L3[block](ni, mj);
            }

            const auto key =
              lk->NormalOrder(m.label, n.label, i.label, j.label);
            for (int k = k_min; k <= k_max; ++k) {
              double *const L_stored = lk->get(k, key);
              if (L_stored == nullptr)
                continue;
              if (only_new && *L_stored != 0.0)
                continue;
              // L1 + L4:
              double Lk_mnij = 0.0;
              for (int J = J_min; J <= J_max; ++J) {
                Lk_mnij += sixj.get(m.twoj, n.twoj, J, j.twoj, i.twoj, k) *
                           LJ[std::size_t(J)];
              }
              Lk_mnij *= Angular::neg1pow(k);
              // L2 + L3:
              for (int K = K_min; K <= K_max; ++K) {
                Lk_mnij += sixj.get(m.twoj, j.twoj, K, n.twoj, i.twoj, k) *
                           YK_mjni[std::size_t(K)];
                Lk_mnij += sixj.get(n.twoj, i.twoj, K, m.twoj, j.twoj, k) *
                           YK_nimj[std::size_t(K)];
              }
              Lk_mnij *= (2 * k + 1);

              const auto L_new =
                (only_new || *L_stored == 0.0 || damp_old == 0.0) ?
                  Lk_mnij :
                  damp_new * Lk_mnij + damp_old * *L_stored;
              max_change = std::max(max_change, std::abs(L_new - *L_stored));
              max_value = std::max(max_value, std::abs(L_new));
              *L_stored = L_new;
            }
          }
        }
      }
    }
  }
  return {max_value, max_change};
}

} // namespace MBPT
