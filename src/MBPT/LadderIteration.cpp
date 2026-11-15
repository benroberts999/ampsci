#include "MBPT/LadderIteration.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "LinAlg/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
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

using Matrix = LinAlg::Matrix<double>;

//==============================================================================
std::pair<int, int> k_minmax_L(int kappa_a, int kappa_b, int kappa_c,
                               int kappa_d) {
  // Each Coulomb vertex shifts the parity of both pair lines together, so
  // only the total parity is constrained
  if ((Angular::l_k(kappa_a) + Angular::l_k(kappa_b) + Angular::l_k(kappa_c) +
       Angular::l_k(kappa_d)) %
        2 !=
      0) {
    return {1, 0};
  }
  const auto [l1, u1] =
    Coulomb::k_minmax_tj(Angular::twoj_k(kappa_a), Angular::twoj_k(kappa_c));
  const auto [l2, u2] =
    Coulomb::k_minmax_tj(Angular::twoj_k(kappa_b), Angular::twoj_k(kappa_d));
  return {std::max(l1, l2), std::min({u1, u2, ladder_max_k})};
}

std::pair<int, int> k_minmax_L(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d) {
  return k_minmax_L(a.kappa(), b.kappa(), c.kappa(), d.kappa());
}

//==============================================================================
bool ExternalLine::is_ij_line(const DiracSpinor &a) const {
  return L_mnxb != nullptr && std::find(x_orbitals.cbegin(), x_orbitals.cend(),
                                        a) != x_orbitals.cend();
}

bool ExternalLine::is_m_line(const DiracSpinor &a) const {
  return L_xnab != nullptr && std::find(x_orbitals.cbegin(), x_orbitals.cend(),
                                        a) != x_orbitals.cend();
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
OrbitalPairIndex::OrbitalPairIndex(const std::vector<DiracSpinor> &first,
                                   const std::vector<DiracSpinor> &second,
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
      const auto parity = pair_parity(first[ix], second[iy]);
      const auto [rank_min, rank_max] =
        Coulomb::k_minmax_tj(first[ix].twoj(), second[iy].twoj());
      for (int rank = rank_min; rank <= std::min(rank_max, max_rank); ++rank) {
        m_index[std::size_t(rank)][ix * m_size_second + iy] =
          int(m_count[block(rank, parity)]);
        ++m_count[block(rank, parity)];
      }
    }
  }
}

//==============================================================================
Multipoles coulomb_multipoles(const DiracSpinor &a, const DiracSpinor &b,
                              const DiracSpinor &c, const DiracSpinor &d,
                              const Coulomb::QkTable &qk) {
  Multipoles Q;
  Q.k_step = 2;
  const auto [k_min, k_max] =
    Coulomb::k_minmax_Q(a.kappa(), b.kappa(), c.kappa(), d.kappa());
  Q.k_min = k_min;
  Q.k_max = std::min(k_max, ladder_max_k);
  if (Q.k_max < Q.k_min)
    return Q;
  const auto key =
    qk.NormalOrder(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
  for (int k = Q.k_min; k <= Q.k_max; k += Q.k_step) {
    Q.value[std::size_t(k)] = qk.Q(k, key);
  }
  return Q;
}

//==============================================================================
template <typename LadderTable>
Multipoles QL_multipoles(const DiracSpinor &a, const DiracSpinor &b,
                         const DiracSpinor &c, const DiracSpinor &d,
                         const Coulomb::QkTable &qk, const LadderTable *lk) {
  Multipoles QL;
  std::tie(QL.k_min, QL.k_max) = k_minmax_L(a, b, c, d);
  if (QL.k_max < QL.k_min)
    return QL;
  const auto [q_min, q_max] =
    Coulomb::k_minmax_Q(a.kappa(), b.kappa(), c.kappa(), d.kappa());
  const auto Q_key =
    qk.NormalOrder(a.nk_index(), b.nk_index(), c.nk_index(), d.nk_index());
  const auto L_key = lk ? lk->NormalOrder(a.nk_index(), b.nk_index(),
                                          c.nk_index(), d.nk_index()) :
                          0ul;
  for (int k = QL.k_min; k <= QL.k_max; ++k) {
    const bool has_Q = k >= q_min && k <= q_max && (k - q_min) % 2 == 0;
    QL.value[std::size_t(k)] =
      (has_Q ? qk.Q(k, Q_key) : 0.0) + (lk ? lk->Q(k, L_key) : 0.0);
  }
  return QL;
}

template Multipoles QL_multipoles(const DiracSpinor &, const DiracSpinor &,
                                  const DiracSpinor &, const DiracSpinor &,
                                  const Coulomb::QkTable &,
                                  const Coulomb::LkTable *);
template Multipoles QL_multipoles(const DiracSpinor &, const DiracSpinor &,
                                  const DiracSpinor &, const DiracSpinor &,
                                  const Coulomb::QkTable &,
                                  const Coulomb::NkTable *);

//==============================================================================
double pair_coupled(int J, const DiracSpinor &a, const DiracSpinor &b,
                    const DiracSpinor &c, const DiracSpinor &d,
                    const Multipoles &L_k, const SixJArray &sixj) {
  double sum = 0.0;
  for (int k = L_k.k_min; k <= L_k.k_max; k += L_k.k_step) {
    sum += Angular::neg1pow(k) *
           sixj.get(a.twoj(), b.twoj(), J, d.twoj(), c.twoj(), k) *
           L_k.value[std::size_t(k)];
  }
  return pair_sign(J, a, c) * sum;
}

//==============================================================================
double multipole_from_pair_coupled(
  int k, const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
  const DiracSpinor &d, const std::array<double, ladder_max_multipoles> &L_J,
  int J_min, int J_max, const SixJArray &sixj) {
  double sum = 0.0;
  for (int J = J_min; J <= J_max; ++J) {
    sum += (2 * J + 1) * pair_sign(J, a, c) *
           sixj.get(a.twoj(), b.twoj(), J, d.twoj(), c.twoj(), k) *
           L_J[std::size_t(J)];
  }
  return Angular::neg1pow(k) * (2 * k + 1) * sum;
}

//==============================================================================
double particle_hole_coupled(int K, const DiracSpinor &a, const DiracSpinor &b,
                             const DiracSpinor &c, const DiracSpinor &d,
                             const Multipoles &L_k, const SixJArray &sixj) {
  double sum = 0.0;
  for (int k = L_k.k_min; k <= L_k.k_max; k += L_k.k_step) {
    sum += sixj.get(a.twoj(), d.twoj(), K, b.twoj(), c.twoj(), k) *
           L_k.value[std::size_t(k)];
  }
  return sum;
}

//==============================================================================
double multipole_from_particle_hole_coupled(
  int k, const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
  const DiracSpinor &d, const std::array<double, ladder_max_multipoles> &L_K,
  int K_min, int K_max, const SixJArray &sixj) {
  double sum = 0.0;
  for (int K = K_min; K <= K_max; ++K) {
    sum += (2 * K + 1) *
           sixj.get(a.twoj(), d.twoj(), K, b.twoj(), c.twoj(), k) *
           L_K[std::size_t(K)];
  }
  return (2 * k + 1) * sum;
}

//==============================================================================
std::vector<DiracSpinor>
ladder_ij_orbitals(const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const std::vector<DiracSpinor> &i_orbs) {
  auto ij_orbitals = core;
  for (const auto &Fi : i_orbs) {
    if (std::find(core.cbegin(), core.cend(), Fi) != core.cend())
      continue;
    const auto in_excited = std::find(excited.cbegin(), excited.cend(), Fi);
    if (in_excited != excited.cend() &&
        std::find(ij_orbitals.cbegin(), ij_orbitals.cend(), Fi) ==
          ij_orbitals.cend()) {
      ij_orbitals.push_back(*in_excited);
    }
  }
  return ij_orbitals;
}

//==============================================================================
std::vector<char> required_ij_pairs(const std::vector<DiracSpinor> &ij_orbitals,
                                    const std::vector<DiracSpinor> &core,
                                    const std::vector<DiracSpinor> &i_orbs) {
  const auto num_ij = ij_orbitals.size();
  std::vector<char> is_core(num_ij, 0);
  std::vector<char> is_i_orbital(num_ij, 0);
  for (std::size_t index = 0; index < num_ij; ++index) {
    const auto &Fij = ij_orbitals[index];
    is_core[index] = std::find(core.cbegin(), core.cend(), Fij) != core.cend();
    is_i_orbital[index] =
      std::find(i_orbs.cbegin(), i_orbs.cend(), Fij) != i_orbs.cend();
  }
  std::vector<char> required(num_ij * num_ij, 0);
  for (std::size_t index_i = 0; index_i < num_ij; ++index_i) {
    for (std::size_t index_j = 0; index_j < num_ij; ++index_j) {
      required[index_i * num_ij + index_j] =
        (is_i_orbital[index_i] && is_core[index_j]) ||
        (is_core[index_i] && is_i_orbital[index_j]);
    }
  }
  return required;
}

//==============================================================================
std::vector<std::vector<std::size_t>>
group_by_kappa(const std::vector<DiracSpinor> &orbitals) {
  std::vector<std::vector<std::size_t>> groups;
  for (std::size_t index = 0; index < orbitals.size(); ++index) {
    const auto group =
      std::find_if(groups.begin(), groups.end(), [&](const auto &g) {
        return orbitals[g.front()].kappa() == orbitals[index].kappa();
      });
    if (group == groups.end()) {
      groups.push_back({index});
    } else {
      group->push_back(index);
    }
  }
  return groups;
}

//==============================================================================
std::vector<Matrix> block_matrices(const OrbitalPairIndex &row_pairs,
                                   const OrbitalPairIndex &column_pairs) {
  const auto max_rank = std::min(row_pairs.max_rank(), column_pairs.max_rank());
  std::vector<Matrix> blocks(row_pairs.num_blocks());
  for (int rank = 0; rank <= max_rank; ++rank) {
    for (int parity = 0; parity < 2; ++parity) {
      const auto rows = row_pairs.size(rank, parity);
      const auto cols = column_pairs.size(rank, parity);
      if (rows != 0 && cols != 0) {
        blocks[OrbitalPairIndex::block(rank, parity)] = Matrix(rows, cols);
      }
    }
  }
  return blocks;
}

//==============================================================================
std::vector<Matrix> QL_particle_particle_over_de(
  const std::vector<DiracSpinor> &excited,
  const std::vector<DiracSpinor> &ij_orbitals,
  const std::vector<char> &required_ij, const OrbitalPairIndex &rs_pairs,
  const OrbitalPairIndex &ij_pairs, const Coulomb::QkTable &qk,
  const Coulomb::LkTable *lk_prev, const SixJArray &sixj,
  const ExternalLine *external) {
  const auto max_J = std::min(rs_pairs.max_rank(), ij_pairs.max_rank());
  const auto num_ij = ij_orbitals.size();
  auto QLJ_rs_ij_over_de = block_matrices(rs_pairs, ij_pairs);

  // The external line carries en_v, and its integrals are in L_mnxb
  std::vector<double> en_ij(num_ij);
  std::vector<char> is_line(num_ij, 0);
  for (std::size_t index = 0; index < num_ij; ++index) {
    is_line[index] = external && external->is_ij_line(ij_orbitals[index]);
    en_ij[index] = is_line[index] ? external->en_v : ij_orbitals[index].en();
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t index_r = 0; index_r < excited.size(); ++index_r) {
    const auto &r = excited[index_r];
    for (std::size_t index_s = 0; index_s < excited.size(); ++index_s) {
      const auto &s = excited[index_s];
      const auto parity = pair_parity(r, s);
      for (std::size_t index_i = 0; index_i < num_ij; ++index_i) {
        const auto &i = ij_orbitals[index_i];
        for (std::size_t index_j = 0; index_j < num_ij; ++index_j) {
          const auto &j = ij_orbitals[index_j];
          if (!required_ij[index_i * num_ij + index_j] ||
              pair_parity(i, j) != parity)
            continue;
          const auto QL = is_line[index_i] || is_line[index_j] ?
                            QL_multipoles(r, s, i, j, qk, external->L_mnxb) :
                            QL_multipoles(r, s, i, j, qk, lk_prev);
          if (QL.k_max < QL.k_min)
            continue;
          const auto de = en_ij[index_i] + en_ij[index_j] - r.en() - s.en();
          const auto [J_min, J_max] = coupled_rank_range(r, s, i, j, max_J);
          for (int J = J_min; J <= J_max; ++J) {
            const auto row = std::size_t(rs_pairs.index(J, index_r, index_s));
            const auto col = std::size_t(ij_pairs.index(J, index_i, index_j));
            QLJ_rs_ij_over_de[OrbitalPairIndex::block(J, parity)](row, col) =
              pair_coupled(J, r, s, i, j, QL, sixj) / de;
          }
        }
      }
    }
  }
  return QLJ_rs_ij_over_de;
}

//==============================================================================
std::vector<Matrix> Q_hole_hole(const std::vector<DiracSpinor> &core,
                                const std::vector<DiracSpinor> &ij_orbitals,
                                const std::vector<char> &required_ij,
                                const OrbitalPairIndex &cd_pairs,
                                const OrbitalPairIndex &ij_pairs,
                                const Coulomb::QkTable &qk,
                                const SixJArray &sixj) {
  const auto max_J = std::min(cd_pairs.max_rank(), ij_pairs.max_rank());
  const auto num_ij = ij_orbitals.size();
  auto QJ_cd_ij = block_matrices(cd_pairs, ij_pairs);

#pragma omp parallel for schedule(dynamic)
  for (std::size_t index_c = 0; index_c < core.size(); ++index_c) {
    const auto &c = core[index_c];
    for (std::size_t index_d = 0; index_d < core.size(); ++index_d) {
      const auto &d = core[index_d];
      const auto parity = pair_parity(c, d);
      for (std::size_t index_i = 0; index_i < num_ij; ++index_i) {
        const auto &i = ij_orbitals[index_i];
        for (std::size_t index_j = 0; index_j < num_ij; ++index_j) {
          const auto &j = ij_orbitals[index_j];
          if (!required_ij[index_i * num_ij + index_j] ||
              pair_parity(i, j) != parity)
            continue;
          const auto Q = coulomb_multipoles(c, d, i, j, qk);
          if (Q.k_max < Q.k_min)
            continue;
          const auto [J_min, J_max] = coupled_rank_range(c, d, i, j, max_J);
          for (int J = J_min; J <= J_max; ++J) {
            const auto row = std::size_t(cd_pairs.index(J, index_c, index_d));
            const auto col = std::size_t(ij_pairs.index(J, index_i, index_j));
            QJ_cd_ij[OrbitalPairIndex::block(J, parity)](row, col) =
              pair_coupled(J, c, d, i, j, Q, sixj);
          }
        }
      }
    }
  }
  return QJ_cd_ij;
}

//==============================================================================
std::vector<Matrix> Q_particle_particle_rows(
  const DiracSpinor &m, const std::vector<std::size_t> &n_group,
  const std::vector<DiracSpinor> &excited, const OrbitalPairIndex &rs_pairs,
  const OrbitalPairIndex &ij_pairs, const Coulomb::QkTable &qk,
  const SixJArray &sixj) {
  const auto max_J = std::min(rs_pairs.max_rank(), ij_pairs.max_rank());
  const auto num_n = n_group.size();
  std::vector<Matrix> QJ_mn_rs(std::size_t(max_J) + 1);
  if (n_group.empty())
    return QJ_mn_rs;

  // All n of the group share kappa: the same J range, parity and blocks
  const auto &n_any = excited[n_group.front()];
  const auto parity_mn = pair_parity(m, n_any);
  const auto [J_min_mn, J_max_mn] =
    coupled_rank_range(m, n_any, m, n_any, max_J);

  bool any_J = false;
  for (int J = J_min_mn; J <= J_max_mn; ++J) {
    const auto num_rs = rs_pairs.size(J, parity_mn);
    if (num_rs != 0 && ij_pairs.size(J, parity_mn) != 0) {
      QJ_mn_rs[std::size_t(J)] = Matrix(num_n, num_rs);
      any_J = true;
    }
  }
  if (!any_J)
    return QJ_mn_rs;

  for (std::size_t index_r = 0; index_r < excited.size(); ++index_r) {
    const auto &r = excited[index_r];
    for (std::size_t index_s = 0; index_s < excited.size(); ++index_s) {
      const auto &s = excited[index_s];
      if (pair_parity(r, s) != parity_mn)
        continue;
      const auto [J_min, J_max] = coupled_rank_range(m, n_any, r, s, max_J);
      if (J_max < J_min)
        continue;
      for (std::size_t row = 0; row < num_n; ++row) {
        const auto &n = excited[n_group[row]];
        const auto Q = coulomb_multipoles(m, n, r, s, qk);
        // The k range depends on kappa only: the same for all n of the group
        if (Q.k_max < Q.k_min)
          break;
        for (int J = J_min; J <= J_max; ++J) {
          auto &QJ = QJ_mn_rs[std::size_t(J)];
          if (QJ.rows() == 0)
            continue;
          const auto col = std::size_t(rs_pairs.index(J, index_r, index_s));
          QJ(row, col) = pair_coupled(J, m, n, r, s, Q, sixj);
        }
      }
    }
  }
  return QJ_mn_rs;
}

//==============================================================================
std::vector<Matrix> QL_hole_hole_rows_over_de(
  const DiracSpinor &m, const std::vector<std::size_t> &n_group,
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const OrbitalPairIndex &cd_pairs, const OrbitalPairIndex &ij_pairs,
  const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
  const SixJArray &sixj, const ExternalLine *external) {
  const auto max_J = std::min(cd_pairs.max_rank(), ij_pairs.max_rank());
  const auto num_n = n_group.size();
  std::vector<Matrix> QLJ_mn_cd_over_de(std::size_t(max_J) + 1);
  // The external line carries en_v, and its integrals are in L_xnab
  const bool m_is_line = external && external->is_m_line(m);
  const auto en_m = m_is_line ? external->en_v : m.en();
  if (n_group.empty())
    return QLJ_mn_cd_over_de;

  const auto &n_any = excited[n_group.front()];
  const auto parity_mn = pair_parity(m, n_any);
  const auto [J_min_mn, J_max_mn] =
    coupled_rank_range(m, n_any, m, n_any, max_J);

  bool any_J = false;
  for (int J = J_min_mn; J <= J_max_mn; ++J) {
    const auto num_cd = cd_pairs.size(J, parity_mn);
    if (num_cd != 0 && ij_pairs.size(J, parity_mn) != 0) {
      QLJ_mn_cd_over_de[std::size_t(J)] = Matrix(num_n, num_cd);
      any_J = true;
    }
  }
  if (!any_J)
    return QLJ_mn_cd_over_de;

  for (std::size_t row = 0; row < num_n; ++row) {
    const auto &n = excited[n_group[row]];
    for (std::size_t index_c = 0; index_c < core.size(); ++index_c) {
      const auto &c = core[index_c];
      for (std::size_t index_d = 0; index_d < core.size(); ++index_d) {
        const auto &d = core[index_d];
        if (pair_parity(c, d) != parity_mn)
          continue;
        const auto [J_min, J_max] = coupled_rank_range(m, n, c, d, max_J);
        if (J_max < J_min)
          continue;
        const auto QL = m_is_line ?
                          QL_multipoles(m, n, c, d, qk, external->L_xnab) :
                          QL_multipoles(m, n, c, d, qk, lk_prev);
        const auto de = c.en() + d.en() - en_m - n.en();
        for (int J = J_min; J <= J_max; ++J) {
          auto &QLJ = QLJ_mn_cd_over_de[std::size_t(J)];
          if (QLJ.rows() == 0)
            continue;
          const auto col = std::size_t(cd_pairs.index(J, index_c, index_d));
          QLJ(row, col) = pair_coupled(J, m, n, c, d, QL, sixj) / de;
        }
      }
    }
  }
  return QLJ_mn_cd_over_de;
}

//==============================================================================
std::vector<Matrix> QL_particle_hole_over_de(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &ij_orbitals, const OrbitalPairIndex &ph_pairs,
  const OrbitalPairIndex &rc_pairs, const Coulomb::QkTable &qk,
  const Coulomb::LkTable *lk_prev, const SixJArray &sixj,
  const ExternalLine *external) {
  const auto max_K = std::min(ph_pairs.max_rank(), rc_pairs.max_rank());
  auto QLK_mj_rc_over_de = block_matrices(ph_pairs, rc_pairs);

  // The external line carries en_v. As j its integrals are in L_mnxb, as m
  // in L_xnab
  std::vector<char> j_is_line(ij_orbitals.size(), 0);
  for (std::size_t index_j = 0; index_j < ij_orbitals.size(); ++index_j) {
    j_is_line[index_j] = external && external->is_ij_line(ij_orbitals[index_j]);
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t index_m = 0; index_m < excited.size(); ++index_m) {
    const auto &m = excited[index_m];
    const bool m_is_line = external && external->is_m_line(m);
    const auto en_m = m_is_line ? external->en_v : m.en();
    for (std::size_t index_j = 0; index_j < ij_orbitals.size(); ++index_j) {
      const auto &j = ij_orbitals[index_j];
      const auto en_j = j_is_line[index_j] ? external->en_v : j.en();
      const auto parity = pair_parity(m, j);
      for (std::size_t index_r = 0; index_r < excited.size(); ++index_r) {
        const auto &r = excited[index_r];
        for (std::size_t index_c = 0; index_c < core.size(); ++index_c) {
          const auto &c = core[index_c];
          if (pair_parity(r, c) != parity)
            continue;
          const auto QL = j_is_line[index_j] ?
                            QL_multipoles(m, r, c, j, qk, external->L_mnxb) :
                          m_is_line ?
                            QL_multipoles(m, r, c, j, qk, external->L_xnab) :
                            QL_multipoles(m, r, c, j, qk, lk_prev);
          if (QL.k_max < QL.k_min)
            continue;
          const auto de = c.en() + en_j - en_m - r.en();
          const auto [K_min, K_max] = coupled_rank_range(m, j, r, c, max_K);
          for (int K = K_min; K <= K_max; ++K) {
            const auto row = std::size_t(ph_pairs.index(K, index_m, index_j));
            const auto col = std::size_t(rc_pairs.index(K, index_r, index_c));
            QLK_mj_rc_over_de[OrbitalPairIndex::block(K, parity)](row, col) =
              particle_hole_coupled(K, m, r, c, j, QL, sixj) / de;
          }
        }
      }
    }
  }
  return QLK_mj_rc_over_de;
}

//==============================================================================
std::vector<Matrix> Q_particle_hole(const std::vector<DiracSpinor> &excited,
                                    const std::vector<DiracSpinor> &core,
                                    const std::vector<DiracSpinor> &n_orbitals,
                                    const std::vector<DiracSpinor> &i_orbitals,
                                    const OrbitalPairIndex &rc_pairs,
                                    const OrbitalPairIndex &ni_pairs,
                                    const Coulomb::QkTable &qk,
                                    const SixJArray &sixj) {
  const auto max_K = std::min(ni_pairs.max_rank(), rc_pairs.max_rank());
  auto QK_cr_ni = block_matrices(rc_pairs, ni_pairs);

#pragma omp parallel for schedule(dynamic)
  for (std::size_t index_r = 0; index_r < excited.size(); ++index_r) {
    const auto &r = excited[index_r];
    for (std::size_t index_c = 0; index_c < core.size(); ++index_c) {
      const auto &c = core[index_c];
      const auto parity = pair_parity(r, c);
      for (std::size_t index_n = 0; index_n < n_orbitals.size(); ++index_n) {
        const auto &n = n_orbitals[index_n];
        for (std::size_t index_i = 0; index_i < i_orbitals.size(); ++index_i) {
          const auto &i = i_orbitals[index_i];
          if (pair_parity(n, i) != parity)
            continue;
          const auto Q = coulomb_multipoles(c, n, i, r, qk);
          if (Q.k_max < Q.k_min)
            continue;
          const auto [K_min, K_max] = coupled_rank_range(r, c, n, i, max_K);
          for (int K = K_min; K <= K_max; ++K) {
            const auto row = std::size_t(rc_pairs.index(K, index_r, index_c));
            const auto col = std::size_t(ni_pairs.index(K, index_n, index_i));
            QK_cr_ni[OrbitalPairIndex::block(K, parity)](row, col) =
              particle_hole_coupled(K, c, n, i, r, Q, sixj);
          }
        }
      }
    }
  }
  return QK_cr_ni;
}

//==============================================================================
std::vector<Matrix> L1_particle_particle(
  const DiracSpinor &m, const std::vector<std::size_t> &n_group,
  const std::vector<DiracSpinor> &excited, const OrbitalPairIndex &rs_pairs,
  const OrbitalPairIndex &ij_pairs,
  const std::vector<Matrix> &QLJ_rs_ij_over_de, const Coulomb::QkTable &qk,
  const SixJArray &sixj) {
  const auto QJ_mn_rs =
    Q_particle_particle_rows(m, n_group, excited, rs_pairs, ij_pairs, qk, sixj);
  std::vector<Matrix> L1J_mn_ij(QJ_mn_rs.size());
  if (n_group.empty())
    return L1J_mn_ij;
  const auto parity_mn = pair_parity(m, excited[n_group.front()]);
  for (std::size_t J = 0; J < QJ_mn_rs.size(); ++J) {
    const auto &left = QJ_mn_rs[J];
    if (left.rows() == 0)
      continue;
    const auto &right =
      QLJ_rs_ij_over_de[OrbitalPairIndex::block(int(J), parity_mn)];
    L1J_mn_ij[J] = Matrix(left.rows(), right.cols());
    LinAlg::GEMM(left, right, &L1J_mn_ij[J]);
  }
  return L1J_mn_ij;
}

//==============================================================================
std::vector<Matrix>
L4_hole_hole(const DiracSpinor &m, const std::vector<std::size_t> &n_group,
             const std::vector<DiracSpinor> &excited,
             const std::vector<DiracSpinor> &core,
             const OrbitalPairIndex &cd_pairs, const OrbitalPairIndex &ij_pairs,
             const std::vector<Matrix> &QJ_cd_ij, const Coulomb::QkTable &qk,
             const Coulomb::LkTable *lk_prev, const SixJArray &sixj,
             const ExternalLine *external) {
  const auto QLJ_mn_cd_over_de = QL_hole_hole_rows_over_de(
    m, n_group, excited, core, cd_pairs, ij_pairs, qk, lk_prev, sixj, external);
  std::vector<Matrix> L4J_mn_ij(QLJ_mn_cd_over_de.size());
  if (n_group.empty())
    return L4J_mn_ij;
  const auto parity_mn = pair_parity(m, excited[n_group.front()]);
  for (std::size_t J = 0; J < QLJ_mn_cd_over_de.size(); ++J) {
    const auto &left = QLJ_mn_cd_over_de[J];
    if (left.rows() == 0)
      continue;
    const auto &right = QJ_cd_ij[OrbitalPairIndex::block(int(J), parity_mn)];
    L4J_mn_ij[J] = Matrix(left.rows(), right.cols());
    LinAlg::GEMM(left, right, &L4J_mn_ij[J]);
  }
  return L4J_mn_ij;
}

//==============================================================================
std::vector<Matrix> L2_particle_hole(
  const std::vector<DiracSpinor> &excited, const std::vector<DiracSpinor> &core,
  const std::vector<DiracSpinor> &j_orbitals, const OrbitalPairIndex &mj_pairs,
  const std::vector<DiracSpinor> &n_orbitals,
  const std::vector<DiracSpinor> &i_orbitals, const OrbitalPairIndex &ni_pairs,
  const Coulomb::QkTable &qk, const Coulomb::LkTable *lk_prev,
  const SixJArray &sixj, const ExternalLine *external) {
  const auto max_K = std::min(mj_pairs.max_rank(), ni_pairs.max_rank());
  const OrbitalPairIndex rc_pairs(excited, core, max_K);

  const auto QLK_mj_rc_over_de = QL_particle_hole_over_de(
    excited, core, j_orbitals, mj_pairs, rc_pairs, qk, lk_prev, sixj, external);
  const auto QK_cr_ni = Q_particle_hole(excited, core, n_orbitals, i_orbitals,
                                        rc_pairs, ni_pairs, qk, sixj);

  // One matrix product for each block. The blocks are independent, so they
  // run in parallel, each with single-threaded BLAS: this does not rely on a
  // multi-threaded BLAS library
  std::vector<Matrix> L2K_mj_ni(mj_pairs.num_blocks());
  const qip::SingleThreadBlas single_thread_blas{};
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (int K = 0; K <= max_K; ++K) {
    for (int parity = 0; parity < 2; ++parity) {
      const auto block = OrbitalPairIndex::block(K, parity);
      // Empty unless rows, columns and the (r,c) pairs all have members
      if (QLK_mj_rc_over_de[block].rows() == 0 || QK_cr_ni[block].rows() == 0)
        continue;
      L2K_mj_ni[block] =
        Matrix(mj_pairs.size(K, parity), ni_pairs.size(K, parity));
      LinAlg::GEMM(QLK_mj_rc_over_de[block], QK_cr_ni[block],
                   &L2K_mj_ni[block]);
      // Net sign of the m-scheme minus sign and the recoupling phases
      L2K_mj_ni[block] *= double(Angular::neg1pow(K));
    }
  }
  return L2K_mj_ni;
}

//==============================================================================
void add_blocks(std::vector<Matrix> *sum, std::vector<Matrix> term) {
  assert(sum->size() == term.size());
  for (std::size_t index = 0; index < term.size(); ++index) {
    if (term[index].rows() == 0)
      continue;
    auto &sum_block = (*sum)[index];
    if (sum_block.rows() == 0) {
      sum_block = std::move(term[index]);
    } else {
      sum_block += term[index];
    }
  }
}

//==============================================================================
std::array<double, ladder_max_multipoles>
read_L1_L4_J(const std::vector<Matrix> &L1_L4_J_mn_ij, std::size_t row,
             const OrbitalPairIndex &ij_pairs, std::size_t index_i,
             std::size_t index_j, int J_min, int J_max) {
  std::array<double, ladder_max_multipoles> L1_L4_J{};
  for (int J = J_min; J <= J_max; ++J) {
    const auto &block = L1_L4_J_mn_ij[std::size_t(J)];
    if (block.rows() == 0)
      continue;
    L1_L4_J[std::size_t(J)] =
      block(row, std::size_t(ij_pairs.index(J, index_i, index_j)));
  }
  return L1_L4_J;
}

//==============================================================================
std::pair<std::array<double, ladder_max_multipoles>,
          std::array<double, ladder_max_multipoles>>
read_L2K(const std::vector<Matrix> &L2K_mj_ni, const OrbitalPairIndex &ph_pairs,
         int parity_mj, std::size_t index_m, std::size_t index_j,
         std::size_t index_n, std::size_t index_i, int K_min, int K_max) {
  // L2^K_{mj,ni} gives L2^k_mnij; L2^K_{ni,mj} gives L3^k_mnij = L2^k_nmji.
  // Both pairs are in the block of parity(m,j) = parity(n,i)
  return {read_L2K_element(L2K_mj_ni, ph_pairs, ph_pairs, parity_mj, index_m,
                           index_j, index_n, index_i, K_min, K_max),
          read_L2K_element(L2K_mj_ni, ph_pairs, ph_pairs, parity_mj, index_n,
                           index_i, index_m, index_j, K_min, K_max)};
}

//==============================================================================
std::array<double, ladder_max_multipoles>
read_L2K_element(const std::vector<Matrix> &L2K_mj_ni,
                 const OrbitalPairIndex &mj_pairs,
                 const OrbitalPairIndex &ni_pairs, int parity_mj,
                 std::size_t index_m, std::size_t index_j, std::size_t index_n,
                 std::size_t index_i, int K_min, int K_max) {
  std::array<double, ladder_max_multipoles> L2K{};
  for (int K = K_min; K <= K_max; ++K) {
    const auto &block = L2K_mj_ni[OrbitalPairIndex::block(K, parity_mj)];
    if (block.rows() == 0)
      continue;
    L2K[std::size_t(K)] =
      block(std::size_t(mj_pairs.index(K, index_m, index_j)),
            std::size_t(ni_pairs.index(K, index_n, index_i)));
  }
  return L2K;
}

//==============================================================================
double damped(double L_old, double L_new, double a_damp) {
  const auto a = std::clamp(a_damp, 0.0, 1.0);
  return (L_old == 0.0 || a == 0.0) ? L_new : (1.0 - a) * L_new + a * L_old;
}

//==============================================================================
void store_Lk_mnij(Coulomb::LkTable *lk, std::size_t index_m,
                   const std::vector<std::size_t> &n_group,
                   const std::vector<DiracSpinor> &excited,
                   const std::vector<DiracSpinor> &ij_orbitals,
                   const std::vector<char> &required_ij,
                   const OrbitalPairIndex &ij_pairs,
                   const OrbitalPairIndex &ph_pairs,
                   const std::vector<Matrix> &L1_L4_J_mn_ij,
                   const std::vector<Matrix> &L2K_mj_ni,
                   const LadderDiagrams &diagrams, const SixJArray &sixj,
                   int max_k, bool only_new, double a_damp) {
  const auto &m = excited[index_m];
  const auto num_ij = ij_orbitals.size();
  const auto max_J = ij_pairs.max_rank();
  const auto max_K = ph_pairs.max_rank();

  for (std::size_t row = 0; row < n_group.size(); ++row) {
    const auto index_n = n_group[row];
    const auto &n = excited[index_n];
    const auto parity_mn = pair_parity(m, n);
    for (std::size_t index_i = 0; index_i < num_ij; ++index_i) {
      const auto &i = ij_orbitals[index_i];
      for (std::size_t index_j = 0; index_j < num_ij; ++index_j) {
        const auto &j = ij_orbitals[index_j];
        if (!required_ij[index_i * num_ij + index_j] ||
            pair_parity(i, j) != parity_mn)
          continue;
        // For m = n the two orders (i,j) and (j,i) are the same integral
        if (index_m == index_n &&
            !lk->is_NormalOrdered(m.nk_index(), n.nk_index(), i.nk_index(),
                                  j.nk_index()))
          continue;
        const auto [k_min, k_max_L] = k_minmax_L(m, n, i, j);
        const auto k_max = std::min(k_max_L, max_k);
        if (k_max < k_min)
          continue;

        const auto [J_min, J_max] = coupled_rank_range(m, n, i, j, max_J);
        const auto L1_L4_J = read_L1_L4_J(L1_L4_J_mn_ij, row, ij_pairs, index_i,
                                          index_j, J_min, J_max);
        const auto [K_min, K_max] = coupled_rank_range(m, j, n, i, max_K);
        const auto [L2K_for_L2, L2K_for_L3] =
          read_L2K(L2K_mj_ni, ph_pairs, pair_parity(m, j), index_m, index_j,
                   index_n, index_i, K_min, K_max);

        // L_mnij and L_nmji are one entry, and L2^k_nmji = L3^k_mnij. With
        // only one of L2, L3, the entry holds that diagram with the later of
        // i, j in the third index (L^k_mnib, b core): for j after i, the two
        // trade places
        const bool as_nmji = index_j > index_i;
        const bool with_L2 = as_nmji ? diagrams.L3 : diagrams.L2;
        const bool with_L3 = as_nmji ? diagrams.L2 : diagrams.L3;

        const auto key = lk->NormalOrder(m.nk_index(), n.nk_index(),
                                         i.nk_index(), j.nk_index());
        for (int k = k_min; k <= k_max; ++k) {
          double *const L_stored = lk->get(k, key);
          if (L_stored == nullptr)
            continue;
          if (only_new && *L_stored != 0.0)
            continue;
          const auto L1_L4 = multipole_from_pair_coupled(k, m, n, i, j, L1_L4_J,
                                                         J_min, J_max, sixj);
          const auto L2 = multipole_from_particle_hole_coupled(
            k, m, n, i, j, L2K_for_L2, K_min, K_max, sixj);
          // L3^k_mnij = L2^k_nmji
          const auto L3 = multipole_from_particle_hole_coupled(
            k, n, m, j, i, L2K_for_L3, K_min, K_max, sixj);
          *L_stored =
            damped(*L_stored,
                   L1_L4 + (with_L2 ? L2 : 0.0) + (with_L3 ? L3 : 0.0), a_damp);
        }
      }
    }
  }
}

//==============================================================================
void add_Lk_mnib_entries(Coulomb::LkTable *lk,
                         const std::vector<DiracSpinor> &excited,
                         const std::vector<DiracSpinor> &core,
                         const std::vector<DiracSpinor> &i_orbs, int max_k) {
  const auto ij_orbitals = ladder_ij_orbitals(core, excited, i_orbs);
  const auto required_ij = required_ij_pairs(ij_orbitals, core, i_orbs);
  const auto num_ij = ij_orbitals.size();

  const auto num_k = std::size_t(max_k) + 1;
  if ((*lk)->size() < num_k) {
    (*lk)->resize(num_k);
  }

  // Each k has its own map, so the entries can be inserted in parallel over k.
  // (m,n,i,j) and (n,m,j,i) are the same entry: only the normal-ordered one
#pragma omp parallel for
  for (std::size_t k = 0; k < num_k; ++k) {
    auto &Lk_map = (*lk)->at(k);
    // First pass counts (to reserve space); second pass inserts
    for (const bool insert : {false, true}) {
      std::size_t count = 0;
      for (std::size_t index_i = 0; index_i < num_ij; ++index_i) {
        const auto &i = ij_orbitals[index_i];
        for (std::size_t index_j = 0; index_j < num_ij; ++index_j) {
          const auto &j = ij_orbitals[index_j];
          if (!required_ij[index_i * num_ij + index_j])
            continue;
          for (const auto &m : excited) {
            for (const auto &n : excited) {
              const auto [k0, kI] = k_minmax_L(m, n, i, j);
              if (int(k) < k0 || int(k) > kI)
                continue;
              if (!lk->is_NormalOrdered(m.nk_index(), n.nk_index(),
                                        i.nk_index(), j.nk_index()))
                continue;
              if (insert) {
                Lk_map.insert({lk->NormalOrder(m.nk_index(), n.nk_index(),
                                               i.nk_index(), j.nk_index()),
                               0.0});
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
void calculate_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &i_orbs,
                       const LadderDiagrams &diagrams,
                       const Coulomb::LkTable *lk_prev, int max_k,
                       bool only_new, double a_damp,
                       const ExternalLine *external) {
  assert(lk != lk_prev);
  // With the line as the m orbital, L2 and L3 need different matrices
  assert(external == nullptr || external->L_xnab == nullptr);
  if (excited.empty() || core.empty() || i_orbs.empty())
    return;

  // Orbitals of the (i,j) indices, and which (i,j) pairs are required
  const auto ij_orbitals = ladder_ij_orbitals(core, excited, i_orbs);
  const auto required_ij = required_ij_pairs(ij_orbitals, core, i_orbs);

  // Largest angular momenta: J of an (i,j) pair, K of a particle-hole pair
  const auto max_tj_excited = DiracSpinor::max_tj(excited);
  const auto max_tj_hole =
    std::max(DiracSpinor::max_tj(core), DiracSpinor::max_tj(ij_orbitals));
  const auto max_tj = std::max(max_tj_excited, max_tj_hole);
  const int max_J = std::min(DiracSpinor::max_tj(ij_orbitals), ladder_max_k);
  const int max_K = std::min((max_tj_excited + max_tj_hole) / 2, ladder_max_k);
  // every k, J and K is at most (2j_a + 2j_b)/2 <= max_tj
  const SixJArray sixj(max_tj, std::min(max_tj, ladder_max_k));

  // Rows and columns of the coupled matrices
  const OrbitalPairIndex rs_pairs(excited, excited, max_J);
  const OrbitalPairIndex cd_pairs(core, core, max_J);
  const OrbitalPairIndex ij_pairs(ij_orbitals, ij_orbitals, max_J, required_ij);
  const OrbitalPairIndex ph_pairs(excited, ij_orbitals, max_K);

  // Factors shared by every (m,n): empty blocks for an excluded diagram
  const auto QLJ_rs_ij_over_de =
    diagrams.L1 ?
      QL_particle_particle_over_de(excited, ij_orbitals, required_ij, rs_pairs,
                                   ij_pairs, qk, lk_prev, sixj, external) :
      std::vector<Matrix>(rs_pairs.num_blocks());
  const auto QJ_cd_ij = diagrams.L4 ?
                          Q_hole_hole(core, ij_orbitals, required_ij, cd_pairs,
                                      ij_pairs, qk, sixj) :
                          std::vector<Matrix>(cd_pairs.num_blocks());
  // L2 and L3 both come from L2^K
  const auto L2K_mj_ni =
    diagrams.L2 || diagrams.L3 ?
      L2_particle_hole(excited, core, ij_orbitals, ph_pairs, excited,
                       ij_orbitals, ph_pairs, qk, lk_prev, sixj, external) :
      std::vector<Matrix>(ph_pairs.num_blocks());

  // The n of one kappa group share all angular factors and block structure:
  // one matrix (rows n) per J. Only n >= m: L^k_mnij = L^k_nmji
  const auto kappa_groups = group_by_kappa(excited);
  const auto num_groups = kappa_groups.size();

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t index_m = 0; index_m < excited.size(); ++index_m) {
    for (std::size_t index_group = 0; index_group < num_groups; ++index_group) {
      const auto &m = excited[index_m];
      std::vector<std::size_t> n_group;
      for (const auto index_n : kappa_groups[index_group]) {
        if (index_n >= index_m) {
          n_group.push_back(index_n);
        }
      }
      if (n_group.empty())
        continue;

      // Pair-coupled (L1 + L4)^J_{mn,ij}: one matrix per J, rows n
      auto L1_L4_J_mn_ij =
        diagrams.L1 ?
          L1_particle_particle(m, n_group, excited, rs_pairs, ij_pairs,
                               QLJ_rs_ij_over_de, qk, sixj) :
          std::vector<Matrix>(std::size_t(max_J) + 1);
      if (diagrams.L4) {
        add_blocks(&L1_L4_J_mn_ij,
                   L4_hole_hole(m, n_group, excited, core, cd_pairs, ij_pairs,
                                QJ_cd_ij, qk, lk_prev, sixj));
      }

      store_Lk_mnij(lk, index_m, n_group, excited, ij_orbitals, required_ij,
                    ij_pairs, ph_pairs, L1_L4_J_mn_ij, L2K_mj_ni, diagrams,
                    sixj, max_k, only_new, a_damp);
    }
  }
}

} // namespace MBPT
