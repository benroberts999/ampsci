#include "MBPT/LadderExternalLine.hpp"
#include "Coulomb/include.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/LadderIteration.hpp"
#include "MBPT/SpinorMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <vector>

namespace MBPT {

using Matrix = LinAlg::Matrix<double>;

//==============================================================================
ExternalLineLadder external_line_start(const DiracSpinor &v,
                                       const std::vector<DiracSpinor> &core,
                                       const std::vector<DiracSpinor> &excited,
                                       const Coulomb::LkTable &lk) {
  assert(std::find(excited.cbegin(), excited.cend(), v) != excited.cend());
  ExternalLineLadder ladder{v.en(),
                            qip::select_if(excited, [&v](const DiracSpinor &x) {
                              return x.kappa() == v.kappa();
                            })};
  const auto max_k = lk.max_k();

  add_Lk_mnib_entries(&ladder.L_mnxb, excited, core, ladder.x_orbitals, max_k);
  for (const auto &m : excited) {
    for (const auto &n : excited) {
      for (const auto &b : core) {
        const auto [k_min, k_max] = k_minmax_L(m, n, v, b);
        for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
          double *const L_v =
            ladder.L_mnxb.get(k, ladder.L_mnxb.NormalOrder(m, n, v, b));
          if (L_v != nullptr) {
            *L_v = lk.Q(k, m, n, v, b);
          }
        }
      }
    }
  }

  for (const auto &x : ladder.x_orbitals) {
    for (const auto &n : excited) {
      for (const auto &a : core) {
        for (const auto &b : core) {
          const auto [k_min, k_max] = k_minmax_L(x, n, a, b);
          for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
            ladder.L_xnab.add(k, x, n, a, b,
                              x == v ? lk.Q(k, v, n, a, b) : 0.0);
          }
        }
      }
    }
  }
  return ladder;
}

//==============================================================================
void calculate_Lk_mnxb(Coulomb::LkTable *L_mnxb, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &x_orbitals, double en_v,
                       const LadderDiagrams &diagrams,
                       const Coulomb::LkTable &lk,
                       const Coulomb::LkTable &L_mnxb_prev, double a_damp) {
  const ExternalLine line{x_orbitals, en_v, &L_mnxb_prev, nullptr};
  calculate_Lk_mnib(L_mnxb, qk, excited, core, x_orbitals, diagrams, &lk,
                    L_mnxb->max_k(), false, a_damp, &line);
}

//==============================================================================
void calculate_Lk_xnab(Coulomb::NkTable *L_xnab, const Coulomb::QkTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &x_orbitals, double en_v,
                       const LadderDiagrams &diagrams,
                       const Coulomb::LkTable &lk,
                       const Coulomb::NkTable &L_xnab_prev, double a_damp) {
  // P with the external line is not derived
  assert(!(diagrams.P1 || diagrams.P2 || diagrams.P3));
  if (excited.empty() || core.empty() || x_orbitals.empty())
    return;
  const ExternalLine line{x_orbitals, en_v, nullptr, &L_xnab_prev};

  // As calculate_Lk_mnib() with i_orbs = core: (i,j) = (a,b), all core pairs
  const auto required_ab = required_ij_pairs(core, core, core);
  const auto max_tj_excited = DiracSpinor::max_tj(excited);
  const auto max_tj_core = DiracSpinor::max_tj(core);
  const auto max_tj = std::max(max_tj_excited, max_tj_core);
  const int max_J = std::min(max_tj_core, ladder_max_k);
  const int max_K = std::min((max_tj_excited + max_tj_core) / 2, ladder_max_k);
  const SixJArray sixj(max_tj, std::min(max_tj, ladder_max_k));
  const OrbitalPairIndex rs_pairs(excited, excited, max_J);
  const OrbitalPairIndex ab_pairs(core, core, max_J);
  const OrbitalPairIndex ph_pairs(excited, core, max_K);

  const auto QLJ_rs_ab_over_de =
    diagrams.L1 ?
      QL_particle_particle_over_de(excited, core, required_ab, rs_pairs,
                                   ab_pairs, qk, &lk, sixj) :
      std::vector<Matrix>(rs_pairs.num_blocks());
  const auto QJ_cd_ab =
    diagrams.L4 ?
      Q_hole_hole(core, core, required_ab, ab_pairs, ab_pairs, qk, sixj) :
      std::vector<Matrix>(ab_pairs.num_blocks());
  // L2^K_{xb,na}: rows (x,b) hold the line. L2^K_{na,xb} (for L3): no line
  const auto L2K_mb_na_line =
    diagrams.L2 ? L2_particle_hole(excited, core, core, ph_pairs, excited, core,
                                   ph_pairs, qk, &lk, sixj, &line) :
                  std::vector<Matrix>(ph_pairs.num_blocks());
  const auto L2K_mb_na =
    diagrams.L3 ? L2_particle_hole(excited, core, core, ph_pairs, excited, core,
                                   ph_pairs, qk, &lk, sixj) :
                  std::vector<Matrix>(ph_pairs.num_blocks());

  // Positions of the x in excited; all n of each kappa group (no symmetry)
  std::vector<std::size_t> x_positions;
  for (std::size_t index = 0; index < excited.size(); ++index) {
    if (line.is_m_line(excited[index])) {
      x_positions.push_back(index);
    }
  }
  const auto kappa_groups = group_by_kappa(excited);

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t ix = 0; ix < x_positions.size(); ++ix) {
    for (std::size_t index_group = 0; index_group < kappa_groups.size();
         ++index_group) {
      const auto &x = excited[x_positions[ix]];
      const auto &n_group = kappa_groups[index_group];
      auto L1_L4_J_xn_ab =
        diagrams.L1 ?
          L1_particle_particle(x, n_group, excited, rs_pairs, ab_pairs,
                               QLJ_rs_ab_over_de, qk, sixj) :
          std::vector<Matrix>(std::size_t(max_J) + 1);
      if (diagrams.L4) {
        add_blocks(&L1_L4_J_xn_ab,
                   L4_hole_hole(x, n_group, excited, core, ab_pairs, ab_pairs,
                                QJ_cd_ab, qk, &lk, sixj, &line));
      }
      store_Lk_xnab(L_xnab, x_positions[ix], n_group, excited, core, ab_pairs,
                    ph_pairs, L1_L4_J_xn_ab, L2K_mb_na_line, L2K_mb_na, sixj,
                    a_damp);
    }
  }
}

//==============================================================================
void store_Lk_xnab(Coulomb::NkTable *L_xnab, std::size_t index_x,
                   const std::vector<std::size_t> &n_group,
                   const std::vector<DiracSpinor> &excited,
                   const std::vector<DiracSpinor> &core,
                   const OrbitalPairIndex &ab_pairs,
                   const OrbitalPairIndex &ph_pairs,
                   const std::vector<Matrix> &L1_L4_J_xn_ab,
                   const std::vector<Matrix> &L2K_mb_na_line,
                   const std::vector<Matrix> &L2K_mb_na, const SixJArray &sixj,
                   double a_damp) {
  const auto &x = excited[index_x];
  for (std::size_t row = 0; row < n_group.size(); ++row) {
    const auto index_n = n_group[row];
    const auto &n = excited[index_n];
    for (std::size_t index_a = 0; index_a < core.size(); ++index_a) {
      const auto &a = core[index_a];
      for (std::size_t index_b = 0; index_b < core.size(); ++index_b) {
        const auto &b = core[index_b];
        const auto [k_min, k_max] = k_minmax_L(x, n, a, b);
        if (k_max < k_min)
          continue;

        const auto [J_min, J_max] =
          coupled_rank_range(x, n, a, b, ab_pairs.max_rank());
        const auto L1_L4_J = read_L1_L4_J(L1_L4_J_xn_ab, row, ab_pairs, index_a,
                                          index_b, J_min, J_max);
        const auto [K_min, K_max] =
          coupled_rank_range(x, b, n, a, ph_pairs.max_rank());
        // Both pairs are in the block of parity(x,b) = parity(n,a)
        const auto parity_xb = pair_parity(x, b);
        const auto L2K_for_L2 =
          read_L2K_element(L2K_mb_na_line, ph_pairs, ph_pairs, parity_xb,
                           index_x, index_b, index_n, index_a, K_min, K_max);
        const auto L2K_for_L3 =
          read_L2K_element(L2K_mb_na, ph_pairs, ph_pairs, parity_xb, index_n,
                           index_a, index_x, index_b, K_min, K_max);

        const auto key = L_xnab->NormalOrder(x, n, a, b);
        for (int k = k_min; k <= k_max; ++k) {
          double *const L_stored = L_xnab->get(k, key);
          if (L_stored == nullptr)
            continue;
          const auto L1_L4 = multipole_from_pair_coupled(k, x, n, a, b, L1_L4_J,
                                                         J_min, J_max, sixj);
          const auto L2 = multipole_from_particle_hole_coupled(
            k, x, n, a, b, L2K_for_L2, K_min, K_max, sixj);
          // L3^k_xnab = L2^k_nxba
          const auto L3 = multipole_from_particle_hole_coupled(
            k, n, x, b, a, L2K_for_L3, K_min, K_max, sixj);
          *L_stored = damped(*L_stored, L1_L4 + L2 + L3, a_damp);
        }
      }
    }
  }
}

//==============================================================================
double external_line_change(ExternalLineLadder &ladder,
                            const ExternalLineLadder &ladder_prev) {
  double max_change = 0.0;
  double largest = 0.0;
  for (std::size_t k = 0; k < ladder.L_mnxb->size(); ++k) {
    for (const auto &[key, value] : ladder.L_mnxb->at(k)) {
      max_change = std::max(
        max_change, std::abs(value - ladder_prev.L_mnxb.Q(int(k), key)));
      largest = std::max(largest, std::abs(value));
    }
  }
  for (std::size_t k = 0; k < ladder.L_xnab->size(); ++k) {
    for (const auto &[key, value] : ladder.L_xnab->at(k)) {
      max_change = std::max(
        max_change, std::abs(value - ladder_prev.L_xnab.Q(int(k), key)));
      largest = std::max(largest, std::abs(value));
    }
  }
  return largest == 0.0 ? 0.0 : max_change / largest;
}

//==============================================================================
ExternalLineLadder
external_line_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable &lk,
                     const LadderDiagrams &diagrams, int max_it,
                     double eps_target, double a_damp, bool print) {
  auto ladder = external_line_start(v, core, excited, lk);
  // Not converged until at least one sweep is done
  double eps = 1.0;
  for (int it = 1; it <= max_it; ++it) {
    const auto ladder_prev = ladder;
    calculate_Lk_mnxb(&ladder.L_mnxb, qk, excited, core, ladder.x_orbitals,
                      ladder.en_v, diagrams, lk, ladder_prev.L_mnxb, a_damp);
    calculate_Lk_xnab(&ladder.L_xnab, qk, excited, core, ladder.x_orbitals,
                      ladder.en_v, diagrams, lk, ladder_prev.L_xnab, a_damp);
    eps = external_line_change(ladder, ladder_prev);
    if (print) {
      fmt::print("{:4s} external line it:{:2}  eps = {:.1e}\n", v.shortSymbol(),
                 it, eps);
      std::cout << std::flush;
    }
    if (eps < eps_target)
      break;
  }
  if (eps >= eps_target) {
    fmt::print("Warning: external line of {} not converged in {} sweeps\n",
               v.shortSymbol(), max_it);
  }
  return ladder;
}

//==============================================================================
template <typename MnxbTable, typename XnabTable>
GMatrix
Sigma_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &x_orbitals,
             const MnxbTable &L_mnxb, const XnabTable &L_xnab,
             const std::vector<DiracSpinor> &core,
             const std::vector<DiracSpinor> &excited, double r0, double rmax,
             std::size_t stride, bool include_G, bool include_exchange) {
  const auto grid = v.grid_sptr();
  const auto i0 = grid->getIndex(r0);
  const auto size = (grid->getIndex(rmax) - i0) / stride + 1;
  GMatrix Sigma{i0, stride, size, include_G, grid};
  std::vector<GMatrix> Sigma_thread(std::size_t(omp_get_max_threads()), Sigma);
  const auto max_k_mnxb = L_mnxb.max_k();
  const auto max_k_xnab = L_xnab.max_k();

#pragma omp parallel for schedule(dynamic)
  for (std::size_t index_n = 0; index_n < excited.size(); ++index_n) {
    const auto &n = excited[index_n];
    auto &Sigma_t = Sigma_thread[std::size_t(omp_get_thread_num())];
    // sum_x L_x |x>, for one set of the other indices
    DiracSpinor L_projected(0, v.kappa(), grid);
    for (const auto &a : core) {

      // |W^k_{.amn}> L^k_{mnxa} <x| / ([k][j_v] de); or Q in place of W
      for (const auto &m : excited) {
        const auto [k_min, k_max] = k_minmax_L(m, n, v, a);
        const auto de = v.en() + a.en() - m.en() - n.en();
        for (int k = k_min; k <= std::min(k_max, max_k_mnxb); ++k) {
          if (!include_exchange && !Coulomb::Qk_abcd_SR(k, v.kappa(), a.kappa(),
                                                        m.kappa(), n.kappa()))
            continue;
          L_projected *= 0.0;
          for (const auto &x : x_orbitals) {
            L_projected += L_mnxb.Q(k, m, n, x, a) * x;
          }
          Sigma_t.add(include_exchange ?
                        Coulomb::Wkv_bcd(k, v.kappa(), a, m, n) :
                        Coulomb::Qkv_bcd(k, v.kappa(), a, m, n),
                      L_projected, 1.0 / ((2 * k + 1) * (v.twoj() + 1) * de));
        }
      }

      // |W^k_{.nab}> L^k_{xnab} <x| / ([k][j_v] de); or Q in place of W
      for (const auto &b : core) {
        const auto [k_min, k_max] = k_minmax_L(v, n, a, b);
        const auto de = v.en() + n.en() - a.en() - b.en();
        for (int k = k_min; k <= std::min(k_max, max_k_xnab); ++k) {
          if (!include_exchange && !Coulomb::Qk_abcd_SR(k, v.kappa(), n.kappa(),
                                                        a.kappa(), b.kappa()))
            continue;
          L_projected *= 0.0;
          for (const auto &x : x_orbitals) {
            L_projected += L_xnab.Q(k, x, n, a, b) * x;
          }
          Sigma_t.add(include_exchange ?
                        Coulomb::Wkv_bcd(k, v.kappa(), n, a, b) :
                        Coulomb::Qkv_bcd(k, v.kappa(), n, a, b),
                      L_projected, 1.0 / ((2 * k + 1) * (v.twoj() + 1) * de));
        }
      }
    }
  }
  for (const auto &Sigma_t : Sigma_thread) {
    Sigma += Sigma_t;
  }
  return Sigma.drj_in_place();
}

template GMatrix
Sigma_ladder(const DiracSpinor &, const std::vector<DiracSpinor> &,
             const Coulomb::LkTable &, const Coulomb::NkTable &,
             const std::vector<DiracSpinor> &, const std::vector<DiracSpinor> &,
             double, double, std::size_t, bool, bool);
template GMatrix
Sigma_ladder(const DiracSpinor &, const std::vector<DiracSpinor> &,
             const Coulomb::NkTable &, const Coulomb::NkTable &,
             const std::vector<DiracSpinor> &, const std::vector<DiracSpinor> &,
             double, double, std::size_t, bool, bool);
template GMatrix
Sigma_ladder(const DiracSpinor &, const std::vector<DiracSpinor> &,
             const Coulomb::LkTable &, const Coulomb::LkTable &,
             const std::vector<DiracSpinor> &, const std::vector<DiracSpinor> &,
             double, double, std::size_t, bool, bool);

//==============================================================================
GMatrix Sigma_ladder(const DiracSpinor &v, const ExternalLineLadder &ladder,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, double r0,
                     double rmax, std::size_t stride, bool include_G) {
  return Sigma_ladder(v, ladder.x_orbitals, ladder.L_mnxb, ladder.L_xnab, core,
                      excited, r0, rmax, stride, include_G, true);
}

//==============================================================================
Coulomb::NkTable
lowest_order_Lk_mnxb(const Coulomb::QkTable &qk,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &x_orbitals,
                     const LadderDiagrams &diagrams, int max_k) {
  // L1 and L3 have e_x in an energy denominator
  assert(!diagrams.L1 && !diagrams.L3);
  // P with the external line is not derived
  assert(!(diagrams.P1 || diagrams.P2 || diagrams.P3));
  Coulomb::NkTable L_mnxb;
  for (const auto &m : excited) {
    for (const auto &n : excited) {
      for (const auto &x : x_orbitals) {
        for (const auto &b : core) {
          const auto [k_min, k_max] = k_minmax_L(m, n, x, b);
          for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
            L_mnxb.add(k, m, n, x, b, 0.0);
          }
        }
      }
    }
  }
  if (L_mnxb.count() == 0)
    return L_mnxb;

  // The (i,j) orbitals: the x, then the core, by position (an x may have the
  // key of a core orbital). Only the pairs (x,b) are required
  const auto num_x = x_orbitals.size();
  const auto xb_orbitals = qip::merge(x_orbitals, core);
  const auto num_xb = xb_orbitals.size();
  std::vector<char> required_xb(num_xb * num_xb, 0);
  for (std::size_t index_x = 0; index_x < num_x; ++index_x) {
    for (std::size_t index_b = 0; index_b < core.size(); ++index_b) {
      required_xb[index_x * num_xb + num_x + index_b] = 1;
    }
  }

  const auto max_tj_excited = DiracSpinor::max_tj(excited);
  const auto max_tj_hole = DiracSpinor::max_tj(xb_orbitals);
  const auto max_tj = std::max(max_tj_excited, max_tj_hole);
  const int max_J = std::min(max_tj_hole, ladder_max_k);
  const int max_K = std::min((max_tj_excited + max_tj_hole) / 2, ladder_max_k);
  const SixJArray sixj(max_tj, std::min(max_tj, ladder_max_k));
  const OrbitalPairIndex cd_pairs(core, core, max_J);
  const OrbitalPairIndex xb_pairs(xb_orbitals, xb_orbitals, max_J, required_xb);
  const OrbitalPairIndex mb_pairs(excited, core, max_K);
  const OrbitalPairIndex nx_pairs(excited, x_orbitals, max_K);

  const auto QJ_cd_xb = diagrams.L4 ?
                          Q_hole_hole(core, xb_orbitals, required_xb, cd_pairs,
                                      xb_pairs, qk, sixj) :
                          std::vector<Matrix>(cd_pairs.num_blocks());
  // L2^K_{mb,nx}: x is the i orbital of the bare Q^K
  const auto L2K_mb_nx =
    diagrams.L2 ? L2_particle_hole(excited, core, core, mb_pairs, excited,
                                   x_orbitals, nx_pairs, qk, nullptr, sixj) :
                  std::vector<Matrix>(mb_pairs.num_blocks());

  // All n of each kappa group: L_mnxb and L_nmbx are different integrals
  const auto kappa_groups = group_by_kappa(excited);

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t index_m = 0; index_m < excited.size(); ++index_m) {
    for (std::size_t index_group = 0; index_group < kappa_groups.size();
         ++index_group) {
      const auto &m = excited[index_m];
      const auto &n_group = kappa_groups[index_group];
      const auto L4J_mn_xb =
        diagrams.L4 ? L4_hole_hole(m, n_group, excited, core, cd_pairs,
                                   xb_pairs, QJ_cd_xb, qk, nullptr, sixj) :
                      std::vector<Matrix>(std::size_t(max_J) + 1);
      for (std::size_t row = 0; row < n_group.size(); ++row) {
        const auto index_n = n_group[row];
        const auto &n = excited[index_n];
        for (std::size_t index_x = 0; index_x < num_x; ++index_x) {
          const auto &x = x_orbitals[index_x];
          for (std::size_t index_b = 0; index_b < core.size(); ++index_b) {
            const auto &b = core[index_b];
            const auto [k_min, k_max_L] = k_minmax_L(m, n, x, b);
            const auto k_max = std::min(k_max_L, max_k);
            if (k_max < k_min)
              continue;
            const auto [J_min, J_max] = coupled_rank_range(m, n, x, b, max_J);
            const auto L4J = read_L1_L4_J(L4J_mn_xb, row, xb_pairs, index_x,
                                          num_x + index_b, J_min, J_max);
            const auto [K_min, K_max] = coupled_rank_range(m, b, n, x, max_K);
            const auto L2K = read_L2K_element(
              L2K_mb_nx, mb_pairs, nx_pairs, pair_parity(m, b), index_m,
              index_b, index_n, index_x, K_min, K_max);
            const auto key = L_mnxb.NormalOrder(m, n, x, b);
            for (int k = k_min; k <= k_max; ++k) {
              *L_mnxb.get(k, key) = multipole_from_pair_coupled(
                                      k, m, n, x, b, L4J, J_min, J_max, sixj) +
                                    multipole_from_particle_hole_coupled(
                                      k, m, n, x, b, L2K, K_min, K_max, sixj);
            }
          }
        }
      }
    }
  }
  return L_mnxb;
}

//==============================================================================
Coulomb::NkTable
lowest_order_Lk_xnab(const Coulomb::QkTable &qk,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &x_orbitals,
                     const LadderDiagrams &diagrams, int max_k) {
  // L2 and L4 have e_x in an energy denominator
  assert(!diagrams.L2 && !diagrams.L4);
  // P with the external line is not derived
  assert(!(diagrams.P1 || diagrams.P2 || diagrams.P3));
  Coulomb::NkTable L_xnab;
  for (const auto &x : x_orbitals) {
    for (const auto &n : excited) {
      for (const auto &a : core) {
        for (const auto &b : core) {
          const auto [k_min, k_max] = k_minmax_L(x, n, a, b);
          for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
            L_xnab.add(k, x, n, a, b, 0.0);
          }
        }
      }
    }
  }
  if (L_xnab.count() == 0)
    return L_xnab;

  // (i,j) = (a,b), all core pairs
  const auto required_ab = required_ij_pairs(core, core, core);
  const auto max_tj_excited = DiracSpinor::max_tj(excited);
  const auto max_tj_core = DiracSpinor::max_tj(core);
  const auto max_tj =
    std::max({max_tj_excited, max_tj_core, DiracSpinor::max_tj(x_orbitals)});
  const int max_J = std::min(max_tj_core, ladder_max_k);
  const int max_K = std::min((max_tj_excited + max_tj_core) / 2, ladder_max_k);
  const SixJArray sixj(max_tj, std::min(max_tj, ladder_max_k));
  const OrbitalPairIndex rs_pairs(excited, excited, max_J);
  const OrbitalPairIndex ab_pairs(core, core, max_J);
  const OrbitalPairIndex na_pairs(excited, core, max_K);
  const OrbitalPairIndex xb_pairs(x_orbitals, core, max_K);

  const auto QLJ_rs_ab_over_de =
    diagrams.L1 ?
      QL_particle_particle_over_de(excited, core, required_ab, rs_pairs,
                                   ab_pairs, qk, nullptr, sixj) :
      std::vector<Matrix>(rs_pairs.num_blocks());
  // L3^k_xnab = L2^k_nxba, from L2^K_{na,xb}: x is the n orbital of the bare
  // Q^K
  const auto L2K_na_xb =
    diagrams.L3 ? L2_particle_hole(excited, core, core, na_pairs, x_orbitals,
                                   core, xb_pairs, qk, nullptr, sixj) :
                  std::vector<Matrix>(na_pairs.num_blocks());

  const auto kappa_groups = group_by_kappa(excited);

  // BLAS must run single-threaded inside the omp region below
  const qip::SingleThreadBlas single_thread_blas{};

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t index_x = 0; index_x < x_orbitals.size(); ++index_x) {
    for (std::size_t index_group = 0; index_group < kappa_groups.size();
         ++index_group) {
      const auto &x = x_orbitals[index_x];
      const auto &n_group = kappa_groups[index_group];
      const auto L1J_xn_ab =
        diagrams.L1 ?
          L1_particle_particle(x, n_group, excited, rs_pairs, ab_pairs,
                               QLJ_rs_ab_over_de, qk, sixj) :
          std::vector<Matrix>(std::size_t(max_J) + 1);
      for (std::size_t row = 0; row < n_group.size(); ++row) {
        const auto index_n = n_group[row];
        const auto &n = excited[index_n];
        for (std::size_t index_a = 0; index_a < core.size(); ++index_a) {
          const auto &a = core[index_a];
          for (std::size_t index_b = 0; index_b < core.size(); ++index_b) {
            const auto &b = core[index_b];
            const auto [k_min, k_max_L] = k_minmax_L(x, n, a, b);
            const auto k_max = std::min(k_max_L, max_k);
            if (k_max < k_min)
              continue;
            const auto [J_min, J_max] = coupled_rank_range(x, n, a, b, max_J);
            const auto L1J = read_L1_L4_J(L1J_xn_ab, row, ab_pairs, index_a,
                                          index_b, J_min, J_max);
            const auto [K_min, K_max] = coupled_rank_range(n, a, x, b, max_K);
            const auto L2K = read_L2K_element(
              L2K_na_xb, na_pairs, xb_pairs, pair_parity(n, a), index_n,
              index_a, index_x, index_b, K_min, K_max);
            const auto key = L_xnab.NormalOrder(x, n, a, b);
            for (int k = k_min; k <= k_max; ++k) {
              // L3^k_xnab = L2^k_nxba
              *L_xnab.get(k, key) = multipole_from_pair_coupled(
                                      k, x, n, a, b, L1J, J_min, J_max, sixj) +
                                    multipole_from_particle_hole_coupled(
                                      k, n, x, b, a, L2K, K_min, K_max, sixj);
            }
          }
        }
      }
    }
  }
  return L_xnab;
}

//==============================================================================
GMatrix Sigma_ladder_mirror(const DiracSpinor &v,
                            const std::vector<DiracSpinor> &x_orbitals,
                            const std::vector<DiracSpinor> &core,
                            const std::vector<DiracSpinor> &excited,
                            const Coulomb::QkTable &qk, int max_k,
                            bool include_L4, double r0, double rmax,
                            std::size_t stride, bool include_G) {
  const LadderDiagrams diagram_L1{true, false, false, false};
  const LadderDiagrams diagram_L2{false, true, false, false};
  const LadderDiagrams diagram_L3{false, false, true, false};
  const LadderDiagrams diagram_L4{false, false, false, true};

  // Closed with Q: L2 with x as i, L3 with x as m
  auto Sigma = Sigma_ladder(
    v, x_orbitals,
    lowest_order_Lk_mnxb(qk, excited, core, x_orbitals, diagram_L2, max_k),
    lowest_order_Lk_xnab(qk, excited, core, x_orbitals, diagram_L3, max_k),
    core, excited, r0, rmax, stride, include_G, false);

  // Closed with W: L4 with x as i, L1 with x as m
  Sigma += Sigma_ladder(
    v, x_orbitals,
    include_L4 ?
      lowest_order_Lk_mnxb(qk, excited, core, x_orbitals, diagram_L4, max_k) :
      Coulomb::NkTable{},
    lowest_order_Lk_xnab(qk, excited, core, x_orbitals, diagram_L1, max_k),
    core, excited, r0, rmax, stride, include_G, true);

  // The mirror image of a term of the correlation potential is its transpose
  return Sigma.transpose_drj();
}

//==============================================================================
std::vector<DiracSpinor>
projection_states(const DiracSpinor &v, const std::vector<DiracSpinor> &basis) {
  return qip::select_if(
    basis, [&v](const DiracSpinor &x) { return x.kappa() == v.kappa(); });
}

//==============================================================================
std::size_t fill_Qk_projection_states(Coulomb::QkTable *qk,
                                      const std::vector<DiracSpinor> &orbitals,
                                      const std::vector<DiracSpinor> &basis,
                                      const std::vector<DiracSpinor> &valence,
                                      int max_k) {
  const auto count_before = qk->count();

  // The projection states of every valence kappa (each kappa once) that are
  // not in orbitals
  std::vector<DiracSpinor> new_states;
  std::vector<int> kappas_done;
  for (const auto &v : valence) {
    if (std::find(kappas_done.cbegin(), kappas_done.cend(), v.kappa()) !=
        kappas_done.cend())
      continue;
    kappas_done.push_back(v.kappa());
    for (const auto &x : projection_states(v, basis)) {
      if (std::find(orbitals.cbegin(), orbitals.cend(), x) == orbitals.cend()) {
        new_states.push_back(x);
      }
    }
  }
  if (new_states.empty())
    return 0;
  const auto orbitals_and_x = qip::merge(orbitals, new_states);
  const Coulomb::YkTable yk(orbitals_and_x);

  // The mirror terms hold a projection state only once in a Coulomb
  // integral, with three of orbitals: only the integrals with exactly one
  // new state are required. One entry for each nk_index: 1 for a new state
  std::vector<char> is_new_state(
    std::size_t(std::numeric_limits<DiracSpinor::Index>::max()) + 1, 0);
  for (const auto &x : new_states) {
    is_new_state[x.nk_index()] = 1;
  }
  const auto one_new_state =
    [&is_new_state](int, const DiracSpinor &a, const DiracSpinor &b,
                    const DiracSpinor &c, const DiracSpinor &d) {
      return is_new_state[a.nk_index()] + is_new_state[b.nk_index()] +
               is_new_state[c.nk_index()] + is_new_state[d.nk_index()] ==
             1;
    };
  qk->fill_if(orbitals_and_x, yk, one_new_state, max_k, false);
  return qk->count() - count_before;
}

} // namespace MBPT
