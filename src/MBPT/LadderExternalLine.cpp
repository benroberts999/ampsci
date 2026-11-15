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
    diagrams.L2_L3 ?
      L2_particle_hole(excited, core, core, ph_pairs, qk, &lk, sixj, &line) :
      std::vector<Matrix>(ph_pairs.num_blocks());
  const auto L2K_mb_na =
    diagrams.L2_L3 ?
      L2_particle_hole(excited, core, core, ph_pairs, qk, &lk, sixj) :
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
        const auto parity_xb = pair_parity(x, b);
        const auto L2K_for_L2 =
          read_L2K(L2K_mb_na_line, ph_pairs, parity_xb, index_x, index_b,
                   index_n, index_a, K_min, K_max)
            .first;
        const auto L2K_for_L3 =
          read_L2K(L2K_mb_na, ph_pairs, parity_xb, index_x, index_b, index_n,
                   index_a, K_min, K_max)
            .second;

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
GMatrix Sigma_ladder(const DiracSpinor &v, const ExternalLineLadder &ladder,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, double r0,
                     double rmax, std::size_t stride, bool include_G) {
  const auto grid = v.grid_sptr();
  const auto i0 = grid->getIndex(r0);
  const auto size = (grid->getIndex(rmax) - i0) / stride + 1;
  GMatrix Sigma{i0, stride, size, include_G, grid};
  std::vector<GMatrix> Sigma_thread(std::size_t(omp_get_max_threads()), Sigma);
  const auto max_k = ladder.L_mnxb.max_k();

#pragma omp parallel for schedule(dynamic)
  for (std::size_t index_n = 0; index_n < excited.size(); ++index_n) {
    const auto &n = excited[index_n];
    auto &Sigma_t = Sigma_thread[std::size_t(omp_get_thread_num())];
    // sum_x L_x |x>, for one set of the other indices
    DiracSpinor L_projected(0, v.kappa(), grid);
    for (const auto &a : core) {

      // |W^k_{.amn}> L^k_{mnxa} <x| / ([k][j_v] de)
      for (const auto &m : excited) {
        const auto [k_min, k_max] = k_minmax_L(m, n, v, a);
        const auto de = v.en() + a.en() - m.en() - n.en();
        for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
          L_projected *= 0.0;
          for (const auto &x : ladder.x_orbitals) {
            L_projected += ladder.L_mnxb.Q(k, m, n, x, a) * x;
          }
          Sigma_t.add(Coulomb::Wkv_bcd(k, v.kappa(), a, m, n), L_projected,
                      1.0 / ((2 * k + 1) * (v.twoj() + 1) * de));
        }
      }

      // |W^k_{.nab}> L^k_{xnab} <x| / ([k][j_v] de)
      for (const auto &b : core) {
        const auto [k_min, k_max] = k_minmax_L(v, n, a, b);
        const auto de = v.en() + n.en() - a.en() - b.en();
        for (int k = k_min; k <= std::min(k_max, max_k); ++k) {
          L_projected *= 0.0;
          for (const auto &x : ladder.x_orbitals) {
            L_projected += ladder.L_xnab.Q(k, x, n, a, b) * x;
          }
          Sigma_t.add(Coulomb::Wkv_bcd(k, v.kappa(), n, a, b), L_projected,
                      1.0 / ((2 * k + 1) * (v.twoj() + 1) * de));
        }
      }
    }
  }
  for (const auto &Sigma_t : Sigma_thread) {
    Sigma += Sigma_t;
  }
  return Sigma.drj_in_place();
}

} // namespace MBPT
