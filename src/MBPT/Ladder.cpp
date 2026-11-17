#include "Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/LadderExternalLine.hpp"
#include "MBPT/LadderIteration.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include "qip/Vector.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace MBPT {

//==============================================================================
SigmaLMethod parseSigmaLMethod(const std::string &method) {
  if (qip::ci_compare(method, "basis") || qip::ci_compare(method, "full"))
    return SigmaLMethod::basis;
  if (qip::ci_compare(method, "Dzuba") || qip::ci_compare(method, "ratio"))
    return SigmaLMethod::ratio;
  if (qip::ci_compare(method, "ratio2"))
    return SigmaLMethod::ratio2;
  std::cout << "Warning: unknown Sigma_L method: " << method
            << " ?? Defaulting to basis\n";
  return SigmaLMethod::basis;
}

std::string parseSigmaLMethod(SigmaLMethod method) {
  switch (method) {
  case SigmaLMethod::basis:
    return "basis";
  case SigmaLMethod::ratio:
    return "ratio";
  case SigmaLMethod::ratio2:
    return "ratio2";
  }
  assert(false);
  return "ratio";
}

//==============================================================================
void fill_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs,
                  const LadderDiagrams &diagrams, int max_k, bool print) {

  // Required integrals: L^k_mnib with m,n excited, i in i_orbs, b in core.
  // nb: k runs over the full ladder range (k_minmax_L): triangle rules plus
  // the combined parity rule only. Unlike Q^k, L^k has components at BOTH
  // parities of k; the wrong-parity (Coulomb-forbidden) components enter the
  // exchange diagrams (via P/W) and feed back into the dressed integrals
  // (Q+L) of the next iteration.

  const auto kmax = (max_k < 0) ? qk.max_k() : std::min(max_k, qk.max_k());

  // Maximum multipolarity of the stored integrals: same convention as
  // CoulombTable::fill (kmax <= 0 means no cut-off)
  const auto basis_max_k = DiracSpinor::max_tj(qip::merge(core, excited));
  const auto max_k_stored =
    (kmax <= 0) ? basis_max_k : std::min(basis_max_k, kmax);

  IO::ChronoTimer timer("", print);

  // Add the required integrals that are not already in the table. Existing
  // entries are kept (e.g., read from file), and only the new ones calculated
  const auto num_existing = lk->count();
  add_Lk_mnib_entries(lk, excited, core, i_orbs, max_k_stored);
  if (lk->count() == num_existing)
    return;

  if (print) {
    fmt::print("Calculating {}{} integrals\n",
               num_existing == 0 ? "" : "(subset of) ", lk->count());
  }
  calculate_Lk_mnib(lk, qk, excited, core, i_orbs, diagrams, nullptr,
                    max_k_stored, true, 0.0);
}
//==============================================================================
GMatrix Sigma_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     double r0, double rmax, std::size_t stride,
                     bool include_G) {

  // Ladder correction to the correlation potential, evaluated at energy e_v,
  // by the ratio method (follows Dzuba, PRA 78, 042502 (2008)). Each term of
  // the regular second-order Sigma (as in Goldstone::Sigma_both) is rescaled
  // by the scalar ratio L/Q, both taken straight from the stored tables (the
  // L entries with i = v are already at the valence energy - nothing computed
  // on-the-fly; fast):
  //   (a+b): |Q^k_{.amn}><W^k_{.amn}| (L^k_{mnva}/Q^k_{mnva}) / ([k][j_v] de)
  //   (c+d): |Q^k_{.nab}><W^k_{.nab}| (L^k_{vnab}/Q^k_{vnab}) / ([k][j_v] de)
  // The diagonal then reproduces the Coulomb-parity part of the ladder
  // energy exactly; the radial shape is approximate (each term keeps its
  // second-order shape). Terms with Q^k = 0 cannot be rescaled: dropped. In
  // particular, the exchange-only k channels (Coulomb-forbidden parity,
  // where W = P and L^k != 0 - see k_minmax_L) are inherently absent, so
  // <v|Sigma_L|v> differs from de_valence_w(v) by those (small) terms.

  const auto kappa_v = v.kappa();
  const auto en_v = v.en();
  const auto tjv = v.twoj();

  // Sub-grid + GMatrix (lower g part included if include_G):
  const auto grid = v.grid_sptr();
  const auto i0 = grid->getIndex(r0);
  const auto size = (grid->getIndex(rmax) - i0) / stride + 1;
  GMatrix Sd{i0, stride, size, include_G, grid};

  if (core.empty() || excited.empty())
    return Sd;

  std::vector<GMatrix> Sd_ts(std::size_t(omp_get_max_threads()), Sd);

  qip::ProgressBar bar(excited.size());
#pragma omp parallel for schedule(dynamic)
  for (auto in = 0ul; in < excited.size(); ++in) {
    const auto &n = excited[in];
    // Per-thread accumulator: must be bound INSIDE the parallel region so
    // each thread writes to its own matrix (else all threads race on one).
    auto &Sd_t = Sd_ts[std::size_t(omp_get_thread_num())];
    for (const auto &a : core) {

      // Diagrams (a)+(b): |Q^k_{.amn}><W^k_{.amn}| L^k_{mnva}/Q^k_{mnva}
      // The second-order |Q><W| terms require Q^k != 0 (Coulomb parity)
      const auto [kmin_na, kmax_na] = Coulomb::k_minmax_tj(n.twoj(), a.twoj());
      for (int k = kmin_na; k <= kmax_na; ++k) {

        const auto f_kkjj = (2 * k + 1) * (tjv + 1);
        if (!Angular::Ck_kk_SR(k, n.kappa(), a.kappa()))
          continue;
        for (const auto &m : excited) {
          if (!Angular::Ck_kk_SR(k, kappa_v, m.kappa()))
            continue;
          const auto dele = en_v + a.en() - m.en() - n.en();
          const auto Qkmnva = qk.Q(k, m, n, v, a);
          if (Qkmnva == 0.0)
            continue;
          const auto Lkmnva = lk ? lk->Q(k, m, n, v, a) : 0.0;
          if (Lkmnva == 0.0)
            continue;
          const auto Qkv = Coulomb::Qkv_bcd(k, kappa_v, a, m, n);
          const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, a, m, n);
          Sd_t.add(Qkv, Wkv, (Lkmnva / Qkmnva) / (f_kkjj * dele));
        }
      }

      // Diagrams (c)+(d): |Q^k_{.nab}><W^k_{.nab}| L^k_{vnab}/Q^k_{vnab}
      for (const auto &b : core) {
        const auto [kmin_nb, kmax_nb] =
          Coulomb::k_minmax_tj(n.twoj(), b.twoj());
        for (int k = kmin_nb; k <= kmax_nb; ++k) {

          const auto f_kkjj = (2 * k + 1) * (tjv + 1);
          if (!Angular::Ck_kk_SR(k, n.kappa(), b.kappa()) ||
              !Angular::Ck_kk_SR(k, kappa_v, a.kappa()))
            continue;

          const auto dele = en_v + n.en() - a.en() - b.en();
          const auto Qkvnab = qk.Q(k, v, n, a, b);
          if (Qkvnab == 0.0)
            continue;
          const auto Lkvnab = lk ? lk->Q(k, v, n, a, b) : 0.0;
          if (Lkvnab == 0.0)
            continue;
          const auto Qkv = Coulomb::Qkv_bcd(k, kappa_v, n, a, b);
          const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, n, a, b);
          Sd_t.add(Qkv, Wkv, (Lkvnab / Qkvnab) / (f_kkjj * dele));
        }
      }
    }
    bar.update();
  }

  for (const auto &Sd_t : Sd_ts) {
    Sd += Sd_t;
  }

  return Sd.drj_in_place();
}

//==============================================================================

GMatrix Sigma_ladder_ratio2(const DiracSpinor &v,
                            const std::vector<DiracSpinor> &core,
                            const std::vector<DiracSpinor> &excited,
                            const Coulomb::QkTable &qk,
                            const Coulomb::LkTable &lk, double r0, double rmax,
                            std::size_t stride, bool include_G) {
  // Ratio method, made symmetric
  const auto Sigma_ratio =
    Sigma_ladder(v, core, excited, qk, &lk, r0, rmax, stride, include_G);
  auto Sigma = 0.5 * (Sigma_ratio + Sigma_ratio.transpose_drj());

  // Exact column Sigma_L|v>: the x = v term of the projection onto the basis,
  // whose ladder integrals L^k_mnva and L^k_vnab are in the ladder table
  const auto column = Sigma_ladder(v, std::vector{v}, lk, lk, core, excited, r0,
                                   rmax, stride, include_G, true) *
                      v;

  // Symmetric rank-2 correction, so that Sigma_L|v> is the exact column
  const auto c = column - Sigma * v;
  GMatrix correction{Sigma.i0(), stride, Sigma.size(), include_G,
                     v.grid_sptr()};
  correction.add(c, v);
  correction.add(v, c);
  correction.add(v, v, -(v * c));
  Sigma += correction.drj_in_place();
  return Sigma;
}

//==============================================================================
void update_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                    const std::vector<DiracSpinor> &excited,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &update_i,
                    const LadderDiagrams &diagrams,
                    const Coulomb::LkTable *const lk_prev, double a_damp,
                    bool print, const LinAlg::Matrix<double> *rho_ma) {

  // Restrict the re-iteration to entries whose i-slot is in 'update_i' (b is
  // always core; by the L^k_{mnib}=L^k_{nmbi} symmetry the i index sits in
  // slot 3 or 4). Lets us converge core (update_i=core) before valence
  // (update_i=valence), the latter using the now-converged core integrals.
  // All these integrals are calculated at once [see calculate_Lk_mnib()].
  assert(a_damp != 1.0);
  IO::ChronoTimer timer("", print);
  if (print) {
    fmt::print("Updating (subset of) {} integrals\n", lk->count());
  }
  calculate_Lk_mnib(lk, qk, excited, core, update_i, diagrams, lk_prev,
                    lk->max_k(), false, a_damp, nullptr, rho_ma);
}

//==============================================================================

RelaxationDensity relaxation_density(const Coulomb::QkTable &qk,
                                     const std::vector<DiracSpinor> &core,
                                     const std::vector<DiracSpinor> &excited) {
  RelaxationDensity density{
    LinAlg::Matrix<double>(core.size(), core.size()),
    LinAlg::Matrix<double>(excited.size(), excited.size())};
  if (core.empty() || excited.empty())
    return density;
  const auto max_tj =
    std::max(DiracSpinor::max_tj(excited), DiracSpinor::max_tj(core));
  const SixJArray sixj(max_tj, std::min(max_tj, ladder_max_k));

  // Core: gamma_cb, the core orbitals partly emptied
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t index_b = 0; index_b < core.size(); ++index_b) {
    for (std::size_t index_c = 0; index_c < core.size(); ++index_c) {
      const auto &b = core[index_b];
      const auto &c = core[index_c];
      if (c.kappa() != b.kappa())
        continue;
      double sum = 0.0;
      for (const auto &a : core) {
        for (const auto &m : excited) {
          for (const auto &n : excited) {
            const auto Q_abmn = coulomb_multipoles(a, b, m, n, qk);
            if (Q_abmn.k_max < Q_abmn.k_min)
              continue;
            const auto W_mnac = W_multipoles(m, n, a, c, qk, sixj);
            const auto de = (a.en() + b.en() - m.en() - n.en()) *
                            (a.en() + c.en() - m.en() - n.en());
            for (int k = Q_abmn.k_min; k <= Q_abmn.k_max; k += 2) {
              sum -= Q_abmn.value[std::size_t(k)] *
                     W_mnac.value[std::size_t(k)] / ((2 * k + 1) * de);
            }
          }
        }
      }
      density.gamma_cb(index_b, index_c) = sum / b.twojp1();
    }
  }

  // Excited: gamma_ms, the excited orbitals partly filled
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (std::size_t index_m = 0; index_m < excited.size(); ++index_m) {
    for (std::size_t index_s = 0; index_s < excited.size(); ++index_s) {
      const auto &m = excited[index_m];
      const auto &s = excited[index_s];
      if (s.kappa() != m.kappa())
        continue;
      double sum = 0.0;
      for (const auto &a : core) {
        for (const auto &b : core) {
          for (const auto &n : excited) {
            const auto Q_abmn = coulomb_multipoles(a, b, m, n, qk);
            if (Q_abmn.k_max < Q_abmn.k_min)
              continue;
            const auto W_nsba = W_multipoles(n, s, b, a, qk, sixj);
            const auto de = (a.en() + b.en() - m.en() - n.en()) *
                            (a.en() + b.en() - n.en() - s.en());
            for (int k = Q_abmn.k_min; k <= Q_abmn.k_max; k += 2) {
              sum += Q_abmn.value[std::size_t(k)] *
                     W_nsba.value[std::size_t(k)] / ((2 * k + 1) * de);
            }
          }
        }
      }
      density.gamma_ms(index_m, index_s) = sum / m.twojp1();
    }
  }
  return density;
}

//==============================================================================

std::array<double, 2>
de_valence_relaxation(const DiracSpinor &v, const Coulomb::QkTable &qk,
                      const LinAlg::Matrix<double> &rho_ma,
                      const RelaxationDensity &density,
                      const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &excited) {
  assert(rho_ma.rows() == excited.size() && rho_ma.cols() == core.size());
  assert(density.gamma_cb.rows() == core.size() &&
         density.gamma_ms.rows() == excited.size());
  const auto max_tj = std::max(
    {DiracSpinor::max_tj(excited), DiracSpinor::max_tj(core), v.twoj()});
  const SixJArray sixj(max_tj, std::min(max_tj, ladder_max_k));

  // Each term: f_vp W^0_vpvq gamma, over p, q of one kappa
  double rho_part = 0.0;
  for (std::size_t index_a = 0; index_a < core.size(); ++index_a) {
    const auto &a = core[index_a];
    const auto f_va = closed_shell_factor(v, a);
    for (std::size_t index_m = 0; index_m < excited.size(); ++index_m) {
      const auto &m = excited[index_m];
      if (m.kappa() != a.kappa())
        continue;
      rho_part += 2.0 * f_va * W_multipoles(v, a, v, m, qk, sixj).value[0] *
                  rho_ma(index_m, index_a);
    }
  }

  double occupation_part = 0.0;
  for (std::size_t index_b = 0; index_b < core.size(); ++index_b) {
    const auto &b = core[index_b];
    const auto f_vb = closed_shell_factor(v, b);
    for (std::size_t index_c = 0; index_c < core.size(); ++index_c) {
      const auto &c = core[index_c];
      if (c.kappa() != b.kappa())
        continue;
      occupation_part += f_vb * W_multipoles(v, b, v, c, qk, sixj).value[0] *
                         density.gamma_cb(index_b, index_c);
    }
  }
  for (std::size_t index_m = 0; index_m < excited.size(); ++index_m) {
    const auto &m = excited[index_m];
    const auto f_vm = closed_shell_factor(v, m);
    for (std::size_t index_s = 0; index_s < excited.size(); ++index_s) {
      const auto &s = excited[index_s];
      if (s.kappa() != m.kappa())
        continue;
      occupation_part += f_vm * W_multipoles(v, m, v, s, qk, sixj).value[0] *
                         density.gamma_ms(index_m, index_s);
    }
  }
  return {rho_part, occupation_part};
}

//==============================================================================
bool write_SigmaL(const std::string &fname, const std::vector<SigmaLData> &SLs,
                  const Grid &grid) {
  if (fname.empty() || fname == "false" || SLs.empty())
    return false;

  std::cout << "\nWriting Sigma_L (ladder) to file: " << fname << " ... "
            << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, IO::FRW::write);
  const auto rw = IO::FRW::write;

  // Full-grid parameters - just to check on read:
  auto r0 = grid.r0();
  auto rmax = grid.rmax();
  auto b = grid.loglin_b();
  auto pts = grid.num_points();
  rw_binary(iofs, rw, r0, rmax, b, pts);

  auto num = SLs.size();
  rw_binary(iofs, rw, num);

  for (const auto &sig : SLs) {
    auto kappa = sig.kappa;
    auto n = sig.n;
    auto en = sig.en;
    auto i0 = sig.SL.i0();
    auto stride = sig.SL.stride();
    auto size = sig.SL.size();
    auto incl_g = sig.SL.includes_g();
    rw_binary(iofs, rw, kappa, n, en, i0, stride, size, incl_g);
    for (auto i = 0ul; i < size; ++i) {
      for (auto j = 0ul; j < size; ++j) {
        auto ff = sig.SL.ff(i, j);
        rw_binary(iofs, rw, ff);
        if (incl_g) {
          auto fg = sig.SL.fg(i, j);
          auto gf = sig.SL.gf(i, j);
          auto gg = sig.SL.gg(i, j);
          rw_binary(iofs, rw, fg, gf, gg);
        }
      }
    }
  }
  std::cout << "done.\n";
  return true;
}

//==============================================================================
std::vector<SigmaLData> read_SigmaL(const std::string &fname,
                                    const std::shared_ptr<const Grid> &grid) {
  std::vector<SigmaLData> SLs;
  if (fname.empty() || fname == "false" || !grid)
    return SLs;
  if (!IO::FRW::file_exists(fname)) {
    std::cout << "\nNo Sigma_L (ladder) file: " << fname << "\n";
    return SLs;
  }

  std::cout << "\nReading Sigma_L (ladder) from file: " << fname << " ... "
            << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, IO::FRW::read);
  const auto rw = IO::FRW::read;

  // Full-grid parameters - must match current grid:
  double r0{}, rmax{}, b{};
  std::size_t pts{};
  rw_binary(iofs, rw, r0, rmax, b, pts);
  const bool grid_ok = std::abs((r0 - grid->r0()) / r0) < 1.0e-6 &&
                       std::abs(rmax - grid->rmax()) < 0.001 &&
                       std::abs(b - grid->loglin_b()) < 0.001 &&
                       pts == grid->num_points();
  if (!grid_ok) {
    std::cout << "\nCannot read from: " << fname << ". Grid mismatch.\n";
    return SLs;
  }

  std::size_t num{};
  rw_binary(iofs, rw, num);

  for (std::size_t iS = 0; iS < num; ++iS) {
    int kappa{}, n{};
    double en{};
    std::size_t i0{}, stride{}, size{};
    bool incl_g{};
    rw_binary(iofs, rw, kappa, n, en, i0, stride, size, incl_g);
    GMatrix SL{i0, stride, size, incl_g, grid};
    for (auto i = 0ul; i < size; ++i) {
      for (auto j = 0ul; j < size; ++j) {
        rw_binary(iofs, rw, SL.ff(i, j));
        if (incl_g) {
          rw_binary(iofs, rw, SL.fg(i, j), SL.gf(i, j), SL.gg(i, j));
        }
      }
    }
    SLs.push_back({kappa, n, en, std::move(SL)});
  }
  std::cout << "done.\n";
  return SLs;
}

} // namespace MBPT
