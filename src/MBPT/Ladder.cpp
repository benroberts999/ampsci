#include "Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <vector>

namespace MBPT {

// Size of the stack arrays used to cache the k-dependent (Q,L) integrals in
// the L1/L2/L4 inner loops, indexed directly by multipolarity k. Comfortably
// larger than any physical k = k_minmax_Q(...) for realistic bases.
constexpr std::size_t sk_array_size = 32;

//==============================================================================
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk) {

  const auto L123 = L1(k, m, n, i, j, qk, excited, SJ, Lk) +
                    L2(k, m, n, i, j, qk, core, excited, SJ, Lk) +
                    L3(k, m, n, i, j, qk, core, excited, SJ, Lk);
  // Optionally include "4th" ladder diagram
  // nb: L4 not fully checked!
  if (include_L4)
    return L123 + L4(k, m, n, i, j, qk, core, SJ, Lk);
  else
    return L123;
}

//------------------------------------------------------------------------------
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk) {

  // m (and n) must be excited states, as should 'excited'
  // Therefore, can test:
  // Ensured 'excited' is actually the excited orbitals
  // and that m and n are excited orbitals
  // (Still possible BOTH wrong at the same time..)
  // assert(std::find(excited.cbegin(), excited.cend(), m) != excited.cend());
  // assert(std::find(excited.cbegin(), excited.cend(), n) != excited.cend());

  double l1 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnij1 =
    Angular::neg1pow_2(2 + m.twoj() + n.twoj() + i.twoj() + j.twoj());

  //  6j(r) Triads: {m,i,k}, {k,l,u}, {i,l,r}, {u,r,m}
  //  6j(s) Triads: {n,b,k}, {k,l,u}, {b,l,s}, {u,s,n}
  // if (!Coulomb::sixjTriads(m, i, k, {}, {}, {}) ||
  //     !Coulomb::sixjTriads(n, j, k, {}, {}, {})) {
  //   std::cout << "!!\n";
  //   return 0.0;
  // }

  // Local dense cache of the recoupling 6j symbols. sj_r = {m,i,k;l,u,r} and
  // sj_s = {n,j,k;l,u,s} depend on the intermediate orbital (r or s) only
  // through its 2j, so two orbitals of the same kappa give the same value.
  // Precompute once per call, indexed by (2j, l, u), turning the inner-loop
  // SixJTable hash lookup into a direct array read. Bit-identical.
  // (kmax bounds every l,u from k_minmax_Q: l,u <= (max_ext_2j + max_2j)/2.)
  const int max_2j = DiracSpinor::max_tj(excited);
  const int kmax =
    (std::max({m.twoj(), n.twoj(), i.twoj(), j.twoj()}) + max_2j) / 2;
  // sj_dim = kmax+1; sj_cache indexed as [2j][l][u]
  const auto sj_dim = std::size_t(kmax) + 1;
  const auto sj_cache_index = [sj_dim](int t2, int l, int u) {
    return (std::size_t(t2) * sj_dim + std::size_t(l)) * sj_dim +
           std::size_t(u);
  };
  static thread_local std::vector<double> sjr_cache, sjs_cache;
  sjr_cache.resize((std::size_t(max_2j) + 1) * sj_dim * sj_dim);
  sjs_cache.resize((std::size_t(max_2j) + 1) * sj_dim * sj_dim);
  for (int t2 = 1; t2 <= max_2j; t2 += 2) {
    for (int l = 0; l <= kmax; ++l) {
      for (int u = 0; u <= kmax; ++u) {
        sjr_cache[sj_cache_index(t2, l, u)] =
          SJ.get_2(m.twoj(), i.twoj(), 2 * k, 2 * l, 2 * u, t2);
        sjs_cache[sj_cache_index(t2, l, u)] =
          SJ.get_2(n.twoj(), j.twoj(), 2 * k, 2 * l, 2 * u, t2);
      }
    }
  }

  // Thread-local cache of Q^u_{mnrs} for the current (m,n), indexed by
  // (r,s,u). Q^u_{mnrs} depends only on (m,n,r,s) - not on the external (i,j) -
  // but L1 is called once per (m,n,i,j). fill()/update() hold (m,n)=(a,b) fixed
  // while sweeping the inner (i,j)=(c,d), so caching here turns the dominant
  // cross-call Qk hash lookups into array reads (rebuilt only when (m,n) or the
  // excited set changes). The u-dimension bound depends only on (m,n), so the
  // cache stays valid across the whole (i,j) sweep. Bit-identical.
  const auto excited_size = excited.size();
  const int u_max = (std::max(m.twoj(), n.twoj()) + max_2j) / 2;
  // u_stride = u_max+1; upper bound on u from k_minmax_Q(m,n,r,s)
  const auto u_stride = std::size_t(u_max) + 1;
  const auto Qu_mnrs_index = [excited_size, u_stride](std::size_t ir,
                                                      std::size_t is, int u) {
    return (ir * excited_size + is) * u_stride + std::size_t(u);
  };
  static thread_local std::vector<double> Qu_mnrs;
  // Pointer comparisons detect if the basis vector is reallocated (basis change).
  static thread_local const DiracSpinor *prev_excited = nullptr;
  static thread_local std::size_t prev_excited_size = 0;
  static thread_local int prev_u_max = -1;
  static thread_local DiracSpinor::Index prev_m = 0, prev_n = 0;
  if (prev_excited != excited.data() || prev_excited_size != excited_size ||
      prev_u_max != u_max || prev_m != m.nk_index() || prev_n != n.nk_index()) {
    // resize (not assign): every cell the main loop reads (u in [u0,uI]) is
    // written below for the same (r,s); cells outside are never read.
    Qu_mnrs.resize(excited_size * excited_size * u_stride);
    for (auto ir = 0ul; ir < excited_size; ++ir) {
      for (auto is = 0ul; is < excited_size; ++is) {
        const auto [u0, uI] =
          Coulomb::k_minmax_Q(m, n, excited[ir], excited[is]);
        if (uI < u0)
          continue;
        const auto Qkey_mnrs = qk.NormalOrder(m, n, excited[ir], excited[is]);
        for (auto u = u0; u <= uI; u += 2) {
          Qu_mnrs[Qu_mnrs_index(ir, is, u)] = qk.Q(u, Qkey_mnrs);
        }
      }
    }
    prev_excited = excited.data();
    prev_excited_size = excited_size;
    prev_u_max = u_max;
    prev_m = m.nk_index();
    prev_n = n.nk_index();
  }

  for (auto ir = 0ul; ir < excited_size; ++ir) {
    const auto &r = excited[ir];
    for (auto is = 0ul; is < excited_size; ++is) {
      const auto &s = excited[is];

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_e_ijrs = 1.0 / (i.en() + j.en() - r.en() - s.en());
      const auto Qkey_rsij = qk.NormalOrder(r, s, i, j);
      const auto Lkey_rsij = Lk ? Lk->NormalOrder(r, s, i, j) : 0ul;

      // Cache (Q+L)^l_rsij: depends only on l, avoid N_u lookups
      assert(lI < int(sk_array_size));
      std::array<double, sk_array_size> QLl_rsij{};
      for (auto l = l0; l <= lI; l += 2) {
        QLl_rsij[std::size_t(l)] =
          qk.Q(l, Qkey_rsij) + (Lk ? Lk->Q(l, Lkey_rsij) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = Qu_mnrs[Qu_mnrs_index(ir, is, u)];
        // Zero when parity or triangle selection rules forbid this u.
        if (Q_umnrs == 0.0)
          continue;

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(u, r, m) == 0 || Coulomb::triangle(u, s, n) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto sj_r = sjr_cache[sj_cache_index(r.twoj(), l, u)];
          const auto sj_s = sjs_cache[sj_cache_index(s.twoj(), l, u)];

          const auto QL_lrsij = QLl_rsij[std::size_t(l)];

          l1 += (s_rs * sj_r * sj_s) * Q_umnrs * QL_lrsij * inv_e_ijrs;
        }
      }
    }
  }
  l1 *= s_mnij1 * tkp1;
  return l1;
}

//------------------------------------------------------------------------------
double L4(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk) {

  // m (and n) must be excited states, as should 'excited'
  // Therefore, can test:
  // Ensured 'excited' is actually the excited orbitals
  // and that m and n are excited orbitals
  // assert(std::find(core.cbegin(), core.cend(), m) == core.cend());
  // assert(std::find(core.cbegin(), core.cend(), n) == core.cend());

  double l4 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnij1 =
    Angular::neg1pow_2(2 + m.twoj() + n.twoj() + i.twoj() + j.twoj());

  //  6j(r) Triads: {m,i,k}, {k,u,l}, {i,u,c}, {l,c,m}
  //  6j(s) Triads: {n,b,k}, {k,u,l}, {b,u,d}, {l,d,n}

  // Local dense cache of the recoupling 6j symbols (see L1). sj_c = {m,i,k;u,l,c}
  // and sj_d = {n,j,k;u,l,d} depend on the intermediate orbital only through its
  // 2j; precompute once per call, indexed by (2j, u, l). c, d are core orbitals.
  // Bit-identical.
  const int max_2j = DiracSpinor::max_tj(core);
  const int kmax =
    (std::max({m.twoj(), n.twoj(), i.twoj(), j.twoj()}) + max_2j) / 2;
  // sj_dim = kmax+1; sj_cache indexed as [2j][u][l]
  const auto sj_dim = std::size_t(kmax) + 1;
  const auto sj_cache_index = [sj_dim](int t2, int u, int l) {
    return (std::size_t(t2) * sj_dim + std::size_t(u)) * sj_dim +
           std::size_t(l);
  };
  static thread_local std::vector<double> sjc_cache, sjd_cache;
  sjc_cache.resize((std::size_t(max_2j) + 1) * sj_dim * sj_dim);
  sjd_cache.resize((std::size_t(max_2j) + 1) * sj_dim * sj_dim);
  for (int t2 = 1; t2 <= max_2j; t2 += 2) {
    for (int u = 0; u <= kmax; ++u) {
      for (int l = 0; l <= kmax; ++l) {
        sjc_cache[sj_cache_index(t2, u, l)] =
          SJ.get_2(m.twoj(), i.twoj(), 2 * k, 2 * u, 2 * l, t2);
        sjd_cache[sj_cache_index(t2, u, l)] =
          SJ.get_2(n.twoj(), j.twoj(), 2 * k, 2 * u, 2 * l, t2);
      }
    }
  }

  // Thread-local cache of Q^u_{c,d,i,j} for fixed (i,j): avoids repeated hash
  // lookups when L4 is called repeatedly with the same (i,j).
  // Indexed as [ic * core_size + id][u].
  const auto core_size = core.size();
  const int u_max = (std::max(i.twoj(), j.twoj()) + max_2j) / 2;
  // u_stride = u_max+1; upper bound on u from k_minmax_Q(c,d,i,j)
  const auto u_stride = std::size_t(u_max) + 1;
  const auto Qu_cdij_index = [core_size, u_stride](std::size_t ic,
                                                   std::size_t id, int u) {
    return (ic * core_size + id) * u_stride + std::size_t(u);
  };
  static thread_local std::vector<double> Qu_cdij;
  // Pointer comparison detects if the core basis is reallocated (basis change).
  static thread_local const DiracSpinor *prev_core = nullptr;
  static thread_local DiracSpinor::Index prev_i = 0, prev_j = 0;
  static thread_local int prev_u_max = -1;
  if (prev_core != core.data() || prev_i != i.nk_index() ||
      prev_j != j.nk_index() || prev_u_max != u_max) {
    Qu_cdij.assign(core_size * core_size * u_stride, 0.0);
    for (auto ic = 0ul; ic < core_size; ++ic) {
      for (auto id = 0ul; id < core_size; ++id) {
        const auto [u0, uI] = Coulomb::k_minmax_Q(core[ic], core[id], i, j);
        if (uI < u0)
          continue;
        const auto Qkey_cdij = qk.NormalOrder(core[ic], core[id], i, j);
        for (auto u = u0; u <= uI; u += 2) {
          Qu_cdij[Qu_cdij_index(ic, id, u)] = qk.Q(u, Qkey_cdij);
        }
      }
    }
    prev_core = core.data();
    prev_i = i.nk_index();
    prev_j = j.nk_index();
    prev_u_max = u_max;
  }

  for (auto ic = 0ul; ic < core_size; ++ic) {
    const auto &c = core[ic];
    for (auto id = 0ul; id < core_size; ++id) {
      const auto &d = core[id];

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, d, i, j);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, c, d);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
      const auto inv_e_cdmn = 1.0 / (c.en() + d.en() - m.en() - n.en());
      const auto Qkey_mncd = qk.NormalOrder(m, n, c, d);
      const auto Lkey_mncd = Lk ? Lk->NormalOrder(m, n, c, d) : 0ul;

      // Cache (Q+L)^l_mncd: depends only on l, used inside the u loop
      assert(lI < int(sk_array_size));
      std::array<double, sk_array_size> QLl_mncd{};
      for (auto l = l0; l <= lI; l += 2) {
        QLl_mncd[std::size_t(l)] =
          qk.Q(l, Qkey_mncd) + (Lk ? Lk->Q(l, Lkey_mncd) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucdij = Qu_cdij[Qu_cdij_index(ic, id, u)];
        // Zero when parity or triangle selection rules forbid this u.
        if (Q_ucdij == 0.0)
          continue;

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(i, u, c) == 0 || Coulomb::triangle(j, u, d) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, u, l) == 0)
            continue;

          const auto sj_c = sjc_cache[sj_cache_index(c.twoj(), u, l)];
          const auto sj_d = sjd_cache[sj_cache_index(d.twoj(), u, l)];

          const auto QL_lmncd = QLl_mncd[std::size_t(l)];

          l4 += (s_cd * sj_c * sj_d) * Q_ucdij * QL_lmncd * inv_e_cdmn;
        }
      }
    }
  }
  l4 *= s_mnij1 * tkp1;
  return l4;
}

//------------------------------------------------------------------------------
double L2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk) {

  // m (and n) must be excited states, as should 'excited'
  // Therefore, can test:
  // Ensured 'excited' is actually the excited orbitals
  // and that m and n are excited orbitals
  // assert(std::find(excited.cbegin(), excited.cend(), m) != excited.cend());
  // assert(std::find(excited.cbegin(), excited.cend(), n) != excited.cend());
  // assert(std::find(core.cbegin(), core.cend(), m) == core.cend());

  double l2 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnijk =
    Angular::neg1pow_2(2 * k + m.twoj() + n.twoj() + i.twoj() + j.twoj());
  const auto ejm = j.en() - m.en();

  // Local dense cache of the recoupling 6j symbols (see L1). sj_c = {m,i,k;u,l,c}
  // and sj_r = {j,n,k;u,l,r} depend on the intermediate orbital only through its
  // 2j; precompute once per call, indexed by (2j, u, l). c is a core, r an
  // excited orbital, so max_2j spans both. Bit-identical.
  const int max_2j =
    std::max(DiracSpinor::max_tj(core), DiracSpinor::max_tj(excited));
  const int kmax =
    (std::max({m.twoj(), n.twoj(), i.twoj(), j.twoj()}) + max_2j) / 2;
  // sj_dim = kmax+1; sj_cache indexed as [2j][u][l]
  const auto sj_dim = std::size_t(kmax) + 1;
  const auto sj_cache_index = [sj_dim](int t2, int u, int l) {
    return (std::size_t(t2) * sj_dim + std::size_t(u)) * sj_dim +
           std::size_t(l);
  };
  static thread_local std::vector<double> sjc_cache, sjr_cache;
  sjc_cache.resize((std::size_t(max_2j) + 1) * sj_dim * sj_dim);
  sjr_cache.resize((std::size_t(max_2j) + 1) * sj_dim * sj_dim);
  for (int t2 = 1; t2 <= max_2j; t2 += 2) {
    for (int u = 0; u <= kmax; ++u) {
      for (int l = 0; l <= kmax; ++l) {
        sjc_cache[sj_cache_index(t2, u, l)] =
          SJ.get_2(m.twoj(), i.twoj(), 2 * k, 2 * u, 2 * l, t2);
        sjr_cache[sj_cache_index(t2, u, l)] =
          SJ.get_2(j.twoj(), n.twoj(), 2 * k, 2 * u, 2 * l, t2);
      }
    }
  }

  // (n,i) change on every call so a cross-call cache for Q^u_{cnir} would never
  // hit; look up inline (only for pairs surviving selection rules).
  const auto core_size = core.size();
  const auto excited_size = excited.size();

  for (auto ir = 0ul; ir < excited_size; ++ir) {
    const auto &r = excited[ir];
    for (auto ic = 0ul; ic < core_size; ++ic) {
      const auto &c = core[ic];

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, i, r);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
      const auto inv_e_cjmr = 1.0 / (c.en() + ejm - r.en());
      const auto Qkey_cnir = qk.NormalOrder(c, n, i, r);
      const auto Qkey_mrcj = qk.NormalOrder(m, r, c, j);
      const auto Lkey_mrcj = Lk ? Lk->NormalOrder(m, r, c, j) : 0ul;

      // Cache (Q+L)^l_mrcj: depends only on l, used inside the u loop (see L1).
      assert(lI < int(sk_array_size));
      std::array<double, sk_array_size> QLl_mrcj{};
      for (auto l = l0; l <= lI; l += 2) {
        QLl_mrcj[std::size_t(l)] =
          qk.Q(l, Qkey_mrcj) + (Lk ? Lk->Q(l, Lkey_mrcj) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucnir = qk.Q(u, Qkey_cnir);
        // Zero when parity or triangle selection rules forbid this u.
        if (Q_ucnir == 0.0)
          continue;

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(i, u, c) == 0 || Coulomb::triangle(n, u, r) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto s_ul = Angular::neg1pow(u + l);
          const auto sj_c = sjc_cache[sj_cache_index(c.twoj(), u, l)];
          const auto sj_r = sjr_cache[sj_cache_index(r.twoj(), u, l)];
          const auto QL_lmrcj = QLl_mrcj[std::size_t(l)];

          l2 += (s_ul * s_rc * sj_c * sj_r) * Q_ucnir * QL_lmrcj * inv_e_cjmr;
        }
      }
    }
  }
  l2 *= s_mnijk * tkp1;
  return l2;
}

//==============================================================================
void fill_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                  const Angular::SixJTable &sjt, bool print) {

  // const double a_damp = 0.35;
  // const double b_damp = 1.0 - a_damp;

  // Build combined basis: excited + core + any extra i_orbs (e.g. valence)
  const auto basis = qip::merge(core, excited);

  // Flat arrays for O(1) membership checks (no hash overhead).
  // Indexed directly by nk_index (uint16_t), so size is max_idx+1.
  DiracSpinor::Index max_idx = 0;
  for (const auto &x : basis)
    max_idx = std::max(max_idx, x.nk_index());
  // i_orbs is almost always a subset of basis - but not always
  // For example, basis might restict n_min_core
  // whereas i_orbs may be full basis, eg.g., L|i><i| = L_mnib|i>
  for (const auto &x : i_orbs)
    max_idx = std::max(max_idx, x.nk_index());
  std::vector<uint8_t> is_excited(std::size_t(max_idx) + 1, 0);
  std::vector<uint8_t> is_core(std::size_t(max_idx) + 1, 0);
  std::vector<uint8_t> is_i_orb(std::size_t(max_idx) + 1, 0);
  for (const auto &n : excited)
    is_excited[n.nk_index()] = 1;
  for (const auto &a : core)
    is_core[a.nk_index()] = 1;
  for (const auto &x : i_orbs)
    is_i_orb[x.nk_index()] = 1;

  const auto kmax = qk.max_k();

  // Base selection rule: m,n in excited, i in i_orbs, b in core + angular SR.
  // Determines which L^k_mnib integrals we actually need.
  const auto Lk_SR_one = [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                             const DiracSpinor &i,
                             const DiracSpinor &b) -> bool {
    // Require m and n to be excited
    if (!is_excited[m.nk_index()] || !is_excited[n.nk_index()])
      return false;
    // Require i to be in {i}, and b to be in core
    if (!is_i_orb[i.nk_index()] || !is_core[b.nk_index()])
      return false;
    const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, i, b);
    return k >= k0 && k <= kI;
  };

  // Selection rule passed to fill(): must be invariant under the Lk table
  // symmetry L^k_{abcd} = L^k_{badc}. fill() only stores/computes the canonical
  // tuple, so an entry must be accepted if EITHER it OR its symmetry partner
  // (b,a,d,c) passes the base rule - otherwise needed entries whose canonical
  // form puts a non-core hole in the b-slot (e.g. L^k_{mnva} -> canonical
  // (n,m,a,v) with valence v in slot b) are silently dropped.
  const auto Lk_SR = [&](int k, const DiracSpinor &a, const DiracSpinor &b,
                         const DiracSpinor &c, const DiracSpinor &d) -> bool {
    return Lk_SR_one(k, a, b, c, d) || Lk_SR_one(k, b, a, d, c);
  };

  // Lk integral
  const auto Lk_function = [&](int k, const DiracSpinor &m,
                               const DiracSpinor &n, const DiracSpinor &i,
                               const DiracSpinor &b) -> double {
    return Lkmnij(k, m, n, i, b, qk, core, excited, include_L4, sjt, nullptr);
  };

  lk->fill(basis, Lk_function, Lk_SR, kmax, print);
}

//==============================================================================
GMatrix Sigma_ladder(int kappa_v, double en_v,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &basis,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     const Angular::SixJTable &sjt, bool include_L4, double r0,
                     double rmax, std::size_t stride) {

  // Ladder correction to the correlation potential, evaluated at energy en_v.
  // The exchange is folded into the Coulomb vertex via W = Q + P (mirrors
  // de_valence_w), so no ladder-P is needed. The bra index i runs over the
  // excited basis states of kappa_v (projection basis, approximating
  // completeness).
  //
  // Diagrams (a)+(b)  [particle-particle]:
  //   Sigma_L += sum_{i,amn,k} |W^k_{.amn}> w <i| ,  w = L^k_{mn,i,a}/([k][j_v]de)
  //   de = en_v + e_a - e_m - e_n
  // Diagrams (c)+(d)  [particle-hole, external line in the m-slot of L]:
  //   Sigma_L += sum_{i,nab,k} |W^k_{.nab}> w <i| ,  w = L^k_{i,n,a,b}/([k][j_v]de)
  //   de = en_v + e_n - e_a - e_b
  //
  // The ladder integrals are always computed on-the-fly via Lkmnij() (using the
  // converged ladder table lk as the internal rung), evaluated at the fixed
  // external energy en_v. We cannot reuse the stored lk entries: those were
  // evaluated at the orbital energy e_i, not en_v (equal only for i = the
  // valence state at en_v = e_valence). The energy fix is applied by passing a
  // copy of the projection orbital with its energy set to en_v: the energy only
  // enters the L denominators, never the integral lookups, so this overrides
  // every denominator uniformly regardless of which slot the external line is
  // in (i-slot for a+b, m-slot for c+d).

  // 'lk' is the internal-rung ladder table passed straight to Lkmnij: pass
  // nullptr for L(Q,Q) = L^(1) (matches an un-iterated table), or a converged
  // table (its fixed point) for the full ladder.

  const auto tjv = Angular::twoj_k(kappa_v);

  // Sub-grid + GMatrix (no lower g part):
  const auto grid = excited.empty() ?
                      (core.empty() ? nullptr : core.front().grid_sptr()) :
                      excited.front().grid_sptr();
  const auto i0 = grid ? grid->getIndex(r0) : 0ul;
  const auto size = grid ? (grid->getIndex(rmax) - i0) / stride + 1 : 0ul;
  GMatrix Sd{i0, stride, size, true, grid};

  if (core.empty() || excited.empty())
    return Sd;

  // Projection basis: states with kappa == kappa_v
  std::vector<const DiracSpinor *> proj;
  for (const auto &x : basis) {
    if (x.kappa() == kappa_v)
      proj.push_back(&x);
  }
  if (proj.empty())
    return Sd;

  std::vector<GMatrix> Sd_ts(std::size_t(omp_get_max_threads()), Sd);

  for (auto ii = 0ul; ii < proj.size(); ++ii) {
    const auto &i = *proj[ii];

    // Copy of the projection orbital with energy fixed to en_v, used for the
    // external line in the ladder integrals (the bra <i| keeps the true energy).
    auto i_ev = i;
    i_ev.en() = en_v;

    qip::ProgressBar bar(core.size());
#pragma omp parallel for schedule(dynamic)
    for (auto ia = 0ul; ia < core.size(); ++ia) {
      const auto &a = core[ia];
      // Per-thread accumulator: must be bound INSIDE the parallel region so
      // each thread writes to its own matrix (else all threads race on one).
      auto &Sd_t = Sd_ts[std::size_t(omp_get_thread_num())];
      for (const auto &n : excited) {

        // Diagrams (a)+(b): W^k_{.amn} L^k_{mn,i,a}
        const auto [kmin_na, kmax_na] = Coulomb::k_minmax_Ck(n, a);
        for (int k = kmin_na; k <= kmax_na; ++k) {

          const auto f_kkjj = (2 * k + 1) * (tjv + 1);
          // Enforce the (a,n) Coulomb parity: only k with Q^k_{vamn} != 0
          // contribute. The W=Q+P ket does NOT self-gate parity the way the
          // bare Q ket does, so without this its P-part adds spurious
          // wrong-parity terms (not present in de_valence_w).
          if (!Angular::Ck_kk_SR(k, n.kappa(), a.kappa()))
            continue;

          for (const auto &m : excited) {
            if (!Angular::Ck_kk_SR(k, kappa_v, m.kappa()))
              continue;
            const auto Lkmnia =
              Lkmnij(k, m, n, i_ev, a, qk, core, excited, include_L4, sjt, lk);
            if (Lkmnia == 0.0)
              continue;
            const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, a, m, n);
            const auto dele = en_v + a.en() - m.en() - n.en();
            Sd_t.add(Wkv, i, Lkmnia / (f_kkjj * dele));
          }
        }

        // Diagrams (c)+(d): W^k_{.nab} L^k_{i,n,a,b}
        for (const auto &b : core) {
          const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax_Ck(n, b);
          for (int k = kmin_nb; k <= kmax_nb; ++k) {

            const auto f_kkjj = (2 * k + 1) * (tjv + 1);
            if (!Angular::Ck_kk_SR(k, kappa_v, a.kappa()))
              continue;
            // Enforce the (n,b) Coulomb parity (only k with Q^k_{vnab} != 0);
            // see note in the (a)+(b) block.
            if (!Angular::Ck_kk_SR(k, n.kappa(), b.kappa()))
              continue;

            // Always calculate? For some values of i_ev, might already have!
            // Energy denominator??
            const auto Lkinab =
              Lkmnij(k, i_ev, n, a, b, qk, core, excited, include_L4, sjt, lk);
            if (Lkinab == 0.0)
              continue;
            const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, n, a, b);
            const auto dele = en_v + n.en() - a.en() - b.en();
            Sd_t.add(Wkv, i, Lkinab / (f_kkjj * dele));
          }
        }
      }
      bar.update();
    }
  }

  for (const auto &Sd_t : Sd_ts) {
    Sd += Sd_t;
  }

  return Sd.drj_in_place();
}

//==============================================================================
void update_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                    const std::vector<DiracSpinor> &excited,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &i_orbs, bool include_L4,
                    const Angular::SixJTable &sjt,
                    const Coulomb::LkTable *const lk_prev, double a_damp,
                    bool print) {

  // Build combined basis: excited + core + any extra i_orbs (e.g. valence)
  // i_orbs not required (either core or subset of excited)
  // Used previously for "selection rules" - now SR is if we already calc'd it!
  (void)i_orbs;
  const auto basis = qip::merge(core, excited);

  // Lk integral with damping folded in
  const auto Lk_function = [&](int k, const DiracSpinor &m,
                               const DiracSpinor &n, const DiracSpinor &i,
                               const DiracSpinor &b) -> double {
    return Lkmnij(k, m, n, i, b, qk, core, excited, include_L4, sjt, lk_prev);
  };

  lk->update(basis, Lk_function, a_damp, print);
}

} // namespace MBPT
