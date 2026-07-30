#include "Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/String.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

namespace MBPT {

//==============================================================================
SigmaLMethod parseSigmaLMethod(const std::string &method) {
  if (qip::ci_compare(method, "basis") || qip::ci_compare(method, "full"))
    return SigmaLMethod::basis;
  if (qip::ci_compare(method, "Dzuba") || qip::ci_compare(method, "ratio"))
    return SigmaLMethod::ratio;
  if (qip::ci_compare(method, "direct"))
    return SigmaLMethod::direct;
  std::cout << "Warning: unknown Sigma_L method: " << method
            << " ?? Defaulting to ratio\n";
  return SigmaLMethod::ratio;
}

std::string parseSigmaLMethod(SigmaLMethod method) {
  switch (method) {
  case SigmaLMethod::basis:
    return "basis";
  case SigmaLMethod::ratio:
    return "ratio";
  case SigmaLMethod::direct:
    return "direct";
  }
  assert(false);
  return "ratio";
}

//==============================================================================
// Size of the stack arrays used to cache the k-dependent (Q,L) integrals in
// the L1/L2/L4 inner loops, indexed directly by multipolarity k. Comfortably
// larger than any physical k = k_minmax_Q(...) for realistic bases.
constexpr std::size_t sk_array_size = 32;

namespace {
/*
Dense cache of the recoupling 6j symbols {tj1/2, tj2/2, k; a, b, t2/2},
indexed by (t2, a, b) with a,b in [0, kmax] - see L1/L2/L4. The values depend
on the external orbitals only through (tj1, tj2) and on k; the dimensions
depend on (kmax, max_2j). The cache is refilled only when this key changes, so
it is re-used across calls that vary only the other externals - e.g. the
(i,j) sweep in fill/update, or the projection-state sweep in Sigma_ladder.
Bit-identical to an unconditional refill.
Intended to be used as a (static) thread_local (one per call site).
*/
class SixJCache {
  std::vector<double> m_cache{};
  int m_tj1{-1}, m_tj2{-1}, m_k{-1}, m_max_2j{-1}, m_kmax{-1};
  const Angular::SixJTable *m_SJ{nullptr};
  std::size_t m_dim{0};

  std::size_t index(int t2, int a, int b) const {
    return (std::size_t(t2) * m_dim + std::size_t(a)) * m_dim + std::size_t(b);
  }

public:
  // Fills cache with {tj1/2, tj2/2, k; a, b, t2/2} for all (t2, a, b).
  // No-op if the key is unchanged since the last call.
  void update(int tj1, int tj2, int k, int max_2j, int kmax,
              const Angular::SixJTable &SJ) {
    if (tj1 == m_tj1 && tj2 == m_tj2 && k == m_k && max_2j == m_max_2j &&
        kmax == m_kmax && &SJ == m_SJ)
      return;
    m_tj1 = tj1;
    m_tj2 = tj2;
    m_k = k;
    m_max_2j = max_2j;
    m_kmax = kmax;
    m_SJ = &SJ;
    m_dim = std::size_t(kmax) + 1;
    m_cache.resize((std::size_t(max_2j) + 1) * m_dim * m_dim);
    for (int t2 = 1; t2 <= max_2j; t2 += 2) {
      for (int a = 0; a <= kmax; ++a) {
        for (int b = 0; b <= kmax; ++b) {
          m_cache[index(t2, a, b)] =
            SJ.get_2(tj1, tj2, 2 * k, 2 * a, 2 * b, t2);
        }
      }
    }
  }

  double get(int t2, int a, int b) const { return m_cache[index(t2, a, b)]; }
};
} // namespace

//==============================================================================
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk,
              std::optional<double> e_i, std::optional<double> e_m) {

  // nb: i energy enters only L1 and L3; m energy only L2 and L4
  const auto L123 = L1(k, m, n, i, j, qk, excited, SJ, Lk, e_i) +
                    L2(k, m, n, i, j, qk, core, excited, SJ, Lk, {}, e_m) +
                    L3(k, m, n, i, j, qk, core, excited, SJ, Lk, e_i);
  // Optionally include the hole-hole (L4) ladder diagram
  if (include_L4)
    return L123 + L4(k, m, n, i, j, qk, core, SJ, Lk, e_m);
  else
    return L123;
}

//------------------------------------------------------------------------------
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk,
          std::optional<double> e_i) {

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

  // Cached recoupling 6j symbols (see SixJCache). sj_r = {m,i,k;l,u,r} and
  // sj_s = {n,j,k;l,u,s} depend on the intermediate orbital (r or s) only
  // through its 2j, so two orbitals of the same kappa give the same value.
  // Indexed by (2j, l, u), turning the inner-loop SixJTable hash lookup into
  // a direct array read; refilled only when (2j pair, k, dims) change.
  // (kmax bounds every l,u from k_minmax_Q: l,u <= (max_ext_2j + max_2j)/2.)
  const int max_2j = DiracSpinor::max_tj(excited);
  const int kmax =
    (std::max({m.twoj(), n.twoj(), i.twoj(), j.twoj()}) + max_2j) / 2;
  static thread_local SixJCache sjr_cache, sjs_cache;
  sjr_cache.update(m.twoj(), i.twoj(), k, max_2j, kmax, SJ);
  sjs_cache.update(n.twoj(), j.twoj(), k, max_2j, kmax, SJ);

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
      // Rung range: the bare Coulomb part exists only in the parity range
      // [l0,lI]; the L^l part exists at both parities of l (see k_minmax_L)
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, j);
      const auto [w0, wI] = k_minmax_L(r, s, i, j);
      if (uI < u0 || wI < w0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_e_ijrs =
        1.0 / (e_i.value_or(i.en()) + j.en() - r.en() - s.en());
      const auto Qkey_rsij = qk.NormalOrder(r, s, i, j);
      const auto Lkey_rsij = Lk ? Lk->NormalOrder(r, s, i, j) : 0ul;

      // Cache (Q+L)^l_rsij: depends only on l, avoid N_u lookups
      assert(wI < int(sk_array_size));
      std::array<double, sk_array_size> QLl_rsij{};
      for (auto l = w0; l <= wI; ++l) {
        const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
        QLl_rsij[std::size_t(l)] =
          (has_Q ? qk.Q(l, Qkey_rsij) : 0.0) + (Lk ? Lk->Q(l, Lkey_rsij) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = Qu_mnrs[Qu_mnrs_index(ir, is, u)];
        // Zero when parity or triangle selection rules forbid this u.
        if (Q_umnrs == 0.0)
          continue;

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(u, r, m) == 0 || Coulomb::triangle(u, s, n) == 0)
          continue;

        for (auto l = w0; l <= wI; ++l) {

          const auto QL_lrsij = QLl_rsij[std::size_t(l)];
          if (QL_lrsij == 0.0)
            continue;

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto sj_r = sjr_cache.get(r.twoj(), l, u);
          const auto sj_s = sjs_cache.get(s.twoj(), l, u);

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
          const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk,
          std::optional<double> e_m) {

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

  // Cached recoupling 6j symbols (see SixJCache). sj_c = {m,i,k;u,l,c}
  // and sj_d = {n,j,k;u,l,d} depend on the intermediate orbital only through
  // its 2j; indexed by (2j, u, l); refilled only when (2j pair, k, dims)
  // change. c, d are core orbitals.
  const int max_2j = DiracSpinor::max_tj(core);
  const int kmax =
    (std::max({m.twoj(), n.twoj(), i.twoj(), j.twoj()}) + max_2j) / 2;
  static thread_local SixJCache sjc_cache, sjd_cache;
  sjc_cache.update(m.twoj(), i.twoj(), k, max_2j, kmax, SJ);
  sjd_cache.update(n.twoj(), j.twoj(), k, max_2j, kmax, SJ);

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
      // Rung range: L^l part exists at both parities of l (see L1)
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, c, d);
      const auto [w0, wI] = k_minmax_L(m, n, c, d);
      if (uI < u0 || wI < w0)
        continue;

      const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
      const auto inv_e_cdmn =
        1.0 / (c.en() + d.en() - e_m.value_or(m.en()) - n.en());
      const auto Qkey_mncd = qk.NormalOrder(m, n, c, d);
      const auto Lkey_mncd = Lk ? Lk->NormalOrder(m, n, c, d) : 0ul;

      // Cache (Q+L)^l_mncd: depends only on l, used inside the u loop
      assert(wI < int(sk_array_size));
      std::array<double, sk_array_size> QLl_mncd{};
      for (auto l = w0; l <= wI; ++l) {
        const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
        QLl_mncd[std::size_t(l)] =
          (has_Q ? qk.Q(l, Qkey_mncd) : 0.0) + (Lk ? Lk->Q(l, Lkey_mncd) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucdij = Qu_cdij[Qu_cdij_index(ic, id, u)];
        // Zero when parity or triangle selection rules forbid this u.
        if (Q_ucdij == 0.0)
          continue;

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(i, u, c) == 0 || Coulomb::triangle(j, u, d) == 0)
          continue;

        for (auto l = w0; l <= wI; ++l) {

          const auto QL_lmncd = QLl_mncd[std::size_t(l)];
          if (QL_lmncd == 0.0)
            continue;

          // 6j triad:
          if (Angular::triangle(k, u, l) == 0)
            continue;

          const auto sj_c = sjc_cache.get(c.twoj(), u, l);
          const auto sj_d = sjd_cache.get(d.twoj(), u, l);

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
          const Coulomb::LkTable *const Lk, std::optional<double> e_j,
          std::optional<double> e_m) {

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
  const auto ejm = e_j.value_or(j.en()) - e_m.value_or(m.en());

  // Cached recoupling 6j symbols (see SixJCache). sj_c = {m,i,k;u,l,c}
  // and sj_r = {j,n,k;u,l,r} depend on the intermediate orbital only through
  // its 2j; indexed by (2j, u, l); refilled only when (2j pair, k, dims)
  // change. c is a core, r an excited orbital, so max_2j spans both.
  const int max_2j =
    std::max(DiracSpinor::max_tj(core), DiracSpinor::max_tj(excited));
  const int kmax =
    (std::max({m.twoj(), n.twoj(), i.twoj(), j.twoj()}) + max_2j) / 2;
  static thread_local SixJCache sjc_cache, sjr_cache;
  sjc_cache.update(m.twoj(), i.twoj(), k, max_2j, kmax, SJ);
  sjr_cache.update(j.twoj(), n.twoj(), k, max_2j, kmax, SJ);

  // (n,i) change on every call so a cross-call cache for Q^u_{cnir} would never
  // hit; look up inline (only for pairs surviving selection rules).
  const auto core_size = core.size();
  const auto excited_size = excited.size();

  for (auto ir = 0ul; ir < excited_size; ++ir) {
    const auto &r = excited[ir];
    for (auto ic = 0ul; ic < core_size; ++ic) {
      const auto &c = core[ic];

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, i, r);
      // Rung range: L^l part exists at both parities of l (see L1)
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, j);
      const auto [w0, wI] = k_minmax_L(m, r, c, j);
      if (uI < u0 || wI < w0)
        continue;

      const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
      const auto inv_e_cjmr = 1.0 / (c.en() + ejm - r.en());
      const auto Qkey_cnir = qk.NormalOrder(c, n, i, r);
      const auto Qkey_mrcj = qk.NormalOrder(m, r, c, j);
      const auto Lkey_mrcj = Lk ? Lk->NormalOrder(m, r, c, j) : 0ul;

      // Cache (Q+L)^l_mrcj: depends only on l, used inside the u loop (see L1).
      assert(wI < int(sk_array_size));
      std::array<double, sk_array_size> QLl_mrcj{};
      for (auto l = w0; l <= wI; ++l) {
        const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
        QLl_mrcj[std::size_t(l)] =
          (has_Q ? qk.Q(l, Qkey_mrcj) : 0.0) + (Lk ? Lk->Q(l, Lkey_mrcj) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucnir = qk.Q(u, Qkey_cnir);
        // Zero when parity or triangle selection rules forbid this u.
        if (Q_ucnir == 0.0)
          continue;

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(i, u, c) == 0 || Coulomb::triangle(n, u, r) == 0)
          continue;

        for (auto l = w0; l <= wI; ++l) {

          const auto QL_lmrcj = QLl_mrcj[std::size_t(l)];
          if (QL_lmrcj == 0.0)
            continue;

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto s_ul = Angular::neg1pow(u + l);
          const auto sj_c = sjc_cache.get(c.twoj(), u, l);
          const auto sj_r = sjr_cache.get(r.twoj(), u, l);

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
                  const Angular::SixJTable &sjt, int max_k, bool print) {

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

  const auto kmax = (max_k < 0) ? qk.max_k() : std::min(max_k, qk.max_k());

  // Base selection rule: m,n in excited, i in i_orbs, b in core + angular SR.
  // Determines which L^k_mnib integrals we actually need.
  // nb: k runs over the full ladder range (k_minmax_L): triangle rules plus
  // the combined parity rule only. Unlike Q^k, L^k has components at BOTH
  // parities of k; the wrong-parity (Coulomb-forbidden) components enter the
  // exchange diagrams (via P/W) and feed back into the iterated rungs.
  const auto Lk_SR_one = [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                             const DiracSpinor &i,
                             const DiracSpinor &b) -> bool {
    // Require m and n to be excited
    if (!is_excited[m.nk_index()] || !is_excited[n.nk_index()])
      return false;
    // Require i to be in {i}, and b to be in core
    if (!is_i_orb[i.nk_index()] || !is_core[b.nk_index()])
      return false;
    const auto [k0, kI] = k_minmax_L(m, n, i, b);
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
GMatrix Sigma_ladder(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &projection,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     const Angular::SixJTable &sjt, bool include_L4, double r0,
                     double rmax, std::size_t stride, bool include_G) {

  // Ladder correction to the correlation potential, evaluated at energy e_v.
  // The exchange is folded into the Coulomb vertex via W = Q + P (mirrors
  // de_valence_w), so no ladder-P is needed. The bra index i runs over the
  // projection states of kappa_v (approximating completeness).
  //
  // Diagrams (a)+(b)  [particle-particle]:
  //   Sigma_L += sum_{i,amn,k} |W^k_{.amn}> w <i| ,  w = L^k_{mn,i,a}/([k][j_v]de)
  //   de = en_v + e_a - e_m - e_n
  // Diagrams (c)+(d)  [particle-hole, external line in the m-slot of L]:
  //   Sigma_L += sum_{i,nab,k} |W^k_{.nab}> w <i| ,  w = L^k_{i,n,a,b}/([k][j_v]de)
  //   de = en_v + e_n - e_a - e_b
  //
  // The ladder integrals are computed on-the-fly via Lkmnij() (using the
  // converged ladder table lk as the internal rung), evaluated at the fixed
  // external energy en_v via the e_i/e_m energy overrides (the energy only
  // enters the L denominators, never the integral lookups; the override
  // applies to whichever slot holds the external line: i-slot for a+b,
  // m-slot for c+d).
  // Exception: for the valence state itself (the projection state whose
  // orbital energy equals en_v), the stored (converged) table entries ARE the
  // required integrals - use them directly (makes the i = v term essentially
  // free, and consistent with de_valence). For other projection
  // states the stored entries cannot be used: they are not in the table (and
  // would be at the wrong energy).
  //
  // 'lk' is the internal-rung ladder table passed straight to Lkmnij: pass
  // nullptr for L(Q,Q) = L^(1) (matches an un-iterated table), or a converged
  // table (its fixed point) for the full ladder.
  //
  // If 'projection' is empty: ratio method instead (no projection; follows
  // Dzuba, PRA 78, 042502 (2008)). Each term of the regular second-order
  // Sigma (as in Goldstone::Sigma_both) is rescaled by the scalar ratio L/Q,
  // both taken straight from the stored tables (the L entries with i = v are
  // already at the valence energy - nothing computed on-the-fly; fast):
  //   (a+b): |Q^k_{.amn}><W^k_{.amn}| (L^k_{mnva}/Q^k_{mnva}) / ([k][j_v] de)
  //   (c+d): |Q^k_{.nab}><W^k_{.nab}| (L^k_{vnab}/Q^k_{vnab}) / ([k][j_v] de)
  // The diagonal then reproduces the Coulomb-parity part of the ladder
  // energy exactly; the radial shape is approximate (each term keeps its
  // second-order shape). Terms with Q^k = 0 cannot be rescaled: dropped. In
  // particular, the exchange-only k channels (Coulomb-forbidden parity,
  // where W = P and L^k != 0 - see k_minmax_L) are inherently absent, so
  // <v|Sigma_L|v> differs from de_valence_w(v) by those (small) terms.

  const auto ratio_method = projection.empty();

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

  // Projection basis: states with kappa == kappa_v
  std::vector<const DiracSpinor *> proj;
  for (const auto &x : projection) {
    if (x.kappa() == kappa_v)
      proj.push_back(&x);
  }
  if (proj.empty() && !ratio_method)
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

      // Diagrams (a)+(b): W^k_{.amn} L^k_{mn,i,a}
      // k runs over the full ladder range, both parities (see k_minmax_L):
      // at Coulomb-forbidden ("wrong"-parity) k the Q part of W vanishes,
      // but the P part and L^k do not - these exchange-only channels are
      // genuine (also present in de_valence_w). Exception: the ratio method
      // rescales second-order |Q><W| terms, so requires Q^k != 0 (Coulomb
      // parity); the exchange-only channels cannot be represented there.
      const auto [kmin_na, kmax_na] = Coulomb::k_minmax_tj(n.twoj(), a.twoj());
      for (int k = kmin_na; k <= kmax_na; ++k) {

        const auto f_kkjj = (2 * k + 1) * (tjv + 1);
        if (ratio_method && !Angular::Ck_kk_SR(k, n.kappa(), a.kappa()))
          continue;
        for (const auto &m : excited) {
          if (ratio_method) {
            if (!Angular::Ck_kk_SR(k, kappa_v, m.kappa()))
              continue;
          } else {
            const auto [k0, kI] = k_minmax_L(m, n, v, a);
            if (k < k0 || k > kI)
              continue;
          }
          const auto dele = en_v + a.en() - m.en() - n.en();

          if (ratio_method) {
            // Rescale the Sigma(2) term |Q><W| by L/Q (tables only):
            const auto Qkmnva = qk.Q(k, m, n, v, a);
            if (Qkmnva == 0.0)
              continue;
            const auto Lkmnva = lk ? lk->Q(k, m, n, v, a) : 0.0;
            if (Lkmnva == 0.0)
              continue;
            const auto Qkv = Coulomb::Qkv_bcd(k, kappa_v, a, m, n);
            const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, a, m, n);
            Sd_t.add(Qkv, Wkv, (Lkmnva / Qkmnva) / (f_kkjj * dele));
            continue;
          }

          // Wkv, dele don't depend on i: build once, sweep i below.
          const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, a, m, n);
          for (auto ii = 0ul; ii < proj.size(); ++ii) {
            const auto &Fi = *proj[ii];
            // Use stored table entry when we have it (see note above):
            const bool in_table =
              lk && Fi.en() == en_v && lk->contains(k, m, n, Fi, a);
            const auto Lkmnia = in_table ?
                                  lk->Q(k, m, n, Fi, a) :
                                  Lkmnij(k, m, n, Fi, a, qk, core, excited,
                                         include_L4, sjt, lk, en_v);
            if (Lkmnia == 0.0)
              continue;
            Sd_t.add(Wkv, Fi, Lkmnia / (f_kkjj * dele));
          }
        }
      }

      // Diagrams (c)+(d): W^k_{.nab} L^k_{i,n,a,b}
      for (const auto &b : core) {
        // k range: full ladder range, both parities; see (a)+(b) note.
        const auto [kmin_nb, kmax_nb] =
          Coulomb::k_minmax_tj(n.twoj(), b.twoj());
        for (int k = kmin_nb; k <= kmax_nb; ++k) {

          const auto f_kkjj = (2 * k + 1) * (tjv + 1);
          if (ratio_method) {
            if (!Angular::Ck_kk_SR(k, n.kappa(), b.kappa()) ||
                !Angular::Ck_kk_SR(k, kappa_v, a.kappa()))
              continue;
          } else {
            const auto [k0, kI] = k_minmax_L(v, n, a, b);
            if (k < k0 || k > kI)
              continue;
          }

          const auto dele = en_v + n.en() - a.en() - b.en();

          if (ratio_method) {
            // Rescale the Sigma(2) term |Q><W| by L/Q (tables only):
            const auto Qkvnab = qk.Q(k, v, n, a, b);
            if (Qkvnab == 0.0)
              continue;
            const auto Lkvnab = lk ? lk->Q(k, v, n, a, b) : 0.0;
            if (Lkvnab == 0.0)
              continue;
            const auto Qkv = Coulomb::Qkv_bcd(k, kappa_v, n, a, b);
            const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, n, a, b);
            Sd_t.add(Qkv, Wkv, (Lkvnab / Qkvnab) / (f_kkjj * dele));
            continue;
          }

          const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, n, a, b);
          for (auto ii = 0ul; ii < proj.size(); ++ii) {
            const auto &Fi = *proj[ii];
            // Use stored table entry when we have it (see note above):
            const bool in_table =
              lk && Fi.en() == en_v && lk->contains(k, Fi, n, a, b);
            const auto Lkinab = in_table ?
                                  lk->Q(k, Fi, n, a, b) :
                                  Lkmnij(k, Fi, n, a, b, qk, core, excited,
                                         include_L4, sjt, lk, {}, en_v);
            if (Lkinab == 0.0)
              continue;
            Sd_t.add(Wkv, Fi, Lkinab / (f_kkjj * dele));
          }
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
DiracSpinor Lkv_mnia(int k, const DiracSpinor &v, const DiracSpinor &m,
                     const DiracSpinor &n, const DiracSpinor &a,
                     const Coulomb::QkTable &qk, const Coulomb::YkTable &yk,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, bool include_L4,
                     const Angular::SixJTable &SJ,
                     const Coulomb::LkTable *const Lk) {

  // Vertex (ket) form of L^k_{mnia} over the external index i:
  // <x|Lkv_mnia> = Lkmnij(k, m, n, x, a; e_i = v.en()) for any x with
  // kappa_x = kappa_v. Mirrors L1/L2/L3(/L4) exactly, but the single Coulomb
  // line containing the external index is opened as a radial function
  // (yk.Qkv_bcd) instead of contracted with an orbital. The L-part of the
  // internal (Q+L) rung in L1/L3 contains the external line (attached to a
  // dressed rung): for that piece set i = v (scalar dv, one |v> add at the
  // end) - exact at x = v; identical to what the scalar table provides.

  const auto kappa_v = v.kappa();
  const auto tjv = v.twoj();
  const auto en_v = v.en();
  const double tkp1 = 2.0 * k + 1.0;

  DiracSpinor Lkv{0, kappa_v, v.grid_sptr()};
  double dv = 0.0; // coefficient of |v>: dressed-rung (internal-L) piece

  const int max_2j =
    std::max(DiracSpinor::max_tj(excited), DiracSpinor::max_tj(core));
  const int kmax = (std::max({m.twoj(), n.twoj(), tjv, a.twoj()}) + max_2j) / 2;

  //--------------------------------------------------------------------------
  // L1: sum_{rs,ul} Q^u_{mnrs} (Q+L)^l_{rs,x,a} / (en_v + e_a - e_r - e_s)
  // External: Q^l_{rs,x,a} = <x|Q^l(v)_{sra}>; L^l_{rs,x,a} -> i=v.
  {
    const auto s1 =
      Angular::neg1pow_2(2 + m.twoj() + n.twoj() + tjv + a.twoj());
    static thread_local SixJCache sjr_cache, sjs_cache;
    sjr_cache.update(m.twoj(), tjv, k, max_2j, kmax, SJ);
    sjs_cache.update(n.twoj(), a.twoj(), k, max_2j, kmax, SJ);

    for (const auto &r : excited) {
      for (const auto &s : excited) {
        const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
        // Opened Coulomb line exists only in the parity range [l0,lI]; the
        // L^l part (dv) exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] =
          Coulomb::k_minmax_Q(r.kappa(), s.kappa(), kappa_v, a.kappa());
        const auto [w0, wI] =
          k_minmax_L(r.kappa(), s.kappa(), kappa_v, a.kappa());
        if (uI < u0 || wI < w0)
          continue;
        const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
        const auto inv_e = 1.0 / (en_v + a.en() - r.en() - s.en());
        const auto Qkey_mnrs = qk.NormalOrder(m, n, r, s);

        // cl[l] = sum_u Q^u_{mnrs} sj_r sj_s
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> cl{};
        bool non_zero = false;
        for (int u = u0; u <= uI; u += 2) {
          const auto Q_umnrs = qk.Q(u, Qkey_mnrs);
          if (Q_umnrs == 0.0)
            continue;
          if (!Coulomb::triangle(u, r, m) || !Coulomb::triangle(u, s, n))
            continue;
          for (int l = w0; l <= wI; ++l) {
            if (Angular::triangle(k, l, u) == 0)
              continue;
            const auto sjsj =
              sjr_cache.get(r.twoj(), l, u) * sjs_cache.get(s.twoj(), l, u);
            if (sjsj == 0.0)
              continue;
            cl[std::size_t(l)] += sjsj * Q_umnrs;
            non_zero = true;
          }
        }
        if (!non_zero)
          continue;

        const auto factor = s1 * tkp1 * s_rs * inv_e;
        for (int l = w0; l <= wI; ++l) {
          const auto coef = factor * cl[std::size_t(l)];
          if (coef == 0.0)
            continue;
          // ket: bare Coulomb line - Coulomb parity range only
          if (l >= l0 && l <= lI && (l - l0) % 2 == 0) {
            Lkv += coef * yk.Qkv_bcd(l, kappa_v, s, r, a);
          }
          if (Lk) {
            dv += coef * Lk->Q(l, r, s, v, a);
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // L2: sum_{rc,ul} Q^u_{c,n,x,r} (Q+L)^l_{mrca} / (e_c + e_a - e_m - e_r)
  // External: Q^u_{c,n,x,r} = <x|Q^u(v)_{ncr}>; internal rung fully in tables.
  {
    const auto s2 =
      Angular::neg1pow_2(2 * k + m.twoj() + n.twoj() + tjv + a.twoj());
    static thread_local SixJCache sjc_cache, sjr_cache;
    sjc_cache.update(m.twoj(), tjv, k, max_2j, kmax, SJ);
    sjr_cache.update(a.twoj(), n.twoj(), k, max_2j, kmax, SJ);

    for (const auto &r : excited) {
      for (const auto &c : core) {
        const auto [u0, uI] =
          Coulomb::k_minmax_Q(c.kappa(), n.kappa(), kappa_v, r.kappa());
        // Rung range: L^l part exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, a);
        const auto [w0, wI] = k_minmax_L(m, r, c, a);
        if (uI < u0 || wI < w0)
          continue;
        const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
        const auto inv_e = 1.0 / (c.en() + a.en() - m.en() - r.en());
        const auto Qkey_mrca = qk.NormalOrder(m, r, c, a);
        const auto Lkey_mrca = Lk ? Lk->NormalOrder(m, r, c, a) : 0ul;

        // (Q+L)^l_{mrca}: fully internal - from tables
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> QLl{};
        for (int l = w0; l <= wI; ++l) {
          const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
          QLl[std::size_t(l)] = (has_Q ? qk.Q(l, Qkey_mrca) : 0.0) +
                                (Lk ? Lk->Q(l, Lkey_mrca) : 0.0);
        }

        const auto factor = s2 * tkp1 * s_rc * inv_e;
        for (int u = u0; u <= uI; u += 2) {
          if (Angular::triangle(tjv, 2 * u, c.twoj()) == 0 ||
              Angular::triangle(n.twoj(), 2 * u, r.twoj()) == 0)
            continue;
          double cu = 0.0;
          for (int l = w0; l <= wI; ++l) {
            if (QLl[std::size_t(l)] == 0.0)
              continue;
            if (Angular::triangle(k, l, u) == 0)
              continue;
            const auto s_ul = Angular::neg1pow(u + l);
            cu += s_ul * sjc_cache.get(c.twoj(), u, l) *
                  sjr_cache.get(r.twoj(), u, l) * QLl[std::size_t(l)];
          }
          if (cu == 0.0)
            continue;
          Lkv += (factor * cu) * yk.Qkv_bcd(u, kappa_v, n, c, r);
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // L3 = L2(k, n, m, a, x; e_j = en_v):
  // sum_{rc,ul} Q^u_{c,m,a,r} (Q+L)^l_{n,r,c,x} / (e_c + en_v - e_n - e_r)
  // External: Q^l_{n,r,c,x} = <x|Q^l(v)_{nrc}>;
  // L^l_{n,r,c,x} = L^l_{r,n,x,c} -> i=v.
  {
    const auto s3 =
      Angular::neg1pow_2(2 * k + n.twoj() + m.twoj() + a.twoj() + tjv);
    static thread_local SixJCache sjc_cache, sjr_cache;
    sjc_cache.update(n.twoj(), a.twoj(), k, max_2j, kmax, SJ);
    sjr_cache.update(tjv, m.twoj(), k, max_2j, kmax, SJ);

    for (const auto &r : excited) {
      for (const auto &c : core) {
        const auto [u0, uI] = Coulomb::k_minmax_Q(c, m, a, r);
        // Opened Coulomb line exists only in the parity range [l0,lI]; the
        // L^l part (dv) exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] =
          Coulomb::k_minmax_Q(n.kappa(), r.kappa(), c.kappa(), kappa_v);
        const auto [w0, wI] =
          k_minmax_L(n.kappa(), r.kappa(), c.kappa(), kappa_v);
        if (uI < u0 || wI < w0)
          continue;
        const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
        const auto inv_e = 1.0 / (c.en() + en_v - n.en() - r.en());
        const auto Qkey_cmar = qk.NormalOrder(c, m, a, r);

        // dl[l] = sum_u (-1)^{u+l} Q^u_{cmar} sj_c sj_r
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> dl{};
        bool non_zero = false;
        for (int u = u0; u <= uI; u += 2) {
          const auto Q_ucmar = qk.Q(u, Qkey_cmar);
          if (Q_ucmar == 0.0)
            continue;
          if (Angular::triangle(a.twoj(), 2 * u, c.twoj()) == 0 ||
              Angular::triangle(m.twoj(), 2 * u, r.twoj()) == 0)
            continue;
          for (int l = w0; l <= wI; ++l) {
            if (Angular::triangle(k, l, u) == 0)
              continue;
            const auto s_ul = Angular::neg1pow(u + l);
            const auto sjsj =
              sjc_cache.get(c.twoj(), u, l) * sjr_cache.get(r.twoj(), u, l);
            if (sjsj == 0.0)
              continue;
            dl[std::size_t(l)] += s_ul * sjsj * Q_ucmar;
            non_zero = true;
          }
        }
        if (!non_zero)
          continue;

        const auto factor = s3 * tkp1 * s_rc * inv_e;
        for (int l = w0; l <= wI; ++l) {
          const auto coef = factor * dl[std::size_t(l)];
          if (coef == 0.0)
            continue;
          // ket: bare Coulomb line - Coulomb parity range only
          if (l >= l0 && l <= lI && (l - l0) % 2 == 0) {
            Lkv += coef * yk.Qkv_bcd(l, kappa_v, n, r, c);
          }
          if (Lk) {
            dv += coef * Lk->Q(l, r, n, v, c);
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // L4: sum_{cd,ul} Q^u_{c,d,x,a} (Q+L)^l_{mncd} / (e_c + e_d - e_m - e_n)
  // External: Q^u_{c,d,x,a} = <x|Q^u(v)_{dca}>; internal rung fully in tables.
  if (include_L4) {
    const auto s4 =
      Angular::neg1pow_2(2 + m.twoj() + n.twoj() + tjv + a.twoj());
    static thread_local SixJCache sjc_cache, sjd_cache;
    sjc_cache.update(m.twoj(), tjv, k, max_2j, kmax, SJ);
    sjd_cache.update(n.twoj(), a.twoj(), k, max_2j, kmax, SJ);

    for (const auto &c : core) {
      for (const auto &d : core) {
        const auto [u0, uI] =
          Coulomb::k_minmax_Q(c.kappa(), d.kappa(), kappa_v, a.kappa());
        // Rung range: L^l part exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, c, d);
        const auto [w0, wI] = k_minmax_L(m, n, c, d);
        if (uI < u0 || wI < w0)
          continue;
        const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
        const auto inv_e = 1.0 / (c.en() + d.en() - m.en() - n.en());
        const auto Qkey_mncd = qk.NormalOrder(m, n, c, d);
        const auto Lkey_mncd = Lk ? Lk->NormalOrder(m, n, c, d) : 0ul;

        // (Q+L)^l_{mncd}: fully internal - from tables
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> QLl{};
        for (int l = w0; l <= wI; ++l) {
          const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
          QLl[std::size_t(l)] = (has_Q ? qk.Q(l, Qkey_mncd) : 0.0) +
                                (Lk ? Lk->Q(l, Lkey_mncd) : 0.0);
        }

        const auto factor = s4 * tkp1 * s_cd * inv_e;
        for (int u = u0; u <= uI; u += 2) {
          if (Angular::triangle(tjv, 2 * u, c.twoj()) == 0 ||
              Angular::triangle(a.twoj(), 2 * u, d.twoj()) == 0)
            continue;
          double cu = 0.0;
          for (int l = w0; l <= wI; ++l) {
            if (QLl[std::size_t(l)] == 0.0)
              continue;
            if (Angular::triangle(k, u, l) == 0)
              continue;
            cu += sjc_cache.get(c.twoj(), u, l) *
                  sjd_cache.get(d.twoj(), u, l) * QLl[std::size_t(l)];
          }
          if (cu == 0.0)
            continue;
          Lkv += (factor * cu) * yk.Qkv_bcd(u, kappa_v, d, c, a);
        }
      }
    }
  }

  if (dv != 0.0) {
    Lkv += dv * v;
  }
  return Lkv;
}

//==============================================================================
DiracSpinor Lkv_inab(int k, const DiracSpinor &v, const DiracSpinor &n,
                     const DiracSpinor &a, const DiracSpinor &b,
                     const Coulomb::QkTable &qk, const Coulomb::YkTable &yk,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited, bool include_L4,
                     const Angular::SixJTable &SJ,
                     const Coulomb::LkTable *const Lk) {

  // Vertex (ket) form of L^k_{inab} over the external index i (the m-slot):
  // <x|Lkv_inab> = Lkmnij(k, x, n, a, b; e_m = v.en()). Mirror of Lkv_mnia
  // for the particle-hole (c+d) diagrams: here L1/L3 have the external on a
  // bare Coulomb line (opened exactly), while the internal (Q+L) rung of
  // L2/L4 contains it (L-part: i = v, as in Lkv_mnia).

  const auto kappa_v = v.kappa();
  const auto tjv = v.twoj();
  const auto en_v = v.en();
  const double tkp1 = 2.0 * k + 1.0;

  DiracSpinor Lkv{0, kappa_v, v.grid_sptr()};
  double dv = 0.0; // coefficient of |v>: dressed-rung (internal-L) piece

  const int max_2j =
    std::max(DiracSpinor::max_tj(excited), DiracSpinor::max_tj(core));
  const int kmax = (std::max({tjv, n.twoj(), a.twoj(), b.twoj()}) + max_2j) / 2;

  //--------------------------------------------------------------------------
  // L1: sum_{rs,ul} Q^u_{x,n,r,s} (Q+L)^l_{rsab} / (e_a + e_b - e_r - e_s)
  // External: Q^u_{x,n,r,s} = <x|Q^u(v)_{nrs}>; internal rung fully in tables.
  {
    const auto s1 =
      Angular::neg1pow_2(2 + tjv + n.twoj() + a.twoj() + b.twoj());
    static thread_local SixJCache sjr_cache, sjs_cache;
    sjr_cache.update(tjv, a.twoj(), k, max_2j, kmax, SJ);
    sjs_cache.update(n.twoj(), b.twoj(), k, max_2j, kmax, SJ);

    for (const auto &r : excited) {
      for (const auto &s : excited) {
        const auto [u0, uI] =
          Coulomb::k_minmax_Q(kappa_v, n.kappa(), r.kappa(), s.kappa());
        // Rung range: L^l part exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, a, b);
        const auto [w0, wI] = k_minmax_L(r, s, a, b);
        if (uI < u0 || wI < w0)
          continue;
        const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
        const auto inv_e = 1.0 / (a.en() + b.en() - r.en() - s.en());
        const auto Qkey_rsab = qk.NormalOrder(r, s, a, b);
        const auto Lkey_rsab = Lk ? Lk->NormalOrder(r, s, a, b) : 0ul;

        // (Q+L)^l_{rsab}: fully internal - from tables
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> QLl{};
        for (int l = w0; l <= wI; ++l) {
          const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
          QLl[std::size_t(l)] = (has_Q ? qk.Q(l, Qkey_rsab) : 0.0) +
                                (Lk ? Lk->Q(l, Lkey_rsab) : 0.0);
        }

        const auto factor = s1 * tkp1 * s_rs * inv_e;
        for (int u = u0; u <= uI; u += 2) {
          if (Angular::triangle(2 * u, r.twoj(), tjv) == 0 ||
              Angular::triangle(2 * u, s.twoj(), n.twoj()) == 0)
            continue;
          double cu = 0.0;
          for (int l = w0; l <= wI; ++l) {
            if (QLl[std::size_t(l)] == 0.0)
              continue;
            if (Angular::triangle(k, l, u) == 0)
              continue;
            cu += sjr_cache.get(r.twoj(), l, u) *
                  sjs_cache.get(s.twoj(), l, u) * QLl[std::size_t(l)];
          }
          if (cu == 0.0)
            continue;
          Lkv += (factor * cu) * yk.Qkv_bcd(u, kappa_v, n, r, s);
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // L2: sum_{rc,ul} Q^u_{cnar} (Q+L)^l_{x,r,c,b} / (e_c + e_b - en_v - e_r)
  // External: Q^l_{x,r,c,b} = <x|Q^l(v)_{rcb}>; L^l_{x,r,c,b} -> i=v.
  {
    const auto s2 =
      Angular::neg1pow_2(2 * k + tjv + n.twoj() + a.twoj() + b.twoj());
    static thread_local SixJCache sjc_cache, sjr_cache;
    sjc_cache.update(tjv, a.twoj(), k, max_2j, kmax, SJ);
    sjr_cache.update(b.twoj(), n.twoj(), k, max_2j, kmax, SJ);

    for (const auto &r : excited) {
      for (const auto &c : core) {
        const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, a, r);
        // Opened Coulomb line exists only in the parity range [l0,lI]; the
        // L^l part (dv) exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] =
          Coulomb::k_minmax_Q(kappa_v, r.kappa(), c.kappa(), b.kappa());
        const auto [w0, wI] =
          k_minmax_L(kappa_v, r.kappa(), c.kappa(), b.kappa());
        if (uI < u0 || wI < w0)
          continue;
        const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
        const auto inv_e = 1.0 / (c.en() + b.en() - en_v - r.en());
        const auto Qkey_cnar = qk.NormalOrder(c, n, a, r);

        // dl[l] = sum_u (-1)^{u+l} Q^u_{cnar} sj_c sj_r
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> dl{};
        bool non_zero = false;
        for (int u = u0; u <= uI; u += 2) {
          const auto Q_ucnar = qk.Q(u, Qkey_cnar);
          if (Q_ucnar == 0.0)
            continue;
          if (Angular::triangle(a.twoj(), 2 * u, c.twoj()) == 0 ||
              Angular::triangle(n.twoj(), 2 * u, r.twoj()) == 0)
            continue;
          for (int l = w0; l <= wI; ++l) {
            if (Angular::triangle(k, l, u) == 0)
              continue;
            const auto s_ul = Angular::neg1pow(u + l);
            const auto sjsj =
              sjc_cache.get(c.twoj(), u, l) * sjr_cache.get(r.twoj(), u, l);
            if (sjsj == 0.0)
              continue;
            dl[std::size_t(l)] += s_ul * sjsj * Q_ucnar;
            non_zero = true;
          }
        }
        if (!non_zero)
          continue;

        const auto factor = s2 * tkp1 * s_rc * inv_e;
        for (int l = w0; l <= wI; ++l) {
          const auto coef = factor * dl[std::size_t(l)];
          if (coef == 0.0)
            continue;
          // ket: bare Coulomb line - Coulomb parity range only
          if (l >= l0 && l <= lI && (l - l0) % 2 == 0) {
            Lkv += coef * yk.Qkv_bcd(l, kappa_v, r, c, b);
          }
          if (Lk) {
            dv += coef * Lk->Q(l, v, r, c, b);
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // L3 = L2(k, n, x, b, a):
  // sum_{rc,ul} Q^u_{c,x,b,r} (Q+L)^l_{nrca} / (e_c + e_a - e_n - e_r)
  // External: Q^u_{c,x,b,r} = <x|Q^u(v)_{crb}>; internal rung fully in tables.
  {
    const auto s3 =
      Angular::neg1pow_2(2 * k + n.twoj() + tjv + b.twoj() + a.twoj());
    static thread_local SixJCache sjc_cache, sjr_cache;
    sjc_cache.update(n.twoj(), b.twoj(), k, max_2j, kmax, SJ);
    sjr_cache.update(a.twoj(), tjv, k, max_2j, kmax, SJ);

    for (const auto &r : excited) {
      for (const auto &c : core) {
        const auto [u0, uI] =
          Coulomb::k_minmax_Q(c.kappa(), kappa_v, b.kappa(), r.kappa());
        // Rung range: L^l part exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] = Coulomb::k_minmax_Q(n, r, c, a);
        const auto [w0, wI] = k_minmax_L(n, r, c, a);
        if (uI < u0 || wI < w0)
          continue;
        const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
        const auto inv_e = 1.0 / (c.en() + a.en() - n.en() - r.en());
        const auto Qkey_nrca = qk.NormalOrder(n, r, c, a);
        const auto Lkey_nrca = Lk ? Lk->NormalOrder(n, r, c, a) : 0ul;

        // (Q+L)^l_{nrca}: fully internal - from tables
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> QLl{};
        for (int l = w0; l <= wI; ++l) {
          const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
          QLl[std::size_t(l)] = (has_Q ? qk.Q(l, Qkey_nrca) : 0.0) +
                                (Lk ? Lk->Q(l, Lkey_nrca) : 0.0);
        }

        const auto factor = s3 * tkp1 * s_rc * inv_e;
        for (int u = u0; u <= uI; u += 2) {
          if (Angular::triangle(b.twoj(), 2 * u, c.twoj()) == 0 ||
              Angular::triangle(tjv, 2 * u, r.twoj()) == 0)
            continue;
          double cu = 0.0;
          for (int l = w0; l <= wI; ++l) {
            if (QLl[std::size_t(l)] == 0.0)
              continue;
            if (Angular::triangle(k, l, u) == 0)
              continue;
            const auto s_ul = Angular::neg1pow(u + l);
            cu += s_ul * sjc_cache.get(c.twoj(), u, l) *
                  sjr_cache.get(r.twoj(), u, l) * QLl[std::size_t(l)];
          }
          if (cu == 0.0)
            continue;
          Lkv += (factor * cu) * yk.Qkv_bcd(u, kappa_v, c, r, b);
        }
      }
    }
  }

  //--------------------------------------------------------------------------
  // L4: sum_{cd,ul} Q^u_{cdab} (Q+L)^l_{x,n,c,d} / (e_c + e_d - en_v - e_n)
  // External: Q^l_{x,n,c,d} = <x|Q^l(v)_{ncd}>; L^l_{x,n,c,d} -> i=v.
  if (include_L4) {
    const auto s4 =
      Angular::neg1pow_2(2 + tjv + n.twoj() + a.twoj() + b.twoj());
    static thread_local SixJCache sjc_cache, sjd_cache;
    sjc_cache.update(tjv, a.twoj(), k, max_2j, kmax, SJ);
    sjd_cache.update(n.twoj(), b.twoj(), k, max_2j, kmax, SJ);

    for (const auto &c : core) {
      for (const auto &d : core) {
        const auto [u0, uI] = Coulomb::k_minmax_Q(c, d, a, b);
        // Opened Coulomb line exists only in the parity range [l0,lI]; the
        // L^l part (dv) exists at both parities of l (see k_minmax_L)
        const auto [l0, lI] =
          Coulomb::k_minmax_Q(kappa_v, n.kappa(), c.kappa(), d.kappa());
        const auto [w0, wI] =
          k_minmax_L(kappa_v, n.kappa(), c.kappa(), d.kappa());
        if (uI < u0 || wI < w0)
          continue;
        const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
        const auto inv_e = 1.0 / (c.en() + d.en() - en_v - n.en());
        const auto Qkey_cdab = qk.NormalOrder(c, d, a, b);

        // dl[l] = sum_u Q^u_{cdab} sj_c sj_d
        assert(wI < int(sk_array_size));
        std::array<double, sk_array_size> dl{};
        bool non_zero = false;
        for (int u = u0; u <= uI; u += 2) {
          const auto Q_ucdab = qk.Q(u, Qkey_cdab);
          if (Q_ucdab == 0.0)
            continue;
          if (Angular::triangle(a.twoj(), 2 * u, c.twoj()) == 0 ||
              Angular::triangle(b.twoj(), 2 * u, d.twoj()) == 0)
            continue;
          for (int l = w0; l <= wI; ++l) {
            if (Angular::triangle(k, u, l) == 0)
              continue;
            const auto sjsj =
              sjc_cache.get(c.twoj(), u, l) * sjd_cache.get(d.twoj(), u, l);
            if (sjsj == 0.0)
              continue;
            dl[std::size_t(l)] += sjsj * Q_ucdab;
            non_zero = true;
          }
        }
        if (!non_zero)
          continue;

        const auto factor = s4 * tkp1 * s_cd * inv_e;
        for (int l = w0; l <= wI; ++l) {
          const auto coef = factor * dl[std::size_t(l)];
          if (coef == 0.0)
            continue;
          // ket: bare Coulomb line - Coulomb parity range only
          if (l >= l0 && l <= lI && (l - l0) % 2 == 0) {
            Lkv += coef * yk.Qkv_bcd(l, kappa_v, n, c, d);
          }
          if (Lk) {
            dv += coef * Lk->Q(l, v, n, c, d);
          }
        }
      }
    }
  }

  if (dv != 0.0) {
    Lkv += dv * v;
  }
  return Lkv;
}

//==============================================================================
GMatrix
Sigma_ladder_direct(const DiracSpinor &v, const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &excited,
                    const Coulomb::QkTable &qk, const Coulomb::YkTable &yk,
                    const Coulomb::LkTable *lk, const Angular::SixJTable &sjt,
                    bool include_L4, double r0, double rmax, std::size_t stride,
                    bool include_G) {

  // Ladder correction to the correlation potential via the direct (open
  // external line) method. Same diagram structure as Sigma_ladder, but the
  // bra side is the ladder vertex |L^k> (Lkv_mnia / Lkv_inab) rather than a
  // projection sum_i L^k_i <i|:
  //
  // Diagrams (a)+(b): Sigma_L += |W^k_{.amn}> <L^k_{mn.a}| / ([k][j_v] de)
  //   de = en_v + e_a - e_m - e_n
  // Diagrams (c)+(d): Sigma_L += |W^k_{.nab}> <L^k_{.nab}| / ([k][j_v] de)
  //   de = en_v + e_n - e_a - e_b

  const auto kappa_v = v.kappa();
  const auto en_v = v.en();
  const auto tjv = v.twoj();

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

      // Diagrams (a)+(b): |W^k_{.amn}> <L^k_{mn.a}|
      // k runs over the full ladder range, both parities (see k_minmax_L):
      // at Coulomb-forbidden k the Q part of W vanishes, but the P part and
      // L^k do not (as in de_valence_w and Sigma_ladder)
      const auto [kmin_na, kmax_na] = Coulomb::k_minmax_tj(n.twoj(), a.twoj());
      for (int k = kmin_na; k <= kmax_na; ++k) {
        const auto f_kkjj = (2 * k + 1) * (tjv + 1);
        for (const auto &m : excited) {
          const auto [k0, kI] = k_minmax_L(m, n, v, a);
          if (k < k0 || k > kI)
            continue;
          const auto dele = en_v + a.en() - m.en() - n.en();
          const auto Lv =
            Lkv_mnia(k, v, m, n, a, qk, yk, core, excited, include_L4, sjt, lk);
          const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, a, m, n);
          Sd_t.add(Wkv, Lv, 1.0 / (f_kkjj * dele));
        }
      }

      // Diagrams (c)+(d): |W^k_{.nab}> <L^k_{.nab}|
      for (const auto &b : core) {
        // k range: full ladder range, both parities; see (a)+(b) note
        const auto [kmin_nb, kmax_nb] =
          Coulomb::k_minmax_tj(n.twoj(), b.twoj());
        for (int k = kmin_nb; k <= kmax_nb; ++k) {
          const auto [k0, kI] = k_minmax_L(v, n, a, b);
          if (k < k0 || k > kI)
            continue;
          const auto f_kkjj = (2 * k + 1) * (tjv + 1);
          const auto dele = en_v + n.en() - a.en() - b.en();
          const auto Lv =
            Lkv_inab(k, v, n, a, b, qk, yk, core, excited, include_L4, sjt, lk);
          const auto Wkv = Coulomb::Wkv_bcd(k, kappa_v, n, a, b);
          Sd_t.add(Wkv, Lv, 1.0 / (f_kkjj * dele));
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
void update_Lk_mnib(Coulomb::LkTable *lk, const Coulomb::QkTable &qk,
                    const std::vector<DiracSpinor> &excited,
                    const std::vector<DiracSpinor> &core,
                    const std::vector<DiracSpinor> &update_i, bool include_L4,
                    const Angular::SixJTable &sjt,
                    const Coulomb::LkTable *const lk_prev, double a_damp,
                    bool print) {

  const auto basis = qip::merge(core, excited);

  // Lk integral with damping folded in
  const auto Lk_function = [&](int k, const DiracSpinor &m,
                               const DiracSpinor &n, const DiracSpinor &i,
                               const DiracSpinor &b) -> double {
    return Lkmnij(k, m, n, i, b, qk, core, excited, include_L4, sjt, lk_prev);
  };

  // Empty update_i => update everything.
  if (update_i.empty()) {
    lk->update(basis, Lk_function, a_damp, print);
    return;
  }

  // Otherwise restrict the re-iteration to entries whose i-slot is in
  // 'update_i' (b is always core; by the L^k_{mnib}=L^k_{nmbi} symmetry the i
  // index sits in slot 3 or 4). The filter is applied in update() before any
  // per-k lookup, so skipped entries cost (almost) nothing. Lets us converge
  // core (update_i=core) before valence (update_i=valence), the latter using
  // the now-frozen core rungs.
  std::unordered_set<DiracSpinor::Index> core_set, update_set;
  for (const auto &x : core) {
    core_set.insert(x.nk_index());
  }
  for (const auto &x : update_i) {
    update_set.insert(x.nk_index());
  }
  const auto in = [](const std::unordered_set<DiracSpinor::Index> &set,
                     const DiracSpinor &x) {
    return set.find(x.nk_index()) != set.cend();
  };
  const auto filter = [&](const DiracSpinor &, const DiracSpinor &,
                          const DiracSpinor &i, const DiracSpinor &b) {
    return (in(update_set, i) && in(core_set, b)) ||
           (in(update_set, b) && in(core_set, i));
  };

  lk->update(basis, Lk_function, a_damp, print, filter);
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
