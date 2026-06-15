#include "Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/omp.hpp"
#include <numeric>
#include <unordered_set>

namespace MBPT {

//==============================================================================
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited, bool include_L4,
              const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk,
              const std::vector<double> &fk) {

  const auto L123 = L1(k, m, n, i, j, qk, excited, SJ, Lk, fk) +
                    L2(k, m, n, i, j, qk, core, excited, SJ, Lk, fk) +
                    L3(k, m, n, i, j, qk, core, excited, SJ, Lk, fk);
  // Optionally include "4th" ladder diagram
  // nb: L4 not fully checked!
  if (include_L4)
    return L123 + L4(k, m, n, i, j, qk, core, SJ, Lk, fk);
  else
    return L123;
}

//------------------------------------------------------------------------------
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk,
          const std::vector<double> &fk) {

  // m (and n) must be excited states, as should 'excited'
  // Therefore, can test:
  // Ensured 'excited' is actually the excited orbitals
  // and that m and n are excited orbitals
  // (Still possible BOTH wrong at the same time..)
  // assert(std::find(excited.cbegin(), excited.cend(), m) != excited.cend());
  // assert(std::find(excited.cbegin(), excited.cend(), n) != excited.cend());

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

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

  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &s : excited) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_e_ijrs = 1.0 / (i.en() + j.en() - r.en() - s.en());
      const auto key_mnrs = qk.NormalOrder(m, n, r, s);
      const auto key_rsij = qk.NormalOrder(r, s, i, j);
      const auto lkey_rsij = Lk ? Lk->NormalOrder(r, s, i, j) : 0ul;

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = Fk(u) * qk.Q(u, key_mnrs);
        if (Q_umnrs == 0.0)
          continue; // never? Unless have k_cut

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(u, r, m) == 0 || Coulomb::triangle(u, s, n) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto sj_r = SJ.get(m, i, k, l, u, r);
          const auto sj_s = SJ.get(n, j, k, l, u, s);

          const auto Q_lrsij = Fk(l) * qk.Q(l, key_rsij);
          const auto L_lrsij = Lk ? Lk->Q(l, lkey_rsij) : 0.0;

          l1 +=
            (s_rs * sj_r * sj_s) * Q_umnrs * (Q_lrsij + L_lrsij) * inv_e_ijrs;
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
          const std::vector<double> &fk) {

  // m (and n) must be excited states, as should 'excited'
  // Therefore, can test:
  // Ensured 'excited' is actually the excited orbitals
  // and that m and n are excited orbitals
  // assert(std::find(core.cbegin(), core.cend(), m) == core.cend());
  // assert(std::find(core.cbegin(), core.cend(), n) == core.cend());

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  double l4 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnij1 =
    Angular::neg1pow_2(2 + m.twoj() + n.twoj() + i.twoj() + j.twoj());

  //  6j(r) Triads: {m,i,k}, {k,u,l}, {i,u,c}, {l,c,m}
  //  6j(s) Triads: {n,b,k}, {k,u,l}, {b,u,d}, {l,d,n}

  for (auto c_index = 0ul; c_index < core.size(); ++c_index) {
    const auto &c = core[c_index];
    for (const auto &d : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, d, i, j);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, c, d);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
      const auto inv_e_cdmn = 1.0 / (c.en() + d.en() - m.en() - n.en());
      const auto key_cdij = qk.NormalOrder(c, d, i, j);
      const auto key_mncd = qk.NormalOrder(m, n, c, d);
      const auto lkey_mncd = Lk ? Lk->NormalOrder(m, n, c, d) : 0ul;

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucdij = Fk(u) * qk.Q(u, key_cdij);
        if (Q_ucdij == 0.0)
          continue; // never? Unless have k_cut

        // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(i, u, c) == 0 || Coulomb::triangle(j, u, d) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, u, l) == 0)
            continue;

          const auto sj_c = SJ.get(m, i, k, u, l, c);
          const auto sj_d = SJ.get(n, j, k, u, l, d);

          const auto Q_lmncd = Fk(l) * qk.Q(l, key_mncd);
          const auto L_lmncd = Lk ? Lk->Q(l, lkey_mncd) : 0.0;

          l4 +=
            (s_cd * sj_c * sj_d) * Q_ucdij * (Q_lmncd + L_lmncd) * inv_e_cdmn;
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
          const Coulomb::LkTable *const Lk, const std::vector<double> &fk) {

  // m (and n) must be excited states, as should 'excited'
  // Therefore, can test:
  // Ensured 'excited' is actually the excited orbitals
  // and that m and n are excited orbitals
  // assert(std::find(excited.cbegin(), excited.cend(), m) != excited.cend());
  // assert(std::find(excited.cbegin(), excited.cend(), n) != excited.cend());
  // assert(std::find(core.cbegin(), core.cend(), m) == core.cend());

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  double l2 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnijk =
    Angular::neg1pow_2(2 * k + m.twoj() + n.twoj() + i.twoj() + j.twoj());
  const auto ejm = j.en() - m.en();

  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, i, r);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
      const auto inv_e_cjmr = 1.0 / (c.en() + ejm - r.en());
      const auto key_cnir = qk.NormalOrder(c, n, i, r);
      const auto key_mrcj = qk.NormalOrder(m, r, c, j);
      const auto lkey_mrcj = Lk ? Lk->NormalOrder(m, r, c, j) : 0ul;

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucnir = Fk(u) * qk.Q(u, key_cnir);
        if (Q_ucnir == 0.0)
          continue; // never? Unless have k_cut

        // // From 6J triads (this makes 1.5x speedup):
        if (Coulomb::triangle(i, u, c) == 0 || Coulomb::triangle(n, u, r) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto s_ul = Angular::neg1pow(u + l);

          const auto sj_c = SJ.get(m, i, k, u, l, c);
          const auto sj_r = SJ.get(j, n, k, u, l, r);

          const auto Q_lmrcj = Fk(l) * qk.Q(l, key_mrcj);
          const auto L_lmrcj = Lk ? Lk->Q(l, lkey_mrcj) : 0.0;

          l2 += (s_ul * s_rc * sj_c * sj_r) * Q_ucnir * (Q_lmrcj + L_lmrcj) *
                inv_e_cjmr;
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
                  const Angular::SixJTable &sjt, bool print,
                  const std::vector<double> &fk) {

  // const double a_damp = 0.35;
  // const double b_damp = 1.0 - a_damp;

  // Build combined basis: excited + core + any extra i_orbs (e.g. valence)
  const auto basis = qip::merge(core, excited);

  // NB: This is dumb - probably a much faster way!
  // Sets for O(1) membership checks for "selection rules"
  // (m,n) in excited, b in core, i in i_orbs
  std::unordered_set<DiracSpinor::Index> excited_set, core_set, i_orbs_set;
  for (const auto &n : excited) {
    excited_set.insert(n.nk_index());
  }
  for (const auto &a : core) {
    core_set.insert(a.nk_index());
  }
  for (const auto &i : i_orbs) {
    i_orbs_set.insert(i.nk_index());
  }

  const auto kmax = qk.max_k();

  // Selection rule: m,n in excited, i in i_orbs, b in core + angular SR
  const auto Lk_SR = [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                         const DiracSpinor &i, const DiracSpinor &b) -> bool {
    // XXX This should be checked! XXX
    // It determines which L^k_mnib integrals are calculated

    // Require m and n to be excited
    if (!excited_set.count(m.nk_index()) || !excited_set.count(n.nk_index()))
      return false;
    // Require i to be in {i}, and b to be in core
    if (!i_orbs_set.count(i.nk_index()) || !core_set.count(b.nk_index()))
      return false;
    const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, i, b); // correct??
    return k >= k0 && k <= kI;
  };

  // Lk integral
  const auto Lk_function = [&](int k, const DiracSpinor &m,
                               const DiracSpinor &n, const DiracSpinor &i,
                               const DiracSpinor &b) -> double {
    return Lkmnij(k, m, n, i, b, qk, core, excited, include_L4, sjt, nullptr,
                  fk);
  };

  lk->fill(basis, Lk_function, Lk_SR, kmax, print);
}

//==============================================================================
GMatrix Sigma_ladder(int kappa_v, double en_v,
                     const std::vector<DiracSpinor> &core,
                     const std::vector<DiracSpinor> &excited,
                     const std::vector<DiracSpinor> &basis,
                     const Coulomb::QkTable &qk, const Coulomb::LkTable *lk,
                     const Angular::SixJTable &sjt, bool include_L4,
                     const std::vector<double> &, const std::vector<double> &,
                     double r0, double rmax, std::size_t stride) {

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

  // get_k helper (screening / enhancement factors):
  // const auto get_k = [](int k, const std::vector<double> &f) {
  //   const auto sk = std::size_t(k);
  //   return sk < f.size() ? f[sk] : 1.0;
  // };

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
          // if (get_k(k, fk) == 0.0 || get_k(k, etak) == 0.0)
          //   continue;
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
            // if (get_k(k, fk) == 0.0 || get_k(k, etak) == 0.0)
            //   continue;
            if (!Angular::Ck_kk_SR(k, kappa_v, a.kappa()))
              continue;
            // Enforce the (n,b) Coulomb parity (only k with Q^k_{vnab} != 0);
            // see note in the (a)+(b) block.
            if (!Angular::Ck_kk_SR(k, n.kappa(), b.kappa()))
              continue;

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
                    bool print, const std::vector<double> &fk) {

  // Build combined basis: excited + core + any extra i_orbs (e.g. valence)
  // i_orbs not required (either core or subset of excited)
  // Used previously for "selection rules" - now SR is if we already calc'd it!
  (void)i_orbs;
  const auto basis = qip::merge(core, excited);

  // Lk integral with damping folded in
  const auto Lk_function = [&](int k, const DiracSpinor &m,
                               const DiracSpinor &n, const DiracSpinor &i,
                               const DiracSpinor &b) -> double {
    return Lkmnij(k, m, n, i, b, qk, core, excited, include_L4, sjt, lk_prev,
                  fk);
  };

  lk->update(basis, Lk_function, a_damp, print);
}

} // namespace MBPT
