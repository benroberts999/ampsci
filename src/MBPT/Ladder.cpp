#include "Ladder.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
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
  assert(std::find(excited.cbegin(), excited.cend(), m) != excited.cend());
  assert(std::find(excited.cbegin(), excited.cend(), n) != excited.cend());

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
  assert(std::find(core.cbegin(), core.cend(), m) == core.cend());
  assert(std::find(core.cbegin(), core.cend(), n) == core.cend());

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
  assert(std::find(excited.cbegin(), excited.cend(), m) != excited.cend());
  assert(std::find(excited.cbegin(), excited.cend(), n) != excited.cend());
  assert(std::find(core.cbegin(), core.cend(), m) == core.cend());

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  double l2 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnijk =
    Angular::neg1pow_2(2 * k + m.twoj() + n.twoj() + i.twoj() + j.twoj());

  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, i, r);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
      const auto inv_e_cjmr = 1.0 / (c.en() + j.en() - m.en() - r.en());
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
                  const Angular::SixJTable &sjt,
                  const Coulomb::LkTable *const lk_prev, bool print_progbar,
                  const std::vector<double> &fk) {

  const double a_damp = 0.35;
  const double b_damp = 1.0 - a_damp;

  if (!core.empty() && !excited.empty())
    assert(core.front().en() < excited.front().en());

  // Build combined basis: excited + core + any extra i_orbs (e.g. valence)
  std::vector<DiracSpinor> basis = excited;
  for (const auto &x : core)
    if (std::find(basis.begin(), basis.end(), x) == basis.end())
      basis.push_back(x);
  for (const auto &x : i_orbs)
    if (std::find(basis.begin(), basis.end(), x) == basis.end())
      basis.push_back(x);

  // Sets for O(1) membership checks inside lambdas
  using Index = DiracSpinor::Index;
  std::unordered_set<Index> excited_set, core_set, i_orbs_set;
  for (const auto &x : excited)
    excited_set.insert(x.nk_index());
  for (const auto &x : core)
    core_set.insert(x.nk_index());
  for (const auto &x : i_orbs)
    i_orbs_set.insert(x.nk_index());

  // Selection rule: m,n in excited, i in i_orbs, b in core + angular SR
  const auto Lk_SR = [&](int k, const DiracSpinor &m, const DiracSpinor &n,
                         const DiracSpinor &i, const DiracSpinor &b) -> bool {
    if (!excited_set.count(m.nk_index()) || !excited_set.count(n.nk_index()))
      return false;
    if (!i_orbs_set.count(i.nk_index()) || !core_set.count(b.nk_index()))
      return false;
    const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, i, b);
    return k >= k0 && k <= kI;
  };

  // Lk integral with damping folded in
  const auto Lk_function = [&](int k, const DiracSpinor &m,
                               const DiracSpinor &n, const DiracSpinor &i,
                               const DiracSpinor &b) -> double {
    auto L_new =
      Lkmnij(k, m, n, i, b, qk, core, excited, include_L4, sjt, lk_prev, fk);
    if (lk_prev) {
      const auto L_prev = lk_prev->Q(k, m, n, i, b);
      if (L_prev != 0.0)
        L_new = b_damp * L_new + a_damp * L_prev;
    }
    return L_new;
  };

  lk->fill(basis, Lk_function, Lk_SR, -1, print_progbar);
}

} // namespace MBPT
