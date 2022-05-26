#include "Ladder.hpp"
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Widgets.hpp"
#include <numeric>

namespace MBPT {

//==============================================================================
double
Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
       const DiracSpinor &j, const Coulomb::CoulombTable &qk,
       const std::vector<DiracSpinor> &core,
       const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
       const Coulomb::CoulombTable *const Lk, const std::vector<double> &fk) {

  return L1(k, m, n, i, j, qk, excited, SJ, Lk, fk) +
         L2(k, m, n, i, j, qk, core, excited, SJ, Lk, fk) +
         // L2(k, n, m, j, i, qk, core, excited, SJ, Lk, fk) +
         L3(k, m, n, i, j, qk, core, excited, SJ, Lk, fk) +
         L4(k, m, n, i, j, qk, core, SJ, Lk, fk);
}

//------------------------------------------------------------------------------
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::CoulombTable &qk,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::CoulombTable *const Lk,
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

#pragma omp parallel for reduction(+ : l1)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &s : excited) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_e_ijrs = 1.0 / (i.en() + j.en() - r.en() - s.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = Fk(u) * qk.Q(u, m, n, r, s);
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

          const auto Q_lrsij = Fk(l) * qk.Q(l, r, s, i, j);
          const auto L_lrsij = Lk ? Lk->Q(l, r, s, i, j) : 0.0;

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
          const Coulomb::CoulombTable &qk, const std::vector<DiracSpinor> &core,
          const Angular::SixJTable &SJ, const Coulomb::CoulombTable *const Lk,
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

#pragma omp parallel for reduction(+ : l4)
  for (auto c_index = 0ul; c_index < core.size(); ++c_index) {
    const auto &c = core[c_index];
    for (const auto &d : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, d, i, j);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, c, d);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
      const auto inv_e_cdmn = 1.0 / (c.en() + d.en() - m.en() - n.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucdij = Fk(u) * qk.Q(u, c, d, i, j);
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

          const auto Q_lmncd = Fk(l) * qk.Q(l, m, n, c, d);
          const auto L_lmncd = Lk ? Lk->Q(l, m, n, c, d) : 0.0;

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
          const Coulomb::CoulombTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::CoulombTable *const Lk,
          const std::vector<double> &fk) {

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

#pragma omp parallel for reduction(+ : l2)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, i, r);
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, j);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
      const auto inv_e_cjmr = 1.0 / (c.en() + j.en() - m.en() - r.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucnir = Fk(u) * qk.Q(u, c, n, i, r);
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

          const auto Q_lmrcj = Fk(l) * qk.Q(l, m, r, c, j);
          const auto L_lmrcj = Lk ? Lk->Q(l, m, r, c, j) : 0.0;

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
void fill_Lk_mnib(Coulomb::CoulombTable *lk, const Coulomb::CoulombTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs,
                  const Angular::SixJTable &sjt,
                  const Coulomb::CoulombTable *const lk_prev,
                  bool print_progbar, const std::vector<double> &fk) {

  const double a_damp = 0.35; // 0 means no damping
  const double b_damp = 1.0 - a_damp;

  // Check core/excited orbitals around right way:
  if (!core.empty() && !excited.empty())
    assert(core.front().en() < excited.front().en());

  // nb: we can apply the symmetry condition, and only calculate the unique
  // integrals by only calculating those when {mnib} is already NormalOrdered

  int count = 0; // for prog bar
  for (const auto &m : excited) {
    if (print_progbar)
      qip::progbar50(count++, int(excited.size()));

    for (const auto &n : excited) {
      for (const auto &i : i_orbs) {
        for (const auto &b : core) {
          if (!lk->is_NormalOrdered(m, n, i, b))
            continue;
          // we only need L's when there are non-zero Q's
          const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, i, b);
          for (int k = k0; k <= kI; k += 2) {

            // a) only need L's with non-zero Q's [prob already case], b) may
            // have k_max cutoff in Qk
            if (qk.Q(k, m, n, i, b) == 0.0)
              continue;

            auto L_kmnib = MBPT::Lkmnij(k, m, n, i, b, qk, core, excited, sjt,
                                        lk_prev, fk);

            // If we have old value, 'damp' new value
            const auto L_prev = lk_prev ? lk_prev->Q(k, m, n, i, b) : 0.0;
            if (L_prev != 0.0) {
              L_kmnib = b_damp * L_kmnib + a_damp * L_prev;
            }

            // store new value in table
            lk->update(k, m, n, i, b, L_kmnib);
          }
        }
      }
    }
  }
}

} // namespace MBPT
