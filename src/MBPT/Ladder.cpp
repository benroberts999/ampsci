#include "Ladder.hpp"
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Widgets.hpp"
#include <numeric>

namespace MBPT {

//******************************************************************************
double
Lkmnab(int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
       const DiracSpinor &b, const Coulomb::CoulombTable &qk,
       const std::vector<DiracSpinor> &core,
       const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
       const Coulomb::CoulombTable *const Lk, const std::vector<double> &fk) {

  return L1(k, m, n, i, b, qk, excited, SJ, Lk, fk) +
         L2_v2(k, m, n, i, b, qk, core, excited, SJ, Lk, fk) +
         L3_v2(k, m, n, i, b, qk, core, excited, SJ, Lk, fk);

  // return L23(k, m, n, i, b, qk, core, excited, SJ, Lk, fk);
  //
  // return L1(k, m, n, i, b, qk, excited, SJ, Lk, fk);
  //
  // return L1(k, m, n, i, b, qk, excited, SJ, Lk, fk) +
  //        L23(k, m, n, i, b, qk, core, excited, SJ, Lk, fk);
}

//------------------------------------------------------------------------------
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &b,
          const Coulomb::CoulombTable &qk,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::CoulombTable *const Lk,
          const std::vector<double> &fk) {

  // L1^k_mnib
  //   = sum_{rs,ul} A^{kul}_mnrsib * Q^u_mnrs * (Q+L)^l_rsib / (e_ib - e_rs)
  // A^{kul}_mnrsib
  //   = (-1)^{m+n+r+s+i+b+1} * [k] * {m,i,k;l,u,r} * {n,b,k;l,u,s}

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  double l1 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnab1 =
      Angular::neg1pow_2(2 + m.twoj() + n.twoj() + i.twoj() + b.twoj());

  // nb: including early triad checks cuts comp time down by ~5x!!
  // // 6j(r) Triads: {m,i,k}, {k,l,u}, {i,l,r}, {u,r,m}
  // // 6j(s) Triads: {n,b,k}, {k,l,u}, {b,l,s}, {u,s,n}
  if (!Coulomb::sixjTriads(m, i, k, {}, {}, {})) // {m,i,k;l,u,r}
    return 0.0;
  if (!Coulomb::sixjTriads(n, b, k, {}, {}, {})) // {n,b,k;l,u,s}
    return 0.0;

#pragma omp parallel for reduction(+ : l1)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &s : excited) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, b);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_de = 1.0 / (i.en() + b.en() - r.en() - s.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = Fk(u) * qk.Q(u, m, n, r, s);
        if (Q_umnrs == 0.0)
          continue; // never? Unless have k_cut

        // From 6J triads:
        if (Coulomb::triangle(u, r, m) == 0 || Coulomb::triangle(u, s, n) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto sj_r = SJ.get(m, i, k, l, u, r);
          const auto sj_s = SJ.get(n, b, k, l, u, s);

          const auto Q_lrsab = Fk(l) * qk.Q(l, r, s, i, b);
          const auto L_lrsab = Lk ? Lk->Q(l, r, s, i, b) : 0.0;

          l1 += s_rs * sj_r * sj_s * Q_umnrs * (Q_lrsab + L_lrsab) * inv_de;
        }
      }
    }
  }
  l1 *= s_mnab1 * tkp1;
  return l1;
}

//------------------------------------------------------------------------------
double L23(int k, const DiracSpinor &m, const DiracSpinor &n,
           const DiracSpinor &i, const DiracSpinor &b,
           const Coulomb::CoulombTable &qk,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited,
           const Angular::SixJTable &SJ, const Coulomb::CoulombTable *const Lk,
           const std::vector<double> &fk) {
  // L2^k_mnib
  //   = sum_{rc} ((-1)^k / [k]) * P^k_cnrb * (P+Lambda)^k_mric / (e_ic - e_mr)
  // L3^k_mnib
  //   = sum_{rc} ((-1)^k / [k]) * P^k_cmri * (P+Lambda)^k_nrbc / (e_bc - e_nr)

  double l23 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_k = Angular::neg1pow(k);

  // 6J triads (common to L2 and L3)
  if (Coulomb::triangle(m, i, k) == 0 || Coulomb::triangle(n, b, k) == 0)
    return 0.0;

#pragma omp parallel for reduction(+ : l23)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      // 6J triads (common to L2 and L3)
      if (Coulomb::triangle(c, r, k) == 0)
        continue;

      const auto inv_de_acmr = 1.0 / (i.en() + c.en() - m.en() - r.en());

      const auto P_kcnrb = qk.P2(k, c, n, r, b, SJ, fk);
      if (P_kcnrb != 0.0) {
        const auto P_kmrac = qk.P2(k, m, r, i, c, SJ, fk);
        const auto Lambda_kmrac = Lk ? Lk->P(k, m, r, i, c, &SJ) : 0.0;
        l23 += P_kcnrb * (P_kmrac + Lambda_kmrac) * inv_de_acmr;
      }

      //-------------------------------------------------------

      const auto inv_de_bcnr = 1.0 / (b.en() + c.en() - n.en() - r.en());

      const auto P_kcmra = qk.P2(k, c, m, r, i, SJ, fk);
      if (P_kcmra != 0.0) {
        const auto P_knrbc = qk.P2(k, n, r, b, c, SJ, fk);
        const auto Lambda_knrbc = Lk ? Lk->P(k, n, r, b, c, &SJ) : 0.0;
        l23 += P_kcmra * (P_knrbc + Lambda_knrbc) * inv_de_bcnr;
      }
    }
  }
  l23 *= (s_k / tkp1);
  return l23;
}

//------------------------------------------------------------------------------
double
L2_v2(int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
      const DiracSpinor &b, const Coulomb::CoulombTable &qk,
      const std::vector<DiracSpinor> &core,
      const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
      const Coulomb::CoulombTable *const Lk, const std::vector<double> &fk) {

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  double l2 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnibk =
      Angular::neg1pow_2(m.twoj() + n.twoj() + i.twoj() + b.twoj() + 2 * k);

  // nb: including early triad checks cuts comp time down by ~5x!!
  // // 6j(r) Triads: {m,i,k}, {k,l,u}, {i,l,r}, {u,r,m}
  // // 6j(s) Triads: {n,b,k}, {k,l,u}, {b,l,s}, {u,s,n}
  // if (!Coulomb::sixjTriads(m, i, k, {}, {}, {})) // {m,i,k;l,u,r}
  //   return 0.0;
  // if (!Coulomb::sixjTriads(n, b, k, {}, {}, {})) // {n,b,k;l,u,s}
  //   return 0.0;

#pragma omp parallel for reduction(+ : l2)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, r, c, b);
      const auto [l0, lI] = Coulomb::k_minmax_Q(c, n, i, r);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rc = Angular::neg1pow_2(r.twoj() + c.twoj());
      const auto inv_de = 1.0 / (c.en() + b.en() - m.en() - r.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umrcb = Fk(u) * qk.Q(u, m, r, c, b);
        if (Q_umrcb == 0.0)
          continue; // never? Unless have k_cut

        // From 6J triads:
        // if (Coulomb::triangle(u, r, m) == 0 || Coulomb::triangle(u, s, n) ==
        // 0)
        //   continue;

        for (auto l = l0; l <= lI; l += 2) {

          // // 6j triad:
          // if (Angular::triangle(k, l, u) == 0)
          //   continue;

          const auto sj_c = SJ.get(m, i, k, l, u, c);
          const auto sj_r = SJ.get(b, n, k, l, u, r);

          const auto Q_lcnir = Fk(l) * qk.Q(l, c, n, i, r);
          const auto L_lcnir = Lk ? Lk->Q(l, c, n, i, r) : 0.0;

          const auto s_ul = Angular::neg1pow(l + u);

          l2 += s_rc * s_ul * sj_c * sj_r * Q_umrcb * (Q_lcnir + L_lcnir) *
                inv_de;
        }
      }
    }
  }
  l2 *= s_mnibk * tkp1;
  return l2;
}
//------------------------------------------------------------------------------
double
L3_v2(int k, const DiracSpinor &m, const DiracSpinor &n, const DiracSpinor &i,
      const DiracSpinor &b, const Coulomb::CoulombTable &qk,
      const std::vector<DiracSpinor> &core,
      const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
      const Coulomb::CoulombTable *const Lk, const std::vector<double> &fk) {

  // screening factors:
  auto Fk = [&fk](int l) {
    return l < (int)fk.size() ? fk[std::size_t(l)] : 1.0;
  };

  double l3 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnibk =
      Angular::neg1pow_2(m.twoj() + n.twoj() + i.twoj() + b.twoj() + 2 * k);

  // nb: including early triad checks cuts comp time down by ~5x!!
  // // 6j(r) Triads: {m,i,k}, {k,l,u}, {i,l,r}, {u,r,m}
  // // 6j(s) Triads: {n,b,k}, {k,l,u}, {b,l,s}, {u,s,n}
  // if (!Coulomb::sixjTriads(m, i, k, {}, {}, {})) // {m,i,k;l,u,r}
  //   return 0.0;
  // if (!Coulomb::sixjTriads(n, b, k, {}, {}, {})) // {n,b,k;l,u,s}
  //   return 0.0;

#pragma omp parallel for reduction(+ : l3)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, c, r, b);
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, n, i, c);
      if (uI < u0 || lI < l0)
        continue;

      const auto s_rc = Angular::neg1pow_2(c.twoj() + r.twoj());
      const auto inv_de = 1.0 / (c.en() + i.en() - n.en() - r.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umcrb = Fk(u) * qk.Q(u, m, c, r, b);
        if (Q_umcrb == 0.0)
          continue; // never? Unless have k_cut

        // From 6J triads:
        // if (Coulomb::triangle(u, r, m) == 0 || Coulomb::triangle(u, s, n) ==
        // 0)
        //   continue;

        for (auto l = l0; l <= lI; l += 2) {

          // // 6j triad:
          // if (Angular::triangle(k, l, u) == 0)
          //   continue;

          const auto sj_r = SJ.get(m, i, k, l, u, r);
          const auto sj_c = SJ.get(b, n, k, l, u, c);

          const auto Q_lrnic = Fk(l) * qk.Q(l, r, n, i, c);
          const auto L_lrnic = Lk ? Lk->Q(l, r, n, i, c) : 0.0;

          const auto s_ul = Angular::neg1pow(l + u);

          l3 += s_rc * s_ul * sj_c * sj_r * Q_umcrb * (Q_lrnic + L_lrnic) *
                inv_de;
        }
      }
    }
  }
  l3 *= s_mnibk * tkp1;
  return l3;
}

//******************************************************************************
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
  assert(core.front().en() < excited.front().en());

  // nb: we can apply the symmetry condition, and only calculate the unique
  // integrals, but ONLY if {i}={core}. When {i}={valence}, they are all unique!
  // Due to L_abcd = L_badc symmetry. ~2x speedup
  // NOTE: When {i}!=core, we need to calculate L_nnib AND L_nnbi - ie., in the
  // case where n=m, we need to add the {ib,bi case!}
  const auto apply_symmetry = &core == &i_orbs;

  int count = 0; // for prog bar
  for (const auto &m : excited) {
    if (print_progbar)
      qip::progbar50(count++, int(excited.size()));

    for (const auto &n : excited) {
      if (apply_symmetry && n < m)
        continue;
      for (const auto &i : i_orbs) {
        for (const auto &b : core) {
          // we only need L's when there are non-zero Q's
          const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, i, b);
          for (int k = k0; k <= kI; k += 2) {

            // a) only need L's with non-zero Q's [prob already case], b) may
            // have k_max cutoff in Qk
            if (qk.Q(k, m, n, i, b) == 0.0)
              continue;

            auto L_kmnib = MBPT::Lkmnab(k, m, n, i, b, qk, core, excited, sjt,
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

    // For case when {i}!=core, we have skipped L_mmbi terms; add them here
    if (!apply_symmetry) {
      for (const auto &b : core) {
        for (const auto &i : i_orbs) {
          const auto [k0, kI] = Coulomb::k_minmax_Q(m, m, b, i);
          for (int k = k0; k <= kI; k += 2) {

            auto L_kmmbi = MBPT::Lkmnab(k, m, m, b, i, qk, core, excited, sjt,
                                        lk_prev, fk);

            // If we have old value, 'damp' new value
            const auto L_prev = lk_prev ? lk_prev->Q(k, m, m, b, i) : 0.0;
            if (L_prev != 0.0) {
              L_kmmbi = b_damp * L_kmmbi + a_damp * L_prev;
            }

            // store new value in table
            lk->update(k, m, m, b, i, L_kmmbi);
          }
        }
      }
    }
  }
}

} // namespace MBPT
