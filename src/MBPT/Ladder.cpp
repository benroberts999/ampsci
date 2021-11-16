#include "Ladder.hpp"
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Widgets.hpp"
#include <numeric>

namespace MBPT {

//******************************************************************************
double Lkmnab(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &b,
              const Coulomb::CoulombTable &qk,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited,
              const Angular::SixJTable *const SJ,
              const Coulomb::CoulombTable *const Lk) {
  return L1(k, m, n, i, b, qk, excited, SJ, Lk) +
         L23(k, m, n, i, b, qk, core, excited, SJ, Lk);
}

//------------------------------------------------------------------------------
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &b,
          const Coulomb::CoulombTable &qk,
          const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable *const SJ,
          const Coulomb::CoulombTable *const Lk) {

  // L1^k_mnib
  //   = sum_{rs,ul} A^{kul}_mnrsib * Q^u_mnrs * (Q+L)^l_rsib / (e_ib - e_rs)
  // A^{kul}_mnrsib
  //   = (-1)^{m+n+r+s+i+b+1} * [k] * {m,i,k;l,u,r} * {n,b,k;l,u,s}

  double l1 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_mnab1 =
      Angular::neg1pow_2(2 + m.twoj() + n.twoj() + i.twoj() + b.twoj());

  // nb: including early triad checks cuts comp time down by ~5x!!
  if (!Coulomb::sixjTriads(m, i, k, {}, {}, {})) // {m,i,k;l,u,r}
    return 0.0;
  if (!Coulomb::sixjTriads(n, b, k, {}, {}, {})) // {n,b,k;l,u,s}
    return 0.0;

    // // nb: including early triad checks cuts comp time down by ~5x!!
    // // 6j(r) Triads: {m,i,k}, {k,l,u}, {i,l,r}, {u,r,m}
    // // 6j(s) Triads: {n,b,k}, {k,l,u}, {b,l,s}, {u,s,n}
    // if (Angular::triangle(m.twoj(), i.twoj(), 2 * k) == 0 ||
    //     Angular::triangle(n.twoj(), b.twoj(), 2 * k) == 0)
    //   return 0.0;

#pragma omp parallel for reduction(+ : l1)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &s : excited) {

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, b);
      if (lI < l0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_de = 1.0 / (i.en() + b.en() - r.en() - s.en());

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = qk.Q(u, m, n, r, s);
        if (Q_umnrs == 0.0)
          continue; // never?

        // From 6J triads:
        if (Angular::triangle(2 * u, r.twoj(), m.twoj()) == 0 ||
            Angular::triangle(2 * u, s.twoj(), n.twoj()) == 0)
          continue;

        for (auto l = l0; l <= lI; l += 2) {

          // 6j triad:
          if (Angular::triangle(k, l, u) == 0)
            continue;

          const auto sj_r = SJ->get(m, i, k, l, u, r);
          const auto sj_s = SJ->get(n, b, k, l, u, s);

          const auto Q_lrsab = qk.Q(l, r, s, i, b);
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
           const Angular::SixJTable *const SJ,
           const Coulomb::CoulombTable *const Lk) {
  // L2^k_mnib
  //   = sum_{rc} ((-1)^k / [k]) * P^k_cnrb * (P+Lambda)^k_mric / (e_ic - e_mr)
  // L3^k_mnib
  //   = sum_{rc} ((-1)^k / [k]) * P^k_cmri * (P+Lambda)^k_nrbc / (e_bc - e_nr)

  double l23 = 0.0;
  const double tkp1 = 2.0 * k + 1.0;
  const auto s_k = Angular::neg1pow(k);

  // // 6J triad (from Pmric)
  // // Valid?
  // if (Angular::triangle(m.twoj(), i.twoj(), 2 * k) == 0 ||
  //     Angular::triangle(n.twoj(), b.twoj(), 2 * k) == 0)
  //   return 0.0;

#pragma omp parallel for reduction(+ : l23)
  for (auto r_index = 0ul; r_index < excited.size(); ++r_index) {
    const auto &r = excited[r_index];
    for (const auto &c : core) {

      const auto inv_de_acmr = 1.0 / (i.en() + c.en() - m.en() - r.en());

      const auto P_kcnrb = qk.P(k, c, n, r, b);
      if (P_kcnrb != 0.0) {
        const auto P_kmrac = qk.P(k, m, r, i, c, SJ);
        const auto Lambda_kmrac = Lk ? Lk->P(k, m, r, i, c, SJ) : 0.0;
        l23 += P_kcnrb * (P_kmrac + Lambda_kmrac) * inv_de_acmr;
      }

      //-------------------------------------------------------

      const auto inv_de_bcnr = 1.0 / (b.en() + c.en() - n.en() - r.en());

      const auto P_kcmra = qk.P(k, c, m, r, i);
      if (P_kcmra != 0.0) {
        const auto P_knrbc = qk.P(k, n, r, b, c, SJ);
        const auto Lambda_knrbc = Lk ? Lk->P(k, n, r, b, c, SJ) : 0.0;
        l23 += P_kcmra * (P_knrbc + Lambda_knrbc) * inv_de_bcnr;
      }
    }
  }
  l23 *= (s_k / tkp1);
  return l23;
}

//******************************************************************************
void calculate_Lk_mnib(Coulomb::CoulombTable *lk,
                       const Coulomb::CoulombTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &i_orbs,
                       const Angular::SixJTable *const sjt,
                       const Coulomb::CoulombTable *const lk_prev,
                       bool print_progbar) {

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
      qip::progbar(++count, int(excited.size()));

    for (const auto &n : excited) {
      if (apply_symmetry && n < m)
        continue;
      for (const auto &i : i_orbs) {
        for (const auto &b : core) {
          // we only need L's when there are non-zero Q's
          const auto [k0, kI] = Coulomb::k_minmax_Q(m, n, i, b);
          for (int k = k0; k <= kI; k += 2) {

            auto L_kmnib =
                MBPT::Lkmnab(k, m, n, i, b, qk, core, excited, sjt, lk_prev);

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

            auto L_kmmbi =
                MBPT::Lkmnab(k, m, m, b, i, qk, core, excited, sjt, lk_prev);

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

  if (print_progbar)
    std::cout << "\n" << std::flush; // prog bar
}

} // namespace MBPT
