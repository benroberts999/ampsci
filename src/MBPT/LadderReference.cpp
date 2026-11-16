#include "MBPT/LadderReference.hpp"
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "MBPT/LadderIteration.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <array>
#include <optional>
#include <vector>

namespace MBPT {

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

  const auto excited_size = excited.size();

  for (auto ir = 0ul; ir < excited_size; ++ir) {
    const auto &r = excited[ir];
    for (auto is = 0ul; is < excited_size; ++is) {
      const auto &s = excited[is];

      const auto [u0, uI] = Coulomb::k_minmax_Q(m, n, r, s);
      // Range of l of (Q+L)^l: the bare Coulomb part exists only in the
      // parity range [l0,lI]; the L^l part exists at both parities of l (see
      // k_minmax_L)
      const auto [l0, lI] = Coulomb::k_minmax_Q(r, s, i, j);
      const auto [w0, wI] = k_minmax_L(r, s, i, j);
      if (uI < u0 || wI < w0)
        continue;

      const auto s_rs = Angular::neg1pow_2(r.twoj() + s.twoj());
      const auto inv_e_ijrs =
        1.0 / (e_i.value_or(i.en()) + j.en() - r.en() - s.en());
      const auto Qkey_mnrs = qk.NormalOrder(m, n, r, s);
      const auto Qkey_rsij = qk.NormalOrder(r, s, i, j);
      const auto Lkey_rsij = Lk ? Lk->NormalOrder(r, s, i, j) : 0ul;

      // Cache (Q+L)^l_rsij: depends only on l, avoid N_u lookups
      std::array<double, ladder_max_multipoles> QLl_rsij{};
      for (auto l = w0; l <= wI; ++l) {
        const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
        QLl_rsij[std::size_t(l)] =
          (has_Q ? qk.Q(l, Qkey_rsij) : 0.0) + (Lk ? Lk->Q(l, Lkey_rsij) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_umnrs = qk.Q(u, Qkey_mnrs);
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

          const auto sj_r =
            SJ.get_2(m.twoj(), i.twoj(), 2 * k, 2 * l, 2 * u, r.twoj());
          const auto sj_s =
            SJ.get_2(n.twoj(), j.twoj(), 2 * k, 2 * l, 2 * u, s.twoj());

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

  const auto core_size = core.size();

  for (auto ic = 0ul; ic < core_size; ++ic) {
    const auto &c = core[ic];
    for (auto id = 0ul; id < core_size; ++id) {
      const auto &d = core[id];

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, d, i, j);
      // Range of l of (Q+L)^l: L^l part exists at both parities of l (see L1)
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, c, d);
      const auto [w0, wI] = k_minmax_L(m, n, c, d);
      if (uI < u0 || wI < w0)
        continue;

      const auto s_cd = Angular::neg1pow_2(c.twoj() + d.twoj());
      const auto inv_e_cdmn =
        1.0 / (c.en() + d.en() - e_m.value_or(m.en()) - n.en());
      const auto Qkey_cdij = qk.NormalOrder(c, d, i, j);
      const auto Qkey_mncd = qk.NormalOrder(m, n, c, d);
      const auto Lkey_mncd = Lk ? Lk->NormalOrder(m, n, c, d) : 0ul;

      // Cache (Q+L)^l_mncd: depends only on l, used inside the u loop
      std::array<double, ladder_max_multipoles> QLl_mncd{};
      for (auto l = w0; l <= wI; ++l) {
        const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
        QLl_mncd[std::size_t(l)] =
          (has_Q ? qk.Q(l, Qkey_mncd) : 0.0) + (Lk ? Lk->Q(l, Lkey_mncd) : 0.0);
      }

      for (auto u = u0; u <= uI; u += 2) {
        const auto Q_ucdij = qk.Q(u, Qkey_cdij);
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

          const auto sj_c =
            SJ.get_2(m.twoj(), i.twoj(), 2 * k, 2 * u, 2 * l, c.twoj());
          const auto sj_d =
            SJ.get_2(n.twoj(), j.twoj(), 2 * k, 2 * u, 2 * l, d.twoj());

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

  const auto core_size = core.size();
  const auto excited_size = excited.size();

  for (auto ir = 0ul; ir < excited_size; ++ir) {
    const auto &r = excited[ir];
    for (auto ic = 0ul; ic < core_size; ++ic) {
      const auto &c = core[ic];

      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, i, r);
      // Range of l of (Q+L)^l: L^l part exists at both parities of l (see L1)
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
      std::array<double, ladder_max_multipoles> QLl_mrcj{};
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
          const auto sj_c =
            SJ.get_2(m.twoj(), i.twoj(), 2 * k, 2 * u, 2 * l, c.twoj());
          const auto sj_r =
            SJ.get_2(j.twoj(), n.twoj(), 2 * k, 2 * u, 2 * l, r.twoj());

          l2 += (s_ul * s_rc * sj_c * sj_r) * Q_ucnir * QL_lmrcj * inv_e_cjmr;
        }
      }
    }
  }
  l2 *= s_mnijk * tkp1;
  return l2;
}

//------------------------------------------------------------------------------

double Pkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited,
              const Angular::SixJTable &SJ, const Coulomb::LkTable *const Lk) {
  // Each piece, and its partner [m <-> n, i <-> j] at the same k
  const auto P_mnij = P1(k, m, n, i, j, qk, core, excited, Lk) +
                      P2(k, m, n, i, j, qk, core, excited, SJ, Lk) +
                      P3(k, m, n, i, j, qk, core, excited, SJ, Lk);
  const auto P_nmji = P1(k, n, m, j, i, qk, core, excited, Lk) +
                      P2(k, n, m, j, i, qk, core, excited, SJ, Lk) +
                      P3(k, n, m, j, i, qk, core, excited, SJ, Lk);
  return P_mnij + P_nmji;
}

//------------------------------------------------------------------------------

double P1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited,
          const Coulomb::LkTable *const Lk) {

  double p1 = 0.0;
  for (const auto &r : excited) {
    for (const auto &c : core) {
      // Zero unless k is Coulomb-allowed for (c,n,r,j)
      const auto Q_kcnrj = qk.Q(k, c, n, r, j);
      if (Q_kcnrj == 0.0)
        continue;
      const auto QL_kmric =
        qk.Q(k, m, r, i, c) + (Lk ? Lk->Q(k, m, r, i, c) : 0.0);
      p1 += Q_kcnrj * QL_kmric / (i.en() + c.en() - m.en() - r.en());
    }
  }
  return Angular::neg1pow(k) * p1 / (2.0 * k + 1.0);
}

//------------------------------------------------------------------------------

double P2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk) {

  double p2 = 0.0;
  for (const auto &r : excited) {
    for (const auto &c : core) {
      // (Q+L)^k_mric: the L^k part exists at both parities of k (see L1)
      const auto QL_kmric =
        qk.Q(k, m, r, i, c) + (Lk ? Lk->Q(k, m, r, i, c) : 0.0);
      if (QL_kmric == 0.0)
        continue;

      // Exchange multipole of (c,n,r,j), over [k]: the sum over u of Q^u_cnjr
      const auto [u0, uI] = Coulomb::k_minmax_Q(c, n, j, r);
      double sum_u = 0.0;
      for (auto u = u0; u <= uI; u += 2) {
        sum_u +=
          SJ.get_2(c.twoj(), r.twoj(), 2 * k, n.twoj(), j.twoj(), 2 * u) *
          qk.Q(u, c, n, j, r);
      }
      p2 += sum_u * QL_kmric / (i.en() + c.en() - m.en() - r.en());
    }
  }
  return Angular::neg1pow(k) * p2;
}

//------------------------------------------------------------------------------

double P3(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::LkTable *const Lk) {

  double p3 = 0.0;
  for (const auto &r : excited) {
    for (const auto &c : core) {
      // Zero unless k is Coulomb-allowed for (c,n,r,j)
      const auto Q_kcnrj = qk.Q(k, c, n, r, j);
      if (Q_kcnrj == 0.0)
        continue;

      // Exchange multipole of (m,r,i,c), over [k]: the sum over l of
      // (Q+L)^l_mrci. Range of l: the bare Coulomb part exists only in the
      // parity range [l0,lI]; the L^l part exists at both parities of l (see
      // L1)
      const auto [l0, lI] = Coulomb::k_minmax_Q(m, r, c, i);
      const auto [w0, wI] = k_minmax_L(m, r, c, i);
      const auto Qkey_mrci = qk.NormalOrder(m, r, c, i);
      const auto Lkey_mrci = Lk ? Lk->NormalOrder(m, r, c, i) : 0ul;
      double sum_l = 0.0;
      for (auto l = w0; l <= wI; ++l) {
        const bool has_Q = l >= l0 && l <= lI && (l - l0) % 2 == 0;
        const auto QL_lmrci =
          (has_Q ? qk.Q(l, Qkey_mrci) : 0.0) + (Lk ? Lk->Q(l, Lkey_mrci) : 0.0);
        sum_l +=
          SJ.get_2(m.twoj(), i.twoj(), 2 * k, r.twoj(), c.twoj(), 2 * l) *
          QL_lmrci;
      }
      p3 += Q_kcnrj * sum_l / (i.en() + c.en() - m.en() - r.en());
    }
  }
  return Angular::neg1pow(k) * p3;
}

//------------------------------------------------------------------------------

double Vhp(int k, const DiracSpinor &m, const DiracSpinor &n,
           const DiracSpinor &i, const DiracSpinor &j,
           const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited,
           const Coulomb::LkTable *const Lk) {
  // Only a pair whose hole is a core orbital is a polarisation loop
  if (std::find(core.cbegin(), core.cend(), j) == core.cend())
    return 0.0;

  double h = 0.0;
  for (const auto &r : excited) {
    if (r.kappa() != n.kappa())
      continue;
    // y^(j)_nr = <n|y^0_jj|r>: the monopole of the interaction with the hole
    const auto y_nr = qk.R(0, n, j, r, j);
    const auto QL_kmrij =
      qk.Q(k, m, r, i, j) + (Lk ? Lk->Q(k, m, r, i, j) : 0.0);
    h += y_nr * QL_kmrij / (i.en() + j.en() - m.en() - r.en());
  }
  return -h;
}

//------------------------------------------------------------------------------

double Lk_vnab_F(int k, const DiracSpinor &v, const DiracSpinor &n,
                 const DiracSpinor &a, const DiracSpinor &b,
                 const Coulomb::QkTable &qk,
                 const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &excited,
                 const LadderDiagrams &diagrams,
                 const Coulomb::LkTable *const lk_F,
                 const Coulomb::NkTable *const L_vnab) {

  double L = 0.0;
  if (diagrams.P1) {
    // Polarisation loop on (n,b), acting on the closing amplitude
    double loop_nb = 0.0;
    for (const auto &r : excited) {
      for (const auto &c : core) {
        const auto Q_kcnrb = qk.Q(k, c, n, r, b);
        if (Q_kcnrb == 0.0)
          continue;
        const auto QL_kvrac =
          qk.Q(k, v, r, a, c) + (L_vnab ? L_vnab->Q(k, v, r, a, c) : 0.0);
        loop_nb += Q_kcnrb * QL_kvrac / (a.en() + c.en() - v.en() - r.en());
      }
    }
    L += Angular::neg1pow(k) * loop_nb / (2.0 * k + 1.0);
    // Polarisation loop on (v,a), acting on the core amplitude x_F_nrbc
    L += P1(k, n, v, b, a, qk, core, excited, lk_F);
  }
  if (diagrams.hole_particle) {
    // Hole-particle term on (n,b) only
    double h = 0.0;
    for (const auto &r : excited) {
      if (r.kappa() != n.kappa())
        continue;
      const auto y_nr = qk.R(0, n, b, r, b);
      const auto QL_kvrab =
        qk.Q(k, v, r, a, b) + (L_vnab ? L_vnab->Q(k, v, r, a, b) : 0.0);
      h += y_nr * QL_kvrab / (a.en() + b.en() - v.en() - r.en());
    }
    L -= h;
  }
  return L;
}

} // namespace MBPT
