#include "Goldstone.hpp"

namespace MBPT {

//==============================================================================
Goldstone::Goldstone(const std::vector<DiracSpinor> &basis, double e_Fermi,
                     MBPT::rgrid_params subgrid, int n_min_core, bool include_G)
    : m_grid(basis.front().grid_sptr()),
      m_basis(DiracSpinor::split_by_energy(basis, e_Fermi, n_min_core)),
      m_Yeh(m_basis.second, m_basis.first),
      m_i0(m_grid->getIndex(subgrid.r0)),
      m_imax(m_grid->getIndex(subgrid.rmax)),
      m_stride(subgrid.stride),
      m_subgrid_points((m_imax - m_i0) / m_stride + 1),
      m_include_G(include_G) {}

//==============================================================================
GMatrix Goldstone::Sigma_direct(int kappa_v, double en_v,
                                const std::vector<double> &fks,
                                const std::vector<double> &etaks) const {

  // Four second-order diagrams:
  // Diagram (a):
  // |Q^k_amn><Q^k_amn| / de_amn / [k][j]
  // Diagram (b) (exchange):
  // |Q^k_amn><P^k_amn| / de_amn / [k][j]
  // Diagram (c):
  // |Q^k_nba><Q^k_nba| / de_nba / [k][j]
  // Diagram (d) (exchange):
  // |Q^k_nba><P^k_nba| / de_nba / [k][j]
  // where:
  // All indeces are summed over,
  // a & b are core states, n & m are virtual excited states,
  // k is multipolarity [Coloulmb expansion]
  // de_xyz = e_v + e_x - e_y - e_z

  const auto &[holes, excited] = m_basis;
  const auto tjv = Angular::twoj_k(kappa_v);

  GMatrix Sd{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};

  if (holes.empty() || excited.empty())
    return Sd;

#pragma omp parallel for collapse(2)
  for (auto ia = 0ul; ia < holes.size(); ia++) {
    for (auto in = 0ul; in < excited.size(); in++) {
      const auto &a = holes[ia];
      const auto &n = excited[in];

      GMatrix Sd_an{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};

      const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax_Ck(n, a);
      for (int k = kmin_nb; k <= kmax_nb; ++k) {

        const auto f_kkjj = (2 * k + 1) * (tjv + 1);

        // Effective screening parameter:
        const auto fk = get_k(k, fks);     // screening
        const auto etak = get_k(k, etaks); // hole-particle
        if (fk == 0.0 || etak == 0.0)
          continue;

        // Diagram (a) [direct]
        for (const auto &m : excited) {
          if (!Angular::Ck_kk_SR(k, kappa_v, m.kappa()))
            continue;
          const auto Qkv = m_Yeh.Qkv_bcd(k, kappa_v, a, m, n);
          const auto dele = en_v + a.en() - m.en() - n.en();
          const auto factor = etak * fk / (f_kkjj * dele);
          Sd_an.add(Qkv, Qkv, factor);
        }

        // Diagram (c) [direct]
        for (const auto &b : holes) {
          if (!Angular::Ck_kk_SR(k, kappa_v, b.kappa()))
            continue;
          const auto Qkv = m_Yeh.Qkv_bcd(k, kappa_v, n, b, a);
          const auto dele = en_v + n.en() - b.en() - a.en();
          const auto factor = etak * fk / (f_kkjj * dele);
          Sd_an.add(Qkv, Qkv, factor);
        } // b

      } // k

#pragma omp critical(sum_Sd)
      { Sd += Sd_an; }

    } // n
  }   // a

  return Sd.drj_in_place();
}

//==============================================================================
GMatrix Goldstone::Sigma_exchange(int kappa_v, double en_v,
                                  const std::vector<double> &fks) const {

  // Diagram (b) (exchange):
  // |Q^k_amn><P^k_amn| / de_amn / [k][j]
  // Diagram (d) (exchange):
  // |Q^k_nba><P^k_nba| / de_nba / [k][j]

  const auto &[holes, excited] = m_basis;
  const auto tjv = Angular::twoj_k(kappa_v);

  GMatrix Sx{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};

  if (holes.empty() || excited.empty())
    return Sx;

#pragma omp parallel for collapse(2)
  for (auto ia = 0ul; ia < holes.size(); ia++) {
    for (auto in = 0ul; in < excited.size(); in++) {
      const auto &a = holes[ia];
      const auto &n = excited[in];

      GMatrix Sx_an{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};

      const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax_Ck(n, a);
      for (int k = kmin_nb; k <= kmax_nb; ++k) {

        const auto f_kkjj = (2 * k + 1) * (tjv + 1);

        // Effective screening parameter:
        const auto fk = get_k(k, fks); // screening
        if (fk == 0.0)
          continue;

        // Diagram (b) [exchange]
        for (const auto &m : excited) {
          if (!Angular::Ck_kk_SR(k, kappa_v, m.kappa()))
            continue;
          const auto Qkv = m_Yeh.Qkv_bcd(k, kappa_v, a, m, n);
          const auto Pkv = m_Yeh.Pkv_bcd(k, kappa_v, a, m, n, fks);
          const auto dele = en_v + a.en() - m.en() - n.en();
          const auto factor = 1.0 / (f_kkjj * dele);
          Sx_an.add(Qkv, Pkv, fk * factor);
        }

        // Diagram (d) [exchange]
        for (const auto &b : holes) {
          if (!Angular::Ck_kk_SR(k, kappa_v, b.kappa()))
            continue;
          const auto Qkv = m_Yeh.Qkv_bcd(k, kappa_v, n, b, a);
          const auto Pkv = m_Yeh.Pkv_bcd(k, kappa_v, n, b, a, fks);
          const auto dele = en_v + n.en() - b.en() - a.en();
          const auto factor = 1.0 / (f_kkjj * dele);
          Sx_an.add(Qkv, Pkv, fk * factor);
        } // b

      } // k

#pragma omp critical(sum_Sd)
      { Sx += Sx_an; }

    } // n
  }   // a

  return Sx.drj_in_place();
}

} // namespace MBPT