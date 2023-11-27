#include "Goldstone.hpp"
#include "Coulomb/CoulombIntegrals.hpp"

namespace MBPT {

//==============================================================================
Goldstone::Goldstone(const std::vector<DiracSpinor> &basis,
                     const std::vector<DiracSpinor> &core, std::size_t i0,
                     std::size_t stride, std::size_t size, int n_min_core,
                     bool include_G, const HF::Breit *Br)
    : m_grid(basis.front().grid_sptr()),
      m_basis(DiracSpinor::split_by_core(basis, core, n_min_core)),
      m_Yeh(m_basis.second, m_basis.first),
      m_i0(i0),
      // m_imax(m_grid->getIndex(subgrid.rmax)),
      m_stride(stride),
      m_subgrid_points(size),
      m_n_min_core(n_min_core),
      m_include_G(include_G) {

  if (Br != nullptr) {
    m_Br = *Br; // copy is OK, small
  }

  bool verbose = true;
  if (verbose) {
    std::cout << "\nGoldstone diagrams:\n";
    std::cout << "Basis: " << DiracSpinor::state_config(basis) << "\n";
    fmt::print("Including n ≥ {} in internal hole states\n", n_min_core);
    if (m_include_G) {
      std::cout << "Including G parts of matrix\n";
    }
    if (m_Br) {
      std::cout << "Including Breit inside matrix\n";
    }
    printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
           "stride=%i]\n",
           m_grid->r(m_i0), m_grid->r(m_i0 + m_stride * (m_subgrid_points - 1)),
           int(m_subgrid_points), int(m_i0), int(m_stride));
  }
}

//==============================================================================
GMatrix Goldstone::Sigma_direct(int kappa_v, double en_v,
                                const std::vector<double> &fks,
                                const std::vector<double> &etaks,
                                int n_max_core) const {

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

  const auto &[t_holes, t_excited] = m_basis;

  // This is required by clang++14??
  // "reference to local binding 'holes' declared in enclosing function"
  const auto &holes = t_holes;
  const auto &excited = t_excited;

  const auto tjv = Angular::twoj_k(kappa_v);

  GMatrix Sd{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};

  if (holes.empty() || excited.empty())
    return Sd;

#pragma omp parallel for collapse(2)
  for (auto ia = 0ul; ia < holes.size(); ia++) {
    for (auto in = 0ul; in < excited.size(); in++) {
      const auto &a = holes[ia];
      const auto &n = excited[in];
      if (n_max_core > 0 && a.n() > n_max_core)
        continue;

      GMatrix Sd_an{m_i0, m_stride, m_subgrid_points, m_include_G, m_grid};

      const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax_Ck(n, a);
      for (int k = kmin_nb; k <= kmax_nb; k += 2) {

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

          if (m_Br) {
            const auto Bkv = m_Br->Bkv_bcd(k, kappa_v, a, m, n);
            Sd_an.add(Bkv, Qkv, factor);
          }
        }

        // Diagram (c) [direct]
        for (const auto &b : holes) {
          if (!Angular::Ck_kk_SR(k, kappa_v, b.kappa()))
            continue;
          const auto Qkv = m_Yeh.Qkv_bcd(k, kappa_v, n, b, a);
          const auto dele = en_v + n.en() - b.en() - a.en();
          const auto factor = etak * fk / (f_kkjj * dele);
          Sd_an.add(Qkv, Qkv, factor);

          if (m_Br) {
            const auto Bkv = m_Br->Bkv_bcd(k, kappa_v, n, b, a);
            Sd_an.add(Bkv, Qkv, factor);
          }
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

  const auto &[t_holes, t_excited] = m_basis;

  // This is required by clang++14??
  // "reference to local binding 'holes' declared in enclosing function"
  const auto &holes = t_holes;
  const auto &excited = t_excited;

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
      for (int k = kmin_nb; k <= kmax_nb; k += 2) {

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
          // screen both Coulomb lines??
          // const auto Pkv = m_Yeh.Pkv_bcd(k, kappa_v, a, m, n, fks);
          const auto Pkv = m_Yeh.Pkv_bcd(k, kappa_v, a, m, n);
          const auto dele = en_v + a.en() - m.en() - n.en();
          const auto factor = fk / (f_kkjj * dele);
          Sx_an.add(Qkv, Pkv, factor);

          if (m_Br) {
            const auto Bkv = m_Br->Bkv_bcd(k, kappa_v, a, m, n);
            Sx_an.add(Bkv, Pkv, factor);
          }
        }

        // Diagram (d) [exchange]
        for (const auto &b : holes) {
          if (!Angular::Ck_kk_SR(k, kappa_v, b.kappa()))
            continue;
          const auto Qkv = m_Yeh.Qkv_bcd(k, kappa_v, n, b, a);
          // screen both Coulomb lines??
          // const auto Pkv = m_Yeh.Pkv_bcd(k, kappa_v, n, b, a, fks);
          const auto Pkv = m_Yeh.Pkv_bcd(k, kappa_v, n, b, a);
          const auto dele = en_v + n.en() - b.en() - a.en();
          const auto factor = fk / (f_kkjj * dele);
          Sx_an.add(Qkv, Pkv, factor);

          if (m_Br) {
            const auto Bkv = m_Br->Bkv_bcd(k, kappa_v, n, b, a);
            Sx_an.add(Bkv, Pkv, factor);
          }
        } // b

      } // k

#pragma omp critical(sum_Sx)
      { Sx += Sx_an; }

    } // n
  }   // a

  // std::cin.get();

  return Sx.drj_in_place();
}

} // namespace MBPT