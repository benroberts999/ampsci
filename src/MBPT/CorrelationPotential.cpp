#include "MBPT/CorrelationPotential.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

static int find_max_tj(const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &excited) {
  auto maxtj1 =
      std::max_element(core.cbegin(), core.cend(), DiracSpinor::comp_j)->twoj();
  auto maxtj2 =
      std::max_element(excited.cbegin(), excited.cend(), DiracSpinor::comp_j)
          ->twoj();
  return std::max(maxtj1, maxtj2);
}

namespace MBPT {

// Good idea to limit k? Doesn't seem to help much>
CorrelationPotential::CorrelationPotential(
    const std::vector<DiracSpinor> &core,
    const std::vector<DiracSpinor> &excited)
    : m_core(core), m_excited(excited),
      m_yec(m_excited.front().p_rgrid, &m_excited, &m_core), //
      m_maxk(find_max_tj(core, excited)),                    //
      m_6j(m_maxk, m_maxk)                                   //
{}

//******************************************************************************
DiracSpinor CorrelationPotential::Sigma2Fv(const DiracSpinor &v) const {

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

  const auto &Ck = m_yec.Ck();

  std::vector<DiracSpinor> SFv_a(m_core.size(),
                                 DiracSpinor(0, v.k, *(v.p_rgrid)));

  // Note: get_yk_ab() must only be called with k for which y^k_ab exists!
  // Therefore, must use the k_minmax() function provided

#pragma omp parallel for
  for (auto ia = 0ul; ia < m_core.size(); ia++) {
    const auto &a = m_core[ia];
    auto &sigma_a = SFv_a[ia];
    auto Qkv = DiracSpinor(0, v.k, *(v.p_rgrid)); // re-use to reduce alloc'ns
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = m_yec.k_minmax(n, a);
      const auto max_k = std::min(m_maxk, kmax_nb);
      for (int k = kmin_nb; k <= max_k; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * v.twojp1();
        const auto &yknb = m_yec(k, n, a);

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          if (Ck(k, v.k, m.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, a, m, n, k, yknb, Ck);
          const auto Qk = v * Qkv;
          const auto Pk =
              Coulomb::Pk_abcd(v, a, m, n, k, m_yec(m, a), Ck, m_6j);
          const auto dele = v.en + a.en - m.en - n.en;
          // nb: *= here "ruins" Q, but faster and we don't need Q any more
          sigma_a += (Qkv *= ((1.0 / dele / f_kkjj) * (Qk + Pk)));
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_core) {
          if (Ck(k, v.k, b.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, n, b, a, k, yknb, Ck);
          const auto Qk = v * Qkv;
          const auto Pk =
              Coulomb::Pk_abcd(v, n, b, a, k, m_yec(n, b), Ck, m_6j);
          const auto dele = v.en + n.en - b.en - a.en;
          // nb: *= here "ruins" Q, but faster and we don't need Q any more
          sigma_a += (Qkv *= ((1.0 / dele / f_kkjj) * (Qk + Pk)));
        } // b

      } // k
    }   // n
  }     // a

  return std::accumulate(SFv_a.begin(), SFv_a.end(),
                         DiracSpinor(0, v.k, *(v.p_rgrid)));
}

//******************************************************************************
double CorrelationPotential::Sigma2vw(const DiracSpinor &v,
                                      const DiracSpinor &w) const {
  if (v.k != w.k)
    return 0;

  const auto &Ck = m_yec.Ck();

  std::vector<double> delta_a(m_core.size());
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_core.size(); ia++) {
    const auto &a = m_core[ia];
    auto &del_a = delta_a[ia];
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = m_yec.k_minmax(n, a);
      const auto max_k = std::min(m_maxk, kmax_nb);
      for (int k = kmin_nb; k <= max_k; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * v.twojp1();
        const auto &yknb = m_yec.get_yk_ab(k, n, a);

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          const auto Qk = Coulomb::Qk_abcd(v, a, m, n, k, yknb, Ck);
          if (Qk == 0.0)
            continue;
          const auto &ybm = m_yec.get_y_ab(m, a);
          const auto Pk = Coulomb::Pk_abcd(v, a, m, n, k, ybm, Ck, m_6j);
          const auto dele = v.en + a.en - m.en - n.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qk + Pk)) * Qk;
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_core) {
          const auto Qk = Coulomb::Qk_abcd(v, n, b, a, k, yknb, Ck);
          if (Qk == 0.0)
            continue;
          const auto &yna = m_yec.get_y_ab(n, b);
          const auto Pk = Coulomb::Pk_abcd(v, n, b, a, k, yna, Ck, m_6j);
          const auto dele = v.en + n.en - b.en - a.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qk + Pk)) * Qk;
        } // b

      } // k
    }   // n
  }     // a

  return std::accumulate(delta_a.begin(), delta_a.end(), 0.0);
}

} // namespace MBPT
