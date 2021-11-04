#pragma once
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <type_traits>
// namespace MBPT

//******************************************************************************
// Calculate energy shift (either ladder, or sigma2)
// nb: set lk to qk to det de(2)
template <typename Qintegrals, typename QorLintegrals>
double de_valence(const DiracSpinor &v, const Qintegrals &qk,
                  const QorLintegrals &lk, const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited) {

  static_assert(std::is_same_v<Qintegrals, Coulomb::YkTable> ||
                    std::is_base_of_v<Coulomb::CoulombTable, Qintegrals>,
                "Qintegrals must be YkTable or CoulombTable (QkTable/LkTable)");
  static_assert(
      std::is_same_v<QorLintegrals, Coulomb::YkTable> ||
          std::is_base_of_v<Coulomb::CoulombTable, QorLintegrals>,
      "QorLintegrals must be YkTable or CoulombTable (QkTable/LkTable)");

  const bool LisQ = [&]() {
    if constexpr (std::is_same_v<Qintegrals, QorLintegrals>)
      return &qk == &lk;
    else
      return false;
  }();

  auto de_v = 0.0;

  // nb: always put 'v' in correct* position. This ensures that actual state v
  // is used when evaluating Q/W integrals [as opposed to the excited basis
  // version of v]. Note: only matters when using YkTable for integrals
  // *correct == 1st or 3rd in Q, 1st or 4th in P

#pragma omp parallel for reduction(+ : de_v)
  for (auto in = 0ul; in < excited.size(); ++in) {
    const auto &n = excited[in];
    for (const auto &a : core) {

      // Diagrams (a) + (b)
      for (const auto &m : excited) {
        const auto inv_de = 1.0 / (v.en() + a.en() - m.en() - n.en());
        const auto [k0, kI] = Coulomb::k_minmax_Q(v, a, m, n);
        for (int k = k0; k <= kI; k += 2) {
          const auto Q_kvamn = qk.Q(k, v, a, m, n);

          const auto L_kvamn = LisQ ? Q_kvamn : lk.Q(k, m, n, v, a);
          const auto P_kvamn = lk.P(k, n, m, a, v); // \propto Q_mnva
          de_v += Q_kvamn * (L_kvamn + P_kvamn) * inv_de / (2 * k + 1);
        } // k
      }   // m

      // Diagrams (c) + (d)
      for (const auto &b : core) {
        const auto inv_de = 1.0 / (v.en() + n.en() - a.en() - b.en());
        const auto [k0, kI] = Coulomb::k_minmax_Q(v, n, a, b);
        for (int k = k0; k <= kI; k += 2) {
          const auto Q_kvnab = qk.Q(k, v, n, a, b);

          const auto L_kvnab = LisQ ? Q_kvnab : lk.Q(k, v, n, a, b);
          const auto P_kvnab = lk.P(k, v, n, a, b);
          de_v += Q_kvnab * (L_kvnab + P_kvnab) * inv_de / (2 * k + 1);
        } // k
      }   // b

      //
    } // a
  }   // n

  return de_v / v.twojp1();
}

//------------------------------------------------------------------------------
// Calculate energy shift (either ladder, or sigma2) for CORE
// lk may be regular Coulomb integrals [in which case this returns MBPT(2)
// correction], or Ladder diagrams [in which case this returns the ladder
// diagram correction]
template <typename Qintegrals, typename QorLintegrals>
double de_core(const Qintegrals &qk, const QorLintegrals &lk,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited) {

  static_assert(std::is_same_v<Qintegrals, Coulomb::YkTable> ||
                    std::is_base_of_v<Coulomb::CoulombTable, QorLintegrals>,
                "Qintegrals must be YkTable or CoulombTable (QkTable/LkTable)");

  const bool LisQ = [&]() {
    if constexpr (std::is_same_v<Qintegrals, QorLintegrals>)
      return &qk == &lk;
    else
      return false;
  }();

  auto de_c = 0.0;
#pragma omp parallel for reduction(+ : de_c)
  for (auto in = 0ul; in < excited.size(); ++in) {
    const auto &n = excited[in];
    for (const auto &m : excited) {
      for (const auto &a : core) {
        for (const auto &b : core) {
          const auto inv_de = 1.0 / (a.en() + b.en() - m.en() - n.en());
          const auto [k0, kI] = Coulomb::k_minmax_Q(a, b, m, n);
          for (int k = k0; k <= kI; k += 2) {
            const auto Q_kabmn = qk.Q(k, a, b, m, n);
            const auto L_kmnab = LisQ ? Q_kabmn : lk.Q(k, m, n, a, b);
            const auto P_kmnab = lk.P(k, m, n, a, b);
            de_c += 0.5 * Q_kabmn * (L_kmnab + P_kmnab) * inv_de / (2 * k + 1);
          }
        }
      }
    }
  }
  return de_c;
}
