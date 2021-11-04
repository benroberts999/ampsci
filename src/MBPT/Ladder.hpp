#pragma once
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

//! Calculates ladder integral, L^k_mnab
double Lkmnab(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &a, const DiracSpinor &b,
              const Coulomb::CoulombTable &qk,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited,
              const Angular::SixJTable *const SJ = nullptr,
              const Coulomb::CoulombTable *const Lk = nullptr);

//! Ladder integral, L^k_mnab := L1 + L2
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b,
          const Coulomb::CoulombTable &qk,
          const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable *const SJ = nullptr,
          const Coulomb::CoulombTable *const Lk = nullptr);

//! Ladder integral, L^k_mnab := L1 + L2
double L23(int k, const DiracSpinor &m, const DiracSpinor &n,
           const DiracSpinor &a, const DiracSpinor &b,
           const Coulomb::CoulombTable &qk,
           const std::vector<DiracSpinor> &core,
           const std::vector<DiracSpinor> &excited,
           const Angular::SixJTable *const SJ = nullptr,
           const Coulomb::CoulombTable *const Lk = nullptr);

//! Fills Lk matrix
void calculate_Lk_mnib(Coulomb::CoulombTable *lk,
                       const Coulomb::CoulombTable &qk,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &i_orbs,
                       const Angular::SixJTable *const sjt = nullptr,
                       const Coulomb::CoulombTable *const lk_prev = nullptr,
                       bool print_progbar = true);

//! Calculate energy shift (either ladder, or sigma2) for valence
/*! @details
   lk may be regular Coulomb integrals [in which case this returns
   MBPT(2) correction], or Ladder diagrams [in which case this returns the
   ladder diagram correction]
*/
template <typename Qintegrals, typename QorLintegrals>
double de_valence(const DiracSpinor &v, const Qintegrals &qk,
                  const QorLintegrals &lk, const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited);

//! Calculate energy shift (either ladder, or sigma2) for CORE
/*! @details
lk may be regular Coulomb integrals [in which case this returns MBPT(2)
correction], or Ladder diagrams [in which case this returns the ladder
diagram correction]
*/
template <typename Qintegrals, typename QorLintegrals>
double de_core(const Qintegrals &qk, const QorLintegrals &lk,
               const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited);

// template implementations:
#include "Ladder.ipp"

} // namespace MBPT
