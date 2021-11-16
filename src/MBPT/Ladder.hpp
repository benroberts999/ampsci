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

//! Ladder integral, L^k_mnab := L1 + L23
/*! @details
L1^k_mnib
  = sum_{rs,ul} A^{kul}_mnrsib * Q^u_mnrs * (Q+L)^l_rsib / (e_ib - e_rs)

A^{kul}_mnrsib
  = (-1)^{m+n+r+s+i+b+1} * [k] * {m,i,k;l,u,r} * {n,b,k;l,u,s}
*/
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &a, const DiracSpinor &b,
          const Coulomb::CoulombTable &qk,
          const std::vector<DiracSpinor> &excited,
          const Angular::SixJTable *const SJ = nullptr,
          const Coulomb::CoulombTable *const Lk = nullptr);

//! Ladder integral, L^k_mnab := L1 + L23; L23 = L2 + L3
/*! @details
L2^k_mnib
  = sum_{rc} ((-1)^k / [k]) * P^k_cnrb * (P+Lambda)^k_mric / (e_ic - e_mr)

L3^k_mnib
  = sum_{rc} ((-1)^k / [k]) * P^k_cmri * (P+Lambda)^k_nrbc / (e_bc - e_nr)
*/
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
