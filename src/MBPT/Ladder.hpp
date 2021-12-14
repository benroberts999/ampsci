#pragma once
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

// struct Fk {
//   const std::vector<double> fk;
//   double scr(int k) const {
//     return k < (int)fk.size() ? fk[std::size_t(k)] : 1.0;
//   }
// };
// // XXX Of course, this is temporary!
// const Fk global_fk{{0.654, 0.514, 0.792, 0.857, 0.930, 0.967, 0.986, 0.994}};
// const Fk
// global_eta{{1.258, 1.511, 1.281, 1.175, 1.114, 1.073, 1.062, 1.066}};

// const Fk global_fk{{1.0}};
// const Fk global_eta{{1.0}};

// // Weighted average:
// fk = 0.714, 0.617, 0.832, 0.885, 0.941, 0.970, 0.987, 0.994
// fk = 0.654, 0.514, 0.792, 0.857, 0.930, 0.967, 0.986, 0.994 // w/ hp
// eta= 1.258, 1.511, 1.281, 1.175, 1.114, 1.073, 1.062, 1.066

//! Calculates ladder integral, L^k_mnab
/*! @details Lk pointer is pointer to previous iteration of Lk. fk is optional
 * vector of screening factors
 */
double Lkmnab(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &a, const DiracSpinor &b,
              const Coulomb::CoulombTable &qk,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited,
              const Angular::SixJTable &SJ,
              const Coulomb::CoulombTable *const Lk = nullptr,
              const std::vector<double> &fk = {});

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
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::CoulombTable *const Lk = nullptr,
          const std::vector<double> &fk = {});

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
           const Angular::SixJTable &SJ,
           const Coulomb::CoulombTable *const Lk = nullptr,
           const std::vector<double> &fk = {});

//! Fills Lk matrix
void fill_Lk_mnib(Coulomb::CoulombTable *lk, const Coulomb::CoulombTable &qk,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &i_orbs,
                  const Angular::SixJTable &sjt,
                  const Coulomb::CoulombTable *const lk_prev = nullptr,
                  bool print_progbar = true,
                  const std::vector<double> &fk = {});

//! Calculate energy shift (either ladder, or sigma2) for valence
/*! @details
   lk may be regular Coulomb integrals [in which case this returns
   MBPT(2) correction], or Ladder diagrams [in which case this returns the
   ladder diagram correction]
*/
template <typename Qintegrals, typename QorLintegrals>
double de_valence(const DiracSpinor &v, const Qintegrals &qk,
                  const QorLintegrals &lk, const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const std::vector<double> &fk = {},
                  const std::vector<double> &etak = {});

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
