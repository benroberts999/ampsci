#pragma once
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "MBPT/RDMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"

namespace MBPT {

// // Weighted average:
// fk = 0.714, 0.617, 0.832, 0.885, 0.941, 0.970, 0.987, 0.994
// fk = 0.654, 0.514, 0.792, 0.857, 0.930, 0.967, 0.986, 0.994 // w/ hp
// eta= 1.258, 1.511, 1.281, 1.175, 1.114, 1.073, 1.062, 1.066

//! Calculates ladder integral, L^k_mnab
/*! @details Lk pointer is pointer to previous iteration of Lk. fk is optional
 * vector of screening factors
 */
double Lkmnij(int k, const DiracSpinor &m, const DiracSpinor &n,
              const DiracSpinor &i, const DiracSpinor &j,
              const Coulomb::CoulombTable &qk,
              const std::vector<DiracSpinor> &core,
              const std::vector<DiracSpinor> &excited,
              const Angular::SixJTable &SJ,
              const Coulomb::CoulombTable *const Lk = nullptr,
              const std::vector<double> &fk = {});

//! Ladder integral, L^k_mnij := L1_mnij + L2_mnij + L2_nmji
/*! @details
L1^k_mnij
  = sum_{rs,ul} A^{kul}_mnrsij * Q^u_mnrs * (Q+L)^l_rsij / (e_ij - e_rs)

A^{kul}_mnrsij
  = (-1)^{m+n+r+s+i+j+1} * [k] * {m,i,k;l,u,r} * {n,j,k;l,u,s}
*/
double L1(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::CoulombTable &qk,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::CoulombTable *const Lk = nullptr,
          const std::vector<double> &fk = {});

//! Ladder integral, L^k_mnab := L1_mnij + L2_mnij + L3_nmji,  L3_mnij = L2_nmji
/*! @details
L2^k_mnij
  = sum_{rc,ul} (-1)^{k+u+l+1} A^{klu}_mjcrin
    * Q^u_cnir * (Q+L)^l_mrcj / (e_cj - e_mr)
 */
double L2(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::CoulombTable &qk, const std::vector<DiracSpinor> &core,
          const std::vector<DiracSpinor> &excited, const Angular::SixJTable &SJ,
          const Coulomb::CoulombTable *const Lk = nullptr,
          const std::vector<double> &fk = {});

inline double L3(int k, const DiracSpinor &m, const DiracSpinor &n,
                 const DiracSpinor &i, const DiracSpinor &j,
                 const Coulomb::CoulombTable &qk,
                 const std::vector<DiracSpinor> &core,
                 const std::vector<DiracSpinor> &excited,
                 const Angular::SixJTable &SJ,
                 const Coulomb::CoulombTable *const Lk = nullptr,
                 const std::vector<double> &fk = {}) {
  return L2(k, n, m, j, i, qk, core, excited, SJ, Lk, fk);
}

double L4(int k, const DiracSpinor &m, const DiracSpinor &n,
          const DiracSpinor &i, const DiracSpinor &j,
          const Coulomb::CoulombTable &qk, const std::vector<DiracSpinor> &core,
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
