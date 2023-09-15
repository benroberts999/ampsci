#pragma once
#include "Angular/Angular.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

namespace MBPT {

//! Type of energy demoninators: Rayleigh-Schrodinger (RS),
//! Brillouin-Wigner (BW). (not exact)
enum class Denominators { RS, BW };

//! Splits the basis into the core (holes) and excited states.
/*! @details 
 States with energy below E_Fermi are considered core/holes; only core states
 with n>=min_n_core, and excited states with n<=max_n_excited are kept.
*/
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>>
split_basis(const std::vector<DiracSpinor> &basis, double E_Fermi,
            int min_n_core = 1, int max_n_excited = 999);

//! Reduced two-body Sigma (2nd order correlation) operator. Sum of 6 diagrams
/*! @details 
  Note: these have fewer symmetries to Q^k; S_vwxy = S_xyvw
*/
double Sk_vwxy(int k, const DiracSpinor &v, const DiracSpinor &w,
               const DiracSpinor &x, const DiracSpinor &y,
               const Coulomb::QkTable &qk, const std::vector<DiracSpinor> &core,
               const std::vector<DiracSpinor> &excited,
               const Angular::SixJTable &SixJ,
               Denominators denominators = Denominators::BW);

//! Selection rule for Sk_vwxy (differs from Qk_vwxy due to parity)
bool Sk_vwxy_SR(int k, const DiracSpinor &v, const DiracSpinor &w,
                const DiracSpinor &x, const DiracSpinor &y);

//! Minimum/maximum k allowed by selectrion rules for Sk_vwxy. Cannot +=2.
std::pair<int, int> k_minmax_S(const DiracSpinor &v, const DiracSpinor &w,
                               const DiracSpinor &x, const DiracSpinor &y);

//! Matrix element of 1-body Sigma (2nd-order correlation) operator;
//! de_v = <v|Sigma|v>.
/*! @details 
 - ev is optional energy at which Sigma is calculated.
 - If ev is not given, uses 0.5*(ev+ew)
 - max_l_internal is largest angular momentum l to include in summation in
   internal lines; only used for tests.
 - qk (CoulombIntegral) may be YkTable or QkTable.
*/
template <class CoulombIntegral> // CoulombIntegral may be YkTable or QkTable
double Sigma_vw(const DiracSpinor &v, const DiracSpinor &w,
                const CoulombIntegral &qk, const std::vector<DiracSpinor> &core,
                const std::vector<DiracSpinor> &excited,
                //const Angular::SixJTable *const SixJ = nullptr,
                int max_l_internal = 99,
                std::optional<double> ev = std::nullopt);

//==============================================================================
//==============================================================================

namespace InternalSigma {
//! Diagram a for reduced two-body Sigma (sum of six Goldstone diagrams)
double S_Sigma2_a(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable &SixJ, Denominators denominators);

//! Diagram b for reduced two-body Sigma
double S_Sigma2_b(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable &SixJ, Denominators denominators);

//! Diagram c for reduced two-body Sigma
double S_Sigma2_c(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable &SixJ, Denominators denominators);

//! Diagram d for reduced two-body Sigma (no choice for denominators)
double S_Sigma2_d(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable &SixJ);
} // namespace InternalSigma

} // namespace MBPT

//==============================================================================
#include "Sigma2.ipp"