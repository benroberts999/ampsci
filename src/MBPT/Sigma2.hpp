#pragma once
#include "Angular/include.hpp"
#include "Coulomb/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/String.hpp"
#include <string>
#include <vector>

namespace MBPT {

/*! @brief Type of energy demoninators: RS, Fermi, Fermi0

 - RS     : Use energy denominator ascociated with actual orbital for external legs. Symmeterised so that diagram is symmetric. May be danger of accidental enhancement.
 - Fermi  : Use energy from the "Fermi" level (i.e., lowest state for given kappa in excited spectrum).
 - Fermi0 : As above, but assume Fermi level for all kappas the same. These often cancel, so there is no (excited-excited) term in denominator (except diagram d). Fine, since the remaining core-excited always dominates.

Energy for internal legs (hole-particle) always actual orbtials.
*/
enum class Denominators { RS, Fermi, Fermi0 };

//! Returns string representation of Denominators enum
std::string parse_Denominators(Denominators d);

//! Parses string to Denominators enum (case-insensitive); returns Fermi0 if unrecognised
Denominators parse_Denominators(std::string_view s);

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
               Denominators denominators = Denominators::Fermi0);

//! Selection rule for Sk_vwxy (differs from Qk_vwxy due to parity)
bool Sk_vwxy_SR(int k, const DiracSpinor &v, const DiracSpinor &w,
                const DiracSpinor &x, const DiracSpinor &y);

//! Minimum/maximum k allowed by selectrion rules for Sk_vwxy. Cannot +=2.
std::pair<int, int> k_minmax_S(const DiracSpinor &v, const DiracSpinor &w,
                               const DiracSpinor &x, const DiracSpinor &y);

std::pair<int, int> k_minmax_S(int twoj_v, int twoj_w, int twoj_x, int twoj_y);

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
                int max_l_internal = 99,
                std::optional<double> ev = std::nullopt);

//! Returns energy of first state in excited that matches given kappa
double e_bar(int kappa_v, const std::vector<DiracSpinor> &excited);

//==============================================================================
//==============================================================================

namespace InternalSigma {
//! Diagrams a+b for reduced two-body Sigma (sum of six Goldstone diagrams)
double S_Sigma2_ab(int k, const DiracSpinor &v, const DiracSpinor &w,
                   const DiracSpinor &x, const DiracSpinor &y,
                   const Coulomb::QkTable &qk,
                   const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Angular::SixJTable &SixJ, Denominators denominators);

//! Diagram b for reduced two-body Sigma
double S_Sigma2_c1(int k, const DiracSpinor &v, const DiracSpinor &w,
                   const DiracSpinor &x, const DiracSpinor &y,
                   const Coulomb::QkTable &qk,
                   const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Angular::SixJTable &SixJ, Denominators denominators);

//! Diagram c for reduced two-body Sigma
double S_Sigma2_c2(int k, const DiracSpinor &v, const DiracSpinor &w,
                   const DiracSpinor &x, const DiracSpinor &y,
                   const Coulomb::QkTable &qk,
                   const std::vector<DiracSpinor> &core,
                   const std::vector<DiracSpinor> &excited,
                   const Angular::SixJTable &SixJ, Denominators denominators);

//! Diagram d for reduced two-body Sigma
double S_Sigma2_d(int k, const DiracSpinor &v, const DiracSpinor &w,
                  const DiracSpinor &x, const DiracSpinor &y,
                  const Coulomb::QkTable &qk,
                  const std::vector<DiracSpinor> &core,
                  const std::vector<DiracSpinor> &excited,
                  const Angular::SixJTable &SixJ, Denominators denominators);
} // namespace InternalSigma

} // namespace MBPT

//==============================================================================
#include "Sigma2.ipp"