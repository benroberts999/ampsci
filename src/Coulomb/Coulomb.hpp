#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>
namespace Angular {
class Ck_ab;
class SixJ;
} // namespace Angular

//! Functions (+classes) for computing Coulomb integrals
namespace Coulomb {

//! Calculates Hartree Screening functions \f$y^k_{ab}(r)\f$
//! @details maxi is max point to calculate; blank or zero means all the way
std::vector<double> yk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb,
                          const int k, const std::size_t maxi = 0);
//! Overload: does not allocate ykab
void yk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &ykab, const std::size_t maxi = 0);

//! Calculates R^k_abcd for given k. From scratch (calculates y)
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);
//! Overload for when y^k_bd already exists [much faster]
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fc,
               const std::vector<double> &ykbd);

//! "Right-hand-side" R^k_{a}bcd [i.e., without Fa integral]
//! @details Fa only used for kappa: returns spinor with k = Fa.k
DiracSpinor Rk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                        const DiracSpinor &Fc, const DiracSpinor &Fd,
                        const int k);
//! Overload for when y^k_bd already exists [much faster]
DiracSpinor Rk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fc,
                        const std::vector<double> &ykbd);
//! Overload: Rk_rhs depends on kappa_a, not Fa
DiracSpinor Rk_abcd_rhs(const int kappa_a, const DiracSpinor &Fc,
                        const std::vector<double> &ykbd);
//! Overload for when spinor exists. Rkv is overwritten
void Rk_abcd_rhs(DiracSpinor *const Rkv, const DiracSpinor &Fc,
                 const std::vector<double> &ykbd);

//! Calculates Q^k_abcd for given k. From scratch (calculates y)
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Overload for when Coulomb int + angular is known; faster
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<double> &ykbd, const Angular::Ck_ab &Ck);

//! Calculates Qkv_bcd (radial spinor).
DiracSpinor Qkv_bcd(const int kappa_v, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
                    const std::vector<double> &ykbd, const Angular::Ck_ab &Ck);

//! Overload for when spinor exists. Qkv is overwritten
void Qkv_bcd(DiracSpinor *const Qkv, const DiracSpinor &Fb,
             const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
             const std::vector<double> &ykbd, const Angular::Ck_ab &Ck);

//! Calculates Z^k_abcd for given k. From scratch (calculates y)
//! @details
//! \f[ Z^k_abcd = (-1)^{ja+jb+1} * ( Q^k_abcd + \sum_l [k] 6j * Q^l_abdc)
//! \f]
double Zk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Calculates W^k_abcd for given k. From scratch (calculates y)
//! @details
//! \f[ W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc \f]
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Exchange only version of W (W-Q): W = Q + P
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj);

//! Exchange only version of W (W-Q): W = Q + P
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Exchange only version of W (W-Q): W = Q + P
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj);

} // namespace Coulomb
