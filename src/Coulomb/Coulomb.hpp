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

//! Breit b^k function: (0,r) and (r,inf) part stored sepperately (in/out)
void bk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &b0, std::vector<double> &binf,
           const std::size_t maxi = 0);

//! Breit g^k function: (0,r) + (r,inf) part stored together (in/out)
void gk_ab(const DiracSpinor &Fa, const DiracSpinor &Fb, const int k,
           std::vector<double> &g0, std::vector<double> &ginf,
           const std::size_t maxi = 0);

//******************************************************************************

//! Calculates R^k_abcd for given k. From scratch (calculates y)
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Overload for when y^k_bd already exists [much faster]
double Rk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fc,
               const std::vector<double> &ykbd);

//! "Right-hand-side" R^k{v}_bcd [i.e., without Fv integral]
DiracSpinor Rkv_bcd(const int kappa_v, const DiracSpinor &Fb,
                    const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);
//! Overload for when y^k_bd already exists [much faster]
DiracSpinor Rkv_bcd(const int kappa_v, const DiracSpinor &Fc,
                    const std::vector<double> &ykbd);
//! Overload for when spinor exists. Rkv is overwritten
void Rkv_bcd(DiracSpinor *const Rkv, const DiracSpinor &Fc,
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

// void Qkv_bcd(DiracSpinor *Qkv, const DiracSpinor &Fb, const DiracSpinor &Fc,
//              const DiracSpinor &Fd, const int k);
DiracSpinor Qkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k);

//! Overload for when spinor exists. Qkv is overwritten
void Qkv_bcd(DiracSpinor *const Qkv, const DiracSpinor &Fb,
             const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
             const std::vector<double> &ykbd, const Angular::Ck_ab &Ck);

//! Exchange only version of W (W-Q): W = Q + P
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Exchange only version of W (W-Q): W = Q + P
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj);

//! Exchange only version of W (W-Q): W = Q + P
void Pkv_bcd(DiracSpinor *const Pkv, const DiracSpinor &Fb,
             const DiracSpinor &Fc, const DiracSpinor &Fd, const int k,
             const std::vector<std::vector<double>> &ybc,
             const Angular::Ck_ab &Ck, const Angular::SixJ &sixj);

DiracSpinor Pkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k,
                    const std::vector<std::vector<double>> &ybc,
                    const Angular::Ck_ab &Ck, const Angular::SixJ &sixj);

DiracSpinor Pkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k);

// Includes screening factor in a dumb way..
void Pkv_bcd_2(DiracSpinor *Pkv, const DiracSpinor &Fb, const DiracSpinor &Fc,
               const DiracSpinor &Fd, const int k,
               const std::vector<std::vector<double>> &ybc,
               const Angular::Ck_ab &Ck, const Angular::SixJ &sixj,
               const std::vector<double> &f2k);

//! Calculates W^k_abcd for given k. From scratch (calculates y)
//! @details
//! \f[ W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc \f]
//! \f[ W^k_abcd = Q^k_abcd + \P^k_abcd  \f]
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

DiracSpinor Wkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k);

//! Calculates Z^k_abcd for given k. From scratch (calculates y)
//! @details
//! \f[ Z^k_abcd = (-1)^{ja+jb+1} * ( Q^k_abcd + \sum_l [k] 6j * Q^l_abdc)
//! \f]
double Zk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Returns min and max k (multipolarity) allowed for C^k_ab
std::pair<int, int> k_minmax(const DiracSpinor &a, const DiracSpinor &b);
std::pair<int, int> k_minmax_tj(int tja, int tjb);

//! Returns min and max k (multipolarity) allowed for Q^k_abcd, Wk, Pk. For Q
//! (and Q only), parity rule is included, so you may safely call k+=2
std::pair<int, int> k_minmax_Q(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d);
//! DOES NOT contain parity rules (6j only) - so NOT safe to call k+=2
std::pair<int, int> k_minmax_P(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d);
std::pair<int, int> k_minmax_P(int kappa_a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d);
//! DOES NOT contain parity rules (6j only) - so NOT safe to call k+=2
std::pair<int, int> k_minmax_W(const DiracSpinor &a, const DiracSpinor &b,
                               const DiracSpinor &c, const DiracSpinor &d);

} // namespace Coulomb
