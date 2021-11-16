#pragma once
#include "Angular/Angular.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <optional>
#include <vector>

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

//******************************************************************************

//! Calculates Q^k_abcd for given k. From scratch (calculates y) [see YkTable
//! version if already have YkTable]
double Qk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Calculates Q^k(v)_bcd for given k,kappa_v. From scratch (calculates y) [see
//! YkTable version if already have YkTable]
DiracSpinor Qkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k);

//******************************************************************************

//! Exchange only version of W (W-Q): W = Q + P [see Qk above]
double Pk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

//! Exchange only version of W (W-Q): W = Q + P [see Qk above]
DiracSpinor Pkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k);

//******************************************************************************

//! Calculates W^k_abcd for given k. From scratch (calculates y)
//! @details
//! \f[ W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc \f]
//! \f[ W^k_abcd = Q^k_abcd + \P^k_abcd  \f]
double Wk_abcd(const DiracSpinor &Fa, const DiracSpinor &Fb,
               const DiracSpinor &Fc, const DiracSpinor &Fd, const int k);

DiracSpinor Wkv_bcd(int kappa_v, const DiracSpinor &Fb, const DiracSpinor &Fc,
                    const DiracSpinor &Fd, const int k);

//******************************************************************************

//! Returns min and max k (multipolarity) allowed for C^k_ab, accounting for
//! parity
std::pair<int, int> k_minmax(const DiracSpinor &a, const DiracSpinor &b);
//! Returns min and max k (multipolarity) allowed for C^k_ab, NOT accounting for
//! parity (2j only, not kappa/l)
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

//******************************************************************************

template <class A> static int twojk(const A &a) {
  if constexpr (std::is_same_v<A, DiracSpinor>) {
    return a.twoj();
  } else {
    static_assert(std::is_same_v<A, int>);
    return 2 * a;
  }
}

template <class A> static std::optional<int> twojknull(const A &a) {
  if constexpr (std::is_same_v<A, DiracSpinor>) {
    return a.twoj();
  } else if constexpr (std::is_same_v<A, int>) {
    static_assert(std::is_same_v<A, int>);
    return 2 * a;
  } else {
    return std::nullopt;
  }
}

template <class A, class B, class C, class D, class E, class F>
static double sixj(const A &a, const B &b, const C &c, const D &d, const E &e,
                   const F &f) {
  return Angular::sixj_2(twojk(a), twojk(b), twojk(c), twojk(d), twojk(e),
                         twojk(f));
}

template <class A = std::optional<int>, class B = std::optional<int>,
          class C = std::optional<int>, class D = std::optional<int>,
          class E = std::optional<int>, class F = std::optional<int>>
static bool sixjTriads(const A &a, const B &b, const C &c, const D &d,
                       const E &e, const F &f) {
  return Angular::sixjTriads(twojknull(a), twojknull(b), twojknull(c),
                             twojknull(d), twojknull(e), twojknull(f));
}

} // namespace Coulomb
