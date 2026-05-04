#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <vector>

// Functions for computing Breit-related Coulomb integrals

namespace Coulomb {

//! Breit b^k function: (0,r) and (r,inf) part stored sepperately (in/out)
void bk_ab(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
           std::vector<double> &b0, std::vector<double> &binf,
           const std::size_t maxi = 0);

//! Breit g^k function: (0,r) + (r,inf) part stored together (in/out)
void gk_ab(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
           std::vector<double> &g0, std::vector<double> &ginf,
           const std::size_t maxi = 0);

/*!
  @brief Frequency-dependent Breit screening function g^k_ab(r,w), X_ab density.
  @details
  Computes the frequency-dependent Breit screening function using the X_ab
  radial density [arXiv:2602.17129](https://arxiv.org/abs/2602.17129):
  \f[ X_{ab}(r) = f_a(r)\,g_b(r) + g_a(r)\,f_b(r). \f]
  The static radial kernel \f$r_<^k/r_>^{k+1}\f$ is replaced by the
  frequency-dependent form:
  \f[ \frac{r_<^k}{r_>^{k+1}} \to -\omega(2k+1)\,j_k(\omega r_<)\,y_k(\omega r_>), \f]
  where \f$j_k\f$ and \f$y_k\f$ are spherical Bessel functions of the first and
  second kind.
  This screening function enters the frequency-dependent \f$u^k_{abcd}\f$
  integral.
  The \f$(0,r)\f$ and \f$(r,\infty)\f$ parts of the inner integral are stored
  separately in @p g0 and @p ginf.

  Refer to [arXiv:2602.17129](https://arxiv.org/abs/2602.17129) for details.
*/
void gk_ab_freqw(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 std::vector<double> &g0, std::vector<double> &ginf,
                 const std::size_t maxi = 0, const double w = 0.0);

/*!
  @brief Frequency-dependent Breit screening function h^k_ab(r,w), Y_ab density.
  @details
  Same as gk_ab_freqw but uses the Y_ab radial density
  [arXiv:2602.17129](https://arxiv.org/abs/2602.17129):
  \f[ Y_{ab}(r) = f_a(r)\,g_b(r) - g_a(r)\,f_b(r), \f]
  with the same frequency-dependent kernel replacement:
  \f[ \frac{r_<^k}{r_>^{k+1}} \to -\omega(2k+1)\,j_k(\omega r_<)\,y_k(\omega r_>). \f]
  Together with gk_ab_freqw, this is used to construct the \f$P^k\f$ and
  \f$Q^k\f$ screening functions that enter the
  \f$s^k_{abcd}\f$ and \f$t^k_{abcd}\f$ radial integrals.
  The \f$(0,r)\f$ and \f$(r,\infty)\f$ parts are stored separately in @p g0 and @p ginf.

  Refer to [arXiv:2602.17129](https://arxiv.org/abs/2602.17129) for details.
*/
void hk_ab_freqw(const int k, const DiracSpinor &Fa, const DiracSpinor &Fb,
                 std::vector<double> &g0, std::vector<double> &ginf,
                 const std::size_t maxi = 0, const double w = 0.0);

/*!
  @brief Frequency-dependent Breit v^k_ab(r,w): four partial screening integrals.
  @details
  Computes the frequency-dependent replacement of the \f$v^k_{abcd}\f$ radial
  integral [arXiv:2602.17129](https://arxiv.org/abs/2602.17129).
  The \f$P^k_{ij}\f$ and \f$Q^k_{ij}\f$ densities are
  constructed internally from @p Fi and @p Fj:
  \f[
    P^k_{ij}(r) = \frac{\kappa_i-\kappa_j}{k}\,X_{ij}(r) - Y_{ij}(r), \qquad
    Q^k_{ij}(r) = \frac{\kappa_i-\kappa_j}{k+1}\,X_{ij}(r) + Y_{ij}(r).
  \f]

  The four outputs correspond to the four double-integral terms:
  \f{align*}{
    &-2\!\int_0^\infty\!\!dr_1\!\int_0^{r_1}\!\!dr_2
      \left\{\left[\omega j_{k-1}(\omega r_<)y_{k+1}(\omega r_>)
              + \tfrac{2k+1}{\omega^2}\tfrac{r_<^{k-1}}{r_>^{k+2}}\right]
             Q^k_{ac}(r_1)P^k_{ij}(r_2)
             + \omega j_{k+1}(\omega r_<)y_{k-1}(\omega r_>)
               P^k_{ac}(r_1)Q^k_{ij}(r_2)\right\}\\
    &-2\!\int_0^\infty\!\!dr_1\!\int_{r_1}^{\infty}\!\!dr_2
      \left\{\left[\omega j_{k-1}(\omega r_<)y_{k+1}(\omega r_>)
              + \tfrac{2k+1}{\omega^2}\tfrac{r_<^{k-1}}{r_>^{k+2}}\right]
             P^k_{ac}(r_1)Q^k_{ij}(r_2)
             + \omega j_{k+1}(\omega r_<)y_{k-1}(\omega r_>)
               Q^k_{ac}(r_1)P^k_{ij}(r_2)\right\},
  \f}
  split as screening functions of \f$r_1\f$ (with \f$r_2\f$ integrated out):
  - @p v1: \f$(0,r)\f$ integral over \f$P^k_{ij}\f$; contracted externally with \f$Q^k_{ac}(r)\f$
  - @p v2: \f$(0,r)\f$ integral over \f$Q^k_{ij}\f$; contracted externally with \f$P^k_{ac}(r)\f$
  - @p v3: \f$(r,\infty)\f$ integral over \f$Q^k_{ij}\f$; contracted externally with \f$P^k_{ac}(r)\f$
  - @p v4: \f$(r,\infty)\f$ integral over \f$P^k_{ij}\f$; contracted externally with \f$Q^k_{ac}(r)\f$

  Refer to [arXiv:2602.17129](https://arxiv.org/abs/2602.17129) for details.

  @note For \f$\omega \leq \alpha\f$, a low-\f$\omega\f$ Taylor expansion of the
        Bessel functions is used to avoid numerical cancellation from divergent
        leading-order Bessel behaviour as \f$\omega\to 0\f$.
*/
void vk_ab_freqw(const int k, const DiracSpinor &Fi, const DiracSpinor &Fj,
                 const Grid &gr, std::vector<double> &v1,
                 std::vector<double> &v2, std::vector<double> &v3,
                 std::vector<double> &v4, std::size_t maxi, const double w);

} // namespace Coulomb
