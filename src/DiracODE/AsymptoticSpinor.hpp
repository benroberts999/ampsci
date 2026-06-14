#pragma once
#include "Physics/PhysConst_constants.hpp"
#include "qip/Maths.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <utility>

namespace DiracODE {

/*!
  @brief Performs asymptotic expansion for f and g at large r, up to order Nx in (1/r).
*/
template <std::size_t Nx = 15>
class AsymptoticSpinor {
private:
  int kappa;
  double Zeff, en, alpha, m_mass, eps_target;
  double kappa2, alpha2, c, lambda, sigma;
  std::array<double, Nx> bx; // bx must be first
  std::array<double, Nx> ax; // ax depends on bx

public:
  AsymptoticSpinor(int in_kappa, double in_Zeff, double in_en,
                   double in_alpha = PhysConst::alpha,
                   double in_eps_target = 1.0e-14, double m = 1.0)
    : kappa(in_kappa),
      Zeff(in_Zeff),
      en(in_en),
      alpha(in_alpha),
      m_mass(m),
      eps_target(in_eps_target),
      kappa2(double(kappa * kappa)),
      alpha2(alpha * alpha),
      c(1.0 / alpha),
      lambda(std::sqrt(-en * (2.0 * m_mass + en * alpha2))),
      sigma((m + en * alpha2) * (Zeff / lambda)),
      // Ren(en + m * c2),
      bx(make_bx()),
      ax(make_ax()) {
    // assert(en < 0.0 && "Must have en<0 in AsymptoticSpinor");
  }

  /*!
    @brief Returns {f(r), g(r)} via asymptotic expansion at large r.
    @details
    Large-r expansion of upper/lower radial components of the Dirac solution,
    see Johnson (2007), Eqs. (2.170) -- (2.171).

    f(r) = r^s exp(-yr) * { A(1 + O(1/r) + ...) + B(O(1/r) + ...)},

    g(r) = r^s exp(-yr) * { -B(1 + O(1/r) + ...) + A(O(1/r) + ...)},

    where s~1, y~1, A~1, B<<1.

    The 1/r expansion inside the braces is truncated at order Nx. The series is
    terminated early if the relative change drops below eps_target (typically
    around order ~5).
  */
  std::pair<double, double> fg(double r) const {
    // See Johnson (2007), Eqs. (2.170) -- (2.171)
    // Notation difference:
    // P(r) = f(r)
    // Q(r) = -g(r)
    // There appears to by typo in Eq. (2.171)

    const double A_large = std::sqrt(1.0 + 0.5 * en * alpha2 / m_mass);
    const double A_small = std::sqrt(-0.5 * en / m_mass) * alpha;

    const double rfac = /*2.0 * */ std::pow(r, sigma) * std::exp(-lambda * r);
    double fs = 1.0;
    double gs = 0.0;
    // Continue the expansion until reach eps, or Nx
    for (std::size_t k = 0; k < Nx; k++) {
      const auto rkp1 = qip::pow(r, int(k) + 1);
      const auto df = ax[k] / rkp1;
      const auto dg = bx[k] / rkp1;
      fs += df;
      gs += dg;
      const auto eps = std::max(std::abs(df / fs), std::abs(dg / gs));
      if (eps < eps_target) {
        break;
      }
    }
    // here: typo in Johnson, or not? Both work
    return {rfac * (A_large * fs + A_small * gs),
            rfac * (A_large * gs - A_small * fs)};
    // -rfac * (A_large * fs - A_small * gs)};
  }

private:
  std::array<double, Nx> make_bx() const {
    // See Johnson (2007), Eqs. (2.172) -- (2.173)
    std::array<double, Nx> tbx;
    const auto Zalpha2 = Zeff * Zeff * alpha2;
    tbx[0] = (kappa / m_mass + (Zeff / lambda)) * (0.5 * alpha);
    for (std::size_t i = 1; i < Nx; i++) {
      tbx[i] = (kappa2 - qip::pow<2>((double(i) - sigma)) - Zalpha2) *
               tbx[i - 1] / (double(2 * i) * lambda);
    }
    return tbx;
  }

  std::array<double, Nx> make_ax() const {
    // See Johnson (2007), Eq. (2.174)
    // bx must already be initialised
    std::array<double, Nx> tax;
    const auto RenAlpha2 = m_mass + en * alpha2;
    for (std::size_t i = 0; i < Nx; i++) {
      tax[i] = (kappa * m_mass + (double(i + 1) - sigma) * RenAlpha2 -
                Zeff * lambda * alpha2) *
               (bx[i] * c) / (double(i + 1) * lambda);
    }
    return tax;
  }
};

//==============================================================================

/*!
  @brief
  Pair of continuum spinors at a single r ('tail'), both regular (F) and irregular (G) flavours.
  @details
  Each spinor is stored as {f, g} (large, small radial components).
  - F^C ("regular" Coulomb): large component ~ cos(theta_kappa)
  - G^C ("irregular" Coulomb): large component ~ sin(theta_kappa)
  i.e. the two are a 90-degree (oscilating) pair.
*/
struct ContinuumTailSpinors {
  //! F^C = {fC, gC}, large ~ cos
  double fC, gC;
  //! G^C = {fG, gG}, large ~ sin
  double fG, gG;
};

/*!
  @brief Large-r Dirac-Coulomb oscilating tail spinors for a continuum
  (en>0) state, energy-normalised.
  @details
  This is the en>0 continuation of @ref AsymptoticSpinor. For a bound state the
  decay constant \f$ \lambda = \sqrt{(-\en(2 + \en*\alpha^2/m))}\f$
  is real and the spinor decays as r^sigma exp(-lambda r). 
  For a continuum state en>0, lambda becomes
  imaginary, lambda -> i*p with the relativistic momentum

  \f[ 
      p = \sqrt{\en(2 + \en\alpha^2/m)}
        = \frac{\sqrt{\en(\en + 2c^2)}}{c}, 
  \f]

  so r^sigma exp(-lambda r) -> exp(-i(pr + nu\ln r)) becomes oscillatory, with
  \f$ \sigma = -i\nu \f$, \f$ \nu = (m + \en\alpha^2)Z_{\rm ion}/p \f$. The
  expansion coefficients ax, bx and the small-component amplitude A_small all
  acquire imaginary parts, but obey the *same* recurrence as the bound case.
  The real and imaginary parts of the
  resulting complex spinor are the two real, linearly-independent oscilating
  solutions

  \f[
    F^C \to \begin{pmatrix} A_L\cos\theta \\ A_S\sin\theta\end{pmatrix},
    \qquad
    G^C \to \begin{pmatrix} -A_L\sin\theta \\ A_S\cos\theta\end{pmatrix},
    \qquad A_S = \beta A_L,
  \f]
  (up to an overall sign/phase convention) with the large/small ratio
  \f$ \beta = \sqrt{\en/(\en+2c^2)} \f$. They are scaled to the ampsci energy
  normalisation \f$ A_L = \sqrt{\alpha/(\pi\beta)} \f$ so that the conserved
  Wronskian is \f$ W[F^C,G^C] = f^C g^G - f^G g^C = -\alpha/\pi \f$, exactly
  (the sign reflects the chosen orientation of the pair).

  These are used to seed the inward integration of the irregular partner
  F_irr (method B: Coulomb match), exactly as the bound @ref AsymptoticSpinor
  seeds the decaying solution.

  There might be a way to do this without explict complex numbers?

  @note Energy-normalised assuming electron mass m=1 (the physical continuum
        case). The overall normalisation cancels in the method-B linear
        projection, so it is only essential for the W = -alpha/pi check.
  @warning Valid only for en>0 and at large r (where the 1/r expansion has
           converged); for en<0 use @ref AsymptoticSpinor.
*/
template <std::size_t Nx = 15>
class AsymptoticSpinorContinuum {
private:
  using Complex = std::complex<double>;
  int kappa;
  double Zeff, en, alpha, eps_target, m_mass;
  double kappa2, alpha2, c, c2;
  // relativistic momentum
  double p;
  // raw large-component amplitude
  double A_large;
  // |raw small-component amplitude|
  double As_mag;
  // energy-normalisation scale factor
  double m_scale;
  // = i*p
  Complex lambda;
  // = -i*nu
  Complex sigma;
  // bx must be first
  std::array<Complex, Nx> bx;
  // ax depends on bx
  std::array<Complex, Nx> ax;

public:
  AsymptoticSpinorContinuum(int in_kappa, double in_Zeff, double in_en,
                            double in_alpha = PhysConst::alpha,
                            double in_eps_target = 1.0e-14, double m = 1.0)
    : kappa(in_kappa),
      Zeff(in_Zeff),
      en(in_en),
      alpha(in_alpha),
      eps_target(in_eps_target),
      m_mass(m),
      kappa2(double(kappa * kappa)),
      alpha2(alpha * alpha),
      c(1.0 / alpha),
      c2(c * c),
      p(std::sqrt(en * (2.0 + en * alpha2 / m))),
      A_large(std::sqrt(1.0 + 0.5 * en * alpha2 / m)),
      As_mag(std::sqrt(0.5 * en / m) * alpha),
      // Energy-normalisation: scale so large-component envelope = A_L
      // A_L = sqrt(alpha/(pi*beta)) = 1/sqrt(pi*ceps); |f| envelope = 2*A_large.
      m_scale((1.0 / std::sqrt(M_PI * std::sqrt(en / (en * alpha2 + 2.0)))) /
              (2.0 * A_large)),
      lambda(0.0, p),
      sigma((m + en * alpha2) * (Zeff / lambda)),
      bx(make_bx()),
      ax(make_ax()) {
    assert(en > 0.0 && "Must have en>0 in AsymptoticSpinorContinuum");
  }

  //! Relativistic momentum p = sqrt(en(en+2c^2))/c.
  double momentum() const { return p; }
  //! Energy-normalised large-component amplitude A_L = sqrt(alpha/(pi*beta)).
  double amplitude_large() const { return 2.0 * A_large * m_scale; }
  //! Small/large amplitude ratio beta = sqrt(en/(en+2c^2)).
  double beta() const { return As_mag / A_large; }

  /*!
    @brief Returns the two real oscilating tail spinors {F^C, G^C} at r.
    @details
    F^C (large ~ cos) and G^C (large ~ sin) are the real and imaginary parts
    of the complex en>0 asymptotic spinor, energy-normalised. The 1/r series
    is truncated at order Nx or when the relative change drops below the
    eps_target supplied at construction.
  */
  ContinuumTailSpinors fg(double r) const {
    const Complex rfac = 2.0 * std::exp(sigma * std::log(r) - lambda * r);
    Complex fs{1.0, 0.0};
    Complex gs{0.0, 0.0};
    double rk = 1.0;
    for (std::size_t k = 0; k < Nx; k++) {
      rk *= r;
      fs += ax[k] / rk;
      gs += bx[k] / rk;
      const auto eps =
        std::max(std::abs(ax[k] / fs), std::abs(bx[k] / gs)) / rk;
      if (eps < eps_target) {
        break;
      }
    }
    // A_small = i*As_mag for en>0 (continuation of sqrt(-0.5 en/m)*alpha)
    const Complex As{0.0, As_mag};
    const Complex f = rfac * (A_large * fs + As * gs);
    const Complex g = rfac * (A_large * gs - As * fs);
    return {m_scale * f.real(), m_scale * g.real(), m_scale * f.imag(),
            m_scale * g.imag()};
  }

private:
  std::array<Complex, Nx> make_bx() const {
    std::array<Complex, Nx> tbx;
    const auto Zalpha2 = Zeff * Zeff * alpha2;
    tbx[0] = (double(kappa) + (Zeff / lambda)) * (0.5 * alpha);
    for (std::size_t i = 1; i < Nx; i++) {
      const auto im_s = double(i) - sigma;
      tbx[i] = (kappa2 - im_s * im_s - Zalpha2) * tbx[i - 1] /
               (double(2 * i) * lambda);
    }
    return tbx;
  }

  std::array<Complex, Nx> make_ax() const {
    std::array<Complex, Nx> tax;
    const auto RenAlpha2 = 1.0 + en * alpha2;
    for (std::size_t i = 0; i < Nx; i++) {
      tax[i] = (double(kappa) + (double(i + 1) - sigma) * RenAlpha2 -
                Zeff * lambda * alpha2) *
               (bx[i] * c) / (double(i + 1) * lambda);
    }
    return tax;
  }
};

} // namespace DiracODE