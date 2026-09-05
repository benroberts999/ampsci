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
  @brief 
  Performs asymptotic expansion for f and g at large r, up to order Nx in (1/r).
  
  @details 
  Templated on energy type T: T=double for bound states, or
  T=std::complex<double> for Green's function solutions at complex energy.
  The expansion coefficients, lambda, and sigma extend analytically; the
  principal branch of sqrt gives Re(lambda)>0, i.e. the decaying solution.

  The branch may instead be chosen explicitly via @ref with_lambda 
  (e.g. lambda = i p for the oscillating en > 0 tail, see
  @ref AsymptoticSpinorContinuum).
  Everything else (sigma, the small-component amplitude, 
  and the 1/r coefficients) is derived from lambda.
*/
template <typename T = double, std::size_t Nx = 15>
class AsymptoticSpinor {
private:
  // Selects the constructor that takes lambda explicitly (see with_lambda)
  struct ExplicitLambda {};

  int kappa;
  double Zeff;
  T en;
  double alpha, m_mass, eps_target;
  double kappa2, alpha2, c;
  T lambda, sigma;
  // Large/small component amplitudes: A_large = sqrt(1 + en alpha^2/(2m)),
  // A_small = sqrt(-en/(2m)) alpha. Since lambda = sqrt(-2 m en) A_large,
  // A_small is written in terms of lambda so both share one branch of sqrt(-en)
  T A_large, A_small;
  // bx must be first (ax depends on bx in initialisation), make_ax()
  std::array<T, Nx> bx;
  std::array<T, Nx> ax;

  AsymptoticSpinor(ExplicitLambda, int in_kappa, double in_Zeff, T in_en,
                   T in_lambda, double in_alpha, double in_eps_target, double m)
    : kappa(in_kappa),
      Zeff(in_Zeff),
      en(in_en),
      alpha(in_alpha),
      m_mass(m),
      eps_target(in_eps_target),
      kappa2(double(kappa * kappa)),
      alpha2(alpha * alpha),
      c(1.0 / alpha),
      lambda(in_lambda),
      sigma((m + en * alpha2) * (Zeff / lambda)),
      A_large(std::sqrt(1.0 + 0.5 * en * alpha2 / m_mass)),
      A_small(lambda * alpha / (2.0 * m_mass * A_large)),
      bx(make_bx()),
      ax(make_ax()) {}

public:
  AsymptoticSpinor(int in_kappa, double in_Zeff, T in_en,
                   double in_alpha = PhysConst::alpha,
                   double in_eps_target = 1.0e-14, double m = 1.0)
    : AsymptoticSpinor(
        ExplicitLambda{}, in_kappa, in_Zeff, in_en,
        std::sqrt(-in_en * (2.0 * m + in_en * in_alpha * in_alpha)), in_alpha,
        in_eps_target, m) {
    // assert(en < 0.0 && "Must have en<0 in AsymptoticSpinor");
  }

  /*!
    @brief Constructs with an explicitly chosen lambda (exponent in
    exp(-lambda r)), i.e. with a chosen branch of sqrt(-en(2m + en alpha^2)).
    @details in_lambda must satisfy lambda^2 = -en(2m + en alpha^2); only its
    sign is free. Used for the en > 0 tail, where lambda = +/- i p and the
    principal branch is ambiguous (it sits on the branch cut).
  */
  static AsymptoticSpinor with_lambda(int in_kappa, double in_Zeff, T in_en,
                                      T in_lambda,
                                      double in_alpha = PhysConst::alpha,
                                      double in_eps_target = 1.0e-14,
                                      double m = 1.0) {
    return AsymptoticSpinor(ExplicitLambda{}, in_kappa, in_Zeff, in_en,
                            in_lambda, in_alpha, in_eps_target, m);
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
  std::pair<T, T> fg(double r) const {
    // See Johnson (2007), Eqs. (2.170) -- (2.171)
    // Notation difference:
    // P(r) = f(r)
    // Q(r) = -g(r)
    // There appears to by typo in Eq. (2.171)

    const T rfac = /*2.0 * */ std::pow(r, sigma) * std::exp(-lambda * r);
    T fs{1.0};
    T gs{0.0};
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
  std::array<T, Nx> make_bx() const {
    // See Johnson (2007), Eqs. (2.172) -- (2.173)
    std::array<T, Nx> tbx;
    const auto Zalpha2 = Zeff * Zeff * alpha2;
    tbx[0] = (kappa / m_mass + (Zeff / lambda)) * (0.5 * alpha);
    for (std::size_t i = 1; i < Nx; i++) {
      tbx[i] = (kappa2 - qip::pow<2>((double(i) - sigma)) - Zalpha2) *
               tbx[i - 1] / (double(2 * i) * lambda);
    }
    return tbx;
  }

  std::array<T, Nx> make_ax() const {
    // See Johnson (2007), Eq. (2.174)
    // bx must already be initialised
    std::array<T, Nx> tax;
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
  @brief Pair of oscillating continuum spinor values at one radius: the regular
  (F) and irregular (G) large-r Dirac-Coulomb solutions.
  @details Each spinor is stored as {f, g} (large, small components).
  F^C has large component ~ cos(theta), G^C has large component ~ sin(theta):
  the two are a quarter-wave (90-degree) pair.
*/
struct ContinuumTailSpinors {
  //! F^C = {fC, gC}, large ~ cos
  double fC, gC;
  //! G^C = {fG, gG}, large ~ sin
  double fG, gG;
};

/*!
  @brief Large-r Dirac-Coulomb oscillating tail spinors of a continuum (en > 0)
  state, energy-normalised.
  @details
  The en > 0 continuation of @ref AsymptoticSpinor. For a bound state the decay
  constant lambda = sqrt(-en(2 + en alpha^2/m)) is real and the spinor decays
  as r^sigma exp(-lambda r). For en > 0, lambda -> i p, with the relativistic
  momentum
  \f[ p = \sqrt{\en(2 + \en\alpha^2/m)} = \sqrt{\en(\en + 2c^2)}/c , \f]
  so the solution oscillates as exp(-i(pr + nu ln r)), with sigma = -i nu and
  nu = (m + en alpha^2) Z_ion / p. The expansion coefficients obey the same
  recurrence as the bound case, with complex values. The real and imaginary
  parts of the resulting complex spinor are the two real, linearly independent
  oscillating solutions
  \f[
    F^C \to \begin{pmatrix} A_L\cos\theta \\ -A_S\sin\theta\end{pmatrix},
    \qquad
    G^C \to \begin{pmatrix} -A_L\sin\theta \\ -A_S\cos\theta\end{pmatrix},
    \qquad A_S = \beta A_L,
  \f]
  with beta = sqrt(en/(en + 2c^2)). They are scaled to the energy
  normalisation A_L = sqrt(alpha/(pi beta)), so that the Wronskian is
  W[F^C, G^C] = f^C g^G - f^G g^C = -alpha/pi exactly.

  Implemented as the complex-energy @ref AsymptoticSpinor on the branch
  lambda = +i p (chosen explicitly: the principal sqrt is ambiguous on the
  branch cut), with the energy-normalisation scale applied on output.

  Used to seed the inward integration of the irregular continuum solution
  (@ref solveContinuumIrregular), as @ref AsymptoticSpinor seeds the decaying
  bound solution.

  @warning Valid only for en > 0 and at large r (where the 1/r expansion has
           converged); for en < 0 use @ref AsymptoticSpinor.
*/
template <std::size_t Nx = 15>
class AsymptoticSpinorContinuum {
private:
  using Complex = std::complex<double>;
  // relativistic momentum p = sqrt(en(2m + en alpha^2))
  double p;
  // small/large amplitude ratio beta = sqrt(en/(en + 2mc^2))
  double m_beta;
  // energy-normalised large-component amplitude A_L = sqrt(alpha/(pi beta))
  double m_amplitude_large;
  // Scale taking the raw expansion (large-component envelope
  // sqrt(1 + en alpha^2/(2m))) to the energy-normalised amplitude A_L
  double m_scale;
  // Complex expansion on the branch lambda = +i p, oscillating as
  // exp(-i(pr + nu ln r)); F^C and G^C are its real and imaginary parts
  AsymptoticSpinor<Complex, Nx> m_expansion;

public:
  AsymptoticSpinorContinuum(int kappa, double Zeff, double en,
                            double alpha = PhysConst::alpha,
                            double eps_target = 1.0e-14, double m = 1.0)
    : p(std::sqrt(en * (2.0 * m + en * alpha * alpha))),
      m_beta(std::sqrt(en * alpha * alpha / (2.0 * m + en * alpha * alpha))),
      m_amplitude_large(std::sqrt(alpha / (M_PI * m_beta))),
      m_scale(m_amplitude_large /
              std::sqrt(1.0 + 0.5 * en * alpha * alpha / m)),
      m_expansion(AsymptoticSpinor<Complex, Nx>::with_lambda(
        kappa, Zeff, Complex{en}, Complex{0.0, p}, alpha, eps_target, m)) {
    assert(en > 0.0 && "Must have en>0 in AsymptoticSpinorContinuum");
  }

  //! Relativistic momentum p = sqrt(en(en+2c^2))/c.
  double momentum() const { return p; }
  //! Energy-normalised large-component amplitude A_L = sqrt(alpha/(pi*beta)).
  double amplitude_large() const { return m_amplitude_large; }
  //! Small/large amplitude ratio beta = sqrt(en/(en+2c^2)).
  double beta() const { return m_beta; }

  /*!
    @brief Returns the two real oscillating tail spinors {F^C, G^C} at r.
    @details
    F^C (large ~ cos) and G^C (large ~ sin) are the real and imaginary parts
    of the complex en > 0 asymptotic spinor, energy-normalised. The 1/r series
    is truncated at order Nx, or when the relative change drops below the
    eps_target supplied at construction.
  */
  ContinuumTailSpinors fg(double r) const {
    const auto [f, g] = m_expansion.fg(r);
    return {m_scale * f.real(), m_scale * g.real(), m_scale * f.imag(),
            m_scale * g.imag()};
  }
};

} // namespace DiracODE
