#pragma once
#include "Physics/PhysConst_constants.hpp"
#include "qip/Maths.hpp"
#include <array>
#include <cassert>
#include <cmath>

namespace DiracODE {

//! Performs asymptotic expansion for f and g at large r, up to order Nx in (1/r)
template <std::size_t Nx = 15> class AsymptoticSpinor {
private:
  int kappa;
  double Zeff, en, alpha, eps_target;
  double kappa2, alpha2, c, c2, lambda, sigma, Ren;
  double m_mass;
  std::array<double, Nx> bx; // bx must be first
  std::array<double, Nx> ax; // ax depends on bx

public:
  AsymptoticSpinor(int in_kappa, double in_Zeff, double in_en,
                   double in_alpha = PhysConst::alpha,
                   double in_eps_target = 1.0e-14, double m = 1.0)
      : kappa(in_kappa),
        Zeff(std::max(in_Zeff, 1.0)),
        en(in_en),
        alpha(in_alpha),
        eps_target(in_eps_target),
        kappa2(double(kappa * kappa)),
        alpha2(alpha * alpha),
        c(1.0 / alpha),
        c2(c * c),
        lambda(std::sqrt(-en * (2.0 + en * alpha2 / m))),
        sigma((m + en * alpha2) * (Zeff / lambda)),
        Ren(en + m * c2),
        m_mass(m),
        bx(make_bx()),
        ax(make_ax()) {
    // assert(en < 0.0 && "Must have en<0 in AsymptoticSpinor");
  }

  //! Performs asymptotic expansion for f and g at large r, up to order Nx in (1/r)
  /*!
   Large-r expansion of upper/lower radial components of the Dirac solution,
   see Johnson (2007), Eqs. (2.170) -- (2.171).

   f(r) = r^s exp(-yr) * { A(1 + O(1/r) + ...) + B(O(1/r) + ...)},

   g(r) = r^s exp(-yr) * { -B(1 + O(1/r) + ...) + A(O(1/r) + ...)},

   where s~1, y~1, A~1, B<<1.

   It's the 1/r expansion inside the {} brackets that is truncated at order Nx.
   The series is terminated if relative change drops below eps_target.
   This usually hapens around order ~5.
  */
  std::pair<double, double> fg(double r) const {
    // See Johnson (2007), Eqs. (2.170) -- (2.171)
    // Notation difference:
    // P(r) = f(r)
    // Q(r) = -g(r)
    // There appears to by typo in Eq. (2.171)

    const double A_large = std::sqrt(1.0 + 0.5 * en * alpha2 / m_mass);
    const double A_small = std::sqrt(-0.5 * en) * alpha;

    const double rfac = 2.0 * std::pow(r, sigma) * std::exp(-lambda * r);
    double fs = 1.0;
    double gs = 0.0;
    double rk = 1.0;
    // Continue the expansion until reach eps, or Nx
    for (std::size_t k = 0; k < Nx; k++) {
      rk *= r;
      fs += (ax[k] / rk);
      gs += (bx[k] / rk);
      const auto eps =
          std::max(std::abs(ax[k] / fs), std::abs(bx[k] / gs)) / rk;
      if (eps < eps_target) {
        break;
      }
    }
    return {rfac * (A_large * fs + A_small * gs),
            rfac * (A_large * gs - A_small * fs)};
  }

private:
  std::array<double, Nx> make_bx() const {
    // See Johnson (2007), Eqs. (2.172) -- (2.173)
    std::array<double, Nx> tbx;
    const auto Zalpha2 = Zeff * Zeff * alpha2;
    tbx[0] = (kappa + (Zeff / lambda)) * (0.5 * alpha);
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
    const auto RenAlpha2 = 1.0 + en * alpha2;
    for (std::size_t i = 0; i < Nx; i++) {
      tax[i] = (kappa + (double(i + 1) - sigma) * RenAlpha2 -
                Zeff * lambda * alpha2) *
               (bx[i] * c) / (double(i + 1) * lambda);
    }
    return tax;
  }
};

} // namespace DiracODE