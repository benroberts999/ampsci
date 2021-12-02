#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

//! Flambaum-Ginges Radiative potential
/*!@details

First introduced in:
 * V. V. Flambaum and J. S. M. Ginges, Phys. Rev. A 72, 052115 (2005).

With finite-nuclear-size corrections, and updated fitting factors, from;
 * J. S. M. Ginges and J. C. Berengut, Phys. Rev. A 93, 052509 (2016).
 * J. S. M. Ginges and J. C. Berengut, J. Phys. B 49, 095001 (2016).

Note:
 * High- and low- frequency V(r) do NOT include A or B fitting coeficients.
 * All radii (r and rN) are given in atomic units
 * Sign convention: H -> H - V_rad
 * V_rad = V_Uehling + V_SEh + V_SEl + V_WK + i(g.n)H_Magnetic
*/

namespace FGRP {

// Fine structure constant:
constexpr double alpha = 1.0 / 137.035999084;
// electron reduced compton wavelength (atomic units)
constexpr double lam_c = alpha;

//******************************************************************************
//! Uehling potential (r, rN in atomic units)
double V_Uehling(double Z, double r, double rN = 0.0);

//! Magnetic form-factor (r, rN in atomic units)
double H_Magnetic(double Z, double r, double rN = 0.0);

//! High-freq electric SE (NOT including Al) (r, rN in atomic units)
double V_SEh(double Z, double r, double rN = 0.0);

//! Low-freq electric SE (NOT including Bl) (r, rN in atomic units)
double V_SEl(double Z, double r, double rN = 0.0);

//! Effective Wickman-Kroll; not including FNS corrections
double V_WK(double Z, double r, double);

//******************************************************************************
//! Fitting factors from Ginges Berengut, Phys. Rev. A 93, 052509 (2016).
namespace Fit {
// static const auto a01 = std::vector{1.071, 0.0, -1.976, -2.128, 0.169};
static const auto a01 = std::vector{1.071, 0.0, -1.976, -2.128, 0.169};
static const auto b01 = std::vector{0.074, 0.35};
static const auto b2 = std::vector{0.056, 0.050, 0.195};

//! Al(Z) fitting function [PRA 93, 052509 (2016)]
double A(double Z, int l = 0);
//! Bl(Z) fitting function [PRA 93, 052509 (2016)]
double B(double Z, int l = 0);

} // namespace Fit

//******************************************************************************
//! Function that performs the t integral over [1,infinity)
double t_integral(double (*f)(double, void *), std::vector<double> params,
                  double eps = 1.0e-6);

enum class IntType { Linear, Log };

//! Function that performs the r integral over [a,b] using mixed-quadrature
template <IntType IT = IntType::Linear>
double r_integral(std::function<double(double)> f, double a, double b,
                  long unsigned n_pts = 1000);

//******************************************************************************
// Helper functions:

namespace Uehling {
double G_Ueh(double xn, double x);
// Form of function required by GSL integrations; p is {r, Rn}
double J_Ueh_gsl(double t, void *p);
} // namespace Uehling

namespace Magnetic {
double G_mag(double xn, double x);
// Form of function required by GSL integrations; p is {r, Rn}
double J_mag_gsl(double t, void *p);
} // namespace Magnetic

namespace SE {
double xi(double x);
double F_SEl(double Z, double r, double rN);
double I1(double Z, double t);
double I2(double Z, double r, double rn, double t);
// Form of function required by GSL integrations; p is {r, Rn, Z}
double J_SE_gsl(double t, void *p);
} // namespace SE

//******************************************************************************
// Implementation:
template <IntType IT>
double r_integral(std::function<double(double)> f, double a, double b,
                  long unsigned n_pts) {

  static constexpr std::array w{475.0 / 1440, 1902.0 / 1440, 1104.0 / 1440,
                                1586.0 / 1440, 1413.0 / 1440};

  const auto dt = (IT == IntType::Linear) ? (b - a) / double(n_pts - 1) :
                                            std::log(b / a) / double(n_pts - 1);

  const auto x = [=](auto i) {
    if constexpr (IT == IntType::Linear) {
      return a + double(i) * dt;
    } else {
      assert(a != 0.0);
      return a * std::exp(double(i) * dt);
    }
  };
  const auto dxdt = [=](auto i) {
    if constexpr (IT == IntType::Linear) {
      (void)i; // suppress unused variable warning on older gcc compilers
      return 1.0;
    } else {
      return x(i); // dxdt = x(t) for log
    }
  };

  double Rint = 0.0;
  for (long unsigned i = 0; i < w.size(); i++) {
    Rint += w[i] * (f(x(i)) + f(x(n_pts - 1 - i))) * dxdt(i);
  }

  for (auto i = w.size(); i < n_pts - w.size(); i++) {
    Rint += f(x(i)) * dxdt(i);
  }

  return Rint * dt;
}

//******************************************************************************
//******************************************************************************
//! Magnetic-loop vacuum polarisation, includes finite-nuclear size
double Q_MLVP(double r, double rN = 0.0);
namespace Helper {
// Magnetic-loop vacuum polarisation, Q = Int(Qt, {t,1,infty}).
double Qt_MLVP(double t, void *p);
} // namespace Helper

// struct to pass the parameters to the GSL function
struct MLVP_params {
  // simple struct that only stores current point on the radial grid and the
  // nuclear radius
  double r, rN;
};

} // namespace FGRP
