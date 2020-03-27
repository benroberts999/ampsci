#include "RadiativePotential.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_expint.h>
#include <iostream>
#include <vector>

static inline double ExpInt1(double x) {
  // E1(x) = -std::expint(-x)
  // return -std::expint(-x);
  return gsl_sf_expint_E1(x);
}

namespace RadiativePotential {
using namespace Helper;

std::vector<double> form_Hel(const std::vector<double> &r, double x_simple,
                             double x_Ueh, double x_SEe_h, double x_SEe_l,
                             double r_rms_Fermi, double z, double alpha,
                             double rcut) {

  const auto rN_rad = r_rms_Fermi * std::sqrt(5.0 / 3.0) / PhysConst::aB_fm;

  if (rcut == 0.0)
    rcut = r.back();
  const auto imax = r.size();
  std::vector<double> Vel_tmp(r.size());

  if (std::abs(x_simple) > 0) {
    std::cout << "Forming simple exponential potential (Ak=" << x_simple
              << ")\n";
    for (std::size_t i = 0; i < imax; ++i) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      Vel_tmp[i] -=
          x_simple * RadiativePotential::vSimpleExp(ri, rN_rad, z, alpha);
    }
  }

  if (x_Ueh > 0) {
    std::cout << "Forming Uehling potential (scale=" << x_Ueh << ")\n";
    // #pragma omp parallel for
    for (std::size_t i = 0; i < imax; ++i) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      const auto v_Ueh = RadiativePotential::vUehling(ri, rN_rad, z, alpha);
      Vel_tmp[i] -= x_Ueh * v_Ueh;
    }
  }

  if (x_SEe_h > 0) {
    std::cout << "Forming Self-Energy (electric high-f) potential (scale="
              << x_SEe_h << ")\n";
    // The SE high-f part is very slow. We calculate it every few points only,
    // and then interpolate the intermediate points
    const std::size_t stride = 6;
    const auto tmp_max = imax / stride;
    std::vector<double> x(tmp_max), y(tmp_max);
#pragma omp parallel for
    for (std::size_t i = 0; i < tmp_max; ++i) {
      x[i] = r[i * stride];
      if (x[i] > rcut)
        continue;
      y[i] = RadiativePotential::vSEh(x[i], rN_rad, z, alpha);
    }
    const auto vec_SE_h = Interpolator::interpolate(x, y, r);
    for (std::size_t i = 0; i < imax; i++) {
      if (r[i] > rcut)
        break;
      Vel_tmp[i] -= x_SEe_h * vec_SE_h[i];
    }
  }

  if (x_SEe_l > 0) {
    std::cout << "Forming Self-Energy (electric low-f) potential (scale="
              << x_SEe_l << ")\n";
    // #pragma omp parallel for
    for (std::size_t i = 0; i < imax; i++) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      const auto v_SE_l = RadiativePotential::vSEl(ri, rN_rad, z, alpha);
      Vel_tmp[i] -= x_SEe_l * v_SE_l;
    }
  }

  return Vel_tmp;
}
//------------------------------------------------------------------------------
std::vector<double> form_Hmag(const std::vector<double> &r, double x_SEm,
                              double r_rms_Fermi, double z, double alpha,
                              double rcut) {

  const auto rN_rad = r_rms_Fermi * std::sqrt(5.0 / 3.0) / PhysConst::aB_fm;

  if (rcut == 0)
    rcut = r.back();
  const auto imax = r.size();
  std::vector<double> Hmag_tmp(r.size());

  if (x_SEm > 0) {
    std::cout << "Forming Self-Energy (magnetic) potential "
              << "(scale=" << x_SEm << ")\n";
    // #pragma omp parallel for
    for (std::size_t i = 0; i < imax; i++) {
      const auto ri = r[i];
      if (ri > rcut)
        break;
      const auto Hr = RadiativePotential::vSE_Hmag(ri, rN_rad, z, alpha);
      Hmag_tmp[i] += x_SEm * Hr; // nb: plus!
    }
  }
  return Hmag_tmp;
}

//******************************************************************************
//******************************************************************************
double vSimpleExp(double r, double, double z, double alpha) {
  return -z * z * alpha * std::exp(-r / alpha);
}

//******************************************************************************
double vUehling(double r, double rN, double z, double alpha) {
  auto sp1 = SafeProfiler::profile(__func__);

  // Routines return the first approximation which has an absolute error
  // smaller than abs_err_lim or a relative error smaller than rel_err_lim.
  // Note that this is an either-or constraint, not simultaneous. To compute to
  // a specified absolute error, set epsrel to zero (etc.)
  static constexpr double abs_err_lim = 0.0;
  static constexpr double rel_err_lim = 1.0e-6;
  // max_num_subintvls < size(gsl_int_wrk)
  static constexpr unsigned long max_num_subintvls = 750; //?

  gsl_set_error_handler_off(); //?
  if (rN <= 0) {
    rN = 1.0e-7;
  }

  RadPot_params params = {r, rN, z, alpha};

  gsl_function f_gsl;
  f_gsl.function = (r < rN) ? &gslfunc_Ueh_smallr : &gslfunc_Ueh_larger;
  f_gsl.params = &params;

  // This workspace handles the memory for the subinterval ranges, results and
  // error estimates. max_num_subintvls < size(gsl_int_wrk)
  // Allocates a workspace sufficient to hold n double precision
  // intervals, their integration results and error estimates.
  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);
  // nb: i allocate + destroy wrk EACH r... doesn't much matter though...

  double int_result, abs_err;
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  auto pre_factor = z * (alpha / M_PI) / r;

  return pre_factor * int_result;
}

//******************************************************************************
static Fit_AB fit_AB;
double vSEh(double r, double rN, double z, double alpha) {
  auto sp1 = SafeProfiler::profile(__func__);

  static constexpr double abs_err_lim = 0.0;
  static constexpr double rel_err_lim = 1.0e-3;
  static constexpr unsigned long max_num_subintvls = 1000; //?

  gsl_set_error_handler_off(); //?
  if (rN <= 0) {
    rN = 1.0e-7;
  }

  RadPot_params params = {r, rN, z, alpha};

  gsl_function f_gsl;
  f_gsl.function = (r < rN) ? &gslfunc_SEh_smallr : &gslfunc_SEh_larger;
  f_gsl.params = &params;

  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);

  double int_result, abs_err;
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  auto l = 0; // XXX
  auto al = fit_AB.Al(l, z * alpha);

  auto pre_factor = -1.5 * al * (alpha / M_PI) * z / r;

  return pre_factor * int_result;
}

//------------------------------------------------------------------------------
double vSEl(double r, double rN, double z, double alpha) {
  auto sp1 = SafeProfiler::profile(__func__);
  //
  auto l = 0; // XXX
  auto bl = fit_AB.Bl(l, z * alpha);

  auto pre_factor =
      -1.5 * bl * z * z * alpha * alpha * alpha / (rN * rN * rN * r);

  auto f = [=](double x) {
    auto arg_neg = z * std::abs(r - x);
    auto arg_pos = z * (r + x);
    auto a = (arg_neg + 1.0) * std::exp(-arg_neg);
    auto b = (arg_pos + 1.0) * std::exp(-arg_pos);
    return x * (a - b);
  };

  return pre_factor * NumCalc::num_integrate(f, 0.0, rN, 1000, NumCalc::linear);

  // This is the "point-nucleus" version. Gives same result!
  // auto pre_factor2 = -bl * z * z * z * z * alpha * alpha * alpha;
  // return pre_factor2 * std::exp(-z * r);
}

//******************************************************************************
double vSE_Hmag(double r, double rN, double z, double alpha) {
  auto sp1 = SafeProfiler::profile(__func__);

  static constexpr double abs_err_lim = 0.0;
  static constexpr double rel_err_lim = 1.0e-6;
  static constexpr unsigned long max_num_subintvls = 750; //?

  gsl_set_error_handler_off(); //?
  if (rN <= 0) {
    rN = 1.0e-7;
  }

  RadPot_params params = {r, rN, z, alpha};

  gsl_function f_gsl;
  f_gsl.function = &gslfunc_SEmag;
  f_gsl.params = &params;

  gsl_integration_workspace *gsl_int_wrk =
      gsl_integration_workspace_alloc(max_num_subintvls + 1);

  double int_result, abs_err;
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  const auto pre_factor = 3.0 * z * alpha / M_PI; // XXX check!
  return pre_factor * int_result;
}

//******************************************************************************
//******************************************************************************

namespace Helper {

//------------------------------------------------------------------------------
double vUehcommon(double t, double chi) {
  auto t2 = t * t;
  auto chi3 = chi * chi * chi;
  auto a = std::sqrt(t2 - 1.0);
  auto b = (1.0 / t2 + 1.0 / (2.0 * t2 * t2));
  auto c = 2.0 / chi3;
  return a * b * c;
}

//------------------------------------------------------------------------------
double vUehf_smallr(double r, double rN, double chi) {
  // nb: a part can be removed and done analyitcally
  auto a = (r / rN) * chi;
  auto b = std::exp(-chi) * (1.0 + chi);
  auto c = std::sinh(a);
  return a - b * c;
}
double vUehf_larger(double r, double rN, double chi) {
  const auto tmp = (r / rN) * chi;
  const auto a = std::exp(-tmp);
  // const auto b = chi * std::cosh(chi) - sinh(chi); // XXX unstable small chi
  // Good to parts in 10^6
  const auto chi3 = chi * chi * chi;
  const auto b = (chi < 0.5) ? chi3 / 3.0 + (chi3 * chi * chi) / 30.0 +
                                   (chi3 * chi3 * chi) / 840.0
                             : chi * std::cosh(chi) - sinh(chi);
  return a * b;
}

//------------------------------------------------------------------------------
double gslfunc_Ueh_smallr(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto g = vUehcommon(t, chi);
  const auto f = vUehf_smallr(r, rN, chi);
  return f * g;
}
double gslfunc_Ueh_larger(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto g = vUehcommon(t, chi);
  const auto f = vUehf_larger(r, rN, chi);
  return f * g;
}

//******************************************************************************
double gb_I1(double t, double z, double alpha) {
  auto t2 = t * t;
  auto za = z * alpha;
  double a = 1.0 / std::sqrt(t2 - 1.0);
  double b = 1.0 - 1.0 / (2.0 * t2);
  double c = std::log(t2 - 1.0) + 4.0 * std::log(1.0 / za + 0.5);
  double d = -1.5 + 1.0 / t2;
  return a * (b * c + d);
}

double gb_I2(double t, double r, double rN, double z, double alpha) {
  auto za = z * alpha;
  auto ttoa = 2.0 * t / alpha;
  double rA = 0.07 * za * za * alpha;
  auto expr = std::exp(rA * ttoa);

  auto f = [=](double x) {
    auto arg1 = (std::abs(r - x) + rA) * ttoa;
    auto arg2 = (r + x + rA) * ttoa;
    auto a = ExpInt1(arg1) - ExpInt1(arg2);
    return x * a;
  };

  return expr * NumCalc::num_integrate(f, 0.0, rN, 1000, NumCalc::linear);
}

//------------------------------------------------------------------------------
double gb_GSEh_larger(double r, double rN, double chi) {
  // Essentially duplicate of Uehling, but swapped large/small r ?
  const auto chi3 = chi * chi * chi;
  auto a = std::exp(-chi * r / rN) * 2.0 / chi3;
  const auto b = (chi < 0.5) ? chi3 / 3.0 + (chi3 * chi * chi) / 30.0 +
                                   (chi3 * chi3 * chi) / 840.0
                             : chi * std::cosh(chi) - sinh(chi);
  return a * b;
}

double gb_GSEh_smallr(double r, double rN, double chi) {
  // Essentially duplicate of Uehling, but swapped large/small r ?
  const auto chi3 = chi * chi * chi;
  auto a = 2.0 / chi3;
  auto b = (chi * r / rN);
  auto c = std::exp(-chi) * (1.0 + chi) * std::sinh(chi * r / rN);
  return a * (b - c);
}

//------------------------------------------------------------------------------
double gslfunc_SEh_smallr(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  const auto chi = 2.0 * t * rN / alpha;
  double rA = 0.07 * z * z * alpha * alpha * alpha;
  return gb_I1(t, z, alpha) *
         (gb_GSEh_smallr(r, rN, chi) -
          (rA / (rN * rN * rN)) * gb_I2(t, r, rN, z, alpha));
}
double gslfunc_SEh_larger(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  const auto chi = 2.0 * t * rN / alpha;
  double rA = 0.07 * z * z * alpha * alpha * alpha;
  return gb_I1(t, z, alpha) *
         (gb_GSEh_larger(r, rN, chi) -
          (rA / (rN * rN * rN)) * gb_I2(t, r, rN, z, alpha));
}

//******************************************************************************
double seJmag(double x, double y) {
  const auto y3 = y * y * y;
  const auto sihncosh =
      (y < 0.5) ? -y3 / 3.0 - (y3 * y * y) / 30.0 - (y3 * y3 * y) / 840.0
                : sinh(y) - y * std::cosh(y);
  return std::exp(-x) * (1.0 + x) * sihncosh + y3 / 3.0;
}

double gslfunc_SEmag(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto eta = 2.0 * t * r / alpha;
  const auto factor =
      1.0 / (chi * chi * chi * eta * eta * std::sqrt(t * t - 1.0));
  return (r <= rN) ? factor * seJmag(chi, eta) : factor * seJmag(eta, chi);
}

} // namespace Helper
} // namespace RadiativePotential
