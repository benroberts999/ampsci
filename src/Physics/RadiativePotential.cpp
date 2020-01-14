#include "RadiativePotential.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <cmath>
#include <gsl/gsl_errno.h> //?
#include <gsl/gsl_integration.h>
#include <iostream>
// E1(x) = -std::expint(-x)

namespace RadiativePotential {

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

//------------------------------------------------------------------------------
double vUehling(double r, double rN, double z, double alpha) {
  auto sp1 = SafeProfiler::profile(__func__);

  // Routines return the first approximation which has an absolute error
  // smaller than abs_err_lim or a relative error smaller than rel_err_lim.
  // Note that this is an either-or constraint, not simultaneous. To compute to
  // a specified absolute error, set epsrel to zero (etc.)
  static constexpr double abs_err_lim = 0;
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
  // E1(x) = -std::expint(-x)
  auto za = z * alpha;
  auto ttoa = 2.0 * t / alpha;
  double rA = 0.07 * za * za * alpha;
  auto expr = std::exp(rA * ttoa);

  auto f = [=](double x) {
    auto arg1 = (std::abs(r - x) + rA) * ttoa;
    auto arg2 = (r + x + rA) * ttoa;
    auto a = -std::expint(-arg1) + std::expint(-arg2);
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

//------------------------------------------------------------------------------
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

  const auto pre_factor = 3.0 * z / M_PI;
  return pre_factor * int_result;
}

} // namespace RadiativePotential
