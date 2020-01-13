#include "RadiativePotential.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <cmath>
#include <gsl/gsl_errno.h> //?
#include <gsl/gsl_integration.h>
#include <iostream>
// E1(x) = -std::expint(-x)

namespace RadiativePotential {

double vEuhcommon(double t, double chi) {
  auto t2 = t * t;
  auto chi3 = chi * chi * chi;
  auto a = std::sqrt(t2 - 1.0);
  auto b = (1.0 / t2 + 1.0 / (2.0 * t2 * t2));
  auto c = 2.0 / chi3;
  return a * b * c;
}

double vEuhf_smallr(double r, double rN, double chi) {
  // nb: a part can be removed and done analyitcally
  auto a = (r / rN) * chi;
  auto b = std::exp(-chi) * (1.0 + chi);
  auto c = std::sinh(a);
  return a - b * c;
}
double vEuhf_larger(double r, double rN, double chi) {
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

double gslfunc_Ueh_smallr(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto g = vEuhcommon(t, chi);
  const auto f = vEuhf_smallr(r, rN, chi);
  return f * g;
}
double gslfunc_Ueh_larger(double t, void *p) {
  const auto [r, rN, z, alpha] = *(static_cast<RadPot_params *>(p));
  (void)z;
  const auto chi = 2.0 * t * rN / alpha;
  const auto g = vEuhcommon(t, chi);
  const auto f = vEuhf_larger(r, rN, chi);
  return f * g;
}

double vEuhling(double r, double rN, double z, double alpha) {

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
      gsl_integration_workspace_alloc(max_num_subintvls + 1); // -1 ?
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
    // std::cout << -std::expint(-arg1) << " " << std::expint(-arg2) << "\n";
    return x * a;
  };

  // std::cout << "f(0)=" << f(0.0) << " " << f(1.0e-7) << " " << f(rN) << "\n";

  return expr * NumCalc::num_integrate(f, 0.0, rN, 1000, NumCalc::linear);
  // return expr *
  //        NumCalc::num_integrate(f, 1.0e-7, rN, 5000, NumCalc::logarithmic);
}

// double gb_GSEh_smallr(double r, double rN, double chi) {
double gb_GSEh_larger(double r, double rN, double chi) {
  // Essentially duplicate of Euhling, but swapped large/small r ?
  const auto chi3 = chi * chi * chi;
  auto a = std::exp(-chi * r / rN) * 2.0 / chi3;
  const auto b = (chi < 0.5) ? chi3 / 3.0 + (chi3 * chi * chi) / 30.0 +
                                   (chi3 * chi3 * chi) / 840.0
                             : chi * std::cosh(chi) - sinh(chi);
  return a * b;
}

// double gb_GSEh_larger(double r, double rN, double chi) {
double gb_GSEh_smallr(double r, double rN, double chi) {
  // Essentially duplicate of Euhling, but swapped large/small r ?
  const auto chi3 = chi * chi * chi;
  auto a = 2.0 / chi3;
  auto b = (chi * r / rN);
  auto c = std::exp(-chi) * (1.0 + chi) * std::sinh(chi * r / rN);
  return a * (b - c);
}

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

static Fit_AB fit_AB;

double vSEh(double r, double rN, double z, double alpha) {

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
      gsl_integration_workspace_alloc(max_num_subintvls + 1); // -1 ?

  double int_result, abs_err;
  gsl_integration_qagiu(&f_gsl, 1.0, abs_err_lim, rel_err_lim,
                        max_num_subintvls, gsl_int_wrk, &int_result, &abs_err);
  gsl_integration_workspace_free(gsl_int_wrk);

  auto l = 0; // XXX
  auto al = fit_AB.Al(l, z * alpha);

  auto pre_factor = -1.5 * al * (alpha / M_PI) * z / r;

  // std::cout << r << " " << pre_factor * int_result
  //           << ", d=" << abs_err / int_result << "\n";

  return pre_factor * int_result;
}

} // namespace RadiativePotential
