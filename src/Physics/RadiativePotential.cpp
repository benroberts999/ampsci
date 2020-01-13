#include "RadiativePotential.hpp"
#include <cmath>
#include <gsl/gsl_errno.h> //?
#include <gsl/gsl_integration.h>
#include <iostream>

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
  /// XXX add catch?
  gsl_integration_workspace_free(gsl_int_wrk);

  auto pre_factor = z * (alpha / M_PI) / r;

  return pre_factor * int_result;
}

} // namespace RadiativePotential
