#pragma once

namespace RadiativePotential {

// Routines return the first approximation which has an absolute error
// smaller than abs_err_lim or a relative error smaller than rel_err_lim.
// Note that this is an either-or constraint, not simultaneous. To compute to
// a specified absolute error, set epsrel to zero (etc.)
static constexpr double abs_err_lim = 1.0e-3;
static constexpr double rel_err_lim = 1.0e-6;
// max_num_subintvls < size(gsl_int_wrk)
static constexpr unsigned long max_num_subintvls = 500; //?

//*******************************
struct RadPot_params {
  double r, rN, z, alpha;
};

double vEuhcommon(double t, double chi);
double vEuhf_smallr(double r, double rN, double chi);
double vEuhf_larger(double r, double rN, double chi);
double gslfunc_Ueh_smallr(double t, void *p);
double gslfunc_Ueh_larger(double t, void *p);

double vEuhling(double r, double rN, double z, double alpha);

} // namespace RadiativePotential
