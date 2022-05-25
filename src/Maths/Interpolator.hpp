#pragma once
#include "IO/SafeProfiler.hpp"
#include <cassert>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <vector>

// https://www.gnu.org/software/gsl/doc/html/interp.html

//! Interpolates functions using cubic splines. Uses GSL
namespace Interpolator {

//! @brief Interpolates {xin,yin} onto {xout}
//! @details Takes set of points {xin, yin}, interpolates and evaluates new y
//! values at positions defined by {x_out}; returns as vector.
//!  - NOTE: Interpolates, but does NOT extrapolate! Everything outseide the
//!  region (xmin,xmax)_in will be zero
inline std::vector<double> interpolate(const std::vector<double> &x_in,
                                       const std::vector<double> &y_in,
                                       const std::vector<double> &x_out) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  gsl_set_error_handler_off(); //?

  const auto dim_in = x_in.size();
  if (dim_in != y_in.size() || dim_in == 0) {
    std::cerr << "FAIL 11 in Interpolate: xin,yin different sizes: "
              << x_in.size() << " " << y_in.size() << "\n";
    return {};
  }

  // Do Interpolation using cubic splines
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, dim_in);

  gsl_spline_init(spline, &x_in[0], &y_in[0], dim_in);

  auto xmin = x_in.front();
  auto xmax = x_in.back();

  std::vector<double> y_out;
  y_out.reserve(x_out.size());
  for (const auto &x : x_out) {
    if (x < xmin || x > xmax)
      y_out.push_back(0.0);
    else
      y_out.push_back(gsl_spline_eval(spline, x, acc));
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return y_out;
}

//==============================================================================
class Interp {

  gsl_interp_accel *acc; // = gsl_interp_accel_alloc();
  gsl_spline *spline;    // = gsl_spline_alloc(gsl_interp_cspline, dim_in);
  double x0, xmax;

public:
  Interp(const std::vector<double> &x_in, const std::vector<double> &y_in)
      : acc(gsl_interp_accel_alloc()),
        spline(gsl_spline_alloc(gsl_interp_cspline, x_in.size())),
        x0(x_in.front()),
        xmax(x_in.back()) {
    assert(x_in.size() == y_in.size());
    gsl_spline_init(spline, x_in.data(), y_in.data(), x_in.size());
  }

  ~Interp() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  Interp(const Interpolator::Interp &) = delete;
  Interp &operator=(const Interpolator::Interp &) = delete;

  double interp(double x) {
    if (x < x0 || x > xmax)
      return 0.0;
    else
      return gsl_spline_eval(spline, x, acc);
  }

  double operator()(double x) { return interp(x); }
};

} // namespace Interpolator
