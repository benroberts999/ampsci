#pragma once
#include <cassert>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_version.h>
#include <iostream>
#include <map>
#include <vector>

// This should make code work with old versions of GSL
#ifdef GSL_MAJOR_VERSION
#if GSL_MAJOR_VERSION == 1
#define GSL_VERSION_1
#endif
#endif

//! Interpolates functions using cubic splines. Uses GSL:
//! https://www.gnu.org/software/gsl/doc/html/interp.html
namespace Interpolator {

//! Method (type) of 1D Interpolation used
/*! @details 
The following interpolation types are provided by GSL
(see https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_type)

linear

 * Linear interpolation. This interpolation method does not require any additional memory.

polynomial

 * Polynomial interpolation. This method should only be used for interpolating small numbers of points because polynomial interpolation introduces large oscillations, even for well-behaved datasets. The number of terms in the interpolating polynomial is equal to the number of points.

cspline

 * Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second derivatives at the supplied data-points. The second derivative is chosen to be zero at the first point and last point.

cspline_periodic

 * Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic on each interval, with matching first and second derivatives at the supplied data-points. The derivatives at the first and last points are also matched. Note that the last point in the data must have the same y-value as the first point, otherwise the resulting periodic interpolation will have a discontinuity at the boundary.

akima

 * Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.

akima_periodic

 * Non-rounded Akima spline with periodic boundary conditions. This method uses the non-rounded corner algorithm of Wodicka.

steffen

 * Steffen’s method guarantees the monotonicity of the interpolating function between the given data points. Therefore, minima and maxima can only occur exactly at the data points, and there can never be spurious oscillations between data points. The interpolated function is piecewise cubic in each interval. The resulting curve and its first derivative are guaranteed to be continuous, but the second derivative may be discontinuous.

*/
enum class Method {
  //! linear interpolation
  linear,
  //! polynomial interpolation
  polynomial,
  //! cubic b-spline interpolation
  cspline,
  //! cubic b-spline interpolation with periodic boundary condition
  cspline_periodic,
  //! akima interpolation
  akima,
  //! akima interpolation with periodic boundary condition
  akima_periodic,
  //! steffen interpolation (ensure monotonicity between points). Only GSLv2+
  steffen
};

//! Maps from enum 'Interpolator::Method' to gsl_interp_type
static const std::map<Method, const gsl_interp_type *> interp_method{
    {Method::linear, gsl_interp_linear},
    {Method::polynomial, gsl_interp_polynomial},
    {Method::cspline, gsl_interp_cspline},
    {Method::cspline_periodic, gsl_interp_cspline_periodic},
    {Method::akima, gsl_interp_akima},
    {Method::akima_periodic, gsl_interp_akima_periodic}
#ifndef GSL_VERSION_1
    ,
    {Method::steffen, gsl_interp_steffen}
#endif
};

//==============================================================================

//! Performs interpolation using GSL (GNU Scientific Library)
/*! @details 
  On construction takes in vectors x and y, to be interpolated.
  These are intepreted as function values y(x) at descrete points x.
  x and y dimensions must match.
  Given new x', will return y(x') based on interpolation.
  NOTE: Interpolates, but does NOT extrapolate! Everything outside the region (xmin,xmax) from initial input will be zero.
  It may use any interpolation method provided by GSL library.
  See Interpolator::Method for list of options.
  By default, cspline (cubic b-spline) method is used.
  Note: gsl function takes *pointers* to x and y data. Therefore, vectors x and y should outlive the Interp object
*/
class Interp {

private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double x0, xmax;

public:
  //! x_in and y_in
  Interp(const std::vector<double> &x, const std::vector<double> &y,
         Method method = Method::cspline)
      : acc(gsl_interp_accel_alloc()),
        spline(gsl_spline_alloc(interp_method.at(method), x.size())),
        x0(x.front()),
        xmax(x.back()) {
    assert(x.size() == y.size() &&
           "In Interp, input vectors x and y must have same size");
    assert(x.size() >= gsl_interp_type_min_size(interp_method.at(method)) &&
           "In Interp, certain interpolation methods require a minimum number "
           "of points");
    gsl_spline_init(spline, x.data(), y.data(), x.size());
  }
  ~Interp() {
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  Interp(const Interpolator::Interp &) = delete;
  Interp &operator=(const Interpolator::Interp &) = delete;

  //! Evaluates interpolation function at point x. Does not extrapolate.
  double interp(double x) const {
    if (x < x0 || x > xmax)
      return 0.0;
    else
      return gsl_spline_eval(spline, x, acc);
  }

  //! Evaluates interpolation function at set of points {x}. Does not extrapolate.
  std::vector<double> interp(const std::vector<double> &x) const {
    std::vector<double> y;
    y.reserve(x.size());
    for (auto xi : x)
      y.push_back(interp(xi));
    return y;
  }

  //! Evaluates interpolation function at point x. Does not extrapolate.
  double operator()(double x) const { return interp(x); }

  //! Evaluates interpolation function at set of points {x}. Does not extrapolate.
  std::vector<double> operator()(const std::vector<double> &x) const {
    return interp(x);
  }
};

//==============================================================================

//! Performs interpolation using GSL (GNU Scientific Library)
/*! @details Takes set of points {xin, yin}, interpolates and evaluates new y
  values at positions defined by {x_out}; returns as vector.
  Just a wrapper for class Interp
*/
inline std::vector<double> interpolate(const std::vector<double> &x_in,
                                       const std::vector<double> &y_in,
                                       const std::vector<double> &x_out,
                                       Method method = Method::cspline) {
  Interp i_func(x_in, y_in, method);
  return i_func(x_out);
}

//==============================================================================

//! Check if steffen method available (only with GSL version 2+)
static constexpr bool has_steffen_method() {
#ifndef GSL_VERSION_1
  return true;
#else
  return false;
#endif
}

} // namespace Interpolator
