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

/*!
  @brief 1D interpolation using GSL splines.
  @details
  Wraps the GSL interpolation library
  (https://www.gnu.org/software/gsl/doc/html/interp.html).
  Provides a stateful `Interp` class and a convenience `interpolate()`
  function. Several interpolation methods are available via `Method`.
  Does not extrapolate: values outside [xmin, xmax] are returned as zero.
*/
namespace Interpolator {

/*!
  @brief Interpolation method.
  @details
  See https://www.gnu.org/software/gsl/doc/html/interp.html for full details.
*/
enum class Method {
  //! Linear interpolation; no additional memory required
  linear,
  //! Polynomial interpolation; only suitable for small numbers of points
  polynomial,
  //! Cubic spline with natural boundary conditions (zero second derivative at endpoints)
  cspline,
  //! Cubic spline with periodic boundary conditions; first and last y-values must match
  cspline_periodic,
  //! Akima spline with natural boundary conditions (Wodicka non-rounded corner algorithm)
  akima,
  //! Akima spline with periodic boundary conditions
  akima_periodic,
  //! Steffen's monotone spline: no spurious oscillations between data points. GSL v2+ only.
  steffen
};

//! Maps Interpolator::Method to the corresponding GSL interpolation type
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

/*!
  @brief Stateful 1D interpolation object using GSL.
  @details
  Constructed from data vectors x and y (interpreted as samples of y(x)).
  Once constructed, evaluates the interpolated function at arbitrary points
  via `interp()` or `operator()`.

  Does not extrapolate: returns 0.0 for x outside [x.front(), x.back()].

  Copy construction and copy assignment are deleted.

  @note x and y are copied into GSL internal storage during construction;
  the input vectors do not need to outlive the `Interp` object.
*/
class Interp {

private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double x0, xmax;

public:
  /*!
    @brief Construct interpolation object from data points.
    @param x      Strictly increasing x values (abscissae).
    @param y      Function values y(x); must satisfy `y.size() == x.size()`.
    @param method Interpolation method (default: cspline).
  */
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

  //! Copy construction deleted
  Interp(const Interpolator::Interp &) = delete;
  //! Copy assignment deleted
  Interp &operator=(const Interpolator::Interp &) = delete;

  //! Returns interpolated y(x). Returns 0.0 if x is outside [x0, xmax].
  double interp(double x) const {
    if (x < x0 || x > xmax)
      return 0.0;
    else
      return gsl_spline_eval(spline, x, acc);
  }

  //! Returns interpolated y(x) for each point in x. Returns 0.0 outside range.
  std::vector<double> interp(const std::vector<double> &x) const {
    std::vector<double> y;
    y.reserve(x.size());
    for (auto xi : x)
      y.push_back(interp(xi));
    return y;
  }

  //! Returns interpolated y(x). Returns 0.0 if x is outside [x0, xmax].
  double operator()(double x) const { return interp(x); }

  //! Returns interpolated y(x) for each point in x. Returns 0.0 outside range.
  std::vector<double> operator()(const std::vector<double> &x) const {
    return interp(x);
  }
};

//==============================================================================

/*!
  @brief Convenience wrapper: interpolates y_in(x_in) and evaluates at x_out.
  @param x_in   Input x values (strictly increasing).
  @param y_in   Input y values; must satisfy `y_in.size() == x_in.size()`.
  @param x_out  Points at which to evaluate the interpolant.
  @param method Interpolation method (default: cspline).
  @returns Vector of interpolated values at each point in x_out.
*/
inline std::vector<double> interpolate(const std::vector<double> &x_in,
                                       const std::vector<double> &y_in,
                                       const std::vector<double> &x_out,
                                       Method method = Method::cspline) {
  Interp i_func(x_in, y_in, method);
  return i_func(x_out);
}

//==============================================================================

//! Returns true if the Steffen interpolation method is available (GSL v2+)
static constexpr bool has_steffen_method() {
#ifndef GSL_VERSION_1
  return true;
#else
  return false;
#endif
}

} // namespace Interpolator
