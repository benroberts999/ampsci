#pragma once
#include "LinAlg/LinAlg.hpp"
#include "qip/Vector.hpp"
#include <gsl/gsl_bspline.h>
#include <iostream>
#include <numeric>
#include <vector>

/*!
Calculates basis of N B-splines of order K (degree K-1), defined over range [0,
xmax]. x0 is first non-zero ("internal") knot point.
@details
 * N - K + 2 break-points "internal knots" placed between [x0,xmax].
   * By default: on a _logarithmic_ scale
   * May also be placed on a linear or log-linear scale
   * If log-linear scale, b=(xmax/2) (see qip::loglinear_range)
 * knots are the same as breakpoints, but the ends are repeated K times.
 * Uses the GSL Bspline library, gsl/gsl_bspline.h
   * https://www.gnu.org/software/gsl/doc/html/bspline.html
 * See also: Bachau et al.,  Reports Prog. Phys. 64, 1815 (2001).
*/
class BSpline {

  // K is _order_ (degree is K-1)
  std::size_t m_K;
  // Total number of splines (i=0,1,...,N-1)
  std::size_t m_N;
  // Splines defined over range [0,x_max]
  double m_xmax;
  // Worskspace used by GSL library
  gsl_bspline_workspace *gsl_bspl_work{nullptr};

public:
  enum class KnotDistro { logarithmic, linear, loglinear };
  //****************************************************************************
  //! Constructs basis of n splines of order k, defined over [0,xmax]. x0 is
  //! first non-zero knot. kd = { logarithmic, linear, loglinear }
  BSpline(std::size_t n, std::size_t k, double x0, double xmax,
          KnotDistro kd = KnotDistro::logarithmic)
      : m_K(k), m_N(n), m_xmax(xmax) {
    assert(m_N > m_K && "Require N>K for B-splines");
    assert(xmax > x0 && x0 > 0.0 && xmax > 0.0 && "xmax>x0 and both positive");

    // Allocate memory for the GSL workspace:
    const std::size_t n_break = m_N - m_K + 2;
    gsl_bspl_work = gsl_bspline_alloc(m_K, n_break);

    set_knots(x0, xmax, n_break, kd);

    // Consistancy check;
    assert(gsl_bspline_ncoeffs(gsl_bspl_work) == m_N);
  }

  // delete copy assignment (no simple way to copy gsl_bspl_work)
  BSpline &operator=(const BSpline &) = delete;
  // delete copy constructor (no simple way to copy gsl_bspl_work)
  BSpline(const BSpline &) = delete;

  // Desctructor
  ~BSpline() { gsl_bspline_free(gsl_bspl_work); }

  //****************************************************************************
  //! Order of the splines, K
  std::size_t K() { return m_K; }
  //! Degree of the splines, d:=K-1
  std::size_t d() { return m_K - 1; }
  //! Number of splines, N
  std::size_t N() { return m_N; }

  //****************************************************************************
  //! Returns the first nonzero spline index i0; Note that i runs [0,N).
  //! The last non-zero index is min(i0+K-1, N-1).
  std::size_t find_i0(double x) {
    for (std::size_t i = 0; i < m_N; ++i) {
      if (gsl_vector_get(gsl_bspl_work->knots, i + m_K) > x) {
        return i;
      }
    }
    return m_N;
  }

  //****************************************************************************
  //! Returns a std::vector of all splines {b0, b1, ..., b_N-1} evaluated at x
  std::vector<double> get(double x) {
    std::vector<double> b(m_N);
    gsl_vector_view b_gsl = gsl_vector_view_array(b.data(), b.size());
    if (x >= 0.0 && x <= m_xmax) {
      gsl_bspline_eval(x, &b_gsl.vector, gsl_bspl_work);
    }
    return b;
  }

  //****************************************************************************
  //! Returns a pair {i0, M}, where i0 is spline index of first non-zero spline,
  //! and M is a matrix of non-zero splines and their derivatives. M_ij contains
  //! the jth derivative of the spline b[i+i0] evaulated at x.
  /*! @details
    - n_deriv is the maximum derivative calculated (inclusive).
    - Matrix M is has K rows, and (n_deriv+1) columns.
    - If x>xmax, will return a matrix of zeros, and i0+K will be >N
    - If x=xmax, only last spline is nonzero, still i0+K=N [lower derivs may be
    non-zero]
  */
  std::pair<std::size_t, LinAlg::Matrix<double>>
  get_nonzero(double x, std::size_t n_deriv) {
    std::pair<std::size_t, LinAlg::Matrix<double>> out{0, {m_K, n_deriv + 1}};
    auto &i0 = out.first;
    auto &bij = out.second;

    i0 = find_i0(x);

    if (x >= 0.0 && x <= m_xmax) {
      // outside this range, spline set to zero by default. Invalid to call
      // gsl_bspline_deriv_eval_nonzero() outside spline range
      auto i_end = i0 + m_K - 1;
      gsl_matrix_view b_gsl = bij.as_gsl_view();
      gsl_bspline_deriv_eval_nonzero(x, n_deriv, &b_gsl.matrix, &i0, &i_end,
                                     gsl_bspl_work);
    }
    return out;
  }

  //****************************************************************************
  void print_knots() const {

    std::size_t n_cols_to_print = 4;
    std::cout << "Knots:\n ";
    std::size_t count = 0;
    while (count < m_N + m_K) {
      auto t = gsl_bspl_work->knots->data[count];
      printf("%.5e, ", t);
      ++count;
      if (count % n_cols_to_print == 0)
        std::cout << "\n ";
    }
    if (count % n_cols_to_print != 0)
      std::cout << "\n";
  }

private:
  //****************************************************************************
  void set_knots(double x0, double xmax, std::size_t n_break, KnotDistro kd) {
    // Make version that takes custom (internal) knot sequence?
    // Option for linear/log-spaced knots?

    // n_break break points range from [0,xmax] in n_break steps
    // We define the _internal_ (non-zero) breakpoints
    // these run from [x0, xmax] in n_break - 1 steps
    // Knots are the same as breakpoints, but the first (0) and last (xmax)
    // points are repeated K times
    auto breaks = (kd == KnotDistro::logarithmic) ?
                      qip::logarithmic_range(x0, xmax, n_break - 1) :
                      (kd == KnotDistro::linear) ?
                      qip::uniform_range(x0, xmax, n_break - 1) :
                      (kd == KnotDistro::loglinear) ?
                      qip::loglinear_range(x0, xmax, 0.5 * xmax, n_break - 1) :
                      std::vector<double>{};

    breaks.insert(breaks.begin(), 0.0);

    auto gsl_break_vec = gsl_vector_view_array(breaks.data(), breaks.size());
    gsl_bspline_knots(&gsl_break_vec.vector, gsl_bspl_work);
  }
};
