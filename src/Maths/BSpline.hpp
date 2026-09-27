#pragma once
#include "LinAlg/include.hpp"
#include "qip/Vector.hpp"
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_version.h>
#include <iostream>
#include <numeric>
#include <vector>

// This should make code work with old versions of GSL
// See below for documentation of old GSL function
#ifdef GSL_MAJOR_VERSION
#if GSL_MAJOR_VERSION == 1
#define GSL_VERSION_1
#endif
// GSL 2.8 introduces substantial changes to bspline library.
#if GSL_MAJOR_VERSION <= 2 && GSL_MINOR_VERSION < 8
#define GSL_VERSION_PRIOR_2_8
#endif
#endif

/*!
  @brief Basis of N B-splines of order K (degree K-1), defined over [0, xmax].

  @details
  Constructs and evaluates a B-spline basis using the GSL B-spline library
  (gsl/gsl_bspline.h). See also: Bachau et al., Rep. Prog. Phys. 64, 1815
  (2001).

  **Knot structure:**
  - N - K + 2 breakpoints (internal knots) are placed in [x0, xmax], where x0
    is the first non-zero knot. The point 0 is prepended, giving N - K + 2
    breakpoints total spanning [0, xmax].
  - The full knot vector repeats the endpoints K times.
  - Breakpoints may be placed on a logarithmic (default), linear, or
    log-linear scale; see `KnotDistro`.

  **Compatibility:**
  Supports GSL v2.x (<2.8), and v2.8+. API differences are handled
  internally via compile-time version macros.
  Should also support v1.x, but this is difficult to test.

  Copy construction and copy assignment are deleted (GSL workspace is not
  copyable).
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

#ifdef GSL_VERSION_1
  // Worskspace used by old version of GSL library
  gsl_bspline_deriv_workspace *gsl_bspl_deriv_work{nullptr};
#endif

public:
  //! Distribution of internal knot points
  enum class KnotDistro {
    //! Knots on a logarithmic scale (default; suits bound-state grids)
    logarithmic,
    //! Knots evenly spaced
    linear,
    //! Log-linear blend; midpoint at xmax/2
    loglinear
  };

  //============================================================================
  /*!
    @brief Construct a basis of N B-splines of order K over [0, xmax].
    @param n   Number of splines N. Must satisfy N > K.
    @param k   Order K (degree K-1). E.g. K=4 gives cubic splines.
    @param x0  First non-zero (internal) knot; must satisfy 0 < x0 < xmax.
    @param xmax Upper bound of the spline domain.
    @param kd  Knot distribution: logarithmic (default), linear, or loglinear.
  */
  BSpline(std::size_t n, std::size_t k, double x0, double xmax,
          KnotDistro kd = KnotDistro::logarithmic)
    : m_K(k), m_N(n), m_xmax(xmax) {
    assert(m_N > m_K && "Require N>K for B-splines");
    assert(xmax > x0 && x0 > 0.0 && xmax > 0.0 && "xmax>x0 and both positive");

    // Allocate memory for the GSL workspace:
    const std::size_t n_break = m_N - m_K + 2;
    gsl_bspl_work = gsl_bspline_alloc(m_K, n_break);

#ifdef GSL_VERSION_1
    // Worskspace used by old version of GSL library
    gsl_bspl_deriv_work = gsl_bspline_deriv_alloc(m_K);
#endif

    set_knots(x0, xmax, n_break, kd);

    // Consistancy check:
#ifdef GSL_VERSION_PRIOR_2_8
    // gsl_bspline_ncoeffs will be deprecated in future GSL
    assert(gsl_bspline_ncoeffs(gsl_bspl_work) == m_N);
#else
    assert(gsl_bspline_ncontrol(gsl_bspl_work) == m_N);
#endif
    assert(gsl_bspline_order(gsl_bspl_work) == m_K);
    assert(gsl_bspline_nbreak(gsl_bspl_work) == n_break);
  }

  //! Copy assignment deleted: GSL workspace is not copyable
  BSpline &operator=(const BSpline &) = delete;
  //! Copy construction deleted: GSL workspace is not copyable
  BSpline(const BSpline &) = delete;

  ~BSpline() {
    gsl_bspline_free(gsl_bspl_work);
#ifdef GSL_VERSION_1
    gsl_bspline_deriv_free(gsl_bspl_deriv_work);
#endif
  }

  //============================================================================
  //! Order of the splines, K
  std::size_t K() { return m_K; }
  //! Degree of the splines, d := K-1
  std::size_t d() { return m_K - 1; }
  //! Number of splines, N
  std::size_t N() { return m_N; }

  //============================================================================
  //! Returns the index i0 of the first non-zero spline at x.
  //! Spline indices run [0, N); at most K splines are non-zero at any x,
  //! occupying indices [i0, min(i0+K-1, N-1)].
  std::size_t find_i0(double x) {
    for (std::size_t i = 0; i < m_N; ++i) {
      if (gsl_vector_get(gsl_bspl_work->knots, i + m_K) > x) {
        return i;
      }
    }
    return m_N;
  }

  //============================================================================
  //! Returns all N spline values {b_0(x), b_1(x), ..., b_{N-1}(x)}.
  //! Returns a vector of zeros if x is outside [0, xmax].
  std::vector<double> get(double x) {
    std::vector<double> b(m_N);
    gsl_vector_view b_gsl = gsl_vector_view_array(b.data(), b.size());
    if (x >= 0.0 && x <= m_xmax) {

#ifdef GSL_VERSION_PRIOR_2_8
      // gsl_bspline_eval will be deprecated in future GSL
      gsl_bspline_eval(x, &b_gsl.vector, gsl_bspl_work);
#else
      gsl_bspline_eval_basis(x, &b_gsl.vector, gsl_bspl_work);
#endif
    }
    return b;
  }

  //============================================================================
  /*!
    @brief Returns the non-zero splines and their derivatives at x.
    @details
    Returns `{i0, M}`, where:
    - `i0` is the index of the first non-zero spline,
    - `M` is a K  (n_deriv+1) matrix; `M(i, j)` is the j-th derivative of
      spline `b[i0 + i]` evaluated at x.

    At most K splines are non-zero at any x. This is more efficient than
    `get()` when only the non-zero splines are needed.

    If x > xmax, returns a zero matrix with i0+K > N.
    If x = xmax, only the last spline is non-zero (lower derivatives may not
    be zero).

    @param x       Evaluation point.
    @param n_deriv Maximum derivative order to compute (inclusive).
    @returns Pair `{i0, M}`.
  */
  std::pair<std::size_t, LinAlg::Matrix<double>>
  get_nonzero(double x, std::size_t n_deriv) {
    std::pair<std::size_t, LinAlg::Matrix<double>> out{0, {m_K, n_deriv + 1}};

    if (x >= 0.0 && x <= m_xmax) {
      // outside this range, spline set to zero by default. Invalid to call
      // gsl_bspline_deriv_eval_nonzero() outside spline range

      auto &i0 = out.first;
      auto &bij = out.second;
      gsl_matrix_view b_gsl = bij.as_gsl_view();

#ifdef GSL_VERSION_1
      size_t i_end{};
      gsl_bspline_deriv_eval_nonzero(x, n_deriv, &b_gsl.matrix, &i0, &i_end,
                                     gsl_bspl_work, gsl_bspl_deriv_work);
#elif defined GSL_VERSION_PRIOR_2_8
      size_t i_end{};
      gsl_bspline_deriv_eval_nonzero(x, n_deriv, &b_gsl.matrix, &i0, &i_end,
                                     gsl_bspl_work);
#else
      gsl_bspline_basis_deriv(x, n_deriv, &b_gsl.matrix, &i0, gsl_bspl_work);
#endif
    }

    return out;
  }

  //============================================================================
  //! Prints the full knot vector to stdout (N+K values).
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
  //============================================================================
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

#ifdef GSL_VERSION_PRIOR_2_8
    // gsl_bspline_knots will be deprecated in future GSL
    gsl_bspline_knots(&gsl_break_vec.vector, gsl_bspl_work);
#else
    gsl_bspline_init_augment(&gsl_break_vec.vector, gsl_bspl_work);
#endif
  }
};

// Documentation from OLD GSL v:1.x
//  -- Function: gsl_bspline_deriv_workspace * gsl_bspline_deriv_alloc
//           (const size_t K)
//      This function allocates a workspace for computing the derivatives
//      of a B-spline basis function of order K.  The size of the workspace
//      is O(2k^2).
//
//  -- Function: int gsl_bspline_deriv_eval (const double X, const size_t
//           NDERIV, gsl_matrix * DB, gsl_bspline_workspace * W,
//           gsl_bspline_deriv_workspace * DW)
//      This function evaluates all B-spline basis function derivatives of
//      orders 0 through nderiv (inclusive) at the position X and stores
//      them in the matrix DB.  The (i,j)-th element of DB is
//      d^jB_i(x)/dx^j.  The matrix DB must be of size n = nbreak + k - 2
//      by nderiv + 1.  The value n may also be obtained by calling
//      `gsl_bspline_ncoeffs'.  Note that function evaluations are
//      included as the zeroth order derivatives in DB.  Computing all the
//      basis function derivatives at once is more efficient than
//      computing them individually, due to the nature of the defining
//      recurrence relation.
//
//  -- Function: int gsl_bspline_deriv_eval_nonzero (const double X, const
//           size_t NDERIV, gsl_matrix * DB, size_t * ISTART, size_t *
//           IEND, gsl_bspline_workspace * W, gsl_bspline_deriv_workspace
//           * DW)
//      This function evaluates all potentially nonzero B-spline basis
//      function derivatives of orders 0 through nderiv (inclusive) at the
//      position X and stores them in the matrix DB.  The (i,j)-th element
//      of DB is d^j/dx^j B_(istart+i)(x).  The last row of DB contains
//      d^j/dx^j B_(iend)(x).  The matrix DB must be of size k by at least
//      nderiv + 1.  Note that function evaluations are included as the
//      zeroth order derivatives in DB.  By returning only the nonzero
//      basis functions, this function allows quantities involving linear
//      combinations of the B_i(x) and their derivatives to be computed
//      without unnecessary terms.
