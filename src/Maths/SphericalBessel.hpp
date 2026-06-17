#pragma once
#include "LinAlg/Matrix.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <vector>

/*!
  @brief Spherical Bessel functions j_L(x) and y_L(x), and related utilities.
  @details
  Provides an approximate version (good to ~1e-9) and exact GSL-backed versions,
  plus cell-averaged vectors for integration on coarse grids, and a lookup table.
*/
namespace SphericalBessel {

//==============================================================================
//! Spherical Bessel function j_L(x); uses low-x approximation where appropriate
double jL(int L, double x);

//! Spherical Bessel function of first kind via GSL
double exactGSL_JL(int L, double x);

//! Spherical Bessel functions of second kind
double exactGSL_YL(int L, double x);

//==============================================================================

//! Spherical Bessel function of second kind, y_L(x)
double yL(int L, double x);

//! Phi_L(x) = [(2L+1)!! / x^L] * j_L(x). If tilde=true, returns 1 - Phi_L(x)
double PhiL(int L, double x, bool tilde = false);

//! Psi_L(x) = [x^{L+1} / (2L-1)!!] * y_L(x). If tilde=true, returns 1 - Psi_L(x)
double PsiL(int L, double x, bool tilde = false);

//==============================================================================
/*!
  @brief Fills a vector with the spherical Bessel function j_L sampled on the
  supplied grid; oscillation-resolving cell-average is used where needed.
  @details
  Returns a vector v such that v[i] is the value (or cell-average) of j_L on
  the grid point xvec[i]. The "cell" at index i is the interval bounded by the
  midpoints to its neighbours (asymmetric for the first/last point).

  This vector is intended for integration of the form
  \f$ \int j_L(x)\,\psi(x)\,dx \approx \sum_i w_i\,v_i\,\psi(x_i) \f$
  against a smooth wavefunction \f$\psi\f$ sampled on the same grid. When the
  grid cell at index i resolves the oscillation of j_L (\f$\Delta x_i \le
  \lambda_j / N\f$ with \f$\lambda_j = 2\pi\f$ asymptotically and N = 10), v[i]
  is the point value j_L(xvec[i]). When the cell is too coarse, v[i] is the
  cell average
  \f[
    v_i = \frac{1}{\Delta x_i}\int_{\text{cell}_i} j_L(x)\,dx ,
  \f]
  computed by trapezoidal sub-sampling with enough sub-points to resolve j_L
  on the cell. Using these cell averages in the Riemann-sum quadrature is
  equivalent to a Filon-style integration: the (smooth) factor \f$\psi\f$ is
  treated as piecewise-constant per cell while the oscillation of j_L is
  integrated exactly within each cell. This recovers the correct
  Riemann-Lebesgue (~1/x) suppression at large x and removes the aliasing seen
  with pure point-sampling on a coarse grid.

  @param l             Order of the spherical Bessel function.
  @param xvec          Grid (monotonic) on which j_L is sampled.
  @param cell_average  If true (default), enables the cell averaging described
                       above. If false, every entry is a pure point sample
                       j_L(xvec[i]); use this for inspection / plotting or to
                       reproduce legacy behaviour.
  @return Vector v with v[i] = j_L(xvec[i]) or its cell average; same size as xvec.

  @warning The returned values are NOT pure point samples on cells where
           averaging triggers. They are only meaningful in conjunction with a
           Riemann-sum-style quadrature using the same grid weights.
*/
std::vector<double> fillBesselVec(int l, const std::vector<double> &xvec,
                                  bool cell_average = true);

/*!
  @brief As fillBesselVec, with the argument of j_L being k*r instead of x.
  @details
  Returns a vector v with v[i] = j_L(k*rvec[i]) or, where the cell at index i
  is too coarse for the oscillation of j_L(kr) (cell width
  \f$\Delta r_i > 2\pi / (k\,N)\f$ with N = 10), the cell average
  \f[
    v_i = \frac{1}{\Delta r_i}\int_{\text{cell}_i} j_L(k r)\,dr .
  \f]
  See fillBesselVec for the usage and warning.

  @param l             Order of the spherical Bessel function.
  @param k             Momentum (or other) factor multiplying r in the argument.
  @param rvec          Radial grid on which j_L(kr) is sampled.
  @param cell_average  If true (default), enables cell averaging; if false,
                       point-samples j_L(k*rvec[i]) at every grid point.
  @return Vector v of same length as rvec.
*/
std::vector<double> fillBesselVec_kr(int l, double k,
                                     const std::vector<double> &rvec,
                                     bool cell_average = true);

/*!
  @brief In-place variant of fillBesselVec_kr; writes into a caller-owned vector.
  @details Identical semantics to the returning overload (including the
  cell-averaging behaviour); jl is resized to r.size().
  @param l             Order of the spherical Bessel function.
  @param k             Momentum factor multiplying r in the argument.
  @param r             Radial grid.
  @param jl            Output vector (resized internally); must be non-null.
  @param cell_average  If true (default), enables cell averaging; if false,
                       point-samples j_L(k*r[i]) at every grid point.
*/
void fillBesselVec_kr(int l, double k, const std::vector<double> &r,
                      std::vector<double> *jl, bool cell_average = true);

//==============================================================================
/*!
  @brief Lookup table of spherical Bessel functions: j_L(q*r) = J[L][q][r].
  @details
  Pre-computes and stores j_L(q*r) for a grid of q and r values.
  Provides first-match, nearest, and linear-interpolation lookup.
  All lookups are only accurate if the q grid is dense enough.
  Ideally, use the same r grid to extract values as was used to fill the table.
*/
class JL_table {
public:
  //! Default construct empty table; must call fill() before use
  JL_table() {}
  //! Construct and fill the table; equivalent to default-construct + fill()
  JL_table(int max_L, const std::vector<double> &q,
           const std::vector<double> &r, bool cell_average = true) {
    fill(max_L, q, r, cell_average);
  }

  //! If derivatives required for degree L, max_L = L+1.
  //! cell_average is forwarded to fillBesselVec_kr (see its docs); default true.
  void fill(int max_L, const std::vector<double> &q,
            const std::vector<double> &r, bool cell_average = true) {
    m_q = q;
    m_J_L_q.resize(std::size_t(max_L + 1), q.size());
    m_J_L_q_on_qr.resize(std::size_t(max_L + 1), q.size());
    using namespace qip::overloads;
#pragma omp parallel for schedule(dynamic)
    for (auto iq = 0ul; iq < m_J_L_q.cols(); ++iq) {
      for (auto L = 0ul; L <= std::size_t(max_L); ++L) {
        const auto tq = q[iq];
        // Fill jL(qr) (cell-averaged if cell_average and grid is too coarse)
        m_J_L_q[L][iq] = fillBesselVec_kr(int(L), tq, r, cell_average);
        // Store jL(qr)/qr
        m_J_L_q_on_qr[L][iq] = m_J_L_q[L][iq] / (tq * r);
      }
    }
  }

  //----------------------------------------------------------------------------
  //! Direct access to jL(q*r) for specific q grid index
  const std::vector<double> &at(std::size_t L, std::size_t iq) const {
    return m_J_L_q.atc(L, iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q_i*r) for the first grid point q_i such that q_i >= q
  const std::vector<double> &jL(int L, double q) const {
    // The '-1' here ensures we always get a valid index, even if q > m_q.back()
    const auto it = std::lower_bound(m_q.begin(), m_q.end() - 1, q);
    const auto iq = std::size_t(std::distance(m_q.begin(), it));
    return m_J_L_q.atc(std::size_t(L), iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q_i*r) for the grid point q_i nearest to the requested q.
  const std::vector<double> &jL_nearest(int L, double q) const {

    // Clamp
    if (q <= m_q.front())
      return m_J_L_q.atc(std::size_t(L), 0);
    if (q >= m_q.back())
      return m_J_L_q.atc(std::size_t(L), m_q.size() - 1);

    // Find i such that m_q[i] <= q < m_q[i+1]
    const auto it = std::lower_bound(m_q.begin() + 1, m_q.end() - 1, q);
    const auto i = std::size_t(std::distance(m_q.begin(), it)) - 1;

    // Pick nearest
    const auto iq = (q - m_q[i] < m_q[i + 1] - q) ? i : i + 1;

    return m_J_L_q.atc(std::size_t(L), iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q_i*r)/(q_i*r) for the grid point q_i nearest to the requested q.
  const std::vector<double> &jL_on_qr_nearest(int L, double q) const {

    // Clamp
    if (q <= m_q.front())
      return m_J_L_q_on_qr.atc(std::size_t(L), 0);
    if (q >= m_q.back())
      return m_J_L_q_on_qr.atc(std::size_t(L), m_q.size() - 1);

    // Find i such that m_q[i] <= q < m_q[i+1]
    const auto it = std::lower_bound(m_q.begin() + 1, m_q.end() - 1, q);
    const auto i = std::size_t(std::distance(m_q.begin(), it)) - 1;

    // Pick nearest
    const auto iq = (q - m_q[i] < m_q[i + 1] - q) ? i : i + 1;

    return m_J_L_q_on_qr.atc(std::size_t(L), iq);
  }

  //----------------------------------------------------------------------------
  //! Returns jL(q*r) interpolated linearly between grid points.
  //! nb: assumes q grid dense enough that linear interp is reasonable!
  std::vector<double> jL_interp(int L, double q) const {

    // Clamp q to the tabulated range
    if (q <= m_q.front())
      return m_J_L_q.atc(std::size_t(L), 0);
    if (q >= m_q.back())
      return m_J_L_q.atc(std::size_t(L), m_q.size() - 1);

    // Find i such that m_q[i] <= q < m_q[i+1]
    const auto it = std::lower_bound(m_q.begin() + 1, m_q.end() - 1, q);
    const auto i = std::size_t(std::distance(m_q.begin(), it)) - 1;

    // Linear interpolation weight
    // const double t = (q - m_q[i]) / (m_q[i + 1] - m_q[i]);

    const double logq = std::log10(q);
    const double logq0 = std::log10(m_q[i]);
    const double logq1 = std::log10(m_q[i + 1]);
    const double t = (logq - logq0) / (logq1 - logq0);

    const auto &v0 = m_J_L_q.atc(std::size_t(L), i);
    const auto &v1 = m_J_L_q.atc(std::size_t(L), i + 1);

    using namespace qip::overloads;
    return (1.0 - t) * v0 + t * v1;
  }

  //----------------------------------------------------------------------------
private:
  LinAlg::Matrix<std::vector<double>> m_J_L_q{};
  LinAlg::Matrix<std::vector<double>> m_J_L_q_on_qr{};
  std::vector<double> m_q{};
};

} // namespace SphericalBessel
