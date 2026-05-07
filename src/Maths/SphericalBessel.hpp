#pragma once
#include "LinAlg/Matrix.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <vector>

/*!
@brief Wrappers for returning Spherical Bessel functions.
@details Has an "exact" version, and a faster approx version (good to ~ 1.e-9)
*/
namespace SphericalBessel {

//==============================================================================
//! Spherical Bessel function; uses low-x approximation where appropriate
// Good to parts in _ for x and _ otherwise
double JL(int L, double x);

//! Spherical Bessel function of first kind via GSL
double exactGSL_JL(int L, double x);

//! Spherical Bessel functions of second kind
double exactGSL_YL(int L, double x);

//! Used for frequency-dependent Breit - more generic formula that can be used to evaluate spherical Bessel function of second kind for negative orders
double yL(int L, double x);

//==============================================================================
//! Spherical Bessel function of second kind, y_L(x)
double YL(int L, double x);

//! Phi_L(x) = [(2L+1)!! / x^L] * j_L(x). If tilde=true, returns 1 - Phi_L(x)
double PhiL(int L, double x, bool tilde = false);

//! Psi_L(x) = [x^{L+1} / (2L-1)!!] * y_L(x). If tilde=true, returns 1 - Psi_L(x)
double PsiL(int L, double x, bool tilde = false);

//==============================================================================
//! Creates a vector of Jl(r) for given r
std::vector<double> fillBesselVec(int l, const std::vector<double> &xvec);

std::vector<double> fillBesselVec_kr(int l, double k,
                                     const std::vector<double> &rvec);

void fillBesselVec_kr(int l, double k, const std::vector<double> &r,
                      std::vector<double> *jl);

//==============================================================================
//! Spherical Bessel Lookup table; in the form j_L(qr) = J[L][q][r].
/*!
@details
 - Allows 'first match', 'nearest', and 'interpolation' lookup
 - All of these are only accurate if grid dense enough
 - Ideally, use exact same grid to extract values as used to create them
*/
class JL_table {
public:
  JL_table() {}
  JL_table(int max_L, const std::vector<double> &q,
           const std::vector<double> &r) {
    fill(max_L, q, r);
  }

  //! If derivatives required for degree L, max_L = L+1
  void fill(int max_L, const std::vector<double> &q,
            const std::vector<double> &r) {
    m_q = q;
    m_J_L_q.resize(std::size_t(max_L + 1), q.size());
    m_J_L_q_on_qr.resize(std::size_t(max_L + 1), q.size());
    using namespace qip::overloads;
#pragma omp parallel for collapse(2)
    for (auto L = 0ul; L <= std::size_t(max_L); ++L) {
      for (auto iq = 0ul; iq < m_J_L_q.cols(); ++iq) {
        const auto tq = q[iq];
        // Fill jL(qr)
        m_J_L_q[L][iq] = fillBesselVec_kr(int(L), tq, r);
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
