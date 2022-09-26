#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//! Magnetic dipole operator: <a||M1||b>
/*! @details
\f[ <a||M1||b> = R (k_a + k_b) <-k_a||C^1||k_b> \f]
\f[R = \frac{-3}{\alpha^2\omega} \int (f_ag_b+g_af_b) j_1(kr) \, dr\f]
\f$ k = \omega/c = \omega*\alpha \f$.
Negative sign (and alpha) puts into units |mu_B|.
For k<<1 (static case): j1(kr) -> (r*k)/3,
\f[R = \frac{-1}{\alpha} \int (f_ag_b+g_af_b) r \, dr\f]
*/
class M1 final : public TensorOperator {
public:
  M1(const Grid &gr, const double alpha, const double omega = 0.0)
      : TensorOperator(1, Parity::even, +0.0, {}, 0, Realness::real, true),
        m_r(gr.r()),
        m_alpha(alpha) {
    updateFrequency(omega);
  }
  M1 &operator=(const M1 &) = delete;
  M1(const M1 &) = default;
  ~M1() = default;
  std::string name() const override final { return std::string("M1"); }
  std::string units() const override final { return std::string("|mu_B|"); }

  double angularF(const int ka, const int kb) const override final {
    return (ka + kb) * Angular::Ck_kk(1, -ka, kb);
  }
  double angularCff(int, int) const override final { return 0.0; }
  double angularCgg(int, int) const override final { return 0.0; }
  double angularCfg(int, int) const override final { return 1.0; }
  double angularCgf(int, int) const override final { return 1.0; }

  void updateFrequency(const double omega) override final {
    if (std::abs(omega) > 1.0e-10) {
      m_constant = -3.0 / (m_alpha * m_alpha * omega);
      m_vec = SphericalBessel::fillBesselVec_kr(1, omega * m_alpha, m_r);
    } else {
      // j1(kr) -> (r*k)/3, for k<<1
      m_constant = -1.0 / (m_alpha);
      m_vec = m_r;
    }
  }

private:
  const std::vector<double> m_r; // store radial vector (copy)
  const double m_alpha;
};

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_M1(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"", "No input options"}});
  return std::make_unique<M1>(wf.grid(), wf.alpha(), 0.0);
}

} // namespace DiracOperator
