#pragma once
#include "DiracOperator/TensorOperator.hpp"

namespace DiracOperator {

//! E^k (electric multipole) operator
/*! @details
\f[ h = -|e|r^k = -e r^k \f]
\f[<a||d||b> = R C^k_{ab}\f]
\f[R = -e \int r^k (f_a f_b + g_a g_b) \, dr\f]
*/
class Ek : public TensorOperator {
public:
  Ek(const Grid &gr, const int k)
      : TensorOperator(k, Angular::evenQ(k) ? Parity::even : Parity::odd, -1.0,
                       gr.rpow(k), 0),
        m_k(k) {}
  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_k, ka, kb);
  }
  std::string name() const override {
    return std::string("E") + std::to_string(m_k);
  }
  std::string units() const override {
    return m_k == 1 ? "aB" : std::string("aB^") + std::to_string(m_k);
  }

private:
  int m_k;
};

//******************************************************************************
//! Electric dipole operator: -|e|r = -er
/*! @details
\f[<a||d||b> = R (-1)^{ja+1/2} \sqrt{[ja][jb]} \, tjs(ja,jb,1, -1/2,1/2,0)\f]
\f[<a||d||b> = R <k_a||C^k||k_b>\f]
\f[R = -e \int r(f_a f_b + g_a g_b) \, dr\f]
*/
class E1 final : public Ek {
public:
  E1(const Grid &gr) : Ek(gr, 1) {}
};

//******************************************************************************
//! @brief Electric dipole operator, v-form:
//! \f$ \frac{ie}{\omega \alpha} \vec{\alpha}\f$
/*! @details
\f[  <a||d_v||b> =  R \f]
\f[ R = -\frac{2e}{\omega \alpha}
\int( f_ag_b <ka||s||-kb> - g_af_b <-ka||s||kb>)\,dr \f]
*/
class E1v final : public TensorOperator
// d_v = (ie/w alpha) v{alpha}   [v{a} = g0v{g}]\f$
// <a||dv||b> = -2e/(w alpha) Int[ fagb <ka||s||-kb> - gafb <-ka||s||kb>]
{
public:
  E1v(const double alpha, const double omega = 0.0)
      : TensorOperator(1, Parity::odd, -0.0, {}, 0, Realness::real, true),
        m_alpha(alpha) {
    updateFrequency(omega);
  }
  std::string name() const override final { return "E1v"; }
  std::string units() const override final { return "aB"; }

  double angularF(const int, const int) const override final { return 1.0; }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int ka, int kb) const override final {
    return Angular::S_kk(ka, -kb);
  }
  double angularCgf(int ka, int kb) const override final {
    return -Angular::S_kk(-ka, kb);
  }

  void updateFrequency(const double omega) override final {
    m_constant = std::abs(omega) > 1.0e-10 ? -2.0 / (m_alpha * omega) : 1.0;
  }

private:
  double m_alpha; // (including var-alpha)
};

} // namespace DiracOperator
