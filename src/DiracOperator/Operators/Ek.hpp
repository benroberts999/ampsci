#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//==============================================================================
//! E^k (electric multipole) operator, length form, with qr<<1 (static) approximation
/*! @details
  \f[ 
    h_k = -|e|r^k = -e r^k 
  \f]
  \f[
    \redmatel{a}{h_k}{b} = R C^k_{ab}
  \f]
  \f[
    R = -e \int r^k (f_a f_b + g_a g_b) \, dr
  \f]

  - This has \f$ qr \\ 1\f$ (static) approximation; 
  - That is, no bessel functions
  - See @ref EM_multipole operator for full operators
*/
class Ek : public TensorOperator {
public:
  Ek(const Grid &gr, const int k)
    : TensorOperator(k, Angular::evenQ(k) ? Parity::even : Parity::odd, -1.0,
                     gr.rpow(k)),
      m_k(k) {}

  std::unique_ptr<TensorOperator> clone() const override {
    return std::make_unique<Ek>(*this);
  }

  std::string name() const override {
    return std::string("E") + std::to_string(m_k);
  }
  std::string units() const override {
    return m_k == 1 ? "|e|aB" : std::string("|e|aB^") + std::to_string(m_k);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_k, ka, kb);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"}});
    if (input.has_option("help"))
      return nullptr;
    const auto k = input.get("k", 1);
    return std::make_unique<Ek>(wf.grid(), k);
  }

private:
  int m_k;
};

//==============================================================================
//! Electric dipole operator: -|e|r = -er
/*! @details
\f[<a||d||b> = R (-1)^{ja+1/2} \sqrt{[ja][jb]} \, tjs(ja,jb,1, -1/2,1/2,0)\f]
\f[<a||d||b> = R <k_a||C^k||k_b>\f]
\f[R = -e \int r(f_a f_b + g_a g_b) \, dr\f]
*/
class E1 final : public Ek {
public:
  E1(const Grid &gr) : Ek(gr, 1) {}

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<E1>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"no options", ""}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<E1>(wf.grid());
  }
};

//==============================================================================
//! Odd-parity rank-0 operator with radial function r (sigma_r)
class sigma_r final : public ScalarOperator {
public:
  sigma_r(const Grid &rgrid) : ScalarOperator(Parity::odd, -1.0, rgrid.r()) {}
  std::string name() const override final { return "s.r"; }
  std::string units() const override final { return "aB"; }

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<sigma_r>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"no options", ""}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<sigma_r>(wf.grid());
  }
};

//==============================================================================
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
    : TensorOperator(1, Parity::odd, -0.0, {}, Realness::real, true),
      m_alpha(alpha) {
    updateFrequency(omega);
  }
  std::string name() const override final { return "E1v"; }
  std::string units() const override final { return "|e|aB"; }

  double angularF(const int ka, const int kb) const override final {
    // return 1;
    return Angular::Ck_kk(1, ka, kb);
  }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int ka, int kb) const override final {
    // return Angular::S_kk(ka, -kb);
    return ka - kb - 1;
  }
  double angularCgf(int ka, int kb) const override final {
    // return -Angular::S_kk(-ka, kb);
    return ka - kb + 1;
  }

  //! @note At omega=0 (static limit) falls back to m_constant = -1.0;
  //! the velocity-form operator diverges as 1/omega and is not physically
  //! meaningful at zero frequency.
  void updateFrequency(const double omega) override final {
    m_omega = omega;
    // m_constant = std::abs(omega) > 1.0e-10 ? -2.0 / (m_alpha * omega) : 1.0;
    m_constant = std::abs(omega) > 1.0e-10 ? -1.0 / (m_alpha * omega) : -1.0;
  }

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<E1v>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"no options", ""}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<E1v>(wf.alpha(), 0.0);
  }

private:
  double m_alpha; // (including var-alpha)
};

//==============================================================================
//! @brief i\vec{\alpha} matrix element: propto Electric dipole operator.
//! i factored so that ME is real
class ialpha final : public TensorOperator {
public:
  ialpha() : TensorOperator(1, Parity::odd, 1.0, {}, Realness::real, false) {}

  std::string name() const override final { return "i*alpha"; }
  std::string units() const override final { return "au"; }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(1, ka, kb);
  }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int ka, int kb) const override final { return ka - kb - 1; }
  double angularCgf(int ka, int kb) const override final { return ka - kb + 1; }

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<ialpha>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &) {
    input.check({{"no options", ""}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<ialpha>();
  }
};

//==============================================================================
//! @brief Electric quadrupole operator: -|e|r^2
class E2 final : public Ek {
public:
  E2(const Grid &gr) : Ek(gr, 2) {}

  std::unique_ptr<TensorOperator> clone() const override final {
    return std::make_unique<E2>(*this);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"no options", ""}});
    if (input.has_option("help"))
      return nullptr;
    return std::make_unique<E2>(wf.grid());
  }
};

} // namespace DiracOperator
