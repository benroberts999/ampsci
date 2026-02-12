#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {
//==============================================================================
//! Low qr form of Temporal component of the vector multipole operator: $\Phi_K = t^K(q)$
class Phik_lowq final : public TensorOperator {
public:
  Phik_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K),
        m_r2(gr.rpow(2)) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("Phi_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa())) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    if (m_K == 0) {
      Rab_rhs(+1, m_r2, &dF, Fb, -(m_q * m_q / 6.0));
    } else if (m_K == 1) {
      Rab_rhs(+1, m_vec, &dF, Fb, (m_q / 3.0));
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    if (m_K == 0) {
      return -(m_q * m_q / 6.0) * Rab(+1, m_r2, Fa, Fb);
    }
    if (m_K == 1) {
      return (m_q / 3.0) * Rab(+1, m_vec, Fa, Fb);
    }

    return 0.0;
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
  std::vector<double> m_r2;
};

//==============================================================================
//! @brief Low qr form of Scalar multipole operator: $S_K = t^K(q)\gamma^0$
class Sk_lowq final : public TensorOperator {
public:
  Sk_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("t^S_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
};

//==============================================================================
//==============================================================================
// Gamma^5 versions!

//==============================================================================
//! @brief Low qr form of Axial electric multipole operator: $A^E_K = T^{(+1)}_K(q)\gamma^5$
class AEk_lowq final : public TensorOperator {
public:
  AEk_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real),
        m_K(K) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("T^E5_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
};

//==============================================================================
//! @brief Low qr form of Axial longitudinal multipole operator: $A^L_K = T^{(-1)}_K(q)\gamma^5$
class ALk_lowq final : public TensorOperator {
public:
  ALk_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("T^L5_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
};

//==============================================================================
//! @brief Low qr form of Axial magnetic multipole operator: $A^M_K = T^{(0)}_K(q)\gamma^5$
class AMk_lowq final : public TensorOperator {
public:
  AMk_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    if (omega != 0.0)
      updateFrequency(omega);
  }

  std::string name() const override final {
    return std::string("T^M5_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
};

//==============================================================================
//! @brief Low qr form of Temporal component of the axial vector multipole operator: $\Theta_K = \Phi^5_K = t^K(q)\gamma^5$
class Phi5k_lowq final : public TensorOperator {
public:
  Phi5k_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("Phi5_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
};

//==============================================================================
//! @brief Low qr form of Pseudoscalar multipole operator: $P_K = S^5_K = t^K(q)(i\gamma^0\gamma^5)$
class S5k_lowq final : public TensorOperator {
public:
  S5k_lowq(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("S5_k_lowq") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final;

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final;

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  int m_K;
  double m_q{};
};

} // namespace DiracOperator
