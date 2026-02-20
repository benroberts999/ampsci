#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {

//! Low qr form of Vector electric. Only for K=1 (zero otherwise)
class VEk_lowq final : public TensorOperator {
public:
  VEk_lowq(const Grid &, int K, double)
    : TensorOperator(1, Parity::odd, -std::sqrt(2.0) / 3.0, {}, 0,
                     Realness::imaginary, false),
      m_K(K) {}

  std::string name() const override final {
    return std::string("V_lowq^E_") + std::to_string(m_rank);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(1, ka, kb);
  }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int ka, int kb) const override final {
    return m_K == 1 ? ka - kb - 1 : 0.0;
  }
  double angularCgf(int ka, int kb) const override final {
    return m_K == 1 ? ka - kb + 1 : 0.0;
  }

private:
  int m_K;
};

//! Low qr form of Vector magnetic. Only for K=1 (zero otherwise)
class VMk_lowq final : public TensorOperator {
public:
  VMk_lowq(const Grid &gr, int K, double omega)
    : TensorOperator(1, Parity::even, 0, gr.r(), 0, Realness::imaginary, true),
      m_K(K) {
    updateFrequency(omega);
  }

  std::string name() const override final {
    return std::string("V_lowq^E_") + std::to_string(m_rank);
  }

  double angularF(const int ka, const int kb) const override final {
    return (ka + kb) * Angular::Ck_kk(1, ka, -kb);
  }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int, int) const override final {
    return m_K == 1 ? 1.0 : 0.0;
  }
  double angularCgf(int, int) const override final {
    return m_K == 1 ? 1.0 : 0.0;
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = std::abs(PhysConst::alpha * omega);
    m_constant = -m_q / (3.0 * std::sqrt(2.0));
  }

private:
  int m_K;
  double m_q{};
};

//! Low qr form of Vector longitudanal
class VLk_lowq final : public TensorOperator {
public:
  VLk_lowq(const Grid &, int K, double)
    : TensorOperator(1, Parity::odd, 1.0 / 3.0, {}, 0, Realness::imaginary,
                     false),
      m_K(K) {}

  std::string name() const override final {
    return std::string("V_lowq^E_") + std::to_string(m_rank);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(1, ka, kb);
  }

  double angularCff(int, int) const override final { return 0; }
  double angularCgg(int, int) const override final { return 0; }
  double angularCfg(int ka, int kb) const override final {
    return m_K == 1 ? ka - kb - 1 : 0.0;
  }
  double angularCgf(int ka, int kb) const override final {
    return m_K == 1 ? ka - kb + 1 : 0.0;
  }

private:
  int m_K;
};

//==============================================================================
//! Low qr form of Temporal component of the vector multipole operator:
//! \f$ \Phi_K = t^K(q)\f$
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
    return std::string("Phi_lowq_") + std::to_string(m_K);
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
//! @brief Low qr form of Scalar multipole operator: \f$ S_K = t^K(q)\gamma^0\f$
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
      Gab_rhs(&dF, Fb, -2.0);
    } else if (m_K == 1) {
      Rab_rhs(-1, m_vec, &dF, Fb, (m_q / 3.0));
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    if (m_K == 0) {
      return -2.0 * Gab(Fa, Fb);
    }

    if (m_K == 1) {
      return (m_q / 3.0) * Rab(-1, m_vec, Fa, Fb);
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
};

//==============================================================================
//==============================================================================
// Gamma^5 versions!

//==============================================================================
//! @brief Low qr form of Axial electric multipole operator: \f$ A^E_K = T^{(+1)}_K(q)\gamma^5\f$
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
    return std::string("T^E5_lowq_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa()) || kappa_a == -Fb.kappa()) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto dk = double(kappa_a + Fb.kappa());
    const auto dk_int = kappa_a + Fb.kappa();
    const auto same_kap = kappa_a == Fb.kappa();

    if (m_K == 1) {
      const double cx = std::sqrt(2.0) / 3.0;

      if (same_kap || dk_int == 1) {
        Gab_rhs(&dF, Fb, -2.0 * cx * dk);
      } else {
        Rab_rhs(-1, &dF, Fb, cx * dk);
        Rab_rhs(+1, &dF, Fb, -cx);
      }
    }

    if (m_K == 2) {
      const auto cx2 = (1.0 / 15.0) * std::sqrt(3.0 / 2.0);

      if (dk_int == 2) {
        // F-part cancels!
        Gab_rhs(&dF, Fb, -4.0 * cx2 * m_q);
        dF *= m_vec;
      } else {
        Rab_rhs(-1, m_vec, &dF, Fb, cx2 * dk * m_q);
        Rab_rhs(+1, m_vec, &dF, Fb, -2.0 * cx2 * m_q);
      }
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || Fa.kappa() == -Fb.kappa()) {
      return 0.0;
    }

    const auto dk = double(Fa.kappa() + Fb.kappa());
    const auto dk_int = Fa.kappa() + Fb.kappa();

    const auto same_kap = Fa.kappa() == Fb.kappa();

    if (m_K == 1) {
      const auto cx = std::sqrt(2.0);

      if (same_kap || dk_int == 1) {
        return (-2.0 / 3.0) * cx * dk * Gab(Fa, Fb);
      }
      const auto Rm1 = Rab(-1, Fa, Fb);
      const auto Rp1 = Rab(+1, Fa, Fb);
      return cx * (dk * Rm1 - Rp1) / 3.0;
    }

    if (m_K == 2) {
      const auto cx2 = (1.0 / 15.0) * std::sqrt(3.0 / 2.0);

      if (dk_int == 2) {
        // F-part cancels!
        return cx2 * m_q * (-4.0 * Gab(m_vec, Fa, Fb));
      }

      return cx2 * m_q *
             (dk * Rab(-1, m_vec, Fa, Fb) - 2.0 * Rab(+1, m_vec, Fa, Fb));
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
};

//==============================================================================
//! @brief Low qr form of Axial longitudinal multipole operator: \f$ A^L_K = T^{(-1)}_K(q)\gamma^5\f$
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
      Rab_rhs(+1, m_vec, &dF, Fb, -(m_q / 3.0));
    } else if (m_K == 1) {
      const auto dk_int = kappa_a + Fb.kappa();
      const auto dk = double(dk_int);
      const auto same_kap = kappa_a == Fb.kappa();
      if (dk == 0)
        return dF;
      if (same_kap || dk_int == 1) {
        Gab_rhs(&dF, Fb, (2.0 / 3.0) * dk);
      } else {
        const auto c = -(1.0 / 3.0);
        // -(1.0 / 3.0) * (dk * Rab(-1, Fa, Fb) - Rab(+1, Fa, Fb));
        Rab_rhs(-1, &dF, Fb, c * dk);
        Rab_rhs(+1, &dF, Fb, -c);
      }
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    if (m_K == 0) {
      return -(m_q / 3.0) * Rab(+1, m_vec, Fa, Fb);
    }
    if (m_K == 1) {

      const auto dk_int = Fa.kappa() + Fb.kappa();
      const auto dk = double(dk_int);
      const auto same_kap = Fa.kappa() == Fb.kappa();

      if (same_kap || dk_int == 1)
        return (2.0 / 3.0) * dk * Gab(Fa, Fb);

      return -(1.0 / 3.0) * (dk * Rab(-1, Fa, Fb) - Rab(+1, Fa, Fb));
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
};

//==============================================================================
//! @brief Low qr form of Axial magnetic multipole operator: \f$ A^M_K = T^{(0)}_K(q)\gamma^5\f$
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
    return std::string("T^M5_k_lowq_1") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a == Fb.kappa())) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto sk = double(kappa_a - Fb.kappa());
    if (m_K == 1) {
      Rab_rhs(-1, m_vec, &dF, Fb, -(m_q / (std::sqrt(2.0) * 3.0)) * sk);
    }

    if (m_K == 2) {
      using namespace qip::overloads;
      const auto c = -(m_q * m_q / (std::sqrt(6.0) * 15.0)) * sk;
      Rab_rhs(-1, m_vec * m_vec, &dF, Fb, c);
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto sk = double(Fa.kappa() - Fb.kappa());

    if (m_K == 1) {
      const auto ck = sk / std::sqrt(2.0);
      return -ck * m_q * Rab(-1, m_vec, Fa, Fb) / 3.0;
    }
    if (m_K == 2) {
      using namespace qip::overloads;
      return -(m_q * m_q / (std::sqrt(6.0) * 15.0)) * sk *
             Rab(-1, m_vec * m_vec, Fa, Fb);
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
};

//==============================================================================
//! @brief Low qr form of Temporal component of the axial vector multipole operator: \f$ \Theta_K = \Phi^5_K = t^K(q)\gamma^5\f$ .
//! NOTE: If K=0, omega should be (ea-eb); for K=1 should be q = alpha*omega!
class Phi5k_lowq final : public TensorOperator {
public:
  Phi5k_lowq(const Grid &gr, int K, double omega,
             double alpha = PhysConst::alpha)
    : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                     gr.r(), 0, Realness::real, true),
      m_K(K),
      m_alpha(alpha) {
    if (omega != 0.0)
      updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("Phi5_lowq_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa())) {
      return dF;
    }

    if (m_K == 0) {
      // XXX q here is via _frequency_, not momentum
      const auto w_ab = m_q / m_alpha;
      const auto f_rel = 1.0 / (1.0 + m_alpha * m_alpha * kappa_a * w_ab);
      const auto c = -m_alpha * w_ab * f_rel;
      Rab_rhs(+1, m_vec, &dF, Fb, c);
    } else if (m_K == 1) {
      // XXX q here is _momentum_, not frequency
      Pab_rhs(-1, m_vec, &dF, Fb, (m_q / 3.0));
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    if (m_K == 0) {
      // XXX q here is via _frequency_, not momentum
      const auto w_ab = Fa.en() - Fb.en();
      const auto f_rel = 1.0 / (1.0 + m_alpha * m_alpha * Fa.kappa() * w_ab);
      return -m_alpha * w_ab * Rab(+1, m_vec, Fa, Fb) * f_rel;
      //Pab(-1, Fa, Fb); //
    }

    if (m_K == 1) {
      // XXX q here is _momentum_, not frequency
      return (m_q / 3.0) * Pab(-1, m_vec, Fa, Fb);
    }

    return 0.0;
  }

  //! NOTE: If K=0, omega should be (ea-eb); for K=1 should be q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = PhysConst::alpha * omega;
  }

private:
  int m_K;
  double m_alpha;
  double m_q{};
};

//==============================================================================
//! @brief Low qr form of Pseudoscalar multipole operator:
//! \f$ P_K = S^5_K = t^K(q)(i\gamma^0\gamma^5)\f$
//! NOTE: If K=0, omega should be (ea-eb); for K=1 should be q = alpha*omega!
class S5k_lowq final : public TensorOperator {
public:
  S5k_lowq(const Grid &gr, int K, double omega, double alpha = PhysConst::alpha)
    : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                     gr.r(), 0, Realness::real, true),
      m_K(K),
      m_alpha(alpha) {
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
                         const DiracSpinor &Fb) const override final {
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa())) {
      return dF;
    }

    if (m_K == 0) {
      // XXX Need to think about how to do omega here nicely
      // return dF;
      // XXX omega_ab factored out!! XXX

      const auto w_ab = m_q / m_alpha;
      const auto f_rel = 1.0 / (1.0 + m_alpha * m_alpha * kappa_a * w_ab);
      const auto wab2 = std::pow(w_ab, 2);
      const auto c = 0.5 * m_alpha * m_alpha * m_alpha * wab2 * f_rel;

      Rab_rhs(+1, m_vec, &dF, Fb, c);
    } else if (m_K == 1) {
      Pab_rhs(+1, m_vec, &dF, Fb, -(m_q / 3.0));
    }

    return dF;
  }

  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    if (m_K == 0) {

      // XXX q here is via _frequency_, not momentum
      const auto w_ab = Fa.en() - Fb.en();
      const auto f_rel = 1.0 / (1.0 + m_alpha * m_alpha * Fa.kappa() * w_ab);

      const auto wab2 = std::pow(w_ab, 2);
      return 0.5 * m_alpha * m_alpha * m_alpha * wab2 * Rab(+1, m_vec, Fa, Fb) *
             f_rel;
    }

    if (m_K == 1) {
      return -(m_q / 3.0) * Pab(+1, m_vec, Fa, Fb);
    }

    return 0.0;
  }

  //! NOTE: If K=0, omega should be (ea-eb); for K=1 should be q = alpha*omega!
  void updateFrequency(const double omega) override final {
    m_q = m_alpha * omega;
  }

private:
  int m_K;
  double m_alpha;
  double m_q{};
};

} // namespace DiracOperator
