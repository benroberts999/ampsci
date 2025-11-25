#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {

// NOTE: These are quite ineficiant - calculate jK from scratch each time
// Better to give an option to use a lookup table!

//==============================================================================
//! @brief Electric multipole operator, alpha.t^1, L-form.
/*! @details
- in terms of j(q*r), where q = alpha*omega
- If transition_form is true, will use the "t" (transition form); if false, will use the "moment" form. Nb: E1, E2 etc. is the "moment" form, but alpha exp(iqr) goes to transition form - see Eq. (6.129) in Johnson for factors!
- Note: tk (transition) version, different convention to Johnson. Moment version the same!
*/
class Ek_w_L final : public TensorOperator {
public:
  Ek_w_L(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("t^E_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string("1"); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();
    if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto K = double(m_K);
    const auto c1 = double(kappa_a - Fb.kappa()) / (K + 1.0);
    const auto cx = std::sqrt((K + 1.0) / K);
    Rab_rhs(+1, j1, &dF, Fb, cx);
    Pab_rhs(+1, j2, &dF, Fb, -c1 * cx);
    Pab_rhs(-1, j2, &dF, Fb, -cx);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto cc = double(Fa.kappa() - Fb.kappa()) / (K + 1.0);
    const auto cx = std::sqrt((K + 1.0) / K);

    return cx * (Rab(+1, j1, Fa, Fb) - (cc + 1.0) * Vab(j2, Fa, Fb) +
                 (1.0 - cc) * Wab(j2, Fa, Fb));
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
  }

private:
  int m_K;
  std::vector<double> j1{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Electric multipole operator, V-form, including frequency-dependence.
/*! @details
- in terms of j(q*r), where q = alpha*omega
- If transition_form is true, will use the "t" (transition form); if false, will use the "moment" form. Nb: E1, E2 etc. is the "moment" form, but alpha exp(iqr) goes to transition form - see Eq. (6.129) in Johnson for factors!
- Note: tk (transition) version, different convention to Johnson. Moment version the same!
*/
class Ek_w final : public TensorOperator {
public:
  Ek_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^E_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto K = double(m_K);
    const auto dk = double(kappa_a - Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto cx = std::sqrt((K + 1.0) / K);

    Pab_rhs(+1, j1_on_qr, &dF, Fb, cx * dk);
    Pab_rhs(+1, j2, &dF, Fb, -cx * dk / (K + 1.0));
    Pab_rhs(-1, j1_on_qr, &dF, Fb, -cx * K);

    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
      return 0.0;
    }

    const auto K = double(m_K);
    assert(m_K != 0); // should already be discounted!
    const auto cx = std::sqrt((K + 1.0) / K);
    const auto dk = double(Fa.kappa() - Fb.kappa());

    const auto Pp1 = Pab(+1, j1_on_qr, Fa, Fb);
    const auto Pp2 = Pab(+1, j2, Fa, Fb);
    const auto Pm1 = Pab(-1, j1_on_qr, Fa, Fb);

    return cx * (dk * (Pp1 - Pp2 / (K + 1)) - K * Pm1);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      j1_on_qr[i] /= (q * m_vec[i]);
    }
  }

private:
  int m_K;
  std::vector<double> j1_on_qr{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Longitudanal multipole operator, V-form, including frequency-dependence.
class Lk_w final : public TensorOperator {
public:
  Lk_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^L_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
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

    const auto K = double(m_K);
    const auto dk = double(kappa_a - Fb.kappa());

    Pab_rhs(+1, j1_on_qr, &dF, Fb, -dk);
    Pab_rhs(-1, j1_on_qr, &dF, Fb, K);
    Pab_rhs(-1, j2, &dF, Fb, -1.0);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto dk = double(Fa.kappa() - Fb.kappa());

    return -dk * Pab(+1, j1_on_qr, Fa, Fb) + K * Pab(-1, j1_on_qr, Fa, Fb) -
           Pab(-1, j2, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      j1_on_qr[i] /= (q * m_vec[i]);
    }
  }

private:
  int m_K;
  std::vector<double> j1_on_qr{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Magnetic multipole operator, including frequency-dependence.
class Mk_w final : public TensorOperator {
public:
  Mk_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }

  std::string name() const override final {
    return std::string("t^M_") + std::to_string(m_K);
  }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, -ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a == -Fb.kappa()) || m_K == 0) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto K = double(m_K);
    const auto sk = double(kappa_a + Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto ck = sk / std::sqrt(K * (K + 1.0));

    Pab_rhs(+1, j1, &dF, Fb, -ck);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa()) ||
        m_K == 0) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto sk = double(Fa.kappa() + Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto ck = sk / std::sqrt(K * (K + 1.0));

    return -ck * Pab(+1, j1, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1);
  }

private:
  int m_K;
  std::vector<double> j1{};
};

//==============================================================================
//! @brief Temporal component of vector multipole operator, V-form, including frequency-dependence.
class Phik_w final : public TensorOperator {
public:
  Phik_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^V_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
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

    Rab_rhs(+1, jk, &dF, Fb);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Rab(+1, jk, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  int m_K;
  std::vector<double> jk{};
};

//==============================================================================
//! @brief Scalar multipole operator, e^{iqr}gamma^0, including frequency-dependence.
class Sk_w final : public TensorOperator {
public:
  Sk_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^S_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
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

    Rab_rhs(-1, jk, &dF, Fb);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Rab(-1, jk, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  int m_K;
  std::vector<double> jk{};
};

//==============================================================================
//==============================================================================
// Gamma^5 versions!

//==============================================================================
//! @brief Electric multipole operator, V-form, including frequency-dependence.

class E5k_w final : public TensorOperator {
public:
  E5k_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^E5_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa()) || m_K == 0) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto K = double(m_K);
    const auto dk = double(kappa_a + Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto cx = std::sqrt((K + 1.0) / K);

    Rab_rhs(-1, j1_on_qr, &dF, Fb, cx * dk);
    Rab_rhs(-1, j2, &dF, Fb, -cx * dk / (K + 1.0));
    Rab_rhs(+1, j1_on_qr, &dF, Fb, -cx * K);

    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || m_K == 0) {
      return 0.0;
    }

    const auto K = double(m_K);
    assert(m_K != 0); // should already be discounted!
    const auto cx = std::sqrt((K + 1.0) / K);
    const auto dk = double(Fa.kappa() + Fb.kappa());

    const auto Pp1 = Rab(-1, j1_on_qr, Fa, Fb);
    const auto Pp2 = Rab(-1, j2, Fa, Fb);
    const auto Pm1 = Rab(+1, j1_on_qr, Fa, Fb);

    return cx * (dk * (Pp1 - Pp2 / (K + 1)) - K * Pm1);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      j1_on_qr[i] /= (q * m_vec[i]);
    }
  }

private:
  int m_K;
  std::vector<double> j1_on_qr{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Longitudanal multipole operator, V-form, including frequency-dependence.
class L5k_w final : public TensorOperator {
public:
  L5k_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^5L_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
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

    const auto K = double(m_K);
    const auto dk = double(kappa_a + Fb.kappa());

    Rab_rhs(-1, j1_on_qr, &dF, Fb, -dk);
    Rab_rhs(+1, j1_on_qr, &dF, Fb, K);
    Rab_rhs(+1, j2, &dF, Fb, -1.0);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto dk = double(Fa.kappa() + Fb.kappa());

    return -dk * Rab(-1, j1_on_qr, Fa, Fb) + K * Rab(+1, j1_on_qr, Fa, Fb) -
           Rab(+1, j2, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      j1_on_qr[i] /= (q * m_vec[i]);
    }
  }

private:
  int m_K;
  std::vector<double> j1_on_qr{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Magnetic multipole operator, including frequency-dependence.
class M5k_w final : public TensorOperator {
public:
  M5k_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }

  std::string name() const override final {
    return std::string("t^5M_") + std::to_string(m_K);
  }

  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, kb);
  }

  //--------------
  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());
    dF.min_pt() = Fb.min_pt();
    dF.max_pt() = Fb.max_pt();

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a == -Fb.kappa()) || m_K == 0) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto K = double(m_K);
    const auto sk = double(kappa_a - Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto ck = sk / std::sqrt(K * (K + 1.0));

    Rab_rhs(-1, j1, &dF, Fb, -ck);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa()) ||
        m_K == 0) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto sk = double(Fa.kappa() - Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto ck = sk / std::sqrt(K * (K + 1.0));

    return -ck * Rab(-1, j1, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1);
  }

private:
  int m_K;
  std::vector<double> j1{};
};

//==============================================================================
//! @brief Temporal component of vector multipole operator, V-form, including frequency-dependence.
class Phi5k_w final : public TensorOperator {
public:
  Phi5k_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("t^5_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
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

    Pab_rhs(-1, jk, &dF, Fb);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Pab(-1, jk, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  int m_K;
  std::vector<double> jk{};
};

//==============================================================================
//! @brief Pseudoscalar multipole operator, ~ $e^{iqr}\gamma^0\gamma^5$
class S5k_w final : public TensorOperator {
public:
  S5k_w(const Grid &gr, int K, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 1.0,
                       gr.r(), 0, Realness::real, true),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^5S_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return Angular::Ck_kk(m_K, ka, -kb);
  }

  //--------------
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

    Pab_rhs(+1, jk, &dF, Fb);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    return Pab(+1, jk, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(PhysConst::alpha * omega);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  int m_K;
  std::vector<double> jk{};
};
//==============================================================================

//! Helper functions for the multipole operators
namespace multipole {

//! Convert from "transition form" to "moment form"
inline double moment_factor(int K, double omega) {
  const auto q = std::abs(PhysConst::alpha * omega);
  return qip::double_factorial(2 * K + 1) / qip::pow(q, K) *
         std::sqrt(K / (K + 1.0));
}

//! Scalar(pseudoscalar) multipole operator. Note: q/omega not set
inline std::unique_ptr<DiracOperator::TensorOperator>
S_K(const Grid &grid, int k, bool gamma5 = false) {
  if (gamma5)
    return std::make_unique<S5k_w>(grid, k, 1.0e-4);
  return std::make_unique<Sk_w>(grid, k, 1.0e-4);
}

//! Temporal part of vector(pseudovector) multipole operator. Note: q/omega not set
inline std::unique_ptr<DiracOperator::TensorOperator>
Phi_K(const Grid &grid, int k, bool gamma5 = false) {
  if (gamma5)
    return std::make_unique<Phi5k_w>(grid, k, 1.0e-4);
  return std::make_unique<Phik_w>(grid, k, 1.0e-4);
}

//! Spatial part of vector(pseudovector) multipole operator. Note: q/omega not set
inline std::unique_ptr<DiracOperator::TensorOperator>
V_sigma_K(const Grid &grid, int sigma, int k, bool gamma5 = false) {

  switch (sigma) {

    // "Electric"
  case +1:
    if (gamma5)
      return std::make_unique<E5k_w>(grid, k, 1.0e-4);
    else
      return std::make_unique<Ek_w>(grid, k, 1.0e-4);

    // "Longitudanal"
  case -1:
    if (gamma5)
      return std::make_unique<L5k_w>(grid, k, 1.0e-4);
    else
      return std::make_unique<Lk_w>(grid, k, 1.0e-4);

    // "Magnetic"
  case 0:
    if (gamma5)
      return std::make_unique<M5k_w>(grid, k, 1.0e-4);
    else
      return std::make_unique<Mk_w>(grid, k, 1.0e-4);
  }

  assert(false && "make_op: sigma must be -1, 0, or +1");
  return {}; // never reached when assertions are enabled
}

} // namespace multipole

//==============================================================================
//==============================================================================

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek_w_L(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Ek_w_L>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [0.001]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.001);
  return std::make_unique<Ek_w>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Mk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for M1, =2 for M2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Mk_w>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Lk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"k", "Rank [1]"}, {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Lk_w>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Vk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"k", "Rank [1]"}, {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Phik_w>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Sk_w(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"k", "Rank [1]"}, {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Sk_w>(wf.grid(), k, omega);
}

} // namespace DiracOperator
