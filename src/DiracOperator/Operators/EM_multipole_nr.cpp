#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {

//==============================================================================
//! @brief Electric multipole operator, in small qr limit

class E5k_w_nr final : public TensorOperator {
public:
  E5k_w_nr(const Grid &gr, int K, double omega)
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
    return Angular::Ck_kk(1, ka, -kb);
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

    //const auto K = double(m_K);
    assert(m_K != 0); // should already be discounted!
    const auto cx = std::sqrt(2.0) / 3.0;
    const auto dk = double(Fa.kappa() + Fb.kappa());
    std::vector<double> ones(Fa.grid().size(), 1);

    return cx * (dk * Rab(-1, ones, Fa, Fb) - Rab(+1, ones, Fa, Fb));
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
//! @brief Longitudanal multipole operator, in small qr limit
class L5k_w_nr final : public TensorOperator {
public:
  L5k_w_nr(const Grid &gr, int K, double omega)
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
    return Angular::Ck_kk(1, ka, -kb);
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
    std::vector<double> ones(Fa.grid().size(), 1);

    return -(1.0 / 3.0) *
           (dk * Rab(-1, ones, Fa, Fb) - Rab(+1, ones, Fa, Fb)); // sigma = -1
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
//! @brief Magnetic multipole operator, in small qr limit
class M5k_w_nr final : public TensorOperator {
public:
  M5k_w_nr(const Grid &gr, int K, double omega)
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
    return Angular::Ck_kk(1, ka, kb);
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
    //const auto ck = sk / std::sqrt(K * (K + 1.0));

    return -(1.0 / std::sqrt(2.0)) * sk * Rab(-1, j1, Fa, Fb); // sigma = 0
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
//! @brief Temporal component of vector multipole operator, in small qr limit
class Phi5k_w_nr final : public TensorOperator {
public:
  Phi5k_w_nr(const Grid &gr, int K, double omega)
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
    return Angular::Ck_kk(1, ka, -kb);
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

    SphericalBessel::fillBesselVec_kr(1, q, m_vec, &jk);
  }

private:
  int m_K;
  std::vector<double> jk{};
};

//==============================================================================

//! Helper functions for the multipole operators
namespace multipole {

// For the small qr limit
inline std::unique_ptr<DiracOperator::TensorOperator>
V5_sigma_nr(const Grid &grid, int k, int sigma) {

  switch (sigma) {

  case +1:
    return std::make_unique<E5k_w_nr>(grid, 1, 1.0e-4);

  case -1:
    return std::make_unique<L5k_w_nr>(grid, 1, 1.0e-4);

  case 0:
    return std::make_unique<M5k_w_nr>(grid, 1, 1.0e-4);
  }
}

} // namespace multipole

} // namespace DiracOperator