#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {

  //==============================================================================
//! @brief Electric multipole operator, in small qr limit
class E_nr final : public TensorOperator {
public:
  E_nr(const Grid &gr, double omega)
      : TensorOperator(1, Parity::even, 1.0, gr.r(), 0, Realness::real, false) {
        updateFrequency(omega);
  }
  std::string name() const override final { return std::string("tv^E_nr"); }
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

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a + Fb.kappa() == 1)) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto cx = std::sqrt(2.0) / 3.0;
    const auto dk = double(kappa_a + Fb.kappa());

    Rab_rhs(-1, &dF, Fb, cx * (dk - 1.0));
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || Fa.kappa() + Fb.kappa() == 1) {
      return 0.0;
    }

    const auto cx = std::sqrt(2.0) / 3.0;
    const auto dk = double(Fa.kappa() + Fb.kappa());

    return cx * (dk*Rab(-1, Fa, Fb) - Rab(+1, Fa, Fb));
  }
};

//==============================================================================
//! @brief Longitudanal multipole operator, in small qr limit
class L_nr final : public TensorOperator {
public:
  L_nr(const Grid &gr, double omega)
      : TensorOperator(1, Parity::even, 1.0, gr.r(), 0, Realness::real, false) {
        updateFrequency(omega);
  }
  std::string name() const override final { return std::string("tvL_nr"); }
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

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a + Fb.kappa() == 1)) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto dk = double(kappa_a + Fb.kappa());
    Rab_rhs(-1, &dF, Fb, -(1.0 / 3.0) * (dk - 1.0));

    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || Fa.kappa() + Fb.kappa() == 1) {
      return 0.0;
    }

    const auto dk = double(Fa.kappa() + Fb.kappa());

    return -(1.0 / 3.0) * (dk*Rab(-1, Fa, Fb) - Rab(+1, Fa, Fb));
  }
};

//==============================================================================
//! @brief Magnetic multipole operator, in small qr limit
class M_nr final : public TensorOperator {
public:
  M_nr(const Grid &gr, double omega)
      : TensorOperator(1, Parity::odd, 1.0, gr.r(), 0, Realness::real, true) {
    updateFrequency(omega);
  }

  std::string name() const override final { return std::string("tM_nr"); }

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

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a == -Fb.kappa())) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto sk = double(kappa_a - Fb.kappa());
    const auto ck = -(m_q / std::sqrt(2.0) / 3.0) * sk;

    Rab_rhs(-1, m_vec, &dF, Fb, ck);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa())) {
      return 0.0;
    }

    const auto sk = double(Fa.kappa() - Fb.kappa());
    return -(m_q /(std::sqrt(2.0)*3.0)) * sk * Rab(-1, m_vec, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  double m_q;
};

//==============================================================================
//! @brief Electric multipole operator, in small qr limit
class E5_nr final : public TensorOperator {
public:
  E5_nr(const Grid &gr, double omega)
      : TensorOperator(1, Parity::even, 1.0, gr.r(), 0, Realness::real, false) {
  }
  std::string name() const override final { return std::string("tv^E5_nr"); }
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

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a + Fb.kappa() == 1)) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto cx = std::sqrt(2.0) / 3.0;
    const auto dk = double(kappa_a + Fb.kappa());

    Rab_rhs(-1, &dF, Fb, cx * (dk - 1.0));
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || Fa.kappa() + Fb.kappa() == 1) {
      return 0.0;
    }

    const auto cx = std::sqrt(2.0) / 3.0;
    const auto dk = double(Fa.kappa() + Fb.kappa());

    return cx * (dk*Rab(-1, Fa, Fb) - Rab(+1, Fa, Fb));
  }
};

//==============================================================================
//! @brief Longitudanal multipole operator, in small qr limit
class L5_nr final : public TensorOperator {
public:
  L5_nr(const Grid &gr, double omega)
      : TensorOperator(1, Parity::even, 1.0, gr.r(), 0, Realness::real, false) {
        updateFrequency(omega);
  }
  std::string name() const override final { return std::string("tv^5L_nr"); }
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

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a + Fb.kappa() == 1)) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto dk = double(kappa_a + Fb.kappa());
    Rab_rhs(-1, &dF, Fb, -(1.0 / 3.0) * (dk - 1.0));

    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || Fa.kappa() + Fb.kappa() == 1) {
      return 0.0;
    }

    const auto dk = double(Fa.kappa() + Fb.kappa());

    return -(1.0 / 3.0) * (dk*Rab(-1, Fa, Fb) - Rab(+1, Fa, Fb));
  }
};

//==============================================================================
//! @brief Magnetic multipole operator, in small qr limit
class M5_nr final : public TensorOperator {
public:
  M5_nr(const Grid &gr, double omega)
      : TensorOperator(1, Parity::odd, 1.0, gr.r(), 0, Realness::real, true) {
    updateFrequency(omega);
  }

  std::string name() const override final { return std::string("t^5M_nr"); }

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

    if (isZero(kappa_a, Fb.kappa()) || (kappa_a == -Fb.kappa())) {
      dF.min_pt() = 0;
      dF.max_pt() = 0;
      return dF;
    }

    const auto sk = double(kappa_a - Fb.kappa());
    const auto ck = -(m_q / std::sqrt(2.0) / 3.0) * sk;

    Rab_rhs(-1, m_vec, &dF, Fb, ck);
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa())) {
      return 0.0;
    }

    const auto sk = double(Fa.kappa() - Fb.kappa());
    return -(m_q /(std::sqrt(2.0)*3.0)) * sk * Rab(-1, m_vec, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto m_q = std::abs(PhysConst::alpha * omega);
  }

private:
  double m_q;
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
V5_sigma_nr(const Grid &grid, int sigma) {

  switch (sigma) {

  case +1:
    return std::make_unique<E5_nr>(grid);

  case -1:
    return std::make_unique<L5_nr>(grid);

  case 0:
    return std::make_unique<M5_nr>(grid, 1.0e-4);

  default:
    return nullptr;
  }
}

} // namespace multipole

//==============================================================================
//==============================================================================

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_E_w_nr(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<E5_nr>(wf.grid(), k, omega);
  return std::make_unique<E_nr>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_M_w_nr(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for M1, =2 for M2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<M5_nr>(wf.grid(), k, omega);
  return std::make_unique<M_nr>(wf.grid(), k, omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_L_w_nr(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [1.0e-4]"},
               {"gamma5", "Use gamma^5 variant: true/false [false]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 1.0e-4);
  const auto gamma5 = input.get("gamma5", false);
  if (gamma5)
    return std::make_unique<L5_nr>(wf.grid(), k, omega);
  return std::make_unique<L_nr>(wf.grid(), k, omega);
}

} // namespace DiracOperator