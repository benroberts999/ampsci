#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {

//==============================================================================
//! @brief Electric multipole operator, alpha.t^1, L-form.
/*! @details
- Used, e.g., to check validity of dipole approximation for high frequencies.
- in terms of j(q*r), where q = alpha*omega
- If transition_form is true, will use the "t" (transition form); if false, will use the "moment" form. Nb: E1, E2 etc. is the "moment" form, but alpha exp(iqr) goes to transition form - see Eq. (6.129) in Johnson for factors!
- Note: tk (transition) version, different convention to Johnson. Moment version the same!
*/
class Ek_omega final : public TensorOperator {
public:
  Ek_omega(const Grid &gr, int K, double alpha, double omega,
           bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K),
        m_transition_form(transition_form) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return m_transition_form ? std::string("t^E_") + std::to_string(m_K) :
                               std::string("E(") + std::to_string(m_K) + ")";
  }
  std::string units() const override final {
    return m_transition_form ? std::string("") : std::string("aB^k");
  }

  double angularF(const int ka, const int kb) const override final {
    return m_constant * Angular::Ck_kk(m_K, ka, kb);
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

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto cc = double(Fa.kappa() - Fb.kappa()) / (K + 1.0);
    const auto cx = std::sqrt((K + 1.0) / K);

    return cx * (Rab(+1, j1, Fa, Fb) - (cc + 1.0) * Vab(j2, Fa, Fb) +
                 (1.0 - cc) * Wab(j2, Fa, Fb));
  }

  //! nb: q = alpha*omega!
  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K) *
                         std::sqrt(m_K / (m_K + 1.0));

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
  }

private:
  double m_alpha;
  int m_K;
  bool m_transition_form = true;
  std::vector<double> j1{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Electric multipole operator, V-form, including frequency-dependence.
/*! @details
- Used, e.g., to check validity of dipole approximation for high frequencies.
- in terms of j(q*r), where q = alpha*omega
- If transition_form is true, will use the "t" (transition form); if false, will use the "moment" form. Nb: E1, E2 etc. is the "moment" form, but alpha exp(iqr) goes to transition form - see Eq. (6.129) in Johnson for factors!
- Note: tk (transition) version, different convention to Johnson. Moment version the same!
*/
class Ekv_omega final : public TensorOperator {
public:
  Ekv_omega(const Grid &gr, int K, double alpha, double omega,
            bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K),
        m_transition_form(transition_form) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return m_transition_form ? std::string("tv^E_") + std::to_string(m_K) :
                               std::string("Ev(") + std::to_string(m_K) + ")";
  }
  std::string units() const override final {
    return m_transition_form ? std::string("") : std::string("aB^k");
  }

  double angularF(const int ka, const int kb) const override final {
    return m_constant * Angular::Ck_kk(m_K, ka, kb);
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

    if (isZero(Fa.kappa(), Fb.kappa())) {
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
  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K) *
                         std::sqrt(m_K / (m_K + 1.0));

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      j1_on_qr[i] /= (q * m_vec[i]);
    }
  }

private:
  double m_alpha; // (including var-alpha)
  int m_K;
  double m_q{0.0};
  bool m_transition_form = true;
  std::vector<double> j1_on_qr{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Longitudanal multipole operator, V-form, including frequency-dependence.
class Lk_omega final : public TensorOperator {
public:
  Lk_omega(const Grid &gr, int K, double alpha, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^L_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return m_constant * Angular::Ck_kk(m_K, ka, kb);
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
  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = 1.0;

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1_on_qr);
    SphericalBessel::fillBesselVec_kr(m_K + 1, q, m_vec, &j2);
    for (std::size_t i = 0; i < m_vec.size(); ++i) {
      j1_on_qr[i] /= (q * m_vec[i]);
    }
  }

private:
  double m_alpha; // (including var-alpha)
  int m_K;
  double m_q{0.0};
  std::vector<double> j1_on_qr{};
  std::vector<double> j2{};
};

//==============================================================================
//! @brief Magnetic multipole operator, including frequency-dependence.
class Mk_omega final : public TensorOperator {
public:
  Mk_omega(const Grid &gr, int K, double alpha, double omega,
           bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K),
        m_transition_form(transition_form) {
    updateFrequency(omega);
  }

  std::string name() const override final {
    return m_transition_form ? std::string("t^M_") + std::to_string(m_K) :
                               std::string("M(") + std::to_string(m_K) + ")";
  }

  std::string units() const override final {
    return m_transition_form ? std::string("") : std::string("mu_B*aB^k");
  }

  double angularF(const int ka, const int kb) const override final {
    return -m_constant * (ka + kb) * Angular::Ck_kk(m_K, -ka, kb) / (m_K + 1);
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

    if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto sk = double(Fa.kappa() + Fb.kappa());
    assert(m_K != 0); // should already be discounted!
    const auto ck = sk / std::sqrt(K * (K + 1.0));

    return -ck * Pab(+1, j1, Fa, Fb);
  }

  //! nb: q = alpha*omega!
  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     (1.0 / PhysConst::muB_CGS) *
                         qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K) *
                         std::sqrt(m_K / (m_K + 1.0));

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1);
  }

private:
  double m_alpha;
  int m_K;
  bool m_transition_form = true;
  std::vector<double> j1{};
};

//==============================================================================
//! @brief Temporal component of vector multipole operator, V-form, including frequency-dependence.
class Vk_omega final : public TensorOperator {
public:
  Vk_omega(const Grid &gr, int K, double alpha, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^V_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return m_constant * Angular::Ck_kk(m_K, ka, kb);
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
  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;
    m_constant = 1.0;

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  double m_alpha; // (including var-alpha)
  int m_K;
  double m_q{0.0};
  std::vector<double> jk{};
};

//==============================================================================
//! @brief Scalar multipole operator, e^{iqr}gamma^0, including frequency-dependence.
class Sk_omega final : public TensorOperator {
public:
  Sk_omega(const Grid &gr, int K, double alpha, double omega)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return std::string("tv^S_") + std::to_string(m_K);
  }
  std::string units() const override final { return std::string(""); }

  double angularF(const int ka, const int kb) const override final {
    return m_constant * Angular::Ck_kk(m_K, ka, kb);
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
  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;
    m_constant = 1.0;

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  double m_alpha; // (including var-alpha)
  int m_K;
  double m_q{0.0};
  std::vector<double> jk{};
};

//==============================================================================
//==============================================================================

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ek_omega(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [0]"},
               {"transition_form",
                "Use transition form (true), or moment form (fasle). nb: use "
                "moment form to compare to E1/E2 etc, transition form for "
                "alpha*e^{iqr} expansion [true]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  const auto transition_form = input.get("transition_form", true);
  return std::make_unique<Ek_omega>(wf.grid(), k, wf.alpha(), omega,
                                    transition_form);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Ekv_omega(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for E1, =2 for E2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [0.001]"},
               {"transition_form",
                "Use transition form (true), or moment form (fasle). nb: use "
                "moment form to compare to E1/E2 etc, transition form for "
                "alpha*e^{iqr} expansion [true]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.001);
  const auto transition_form = input.get("transition_form", true);
  return std::make_unique<Ekv_omega>(wf.grid(), k, wf.alpha(), omega,
                                     transition_form);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Mk_omega(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"k", "Rank: k=1 for M1, =2 for M2 etc. [1]"},
               {"omega", "Frequency: nb: q = alpha*omega [0]"},
               {"transition_form",
                "Use transition form (true), or moment form (fasle). nb: use "
                "moment form to compare to E1/E2 etc, transition form for "
                "alpha*e^{iqr} expansion [true]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  const auto transition_form = input.get("transition_form", true);
  return std::make_unique<Mk_omega>(wf.grid(), k, wf.alpha(), omega,
                                    transition_form);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Lk_omega(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"k", "Rank [1]"}, {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Lk_omega>(wf.grid(), k, wf.alpha(), omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Vk_omega(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"k", "Rank [1]"}, {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Vk_omega>(wf.grid(), k, wf.alpha(), omega);
}

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Sk_omega(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"k", "Rank [1]"}, {"omega", "Frequency: nb: q = alpha*omega [0]"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto k = input.get("k", 1);
  const auto omega = input.get("omega", 0.0);
  return std::make_unique<Sk_omega>(wf.grid(), k, wf.alpha(), omega);
}

} // namespace DiracOperator
