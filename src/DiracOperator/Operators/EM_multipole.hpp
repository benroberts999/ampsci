#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Maths.hpp"

namespace DiracOperator {

//==============================================================================
//! @brief Electric multipole operator, L-form, including frequency-dependence.
/*! @details
- Used, e.g., to check validity of dipole approximation for high frequencies.
- in terms of j(q*r), where q = alpha*omega
- If transition_form is true, will use the "t" (transition form); if false, will use the "moment" form. Nb: E1, E2 etc. is the "moment" form, but alpha exp(iqr) goes to transition form - see Eq. (6.129) in Johnson for factors!
*/
class Ek_omega final : public TensorOperator {
public:
  Ek_omega(const Grid &gr, int K, double alpha, double omega,
           bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 1.0,
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

    const auto c1 = double(kappa_a - Fb.kappa()) / (m_K + 1) + 1;
    const auto c2 = c1 - 2;
    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = j1[i] * Fb.f(i) - j2[i] * c1 * Fb.g(i);
      dF.g(i) = j1[i] * Fb.g(i) - j2[i] * c2 * Fb.g(i);
    }
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto cc = double(Fa.kappa() - Fb.kappa()) / (m_K + 1);
    // const auto c2 = c1 - 2;

    const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
    const auto pf = std::min(Fa.max_pt(), Fb.max_pt());

    const auto &drdu = Fb.grid().drdu();

    const auto ff = NumCalc::integrate(1.0, pi, pf, j1, Fa.f(), Fb.f(), drdu);
    const auto gg = NumCalc::integrate(1.0, pi, pf, j1, Fa.g(), Fb.g(), drdu);

    // fg only enters if states not the same. Compare memory address, not quantum numbers, may be spline/non-spline?
    const auto fg_not_gf = &Fa != &Fb;
    const auto fg =
        fg_not_gf ? NumCalc::integrate(1.0, pi, pf, j2, Fa.f(), Fb.g(), drdu) :
                    0.0;
    const auto gf =
        fg_not_gf ? NumCalc::integrate(1.0, pi, pf, j2, Fa.g(), Fb.f(), drdu) :
                    fg;

    return (ff + gg - cc * (fg + gf) - (fg - gf)) * Fb.grid().du();
  }

  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K);

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

    const auto c1 = double(kappa_a - Fb.kappa()) / (m_K + 1);
    const auto c2 = -double(m_K);
    const auto c3 = c1 * (m_K + 1);
    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = ((c3 + c2) * j1_on_qr[i] - c1 * j2[i]) * Fb.g(i);
      dF.g(i) = ((c3 - c2) * j1_on_qr[i] - c1 * j2[i]) * Fb.f(i);
    }
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto K = double(m_K);
    const auto a = double(Fa.kappa() - Fb.kappa()) / (K + 1.0);
    double cExp = std::sqrt(2.0 * K + 1.0);
    double cc = 1.0; // Johnson
    // double cc = std::sqrt((K + 1) / K); // Ben

    const auto Pp1 = Pab(+1, j1_on_qr, Fa, Fb);
    const auto Pp2 = Pab(+1, j2, Fa, Fb);
    const auto Pm1 = Pab(-1, j1_on_qr, Fa, Fb);

    return cc * cExp * (K * Pm1 - a * ((K + 1.0) * Pp1 - Pp2));

    // const auto c1 = double(Fa.kappa() - Fb.kappa()) / (m_K + 1);
    // const auto c2 = -double(m_K);
    // const auto c3 = c1 * (m_K + 1);

    // const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
    // const auto pf = std::min(Fa.max_pt(), Fb.max_pt());

    // const auto &drdu = Fb.grid().drdu();

    // const auto fg1 =
    //     NumCalc::integrate(1.0, pi, pf, j1_on_qr, Fa.f(), Fb.g(), drdu);
    // const auto fg2 = NumCalc::integrate(1.0, pi, pf, j2, Fa.f(), Fb.g(), drdu);

    // const auto gf1 =
    //     NumCalc::integrate(1.0, pi, pf, j1_on_qr, Fa.g(), Fb.f(), drdu);
    // const auto gf2 = NumCalc::integrate(1.0, pi, pf, j2, Fa.g(), Fb.f(), drdu);

    // return ((c3 + c2) * fg1 + (c3 - c2) * gf1 - c1 * (fg2 + gf2)) *
    //        Fb.grid().du();
  }

  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K);

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
  Lk_omega(const Grid &gr, int K, double alpha, double omega,
           bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K),
        m_transition_form(transition_form) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return m_transition_form ? std::string("tv^L_") + std::to_string(m_K) :
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

    const auto c0 = double(m_K);
    const auto c1 = double(kappa_a - Fb.kappa());
    const auto cc = std::sqrt(c0 / (2 * c0 + 1) / (c0 + 1));

    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = cc * ((c0 - c1) * j1_on_qr[i] - j2[i]) * Fb.g(i);
      dF.g(i) = cc * ((-c0 - c1) * j1_on_qr[i] + j2[i]) * Fb.f(i);
    }
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto c0 = double(m_K);
    const auto c1 = double(Fa.kappa() - Fb.kappa());
    const auto cc = std::sqrt(c0 / (2 * c0 + 1) / (c0 + 1));

    const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
    const auto pf = std::min(Fa.max_pt(), Fb.max_pt());

    const auto &drdu = Fb.grid().drdu();

    const auto same = &Fa == &Fb;

    const auto fg1 =
        NumCalc::integrate(1.0, pi, pf, j1_on_qr, Fa.f(), Fb.g(), drdu);

    const auto fg2 = NumCalc::integrate(1.0, pi, pf, j2, Fa.f(), Fb.g(), drdu);

    const auto gf1 =
        same ? fg1 :
               NumCalc::integrate(1.0, pi, pf, j1_on_qr, Fa.g(), Fb.f(), drdu);
    const auto gf2 =
        same ? fg2 : NumCalc::integrate(1.0, pi, pf, j2, Fa.g(), Fb.f(), drdu);

    return cc * ((c0 - c1) * fg1 + (-c0 - c1) * gf1 - fg2 + gf2) *
           Fb.grid().du();
  }

  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K);

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

    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = j1[i] * Fb.g(i);
      dF.g(i) = j1[i] * Fb.f(i);
    }
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa()) || (Fa.kappa() == -Fb.kappa())) {
      return 0.0;
    }

    const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
    const auto pf = std::min(Fa.max_pt(), Fb.max_pt());

    const auto &drdu = Fb.grid().drdu();

    const auto fg = NumCalc::integrate(1.0, pi, pf, j1, Fa.f(), Fb.g(), drdu);
    const auto gf = NumCalc::integrate(1.0, pi, pf, j1, Fa.g(), Fb.f(), drdu);

    return (fg + gf) * Fb.grid().du();
  }

  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     (2.0 / m_alpha) * qip::double_factorial(2 * m_K + 1) /
                         qip::pow(q, m_K);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &j1);
  }

private:
  double m_alpha;
  int m_K;
  bool m_transition_form = true;
  std::vector<double> j1{};
};

//==============================================================================
//! @brief Scalar (temporal) multipole operator, V-form, including frequency-dependence.
class Vk_omega final : public TensorOperator {
public:
  Vk_omega(const Grid &gr, int K, double alpha, double omega,
           bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::even : Parity::odd, 0.0,
                       gr.r(), 0, Realness::real, true),
        m_alpha(alpha),
        m_K(K),
        m_transition_form(transition_form) {
    updateFrequency(omega);
  }
  std::string name() const override final {
    return m_transition_form ? std::string("tv^V_") + std::to_string(m_K) :
                               std::string("Vv(") + std::to_string(m_K) + ")";
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

    const auto c0 = double(m_K);
    const auto cc = std::sqrt(c0 / (2 * c0 + 1) / (c0 + 1));

    for (auto i = Fb.min_pt(); i < Fb.max_pt(); i++) {
      dF.f(i) = cc * jk[i] * Fb.f(i);
      dF.g(i) = cc * jk[i] * Fb.g(i);
    }
    return dF;
  }

  //--------------
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {

    if (isZero(Fa.kappa(), Fb.kappa())) {
      return 0.0;
    }

    const auto c0 = double(m_K);
    const auto cc = std::sqrt(c0 / (2 * c0 + 1) / (c0 + 1));

    const auto pi = std::max(Fa.min_pt(), Fb.min_pt());
    const auto pf = std::min(Fa.max_pt(), Fb.max_pt());

    const auto &drdu = Fb.grid().drdu();

    const auto ff = NumCalc::integrate(1.0, pi, pf, jk, Fa.f(), Fb.f(), drdu);
    const auto gg = NumCalc::integrate(1.0, pi, pf, jk, Fa.g(), Fb.g(), drdu);

    return cc * (ff + gg) * Fb.grid().du();
  }

  void update_q(double q) { updateFrequency(q / m_alpha); }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    m_q = q;

    // should be 1 for "transition" operator, otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     qip::double_factorial(2 * m_K + 1) / qip::pow(q, m_K);

    SphericalBessel::fillBesselVec_kr(m_K, q, m_vec, &jk);
  }

private:
  double m_alpha; // (including var-alpha)
  int m_K;
  double m_q{0.0};
  bool m_transition_form = true;
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
generate_Mk_omega(const IO::InputBlock &input, const Wavefunction &wf) {
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
  return std::make_unique<Mk_omega>(wf.grid(), k, wf.alpha(), omega,
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

} // namespace DiracOperator
