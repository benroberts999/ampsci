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
  double m_alpha; // (including var-alpha)
  int m_K;
  bool m_transition_form = true;
  std::vector<double> j1{};
  std::vector<double> j2{};
};

//==============================================================================
//! Magnetic multipole transition operator
class Mk_omega final : public TensorOperator {
public:
  Mk_omega(const Grid &gr, int K, double alpha, double omega,
           bool transition_form)
      : TensorOperator(K, Angular::evenQ(K) ? Parity::odd : Parity::even, +0.0,
                       {}, 0, Realness::real, true),
        m_r(gr.r()),
        m_K(K),
        m_alpha(alpha),
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
    return -(ka + kb) * Angular::Ck_kk(m_K, -ka, kb) / (m_K + 1);
  }
  double angularCff(int, int) const override final { return 0.0; }
  double angularCgg(int, int) const override final { return 0.0; }
  double angularCfg(int, int) const override final { return 1.0; }
  double angularCgf(int, int) const override final { return 1.0; }

  //! nb: q = alpha*omega!
  void updateFrequency(const double omega) override final {
    const auto q = std::abs(m_alpha * omega);
    // should be 1 for "transition" operator,
    // otherwise factor for the "moment" form
    m_constant = m_transition_form ?
                     1.0 :
                     (2.0 / m_alpha) * qip::double_factorial(2 * m_K + 1) /
                         qip::pow(q, m_K);
    SphericalBessel::fillBesselVec_kr(m_K, q, m_r, &m_vec);
  }

private:
  std::vector<double> m_r; // store radial vector (copy)
  int m_K;
  double m_alpha;
  bool m_transition_form = true;
};

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

} // namespace DiracOperator
