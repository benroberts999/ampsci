#pragma once
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Potentials/BSM_Vee_Potentials.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//==============================================================================

class Vee : public TensorOperator {
public:
  Vee(const std::vector<DiracSpinor> &core_in, const bool contact,
      const double mu_in, const std::string type, const bool eN)
    : TensorOperator(0, Parity::odd),
      m_core(core_in),
      m_contact(contact),
      m_mu(mu_in),
      m_type(type),
      m_eN(eN) {}

  // radialIntegral as defined here is actually the *full* integral
  // angularF needs to cancel out 'factor' in fullME (see TensorOperator.cpp)
  double angularF(const int ka, const int kb) const override final {

    // Inverse of factor in fullME
    const auto twoja = Angular::twoj_k(ka);
    const auto twojb = Angular::twoj_k(kb);
    const auto tma = std::min(twoja, twojb);

    const auto sign = Angular::neg1pow_2(twoja - tma);

    const auto threej = Angular::threej_2(twoja, 0, twojb, -tma, 0, tma);

    if (Angular::zeroQ(threej)) {
      return 0.0;
    } else {
      const auto inv_factor =
        1.0 / (sign * Angular::threej_2(twoja, 0, twojb, -tma, 0, tma));

      return inv_factor;
    }
  }

  std::string name() const override { return std::string("Vee"); }

  DiracSpinor radial_rhs(const int kappa_a,
                         const DiracSpinor &Fb) const override final {

    const double gghc = 1.0;

    return BSM_Vee::V_Fv(m_core, m_eN, Fb, m_type, kappa_a, gghc, m_contact,
                         m_mu);
  }

  // BSM_Vee::V_Fv is the *full* RHS, not just radial.
  double radialIntegral(const DiracSpinor &Fa,
                        const DiracSpinor &Fb) const override final {
    return Fa * radial_rhs(Fa.kappa(), Fb);
  }

  static std::unique_ptr<TensorOperator> generate(const IO::InputBlock &input,
                                                  const Wavefunction &wf) {
    input.check({{"contact", "contact interaction? [false]"},
                 {"mu", "mediator mass in mc/hbar [1]"}});
    if (input.has_option("help"))
      return nullptr;
    const bool contact = input.get("contact", false);
    const double mu = input.get("mu", 1.0);
    // const std::string type = input.get<std::string>("type", "sp");

    if (input.has_option("help")) {
      return nullptr;
    }
    return std::make_unique<Vee>(wf.core(), contact, mu, "sp", false);
  }

private:
  const std::vector<DiracSpinor> m_core;
  const bool m_contact;
  const double m_mu;
  const std::string m_type;
  const bool m_eN;

  //==============================================================================
};
} // namespace DiracOperator
