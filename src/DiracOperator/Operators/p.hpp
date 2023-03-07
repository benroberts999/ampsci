#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"

namespace DiracOperator {

//==============================================================================

class p final : public TensorOperator {
public:
  p() : TensorOperator(1, Parity::odd, 1.0, {}, 0, Realness::imaginary) {}
  std::string name() const override final { return "p"; }
  std::string units() const override final { return "i au"; }

  double angularF(const int ka, const int kb) const override final {
    return -1.0 * Angular::Ck_kk(1, ka, kb);
  }

  double angularCff(int, int) const override final { return 1.0; }
  double angularCgg(int, int) const override final { return 1.0; }
  double angularCfg(int, int) const override final { return 0.0; }
  double angularCgf(int, int) const override final { return 0.0; }

  virtual DiracSpinor radial_rhs(const int kappa_a,
                                 const DiracSpinor &Fb) const override final {
    // {(f' + eta/r f) , (g' + xi/r g)}
    const auto kappa_b = Fb.kappa();
    DiracSpinor dF(0, kappa_a, Fb.grid_sptr());

    // const double eta = (kappa_a == kappa_b + 1) ? kappa_a : -kappa_b;
    // const double xi = (kappa_a == kappa_b + 1) ? kappa_b : -kappa_a;

    const auto la = Angular::l_k(kappa_a);
    const auto lb = Fb.l();
    const auto lta = Angular::l_k(-kappa_a);
    const auto ltb = Angular::l_k(-kappa_b);
    const double eta = (la == lb - 1) ? lb : -lb - 1;
    const double xi = (lta == ltb - 1) ? ltb : -ltb - 1;

    const auto &gr = Fb.grid();
    const auto df = NumCalc::derivative(Fb.f(), gr.drdu(), gr.du(), 1);
    const auto dg = NumCalc::derivative(Fb.g(), gr.drdu(), gr.du(), 1);

    using namespace qip::overloads;
    dF.f() = df + eta * (Fb.f() / gr.r());
    dF.g() = dg + xi * (Fb.g() / gr.r());
    dF.max_pt() = Fb.max_pt(); //?
    return dF;
  }
};

//==============================================================================
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_p(const IO::InputBlock &input, const Wavefunction &) {
  input.check({{"", "no input"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  return std::make_unique<DiracOperator::p>();
}

} // namespace DiracOperator