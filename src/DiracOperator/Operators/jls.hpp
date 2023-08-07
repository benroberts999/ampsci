#pragma once
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//==============================================================================
//! l (orbital angular momentum) operator
class l : public TensorOperator {
public:
  l() : TensorOperator(1, Parity::even) {}

  double angularF(const int ka, const int kb) const override final {
    const auto la = Angular::l_k(ka);
    const auto lb = Angular::l_k(kb);
    const auto tja = Angular::twoj_k(ka);
    const auto tjb = Angular::twoj_k(kb);
    if (la != lb)
      return 0.0;
    const auto sign = Angular::neg1pow_2(tjb + 2 * lb - 1);
    const auto fact =
        std::sqrt(double((tja + 1) * (tjb + 1) * (2 * la + 1) * la * (la + 1)));
    const auto sjs = Angular::sixj_2(tjb, tja, 2, 2 * lb, 2 * lb, 1);
    return sign * fact * sjs;
  }
  std::string name() const override { return std::string("l"); }
  std::string units() const override { return "au"; }
};

//==============================================================================
//! s (spin) operator
class s : public TensorOperator {
public:
  s() : TensorOperator(1, Parity::even, 1.0) {}

  double angularF(const int ka, const int kb) const override final {
    const auto la = Angular::l_k(ka);
    const auto lb = Angular::l_k(kb);
    const auto tja = Angular::twoj_k(ka);
    const auto tjb = Angular::twoj_k(kb);
    if (la != lb)
      return 0.0;
    const auto sign = Angular::neg1pow_2(tja + 2 * lb - 1);
    const auto fact = std::sqrt((tja + 1) * (tjb + 1) * 1.5);
    const auto sjs = Angular::sixj_2(tjb, tja, 2, 1, 1, 2 * la);
    return sign * fact * sjs;
  }
  std::string name() const override { return std::string("s"); }
  std::string units() const override { return "au"; }
};

//==============================================================================
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_l(const IO::InputBlock &input, const Wavefunction &) {
  using namespace DiracOperator;
  input.check({{"no options", ""}});
  if (input.has_option("help")) {
    return nullptr;
  }
  return std::make_unique<l>();
}

inline std::unique_ptr<DiracOperator::TensorOperator>
generate_s(const IO::InputBlock &input, const Wavefunction &) {
  using namespace DiracOperator;
  input.check({{"no options", ""}});
  if (input.has_option("help")) {
    return nullptr;
  }
  return std::make_unique<s>();
}

} // namespace DiracOperator