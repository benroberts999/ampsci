#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//! @brief General function of r, even scalar operator
//! @details Pass only grid to just get r, or either pass a lambda/function
//! [f(r)], or a number, n, (for r^n).
class RadialF final : public ScalarOperator {
  // = some function of r
  // Pass only grid to just get r, or
  // either pass a lambda/function [f(r)], or a number, n, (for r^n)
  // Explicitely even with rank 0, so ...
private:
  std::vector<double> fillVec(const Grid &gr,
                              const std::function<double(double)> &f) {
    std::vector<double> f_r;
    f_r.reserve(gr.num_points());
    for (auto r : gr.r())
      f_r.push_back(f(r));
    return f_r;
  }

public:
  RadialF(const Grid &rgrid, const std::function<double(double)> &f)
      : ScalarOperator(Parity::even, 1.0, fillVec(rgrid, f)) {}
  RadialF(const Grid &rgrid, const double n)
      : ScalarOperator(Parity::even, 1.0, fillVec(rgrid, [n](double r) {
                         return std::pow(r, n);
                       })) {}
  std::string name() const override final { return "RadialFunction"; }
  std::string units() const override final { return "au"; }
};

//! radial derivative operator
class dr final : public ScalarOperator {

public:
  dr()
      : ScalarOperator(Parity::even, 1.0, {}, {1, 0, 0, 1}, 1, Realness::real) {
  }
  std::string name() const override final { return "dr"; }
  std::string units() const override final { return "au"; }
};

//------------------------------------------------------------------------------
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_r(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"power", "Power (real) for r^k"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  const auto power = input.get("power", 1.0);
  std::cout << "r^(" << power << ")\n";
  return std::make_unique<RadialF>(wf.grid(), power);
}

inline std::unique_ptr<DiracOperator::TensorOperator>
generate_dr(const IO::InputBlock &input, const Wavefunction &) {
  input.check({{"", "no input"}});
  if (input.has_option("help")) {
    return nullptr;
  }
  using namespace DiracOperator;
  return std::make_unique<dr>();
}

} // namespace DiracOperator
