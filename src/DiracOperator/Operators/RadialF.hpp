#pragma once
#include "DiracOperator/TensorOperator.hpp"

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

} // namespace DiracOperator
