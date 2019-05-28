#pragma once
#include "DiracSpinor.h"
#include <vector>

//******************************************************************************
struct DiracMatrix {
  // Only really makes sense for g0 and g5 - '1' really means identity
  // Same struct is used for Pauli spin matrices too
  DiracMatrix(int in00 = 1, int in01 = 0, int in10 = 0, int in11 = 1,
              bool in_imag = false)
      : e00(in00), e01(in01), e10(in10), e11(in11), imaginary(in_imag) {}
  const int e00, e01, e10, e11;
  const bool imaginary;
};

//******************************************************************************
namespace PauliSpinMatrix {
const DiracMatrix ident(1, 0, 0, 1);
const DiracMatrix sx(0, 1, 1, 0);
const DiracMatrix sy(0, -1, 1, 0, true);
const DiracMatrix sz(1, 0, 0, -1);
} // namespace PauliSpinMatrix

namespace GammaMatrix {
const DiracMatrix ident(1, 0, 0, 1);
const DiracMatrix g0(1, 0, 0, -1);
const DiracMatrix g5(0, 1, 1, 0);
} // namespace GammaMatrix

//******************************************************************************
class DiracOperator {

public:
  DiracOperator(std::vector<double> in_v, DiracMatrix in_g = GammaMatrix::ident,
                int in_diff = 0, bool in_imag = false)
      : v(in_v), g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  DiracOperator(DiracMatrix in_g = GammaMatrix::ident, int in_diff = 0,
                bool in_imag = false)
      : g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  DiracOperator(int in_diff = 0, bool in_imag = false)
      : diff_order(in_diff), imaginary(in_imag) {}

  DiracOperator(DiracMatrix in_g = GammaMatrix::ident, bool in_imag = false)
      : g(in_g), imaginary(in_imag) {}

private: // Data
  const std::vector<double> v;
  const DiracMatrix g = GammaMatrix::ident;
  const int diff_order = 0;
  const bool imaginary = false;

public:
  DiracSpinor operate(DiracSpinor phi);
};
