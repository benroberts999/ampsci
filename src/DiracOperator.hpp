#pragma once
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <vector>

//******************************************************************************
struct DiracMatrix
// Two-component form Dirac matrices
// Only really makes sense for g0 and g5 - '1' really means identity
// Same struct is used for Pauli spin matrices too
// Notation for elements:
//  (e00  e01)
//  (e10  e11)
// Has a multiplication and addition operators.
// Note: (for now), no way to add/mix real and imaginary matrices!
// (This never happens anyway)
{
  DiracMatrix(int in00 = 1, int in01 = 0, int in10 = 0, int in11 = 1,
              bool in_imag = false)
      : e00(in00), e01(in01), e10(in10), e11(in11), imaginary(in_imag) {}
  const int e00, e01, e10, e11;
  const bool imaginary;

  void print() const {
    // this is really just for testing.
    auto i = (this->imaginary) ? "i" : " ";
    printf(" (%2i %2i)\n", this->e00, this->e01);
    printf("%s(%2i %2i)\n", i, this->e10, this->e11);
  }

  DiracMatrix operator*(const DiracMatrix &other) const {
    bool imag = false;
    int sign = 1;
    if (this->imaginary && other.imaginary) {
      sign = -1;
    } else if (this->imaginary || other.imaginary) {
      imag = true;
    }
    int a = sign * (this->e00 * other.e00 + this->e01 * other.e10);
    int b = sign * (this->e00 * other.e01 + this->e01 * other.e11);
    int c = sign * (this->e10 * other.e00 + this->e11 * other.e10);
    int d = sign * (this->e10 * other.e01 + this->e11 * other.e11);
    return DiracMatrix(a, b, c, d, imag);
  }

  DiracMatrix operator+(const DiracMatrix &other) const {
    bool imag = false;
    if (this->imaginary && other.imaginary) {
      imag = true;
    } else if (this->imaginary || other.imaginary) {
      std::cerr << "FAIL 50 in DiracOperator. Cannot mix real and imaginary "
                   "Dirac Matrices! \n";
      this->print();
      std::cerr << "+\n";
      other.print();
      std::abort();
    }
    int a = this->e00 + other.e00;
    int b = this->e01 + other.e01;
    int c = this->e10 + other.e10;
    int d = this->e11 + other.e11;
    return DiracMatrix(a, b, c, d, imag);
  }
  DiracMatrix operator-(const DiracMatrix &other) const {
    bool imag = false;
    if (this->imaginary && other.imaginary) {
      imag = true;
    } else if (this->imaginary || other.imaginary) {
      std::cerr << "FAIL 50 in DiracOperator. Cannot mix real and imaginary "
                   "Dirac Matrices! \n";
      this->print();
      std::cerr << "-\n";
      other.print();
      std::abort();
    }
    int a = this->e00 - other.e00;
    int b = this->e01 - other.e01;
    int c = this->e10 - other.e10;
    int d = this->e11 - other.e11;
    return DiracMatrix(a, b, c, d, imag);
  }
};

//******************************************************************************
// Make some 'global' constants - common dirac matrices
// Assumes DIRAC BASIS
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
  // XXX at the moment, two ways can be imaginary..
  // it's own, and from Dirac matrix....... OK?? XXX

public: // Constructors
  DiracOperator(double in_coef, const std::vector<double> &in_v,
                const DiracMatrix &in_g = GammaMatrix::ident, int in_diff = 0,
                bool in_imag = false)
      : coef(in_coef), v(in_v), g(in_g), diff_order(in_diff),
        imaginary(in_imag) {}

  DiracOperator(const std::vector<double> &in_v,
                const DiracMatrix &in_g = GammaMatrix::ident, int in_diff = 0,
                bool in_imag = false)
      : v(in_v), g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  DiracOperator(DiracMatrix in_g = GammaMatrix::ident, int in_diff = 0,
                bool in_imag = false)
      : g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  DiracOperator(int in_diff = 0, bool in_imag = false)
      : diff_order(in_diff), imaginary(in_imag) {}

  DiracOperator(DiracMatrix in_g = GammaMatrix::ident, bool in_imag = false)
      : g(in_g), imaginary(in_imag) {}

private: // Data
  double coef = 1;
  const std::vector<double> v;
  const DiracMatrix g = GammaMatrix::ident;
  const int diff_order = 0;
  const bool imaginary = false;
  // const Grid *const rgrid = nullptr;

public: // Methods
  DiracSpinor operate(const DiracSpinor &phi) const;
  DiracOperator HermetianConjugate() const;

public: // Operator overloads
  DiracSpinor operator*(const DiracSpinor &phi) const { return operate(phi); }
};

//******************************************************************************
inline DiracSpinor DiracOperator::operate(const DiracSpinor &phi) const {

  // Note: matrix must be either diagonal or off-diagonal
  // This isn't checked or enforced yet!??!?
  bool off_diag = (g.e01 != 0 || g.e10 != 0) ? true : false;

  auto g_imag = off_diag ? !phi.imaginary_g : phi.imaginary_g;

  DiracSpinor dPhi(phi.n, phi.k, *phi.p_rgrid, g_imag);
  dPhi.pinf = phi.pinf; //?
  dPhi.en = phi.en;     //?

  // if imag, sign of "g" changes (unless f was img, then sign of f changes)
  const auto g_sign = (imaginary && phi.imaginary_g) ? -1 : 1;
  const auto f_sign = (imaginary && !phi.imaginary_g) ? -1 : 1;

  if (off_diag) {
    for (std::size_t i = 0; i < phi.p_rgrid->ngp; i++) {
      dPhi.f[i] = g_sign * g.e01 * phi.g[i];
      dPhi.g[i] = f_sign * g.e10 * phi.f[i];
    }
  } else {
    for (std::size_t i = 0; i < phi.p_rgrid->ngp; i++) {
      dPhi.f[i] = f_sign * g.e00 * phi.f[i];
      dPhi.g[i] = g_sign * g.e11 * phi.g[i];
    }
  }

  // Differentiate:
  if (diff_order > 0) {
    dPhi.f = NumCalc::derivative(dPhi.f, phi.p_rgrid->drdu, phi.p_rgrid->du,
                                 diff_order);
    dPhi.g = NumCalc::derivative(dPhi.g, phi.p_rgrid->drdu, phi.p_rgrid->du,
                                 diff_order);
  }

  // Multiply by radial vector:
  if (v.size() > 0) {
    for (std::size_t i = 0; i < phi.p_rgrid->ngp; i++) {
      dPhi.f[i] *= v[i];
      dPhi.g[i] *= v[i];
    }
  }
  if (coef != 1) {
    for (std::size_t i = 0; i < phi.p_rgrid->ngp; i++) {
      dPhi.f[i] *= coef;
      dPhi.g[i] *= coef;
    }
  }
  return dPhi;
}

//******************************************************************************
inline DiracOperator DiracOperator::HermetianConjugate() const {
  // Transpose the Dirac Matrix:
  const auto gT = DiracMatrix(g.e00, g.e10, g.e01, g.e11, g.imaginary);
  // complex conjugate:
  const auto sign = (imaginary) ? -1 : 1;
  return DiracOperator(sign * coef, v, gT, diff_order, imaginary);
}
