#pragma once
#include "Dirac/DiracSpinor.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/Wigner_369j.hpp"
#include <algorithm>
#include <memory>
#include <string>
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

  // Everything that follows is overloading operators:
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
      std::cerr
          << "FAIL 50 in ScalarOperator_old. Cannot mix real and imaginary "
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
      std::cerr
          << "FAIL 50 in ScalarOperator_old. Cannot mix real and imaginary "
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
const DiracMatrix ig5(0, 1, -1, 0); /*?*/
} // namespace GammaMatrix

//******************************************************************************
enum class OperatorParity { even, odd };

class DiracOperator {
protected:
  DiracOperator(int k, OperatorParity pi, double c,
                const std::vector<double> &inv, int d_order = 0)
      : rank(k), parity(pi), constant(c), vec(inv), diff_order(d_order) //
        {};

public:
  virtual ~DiracOperator() = default;

private:
  const int rank;
  const OperatorParity parity;
  const double constant;
  const std::vector<double> vec; // useful to be able to update this!
  const int diff_order;

public:
  virtual double reducedME(const DiracSpinor &Fa,
                           const DiracSpinor &Fb) const = 0;

  bool isZero(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    return isZero(Fa.k, Fb.k);
  }
  bool isZero(const int ka, int kb) const {
    // checks rank and parity
    if (rank < std::abs(Wigner::twoj_k(ka) - Wigner::twoj_k(kb)) / 2)
      return true;
    if ((parity == OperatorParity::even) !=
        (Wigner::parity_k(ka) == Wigner::parity_k(kb)))
      return true;
    return false; /*may still be zero*/
  }

  const std::vector<double> &getv() const { return vec; }

public:
  //****************************************************************************
  double radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const
  // could implement this directly (as before), to avoid the an unnecisary
  // allocation of extra spinor! Woule be ~10% faster to implement seperately
  // with allocation
  {
    if (isZero(Fa, Fb))
      return 0.0;
    return Fa * radial_rhs(Fa.k, Fb);
  }

  // XXX Now, add "reduced_rhs" ? No, just "angular_factor"
  //****************************************************************************
  virtual DiracSpinor radial_rhs(const int k_lhs, const DiracSpinor &Fb) const
  // psi1 * dPsi2 = h.radialIntegral(psi1,psi2)
  // Because of angular factor, _may_ depend on kappa of 'lhs'
  {
    // Note: n and kappa from original psi, but not meaningful!
    const auto &gr = *(Fb.p_rgrid);
    DiracSpinor dPsi(Fb.n, Fb.k, gr);
    if (isZero(k_lhs, Fb.k))
      return dPsi;

    auto cff = angularCff(k_lhs, Fb.k);
    auto cgg = angularCgg(k_lhs, Fb.k);
    auto cfg = angularCfg(k_lhs, Fb.k);
    auto cgf = angularCgf(k_lhs, Fb.k);

    // Strangeness here to account for possible derivatives
    // (I am trying to avoid doing a copy when no derivative)
    // Copy is unavoidable when calcing derivative, but that's rare!
    const std::vector<double> *rhs_f = &(Fb.f);
    const std::vector<double> *rhs_g = &(Fb.g);
    std::unique_ptr<const std::vector<double>> dummy_df = nullptr;
    std::unique_ptr<const std::vector<double>> dummy_dg = nullptr;
    if (diff_order > 0) {
      // rhs_f is either a pointer to F_input, OR (in case of deriv, F')
      dummy_df = std::make_unique<std::vector<double>>(
          NumCalc::derivative(Fb.f, gr.drdu, gr.du, diff_order));
      dummy_dg = std::make_unique<std::vector<double>>(
          NumCalc::derivative(Fb.g, gr.drdu, gr.du, diff_order));
      rhs_f = dummy_df.get();
      rhs_g = dummy_dg.get();
    }

    for (unsigned i = 0; i < Fb.pinf; i++) {
      dPsi.f[i] = constant * (cff * (*rhs_f)[i] + cfg * (*rhs_g)[i]);
      dPsi.g[i] = constant * (cgf * (*rhs_f)[i] + cgg * (*rhs_g)[i]);
    }
    if (!vec.empty()) {
      dPsi *= vec;
    }

    return dPsi;
  }

  //****************************************************************************
  std::string rme_symbol(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    return std::string("<") + Fa.shortSymbol() + "|h|" + Fb.shortSymbol() + ">";
  }

protected:
  virtual double angularCff(int /*k_a*/, int /*k_b*/) const { return 1.0; }
  virtual double angularCgg(int, int) const { return 1.0; }
  virtual double angularCfg(int, int) const { return 0.0; }
  virtual double angularCgf(int, int) const { return 0.0; }
};

//------------------------------------------------------------------------------
class ScalarOperator : public DiracOperator {
public:
  ScalarOperator(OperatorParity pi, double in_coef,
                 const std::vector<double> &in_v,
                 const DiracMatrix &in_g = GammaMatrix::ident, int in_diff = 0)
      : DiracOperator(0, pi, in_coef, in_v, in_diff), c_ff(in_g.e00),
        c_fg(in_g.e01), c_gf(in_g.e10), c_gg(in_g.e11) {}

  ScalarOperator(const std::vector<double> &in_v)
      : DiracOperator(0, OperatorParity::even, 1.0, in_v, 0), c_ff(1.0),
        c_fg(0.0), c_gf(0.0), c_gg(1.0) {}

  ScalarOperator(double in_coef, const std::vector<double> &in_v)
      : DiracOperator(0, OperatorParity::even, in_coef, in_v, 0), c_ff(1.0),
        c_fg(0.0), c_gf(0.0), c_gg(1.0) {}

public:
  virtual double reducedME(const DiracSpinor &Fa,
                           const DiracSpinor &Fb) const override {
    auto inv_threej = std::sqrt(Fb.twoj() + 1.0);
    return inv_threej * matrixEl(Fa, Fb);
  }
  virtual double matrixEl(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    return radialIntegral(Fa, Fb);
  }

private:
  const double c_ff, c_fg, c_gf, c_gg;

protected:
  double virtual angularCff(int, int) const override { return c_ff; }
  double virtual angularCgg(int, int) const override { return c_gg; }
  double virtual angularCfg(int, int) const override { return c_fg; }
  double virtual angularCgf(int, int) const override { return c_gf; }
};

//------------------------------------------------------------------------------
class NullOperator : public ScalarOperator {
public:
  NullOperator() : ScalarOperator(OperatorParity::even, 0, {}) {}

  double reducedME(const DiracSpinor &, const DiracSpinor &) const override {
    return 0.0;
  }
  double matrixEl(const DiracSpinor &, const DiracSpinor &) const override {
    return 0.0;
  }
};
