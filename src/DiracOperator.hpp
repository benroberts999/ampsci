#pragma once
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "Physics/Wigner_369j.hpp"
#include <algorithm>
// #include <functional>
#include <memory>
#include <vector>
//
#include <iostream>

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
} // namespace GammaMatrix

//******************************************************************************
class ScalarOperator_old {
  // XXX at the moment, two ways can be imaginary..
  // it's own, and from Dirac matrix....... OK?? XXX
  // Various contructors: (all after v are optional)
  // ScalarOperator_old(C, v, DiracMatrix(a, b, c, d), diff_order, imag?)
  // ScalarOperator_old(C, DiracMatrix(a, b, c, d), diff_order, imag?)
  // ScalarOperator_old(v, DiracMatrix(a, b, c, d), diff_order, imag?)
  // ScalarOperator_old(DiracMatrix(a, b, c, d), diff_order, imag?)
  // ScalarOperator_old(diff_order, imag?)
  // ScalarOperator_old(DiracMatrix(a, b, c, d), imag?)

public: // Constructors
  ScalarOperator_old(double in_coef, const std::vector<double> &in_v,
                     const DiracMatrix &in_g = GammaMatrix::ident,
                     int in_diff = 0, bool in_imag = false)
      : coef(in_coef), v(in_v), g(in_g), diff_order(in_diff),
        imaginary(in_imag) {}

  ScalarOperator_old(double in_coef,
                     const DiracMatrix &in_g = GammaMatrix::ident,
                     int in_diff = 0, bool in_imag = false)
      : coef(in_coef), g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  ScalarOperator_old(const std::vector<double> &in_v,
                     const DiracMatrix &in_g = GammaMatrix::ident,
                     int in_diff = 0, bool in_imag = false)
      : v(in_v), g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  ScalarOperator_old(DiracMatrix in_g = GammaMatrix::ident, int in_diff = 0,
                     bool in_imag = false)
      : g(in_g), diff_order(in_diff), imaginary(in_imag) {}

  ScalarOperator_old(int in_diff = 0, bool in_imag = false)
      : diff_order(in_diff), imaginary(in_imag) {}

  ScalarOperator_old(DiracMatrix in_g = GammaMatrix::ident,
                     bool in_imag = false)
      : g(in_g), imaginary(in_imag) {}

public: // Data
  const double coef = 1;
  const std::vector<double> v;
  const DiracMatrix g = GammaMatrix::ident;
  const int diff_order = 0;
  const bool imaginary = false;
  // const rank = 0;   // max delta_j
  // const parity = 1; // gives allowed delta_l (??)

public: // Methods
  DiracSpinor operate(const DiracSpinor &phi) const;
  ScalarOperator_old HermetianConjugate() const;

public: // Operator overloads
  DiracSpinor operator*(const DiracSpinor &phi) const { return operate(phi); }
};

//******************************************************************************
inline DiracSpinor ScalarOperator_old::operate(const DiracSpinor &phi) const {

  // Note: matrix must be either diagonal or off-diagonal
  // This isn't checked or enforced yet!??!?
  bool off_diag = (g.e01 != 0 || g.e10 != 0) ? true : false;

  // Diagonal operator swaps f,g, so swaps which comp is imagingary
  auto g_imag = off_diag ? !phi.imaginary_g : phi.imaginary_g;
  // Also: if operator is itself imaginary, also swaps!
  if (imaginary) {
    g_imag = !g_imag;
  }

  // if imag, sign of "g" changes (unless f was img, then sign of f changes)
  const auto g_sign = (imaginary && phi.imaginary_g) ? -1 : 1;
  const auto f_sign = (imaginary && !phi.imaginary_g) ? -1 : 1;

  DiracSpinor dPhi(phi.n, phi.k, *phi.p_rgrid, g_imag);
  dPhi.pinf = phi.pinf; //?
  dPhi.en = phi.en;     //?

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
  if (!v.empty()) {
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
inline ScalarOperator_old ScalarOperator_old::HermetianConjugate() const {
  // Transpose the Dirac Matrix:
  const auto gT = DiracMatrix(g.e00, g.e10, g.e01, g.e11, g.imaginary);
  // complex conjugate:
  const auto sign = (imaginary) ? -1 : 1;
  return ScalarOperator_old(sign * coef, v, gT, diff_order, imaginary);
}

// static inline auto Rint()
//     -> std::function<double(const std::vector<double> &pa,
//                             const std::vector<double> &pb)> {
//   return [&](const std::vector<double> &pa, const std::vector<double> &pb) {
//     return NumCalc::integrate(pa, pb, gr.drdu, 1.0, 0, irmax);
//   };
// }

//******************************************************************************
enum class OperatorParity { even, odd };

class DiracOperator {
protected:
  DiracOperator(int k, OperatorParity pi, double c,
                const std::vector<double> &inv, int d_order)
      : rank(k), parity(pi), constant(c), vec(inv), diff_order(d_order) //
        {};

private:
  const int rank;
  const OperatorParity parity;
  const double constant;
  const std::vector<double> vec; // useful to be able to update this!
  const int diff_order;
  // private:
  //   const Grid *const p_rgrid; //??

public:
  virtual double reducedME(const DiracSpinor &Fa,
                           const DiracSpinor &Fb) const = 0;
  // virtual double matrixEl(const DiracSpinor &Fa,
  //                         const DiracSpinor &Fb) const = 0;

protected:
  double radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    if (rank < std::abs(Fa.twoj() - Fb.twoj()) / 2)
      return 0.0;
    if ((parity == OperatorParity::even) != (Fa.parity() == Fb.parity()))
      return 0.0;

    auto ka = Fa.k;
    auto kb = Fb.k;
    auto cff = angularCff(ka, kb);
    auto cgg = angularCgg(ka, kb);
    auto cfg = angularCfg(ka, kb);
    auto cgf = angularCgf(ka, kb);
    const auto &gr = *(Fa.p_rgrid);
    const auto irmax = std::min(Fa.pinf, Fb.pinf); // check if slow?

    const std::vector<double> *rhs_f = &(Fb.f);
    const std::vector<double> *rhs_g = &(Fb.g);

    // XXX - what about derivative ?? pointers!
    std::unique_ptr<const std::vector<double>> new_rhs_f = nullptr;
    std::unique_ptr<const std::vector<double>> new_rhs_g = nullptr;
    if (diff_order > 0) {
      new_rhs_f = std::make_unique<std::vector<double>>(
          NumCalc::derivative(Fb.f, gr.drdu, gr.du, diff_order));
      new_rhs_g = std::make_unique<std::vector<double>>(
          NumCalc::derivative(Fb.g, gr.drdu, gr.du, diff_order));
      rhs_f = &(*new_rhs_f); // don't cry
      rhs_g = &(*new_rhs_f);
    }

    // XXX dodgy hack..
    if (!vec.empty()) {
      auto Rff = (cff == 0.0) ? 0.0
                              : NumCalc::integrate(Fa.f, *rhs_f, vec, gr.drdu,
                                                   1.0, 0, irmax);
      auto Rgg = (cgg == 0.0) ? 0.0
                              : NumCalc::integrate(Fa.g, *rhs_g, vec, gr.drdu,
                                                   1.0, 0, irmax);
      auto Rfg = (cfg == 0.0) ? 0.0
                              : NumCalc::integrate(Fa.f, *rhs_g, vec, gr.drdu,
                                                   1.0, 0, irmax);
      auto Rgf = (cgf == 0.0) ? 0.0
                              : NumCalc::integrate(Fa.g, *rhs_f, vec, gr.drdu,
                                                   1.0, 0, irmax);

      return constant * (cff * Rff + cgg * Rgg + cfg * Rfg + cgf * Rgf) * gr.du;
    } else {
      auto Rff = (cff == 0.0)
                     ? 0.0
                     : NumCalc::integrate(Fa.f, *rhs_f, gr.drdu, 1.0, 0, irmax);
      auto Rgg = (cgg == 0.0)
                     ? 0.0
                     : NumCalc::integrate(Fa.g, *rhs_g, gr.drdu, 1.0, 0, irmax);
      auto Rfg = (cfg == 0.0)
                     ? 0.0
                     : NumCalc::integrate(Fa.f, *rhs_g, gr.drdu, 1.0, 0, irmax);
      auto Rgf = (cgf == 0.0)
                     ? 0.0
                     : NumCalc::integrate(Fa.g, *rhs_f, gr.drdu, 1.0, 0, irmax);

      return constant * (cff * Rff + cgg * Rgg + cfg * Rfg + cgf * Rgf) * gr.du;
    }
  }

protected:
  double virtual angularCff(int /*k_a*/, int /*k_b*/) const { return 1.0; }
  double virtual angularCgg(int, int) const { return 1.0; }
  double virtual angularCfg(int, int) const { return 0.0; }
  double virtual angularCgf(int, int) const { return 0.0; }
  // XXX Add lookup table for 3j, 6j etc?
};

//------------------------------------------------------------------------------
class ScalarOperator : public DiracOperator {
public:
  ScalarOperator(OperatorParity pi, double in_coef,
                 const std::vector<double> &in_v,
                 const DiracMatrix &in_g = GammaMatrix::ident, int in_diff = 0)
      : DiracOperator(0, pi, in_coef, in_v, in_diff), //
        c_ff(in_g.e00), c_fg(in_g.e01), c_gf(in_g.e10), c_gg(in_g.e11) {}

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
  // const bool imaginary; // XXX ok??
  const double c_ff, c_fg, c_gf, c_gg;

protected:
  double virtual angularCff(int, int) const override { return c_ff; }
  double virtual angularCgg(int, int) const override { return c_gg; }
  double virtual angularCfg(int, int) const override { return c_fg; }
  double virtual angularCgf(int, int) const override { return c_gf; }
};

//******************************************************************************
class DirectHamiltonian : public ScalarOperator // need say Dirac here?
// H_D = [V(r)         -c(dr - k/r)]
//       [c(dr + k/r)   V(r) - 2c^2]
//     = V(r) + c g5 k/r + c d_r (0,-1,1,0) + c^2 (0,0,0,-2)
//     = sum of 4 scalar operators
// nb: V = v_nuc + v_dir
{
public:
  DirectHamiltonian(
      const Grid &gr,
      const std::initializer_list<const std::vector<double> *const> vs,
      double alpha)                                    //
      : ScalarOperator(OperatorParity::even, 1.0, {}), //
        cl(1.0 / alpha),                               //
        v(std::make_unique<ScalarOperator>(
            ScalarOperator(OperatorParity::even, 1.0, NumCalc::sumVecs(vs),
                           GammaMatrix::ident))),
        cg5or(std::make_unique<ScalarOperator>(
            ScalarOperator(OperatorParity::even, cl, gr.inverse_r(),
                           DiracMatrix(0, 1, 1, 0), 0))),
        cdr(std::make_unique<ScalarOperator>(ScalarOperator(
            OperatorParity::even, cl, {}, DiracMatrix(0, -1, 1, 0), 1))),
        c2(std::make_unique<ScalarOperator>(ScalarOperator(
            OperatorParity::even, cl * cl, {}, DiracMatrix(0, 0, 0, -2))))
  //
  {
    updateV(vs);
  }

  void
  updateV(const std::initializer_list<const std::vector<double> *const> vs) {
    v = std::make_unique<ScalarOperator>(
        ScalarOperator(OperatorParity::even, 1.0, NumCalc::sumVecs(vs),
                       GammaMatrix::ident, 0));
  }

public:
  virtual double reducedME(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    auto inv_threej = std::sqrt(Fb.twoj() + 1.0);
    return inv_threej * matrixEl(Fa, Fb);
  }
  virtual double matrixEl(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    return v->matrixEl(Fa, Fb) + Fb.k * cg5or->matrixEl(Fa, Fb) +
           +2.0 * cdr->matrixEl(Fa, Fb) + c2->matrixEl(Fa, Fb);
    // XXX Why extra 2 here??
  }

private:
  const double cl;
  std::unique_ptr<ScalarOperator> v = nullptr;
  const std::unique_ptr<ScalarOperator> cg5or, cdr, c2;
};
