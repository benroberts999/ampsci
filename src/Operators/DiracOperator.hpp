#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <cmath>
#include <string>
#include <vector>

enum class OperatorParity { even, odd };
enum class OperatorC { real, imaginary };

//******************************************************************************
struct IntM4x4
// 4x4 Integer matrix. Can be Pure real, or pure imag. Not mixed
// Notation for elements:
//  (e00  e01)
//  (e10  e11)
{
  IntM4x4(int in00, int in01, int in10, int in11, bool in_imag = false)
      : e00(in00), e01(in01), e10(in10), e11(in11), imaginary(in_imag) {}

  const int e00, e01, e10, e11;
  const bool imaginary;
  IntM4x4 operator*(const IntM4x4 &other) const;
  IntM4x4 operator+(const IntM4x4 &other) const;
  IntM4x4 operator-(const IntM4x4 &other) const;
};

//******************************************************************************
// Make some 'global' constants - common dirac matrices
// Assumes DIRAC BASIS
namespace PauliSpinMatrix {
const IntM4x4 ident(1, 0, 0, 1);
const IntM4x4 sx(0, 1, 1, 0);
const IntM4x4 sy(0, -1, 1, 0, true);
const IntM4x4 sz(1, 0, 0, -1);
} // namespace PauliSpinMatrix

namespace GammaMatrix {
const IntM4x4 ident(1, 0, 0, 1);
const IntM4x4 g0(1, 0, 0, -1);
const IntM4x4 g5(0, 1, 1, 0);
// const IntM4x4 ig5(0, 1, -1, 0); /*?*/
} // namespace GammaMatrix

//******************************************************************************
class DiracOperator {
protected:
  DiracOperator(int k, OperatorParity pi, double c = 1,
                const std::vector<double> &inv = {}, int d_order = 0,
                OperatorC RorI = OperatorC::real)
      : m_rank(k), m_parity(pi), constant(c), vec(inv), diff_order(d_order),
        opC(RorI) //
        {};

public:
  virtual ~DiracOperator() = default;

private:
  const int m_rank;
  const OperatorParity m_parity;
  const double constant;
  const std::vector<double> vec; // useful to be able to update this! ?
  const int diff_order;
  const OperatorC opC;

public:
  bool isZero(const int ka, int kb) const;

  const std::vector<double> &getv() const { return vec; }
  double getc() const { return constant; }

  bool imaginaryQ() const { return (opC == OperatorC::imaginary); }
  int rank() const { return m_rank; }
  int parity() const { return (m_parity == OperatorParity::even) ? 1 : -1; }

  std::string rme_symbol(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  virtual std::string name() const { return "Operator"; };
  virtual std::string units() const { return "au"; };

protected:
  // These are needed for radial integrals
  // Usually just constants, but can also be functions
  virtual double angularCff(int /*k_a*/, int /*k_b*/) const { return 1.0; }
  virtual double angularCgg(int, int) const { return 1.0; }
  virtual double angularCfg(int, int) const { return 0.0; }
  virtual double angularCgf(int, int) const { return 0.0; }
  virtual double StateDepConst(const DiracSpinor &, const DiracSpinor &) const {
    return 1.0;
  }

public:
  // angularF: links radiation integral to RME.
  // RME = <a||h||b> = angularF(a,b) * radial_int(a,b)
  virtual double angularF(const int, const int) const = 0;
  DiracSpinor radial_rhs(const int kappa_a, const DiracSpinor &Fb) const;
  // ME = rme3js * RME
  double rme3js(const int twoja, const int twojb, int two_mb,
                int two_q = 0) const;

  DiracSpinor reduced_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  DiracSpinor reduced_rhs(const int ka, const DiracSpinor &Fb) const;

  DiracSpinor reduced_lhs(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  DiracSpinor reduced_lhs(const int ka, const DiracSpinor &Fb) const;

  double radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  double reducedME(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
};

//****************************************************************************
//****************************************************************************
class ScalarOperator : public DiracOperator {
public:
  ScalarOperator(OperatorParity pi, double in_coef,
                 const std::vector<double> &in_v = {},
                 const IntM4x4 &in_g = GammaMatrix::ident, int in_diff = 0,
                 OperatorC rori = OperatorC::real)
      : DiracOperator(0, pi, in_coef, in_v, in_diff, rori), c_ff(in_g.e00),
        c_fg(in_g.e01), c_gf(in_g.e10), c_gg(in_g.e11) {}

  ScalarOperator(const std::vector<double> &in_v)
      : DiracOperator(0, OperatorParity::even, 1.0, in_v, 0), c_ff(1.0),
        c_fg(0.0), c_gf(0.0), c_gg(1.0) {}

  ScalarOperator(double in_coef, const std::vector<double> &in_v = {})
      : DiracOperator(0, OperatorParity::even, in_coef, in_v, 0), c_ff(1.0),
        c_fg(0.0), c_gf(0.0), c_gg(1.0) {}

public:
  virtual double angularF(const int ka, const int kb) const override {
    // |k| = j+1/2
    // return (ka == kb) ? std::sqrt(Fb.twoj() + 1.0) : 0.0;
    return (std::abs(ka) == std::abs(kb)) ? std::sqrt(2.0 * std::abs(ka)) : 0.0;
  }

  // virtual double matrixEl(const DiracSpinor &Fa, const DiracSpinor &Fb) const
  // {
  //   return radialIntegral(Fa, Fb);
  // }

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

protected:
  double virtual angularCff(int, int) const override { return 0.0; }
  double virtual angularCgg(int, int) const override { return 0.0; }
  double virtual angularCfg(int, int) const override { return 0.0; }
  double virtual angularCgf(int, int) const override { return 0.0; }
};
