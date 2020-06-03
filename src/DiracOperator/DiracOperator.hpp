#pragma once
#include "Angular/Angular_369j.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cmath>
#include <string>
#include <vector>

//! Dirac Operators: General + derived
namespace DiracOperator {

enum class Parity { even, odd };
enum class Realness { real, imaginary };

//******************************************************************************
//! 4x4 Integer matrix (for Gamma/Pauli). Can be real or imag. Not mixed.
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
//! @brief General operator (virtual base class); operators derive from this.
//! @details
//! k is rank, c is multiplicative constant, d_order is derivative order,
//! pi is parity, may be Parity::even or ::odd.
//! RorI may be Realness::real or Realness::imaginary.
//! Note: You may not construct a TensorOperator. Instead, you must construct
//! one of the derived 'operators' (there are some general ones); see
//! operators.hpp for list of operators. Operators work by overrideing the
//! angularCxx() functions and angularF().
//! c, v, and Cxx are included in radial integral.
class TensorOperator {
protected:
  TensorOperator(int k, Parity pi, double c = 1,
                 const std::vector<double> &inv = {}, int d_order = 0,
                 Realness RorI = Realness::real, bool freq_dep = false)
      : m_rank(k),
        m_parity(pi),
        diff_order(d_order),
        opC(RorI),
        m_constant(c),
        m_vec(inv),
        freqDependantQ(freq_dep){};

public:
  virtual ~TensorOperator() = default;

private:
  const int m_rank;
  const Parity m_parity;
  const int diff_order;
  const Realness opC;

protected:           // these may be updated for frequency-dependant operators
  double m_constant; // included in radial integral
  std::vector<double> m_vec; // useful to be able to update this! ?

public:
  const bool freqDependantQ{false};

public:
  //! If matrix element <a|h|b> is zero, returns true
  bool isZero(const int ka, int kb) const;

  //! Update frequency for frequency-dependant operators.
  virtual void updateFrequency(const double){};

  //! Returns a const ref to vector v
  const std::vector<double> &getv() const { return m_vec; }
  //! Returns a const ref to constant c
  double getc() const { return m_constant; }

  //! returns true if operator is imaginary (has imag MEs)
  bool imaginaryQ() const { return (opC == Realness::imaginary); }
  int rank() const { return m_rank; }
  //! returns parity, as integer (+1 or -1)
  int parity() const { return (m_parity == Parity::even) ? 1 : -1; }

  //! Returns string for outputting to screen, e.g.: "<6s+||h||6p->"
  std::string rme_symbol(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  std::string R_symbol(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  //! Returns "name" of operator (e.g., 'E1')
  virtual std::string name() const { return "Operator"; };
  //! Returns units of operator (usually au, may be MHz, etc.)
  virtual std::string units() const { return "au"; };

protected:
  // These are needed for radial integrals
  // Usually just constants, but can also be functions of kappa
  virtual double angularCff(int /*k_a*/, int /*k_b*/) const { return 1.0; }
  virtual double angularCgg(int, int) const { return 1.0; }
  virtual double angularCfg(int, int) const { return 0.0; }
  virtual double angularCgf(int, int) const { return 0.0; }

public:
  //! @brief angularF: links radiation integral to RME.
  //! RME = <a||h||b> = angularF(a,b) * radial_int(a,b)
  virtual double angularF(const int, const int) const = 0;
  //! radial_int = Fa * radial_rhs(a, Fb) (a needed for angular factor)
  DiracSpinor radial_rhs(const int kappa_a, const DiracSpinor &Fb) const;
  //! ME = rme3js * RME
  double rme3js(const int twoja, const int twojb, int two_mb,
                int two_q = 0) const;

  //! <a||h||b> = Fa * reduced_rhs(a, Fb) (a needed for angular factor)
  DiracSpinor reduced_rhs(const int ka, const DiracSpinor &Fb) const;

  //! <b||h||a>  = Fa * reduced_lhs(a, Fb) (a needed for angular factor)
  DiracSpinor reduced_lhs(const int ka, const DiracSpinor &Fb) const;

  //! Defined via <a||h||b> = angularF(a,b) * radialIntegral(a,b)
  double radialIntegral(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
  double reducedME(const DiracSpinor &Fa, const DiracSpinor &Fb) const;
};

//****************************************************************************
//****************************************************************************
//! Speacial case for scalar operator
class ScalarOperator : public TensorOperator {
public:
  ScalarOperator(Parity pi, double in_coef,
                 const std::vector<double> &in_v = {},
                 const IntM4x4 &in_g = GammaMatrix::ident, int in_diff = 0,
                 Realness rori = Realness::real)
      : TensorOperator(0, pi, in_coef, in_v, in_diff, rori),
        c_ff(in_g.e00),
        c_fg(in_g.e01),
        c_gf(in_g.e10),
        c_gg(in_g.e11) {}

  ScalarOperator(const std::vector<double> &in_v)
      : TensorOperator(0, Parity::even, 1.0, in_v, 0),
        c_ff(1.0),
        c_fg(0.0),
        c_gf(0.0),
        c_gg(1.0) {}

  ScalarOperator(double in_coef, const std::vector<double> &in_v = {})
      : TensorOperator(0, Parity::even, in_coef, in_v, 0),
        c_ff(1.0),
        c_fg(0.0),
        c_gf(0.0),
        c_gg(1.0) {}

public:
  virtual double angularF(const int ka, const int kb) const override {
    // |k| = j+1/2
    return (std::abs(ka) == std::abs(kb)) ? std::sqrt(2.0 * std::abs(ka)) : 0.0;
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
//! Speacial operator: 0
class NullOperator : public ScalarOperator {
public:
  NullOperator() : ScalarOperator(Parity::even, 0, {}) {}

protected:
  double virtual angularCff(int, int) const override { return 0.0; }
  double virtual angularCgg(int, int) const override { return 0.0; }
  double virtual angularCfg(int, int) const override { return 0.0; }
  double virtual angularCgf(int, int) const override { return 0.0; }
};

} // namespace DiracOperator
