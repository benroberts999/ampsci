#pragma once
#include "Angular/Wigner369j.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/Vector.hpp"
#include <cmath>
#include <string>
#include <vector>

//! Dirac Operators: General + derived
namespace DiracOperator {

enum class Conjugate { yes, no };
enum class Parity { even, odd, blank };
enum class Realness { real, imaginary };

//! Swaps conjuagte option: yes<->no
inline Conjugate swap_conj(Conjugate conj) {
  return conj == Conjugate::no ? Conjugate::yes : Conjugate::no;
}

//! Dagger? Returns true iff conjugate = yes
inline bool is_conjQ(Conjugate conj) { return conj == Conjugate::yes; }

//! Returns Conjugate::yes if conj = true
inline Conjugate apply_conj(bool conj) {
  return conj ? Conjugate::yes : Conjugate::no;
}

//==============================================================================
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
  TensorOperator(int rank_k, Parity pi, double constant = 1.0,
                 const std::vector<double> &inv = {}, int diff_order = 0,
                 Realness RorI = Realness::real, bool freq_dep = false)
      : m_rank(rank_k),
        m_parity(pi),
        m_diff_order(diff_order),
        m_realness(RorI),
        m_freqDependantQ(freq_dep),
        m_constant(constant),
        m_vec(inv) {};

public:
  virtual ~TensorOperator() = default;

protected:
  int m_rank;
  Parity m_parity;
  int m_diff_order;
  Realness m_realness;
  bool m_freqDependantQ{false};

protected:
  // these may be updated for frequency-dependant operators
  double m_constant; // included in radial integral
  std::vector<double> m_vec;

public:
  bool freqDependantQ() const { return m_freqDependantQ; }

public:
  //! If matrix element <a|h|b> is zero, returns true
  bool isZero(int ka, int kb) const;
  bool isZero(const DiracSpinor &Fa, const DiracSpinor &Fb) const;

  bool selectrion_rule(int twoJA, int piA, int twoJB, int piB) const {
    if (twoJA == twoJB && twoJA == 0.0)
      return false;

    if (Angular::triangle(twoJA, twoJB, 2 * m_rank) == 0)
      return false;

    return (m_parity == Parity::even) == (piA == piB);
  }

  //! Update frequency for frequency-dependant operators.
  virtual void updateFrequency(double) { return; };

  //! Permanently re-scales the operator by constant, lambda
  void scale(double lambda);

  //! Returns a const ref to vector v
  const std::vector<double> &getv() const { return m_vec; }
  //! Returns a const ref to constant c
  double getc() const { return m_constant; }
  int get_d_order() const { return m_diff_order; }

  //! returns true if operator is imaginary (has imag MEs)
  bool imaginaryQ() const { return (m_realness == Realness::imaginary); }
  int rank() const { return m_rank; }
  //! returns parity, as integer (+1 or -1)
  int parity() const { return (m_parity == Parity::even) ? 1 : -1; }

  //! returns relative sign between <a||x||b> and <b||x||a>
  int symm_sign(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    const auto sra_i = imaginaryQ() ? -1 : 1;
    const auto sra = Angular::neg1pow_2(Fa.twoj() - Fb.twoj());
    return sra_i * sra;
  }

  //! Sign imparted by conjugation. -1 if both conjugated and imaginary
  int conj_sign(Conjugate conj = Conjugate::yes) const {
    return conj == Conjugate::yes && m_realness == Realness::imaginary ? -1 : 1;
  }

  //! Returns "name" of operator (e.g., 'E1')
  virtual std::string name() const { return "Operator"; };
  //! Returns units of operator (usually au, may be MHz, etc.)
  virtual std::string units() const { return "au"; };

public:
  // These are needed for radial integrals
  // Usually just constants, but can also be functions of kappa
  virtual double angularCff(int /*k_a*/, int /*k_b*/) const { return 1.0; }
  virtual double angularCgg(int, int) const { return 1.0; }
  virtual double angularCfg(int, int) const { return 0.0; }
  virtual double angularCgf(int, int) const { return 0.0; }

public:
  //! @brief angularF: links radiation integral to RME.
  //! RME = <a||h||b> = angularF(a,b) * radial_int(a,b)
  virtual double angularF(int, int) const = 0;
  //! radial_int = Fa * radial_rhs(a, Fb) (a needed for angular factor)
  virtual DiracSpinor radial_rhs(int kappa_a, const DiracSpinor &Fb) const;

  //! Defined via <a||h||b> = angularF(a,b) * radialIntegral(a,b)
  //! (Note: if radial_rhs is overridden, then radialIntegral must also be_
  virtual double radialIntegral(const DiracSpinor &Fa,
                                const DiracSpinor &Fb) const;

  //! ME = rme3js * RME
  double rme3js(int twoja, int twojb, int two_mb = 1, int two_q = 0) const;

  //! <a||h||b> = Fa * reduced_rhs(a, Fb) (a needed for angular factor)
  DiracSpinor reduced_rhs(int ka, const DiracSpinor &Fb,
                          Conjugate conj = Conjugate::no) const;

  //! <b||h||a>  = Fa * reduced_lhs(a, Fb) (a needed for angular factor)
  DiracSpinor reduced_lhs(int ka, const DiracSpinor &Fb) const;

  //! The reduced matrix element
  double reducedME(const DiracSpinor &Fa, const DiracSpinor &Fb,
                   Conjugate conj = Conjugate::no) const;

  //! Returns "full" matrix element, for optional (ma, mb, q) [taken as int 2*].
  //! If not specified, returns z-component (q=0), with ma=mb=min(ja,jb)
  double fullME(const DiracSpinor &Fa, const DiracSpinor &Fb,
                std::optional<int> two_ma = std::nullopt,
                std::optional<int> two_mb = std::nullopt,
                std::optional<int> two_q = std::nullopt) const;
};

//============================================================================
//============================================================================
//! Speacial case for scalar operator
class ScalarOperator : public TensorOperator {
public:
  ScalarOperator(Parity pi, double in_coef,
                 const std::vector<double> &in_v = {},
                 const std::array<int, 4> &in_g = {1, 0, 0, 1}, int in_diff = 0,
                 Realness rori = Realness::real)
      : TensorOperator(0, pi, in_coef, in_v, in_diff, rori),
        c_ff(in_g[0]),
        c_fg(in_g[1]),
        c_gf(in_g[2]),
        c_gg(in_g[3]) {}

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
    // For scalar operators, <a||h||b> = RadInt / 3js
    // 3js:= 1/(Sqrt[2j+1]) ... depends on m???
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
class NullOperator final : public ScalarOperator {
public:
  NullOperator() : ScalarOperator(Parity::even, 0, {}) {}

protected:
  double virtual angularCff(int, int) const override final { return 0.0; }
  double virtual angularCgg(int, int) const override final { return 0.0; }
  double virtual angularCfg(int, int) const override final { return 0.0; }
  double virtual angularCgf(int, int) const override final { return 0.0; }
};

} // namespace DiracOperator
