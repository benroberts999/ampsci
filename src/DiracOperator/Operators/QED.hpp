#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "Physics/FGRadPot.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>

namespace DiracOperator {

//******************************************************************************
//! Radiative QED operator, electric part
class Hrad_el final : public ScalarOperator {
public:
  Hrad_el(const std::vector<double> &Hel)
      : ScalarOperator(Parity::even, Hel.empty() ? 0.0 : 1.0, Hel,
                       {1, 0, 0, 1}) {}
  std::string name() const override final { return "Hrad_el"; }
  std::string units() const override final { return "au"; }
};

//******************************************************************************
//! Radiative QED operator, off-diagonal magnetic part
class Hrad_mag final : public ScalarOperator {
public:
  Hrad_mag(const std::vector<double> &Hmag)
      : ScalarOperator(Parity::even, Hmag.empty() ? 0.0 : -1.0, Hmag,
                       {0, 1, 1, 0}) {}
  std::string name() const override final { return "Hrad_mag"; }
  std::string units() const override final { return "au"; }
};

//******************************************************************************
class Hrad final : public ScalarOperator {
public:
  Hrad(const std::vector<double> &Hel, const std::vector<double> &Hmag)
      : ScalarOperator(Parity::even, 1.0), Vel(Hel), Vm(Hmag) {}
  std::string name() const override final { return "Hrad"; }
  std::string units() const override final { return "au"; }

  virtual DiracSpinor radial_rhs(const int kappa_a,
                                 const DiracSpinor &Fb) const override final {
    // XXX Does this work?
    return Vel.radial_rhs(kappa_a, Fb) + Vm.radial_rhs(kappa_a, Fb);
  }

private:
  Hrad_el Vel;
  Hrad_mag Vm;
};

//******************************************************************************
//! @brief Effective VertexQED operator
/*! @details
Takes in any TensorOperator (DiracOperator) h, and forms the corresponding
effective QED vertex operator, defined:

\f[
\hat h_{\rm vertex} = A \alpha \exp(-b r / \lambda_c)
\f]

where

\f[ \lambda_c = 1/ \alpha \approx 137 \f]

A and b are fitting factors; typically b=1
 */
class VertexQED final : public TensorOperator {

public: // constructor
  VertexQED(const TensorOperator *const h0, const Grid &rgrid, double a = 1.0,
            double b = 1.0)
      : TensorOperator(
            h0->rank(), h0->parity() == 1 ? Parity::even : Parity::odd,
            h0->getc(), vertex_func(rgrid, a, b, h0->getv()), h0->get_d_order(),
            h0->imaginaryQ() ? Realness::imaginary : Realness::real,
            h0->freqDependantQ),
        m_h0(h0) {}

  std::string name() const override final {
    return m_h0->name() + "_vertexQED";
  }
  std::string units() const override final { return m_h0->units(); }

  double angularF(const int ka, const int kb) const override final {
    return m_h0->angularF(ka, kb);
  }

  double angularCff(int ka, int kb) const override final {
    return m_h0->angularCff(ka, kb);
  }
  double angularCgg(int ka, int kb) const override final {
    return m_h0->angularCgg(ka, kb);
  }
  double angularCfg(int ka, int kb) const override final {
    return m_h0->angularCfg(ka, kb);
  }
  double angularCgf(int ka, int kb) const override final {
    return m_h0->angularCgf(ka, kb);
  }

  // Have m_h0 pointer, so delete copy/asign constructors
  VertexQED(const DiracOperator::VertexQED &) = delete;
  VertexQED &operator=(const DiracOperator::VertexQED &) = delete;

private:
  const TensorOperator *const m_h0;

public:
  //! Fitting factor for hyperfine. Default a(Z)
  static double a(double z) { return 1.0 + 28.5 / z; }

  //! Takes existing radial vector, multiplies by:
  //! @details
  //!  a(Z) * a0 * exp( - b * r / a0).
  //! a0 = alpha = 1/137.
  //! b=1 by default. A should be fitted.
  //! a(Z) = 1.0 + 28.5/Z
  //! nb: can give it an empty vector, to just get the exponential function
  static std::vector<double> vertex_func(const Grid &rgrid, double a, double b,
                                         std::vector<double> v = {}) {

    const double a0 = PhysConst::alpha;
    if (v.empty()) {
      // If v is empty, means it should be {1,1,1,1,...}
      v.resize(rgrid.num_points, 1.0);
    }

    for (auto i = 0ul; i < rgrid.num_points; ++i) {
      auto exp = a * a0 * std::exp(-b * rgrid.r[i] / a0);
      v[i] *= exp;
    }
    return v;
  }
};

//******************************************************************************
//! Magnetic loop vacuum polarisation (Uehling vertex)
class MLVP final : public TensorOperator {

public:
  //! rN is nuclear radius, in atomic units
  MLVP(const TensorOperator *const h0, const Grid &rgrid, double rN)
      : TensorOperator(
            h0->rank(), h0->parity() == 1 ? Parity::even : Parity::odd,
            h0->getc(), MLVP_func(rgrid, rN, h0->getv()), h0->get_d_order(),
            h0->imaginaryQ() ? Realness::imaginary : Realness::real,
            h0->freqDependantQ),
        m_h0(h0) {}

  std::string name() const override final { return m_h0->name() + "_MLVP"; }
  std::string units() const override final { return m_h0->units(); }

  double angularF(const int ka, const int kb) const override final {
    return m_h0->angularF(ka, kb);
  }

  double angularCff(int ka, int kb) const override final {
    return m_h0->angularCff(ka, kb);
  }
  double angularCgg(int ka, int kb) const override final {
    return m_h0->angularCgg(ka, kb);
  }
  double angularCfg(int ka, int kb) const override final {
    return m_h0->angularCfg(ka, kb);
  }
  double angularCgf(int ka, int kb) const override final {
    return m_h0->angularCgf(ka, kb);
  }

  // Have m_h0 pointer, so delete copy/asign constructors
  MLVP(const DiracOperator::MLVP &) = delete;
  MLVP &operator=(const DiracOperator::MLVP &) = delete;

private:
  const TensorOperator *const m_h0;

public:
  // public since may as well be
  // This multiplies the original operator by Z(r), which is the MLVP correction
  static std::vector<double> MLVP_func(const Grid &rgrid, double rN,
                                       std::vector<double> v = {}) {
    // rN must be in atomic units

    if (v.empty()) {
      // If v is empty, means it should be {1,1,1,1,...}
      v.resize(rgrid.num_points, 1.0);
    }

    // compute the integral at each radial grid point
    for (auto i = 0ul; i < rgrid.num_points; ++i) {
      const auto Z_mvlp = FGRP::Q_MLVP(rgrid.r[i], rN);
      // multiply the operator
      v[i] *= Z_mvlp;
    }

    return v;
  }
};

} // namespace DiracOperator
