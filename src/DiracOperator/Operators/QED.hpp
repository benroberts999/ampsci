#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/RadiativePotential.hpp"
#include <cmath>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_integration.h>
// #include <gsl/gsl_sf_expint.h>

// namespace RadiativePotential {
//
// // Function signature required for GSL integration (impl below):
// inline double Z_MLVP_finite(double t, void *p);
// inline double Z_MLVP_pointlike(double t, void *p);
//
// // struct to pass the parameters to the GSL function
// struct MLVP_params {
//   // simple struct that only stores current point on the radial grid and the
//   // nuclear radius
//   double r, rN;
// };
//
// } // namespace RadiativePotential

namespace DiracOperator {

//******************************************************************************
//! Radiative QED operator, electric part
class Hrad_el final : public ScalarOperator {
public:
  Hrad_el(const std::vector<double> &Hel)
      : ScalarOperator(Parity::even, 1.0, Hel, {1, 0, 0, 1}) {}
  std::string name() const override final { return "Hrad_el"; }
  std::string units() const override final { return "au"; }
};

//******************************************************************************
//! Radiative QED operator, off-diagonal magnetic part
class Hrad_mag final : public ScalarOperator {
public:
  Hrad_mag(const std::vector<double> &Hmag)
      : ScalarOperator(Parity::even, -PhysConst::c, Hmag, {0, 1, 1, 0}) {}
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
        m_h0(h0)
  // , m_a(a), m_b(b)
  {}

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
  // const double m_a;
  // const double m_b;

public:
  //! Takes existing radial vector, multiplies by:
  //! @details
  //!  A * a0 * exp( - b * r / a0).
  //! a0 = alpha = 1/137.
  //! b=1 by default. A should be fitted.
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
      const auto Z_mvlp = RadiativePotential::vMLVP(rgrid.r[i], rN);
      // multiply the operator
      v[i] *= Z_mvlp;
    }

    return v;
  }
};

} // namespace DiracOperator

//******************************************************************************
// namespace RadiativePotential {

// // added magnetic loop function outside the name space
// // Note that this is ball nuclear magnetisation MLVP function which
// multiplies
// // 1/r**2 hence has to be executed with "pointlike" F(r)
// inline double Z_MLVP_finite(double t, void *p) {
//
//   // Main:
//   // A. V. Volotka, D. A. Glazov, I. I. Tupitsyn, N. S. Oreshkina, G.
//   // Plunien, and V. M. Shabaev, Phys. Rev. A 78, 062507 (2008) - Equation
//   (18)
//
//   // also:
//   // P. Sunnergren, H. Persson, S. Salomonson, S. M. Schneider, I. Lindgren,
//   // and G. Soff, Phys. Rev. A 58, 1055 (1998) -- Equation (58)
//   // * note: This has typo compared to above
//
//   // Note: This includes the "sphere" BW effect for hyperfine.
//   // If remove, can use this for _any_ operator??
//
//   // obtain the radial grid parameters and convert to relativistic units
//   const auto [r_au, r_N_au] = *(static_cast<MLVP_params *>(p));
//   const auto r = r_au / PhysConst::alpha;
//   const auto r_N = r_N_au / PhysConst::alpha;
//
//   // form the integrand (all calculations are preformed in relativistic
//   units) constexpr double prefactor =
//       (2.0 / 3.0) * (PhysConst::alpha / M_PI) * 3.0 / 16.0;
//   // const double ball = 1.0 / (r_N * r_N * r_N);
//   // if F(r) included in operator, use below???
//   // Which "finite-size" part is just F(r) vs Uehling??
//   // For f > rN, set rN to zero??? (since correspond to pointlike?)
//   // const double point = 1.0 / (r * r * r);
//
//   // This is to cancel out the "ball" part of F(r), since included already.
//   // Note: This does not make Z align with pointlike case - I guess this
//   // includes finite nucleus (charge) too?
//   // This way, F(r) *should* be included in operator
//   const auto cancel_ball =
//       r < r_N ? 1.0 / (r * r * r) : 1.0 / (r_N * r_N * r_N);
//
//   const double t2 = t * t;
//   const double t_part =
//       std::sqrt(t2 - 1.0) / (t2 * t) * (1.0 + 1.0 / (2.0 * t2));
//   const double plus_r = (4.0 * r * r_N + 2.0 / t * (r_N + r) + 1.0 / t2) *
//                         std::exp(-2.0 * t * (r + r_N));
//   const double minus_r =
//       (4.0 * r * r_N - 2.0 / t * std::abs(r_N - r) - 1.0 / t2) *
//       std::exp(-2.0 * t * std::abs(r_N - r));
//
//   return prefactor * cancel_ball * t_part * (plus_r + minus_r);
// }
//
// //------------------------------------------------------------------------------
// inline double Z_MLVP_pointlike(double t, void *p) {
//
//   // Main:
//   // A. V. Volotka, D. A. Glazov, I. I. Tupitsyn, N. S. Oreshkina, G.
//   // Plunien, and V. M. Shabaev, Phys. Rev. A 78, 062507 (2008) - Equation
//   (18)
//
//   // also:
//   // P. Sunnergren, H. Persson, S. Salomonson, S. M. Schneider, I. Lindgren,
//   // and G. Soff, Phys. Rev. A 58, 1055 (1998) -- Equation (58)
//   // * note: This has typo compared to above
//
//   // obtain the radial grid parameters and convert to relativistic units
//   const auto [r_au, r_N_au] = *(static_cast<MLVP_params *>(p));
//   const auto r = r_au / PhysConst::alpha;
//   (void)r_N_au; // don't use
//   // const auto r_N = r_N_au / PhysConst::alpha;
//
//   constexpr double prefactor = (2.0 / 3.0) * (PhysConst::alpha / M_PI);
//   const double t2 = t * t;
//   const double t_part = std::sqrt(t2 - 1.0) * (1.0 + 2 * t2) *
//                         (1.0 + 2 * r * t) * std::exp(-2.0 * r * t) /
//                         (2.0 * t2 * t2);
//
//   return prefactor * t_part;
// }

// } // namespace RadiativePotential
