#pragma once
#include "DiracOperator/GenerateOperator.hpp"
#include "DiracOperator/Operators/hfs.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Physics/FGRadPot.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <cmath>

namespace DiracOperator {

//==============================================================================
//! Flambaum-ginges radiative potential operator
class Vrad final : public ScalarOperator {
public:
  Vrad(QED::RadPot rad_pot)
      : ScalarOperator(Parity::even, 1.0), m_Vrad(std::move(rad_pot)) {}
  std::string name() const override final { return "Vrad"; }
  std::string units() const override final { return "au"; }

  virtual DiracSpinor radial_rhs(const int kappa_a,
                                 const DiracSpinor &Fb) const override final {
    auto dF = m_Vrad.Vel(Fb.l()) * Fb;
    using namespace qip::overloads;
    const auto &Hmag = m_Vrad.Hmag(Fb.l());
    dF.f() -= Hmag * Fb.g();
    dF.g() -= Hmag * Fb.f();
    if (kappa_a != Fb.kappa())
      return 0.0 * dF;
    return dF;
  }

  virtual double radialIntegral(const DiracSpinor &Fa,
                                const DiracSpinor &Fb) const override final {
    // nb: faster not to do this, but nicer this way
    return Fa * radial_rhs(Fa.kappa(), Fb);
  }

  const QED::RadPot &RadPot() const { return m_Vrad; }

private:
  QED::RadPot m_Vrad;
};

//==============================================================================
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_Vrad(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{{"Ueh", "  Uehling (vacuum pol). [1.0]"},
        {"SE_h", "  self-energy high-freq electric. [1.0]"},
        {"SE_l", "  self-energy low-freq electric. [1.0]"},
        {"SE_m", "  self-energy magnetic. [1.0]"},
        {"WK", "  Wickman-Kroll. [0.0]"},
        {"rcut", "Maximum radius (au) to calculate Rad Pot for [5.0]"},
        {"scale_rN", "Scale factor for Nuclear size. 0 for pointlike, 1 for "
                     "typical [1.0]"},
        {"scale_l", "List of doubles. Extra scaling factor for each l e.g., "
                    "1,0,1 => include for s and d, but not for p [1.0]"},
        {"readwrite", "Read/write potential? [true]"}}});
  if (input.has_option("help")) {
    return nullptr;
  }

  const auto x_Ueh = input.get("Ueh", 1.0);
  const auto x_SEe_h = input.get("SE_h", 1.0);
  const auto x_SEe_l = input.get("SE_l", 1.0);
  const auto x_SEm = input.get("SE_m", 1.0);
  const auto x_wk = input.get("WK", 0.0);
  const auto rcut = input.get("rcut", 5.0);
  const auto scale_rN = input.get("scale_rN", 1.0);
  const auto x_spd = input.get("scale_l", std::vector{1.0});
  const auto readwrite = input.get("readwrite", true);

  const auto r_N_au =
      std::sqrt(5.0 / 3.0) * scale_rN * wf.nucleus().r_rms() / PhysConst::aB_fm;

  auto qed = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {x_Ueh, x_SEe_h, x_SEe_l, x_SEm, x_wk}, x_spd, true,
                         readwrite);

  return std::make_unique<Vrad>(std::move(qed));
}

//==============================================================================
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
            h0->freqDependantQ()),
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
      v.resize(rgrid.num_points(), 1.0);
    }

    for (auto i = 0ul; i < rgrid.num_points(); ++i) {
      auto exp = a * a0 * std::exp(-b * rgrid.r(i) / a0);
      v[i] *= exp;
    }
    return v;
  }
};

//==============================================================================
//! Magnetic loop vacuum polarisation (Uehling vertex)
class MLVP final : public TensorOperator {

public:
  //! rN is nuclear (charge) radius, in atomic units
  MLVP(const DiracOperator::hfs *const h0, const Grid &rgrid, double rN)
      : TensorOperator(
            h0->rank(), h0->parity() == 1 ? Parity::even : Parity::odd,
            h0->getc(), MLVP_func(rgrid, rN, h0->getv()), h0->get_d_order(),
            h0->imaginaryQ() ? Realness::imaginary : Realness::real,
            h0->freqDependantQ()),
        m_h0(*h0) {}

  std::string name() const override final { return "MLVP"; }
  std::string units() const override final { return m_h0.units(); }

  double angularF(const int ka, const int kb) const override final {
    return m_h0.angularF(ka, kb);
  }
  double angularCff(int ka, int kb) const override final {
    return m_h0.angularCff(ka, kb);
  }
  double angularCgg(int ka, int kb) const override final {
    return m_h0.angularCgg(ka, kb);
  }
  double angularCfg(int ka, int kb) const override final {
    return m_h0.angularCfg(ka, kb);
  }
  double angularCgf(int ka, int kb) const override final {
    return m_h0.angularCgf(ka, kb);
  }

public:
  // Store a copy?
  DiracOperator::hfs m_h0;

public:
  // public since may as well be
  // This multiplies the original operator by Z(r), which is the MLVP correction
  static std::vector<double> MLVP_func(const Grid &rgrid, double rN,
                                       std::vector<double> v = {}) {
    // rN must be in atomic units

    if (v.empty()) {
      // If v is empty, means it should be {1,1,1,1,...}
      v.resize(rgrid.num_points(), 1.0);
    }

    // compute the integral at each radial grid point
    for (auto i = 0ul; i < rgrid.num_points(); ++i) {
      const auto Z_mvlp = FGRP::Q_MLVP(rgrid.r(i), rN);
      // multiply the operator
      v[i] *= Z_mvlp;
    }

    return v;
  }
};

//==============================================================================
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_MLVP(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"rN",
        "Nuclear radius (in fm), for finite-nuclear size "
        "correction to Uehling loop. If not given, taken from wavefunction."},
       {"hfs_options{}",
        "Options for hyperfine operator that sits inside the MLVP operator. "
        "Note: should use pointlike F(r)! [see `ampsci -o hfs`]."}});
  if (input.has_option("help")) {
    return nullptr;
  }

  // 1. generate regular hfs operator
  const auto t_options = input.getBlock("hfs_options");
  auto oper_options = t_options ? *t_options : IO::InputBlock{};

  if (!oper_options.get("nuc_mag")) {
    // Typically, default is Ball. But in our case, should be pointlike
    oper_options.add(IO::Option{"nuc_mag", "pointlike"});
  }

  // 2. MLVP
  const auto rN_fm =
      input.get("rN", std::sqrt(5.0 / 3.0) * wf.nucleus().r_rms());

  if (oper_options.get("print", true)) {
    std::cout << "\nGenerate MLVP operator for hfs, with parameters:\n";
    if (rN_fm != 0.0)
      std::cout << "Using finite nuclear charge in Uehling loop, with rN="
                << rN_fm << " fm.\n";
    else
      std::cout << "Using pointlike Uehling loop.\n";
  }

  const auto h = generate_hfs(oper_options, wf);

  const auto r_N_au = rN_fm / PhysConst::aB_fm;

  return std::make_unique<MLVP>(dynamic_cast<DiracOperator::hfs *>(h.get()),
                                wf.grid(), r_N_au);
}

} // namespace DiracOperator