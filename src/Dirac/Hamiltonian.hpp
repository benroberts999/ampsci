#pragma once
#include "Angular/Wigner_369j.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include <cmath>
#include <vector>

// XXX Same as parity-even, scalar operator?
// Derive from one?

//******************************************************************************
class RadialHamiltonian {
  // must be set for each kappa!
  // If called with index(kappa) > max, will return max
private:
  const Grid *const m_gr;
  // const int m_z;
  const double m_alpha, m_c, m_c2;
  std::vector<std::vector<double>> m_Vk;
  std::vector<double> m_v_mag;

public:
  RadialHamiltonian(const Grid &rgrid, const double in_alpha)
      : m_gr(&rgrid), m_alpha(in_alpha), m_c(1.0 / m_alpha), m_c2(m_c * m_c),
        m_Vk({}), m_v_mag({}) {}

  template <typename... Args> //
  void set_v(const int kappa, const Args &... args) {
    auto ki = Wigner::indexFromKappa(kappa);
    if ((int)m_Vk.size() < ki + 1)
      m_Vk.resize(ki + 1); // XXX may set some to zero??
    m_Vk[ki] = NumCalc::add_vectors(args...);
  }
  void set_v_mag(const std::vector<double> &vin) { m_v_mag = vin; }

  const std::vector<double> &get_v(int kappa = -1) const {
    auto ki = Wigner::indexFromKappa(kappa);
    return ((int)m_Vk.size() < ki + 1) ? m_Vk.back() : m_Vk[ki];
  }
  const std::vector<double> &get_v_mag() const { return m_v_mag; }
  double get_alpha() const { return m_alpha; }
  // double get_z() const { return m_z; }

  // d_r F = (c00, c01)F
  //         (c10, c11)
  double derivMat_00(int k, double, std::size_t i) const {
    auto h = m_v_mag.empty() ? 0.0 : m_v_mag[i] * m_gr->drdu[i];
    return (double(-k)) * m_gr->drduor[i] + h * m_gr->drdu[i];
  }
  double derivMat_01(int k, double en, std::size_t i) const {
    const auto &v = get_v(k);
    return m_alpha * (en + 2.0 * m_c2 - v[i]) * m_gr->drdu[i];
  }
  double derivMat_10(int k, double en, std::size_t i) const {
    const auto &v = get_v(k);
    return -m_alpha * (en - v[i]) * m_gr->drdu[i];
  }
  double derivMat_12(int k, double, std::size_t i) const {
    auto h = m_v_mag.empty() ? 0.0 : m_v_mag[i] * m_gr->drdu[i];
    return double(k) * m_gr->drduor[i] - h * m_gr->drdu[i];
  }

  double matrixEl(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
    if (Fa.k != Fb.k)
      return 0.0;
    const auto kappa = Fa.k;
    const auto max = std::min(Fa.pinf, Fb.pinf);
    const auto min = std::max(Fa.p0, Fb.p0);
    const auto &drdu = Fa.p_rgrid->drdu;

    auto dga = NumCalc::derivative(Fa.g, drdu, Fb.p_rgrid->du, 1);
    auto dgb = NumCalc::derivative(Fb.g, drdu, Fb.p_rgrid->du, 1);

    for (std::size_t i = 0; i < max; i++) {
      auto r = Fa.p_rgrid->r[i];
      dga[i] -= (kappa * Fa.g[i] / r);
      dgb[i] -= (kappa * Fb.g[i] / r);
    }

    auto D1m2 = NumCalc::integrate(1.0, min, max, Fa.f, dgb, drdu) +
                NumCalc::integrate(1.0, min, max, Fb.f, dga, drdu);

    auto Sab = NumCalc::integrate(1.0, min, max, Fa.g, Fb.g, drdu);

    const auto &v = get_v(kappa);
    auto Vab = NumCalc::integrate(1.0, min, max, Fa.f, Fb.f, v, drdu) +
               NumCalc::integrate(1.0, min, max, Fa.g, Fb.g, v, drdu);

    auto V_mag = 0.0;
    if (!m_v_mag.empty())
      V_mag = NumCalc::integrate(1.0, min, max, Fa.f, Fb.g, m_v_mag, drdu) +
              NumCalc::integrate(1.0, min, max, Fa.g, Fb.f, m_v_mag, drdu);
    // XXX include MAG!

    return (Vab - m_c * (D1m2 + 2.0 * m_c * Sab + V_mag)) * Fa.p_rgrid->du;
  }
};