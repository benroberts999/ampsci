#include "DMionisation/StandardHaloModel.hpp"
#include <cmath>
#include <vector>

StandardHaloModel::StandardHaloModel(double in_cosphi, double dves, double dv0)
    : cosphi(in_cosphi),
      v0(AstroConsts::v0 + dv0 * AstroConsts::delta_v0),
      v_sun(AstroConsts::v_sun + dv0 * AstroConsts::delta_v0),
      vesc(AstroConsts::v_esc + dves * AstroConsts::delta_v_esc),
      veorb(AstroConsts::v_earth_orb) {

  // NormCost should be 1, but *= more general
  // (means) it will work for _any_ existing m_normConst
  m_normConst *= normfv();
}

//******************************************************************************
double StandardHaloModel::fv(double v) const
// Standard halo model for velocity distribution, in laboratory frame.
//  f ~ v^2 exp(-v^2)
// Note: distribution for DM particles that cross paths with Earth.
//
// Note: normalised for all extras = 0.
// Otherwise, NOT normalised!
//
//  sinphi should be [-1,1]
//  dv0 should be [0,1]..
//  dves = [-1,1]
{
  // local velocity (lab)
  const double vl = v_sun + veorb * AstroConsts::cos_beta * cosphi;

  // Norm const * v^2:
  const double A = m_normConst * std::pow(v, 2);

  const double arg1 = -std::pow((v - vl) / v0, 2);

  if (v <= 0) {
    return 0;
  } else if (v < vesc - vl) { // here
    const double arg2 = -std::pow((v + vl) / v0, 2);
    return A * (std::exp(arg1) - std::exp(arg2));
  } else if (v < vesc + vl) { // here
    const double arg2 = -std::pow(vesc / v0, 2);
    return A * (std::exp(arg1) - std::exp(arg2));
  } else {
    return 0;
  }
}

//******************************************************************************
double StandardHaloModel::normfv() const {
  const int num_vsteps = 2000;
  const double dv = AstroConsts::v_max / num_vsteps;

  // Just use Reinmann sum; accurate enough (integrand is zero at boundaries)
  double v = dv;
  double A = 0.0;
  for (int i = 0; i < num_vsteps; i++) {
    A += fv(v);
    v += dv;
  }
  return 1.0 / (A * dv);
}

//******************************************************************************
std::vector<double>
StandardHaloModel::fv(const std::vector<double> &v_array) const {
  // Fills an (external) vector array with f(v) values.
  // Note: units will be km/s. AND v_array must be in same units
  // Possible to re-scale; must rescale v and fv!
  std::vector<double> fv_array;
  fv_array.reserve(v_array.size());
  for (auto v : v_array) {
    fv_array.push_back(fv(v));
  }
  return fv_array;
}
