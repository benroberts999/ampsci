#pragma once
#include <vector>

namespace AstroConsts {
//! galactic escape velocity [1804.01231]
constexpr double v_esc = 550.;
//! ~Error in VESV: Range: 498 - 608
constexpr double delta_v_esc = 55.;

//! 1804.01231 + RevModPhys.85.1561 OR 235
constexpr double v0 = 220;
//! Error in v0: range 220 - 235
constexpr double delta_v0 = 20.;

//! [1804.01231]
constexpr double v_sun = v0 + 13;

//! [RevModPhys.85.1561]
constexpr double v_earth_orb = 29.8;

//! 1804.01231 (Earth inclination to sun dir)
constexpr double cos_beta = 0.49;

//! earth rotation speed (approx) @ equator
constexpr double v_earth_rot_eq = 0.47;

//! Max v used in integrations; f(v) always 0 above this
constexpr double v_max =
    v_esc + v_sun + v_earth_orb + v_earth_rot_eq + delta_v0 + delta_v_esc;
} // namespace AstroConsts

//******************************************************************************
class StandardHaloModel {
public:
  //! in_cosphi is cosine of Earth's orbital phase [for annual modulation]:
  //! cosphi=0 is average, =1 is when earth+sun velocities add maximally, -1 is
  //! when the add minimally. dves and dv0 are error terms (in standard
  //! deviations): +1 for +1sigma error, -1 for -1sigma error - for testing
  //! sensitivity
  StandardHaloModel(double in_cosphi = 0, double in_dves = 0,
                    double in_dv0 = 0);

  double cosphi;
  double v0, v_sun, vesc, veorb;

  //! Given v [in km/s], returns f(v), [in (km/s)^-1]. f(v) is normalised as
  //! Int[f(v) {v,0,infty}]=1
  double fv(double v) const;

  //! As above, for for an array of {v}, returns array {f(v)}
  std::vector<double> fv(const std::vector<double> &v_array) const;

private:
  double m_normConst{1.0}; // will be calculated on construction
  double normfv() const;
};
