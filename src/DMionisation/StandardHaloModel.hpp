#pragma once
#include <vector>

namespace SHMCONSTS {
constexpr double VESC = 550.; // galactic escape velocity [1804.01231]
constexpr double DEL_VESC = 55.;
// Range: 498 - 608

constexpr double V0 = 220;     // 1804.01231 + RevModPhys.85.1561 OR 235
constexpr double DEL_V0 = 20.; // range 220 - 235

constexpr double VSUN = V0 + 13; // [1804.01231]
// const double DEL_VSUN = 15.; //same as above!

constexpr double VEORB = 29.8; //[RevModPhys.85.1561]

constexpr double COSBETA = 0.49; // 1804.01231 (Earth inclination to sun dir)

constexpr double VEROTEQ = 0.47; // earth rotation speed (approx) @ equator

constexpr double MAXV = VESC + VSUN + VEORB + VEROTEQ + DEL_V0 + DEL_VESC;
} // namespace SHMCONSTS

//******************************************************************************
class StandardHaloModel {
public:
  StandardHaloModel(double in_cosphi = 0, double in_dves = 0,
                    double in_dv0 = 0);
  double cosphi;
  double v0, v_sun, vesc, veorb;
  double fv(double v) const;

  std::vector<double> makeFvArray(const std::vector<double> &v_array) const;

private:
  double m_normConst = 1;
  double normfv() const;
};
