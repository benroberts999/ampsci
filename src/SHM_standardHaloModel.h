#include <vector>
#include <cmath>

namespace SHM{

  const double VESC = 550.; // galactic escape velocity [1804.01231]
  const double DEL_VESC = 55.;
  //Range: 498 - 608

  const double V0  = 220;   //1804.01231 + RevModPhys.85.1561 OR 235
  const double DEL_V0 = 20.; //range 220 - 235

  const double VSUN = V0 + 13; // [1804.01231] - or 220 {older}
  //const double DEL_VSUN = 15.; //same as above!

  const double VEORB = 29.8; //[RevModPhys.85.1561]

  const double COSBETA = 0.49; // 1804.01231 (Earth inclination to sun dir)

  const double VEROTEQ = 0.47; //earth rotation speed (approx) @ equator

  const double MAXV = VESC + VSUN + VEORB + VEROTEQ;

  double fv(double v, double sinphi=0, double dves=0, double dv0=0);
  double normfv(double sinphi=0, double dves=0, double dv0=0);

}
