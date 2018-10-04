#include <vector>
#include <cmath>

namespace SHM{

  const double vesc = 544.; // galactic escape velocity
  const double vL0 = 220.; // local frame velocity, average
  const double vc  = 220; // circular velocity
  const double vearth = 30.; //XXX update! ??

  const double max_v = vesc + vL0 + vearth;

  double fv(double v, double phi);

}
