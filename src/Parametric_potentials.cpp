#include "Parametric_potentials.hpp"
#include "AtomInfo.hpp"
#include "PhysConst_constants.hpp"
#include <cmath>
#include <iostream>
#include <vector>

namespace Parametric {

//******************************************************************************
double green(int Z, double r, double H, double d)
// Have default H and d : scaled to fit Alkali atoms!
// Does not include nuclear part!
{
  if (r == 0)
    return 0;
  return (-1. / r) * (1. + (double(Z) - 1.) / (H * (exp(r / d) - 1) + 1)) +
         double(Z) / r;
}

//******************************************************************************
double tietz(int Z, double r, double t, double g)
// TIETZ parametric potential
// Does NOT include nuclear potential
{
  if (r == 0)
    return 0;
  return (-1. / r) *
             (1. + (double(Z) - 1.) * exp(-g * r) / pow(1. + t * r, 2)) +
         double(Z) / r;
}

//******************************************************************************
int defaultGreenCore(int z, double &H, double &d)
// Fitted to match HF closed-shell core energies - aid convergance!
{
  // double zHd[18][3]
  std::vector<std::vector<double>> zHd = {
      {2, 12.275, 3.552},     {3, 3.248, 2.518},       {5, 6.812, 4.535},
      {11, 3.001, 1.038},     {13, 9.837, 2.99813},    {19, 5.55771, 1.65891},
      {29, 2.5, 0.65},        {31, 2.932, 0.736},      {37, 4.7987, 1.06374},
      {47, 1.95883, 0.52759}, {49, 3.363, 0.75499},    {55, 5.846, 1.145},
      {70, 3.7, 0.66},        {79, 4.45331, 0.73931},  {81, 4.7354, 0.7749},
      {87, 6.173, 0.967},     {102, 4.33716, 0.67213}, {113, 5.10036, 0.75391}};

  for (auto &item : zHd) {
    int iz = (int)item[0];
    if (z == iz) {
      H = item[1];
      d = item[2];
      return 0;
    }
  }

  for (std::size_t i = 1; i < zHd.size(); i++) {
    int izm = (int)zHd[i - 1][0];
    int izp = (int)zHd[i][0];
    if (z > izm && z < izp) {
      double H1 = zHd[i - 1][1];
      double H2 = zHd[i][1];
      double d1 = zHd[i - 1][2];
      double d2 = zHd[i][2];

      H = H1 + (H2 - H1) * (z - izm) / (izp - izm);
      d = d1 + (d2 - d1) * (z - izm) / (izp - izm);
      // std::cout<<H<<" "<<d<<"\n";
      return 0;
    }
  }

  if (z > 113) {
    H = 5.;
    d = 0.75;
  }

  return 0;
}

//******************************************************************************
int defaultGreen(int z, double &H, double &d)
/*
Default values for Green potential parameters.
Determined by fitting the 5 lowest J=1/2 states.
Crude quadratic fit used for other Z values.
*/
{

  if (z == 3) {
    H = 0.65074;
    d = 0.37017;
    return 0;
  } else if (z == 11) {
    H = 1.32493;
    d = 0.48946;
    return 0;
  } else if (z == 19) {
    H = 2.49872;
    d = 0.78867;
    return 0;
  } else if (z == 37) {
    H = 3.60424;
    d = 0.80155;
    // Johnson:
    // H=3.4811;
    // d=0.7855;
    return 0;
  } else if (z == 49) {
    H = 3.16518;
    d = 0.71266;
    return 0;
  } else if (z == 55) {
    H = 4.83108;
    d = 0.93920;
    // Johnson:
    // H=4.4691;
    // d=0.8967;
    return 0;
  } else if (z == 79) {
    // Johnson:
    H = 4.4560;
    d = 0.7160;
    return 0;
  } else if (z == 81) {
    H = 4.44226;
    d = 0.72244;
    // Johnson:
    // H=4.4530;
    // d=0.7234;
    return 0;
  } else if (z == 87) {
    H = 7.02125;
    d = 1.00633;
    return 0;
  }

  // If not one of the defaults, use crude model fit
  H = 0.665356 + 0.0747058 * z - 0.000155671 * pow(z, 2);
  d = 0.433183 + 0.0103617 * z - 0.0000531958 * pow(z, 2);
  return 0;
}

//******************************************************************************
int defaultTietz(int z, double &t, double &g)
/*
Default values for Green potential parameters.
Determined by fitting the 5 lowest J=1/2 states.
Crude quadratic fit used for other Z values.
*/
{
  if (z == 3) {
    t = 0.10000;
    g = 1.94363;
    return 0;
  } else if (z == 5) {
    t = 0.81084;
    g = 0.24086;
    return 0;
  } else if (z == 11) {
    t = 0.60635;
    g = 1.38417;
    return 0;
  } else if (z == 13) {
    t = 1.22936;
    g = 0.26552;
    return 0;
  } else if (z == 19) {
    t = 1.13831;
    g = 0.55698;
    return 0;
  } else if (z == 31) {
    t = 2.05009;
    g = 0.18498;
    return 0;
  } else if (z == 37) {
    t = 1.65346;
    g = 0.46553;
    // Johnson:
    // t=1.9530;
    // g=0.2700;
    return 0;
  } else if (z == 49) {
    t = 2.36372;
    g = 0.09911;
    return 0;
  } else if (z == 55) {
    t = 1.91561;
    g = 0.31948;
    // Johnson:
    // t=2.0453;
    // g=0.2445;
    return 0;
  } else if (z == 79) {
    t = 2.4310;
    g = 0.3500;
    return 0;
  } else if (z == 81) {
    // t=2.90583;
    // g=0.10749;
    // Johnson:
    t = 2.3537;
    g = 0.3895;
    return 0;
  } else if (z == 87) {
    t = 2.03525;
    g = 0.50664;
    return 0;
  }

  // If not one of the defaults, use crude model fit
  t = 0.33991 + 0.0495393 * z - 0.000272772 * pow(z, 2);
  g = 1.14316 - 0.028863 * z + 0.000205713 * pow(z, 2);
  return 0;
}

//******************************************************************************
std::vector<double> GreenPotential(int z, const std::vector<double> &r_array,
                                   double h, double d) {
  // double Gh, Gd; // Green potential parameters
  if (fabs(h * d) < 1.0e-6)
    defaultGreenCore(z, h, d);
  // Fill the the potential, using Greens Parametric
  std::vector<double> v;
  v.reserve(r_array.size());
  for (const auto r : r_array) {
    v.push_back(green(z, r, h, d));
  }
  return v;
}

std::vector<double> TietzPotential(int z, const std::vector<double> &r_array,
                                   double g, double t) {
  // double Gh, Gd; // Green potential parameters
  if (fabs(g * t) < 1.0e-6)
    defaultTietz(z, t, g);
  // Fill the the potential, using Greens Parametric
  std::vector<double> v;
  v.reserve(r_array.size());
  for (const auto r : r_array) {
    v.push_back(tietz(z, r, g, t));
  }
  return v;
}

} // namespace Parametric
