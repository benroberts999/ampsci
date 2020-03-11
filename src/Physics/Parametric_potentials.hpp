#pragma once
#include <vector>

//! Set of two-parameter parametric electron potentials
namespace Parametric {

double green(int Z, double r, double H, double d);
double tietz(int Z, double r, double g, double t);

//! Returns "default" parameters H and d, given Z; tuned for core to match HF
int defaultGreenCore(int z, double &H, double &d);
//! Returns "default" parameters H and d, given Z
int defaultGreen(int z, double &H, double &d);
//! Note: g and t are backwards...
int defaultTietz(int z, double &t, double &g);

std::vector<double> GreenPotential(int z, const std::vector<double> &r_array,
                                   double H = 0, double d = 0);
std::vector<double> TietzPotential(int z, const std::vector<double> &r_array,
                                   double g = 0, double t = 0);
} // namespace Parametric
