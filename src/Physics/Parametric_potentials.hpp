#pragma once
#include <vector>

namespace Parametric {

double green(int Z, double r, double H, double d);
double tietz(int Z, double r, double g, double t);
int defaultGreen(int z, double &H, double &d);
int defaultTietz(int z, double &t, double &g);
int defaultGreenCore(int z, double &H, double &d);

std::vector<double> GreenPotential(int z, const std::vector<double> &r_array,
                                   double H = 0, double d = 0);
std::vector<double> TietzPotential(int z, const std::vector<double> &r_array,
                                   double g = 0, double t = 0);
} // namespace Parametric
