#pragma once

namespace PRM {

double green(int Z, double r, double H, double d);
double tietz(int Z, double r, double g, double t);
int defaultGreen(int z, double &H, double &d);
int defaultTietz(int z, double &t, double &g);
int defaultGreenCore(int z, double &H, double &d);

} // namespace PRM
