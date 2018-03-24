#ifndef _PARAMETRICPOT_H
#define _PARAMETRICPOT_H
#include<cmath>

double PRM_green(int Z, double r, double H, double d);
double PRM_tietz(int Z, double r, double g, double t);
int PRM_defaultGreen(int z, double &H, double &d);
int PRM_defaultTietz(int z, double &t, double &g);

#endif
