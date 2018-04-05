#include "adamsSolveLocalContinuum.h"
/*
2018-04-04.
Program to solve single-electron continuum-state Dirac problem for a (given)
local, central potential.
Uses "outwardAM" (from adamsSolveLocalBS.cpp) to solve dirac equation

Normalises wavefunction by solving all the way out to asymptotic region, where
solution should be sinosoidal. Then, matches amplitude with low-r expansion.
Therefore, need grid to go very far out (v. large r).
Also, need at least ~10 points per half-period. More points for higher energy!

===== To do :: Jan 2018 =====
  * p,q -> f,g!

*/


//******************************************************************************
int solveContinuum(std::vector<double> &p, std::vector<double> &q, double en,
    std::vector<double> v, int Z, int ka,
    std::vector<double> r, std::vector<double> drdt, double h, int NGP,
    double alpha)
{

  //Perform the "outwards integration"
  outwardAM(p,q,en,v,Z,ka,r,drdt,h,NGP,NGP-1,alpha);
  // XXX shuod be Z or Zion here?? XXX


  return 0;
}



//******************************************************************************
double fitQuadratic(double x1, double x2, double x3,
    double y1, double y2, double y3)
/*
Takes in three points, and fits them to a quadratic function.
Returns y-value for vertex of quadratic.
Used for finding the amplitude of a sine/cosine function, given thee points.
i.e., will return amplitude of since function.
Note: the given 3 points _MUST_ be close to maximum, otherwise, fit wont work
*/
{

  if(y1<0) y1=fabs(y1);
  if(y2<0) y2=fabs(y2);
  if(y3<0) y3=fabs(y3);

  double d = (x1 - x2)*(x1 - x3)*(x2 - x3);
  double Ad = x3*(x2*(x2 - x3)*y1 + x1*(-x1 + x3)*y2) + x1*(x1 - x2)*x2*y3;
  double Bd = x3*x3*(y1 - y2) + x1*x1*(y2 - y3) + x2*x2*(-y1 + y3);
  double Cd = x3*(-y1 + y2) + x2*(y1 - y3) + x1*(-y2 + y3);
  double y0 = (Ad/d) - Bd*Bd/(4.*Cd*d);

  //Find largest input y:
  double ymax=y2;
  if(y1>ymax) ymax=y1;
  if(y3>ymax) ymax=y3;

  if(ymax>y0) y0 = ymax; //y0 can't be less than (y1,y2,y3)

  return y0;

}
