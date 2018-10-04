#include "SHM_standardHaloModel.h"

namespace SHM{

  // const double vesc = 544.; // galactic escape velocity
  // const double vL0 = 220.; // local frame velocity, average
  // const double vc  = 220; // circular velocity
  // const double vearth = 30.; //XXX update! ??
  //
  // const double max_v = vesc + vL0 + vearth;

//******************************************************************************
double fv(double v, double phi)
/*
Standard halo model for velocity distribution, in laboratory frame.
 f ~ v^2 exp(-v^2)
Note: distribution for DM particles that cross paths with Earth.
We should have: <v> = 370
XXX - Includes phase - BUT not 100% sure it's correct! XXX XXX XXX
*/
{

  double vl = vL0 + vearth*sin(phi);

  double Knorm1 = 0.91706*0.942099; //for phi=0
  double Kn = Knorm1*sqrt(M_PI)*vc*vl*370.;
  double A = pow(v,2)/Kn;

  // This is a rough way to enforce normalisation thoughout the year!
  if(phi!=0){ //save some evals if phi=0
    A /= 1. + 0.0587669*sin(phi);
    A /= 1.00137 - 0.00137*cos(2*phi) - 0.00013*sin(phi);//0.000123*sin(phi);
  }
  //Probably, this isn't good enough... XXX

  double arg1 = -pow((v-vl)/vc,2);

  if(v<=0){
    return 0; //just for safety - should never be called with v<0
  }else if(v<vesc-vl){
    double arg2 = -pow((v+vl)/vc,2);
    return A*(exp(arg1)-exp(arg2));
  }else if(v<vesc+vl){
    double arg2 = -pow(vesc/vc,2);
    return A*(exp(arg1)-exp(arg2));
  }else{
    return 0;
  }

}



} //namespace
