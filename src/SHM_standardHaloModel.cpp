#include "SHM_standardHaloModel.h"

namespace SHM{


//******************************************************************************
double fv(double v, double sinphi, double dves, double dv0)
/*
Standard halo model for velocity distribution, in laboratory frame.
 f ~ v^2 exp(-v^2)
Note: distribution for DM particles that cross paths with Earth.

Note: normalised for all extras = 0.
Otherwise, NOT normalised!

 sinphi should be [-1,1]
 dv0 should be [0,1]..
 dves = [-1,1]
*/
{

  double v0    = V0   + dv0*DEL_V0;   //account for errors
  double v_sun = VSUN + dv0*DEL_V0;   //account for errors
  double vesc  = VESC + dves*DEL_VESC; //account for errors
  double veorb = VEORB;

  //local velocity (lab)
  double vl = v_sun + veorb*COSBETA*sinphi;

  //approx norm const:
  double Kn = sqrt(M_PI)*v0*vl*370.;
  double A = 1.126550845*pow(v,2)/Kn;

  double arg1 = -pow((v-vl)/v0,2);

  if(v<=0){
    return 0;
  }else if(v<vesc-vl){ //here
    double arg2 = -pow((v+vl)/v0,2);
    return A*(exp(arg1)-exp(arg2));
  }else if(v<vesc+vl){ //here
    double arg2 = -pow(vesc/v0,2);
    return A*(exp(arg1)-exp(arg2));
  }else{
    return 0;
  }

}


//******************************************************************************
double normfv(double sinphi, double dves, double dv0)
{
  int num_vsteps = 2000;
  double dv = MAXV/num_vsteps;

  double v = dv;
  double A = 0;
  for(int i=0; i<num_vsteps; i++){
    A += fv(v,sinphi,dves,dv0);
    v += dv;
  }
  return 1./(A*dv);
}


} //namespace
