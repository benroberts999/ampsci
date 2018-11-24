#include "StandardHaloModel.h"
#include <cmath>
#include <vector>

using namespace SHMCONSTS;


StandardHaloModel::StandardHaloModel(double in_cosphi, double dves, double dv0)
{
  cosphi = in_cosphi;

  v0    = V0   + dv0*DEL_V0;   //account for errors
  v_sun = VSUN + dv0*DEL_V0;   //account for errors
  vesc  = VESC + dves*DEL_VESC; //account for errors
  veorb = VEORB;

  NormConst = 1;
  double tmp = normfv();
  NormConst = tmp;
}


//******************************************************************************
double StandardHaloModel::fv(double v)
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
  //local velocity (lab)
  double vl = v_sun + veorb*COSBETA*cosphi;

  //Norm const * v^2:
  double A = NormConst*pow(v,2);

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
double StandardHaloModel::normfv()
{
  int num_vsteps = 2000;
  double dv = MAXV/num_vsteps;

  double v = dv;
  double A = 0;
  for(int i=0; i<num_vsteps; i++){
    A += fv(v);
    v += dv;
  }
  return 1./(A*dv);
}
