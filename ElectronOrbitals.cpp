//class ElectronOrbitals::
#include "ElectronOrbitals.h"
#include "physicalConstants.h"
#include "atomInfo.h"

ElectronOrbitals(std::string s_in_z, int in_a, int in_ngp)
{
  //Work out Z from given atomic symbol
  int iz;
  for(iz=1; iz<200; iz++){
    if(atinfo_sym[iz]==s_in_z || s_in_z==std::to_string(iz)) break;
  }
  //XXX needs some kind of safety-check! XXX
  ElectronOrbitals(int iz, int in_a, int in_ngp);
}

ElectronOrbitals(int in_z, int in_a, int in_ngp)
{
  ngp=in_ngp;
  z=in_z;
  if(in_a==0) a=atinfo_a[z]; //Use default atomic mass
  else a=in_a;

}

//******************************************************************************
int formRadialGrid()
{
  // XXX NOTE: There are several options for the grids!
  // See Dzuba code! Which is better? Option for either??
  //I _think_ this is the Johnson grid.... check.
  // XXX AND add explanation to comments!

  //XXX put a safety check??

  double r0=1.e-4; // XXX input?? private variable? XXX
  // XXX copied from before. WHY like this???
  double paramRmax=500;
  double h=log(paramRmax/r0)/(ngp-2); //XXX ok??


  for(int i=0; i<ngp; i++){
    double temp_drdt = r0*exp(i*h);
    drdt.push_back(temp_drdt);
  }

  for(int i=0; i<ngp; i++){
    // Is it OK that it starts at 0?? should it be r0? 0.01*r0?
    // XXX Check Johnson book..?
    double temp_r = drdt[i]-r0;
    r.push_back(temp_r);
  }

  // (dr/dt)/r [for convinience]
  dror.push_back(0.); //XXX is this correct?? XXX
  for(int i=1; i<ngp; i++){
    double temp_dror = drdt[i]/r[i];
    dror.push_back(temp_dror);
  }

  return 0;
}




// Form vnuc [few options!]

//re-size vectors?? nah, later!
