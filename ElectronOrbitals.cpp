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
  Z=in_z;
  if(in_a==0) A=atinfo_a[Z]; //Use default atomic mass
  else A=in_a;

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


//******************************************************************************
int sphericalNucleus()
{

//XXX change, so i can input a different charge radius!! XXX

  double rN; //nuclear charge radius:
  //Estimate nuclear charge radius. Only for spherical nuclei.
  //https://www-nds.iaea.org/radii/
  if(A==1) rN = 0.8783;       // 1-H
  else if(A==4) rN = 1.6755;  // 4-He
  else if(A==7) rN = 2.4440;  // 7-Li
  else if(A<10) rN = 1.15*pow(A,0.333);
  else rN = (0.836*pow(A,0.333)+0.570);
  rN/=ABOHR_FM;
  // XXX Add data tables of nuclear radii!

  //Fill the vnuc array with spherical nuclear potantial
  vnuc.push_back(0); //XXX ??
  for(int i=1; i<ngp; i++){
    double temp_v;
    //XXX double checK!!!
    //XXX also: don't need to evaluate rn^3 etc. more than once!!!
    if(r[i]<rN) temp_v = Z*(pow(r[i],2)-3.*pow(rN,2))/(2.*pow(rN,3));
    else temp_v = -Z/r[i];
    vnuc.push_back(temp_v);
  }

  return 0;
}


int fermiNucleus(double t, double c)
/*
Uses a Fermi-Dirac distribution for the nuclear potential.

rho(r) = rho_0 {1 + Exp[(r-c)/a]}^-1
V(r) = -(4 Pi)/r [A+B]
  A = Int[ rho(x) x^2 , {x,0,r}]
  B = r * Int[ rho(x) x , {x,r,infty}]
rho_0 is found by either:
  * V(infinity) = -Z/r , or equivilantly
  * \int rho(r) d^3r = Z

Depends on:
  * t: skin thickness [90 to 10% fall-off range]
    note: t = a[4 ln(3)]
  * c: half-density raius [rho(c)=0.5 rho0]

t and c are input values. In 'fermi' of fm (femto metres)
If provided with 0, will use 'default' values, approx. formula.

V(r) is expressed in terms of Complete Fermi-Dirac intagrals.
These are computed using the GSL libraries.
https://www.gnu.org/software/gsl/manual/html_node/Complete-Fermi_002dDirac-Integrals
*/
{

  //XXX test? clear vnuc!?

  if(t==0) t=2.4; // Default skin-thickness (in fm)
  if(c==0) c=1.1*pow(A,0.3333); //default half-charge radius ????
  // XXX Better approx! +/or data tables!

  double a=0.22756*t; // a = t*[4 ln(3)]
  double coa=c/a;
  // Use GSL for the Complete Fermi-Dirac Integrals:
  double F2 = gsl_sf_fermi_dirac_2(coa);
  double pi2 = pow(M_PI,2);
  vnuc.push_back(0); //XXX ??
  for(int i=1; i<ngp; i++){
    double t_r = r[i];
    double t_v = -Z/t_r;
    if(t_r<10.*a){  // XXX OK??
      double roa = ABOHR_FM*t_r/a; // convert fm <-> atomic
      double coa2= pow(coa,2);
      double xF1 = gsl_sf_fermi_dirac_1(roa-coa);
      double xF2 = gsl_sf_fermi_dirac_2(roa-coa);
      double tX  = -pow(roa,3) - 2*coa*(pi2+coa2) + roa*(pi2+3*coa2)
                 +  6*roa*xF1 - 12*xF2;
      t_v+=t_v*tX/(12.*F2);
    }
    vnuc.push_back(t_v);
  }

  return 0;
}
