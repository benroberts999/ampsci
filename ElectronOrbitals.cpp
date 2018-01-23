//class ElectronOrbitals::
#include "ElectronOrbitals.h"
#include "physicalConstants.h"
#include "atomInfo.h"
//#include "adamsSolveLocalBS.h"

//******************************************************************************
//ElectronOrbitals(std::string s_in_z, int in_a, int in_ngp)
//{
//  //Work out Z from given atomic symbol
//  int iz;
//  for(iz=1; iz<200; iz++){
//    if(atinfo_sym[iz]==s_in_z || s_in_z==std::to_string(iz)) break;
//  }
//  //XXX needs some kind of safety-check! XXX
//  ElectronOrbitals(int iz, int in_a, int in_ngp);
//}
//-----Overloaded---------------------------------------------------------------
ElectronOrbitals::ElectronOrbitals(int in_z, int in_a, int in_ngp, double var_alpha)
{

  ngp=in_ngp;
  formRadialGrid();

  alpha=ALPHA*var_alpha;

  Z=in_z;
  if(in_a==0) A=atinfo_a[Z]; //Use default atomic mass
  else A=in_a;
  zeroNucleus();
  //sphericalNucleus(); //input rnuc?
  //fermiNucleus();

}

//******************************************************************************
int ElectronOrbitals::localBoundState(int in_max_n, int in_max_l)
/*
``Wrapper'' function, to use the adamsSolveLocalBS method!
*/
{

  //double alpha=ALPHA;

  max_n = in_max_n;
  max_l = in_max_l;

  for(int n=1; n<=max_n; n++){
    for(int i=1; i<2*n; i++){ //loop through each kappa state
      int k = pow(-1,i)*ceil(0.5*i);
      int l = fabs(k+1/2)-1/2;
      if(l>max_l) continue;
      nlist.push_back(n);
      klist.push_back(k);
      int pinf,its;
      double eps;
      double en_a = -0.5*pow((double)Z/n,2);
      std::vector<double> p_a(ngp);
      std::vector<double> q_a(ngp);
      solveDBS(p_a,q_a,en_a,vnuc,Z,n,k,r,drdt,h,ngp,pinf,its,eps,alpha);
      p.push_back(p_a);
      q.push_back(q_a);
      en.push_back(en_a);
      //store convergance info:
      pinflist.push_back(pinf);
      itslist.push_back(its);
      epslist.push_back(eps);
    }
  }



  return 0;
}


//******************************************************************************
int ElectronOrbitals::formRadialGrid()
{
  // XXX NOTE: There are several options for the grids!
  // See Dzuba code! Which is better? Option for either??
  //I _think_ this is the Johnson grid.... check.
  // XXX AND add explanation to comments!

  //XXX put a safety check??

  double r0=1.e-5; // XXX input?? private variable? XXX
  // XXX copied from before. WHY like this???
  double paramRmax=500;
  h=log(paramRmax/r0)/(ngp-2); //XXX ok??


  drdt.clear();
  for(int i=0; i<ngp; i++){
    double temp_drdt = r0*exp((i+1)*h); //XXX -1? check!
    drdt.push_back(temp_drdt);
  }

  r.clear();
  //r.push_back(1.e-8);
  for(int i=0; i<ngp; i++){
    // Is it OK that it starts at 0?? should it be r0? 0.01*r0?
    // XXX Check Johnson book..?
    double temp_r = drdt[i]-r0;
    r.push_back(temp_r);
  }

  // (dr/dt)/r [for convinience]
  dror.clear();
  dror.push_back(0.); //XXX is this correct?? XXX
  for(int i=1; i<ngp; i++){
    double temp_dror = drdt[i]/r[i];
    dror.push_back(temp_dror);
  }

  return 0;
}




//******************************************************************************
double ElectronOrbitals::diracen(int z, int n, int k){
//
  // double c2 = 1./pow(alpha,2);
  // double za2 = pow(alpha*z,2);
  // double g=sqrt(k*k-za2);
	// double diracE=c2*(1./sqrt(1+za2/pow((g+n-fabs(k)),2))-1.);
  double a2 = pow(alpha,2);
  double c2 = 1./pow(alpha,2);
  double za2 = pow(alpha*z,2);
  double g=sqrt(k*k-za2);

  double w2 = pow(z,2)/pow(g+n-fabs(k),2);
  double d  = 1.+a2*w2;

  double diracE = -1*w2/(2*d) - (a2*w2/2+1-sqrt(1+a2*w2))*(c2/d);
  //double diracE=c2*(1./sqrt(1+za2/pow((g+n-fabs(k)),2))-1.);
	return diracE;
}




//******************************************************************************
int ElectronOrbitals::zeroNucleus()
/*
infinitesimal nucleus.
1/r potential
*/
{
  vnuc.clear(); //just to be sure..
  vnuc.push_back(0); //XXX ??
  for(int i=1; i<ngp; i++)
    vnuc.push_back(-Z/r[i]);
  return 0;
}

//******************************************************************************
int ElectronOrbitals::sphericalNucleus(double rnuc)
/*
Potential due to a spherical nucleus, with (charge) radius, rnuc.
Note: rnuc must be given in "fermi" (fm, femto-metres).
If rnuc=0 is given [which is default], then will use approximate
formula for rN.
See: https://www-nds.iaea.org/radii/
*/
{
  vnuc.clear(); //just to be sure..

  double rN=rnuc; //nuclear charge radius:
  //Estimate nuclear charge radius. Only for spherical nuclei.
  if(rnuc==0){
    if(A==1) rN = 0.8783;       // 1-H
    else if(A==4) rN = 1.6755;  // 4-He
    else if(A==7) rN = 2.4440;  // 7-Li
    else if(A<10) rN = 1.15*pow(A,0.333);
    else rN = (0.836*pow(A,0.333)+0.570);
  }
  rN/=ABOHR_FM;
  // XXX Add data tables of nuclear radii!

  //Fill the vnuc array with spherical nuclear potantial
  vnuc.push_back(0); //XXX ??
  double rn2=pow(rN,2);
  double rn3=pow(rN,3);
  for(int i=1; i<ngp; i++){
    double temp_v;
    double ri = r[i];
    if(ri<rN) temp_v = Z*(pow(ri,2)-3.*rn2)/(2.*rn3); //XXX 2.*rn3?? check!? XXX
    else temp_v = -Z/ri;
    vnuc.push_back(temp_v);
    //std::cout<<i<<" "<<r[i]<<" "<<rN<<" "<<temp_v<<" "<<-Z/ri<<"\n";
  }

  return 0;
}

//******************************************************************************
int ElectronOrbitals::fermiNucleus(double t, double c)
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
  vnuc.clear(); // clear the array [just in case..]

  if(t==0) t=2.4; // Default skin-thickness (in fm)
  if(c==0) c=1.1*pow(A,0.3333); //default half-charge radius ????
  // XXX Better approx! +/or data tables!

  //std::cout<<c/ABOHR_FM<<"\n";

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
    //std::cout<<i<<" "<<r[i]<<" "<<t_v<<" "<<-Z/t_r<<"\n";
  }

  return 0;
}
