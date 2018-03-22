//class ElectronOrbitals::
#include "ElectronOrbitals.h"


//XXX XXX XXX make routine that solves for a signle state. pushes back n,k etc.!

//******************************************************************************
ElectronOrbitals::ElectronOrbitals(int in_z, int in_a, int in_ngp, double rmin,
  double rmax, double var_alpha)
{

  ngp=in_ngp;
  //JohnsonRadialGrid(rmin,rmax);
  DzubaRadialGrid(rmin,rmax);

  alpha=ALPHA*var_alpha;

  Z=in_z;
  if(in_a==0) A=ATI_a[Z]; //Use default atomic mass
  else A=in_a;
  zeroNucleus();
}
//-----Overloaded---------------------------------------------------------------
ElectronOrbitals::ElectronOrbitals(std::string s_in_z, int in_a, int in_ngp,
  double rmin, double rmax, double var_alpha)
{
 //Work out Z from given atomic symbol
 int iz = ATI_get_z(s_in_z);
 ElectronOrbitals(iz,in_a,in_ngp,var_alpha);
}



//******************************************************************************
int ElectronOrbitals::solveLocalDirac(int n, int k, double e_a)
/*
XXX vnuc !!!
*/
{
  int pinf,its;
  double eps;
  std::vector<double> p_a(ngp);
  std::vector<double> q_a(ngp);

  std::vector<double> v_a = vnuc;
  if(vdir.size()!=0){
    for(int i=0; i<ngp; i++) v_a[i] += vdir[i];
  }

  int i_ret = solveDBS(p_a,q_a,e_a,v_a,Z,n,k,r,drdt,h,ngp,pinf,its,eps,alpha);
  //Store wf + energy
  p.push_back(p_a);
  q.push_back(q_a);
  en.push_back(e_a);
  //Store states:
  nlist.push_back(n);
  klist.push_back(k);
  //store convergance info:
  pinflist.push_back(pinf);
  itslist.push_back(its);
  epslist.push_back(eps);

  return i_ret;
}

//******************************************************************************
int ElectronOrbitals::hydrogenLike(int in_max_n, int in_max_l)
/*
``Wrapper'' function, to use the adamsSolveLocalBS method!
Set up specifically for H-like ions
XXX Kill this one! XXX
*/
{
  max_n = in_max_n;
  max_l = in_max_l;

  for(int n=1; n<=max_n; n++){
    for(int i=1; i<2*n; i++){ //loop through each kappa state
      int k = pow(-1,i)*ceil(0.5*i);
      int l = (abs(2*k+1)-1)/2;
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
int ElectronOrbitals::determineCore(std::vector<std::string> str_core,
  std::vector<int> &core)
/*
Takes in a string list for the core configuration, outputs an int list
Takes in previous closed shell (noble), + 'rest' (or just the rest)
E.g:
  Core of Cs: Xe
  Core of Gold: Xe 4f14 5d10
'rest' is in form nLm : n=n, L=l, m=number of electrons in that nl shell.
NOTE: Only works up to n=9, and l=5 [h]
*/
{
  core.clear();
  if(str_core.size()==0) return 1;

  int ibeg=1;
  std::string ng=str_core[0];
  if     (ng=="He") core=ATI_core_He;
  else if(ng=="Ne") core=ATI_core_Ne;
  else if(ng=="Ar") core=ATI_core_Ar;
  else if(ng=="Kr") core=ATI_core_Kr;
  else if(ng=="Xe") core=ATI_core_Xe;
  else if(ng=="Rn") core=ATI_core_Rn;
  else if(ng=="Og") core=ATI_core_Og;
  else ibeg=0;

  std::vector<int> core_ex;
  for(size_t i=ibeg; i<str_core.size(); i++){
    //Parse string, determine config for this term
    int n = std::stoi(str_core[i].substr(0,1));
    int m = std::stoi(str_core[i].substr(2));
    std::string strl=str_core[i].substr(1,1);
    int l=-1;
    if     (strl=="s") l=0;
    else if(strl=="p") l=1;
    else if(strl=="d") l=2;
    else if(strl=="f") l=3;
    else if(strl=="g") l=4;
    else if(strl=="h") l=5;
    else return 2;

    //Check if this term is valid
    if(m>4*l+2) return 2;
    if(l+1>n)   return 2;

    //Form int list for this term:
    for(int in=0; in<=n; in++){
      for(int il=0; il<in; il++){
        if(in==n && il==l) core_ex.push_back(m);
        else               core_ex.push_back(0);
      }
    }

    //Merge this term with the existing core:
    int size = std::max(core_ex.size(),core.size());
    core_ex.resize(size);
    core.resize(size);
    for(size_t j=0; j<core.size(); j++){
      core[j] += core_ex[j];
      if(core[j] > 4*ATI_core_l[j]+2) return 2; //check if valid
    }

  }



  return 0;
}

//******************************************************************************
int ElectronOrbitals::JohnsonRadialGrid(double r0, double rmax)
/*
Non-uniform r grid, taken from Johnson book (w/ minor modification).
Quite a standard exponential grid.
Uses:
   dr/dt = r0 * exp(t)
=> r = r0 * exp(t)
   t = i*h for i=0,1,2,...
*/
{

  h=log(rmax/r0)/(ngp-1);
  if(h>0.03){
    std::cout<<"Warning: h="<<h<<" in formRadialGrid. Is this too large??\n";
  }

  drdt.clear();
  for(int i=0; i<ngp; i++){
    double temp_drdt = r0*exp(i*h);
    drdt.push_back(temp_drdt);
  }

  r.clear();
  for(int i=0; i<ngp; i++){
    double temp_r = drdt[i];
    r.push_back(temp_r);
  }

  // (dr/dt)/r [for convinience]
  //Note: with current radial grid, this is always 1! (not true in general)
  dror.clear();
  dror.push_back(1.);
  for(int i=1; i<ngp; i++){
    double temp_dror = drdt[i]/r[i];
    dror.push_back(temp_dror);
  }

  return 0;
}

//******************************************************************************
int ElectronOrbitals::DzubaRadialGrid(double r0, double rmax, double b)
/*
Non-uniform r-grid taken from V. Dzuba code.
Uses:
  dr/dt = r / (b + r)
  t = r + b*log(r)
Grid will be ~ exponential below 'b', roughly linear afterwards
Default b=4.
*/
{

  h=(rmax-r0+b*log(rmax/r0))/(ngp-1);

  //clear the vectors, just in case
  r.clear();
  drdt.clear();
  dror.clear();
  //initial points:
  r.push_back(r0);
  dror.push_back(1./(b+r0));
  drdt.push_back(dror[0]*r0);

  // Use method from Dzuba code to calculate r grid
  double t = r0 + b*log(r0); //t is linear/uniform grid
  for(int i=1; i<ngp; i++){
    t+=h;
    double t_r = r[i-1];  //"temp" r
    //Integrate dr/dt to find r:
    double delta_r=1.;
    int ii=0; //to count number of iterations
    while(delta_r>r0*1.e-15){
      double delta_t = t - (t_r + b*log(t_r)); // t = t0+i*h
      double t_drdt = t_r/(t_r + b); // temp dr/dt
      delta_r = delta_t * t_drdt;
      t_r += delta_r;
      ii++;
      if(ii>30) break; //usually converges in ~ 2 or 3 steps!
    }
    //std::cout<<ii<<" "<<eps<<"\n";
    r.push_back(t_r);
    dror.push_back(1./(b+t_r));
    drdt.push_back(dror[i]*t_r);
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

NOTE: Only seems to work for fairly large Z !

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
