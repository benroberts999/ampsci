//class ElectronOrbitals::
#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "FPC_physicalConstants.h"
#include "ATI_atomInfo.h"
#include "ADAMS_solveLocalBS.h"
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <gsl/gsl_sf_fermi_dirac.h> //For fermiNucleus
#include <algorithm> //for sort
/*
 ==== To Do ====
 * Write out to disk

*/

//******************************************************************************
ElectronOrbitals::ElectronOrbitals(int in_z, int in_a, int in_ngp, double rmin,
  double rmax, double var_alpha)
{
  clearAll();
  logLinearRadialGrid(in_ngp,rmin,rmax);

  alpha=FPC::alpha*var_alpha;
  num_core_states=0;

  //XXX update this! Make spherical default.
  Z=in_z;
  if(in_a==0) A=ATI::A[Z]; //Use default atomic mass
  else A=in_a;
  zeroNucleus();
}
//-----Overloaded---------------------------------------------------------------
ElectronOrbitals::ElectronOrbitals(std::string s_in_z, int in_a, int in_ngp,
  double rmin, double rmax, double var_alpha)
{
 //Work out Z from given atomic symbol
 int iz = ATI::get_z(s_in_z);
 ElectronOrbitals(iz,in_a,in_ngp,rmin,rmax,var_alpha);
}


//******************************************************************************
int ElectronOrbitals::solveLocalDirac(int n, int k, double e_a, int log_dele_or)
/*
Uses ADAMS::solveDBS to solve Dirac Eqn for local potential (Vnuc + Vdir)
*/
{
  int pinf,its;
  double eps;
  std::vector<double> f_a(ngp);
  std::vector<double> g_a(ngp);

  //Fill V(r) with nulcear + DIRECT part of electron potential
  //nb: for exchange part, need to use reSolveDirac()
  std::vector<double> v_a = vnuc;
  if(vdir.size()!=0) for(int i=0; i<ngp; i++) v_a[i] += vdir[i];

  //Solve local dirac Eq:
  if(e_a==0) e_a = enGuessVal(n,k);
  int i_ret = ADAMS::solveDBS(f_a,g_a,e_a,v_a,Z,n,k,r,drdt,h,ngp,pinf,its,eps,
    alpha,log_dele_or);
  //Store wf + energy
  f.push_back(f_a);
  g.push_back(g_a);
  en.push_back(e_a);
  //Store states:
  nlist.push_back(n);
  kappa.push_back(k);
  //store convergance info:
  pinflist.push_back(pinf);
  itslist.push_back(its);
  epslist.push_back(eps);

  return i_ret;
}

//******************************************************************************
double ElectronOrbitals::get_alpha()const{
  return alpha;
}

//******************************************************************************
int ElectronOrbitals::reSolveDirac(int i, double e_a, int log_dele_or)
/*Overloaded version; see below
This one doesn't have exchange potential
*/
{
  std::vector<double> dummy_vex;
  return reSolveDirac(i,e_a,dummy_vex,log_dele_or);
}
//******************************************************************************
int ElectronOrbitals::reSolveDirac(int i, double e_a,
  const std::vector<double> &vex, int log_dele_or)
/*
"Re"solves dirac eqaution. Use this to re-solve for same state.
Over-rides existing solution.
If no e_a is given, will use the existing one!
(Usually, a better guess should be given, using P.T.)
Note: optionally takes in exchange potential! (see overloaded above)
Note: Uses the "dodgy" re-scaled exchange potenital:
Vex\psi_a = sum_b vex_a psi_b -> [sum_b vex_a (psi_b/psi_a)] psi_a
so, here, vex = [sum_b vex_a (psi_b/psi_a)]
This is not ideal..
*/
{
  int pinf,its;
  double eps;

  //XXX this is inneficient. Fine, except for HF. THen, slow! ?
  std::vector<double> v_a = vnuc;
  if(vdir.size()!=0) for(int i=0; i<ngp; i++) v_a[i] += vdir[i];
  if(vex.size()!=0) for(int i=0; i<ngp; i++) v_a[i] += vex[i];

  int n=nlist[i];
  int k=kappa[i];
  if(e_a==0) e_a = en[i];
  int i_ret = ADAMS::solveDBS(f[i],g[i],e_a,v_a,Z,n,k,r,drdt,h,ngp,pinf,its,eps,
    alpha,log_dele_or);
  en[i] = e_a; //update w/ new energy

  //store convergance info:
  pinflist[i] = pinf;
  itslist[i]  = its;
  epslist[i]  = eps;

  return i_ret;
}


//******************************************************************************
int ElectronOrbitals::determineCore(std::string str_core_in)
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
  //Move comma-seperated string into an array (vector)
  std::vector<std::string> str_core;
  std::stringstream ss(str_core_in);
  while( ss.good() ){
    std::string substr;
    getline( ss, substr, ',' );
    str_core.push_back( substr );
  }

  num_core_shell.clear();
  if(str_core.size()==0) return 1;

  int ibeg=1;
  std::string ng=str_core[0];
  if     (ng=="He") num_core_shell=ATI::core_He;
  else if(ng=="Ne") num_core_shell=ATI::core_Ne;
  else if(ng=="Ar") num_core_shell=ATI::core_Ar;
  else if(ng=="Kr") num_core_shell=ATI::core_Kr;
  else if(ng=="Xe") num_core_shell=ATI::core_Xe;
  else if(ng=="Rn") num_core_shell=ATI::core_Rn;
  else if(ng=="Og") num_core_shell=ATI::core_Og;
  else if(ng=="Zn") num_core_shell=ATI::core_Zn;
  else if(ng=="Cd") num_core_shell=ATI::core_Cd;
  else if(ng=="Hg") num_core_shell=ATI::core_Hg;
  else ibeg=0;

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

    std::vector<int> core_ex;
    //Form int list for this term:
    for(int in=0; in<=n; in++){
      for(int il=0; il<in; il++){
        if(in==n && il==l) core_ex.push_back(m);
        else               core_ex.push_back(0);
      }
    }

    //Merge this term with the existing core:
    int size = std::max(core_ex.size(),num_core_shell.size());
    core_ex.resize(size);
    num_core_shell.resize(size);

    for(size_t j=0; j<num_core_shell.size(); j++){
      num_core_shell[j] += core_ex[j];
      if(num_core_shell[j] > 4*ATI::core_l[j]+2) return 2; //check if valid
    }

  }

  //Count number of electrons in the core
  num_core_electrons = 0;
  for(int &num : num_core_shell) num_core_electrons += num;

  return 0;
}

//******************************************************************************
bool ElectronOrbitals::isInCore(int n, int k)
/*
Checks if given state is in the core.
NOTE: in some cases, given state may be in and out! Account for this?
XXX Also check occupancy fraction? Do seperately!
*/
{
  for(int i=0; i<num_core_states; i++)
    if(n==nlist[i] && k==kappa[i]) return true;
  return false;
}

//******************************************************************************
int ElectronOrbitals::maxCore_n(void)
//Returns the largest n in the core (used for energy guesses)
{
  int max_n = 0;
  for(int i=0; i<num_core_states; i++)
    if(nlist[i]>max_n) max_n = nlist[i];
  return max_n;
}

//******************************************************************************
int ElectronOrbitals::solveInitialCore(int log_dele_or)
/*
Solves the Dirac eqn for each state in the core
Only for local potential (direct part)
HF_hartreeFock.cpp has routines for Hartree Fock
*/
{
  for(size_t i=0; i<num_core_shell.size(); i++){
    int num = num_core_shell[i];
    if(num==0) continue;
    int n = ATI::core_n[i];
    int l = ATI::core_l[i];
    double en_a = enGuessCore(n,l);
    int k1 = l; //j = l-1/2
    if(k1!=0) {
      solveLocalDirac(n,k1,en_a,log_dele_or);
      en_a = 0.95*en[nlist.size()-1]; //update guess for next same l
      if(en_a>0) en_a = enGuessCore(n,l);
    }
    int k2 = -(l+1); //j=l+1/2
    solveLocalDirac(n,k2,en_a,log_dele_or);
  }
  num_core_states = nlist.size(); //store number of states in core

  //occupancy fraction for each core state (avg of Non-rel states!):
  for(int i=0; i<num_core_states; i++){
    int n = nlist[i];
    int ka = kappa[i];
    int l = ATI::l_k(ka);
    //Find the correct core list index (to determine filling factor):
    int ic=-1;
    for(size_t j=0; j<num_core_shell.size(); j++){
      if(n==ATI::core_n[j] && l==ATI::core_l[j]){
        ic = j;
        break;
      }
    }
    if(ic==-1){
      std::cout<<"FAIL 254 in ElectronOrbitals:solveInitialCore\n";
      return 2;
    }
    occ_frac.push_back(double(num_core_shell[ic])/(4*l+2));
  }

  return 0;
}

//******************************************************************************
void ElectronOrbitals::orthonormaliseOrbitals(int num_its)
/*
Forces ALL orbitals to be orthogonal to each other, and normal
Note: workes best if run twice!
|a> ->  |a> - \sum_{b!=a} |b><b|a>
Then:
|a> -> |a> / <a|a>
c_ba = c_ab = <a|b>
num_its is optional parameter. Repeats that many times!
Note: I force all orthog to each other - i.e. double count.
{force <2|1>=0 and then <1|2>=0}
Would be 2x faster not to do this - but that would treat some orbitals special!
Hence factor of 0.5
*/
{
  int Ns = nlist.size();
  std::vector< std::vector<double> > c_ab(Ns, std::vector<double>(Ns));

  //Calculate c_ab = <a|b>  [for b>a -- symmetric]
  for(int a=0; a<Ns; a++){
    for(int b=a+1; b<Ns; b++){
      if(kappa[a]!=kappa[b]) continue;
      double fab = INT::integrate3(f[a],f[b],drdt);
      double gab = INT::integrate3(g[a],g[b],drdt);
      c_ab[a][b] = 0.5*h*(fab+gab); //0.5 avoids double counting
    }
  }

  //fill in the a>b part:
  for(int b=0; b<Ns; b++){
    for(int a=b+1; a<Ns; a++){
      if(kappa[a]!=kappa[b]) continue;
      c_ab[a][b] = c_ab[b][a];
    }
  }

  //Orthogonalise orbitals:
  for(int a=0; a<Ns; a++){
    for(int b=0; b<Ns; b++){
      if(a==b) continue;
      if(kappa[a]!=kappa[b]) continue;
      double cab = c_ab[a][b];
      if(cab==0) continue; //already?
      for(int ir=0; ir<ngp; ir++){
        f[a][ir] -= cab*f[b][ir];
        g[a][ir] -= cab*g[b][ir];
      }
    }
  }

  //Re-normalise orbitals (nb: doesn't make much difference)
  for(int a=0; a<Ns; a++){
    double faa = INT::integrate3(f[a],f[a],drdt);
    double gaa = INT::integrate3(g[a],g[a],drdt);
    double norm = 1./sqrt(h*(faa+gaa));
    for(int ir=0; ir<ngp; ir++){
      f[a][ir] *= norm;
      g[a][ir] *= norm;
    }
  }

  //If necisary: repeat
  if(num_its>1) orthonormaliseOrbitals(num_its-1);
}

//******************************************************************************
void ElectronOrbitals::orthonormaliseValence(int num_its)
/*
Force valence orbitals to be orthogonal to:
  a) core orbitals
  b) other valence orbitals
After the core is 'frozen', don't touch core orbitals!
|v> --> |v> - sum_c |c><c|v> - sum_{w!=v} |w><w|v>
*/
{
  int Ns_T = nlist.size();
  int Ns_c = num_core_states;
  int Ns_v = Ns_T - Ns_c;

  //Calculate the core-valence coeficients <c|v> = A_cv
  std::vector< std::vector<double> > A_vc(Ns_v, std::vector<double>(Ns_c));
  for(int iv=Ns_c; iv<Ns_T; iv++){
    for(int ic=0; ic<Ns_c; ic++){
      if(kappa[iv]!=kappa[ic]) continue;
      double fvc = INT::integrate3(f[iv],f[ic],drdt);
      double gvc = INT::integrate3(g[iv],g[ic],drdt);
      A_vc[iv-Ns_c][ic] = h*(fvc+gvc); //no 0.5 here - no double counting
    }
  }

  //Calculate the valence-valence coefs, <v|w> [only for w>v]
  std::vector< std::vector<double> > A_vw(Ns_v, std::vector<double>(Ns_v));
  for(int iv=Ns_c; iv<Ns_T; iv++){
    for(int iw=iv+1; iw<Ns_T; iw++){
      if(kappa[iv]!=kappa[iw]) continue;
      double fvw = INT::integrate3(f[iv],f[iw],drdt);
      double gvw = INT::integrate3(g[iv],g[iw],drdt);
      A_vw[iv-Ns_c][iw-Ns_c] = 0.5*h*(fvw+gvw); //0.5 due to double counting
    }
  }
  //fill in the symmetric v>w part:
  for(int iw=Ns_c; iw<Ns_T; iw++){
    for(int iv=iw+1; iv<Ns_T; iv++){
      if(kappa[iv]!=kappa[iw]) continue;
      A_vw[iv-Ns_c][iw-Ns_c] = A_vw[iw-Ns_c][iv-Ns_c];
    }
  }

  //Orthogonalise:
  for(int iv=Ns_c; iv<Ns_T; iv++){
    //Core-valence part:
    for(int ic=0; ic<Ns_c; ic++){
      if(kappa[iv]!=kappa[ic]) continue;
      double Avc = A_vc[iv-Ns_c][ic];
      if(Avc==0) continue; //never?
      for(int ir=0; ir<ngp; ir++){
        f[iv][ir] -= Avc*f[ic][ir];
        g[iv][ir] -= Avc*g[ic][ir];
      }
    }
    //valence-valence part:
    for(int iw=Ns_c; iw<Ns_T; iw++){
      if(iv==iw) continue;
      if(kappa[iv]!=kappa[iw]) continue;
      double Avw = A_vw[iv-Ns_c][iw-Ns_c];
      if(Avw==0) continue; //never?
      for(int ir=0; ir<ngp; ir++){
        f[iv][ir] -= Avw*f[iw][ir];
        g[iv][ir] -= Avw*g[iw][ir];
      }
    }
  }

  //Re-normalise the orbitals:
  for(int iv=Ns_c; iv<Ns_T; iv++){
    double fvv = INT::integrate3(f[iv],f[iv],drdt);
    double gvv = INT::integrate3(g[iv],g[iv],drdt);
    double norm = 1./sqrt(h*(fvv+gvv));
    for(int ir=0; ir<ngp; ir++){
      f[iv][ir] *= norm;
      g[iv][ir] *= norm;
    }
  }

  //If necisary: repeat
  if(num_its>1) orthonormaliseValence(num_its-1);
}

//******************************************************************************
double ElectronOrbitals::enGuessCore(int n, int l)
/*
Private
Energy guess for core states. Not perfect, good enough
tot_el = total electrons BELOW
num = num electrons in THIS shell
*/
{

  int tot_el = 0;
  int num = 0;
  for(size_t i=0; i<num_core_shell.size(); i++){
    if(l==ATI::core_l[i] && n==ATI::core_n[i]){
      num = num_core_shell[i];
      break;
    }
    tot_el += num_core_shell[i];
  }

  //effective Z (for energy guess) -- not perfect!
  double Zeff =  double(Z - tot_el - num);
  if(l==1) Zeff = 1. + double(Z - tot_el - 0.5*num);
  if(l==2) Zeff = 1. + double(Z - tot_el - 0.5*num);
  if(Zeff<1.) Zeff=1.;

  double en_a = -0.5 * pow(Zeff/n,2);
  if(n>1) en_a *= 0.5;

  if(n==maxCore_n()-1) en_a *= 1.25;
  if(n==maxCore_n()) en_a *= 2.5;

  return en_a;
}


//******************************************************************************
double ElectronOrbitals::enGuessVal(int n, int ka)
/*Energy guess for valence states. Not perfect, good enough*/
{
  int maxn = maxCore_n();
  int l = ATI::l_k(ka);
  int dn=n-maxn;
  double neff=1.+dn;
  double x=1;
  if(maxn<4) x=0.25;
  if(l==1) neff+=0.5*x;
  if(l==2) neff+=2.*pow(x,0.5);
  if(l>=3) neff+=4.*x;
  return -0.5/pow(neff,2);
}


//******************************************************************************
int ElectronOrbitals::exponentialRadialGrid(int in_ngp, double r0, double rmax)
/*
Non-uniform r grid, taken from Johnson book (w/ minor modification).
Quite a standard exponential grid.
Uses:
  dr/dt = r0 * exp(t)
  =>  r = r0 * exp(t)
      t = i*h for i=0,1,2,...
*/
{

  ngp = in_ngp;
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
int ElectronOrbitals::getRadialIndex(double r_target)
/*
Finds the radial grid index that corresponds to r=r_target
NOTE: returns index that corresponds to r _lower_ that (or equal to) r_target
*/
{
  if (r.size()==0){
    std::cout<<"ERROR 219 in ElectronOrbitals - no grid!\n";
    return -1;
  }
  if(r_target >= r[ngp-1]) return (ngp-1);
  if(r_target <= r[0]) return 0; //nb: in this case: r[i] > r_target! careful!
  //loop through, find i for r[i] <= r_target < r[i+1]:
  for(size_t i=0; i<r.size()-1; i++){
    if(r_target>=r[i] && r_target<r[i+1]) return i;
  }
  //shouldn't get past this loop! This just for safety?
  std::cout<<"ERROR 227 in ElectronOrbitals - didn't find r??\n";
  return -1; //??
}

//******************************************************************************
int ElectronOrbitals::logLinearRadialGrid(double in_h, double r0, double rmax,
   double b)
/*
OVERLOADED.
So, can give an input h, and let program determine NGP.
*/
{
  double d_ngp = (in_h-r0+rmax+b*log(rmax/r0))/in_h;
  int in_ngp = (int)d_ngp+1;
  return logLinearRadialGrid(in_ngp,r0,rmax,b);
}
//******************************************************************************
int ElectronOrbitals::logLinearRadialGrid(int in_ngp, double r0, double rmax,
   double b)
/*
Non-uniform r-grid taken from V. Dzuba code.
Uses:
  dr/dt = r / (b + r)
  t = r + b*log(r)
Grid will be ~ exponential below 'b', roughly linear afterwards
Default b=4.
*/
{

  ngp = in_ngp;
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
double ElectronOrbitals::diracen(double z, double n, int k){
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
  rN/=FPC::aB_fm;
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

rho(r) = rho_0 {1 + Exf[(r-c)/a]}^-1
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

  //std::cout<<c/FPC::aB_fm<<"\n";

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
      double roa = FPC::aB_fm*t_r/a; // convert fm <-> atomic
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


//------------------------------------------------------------------------------
int ElectronOrbitals::sortedEnergyList(std::vector<int> &sort_list)
/*
Outouts a list of integers corresponding to the states
sorted by energy (lowest energy first)
*/
{
  std::vector< std::vector<double> > t_en;
  for(size_t i=0; i<en.size(); i++){
    t_en.push_back({en[i],(double)i+0.1});
    //+0.1 to prevent rounding error when going from double -> int
  }

  //Function pointer: sort 2D vector by first col!
  bool (*sortCol)(const std::vector<double>& v1,
    const std::vector<double>& v2)
    = [](const std::vector<double>& v1, const std::vector<double>& v2)
  {
        return v1[0] > v2[0];
  };

  std::sort(t_en.rbegin(), t_en.rend(), sortCol);

  for(size_t i=0; i<en.size(); i++){
    sort_list.push_back((int)t_en[i][1]);
  }
  return 0;
}

//******************************************************************************
void ElectronOrbitals::clearAll()
/*
Just clears all the info.
*/
{
  //Atom info
  Z=0;
  A=0;
  atom="0";
  //orbitals:
  f.clear();
  g.clear();
  en.clear();
  //state info:
  nlist.clear();
  kappa.clear();
  //info from solveing DE
  pinflist.clear();
  itslist.clear();
  epslist.clear();
  //grid
  r.clear();
  drdt.clear();
  dror.clear();
  ngp=0;
  h=0;
  //Potentials
  vnuc.clear();
  vdir.clear();
  //core:
  occ_frac.clear();
  num_core_states=0;
  num_core_electrons=0;
  // ** private **
  alpha=0;
  num_core_shell.clear();
}
