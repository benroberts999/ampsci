#include "HF_hartreeFock.h"
#include "ElectronOrbitals.h"
#include "PRM_parametricPotentials.h"
#include "ATI_atomInfo.h"
#include "WIG_369j.h"
#include <cmath>
#include <string>
#include <vector>

/*
Calculates self-consistent Hartree-Fock potential, including exchange.
Solves all core and valence states.
*/
namespace HF{

//******************************************************************************
int hartreeFockCore(ElectronOrbitals &wf, double eps_HF)
/*
Calculates the Hartree-Fock potential (w/ exchange) and core orbitals.
Uses a Greens parametric potential as starting approx.
Then calculated V_dir and V_exch + solves all core electrons
Note: V_HF = V_dir + V_ex    -- note sign convention!
*/
{
  // "Mixing" of new + old Potential:
  // V_n+1 = eta*V_n+1 + (1-eta)*V_n
  //NB: must be less than 0.5!, because x2 after first few its!
  double eta=0.35;
  int ngp = wf.ngp;

  wf.vdir.clear(); //make sure it's empty

  //Starting approx:
  //Fill the electron part of the potential, using Greens PRM for initial approx
  double Gh,Gd;  //Green potential parameters
  PRM::defaultGreenCore(wf.Z,Gh,Gd); //Get default values for Green potential
  for(int i=0; i<ngp; i++) wf.vdir.push_back(PRM::green(wf.Z,wf.r[i],Gh,Gd));

  //First step: Solve each core state using above parameteric potential
  wf.solveInitialCore(3); //3, since don't need high accuray here [1 in 10^3]
  int Ncs = wf.num_core_states;

  //Define exchange potential. Note: not stored at the moment!
  std::vector< std::vector<double> > vex; //into class??
  vex.clear(); //move to beginning ?
  vex.resize(Ncs,std::vector<double>(ngp));

  //Start the HF itterative procedure:
  int hits;
  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta *= 2;

    //Form new v_dir and v_ex:
    //XXX XXX XXX Fastor to form these once outside! ? Check if safe! XXX
    std::vector<double> vdir_old = wf.vdir;
    //std::vector<double> vdir_new;
    formNewVdir(wf,wf.vdir,false);
    std::vector< std::vector<double> > vex_old = vex;
    formVexCore(wf,vex);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*wf.vdir[j] + (1.-eta)*vdir_old[j];
      for(int i=0; i<Ncs; i++){
        vex[i][j] = eta*vex[i][j] + (1.-eta)*vex_old[i][j];
      }
    }

    //Solve Dirac Eq. for each state in core, using Vdir+Vex:
    std::vector<double> en_old = wf.en;
    double t_eps = 0;
    for(int i=0; i<Ncs; i++){
      //calculate de from PT
      double del_e = 0;
      for(int j=0; j<wf.ngp; j++){
        double dv = wf.vdir[j]+vex[i][j]-vdir_old[j]-vex_old[i][j];
        del_e += dv*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
      }
      del_e*=wf.h;
      double en_guess = en_old[i] + del_e;
      if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
      wf.reSolveDirac(i,en_guess,vex[i],3); //only go to 1/10^3 here
      //t_eps: weighted average of (de)/e for each orbital:
      double sfac = 2.*wf.kappa[i]*wf.core_ocf[i]; //|2k|=2j+1
      t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
    }
    t_eps /= (wf.num_core_electrons*eta);

    //Force all core orbitals to be orthogonal to each other
    wf.orthonormaliseOrbitals(1);

    printf("\rHF core        it:%3i eps=%6.1e              ",hits,t_eps);
    std::cout<<std::flush;
    if(t_eps<eps_HF && hits>2) break;
  }
  std::cout<<"\n";

  //Now, re-solve core orbitals with higher precission
  // + re-solve direct potential (higher precission)
  for(int i=0; i<Ncs; i++) wf.reSolveDirac(i,wf.en[i],vex[i],15);
  wf.orthonormaliseOrbitals(1);
  formNewVdir(wf,wf.vdir,false); //+ new v_ex??
  for(int i=0; i<Ncs; i++) wf.reSolveDirac(i,wf.en[i],vex[i],15);
  wf.orthonormaliseOrbitals(2);

  return 0;
}

//******************************************************************************
int hartreeFockValence(ElectronOrbitals &wf, int na, int ka, double eps_HF)
/*
Calculate valence states in frozen Hartree-Fock core
*/
{
  //Mixing parameter (fraction of new pot to use)
  //note: is x2 after 4 its!
  double eta=0.35; //must be <0.5!
  //Get initial value, with no exchange:
  double en_g = wf.enGuessVal(na,ka);
  wf.solveLocalDirac(na,ka,en_g,3);
  int a = wf.nlist.size() - 1;

  std::vector<double> vexa(wf.ngp);
  int hits;
  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta*=2;
    std::vector<double> vexa_old = vexa;
    double en_old = wf.en[a];
    //Form new exch. potential:
    std::vector<double> vexa_new;
    formVexA(wf,a,vexa_new);
    for(int i=0; i<wf.ngp; i++) vexa[i] = eta*vexa_new[i]+(1.-eta)*vexa_old[i];
    //Use P.T. to calculate energy change:
    double en_new = 0;
    for(int i=0; i<wf.pinflist[a]; i++) en_new +=
      (vexa[i]-vexa_old[i])*(pow(wf.p[a][i],2)+pow(wf.q[a][i],2))*wf.drdt[i];
    en_new = wf.en[a] + en_new*wf.h;
    //Solve Dirac using new potential:
    wf.reSolveDirac(a,en_new,vexa,3);
    double eps = fabs((wf.en[a]-en_old)/(eta*en_old));
    //Force valence states to be orthogonal to each other + to core:
    wf.orthonormaliseValence(1);
    printf("\rHF val:%2i %2i %2i | %3i eps=%6.1e  en=%11.8f                  "
    ,a,na,ka,hits,eps,wf.en[a]);
    std::cout<<std::flush;
    if(eps<eps_HF && hits>2) break;
  }
  std::cout<<"\n";
  //Re-solve w/ higher precission (probs not needed)
  wf.reSolveDirac(a,wf.en[a],vexa,15);
  wf.orthonormaliseValence(2);
  return hits;
}

//******************************************************************************
int formNewVdir(const ElectronOrbitals &wf, std::vector<double> &vdir_new,
  bool scale)
/*
This takes in the wavefunctions, and forms the direct (local) part of the
electron potential. Does not include exchange.
If scale==true, scales by (Nc-1)/Nc [accounts for self-interaction]
scale should be true if not doing exchange, but false if doing exchange!
Note: Uses efficient integral method:
  vdir(r) :=  Int_0^inf rho(r')/r_max dr'
           =  Int_0^r rho(r')/r dr' + Int_r^inf rho(r')/r' dr'
           =  A(r)/r + B(r)
  A(r0)    =  0
  B(r0)    =  Int_0^inf rho(r')/r' dr'
  A(r_n)   =  A(r_n-1) + rho(r)dr
  B(r_n)   =  B(r_n-1) - (rho(r)/r)dr
  vdir(r)  =  A(r)/r + B(r)
*/
{

  //Make sure vector is correct size, and clear old potential away
  //vdir_new.clear();
  //vdir_new.resize(wf.ngp);

  //Scaling factor. For hartree only (when FOCK excluded)
  //Because I sum over all electrons, includes self-interaction
  //This part is cancelled from FOCK, but if no fock, needs to be scaled down.
  double f=1;
  if(scale) f = 1. - (1.)/wf.num_core_electrons;

  //Determine the total electron (charge) density of core:
  // 2j+1 is closed-shell occupation number of that orbital.
  std::vector<double> rho(wf.ngp,0);
  for(int i=0; i<wf.num_core_states; i++){
    int ka = wf.kappa[i];
    int twoj = ATI::twoj_k(ka);
    double frac = wf.core_ocf[i]; //avgs over non-rel. configs
    int j_max = wf.pinflist[i];
    for(int j=0; j<j_max; j++){
      rho[j] += frac*(twoj+1)*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2));
    }
  }

  //Use the above determined electron (charge) density with Gauss' law
  //to determine the electric potential (energy). Makes use of spherical symm.
  // vdir(r) = A(r)/r + B(r) [see above]
  //Note: I don't use quadrature formula here.
  double a=0, b = 0;
  for(int i=0; i<wf.ngp; i++) b += wf.drdt[i]*rho[i]/wf.r[i];
  vdir_new[0] = b;
  for(int i=1; i<wf.ngp; i++){
    b = b - wf.drdt[i-1]*rho[i-1]/wf.r[i-1];
    a = a + wf.drdt[i-1]*rho[i-1];
    vdir_new[i] = f*wf.h*(a/wf.r[i] + b);
  }

  return 0;
}

//==============================================================================

//******************************************************************************
int formVexCore(const ElectronOrbitals &wf, std::vector< std::vector<double> > &vex)
/*
Calculates exchange potential for each state in the core
Parallelised for speed
*/
{
  //vex.clear();
  vex.resize(wf.num_core_states);
  #pragma omp parallel for
  for(int a=0; a<wf.num_core_states; a++){
    std::vector<double> vex_a;
    formVexA(wf,a,vex_a);
    vex[a] = vex_a;
  }
  return 0;
}

//******************************************************************************
int formVexA(const ElectronOrbitals &wf, int a, std::vector<double> &vex_a)
/*
Forms the exchange potential for state a
  Vex_a(r)  = -\sum_b [Jb]*vex_ab(r)*fac   [fac = (fafb+gagb))/(fafa+gaga)]
  vex_ab(r) =  \sum_k Lam_ab^k v_ab^k(r)
  v_ab^k(r) =  Int_0^inf  [rmin^k/rmax^(k+1)]*(fafb+gagb)(r') dr'
               [rmin= min(r,r')]
  Uses same efficient method as for Vdir to calc this integral
Note sign convention: V = V_dir + V_ex [-ve inside Vex!]
Note: the "fac" method is not 100% correct. Unstable for small f_a
I get around this by introducing a cut-off. OK, but not ideal.
Parallelised over core states. note: Must use omp critical
NOTE: Doesn't work for core f states! Because of above dodgy??
For now: exclude core f states. Ok, but not ideal
*/
{
  int ngp = wf.ngp;
  int tja = ATI::twoj_k(wf.kappa[a]);
  int la = ATI::l_k(wf.kappa[a]);

  vex_a.clear(); //must be clear!
  vex_a.resize(ngp);
  #pragma omp parallel for
  for(int b=0; b<wf.num_core_states; b++){

    //XXX HERE! Ignores f states in core!!! XXX
    if(wf.kappa[b]==3||wf.kappa[b]==-4) continue;
    //Arrays and values needed:
    std::vector<double> Vxab(ngp);
    int irmax = std::min(wf.pinflist[a],wf.pinflist[b]);
    std::vector<double> L_abk; //Lambda_ab^k (angular factor)
    int tjb = ATI::twoj_k(wf.kappa[b]);
    int lb = ATI::l_k(wf.kappa[b]);
    int k_min = formLambdaABk(L_abk,tja,tjb,la,lb); //fills Lam array
    int k_max = (tja + tjb)/2;
    double stf = wf.core_ocf[b]*(tjb+1); //avg over non-rel configs
    //Form Vex_ab^k, and sum over k
    for(int k = k_min; k<=k_max; k++){
      double Labk = L_abk[k];
      if(Labk==0) continue;
      //Calculate Vex_ab^k
      std::vector<double> A(ngp),B(ngp),Vxabk(ngp);
      A[0] = 0;
      double b0=0;
      //Use INT formulas??
      for(int i=0; i<irmax; i++) b0 += wf.drdt[i]*
        (wf.p[a][i]*wf.p[b][i]+wf.q[a][i]*wf.q[b][i])/pow(wf.r[i],k+1);
      B[0] = b0;
      for(int i=1; i<ngp; i++){
        double Fdr = wf.drdt[i-1]*
          (wf.p[a][i-1]*wf.p[b][i-1]+wf.q[a][i-1]*wf.q[b][i-1]);
        A[i] = A[i-1] + Fdr*pow(wf.r[i-1],k);
        B[i] = B[i-1] - Fdr/pow(wf.r[i-1],k+1);
        Vxabk[i] = wf.h*(A[i]/pow(wf.r[i],k+1) + B[i]*pow(wf.r[i],k));
      }
      for(int i=1; i<ngp; i++){
        Vxab[i] += Vxabk[i]*Labk;
      }
    }//k loop
    //Finally, sum each b into Vex:
    for(int ir=0; ir<ngp; ir++){
      double fac = 1;
      if(a!=b){
        //NB: Cutt-off? OK?
        if(fabs(wf.p[a][ir])<1.e-3) continue;
        double fac_top = wf.p[a][ir]*wf.p[b][ir] + wf.q[a][ir]*wf.q[b][ir];
        double fac_bot = wf.p[a][ir]*wf.p[a][ir] + wf.q[a][ir]*wf.q[a][ir];
        fac = fac_top/fac_bot;
      }
      #pragma omp critical
      {
        vex_a[ir] += -1*stf*Vxab[ir]*fac;
      }
    }//r
  }//b loop
  return 0;
}

//******************************************************************************
int formLambdaABk(std::vector<double> &L_abk, int tja, int tjb, int la, int lb)
/*
Lambda angular factor for exchange potential
L_abk = 3js(ja,jb,k,-1/2,1/2,2)^2 * parity(la+lb+k)
*/
{
  int k_min = fabs(tja - tjb)/2;
  int k_max = (tja + tjb)/2;
  L_abk.clear();
  L_abk.resize(k_max+1,0);
  #pragma omp parallel for
  for(int k = k_min; k<=k_max; k++){
    double tL = pow(WIG::threej_2(tja,tjb,2*k,-1,1,0),2);
    tL *= WIG::parity(la,lb,k);
    L_abk[k] = tL;
  }
  return k_min;
}


//******************************************************************************
int hartreeCore(ElectronOrbitals &wf, double eps_hartree)
/*
Solves the Hartree equations (no exchange)
NOTE: can make this nicer.. but will never use it, so whatever.
*/
{
  // "Mixing" of new + old Potential:
  // V_n+1 = eta*V_n+1 + (1-eta)*V_n
  const double eta=0.50;

  //Fill the electron part of the potential, using Greens PRM for initial approx
  double Gh,Gd;  //Green potential parameters
  PRM::defaultGreen(wf.Z,Gh,Gd); //Get default values for Green potential
  for(int i=0; i<wf.ngp; i++) wf.vdir.push_back(PRM::green(wf.Z,wf.r[i],Gh,Gd));

  //First step: Solve each core state using parameteric potential
  wf.solveInitialCore(1); //1, since don't need high accuray here [1 in 10^1]
  int Ncs = wf.num_core_states;

  //Hartree loop:
  int num_its=0;
  for(int n=0; n<MAX_HART_ITS; n++){

    //Use known orbitals to form new potential:
    std::vector<double> vdir_old = wf.vdir;
    //std::vector<double> vdir_new;
    formNewVdir(wf,wf.vdir,true);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*wf.vdir[j] + (1.-eta)*vdir_old[j];
    }

    //Solve dirac equation for each (Core) orbital in new potential
    double prev_e = 0;
    for(int i=0; i<Ncs; i++) prev_e += wf.en[i]/Ncs;
    for(int i=0; i<Ncs; i++){
      double del_e=0;
      for(int j=0; j<wf.pinflist[i]; j++)
        del_e += (wf.vdir[j]-vdir_old[j])*
        (pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
      del_e*=wf.h;
      double new_e = wf.en[i] + 1*del_e;
      if(new_e>0)new_e=-0.1;
      wf.reSolveDirac(i,new_e,3); //only go to 1/10^3 - do better at end!
    }
    double next_e = 0;
    for(int i=0; i<Ncs; i++) next_e += wf.en[i]/Ncs;

    //check for convergence:
    //NB: eta in denom, otherwise v. small eta will spuriously give small delta
    double delta_hartree = (next_e-prev_e)/(next_e*eta);
    num_its = n+1;
    //Output Hartree it number/eps to screen (overwrite line):
    printf("Hartree it:%3i,  del=%6.0e     \r",num_its,delta_hartree);
    std::cout<<std::flush; //flush cout to update screen output
    if(fabs(delta_hartree)<eps_hartree) break;
  }
  std::cout<<"\n";

  //re-run solve Dirac to higher convergance level after Hart pot. ok
  for(int i=0; i<Ncs; i++) wf.reSolveDirac(i,0,14);
  //Form the total core potential using new wfs
  //This time, solved for case of valence states (different factor)
  formNewVdir(wf,wf.vdir,false);


  // //Now, solve using Hartree
  // //Move this to seperate routine! XXX
  // int hits;
  // double eps_hart = 1.e-3;
  // for(hits=1; hits<MAX_HART_ITS; hits++){
  //   //Use known orbitals to form new potential:
  //   std::vector<double> vdir_old = wf.vdir;
  //   std::vector<double> vdir_new;
  //   formNewVdir(wf,vdir_new,true); //true: means scale!
  //   for(int j=0; j<wf.ngp; j++)
  //     wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
  //
  //   std::vector<double> en_old = wf.en;
  //   double t_eps = 0;
  //   for(int i=0; i<Ncs; i++){
  //     //Use PT to find new energy guess
  //     double del_e=0;
  //     for(int j=0; j<wf.ngp; j++) del_e += (wf.vdir[j]-vdir_old[j])
  //       *(pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
  //     del_e*=wf.h;
  //     double en_guess = en_old[i] + del_e;
  //     if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
  //     wf.reSolveDirac(i,en_guess,3); //only go to 1/10^3 here
  //     double sfac = 2.*wf.kappa[i]*wf.core_ocf[i]; //|2k|=2j+1
  //     t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
  //   }
  //   t_eps /= (wf.num_core_electrons*eta);
  //   //Note: just weighted average of eps for each orbital.. good 'nuff
  //   std::cout<<"Hartree it: "<<hits<<" "<<t_eps<<" "<<wf.en[0]<<"\n";
  //   if(t_eps<eps_hart && hits>1) break;
  // }

  return num_its;
}

}//Namespace
