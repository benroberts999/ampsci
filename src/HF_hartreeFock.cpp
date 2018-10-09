#include "HF_hartreeFock.h"

/*
8 Oct 2018

Works, but not perfectly.
Orthonormality only reached to ~1e4 - strange. ~1e10 for Hartree.

 * Use INT function?
 * Do iterations for each core state sepperately! Might be faster??
 * "Cheat" thing with fac_bot & fac_top ???

 * Needs some more clean-up!

*/



namespace HF{

//******************************************************************************
int hartreeCore(ElectronOrbitals &wf, double eps_hartree)
/*
Solves the Hartree equations (no exchange term yet)
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
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new,true);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
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
      wf.reSolveLocalDirac(i,new_e,3); //only go to 1/10^3 - do better at end!
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
  for(int i=0; i<Ncs; i++) wf.reSolveLocalDirac(i,0,14);
  //Form the total core potential using new wfs
  //This time, solved for case of valence states (different factor)
  formNewVdir(wf,wf.vdir,false);

  return num_its;
}

//******************************************************************************
int hartreeFockCore(ElectronOrbitals &wf, double eps_HF)
{
/*
  1) Green potential for starting approx.
  2) Hartree for the core (vdir)
  3) Add exchange. Loop over all states? Or each state?
*/

  // "Mixing" of new + old Potential:
  // V_n+1 = eta*V_n+1 + (1-eta)*V_n
  //NB: must be less than 0.5!, because x2 after first few its!
  double eta=0.35;
  double eta_x=0.35;
  int ngp = wf.ngp;

  wf.vdir.clear(); //make sure it's empty

  //Fill the electron part of the potential, using Greens PRM for initial approx
  double Gh,Gd;  //Green potential parameters
  PRM::defaultGreen(wf.Z,Gh,Gd); //Get default values for Green potential
  for(int i=0; i<ngp; i++) wf.vdir.push_back(PRM::green(wf.Z,wf.r[i],Gh,Gd));

  //First step: Solve each core state using parameteric potential
  wf.solveInitialCore(3); //3, since don't need high accuray here [1 in 10^3]
  int Ncs = wf.num_core_states;

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
  //     wf.reSolveLocalDirac(i,en_guess,3); //only go to 1/10^3 here
  //     double sfac = 2.*wf.kappa[i]*wf.core_ocf[i]; //|2k|=2j+1
  //     t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
  //   }
  //   t_eps /= (wf.num_core_electrons*eta);
  //   //Note: just weighted average of eps for each orbital.. good 'nuff
  //   std::cout<<"Hartree it: "<<hits<<" "<<t_eps<<" "<<wf.en[0]<<"\n";
  //   if(t_eps<eps_hart && hits>1) break;
  // }

  std::vector< std::vector<double> > vex; //into class?? XXX
  vex.clear(); //move to beginning ?
  vex.resize(Ncs,std::vector<double>(ngp));

  //XXX One reason it's slow: Iterates V_ex for each orbital, even when some ok!
  // SO: do for each orbital sepperately!

  int hits;
  for(hits=1; hits<MAX_HART_ITS; hits++){

    if(hits==4){
      eta *= 2;
      eta_x *= 2;
    }

    //Form v_dir & v_ex
    std::vector<double> vdir_old = wf.vdir;
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new,false);

    std::vector< std::vector<double> > vex_old = vex;
    formVexCore(wf,vex);

    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
      for(int i=0; i<Ncs; i++){
        vex[i][j] = eta_x*vex[i][j] + (1.-eta_x)*vex_old[i][j];
      }
    }

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
      wf.reSolveLocalDirac(i,en_guess,vex[i],2); //only go to 1/10^3 here
      double sfac = 2.*wf.kappa[i]*wf.core_ocf[i]; //|2k|=2j+1
      t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
    }
    t_eps /= (wf.num_core_electrons*eta);
    //Note: just weighted average of eps for each orbital.. good 'nuff

    printf("\rHF core it: %3i, eps= %6.1e              ",hits,t_eps);
    std::cout<<std::flush;
    if(t_eps<eps_HF && hits>2) break;
  }
  std::cout<<"\n";

  for(int i=0; i<Ncs; i++) wf.reSolveLocalDirac(i,wf.en[i],vex[i],15);
  formNewVdir(wf,wf.vdir,false);

  return 0;
}

//******************************************************************************
int hartreeFockValence(ElectronOrbitals &wf, int na, int ka, double eps_HF)
/*
Calculate valence states in frozen core
*/
{

  if(wf.isInCore(na,ka)){
    std::cout<<"\nHF Warning199: "<<na<<" "<<ka<<" is in the core!\n";
    return 0;
  }

  double eta=0.35;

  double en_g = wf.enGuessVal(na,ka);
  wf.solveLocalDirac(na,ka,en_g,3);
  int a = wf.nlist.size() - 1;

  std::vector<double> vexa(wf.ngp);
  int hits;
  for(hits=1; hits<MAX_HART_ITS; hits++){
    if(hits==4) eta*=2;
    std::vector<double> vexa_old = vexa;
    double en_old = wf.en[a];
    std::vector<double> vexa_new;
    formVexA(wf,a,vexa_new);

    for(int i=0; i<wf.ngp; i++) vexa[i] = eta*vexa_new[i]+(1.-eta)*vexa_old[i];
    wf.reSolveLocalDirac(a,wf.en[a],vexa,2);
    double eps = fabs((wf.en[a]-en_old)/(eta*en_old));
    printf("\rHF val:%2i %2i %2i | %3i eps=%6.1e  en=%11.8f                  "
    ,a,na,ka,hits,eps,wf.en[a]);
    std::cout<<std::flush;
    if(eps<eps_HF) break;
  }
  std::cout<<"\n";
  wf.reSolveLocalDirac(a,wf.en[a],vexa,15);
  return hits;
}

//******************************************************************************
int formNewVdir(ElectronOrbitals &wf, std::vector<double> &vdir_new, bool scale)
/*
This takes in the wavefunctions, and forms the direct (local) part of the
electron potential. Does not include exchange.
NOTE: For now, assumes core potential is the same
for each orbital, which is a good approximation.
When core=true, will solve for all core states with (1-1/N) (V^N-1 for core)
When core=false, will solve for core. Good for valence states [N is different!]
core=true by default
*/
{

  //Make sure vector is correct size, and clear old potential away
  vdir_new.clear();
  vdir_new.resize(wf.ngp);

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


/*
Use the above determined electron (charge) density with Gauss' law
to determine the electric potential (energy). Makes use of spherical symm.
Note: Uses efficient integral method:

  vdir(r) = Int_0^inf rho(r')/r_max dr'
          = Int_0^r rho(r')/r dr' + Int_r^inf rho(r')/r' dr'
          = A(r)/r + B(r)
  A(r0)   = 0
  B(r0)   = Int_0^inf rho(r')/r' dr'
  A(r_n)  = A(r_n-1) + rho(r)dr
  B(r_n)  = B(r_n-1) - (rho(r)/r)dr
  vdir(r) = A(r)/r + B(r)

Note: I don't use quadrature formula here.
*/

  std::vector<double> A(wf.ngp),B(wf.ngp);

  //Note: I don't use quadrature formula here.
  double b0 = 0;
  for(int i=0; i<wf.ngp; i++) b0 += wf.drdt[i]*rho[i]/wf.r[i];

  B[0] = b0; //INT::integrate2(rho_on_r,wf.drdt,1.,0,wf.ngp,0,0);
  A[0] = 0.;
  vdir_new[0] = B[0];
  for(int i=1; i<wf.ngp; i++){
    B[i] = B[i-1] - wf.drdt[i-1]*rho[i-1]/wf.r[i-1];
    A[i] = A[i-1] + wf.drdt[i-1]*rho[i-1];
    vdir_new[i] = f*wf.h*(A[i]/wf.r[i] + B[i]);
  }

  return 0;
}

//==============================================================================

//******************************************************************************
int formVexCore(ElectronOrbitals &wf, std::vector< std::vector<double> > &vex){

  vex.clear();
  vex.resize(wf.num_core_states);
  //#pragma omp parallel for // - no!, actually better to //ise below
  for(int a=0; a<wf.num_core_states; a++){
    std::vector<double> vex_a;
    formVexA(wf,a,vex_a);
    vex[a] = vex_a; //double-check this works w/ vectors..
  }
  return 0;
}

//******************************************************************************
int formVexA(ElectronOrbitals &wf, int a, std::vector<double> &vex_a)
/*
XXX pass by reference? No?
for now, include a in the b summation too.
Means we should not scale v_dir?
v_tot(r) should still go to (N-M)/r for large r - check!
*/
{
  int ngp = wf.ngp;

  int tja = ATI::twoj_k(wf.kappa[a]);
  int la = ATI::l_k(wf.kappa[a]);

  vex_a.clear();
  vex_a.resize(ngp);
  for(int b=0; b<wf.num_core_states; b++){
    int tjb = ATI::twoj_k(wf.kappa[b]);
    int lb = ATI::l_k(wf.kappa[b]);
    std::vector<double> L_abk;
    int k_min = formLambdaABk(L_abk,tja,tjb,la,lb);
    double stf = wf.core_ocf[b]*(tjb+1); //avg over non-rel configs
    #pragma omp parallel for
    for(int ir=0; ir<ngp; ir++){
      double fac = 1;
      if(a!=b){
        //NB: Creates bad instability when p[a] is small!
        //--so I added cut-offs. OK? Works!(?)
        double fac_top = 1.e-6 +
          (wf.p[a][ir]*wf.p[b][ir] + wf.q[a][ir]*wf.q[b][ir]);
        double fac_bot = 1.e-6 +
          (wf.p[a][ir]*wf.p[a][ir] + wf.q[a][ir]*wf.q[a][ir]);
        if(fac_bot!=0 && fac_top!=0) fac = fac_top/fac_bot;
        else fac=0;
      }
      if(fac==0) continue; //continue, not break. a)OMP, b) r=0 case.
      double vex_ab_r = vexABr(wf,a,b,ir,L_abk,k_min);
      vex_a[ir] += -1*stf*vex_ab_r*fac; //nb: minus sign
    }
  }
  return 0;
}

//******************************************************************************
int formLambdaABk(std::vector<double> &L_abk, int tja, int tjb, int la, int lb)
/*
L_abk = 3js(ja,jb,k,-1/2,1/2,2)^2 * parity(la+lb+k)
*/
{
  int k_min = fabs(tja - tjb)/2;
  int k_max = (tja + tjb)/2;
  L_abk.clear();
  L_abk.resize(k_max+1,0);
  for(int k = k_min; k<=k_max; k++){
    double tL = pow(WIG::threej_2(tja,tjb,2*k,-1,1,0),2);
    tL *= WIG::parity(la,lb,k);
    L_abk[k] = tL;
  }
  return k_min;
}

//******************************************************************************
double vexABr(ElectronOrbitals &wf, int a, int b, int ir,
  std::vector<double> L_abk, int k_min)
/*
Not sure if this is the best way to do this..
vex_ab(r) = sum_k L_abk * vex_abk(r)   [2jb+1 already above]
v_abk(r) = int [r_min^k/r_max^(k+1)]*[fafb + gagb](r2) dr_2
r_min = min(r,r2) //nb: can test the indices!
*/
{
  int k_max = (int) L_abk.size() - 1;

  double vex_ab_r = 0;
  int jr_max = std::min(wf.pinflist[a],wf.pinflist[b]);

  for(int k = k_min; k<=k_max; k++){
    double Lk = L_abk[k];
    if(Lk==0) continue;
    for(int jr=0; jr<jr_max; jr++){
      double xk = 0;
      if(ir>jr) xk = pow(wf.r[jr]/wf.r[ir],k)/wf.r[ir];
      else      xk = pow(wf.r[ir]/wf.r[jr],k)/wf.r[jr];
      double FF = (wf.p[a][jr]*wf.p[b][jr] + wf.q[a][jr]*wf.q[b][jr]);
      vex_ab_r += xk*Lk*FF*wf.drdt[jr];
    }
  }

  return vex_ab_r*wf.h;
}


}//Namespace
