#include "HF_hartree.h"

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
    formNewVdir(wf,vdir_new);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
    }

    //Solve dirac equation for each (Core) orbital in new potential
    double prev_e = 0;
    for(int i=0; i<Ncs; i++) prev_e += wf.en[i]/Ncs;
    for(int i=0; i<Ncs; i++){
      double del_e=0;
      for(int j=0; j<wf.ngp; j++)
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

  // std::ofstream ofile;
  // ofile.open("pot.txt");
  // ofile<<"r Gr Vh Z/r 1/r\n";
  // for(int i=0; i<wf.ngp; i++){
  //   ofile<<wf.r[i]<<" "
  //     <<-PRM::green(wf.Z,wf.r[i],Gh,Gd)-wf.vnuc[i]<<" "
  //     <<-wf.vdir[i]-wf.vnuc[i]<<" "
  //     <<wf.Z/wf.r[i]<<" "<<1./wf.r[i]<<"\n";
  //
  // }
  // ofile.close();

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
  const double eta=0.50;
  int ngp = wf.ngp;

  wf.vdir.clear(); //make sure it's empty

  //Fill the electron part of the potential, using Greens PRM for initial approx
  double Gh,Gd;  //Green potential parameters
  PRM::defaultGreen(wf.Z,Gh,Gd); //Get default values for Green potential
  for(int i=0; i<ngp; i++) wf.vdir.push_back(PRM::green(wf.Z,wf.r[i],Gh,Gd));

  //First step: Solve each core state using parameteric potential
  wf.solveInitialCore(3); //3, since don't need high accuray here [1 in 10^3]
  int Ncs = wf.num_core_states;

  //Now, solve using Hartree
  //Move this to seperate routine! XXX
  int hits;
  double eps_hart = 1.e-3;
  for(hits=1; hits<MAX_HART_ITS; hits++){
    //Use known orbitals to form new potential:
    std::vector<double> vdir_old = wf.vdir;
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new,true); //true: means scale!
    for(int j=0; j<wf.ngp; j++)
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];

    std::vector<double> en_old = wf.en;
    double t_eps = 0;
    for(int i=0; i<Ncs; i++){
      //Use PT to find new energy guess
      double del_e=0;
      for(int j=0; j<wf.ngp; j++) del_e += (wf.vdir[j]-vdir_old[j])
        *(pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
      del_e*=wf.h;
      double en_guess = en_old[i] + del_e;
      if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
      wf.reSolveLocalDirac(i,en_guess,3); //only go to 1/10^3 here
      double sfac = 2.*wf.kappa[i]*wf.core_ocf[i]; //|2k|=2j+1
      t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
    }
    t_eps /= (wf.num_core_electrons*eta);
    //Note: just weighted average of eps for each orbital.. good 'nuff
    std::cout<<"Hartree it: "<<hits<<" "<<t_eps<<"\n";
    if(t_eps<eps_hart && hits>1) break;
  }


  std::vector< std::vector<double> > vex; //into class?? XXX
  vex.clear(); //move to beginning ?
  vex.resize(Ncs,std::vector<double>(ngp));

  for(hits=1; hits<MAX_HART_ITS; hits++){
    std::vector<double> en_old = wf.en;

    //Form v_dir & v_ex
    std::vector<double> vdir_old = wf.vdir;
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new,true);

    std::vector< std::vector<double> > vex_old = vex; //XX class!
    formVexCore(wf,vex);

    //calculate de from PT
    for(int i=0; i<Ncs; i++){
      double del_e = 0;
      for(int j=0; j<wf.ngp; j++){
        double dv = wf.vdir[j]+vex[i][j]-vdir_old[j]-vex_old[i][j];
        del_e += dv*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
      }
      del_e*=wf.h;
      double en_guess = en_old[i] + del_e;
      if(en_guess>0) en_guess = en_old[i]; //safety, should never happen
      wf.reSolveLocalDirac(i,en_guess,3); //only go to 1/10^3 here
      double sfac = 2.*wf.kappa[i]*wf.core_ocf[i]; //|2k|=2j+1
      t_eps += fabs(sfac*(wf.en[i]-en_old[i])/en_old[i]);
    }
    t_eps /= (wf.num_core_electrons*eta);
    //Note: just weighted average of eps for each orbital.. good 'nuff
    std::cout<<"Hartree it: "<<hits<<" "<<t_eps<<"\n";
    if(t_eps<eps_HF && hits>1) break;
  }

}

//******************************************************************************
int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, bool core)
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

  //Count number of electrons in the core
  int Ncore=0;
  for(size_t i=0; i<wf.core_list.size(); i++) Ncore+=wf.core_list[i];

  //Factor: When solveing for core N=Ncore. For valence, N=Ncore+_
  double f=1;
  if(core) f = 1. - (1.)/Ncore;

  //Determine the total electron (charge) density of core:
  // 2j+1 is closed-shell occupation number of that orbital.
  std::vector<double> rho(wf.ngp);
  for(int i=0; i<wf.num_core_states; i++){
    int ka = wf.kappa[i];
    int twoj = ATI::twoj_k(ka);
    double frac = wf.core_ocf[i]; //avgs over non-rel. configs
    for(int j=0; j<wf.ngp; j++){
      rho[j] += frac*(twoj+1)*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2));
    }
  }

  //This is the slow part:
  //Ise the above determined electron (charge) density with Gauss' law
  //to determine the electric potential (energy). Makes use of spherical symm.
  #pragma omp parallel for
  for(int j=0; j<wf.ngp; j++){
    double r = wf.r[j];
    double v_tmp = 0;
    for(int k=0; k<wf.ngp; k++){
      double rp = wf.r[k];
      double rm = std::max(r,rp);//can be little more clever, slight speedup
      v_tmp += (rho[k]/rm)*wf.drdt[k];
    }
    vdir_new[j] = f*v_tmp*wf.h;
  }

  return 0;
}

//==============================================================================

//******************************************************************************
int formVexCore(ElectronOrbitals &wf, std::vector< std::vector<double> > &vex){

  vex.clear();
  vex.resize(wf.num_core_states);
  #pragma omp parallel for
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
  for(int b=0; b<wf.num_core_states; b++){// Core? or all??
    int tjb = ATI::twoj_k(wf.kappa[b]);
    int lb = ATI::l_k(wf.kappa[b]);
    std::vector<double> L_abk;
    int k_min = formLambdaABk(L_abk,tja,tjb,la,lb);
    double stf = wf.core_ocf[b]*(tjb+1); //avg over non-rel configs
    for(int ir=0; ir<ngp; ir++){
      double fac_bot = wf.p[a][ir]*wf.p[a][ir] + wf.q[a][ir]*wf.q[a][ir]; //?
      double fac_top = wf.p[a][ir]*wf.p[b][ir] + wf.q[a][ir]*wf.q[b][ir];
      double vex_ab_r = vexABr(wf,a,b,ir,L_abk,k_min);
      vex_a[ir] += stf*vex_ab_r*fac_top/fac_bot;
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
  L_abk.resize(k_max);
  for(int k = k_min; k<=k_max; k++){
    double tL = pow(WIG::threej_2(tja,tjb,2*k,-1,1,0),2);
    tL *= WIG::parity(la,lb,k);
    L_abk[k] = tL;
  }
  return k_min;
}

//******************************************************************************
double vexABr(ElectronOrbitals &wf, int a, int b, int ir,
  std::vector<double> &L_abk, int k_min)
/*
Not sure if this is the best way to do this..
vex_ab(r) = sum_k L_abk * vex_abk(r)   [2jb+1 already above]
v_abk(r) = int [r_min^k/r_max^(k+1)]*[fafb + gagb](r2) dr_2
r_min = min(r,r2) //nb: can test the indices!
*/
{
  int k_max = (int) L_abk.size();

  double vex_ab_r = 0;
  for(int jr=0; jr<wf.ngp; jr++){
    double xL = 0;
    for(int k = k_min; k<=k_max; k++){
      double xLk = L_abk[k];
      if(xLk==0) continue;
      if(ir>jr) xLk *= pow(wf.r[jr]/wf.r[ir],k)/wf.r[ir];
      else      xLk *= pow(wf.r[ir]/wf.r[jr],k)/wf.r[jr];
      xL += xLk;
    }
    double FF = (wf.p[a][jr]*wf.p[b][jr] + wf.q[a][jr]*wf.q[b][jr]);
    vex_ab_r += xL*FF*wf.drdt[jr];
  }
  return vex_ab_r*wf.h;
}

// //******************************************************************************
// int vkScreeningOLD(ElectronOrbitals &wf, int a, int b, std::vector<double> vk)
// /*
// Pass by reference? Faster? Worse??
// Averaged exchange potential! Break into functions??
// NOTE: must be multiplied
// */
//
// {
//
//   vk.resize(wf.ngp);
//
//   int ka = wf.kappa[a];
//   int kb = wf.kappa[b];
//   int tja = ATI::twoj_k(ka);
//   int tjb = ATI::twoj_k(kb);
//   int la = ATI::l_k(ka);
//   int lb = ATI::l_k(kb);
//
//   int ja_m_jb = fabs(tja - tjb)/2;
//   int ja_p_jb = (tja + tjb)/2;
//
//   //Stats factors (averaged over non-rel configs)
//   double sfa = (tja+1)*wf.core_ocf[a];
//   double sfb = (tjb+1)*wf.core_ocf[b];
//
//   //Count number of electrons in the core
//   int Ncore=0;
//   for(size_t i=0; i<wf.core_list.size(); i++) Ncore+=wf.core_list[i];
//
//   std::vector<double> Lam_abk;
//   std::vector<int> k_list;
//   for(int k = ja_m_jb; k<=ja_p_jb; k++){
//     double tL = pow(WIG::threej_2(tja,tjb,2*k,-1,1,0),2);
//     tL *= parity(la,lb,k);
//     if(tL != 0){
//       k_list.push_back(k);
//       Lam_abk.push_back(tL);
//     }
//   }
//   if(Lam_abk.size()==0) return 0;
//
//   #pragma omp parallel for
//   for(int i=0; i<wf.ngp; i++){
//     double r = wf.r[i];
//
//     double v_tmp = 0;
//     for(int j=0; j<wf.ngp; j++){
//       double rp = wf.r[j];
//       double t_X = (wf.p[a][j]*wf.p[b][j] + wf.q[a][j]*wf.q[b][j])*wf.drdt[i];
//       for(size_t ik=0; ik<k_list.size(); ik++){
//         int k = k_list[ik];
//         double tL = Lam_abk[ik];
//         double x_k;
//         if(rp<r) x_k = pow(rp/r,k)/r;
//         else     x_k = pow(r/rp,k)/rp;
//         v_tmp += tL*x_k*t_X;
//       }
//     }
//     double sf = sfa*sfb*(1./Ncore)
//     v_tmp *= sf*wf.h*(wf.p[a][i]*wf.p[b][i] + wf.q[a][i]*wf.q[b][i]);
//
//   }
// }


}//Namespace
