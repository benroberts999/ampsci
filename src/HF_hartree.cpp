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

  //Hartree loop:
  int num_its=0;
  for(int n=0; n<max_hartree; n++){

    //Use known orbitals to form new potential:
    std::vector<double> vdir_old = wf.vdir;
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
    }

    //Solve dirac equation for each (Core) orbital in new potential
    double prev_e = 0;
    for(int i=0; i<wf.num_core; i++) prev_e += wf.en[i]/wf.num_core;
    for(int i=0; i<wf.num_core; i++){
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
    for(int i=0; i<wf.num_core; i++) next_e += wf.en[i]/wf.num_core;

    //check for convergence:
    //NB: eta in denom, otherwise v. small eta will spuriously give small delta
    double delta_hartree = (next_e-prev_e)/(next_e*eta);
    num_its = n+1;
    //printf("Hart it:%3i,  del=%6.0e\n",num_its,delta_hartree);
    printf("Hart it:%3i,  del=%6.0e     \r",num_its,delta_hartree);
    if(fabs(delta_hartree)<eps_hartree) break;
  }
  std::cout<<"\n";

  //re-run solve Dirac to higher convergance level after Hart pot. ok
  for(int i=0; i<wf.num_core; i++) wf.reSolveLocalDirac(i,0,14);
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
  for(int i=0; i<wf.num_core; i++){
    int ka = wf.klist[i];
    int twoj = ATI::twoj_k(ka);
    double frac = wf.core_ocf[i];
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





}//Namespace
