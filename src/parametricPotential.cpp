#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "ATI_atomInfo.h"
#include "FPC_physicalConstants.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "ContinuumOrbitals.h"

int main(void){

  clock_t ti,tf;
  ti = clock();

  double varalpha=1; //need same number as used for the fitting!

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;

  double Tf,Tt,Tg;  //Teitz potential parameters
  double Gf,Gh,Gd;  //Green potential parameters

  int n_max,l_max;
  std::string str_core;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("parametricPotential.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    ifs >> str_core;              getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> Gf >> Gh >> Gd;        getline(ifs,jnk);
    ifs >> Tf >> Tt >> Tg;        getline(ifs,jnk);
    ifs >> n_max >> l_max;        getline(ifs,jnk);
    ifs.close();
  }

  int Z = ATI::get_z(Z_str);
  if(Z==0) return 2;
  if(A==0) A=ATI::A[Z]; //if none given, get default A

  //Normalise the Teitz/Green weights:
  if(Gf!=0 || Tf!=0){
    double TG_norm = Gf + Tf;
    Gf /= TG_norm;
    Tf /= TG_norm;
  }

  //If H,d etc are zero, use default values
  if(Gf!=0 && Gh==0) PRM::defaultGreen(Z,Gh,Gd);
  if(Tf!=0 && Tt==0) PRM::defaultTietz(Z,Tt,Tg);


  printf("\nRunning parametric potential for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");
  if(Gf!=0) printf("%3.0f%% Green potential: H=%.4f  d=%.4f\n",Gf*100.,Gh,Gd);
  if(Tf!=0) printf("%3.0f%% Tietz potential: T=%.4f  g=%.4f\n",Tf*100.,Tt,Tg);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  //if(A!=0) wf.sphericalNucleus();

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[wf.ngp-1]);

  //Determine which states are in the core:
  int core_ok = wf.determineCore(str_core);
  if(core_ok==2){
    std::cout<<"Problem with core: "<<str_core<<"\n";
    return 1;
  }

  //Fill the electron part of the potential
  wf.vdir.resize(wf.ngp);
  for(int i=0; i<wf.ngp; i++){
    double tmp = 0;
    if(Gf!=0) tmp += Gf*PRM::green(Z,wf.r[i],Gh,Gd);
    if(Tf!=0) tmp += Tf*PRM::tietz(Z,wf.r[i],Tt,Tg);
    wf.vdir[i] = tmp;
  }

  //Solve for core states
  wf.solveInitialCore();

  //store number of calculated core states:
  int num_core_states = wf.nlist.size();

  //Calculate the valence (and excited) states
  for(int n=1; n<=n_max; n++){
    for(int l=0; l<=l_max; l++){
      if(l+1>n) continue;

      for(int tk=0; tk<2; tk++){
        int k;
        if(tk==0) k=l;
        else      k=-(l+1);
        if(k==0) continue;
        if(wf.isInCore(n,k)) continue;

        double en_a = wf.enGuessVal(n,k);
        wf.solveLocalDirac(n,k,en_a);

      }
    }
  }

  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\n n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    if((int)m==num_core_states){
      std::cout<<" ========= Valence: ======\n";
      printf(" n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
    }
    int n=wf.nlist[i];
    int k=wf.kappa[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %15.3f\n",
        n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*FPC::Hartree_invcm);
  }

  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);


  return 0;
}
