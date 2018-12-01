#include "ElectronOrbitals.h"
#include "HF_hartreeFock.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "FPC_physicalConstants.h"
#include "ATI_atomInfo.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "ChronoTimer.h"

int main(void){

  ChronoTimer sw; //start the stopwatch
  sw.start();

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;
  double varalpha,varalpha2;
  double eps_hart;
  int iHF;

  int n_max,l_max;
  std::string str_core;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("hartreeFock.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    ifs >> str_core;              getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> iHF >> eps_hart;       getline(ifs,jnk);
    ifs >> n_max >> l_max;        getline(ifs,jnk);
    ifs >> varalpha2;             getline(ifs,jnk);
    ifs.close();
  }
  bool doHF = true;
  if(iHF==0) doHF = false;

  //Change varAlph^2 to varalph
  if(varalpha2==0) varalpha2 = 1.e-10;
  varalpha = sqrt(varalpha2);

  int Z = ATI::get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI::A[Z]; //if none given, get default A

  if(!doHF) printf("\nRunning HARTREE for %s, Z=%i A=%i\n",Z_str.c_str(),Z,A);
  else printf("\nRunning HARTREE FOCK for %s, Z=%i A=%i\n",Z_str.c_str(),Z,A);
  printf("*************************************************\n");

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  if(A>0) wf.sphericalNucleus();

  printf("Grid: pts=%i h=%7.5f r0=%.1e Rmax=%5.1f\n\n"
  ,  wf.ngp,wf.h,wf.r[0],wf.r[wf.ngp-1]);

  //Determine which states are in the core:
  int core_ok = wf.determineCore(str_core);
  if(core_ok==2){
    std::cout<<"Problem with core: "<<str_core<<"\n";
    return 1;
  }

  sw.start();
  //Solve Hartree equations for the core:
  if(doHF) HF::hartreeFockCore(wf,eps_hart);
  else     HF::hartreeCore(wf,eps_hart);
  std::cout<<"core: "<<sw.lap_reading_str()<<"\n";

  //Make a list of valence states to solve:
  std::vector< std::vector<int> > lst;
  for(int n=0; n<=n_max; n++){
    for(int l=0; l<n; l++){
      if(l>l_max) continue;
      for(int tk=0; tk<2; tk++){ //loop over k
        int k;
        if(tk==0) k=l;      //j = l - 1/2
        else      k=-(l+1); //j = l + 1/2
        if(k==0) continue;  // no j = l - 1/2 for l=0
        if(wf.isInCore(n,k)) continue;
        lst.push_back({n,k});
      }
    }
  }
  //Solve for the valence states:
  sw.start();
  for(uint i=0; i<lst.size(); i++){
    int n = lst[i][0];
    int k = lst[i][1];
    if(doHF) HF::hartreeFockValence(wf,n,k,eps_hart);
    else{
      //double en_g = wf.enGuessVal(n,k);
      wf.solveLocalDirac(n,k,0,15);
    }
  }
  std::cout<<"Valence: "<<sw.lap_reading_str()<<"\n";

  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\nCore: %s, Z=%i A=%i\n",Z_str.c_str(),Z,A);
  printf(" n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
  bool val=false; double en_lim=0;
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    if((int)m==wf.num_core_states){
      en_lim = fabs(wf.en[i]);
      val = true;
      std::cout<<" ========= Valence: ======\n";
      printf(
      " n l_j    k Rinf its    eps      En (au)        En (/cm)   En (/cm)\n");
    }
    int n=wf.nlist[i];
    int k=wf.kappa[i];
    int twoj = ATI::twoj_k(k);
    int l = ATI::l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %15.3f",
        n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*FPC::Hartree_invcm);
    if(val)printf(" %10.2f\n",(eni+en_lim)*FPC::Hartree_invcm);
    //else std::cout<<" /eng="<<wf.enGuessCore(n,l)<<"\n";
    else std::cout<<"\n";
  }

  std::cout<<"\n Total time: "<<sw.reading_str()<<"\n";

  bool run_test = false;
  if(run_test){
    std::cout<<"Test orthonormality [should all read 0]:\n";
    std::cout<<"       ";
    for(uint b=0; b<wf.nlist.size(); b++)
      printf("   %1i %2i  ",wf.nlist[b],wf.kappa[b]);
    std::cout<<"\n";
    for(uint a=0; a<wf.nlist.size(); a++){
      printf("%1i %2i  ",wf.nlist[a],wf.kappa[a]);
      for(uint b=0; b<wf.nlist.size(); b++){
        if(wf.kappa[a]!=wf.kappa[b]){
          std::cout<<" ------- ";
          continue;
        }
        double xf = INT::integrate3(wf.p[a],wf.p[b],wf.drdt);
        double xg = INT::integrate3(wf.q[a],wf.q[b],wf.drdt);
        double xo = wf.h*(xf+xg);
        if(wf.nlist[a]==wf.nlist[b]) xo -= 1;
        printf(" %7.0e ",xo);
      }
      std::cout<<"\n";
    }
  }

  return 0;
}
