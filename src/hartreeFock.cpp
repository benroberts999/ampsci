#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "HF_hartreeFock.h"
#include <iostream>
#include <sstream>
#include <sys/time.h>

int main(void){

  struct timeval start, end;
  gettimeofday(&start, NULL);

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;
  double varalpha;
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
    ifs >> varalpha;              getline(ifs,jnk);
    ifs.close();
  }
  bool doHF = true;
  if(iHF==0) doHF = false;

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

  //Solve Hartree equations for the core:
  if(doHF) HF::hartreeFockCore(wf,eps_hart);
  else     HF::hartreeCore(wf,eps_hart);

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
  for(uint i=0; i<lst.size(); i++){
    int n = lst[i][0];
    int k = lst[i][1];
    if(doHF) HF::hartreeFockValence(wf,n,k,eps_hart);
    else{
      double en_g = wf.enGuessVal(n,k);
      wf.solveLocalDirac(n,k,en_g,15);
    }
  }

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

  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.3f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.2f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  bool run_test = false;
  if(run_test){
    std::vector<double> ppqq(wf.ngp);
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
        double x1=0;
        for (int i=0; i<wf.ngp; i++){
          double y = wf.p[a][i]*wf.p[b][i]+wf.q[a][i]*wf.q[b][i];
          ppqq[i] = y;
          x1 += y*wf.drdt[i];
        }
        x1 *= wf.h;
        double x2 = INT::integrate(ppqq,wf.drdt,wf.h);
        double xo=x2;
        if(wf.nlist[a]==wf.nlist[b]) xo -= 1;
        printf(" %7.0e ",xo);
      }
      std::cout<<"\n";
    }
  }

  return 0;
}
