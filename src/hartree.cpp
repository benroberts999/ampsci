#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "HF_hartree.h"
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

  int n_max,l_max;
  std::string str_core;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("hartree.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    ifs >> str_core;              getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> eps_hart;              getline(ifs,jnk);
    ifs >> n_max >> l_max;        getline(ifs,jnk);
    ifs >> varalpha;              getline(ifs,jnk);
    ifs.close();
  }

  int Z = ATI::get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI::A[Z]; //if none given, get default A

  printf("\nRunning HARTREE for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
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
  HF::hartreeCore(wf,eps_hart);

  int maxn=0; //max 'n' in the core (used for valence energy guess)
  for(int i=0; i<wf.num_core_states; i++) if(wf.nlist[i]>maxn) maxn=wf.nlist[i];

  //Calculate the valence (and excited) states
  for(int n=0; n<=n_max; n++){
    for(int l=0; l<=l_max; l++){ //loop over l
      if(l+1>n) continue;
      for(int tk=0; tk<2; tk++){ //loop over k (ie j=l +/- 1/2)
        int k;
        if(tk==0) k=l;      //j = l - 1/2
        else      k=-(l+1); //j = l + 1/2
        if(k==0) continue;  // no j = l - 1/2 for l=0
        if(wf.isInCore(n,k)) continue; //skip states already in the core
        //This energy guess works very well for Cs, Fr etc.
        //Works poorly (but still converges) for light atoms
        int dn=n-maxn;
        double neff=1.+dn;
        double x=1;
        if(maxn<4) x=0.25;
        if(l==1) neff+=0.5*x;
        if(l==2) neff+=2.*pow(x,0.5);
        if(l>=3) neff+=4.*x;
        double en_a = -0.5/pow(neff,2);
        wf.solveLocalDirac(n,k,en_a);
      }
    }
  }

  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\nHARTEE results for %s, Z=%i A=%i\n",Z_str.c_str(),Z,A);
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
    int k=wf.klist[i];
    int twoj = ATI::twoj_k(k);
    int l = ATI::l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %15.3f",
        n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*FPC::Hartree_invcm);
    if(val)printf(" %10.2f\n",(eni+en_lim)*FPC::Hartree_invcm);
    else std::cout<<"\n";
  }

  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.3f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.2f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  return 0;
}
