#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "HF_hartree.h"
#include <iostream>
#include <sstream>


int main(void){

  clock_t ti,tf;
  ti = clock();

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;
  double varalpha;

  int n_max,l_max;
  std::vector<std::string> str_core;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("hartree.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    while(true){
      std::string str;
      ifs >> str;
      if(str=="."||str=="|"||str=="!") break;
      str_core.push_back(str);
    }
    getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> n_max >> l_max;        getline(ifs,jnk);
    ifs >> varalpha;              getline(ifs,jnk);
    ifs.close();
  }

  int Z = ATI_get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI_a[Z]; //if none given, get default A


  double eps_hartree=1.e-6;



  printf("\nRunning HARTEE for %s, Z=%i A=%i\n",
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
    std::cout<<"Problem with core: ";
    for(size_t i=0; i<str_core.size(); i++) std::cout<<str_core[i]<<" ";
    std::cout<<"\n";
    return 1;
  }

  //Solve Hartree equations for the core:
  HF_hartreeCore(wf,eps_hartree);

  int maxn=0;
  for(int i=0; i<wf.num_core; i++) if(wf.nlist[i]>maxn) maxn=wf.nlist[i];

  //Calculate the valence (and excited) states
  for(int n=0; n<=n_max; n++){
    for(int l=0; l<=l_max; l++){
      if(l+1>n) continue;
      for(int tk=0; tk<2; tk++){
        int k;
        if(tk==0) k=l;
        else      k=-(l+1);
        if(k==0) continue;
        if(wf.isInCore(n,k)) continue;
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
  printf("\n n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    if((int)m==wf.num_core){
      std::cout<<" ========= Valence: ======\n";
      printf(" n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
    }
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %15.3f\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*HARTREE_ICM);
  }



  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}
