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
  HF_hartreeCore(wf);

  wf.solveLocalDirac(6,-1,-0.13);
  wf.solveLocalDirac(7,-1,-0.06);
  wf.solveLocalDirac(8,-1,-0.03);








  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\n n l_j    k Rinf its    eps      En (au)        En (/cm)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
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


// //******************************************************************************
// double formRho(ElectronOrbitals wf, std::vector<double> rho) //XXX OK? lots memory? Reference?
// {
//   rho.clear();
//   rho.resize(wf.ngp);
//   for(size_t i=0; i<wf.nlist.size(); i++){
//     int ka = wf.klist[i];
//     int la = (abs(2*ka+1)-1)/2;
//     for(int j=0; j<wf.ngp; j++){
//       rho[j] += (4*la+2)*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2));
//     }
//   }
// }
