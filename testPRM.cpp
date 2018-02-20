#include "ElectronOrbitals.h"
#include "atomInfo.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>

int main(void){

  clock_t ti,tf;
  ti = clock();

  int Z=55,A=133;
  //int n_max=7,l_max=1;
  int ngp = 2000;
  double r0=1.e-6,rmax=250.;
  double varalpha=1.;
  // std::ifstream ifile;
  // ifile.open("h-like.in");
  // {
  //   std::string junk;
  //   ifile >> Z >> A;            getline(ifile,junk);
  //   ifile >> n_max >> l_max;    getline(ifile,junk);
  //   ifile >> r0 >> rmax >> ngp; getline(ifile,junk);
  //   ifile >> varalpha;          getline(ifile,junk);
  // }
  // ifile.close();

  printf("\nRunning BLAH potential, Z=%i\n",Z);
  printf("*************************************************\n");

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  wf.sphericalNucleus();
  //wf.fermiNucleus();

  double H=4.4691;
  double d=0.8967;
  for(int i=0; i<wf.ngp; i++){
    double vtmp = PRM_green(Z,wf.r[i],H,d);
    wf.vdir.push_back(vtmp);
  }


  int max_l=1;
  int max_n=8;
  for(int n=6; n<=max_n; n++){
    for(int i=1; i<2*n; i++){ //loop through each kappa state
      int k = pow(-1,i)*ceil(0.5*i);
      int l = (abs(2*k+1)-1)/2;
      if(l>max_l) continue;
      //
      double e_a = -0.5*pow(1./n,2);
      wf.solveLocalDirac(n,k,e_a);
    }
  }



  printf(" n l_j    k  R_inf its eps     En (au)            En (/cm)\n");
  //int num_states =
  for(size_t i=0; i<wf.nlist.size(); i++){
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    double en0 = wf.en[0];
    printf("%2i %s_%i/2 (%2i)  %3.0f %3i  %5.0e  %.15f  %.7f\n",
        n,atinfo_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        wf.en[i],(wf.en[i]-en0)*HARTREE_ICM);
  }

  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}
