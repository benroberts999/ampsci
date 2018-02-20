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

  double H=4.4691;
  double d=0.8967;
  for(int i=0; i<wf.ngp; i++){
    wf.vnuc[i] += PRM_green(Z,wf.r[i],H,d);
  }


  int k =-1;
  int n=6;
  wf.nlist.push_back(n);
  wf.klist.push_back(k);
  int pinf,its;
  double eps;
  double en_a = -0.5*pow(1./n,2);
  std::vector<double> p_a(ngp);
  std::vector<double> q_a(ngp);
  solveDBS(p_a,q_a,en_a,wf.vnuc,Z,n,k,wf.r,wf.drdt,wf.h,wf.ngp,pinf,its,eps,wf.alpha);
  wf.p.push_back(p_a);
  wf.q.push_back(q_a);
  wf.en.push_back(en_a);
  //store convergance info:
  wf.pinflist.push_back(pinf);
  wf.itslist.push_back(its);
  wf.epslist.push_back(eps);



  printf(" n l_j    k  R_inf its eps     En (au)            En (/cm)\n");
  //int num_states =
  for(size_t i=0; i<wf.nlist.size(); i++){
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    printf("%2i %s_%i/2 (%2i)  %3.0f %3i  %5.0e  %.15f  %.7f\n",
        n,atinfo_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        wf.en[i],wf.en[i]*HARTREE_ICM);
  }

  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}
