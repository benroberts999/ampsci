#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>

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

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("parametricPotential.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> Gf >> Gh >> Gd;        getline(ifs,jnk);
    ifs >> Tf >> Tt >> Tg;        getline(ifs,jnk);
    ifs.close();
  }

  //Normalise the Teitz/Green weights:
  {
    double TG_norm = Gf + Tf;
    Gf /= TG_norm;
    Tf /= TG_norm;
  }

  //If H,d etc are zero, call default values??

  int Z = ATI_get_z(Z_str);
  if(Z==0) return 2;

  printf("\nRunning for Parametric potential potential, Z=%i\n",Z);
  printf("*************************************************\n");

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  //if(A!=0) wf.sphericalNucleus();
  A=133;
  //wf.fermiNucleus(); //Blah!

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[wf.ngp-1]);
  // XXX Make a "print grid info" function!

  wf.vdir.resize(wf.ngp);
  for(int i=0; i<wf.ngp; i++){
    double tmp = 0;
    if(Gf!=0) tmp += Gf*PRM_green(Z,wf.r[i],Gh,Gd);
    if(Tf!=0) tmp += Tf*PRM_tietz(Z,wf.r[i],Tt,Tg);
    wf.vdir[i] = tmp;
  }


  // * Form the core
  // * Loop over core, solve
  // * Solve valence list

  std::vector<int> core_list; //should be in the class!
  core_list = ATI_core_Xe; //Cs core is Xe-like!

  // Solve for each core state:
  int tot_el=0; // for working out Z_eff
  for(size_t i=0; i<core_list.size(); i++){
    int num = core_list[i];
    if(num==0) continue;

    int n = ATI_core_n[i];
    int l = ATI_core_l[i];

    tot_el+=num;

    double Zeff = 2. + double(Z - tot_el);
    double en_a = -0.5 * pow(Zeff/n,2); //energy guess

    int k1 = l; //j = l-1/2
    if(k1!=0) wf.solveLocalDirac(n,k1,en_a);
    int k2 = -(l+1); //j=l+1/2
    wf.solveLocalDirac(n,k2,en_a);

  }
  //return 1;

  int num_core = wf.nlist.size();

  printf("\n n l_j    k Rinf its    eps      En (au)\n");
  for(size_t i=0; i<wf.nlist.size(); i++){
    if((int)i==num_core) std::cout<<"========= Valence: ======\n";
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni);
  }

  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}
