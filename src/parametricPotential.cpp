#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include <iostream>
// #include <fstream>

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

  int Z = ATI_get_z(Z_str);
  if(Z==0) return 2;

  printf("\nRunning for Parametric potential potential, Z=%i\n",Z);
  printf("*************************************************\n");

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  //if(A!=0) wf.sphericalNucleus();

  printf("Grid: pts=%i h=%7.5f Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[wf.ngp-1]);
  // XXX Make a "print grid info" function!

  //Solve the Dirac equation for H-like ions:
  // wf.hydrogenLike(n_max,l_max);


  // * Form the core
  // * Loop over core, solve
  // * Solve valence list

  std::vector<int> core_list; //should be in the class!
  core_list = ATI_core_Xe; //Cs core is Xe-like!

  for(size_t i=0; i<core_list.size(); i++){
    int num = core_list[i];
    if(num==0) continue;

    int n = ATI_core_n[i];
    int l = ATI_core_n[i];

    double en_a = -0.5 * pow(n/6.,2); //energy guess

    k = l //j = l-1/2
    if(k!=0) wf.solveLocalDirac(n,k,en_a);
    k = -(l+1); //j=l+1/2
    wf.solveLocalDirac(n,k,en_a);

  }




  tf = clock();
  double total_time = 1000.*double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.3f ms.\n",total_time);

  return 0;
}
