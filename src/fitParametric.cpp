#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <omp.h>
/*
Finds the best-fit parameter values for the Green or Tietz potentials
*/

int main(void){

  clock_t ti,tf;
  ti = clock();

  bool sphere=true; //Finite nucleus? makes no difference??

  //Input parameters:
  std::string Z_str;
  int A;
  int ngp;
  double r0,rmax;
  int igreen;
  double varalpha=1.;
  std::vector<int> in_n,in_k;
  std::vector<double> in_en;
  std::ifstream ifile;
  ifile.open("fitParametric.in"); //input file
  {
    std::string junk;
    ifile >> Z_str >> A;        getline(ifile,junk);
    ifile >> r0 >> rmax >> ngp; getline(ifile,junk);
    ifile >> igreen;            getline(ifile,junk);
    int nstates;
    ifile >> nstates;           getline(ifile,junk);
    for(int ns=0; ns<nstates; ns++){
      int tn,tk;
      double te;
      ifile >> tn >> tk >> te; getline(ifile,junk);
      in_n.push_back(tn);
      in_k.push_back(tk);
      in_en.push_back(te);
    }
  }
  ifile.close();

  int Z = ATI::get_z(Z_str);
  if(Z==0) return 2;

  bool green=true;
  std::string which="Green";
  if(igreen==1){
    green=false;
    which = "Tietz";
  }

  //convergence parameters for finding best-fit H and d (or t and g)
  double eps = 1.e-5;
  int max_its = 100;

  printf("\n Finding best-fit parameters for %s potential, %s Z=%i\n",
    which.c_str(),Z_str.c_str(),Z);
  printf("*********************************************************\n");

  double GHmin=0.1, GHmax=10.;
  double Gdmin=0.05, Gdmax=2.;

  double Hmin=GHmin, Hmax=GHmax;
  double dmin=Gdmin, dmax=Gdmax;

  double best_H=0,best_d=0;

  for(int nit=0; nit<max_its; nit++){
    const int n_array=32;
    const int n_params=3;
    double array[n_array][n_array][n_params];
    double dH = (Hmax-Hmin)/n_array;
    double dd = (dmax-dmin)/n_array;

    //Solve equation for each (H,d), store values
    #pragma omp parallel for
    for(int n=0; n<n_array; n++){
      double H=Hmin+n*dH;
      for(int m=0; m<n_array; m++){
        double d=dmin+m*dd;
        ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
        if(sphere) wf.sphericalNucleus();
        if(green) for(int i=0; i<wf.ngp; i++)
          wf.vdir.push_back(PRM::green(Z,wf.r[i],H,d));
        else for(int i=0; i<wf.ngp; i++)
          wf.vdir.push_back(PRM::tietz(Z,wf.r[i],H,d));
        double fx = 0;
        for(size_t ns=0; ns<in_n.size(); ns++){
          wf.solveLocalDirac(in_n[ns],in_k[ns],in_en[ns]);
          fx += pow((wf.en[ns]-in_en[ns]),2);
          //fx += pow((wf.en[ns]-in_en[ns])/in_en[ns],2);
        }
        array[n][m][0]=fx;
        array[n][m][1]=H;
        array[n][m][2]=d;
      }
    }

    //Find the "so-far" best-fit
    double bH=0,bd=0,bfx=99999999.;
    for(int n=0; n<n_array; n++){
      for(int m=0; m<n_array; m++){
        if(array[n][m][0]<bfx){
          bfx=array[n][m][0];
          bH=array[n][m][1];
          bd=array[n][m][2];
        }
      }
    }

    //Adjust H and d search range (hone-in):
    Hmin = bH - (0.2*n_array)*dH;
    Hmax = bH + (0.2*n_array)*dH;
    dmin = bd - (0.2*n_array)*dd;
    dmax = bd + (0.2*n_array)*dd;
    if(Hmax>GHmax) Hmax=GHmax;
    if(dmax>Gdmax) dmax=Gdmax;
    if(Hmin<GHmin) Hmin=GHmin;
    if(dmin<Gdmin) dmin=Gdmin;
    best_H = bH;
    best_d = bd;

    //Print progress to screen; quit if convergence OK
    printf("%2i %6.4f %6.4f  %.1e\n",nit,bH,bd,fmax(dH,dd));
    if(dH<eps && dd<eps) break;
  }

  double H=best_H;
  double d=best_d;
  printf("\nBest fit parameters for ");
  if(green) printf("Green: \n  H=%7.5f  d=%7.5f\n\n",H,d);
  else printf("Tietz: \n  t=%7.5f  g=%7.5f\n\n",H,d);

  //Now, solve using the above-found best-fit parameters:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  if(sphere) wf.sphericalNucleus();

  if(green)
    for(int i=0; i<wf.ngp; i++) wf.vdir.push_back(PRM::green(Z,wf.r[i],H,d));
  else
    for(int i=0; i<wf.ngp; i++) wf.vdir.push_back(PRM::tietz(Z,wf.r[i],H,d));
  for(size_t ns=0; ns<in_n.size(); ns++)
    wf.solveLocalDirac(in_n[ns],in_k[ns],in_en[ns]);

  printf(" n l_j    k Rinf its  eps     En (au)            En (/cm)\n");
  for(size_t i=0; i<wf.nlist.size(); i++){
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double rinf = wf.r[wf.pinflist[i]];
    double en0 = wf.en[0];
    double eni = wf.en[i];
    double enT = in_en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %.15f  %13.7f  %9.4f%%\n",
        n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni,(eni-en0)*HARTREE_ICM,100.*(enT-eni)/enT);
  }

  tf = clock();
  double total_time = double(tf-ti)/CLOCKS_PER_SEC;
  printf ("\nt=%.1f s.\n",total_time);

  return 0;
}
