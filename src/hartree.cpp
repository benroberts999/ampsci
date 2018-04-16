#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "ContinuumOrbitals.h"

int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, bool core=true);

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

  int max_hartree=100; //Max number of Hartree iterations
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

  //Fill the electron part of the potential, using Greens PRM for initial approx
  double Gh,Gd;  //Green potential parameters
  PRM_defaultGreen(Z,Gh,Gd); //Get default values for Green potential
  for(int i=0; i<wf.ngp; i++) wf.vdir.push_back(PRM_green(Z,wf.r[i],Gh,Gd));

  //First step: Solve each core state using parameteric potential
  wf.solveInitialCore(1);


  //Hartree loop:
  for(int n=0; n<max_hartree; n++){

    double eta=0.50;

    std::vector<double> vdir_old = wf.vdir;
    std::vector<double> vdir_new;
    formNewVdir(wf,vdir_new);
    for(int j=0; j<wf.ngp; j++){
      wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
    }


    double prev_e = 0;
    for(size_t i=0; i<wf.nlist.size(); i++) prev_e += wf.en[i]/wf.nlist.size();

    for(size_t i=0; i<wf.nlist.size(); i++){
      double del_e=0;
      for(int j=0; j<wf.ngp; j++)
        del_e += (wf.vdir[j]-vdir_old[j])*
        (pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
      del_e*=wf.h;
      double new_e = wf.en[i] + 1*del_e;
      if(new_e>0)new_e=-0.1;
      wf.reSolveLocalDirac(i,new_e,3); //only go to 1/10^3 - do better at end!
    }

    double next_e = 0;
    for(size_t i=0; i<wf.nlist.size(); i++) next_e += wf.en[i]/wf.nlist.size();

    double delta_hartree = (next_e-prev_e)/(next_e*eta);
    printf("Hart it:%3i,  del=%6.0e\n",n+1,delta_hartree);

    if(fabs(delta_hartree)<eps_hartree) break;
  }

  //re-run solve Dirac to higher convergance level after Hart pot. ok
  for(size_t i=0; i<wf.nlist.size(); i++) wf.reSolveLocalDirac(i,0,14);

  formNewVdir(wf,wf.vdir,false);

  wf.solveLocalDirac(6,-1,-0.13);
  wf.solveLocalDirac(7,-1,-0.06);
  wf.solveLocalDirac(8,-1,-0.03);



  std::ofstream ofile;
  ofile.open("pot.txt");
  ofile<<"r Gr Vh Z/r 1/r\n";
  for(int i=0; i<wf.ngp; i++){
    ofile<<wf.r[i]<<" "
      <<-PRM_green(Z,wf.r[i],Gh,Gd)-wf.vnuc[i]<<" "
      <<-wf.vdir[i]-wf.vnuc[i]<<" "
      <<Z/wf.r[i]<<" "<<1./wf.r[i]<<"\n";

  }
  ofile.close();





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

//******************************************************************************
int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, bool core)
//XXX OK? lots memory? Reference?
//XXX core list!!! make part of class!
{

  vdir_new.clear();
  vdir_new.resize(wf.ngp);

  //XXX use core_list instead!!! XXX ??
  int Ncore=0;
  for(size_t i=0; i<wf.nlist.size(); i++){
    int ka = wf.klist[i];
    int twoj = 2*abs(ka)-1;
    Ncore += twoj+1;
  }


  // //a=-1 means assume vdir same for all orbitals!
  double f=1;
  if(core) f = 1. - (1.)/Ncore;

  std::vector<double> rho(wf.ngp);
  for(size_t i=0; i<wf.nlist.size(); i++){
    int ka = wf.klist[i];
    int twoj = 2*abs(ka)-1;
    for(int j=0; j<wf.ngp; j++){
      rho[j] += (twoj+1)*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2));
      //XXX assumes closed shell!? ia!
    }
  }

  for(int j=0; j<wf.ngp; j++){
    double r = wf.r[j];
    double v_tmp = 0;
    for(int k=0; k<wf.ngp; k++){
      double rp = wf.r[k];
      double rm = std::max(r,rp);//can be little more clever, slight speedup
      v_tmp += (rho[k]/rm)*wf.drdt[k];
    }
    vdir_new[j] = f*v_tmp*wf.h;
  }

  return 0;
}
