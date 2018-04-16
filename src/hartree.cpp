#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "ContinuumOrbitals.h"

double enGuess(int Z, int n, int l, int tot_el, int num);
int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, int a=-1);

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
  if(A==0) A=ATI_a[Z]; //if none given, get default A

  //Get default values for Green potential
  //NB: scaled to fit alkali _valence_ states - better if scaled for core!
  double Gh,Gd;  //Green potential parameters
  PRM_defaultGreen(Z,Gh,Gd);

  printf("\nRunning HART for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");
  printf("Initial: Green potential: H=%.4f  d=%.4f\n",Gh,Gd);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  if(A>0) wf.sphericalNucleus();

  printf("Grid: pts=%i h=%7.5f r0=%.1e Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[0],wf.r[wf.ngp-1]);

  //Determine which states are in the core:
  std::vector<int> core_list; //should be in the class!
  int core_ok = wf.determineCore(str_core,core_list);
  if(core_ok==2){
    std::cout<<"Problem with core: ";
    for(size_t i=0; i<str_core.size(); i++) std::cout<<str_core[i]<<" ";
    std::cout<<"\n";
    return 1;
  }

  //Fill the electron part of the potential, using Greens PRM for initial approx
  wf.vdir.resize(wf.ngp);
  for(int i=0; i<wf.ngp; i++){
    wf.vdir[i] = PRM_green(Z,wf.r[i],Gh,Gd);
  }

  //First step: Solve each core state using parameteric potential
  // XXX clean, put in other function
  int tot_el=0; // for working out Z_eff
  for(size_t i=0; i<core_list.size(); i++){
    int num = core_list[i];
    if(num==0) continue;
    int n = ATI_core_n[i];
    int l = ATI_core_l[i];
    double en_a = enGuess(Z,n,l,tot_el,num);
    tot_el+=num;
    int k1 = l; //j = l-1/2
    if(k1!=0) {
      wf.solveLocalDirac(n,k1,en_a,4);
      en_a = 0.95*wf.en[wf.nlist.size()-1]; //update guess for next same l
    }
    int k2 = -(l+1); //j=l+1/2
    if(num>2*l) wf.solveLocalDirac(n,k2,en_a,4);
  }


  // for(int i=0; i<wf.nlist.size(); i++){
  //   double norm=0;
  //   for(int j=0; j<wf.ngp; j++){
  //     norm += (pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j]*wf.h;
  //   }
  //   std::cout<<norm<<"\n";
  // }
  // return 1;





for(int n=0; n<25; n++){

  double eta=0.3;

  std::vector<double> vdir_new;
  formNewVdir(wf,vdir_new);
  std::vector<double> vdir_old = wf.vdir;
  for(int j=0; j<wf.ngp; j++){
    wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
  }


  double prev_e = wf.en[0];

  for(size_t i=0; i<wf.nlist.size(); i++){

    // std::vector<double> vdir_new;
    // formNewVdir(wf,vdir_new,i);
    // std::vector<double> vdir_old = wf.vdir;
    // for(int j=0; j<wf.ngp; j++){
    //   wf.vdir[j] = eta*vdir_new[j] + (1.-eta)*vdir_old[j];
    // }


    double del_e=0;
    for(int j=0; j<wf.ngp; j++)
      del_e += (wf.vdir[j]-vdir_old[j])*
      (pow(wf.p[i][j],2) + pow(wf.q[i][j],2))*wf.drdt[j];
    del_e*=wf.h;
    double new_e = wf.en[i] + 1*del_e;
    if(new_e>0)new_e=-0.1;
    //std::cout<<wf.en[i]<<" "<<del_e<<" "<<wf.en[i]+del_e<<"\n";
    wf.reSolveLocalDirac(i,new_e,10);
  }

  double next_e = wf.en[0];

  std::cout<<n+1<<" "<<(next_e-prev_e)/(next_e*eta)<<"\n";
  if(fabs((next_e-prev_e)/(next_e*eta))<1.e-5) break;

}



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
int formNewVdir(ElectronOrbitals wf, std::vector<double> &vdir_new, int a)
//XXX OK? lots memory? Reference?
//XXX core list!!! make part of class!
{

  vdir_new.clear();
  vdir_new.resize(wf.ngp);

  //XXX use core_list instead!!! XXX
  int Ncore=0;
  for(size_t i=0; i<wf.nlist.size(); i++){
    int ka = wf.klist[i];
    int twoj = 2*abs(ka)-1;
    Ncore += twoj+1;
  }


  // //a=-1 means assume vdir same for all orbitals!
  double f=1;
  if(a==-1) f = 1. - 1./Ncore;

  std::vector<double> rho(wf.ngp);
  for(size_t i=0; i<wf.nlist.size(); i++){
    int ia=0;
    //if ((int)i==a) ia=0; //number of states to 'skip'
    int ka = wf.klist[i];
    //int la = (abs(2*ka+1)-1)/2;
    int twoj = 2*abs(ka)-1;
    for(int j=0; j<wf.ngp; j++){
      rho[j] += (twoj+1-ia)*(pow(wf.p[i][j],2) + pow(wf.q[i][j],2));
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

//******************************************************************************
double enGuess(int Z, int n, int l, int tot_el, int num)
{
  //effective Z (for energy guess) -- not perfect!
  double Zeff =  double(Z - tot_el - num);
  if(l==1) Zeff = 1. + double(Z - tot_el - 0.5*num);
  if(l==2) Zeff = 1. + double(Z - tot_el - 0.5*num);
  if(Zeff<1.) Zeff=1.;

  double en_a = -0.5 * pow(Zeff/n,2);
  if(n>1) en_a *= 0.5;
  return en_a;
}
