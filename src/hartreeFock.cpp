#include "ElectronOrbitals.h"
#include "HartreeFockClass.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include "FPC_physicalConstants.h"
#include "ATI_atomInfo.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include "ChronoTimer.h"

int main(){

  ChronoTimer timer; //start the stopwatch
  timer.start();

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;
  double varalpha,varalpha2;
  double eps_HF;

  int num_val,l_max;
  std::string str_core;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("hartreeFock.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    ifs >> str_core;              getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> eps_HF;                getline(ifs,jnk);
    ifs >> num_val >> l_max;      getline(ifs,jnk);
    ifs >> varalpha2;             getline(ifs,jnk);
    ifs.close();
  }

  //Change varAlph^2 to varalph
  if(varalpha2==0) varalpha2 = 1.e-10;
  varalpha = sqrt(varalpha2);

  int Z = ATI::get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI::A[Z]; //if none given, get default A

  printf("\nRunning HARTREE FOCK for %s, Z=%i A=%i\n",Z_str.c_str(),Z,A);
  printf("*************************************************\n");

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);

  printf("Grid: pts=%i h=%7.5f r0=%.1e Rmax=%5.1f\n\n"
  ,  wf.ngp,wf.h,wf.r[0],wf.r[wf.ngp-1]);

  timer.start(); //start the timer

  //Solve Hartree equations for the core:
  HartreeFock hf(wf,str_core,eps_HF);
  double core_energy = hf.calculateCoreEnergy();
  std::cout<<"core: "<<timer.lap_reading_str()<<"\n";

  //Make a list of valence states to solve:
  //Goes from n=0 -> n_max (inclusive), same for l
  //Adds any state to 'list to calculate' that isn't already in core
  // if((int)wf.num_core_electrons >= wf.Znuc()) n_max = 0;
  // std::vector< std::vector<int> > lst;
  // for(int n=0; n<=n_max; n++){
  //   for(int l=0; l<n; l++){
  //     if(l>l_max) continue;
  //     for(int tk=0; tk<2; tk++){ //loop over k
  //       int k;
  //       if(tk==0) k=l;      //j = l - 1/2
  //       else      k=-(l+1); //j = l + 1/2
  //       if(k==0) continue;  // no j = l - 1/2 for l=0
  //       if(wf.isInCore(n,k)) continue;
  //       lst.push_back({n,k});
  //     }
  //   }
  // }

  //Create list of valence states to solve for
  // Solves for lowest num_val states with given l
  std::vector< std::vector<int> > lst;
  if((int)wf.num_core_electrons >= wf.Znuc()) num_val = 0;
  for(int l=0; l<=l_max; l++){
    int n0 = wf.maxCore_n(l) + 1;
    if(n0==1) n0 += l;
    for(int nv=0; nv<num_val; nv++){
      for(int tk=0; tk<2; tk++){ //loop over k
        int ka = (tk==0)? l : -(l+1); //kappa
        if(ka==0) continue; //no j=l-1/2 for l=s
        lst.push_back({n0+nv,ka});
      }
    }
  }

  //Solve for the valence states:
  timer.start();
  for(const auto &nk : lst){
    int n = nk[0];
    int k = nk[1];
    hf.solveValence(n,k);
  }
  if(lst.size()>0) std::cout<<"Valence: "<<timer.lap_reading_str()<<"\n";

  //make list of energy indices in sorted order:
  std::vector<int> sorted_by_energy_list;
  wf.sortedEnergyList(sorted_by_energy_list,true);

  //Output results:
  printf("\nCore: %s, Z=%i A=%i\n",Z_str.c_str(),Z,A);
  printf("     state   k   Rinf its    eps       En (au)      En (/cm)\n");
  bool val=false;
  double en_lim=0;
  for(int i : sorted_by_energy_list){
    if(val && en_lim==0) en_lim = fabs(wf.en[i]); //give energies wrt core
    int k=wf.kappa[i];
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    std::string symb = wf.seTermSymbol(i);
    printf("%2i) %7s %2i  %5.1f %3i  %5.0e %13.7f %13.1f",
        i,symb.c_str(),k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*FPC::Hartree_invcm);
    if(val)printf(" %10.2f\n",(eni+en_lim)*FPC::Hartree_invcm);
    else std::cout<<"\n";
    if(i==(int)wf.num_core_states-1){
      printf("E_core = %.5f au\n",core_energy);
      if(wf.num_core_states==wf.nlist.size()) break;
      std::cout<<"Val: state   "
        <<"k   Rinf its    eps       En (au)      En (/cm)   En (/cm)\n";
      val = true;
    }
  }

  std::cout<<"\n Total time: "<<timer.reading_str()<<"\n";

  bool run_test = false;
  if(run_test){
    std::cout<<"Test orthonormality [should all read 0]:\n";
    std::cout<<"       ";
    for(uint b=0; b<wf.nlist.size(); b++)
      printf("   %1i %2i  ",wf.nlist[b],wf.kappa[b]);
    std::cout<<"\n";
    for(uint a=0; a<wf.nlist.size(); a++){
      printf("%1i %2i  ",wf.nlist[a],wf.kappa[a]);
      for(uint b=0; b<wf.nlist.size(); b++){
        if(wf.kappa[a]!=wf.kappa[b]){
          std::cout<<" ------- ";
          continue;
        }
        double xf = INT::integrate3(wf.f[a],wf.f[b],wf.drdt);
        double xg = INT::integrate3(wf.g[a],wf.g[b],wf.drdt);
        double xo = wf.h*(xf+xg);
        if(wf.nlist[a]==wf.nlist[b]) xo -= 1;
        printf(" %7.0e ",xo);
      }
      std::cout<<"\n";
    }
    std::cout<<"\n Total time: "<<timer.reading_str()<<"\n";
  }

  bool print_wfs = false;
  if(print_wfs){
    std::ofstream of("hf-orbitals.txt");
    of<<"r ";
    for(size_t a=0; a<wf.nlist.size(); a++){
      of<<"\""<<wf.seTermSymbol(a)<<"\" ";
    }
    of<<"\n";
    for(int i=0; i<wf.ngp; i++){
      of<<wf.r[i]<<" ";
      for(size_t a=0; a<wf.nlist.size(); a++){
        of<<wf.f[a][i]<<" ";
      }
      of<<"\n";
    }
    std::cout<<"\n Total time: "<<timer.reading_str()<<"\n";
  }


  bool testpnc = false;
  if(testpnc){
    double a=0.22756*2.3;
    double c=5.6710; //1.1*pow(wf.A(),0.3333);

    double rho0 = 0;
    std::vector<double> rho(wf.ngp);
    for(size_t i=0; i<rho.size(); i++){
      rho[i] = 1./(1. + exp((wf.r[i]*FPC::aB_fm-c)/a));
      rho0 += pow(wf.r[i],2)*rho[i]*wf.drdt[i]*4*M_PI*wf.h;
    }
    for(size_t i=0; i<rho.size(); i++){
      rho[i]/=rho0;
    }

    double Gf = FPC::GFe11;
    double Cc = (Gf/sqrt(8))*(-133+55); // Qw/(-N)
    double Ac = 2./6;

    int a6s = wf.getStateIndex(6,-1);
    int a7s = wf.getStateIndex(7,-1);
    std::cout<<a6s<<","<<a7s<<"\n";
    double pnc=0;
    for(auto np=0ul; np<wf.kappa.size(); np++){
      if(wf.kappa[np]!=1) continue; //p_1/2 only
      int n = wf.nlist[np];
      // <7s|d|np><np|hw|6s>/dE6s + <7s|hw|np><np|d|6s>/dE7s
      double d7s = wf.radialIntegral(a7s,np,wf.r);
      double w6s = wf.radialIntegral(np,a6s,rho,Operator::gamma5);
      double dE6s = wf.en[a6s]-wf.en[np];
      double d6s = wf.radialIntegral(np,a6s,wf.r);
      double w7s = wf.radialIntegral(a7s,np,rho,Operator::gamma5);
      double dE7s = wf.en[a7s]-wf.en[np];
      double pnc1 = Cc*Ac*d7s*w6s/dE6s;
      double pnc2 = Cc*Ac*d6s*w7s/dE7s;
      std::cout<<"n="<<n<<" pnc= "<<pnc1<<" + "<<pnc2<<" = "<<pnc1+pnc2<<"\n";
      // std::cout<<" d: "<<d6s*2*0.408<<" "<<d7s*2*0.408<<"\n";
      // std::cout<<" w: "<<w6s*sqrt(2)*Cc<<" "<<w7s*sqrt(2)*Cc<<"\n";
      pnc += pnc1+pnc2;
    }
    std::cout<<"Total= "<<pnc<<"\n";
    std::cout<<"\n Total time: "<<timer.reading_str()<<"\n";
  }

  return 0;
}
