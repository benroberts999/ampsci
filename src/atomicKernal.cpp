#include "AKF_akFunctions.h"
#include "ElectronOrbitals.h"
#include "ContinuumOrbitals.h"
#include "ATI_atomInfo.h"
#include "PRM_parametricPotentials.h"
#include "HartreeFockClass.h"
#include "FPC_physicalConstants.h"
#include "ChronoTimer.h"
#include "ExponentialGrid.h"
#include <iostream>
#include <fstream>
#include <cmath>

//******************************************************************************
int main(){
  ChronoTimer sw(true); //start stopwatch

  //Input options
  std::string Z_str;
  int A;
  double r0,rmax;
  int ngp;
  double varalpha;// for non-relativistic approx
  double hart_del;//HART convergance
  std::string str_core; //States for the core
  //Green potential parameters
  int Gf;
  double Gh,Gd;
  //q and dE grids:
  double qmin,qmax,demin,demax;
  int qsteps,desteps;
  // Max anglular momentums
  int max_l,max_L;
  int iout_format;
  std::string label; //label for output file

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("atomicKernal.in");
    std::string jnk;
    // read in the input parameters:
    ifs >> Z_str >> A;            getline(ifs,jnk);
    ifs >> str_core;              getline(ifs,jnk);
    ifs >> r0 >> rmax >> ngp;     getline(ifs,jnk);
    ifs >> Gf >> Gh >> Gd;        getline(ifs,jnk);
    ifs >> hart_del;              getline(ifs,jnk);
    ifs >> varalpha;              getline(ifs,jnk);
    ifs >>demin>>demax>>desteps;  getline(ifs,jnk);
    ifs >>qmin>> qmax >> qsteps;  getline(ifs,jnk);
    ifs >> max_l;                 getline(ifs,jnk);
    ifs >> max_L;                 getline(ifs,jnk);
    ifs >> iout_format;           getline(ifs,jnk);
    ifs >> label;                 getline(ifs,jnk);
    ifs.close();
  }
  if(Gf<0){
    Gf = 0;
  }

  // If L<0, will use plane-waves (instead of cntm fns)
  bool plane_wave = (max_L<0) ? true : false;

  //default Hartree convergance goal:
  if(hart_del==0) hart_del=1.e-6;

  //allow for single-step in dE or q grid
  if(desteps==1) demax=demin;
  if(qsteps==1) qmax=qmin;
  //Set up the E and q grids
  ExpGrid Egrid(desteps, demin, demax);
  ExpGrid qgrid(qsteps, qmin, qmax);

  // Fix maximum angular momentum values:
  if(max_l<0 || max_l>3) max_l=3; //default: all core states (no >f)
  if(plane_wave) max_L = max_l; //for spherical bessel.

  //alpha can't be zero, just make v. small
  if(varalpha==0) varalpha=1.e-25;

  //Convert units for input q and dE range into atomic units
  double keV = (1.e3/FPC::Hartree_eV);
  demin*=keV;
  demax*=keV;
  double qMeV = (1.e6/(FPC::Hartree_eV*FPC::c));
  qmin*=qMeV;
  qmax*=qMeV;

  //Look-up atomic number, Z, and also A
  int Z = ATI::get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI::A[Z]; //if none given, get default A

  //outut file name (excluding extension):
  std::string fname = "ak-"+Z_str+"_"+label;

  //Write out as text and/or binary file
  bool text_out = (iout_format==1)? false : true;
  bool bin_out  = (iout_format>0) ? true : false;

  //Print some info to screen:
  printf("\nRunning Atomic Kernal for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");
  if(Gf!=0) printf("Using Green potential: H=%.4f  d=%.4f\n",Gh,Gd);
  else printf("Using Hartree Fock (converge to %.0e)\n",hart_del);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);

  //Make sure h (large-r step size) is small enough to
  //calculate (normalise) cntm functions with energy = demax
  double h_target = (M_PI/20.)/sqrt(2.*demax);
  if(wf.h > h_target){
    int old_ngp = ngp;
    wf.logLinearRadialGrid(h_target,r0,rmax);
    ngp = wf.ngp;
    std::cout<<"\nWARNING 101: Grid not dense enough for contimuum state with "
         <<"ec="<<demax<<"au\n";
    std::cout<<"Updateing ngp: "<<old_ngp<<" --> "<<ngp<<"\n";
  }
  printf("Grid: pts=%i h=%6.4f r0=%.0e Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[0]
    ,wf.r[wf.ngp-1]);

  //Do Hartree-fock (or parametric potential) for Core
  if(Gf==0){
    HartreeFock hf(wf,str_core,hart_del);
  }else{
    //Use Green (local parametric) potential
    //Fill the electron part of the (local/direct) potential
    for(int i=0; i<wf.ngp; i++)
      wf.vdir.push_back(PRM::green(Z,wf.r[i],Gh,Gd));
    wf.solveInitialCore(str_core);//solves w/ Green
  }

  //make list of energy indices in sorted order:
  std::vector<int> sorted_by_energy_list;
  wf.sortedEnergyList(sorted_by_energy_list);

  //Output results:
  printf("\n     n l_j    k Rinf its    eps     En (au)     En (/cm)    En (eV)    Oc.Frac.\n");
  for(int i : sorted_by_energy_list){
    int n=wf.nlist[i];
    int k=wf.kappa[i];
    int twoj = ATI::twoj_k(k);
    int l = ATI::l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    double x = wf.occ_frac[i];
    printf("%2i) %2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f   (%.2f)\n",
      i,n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
      eni, eni*FPC::Hartree_invcm, eni*FPC::Hartree_eV,x);
  }

  //////////////////////////////////////////////////

  //Arrays to store results for outputting later:
  std::vector< std::vector< std::vector<float> > > AK; //float ok?
  int num_states = (int) wf.nlist.size();
  AK.resize(desteps, std::vector< std::vector<float> >
    (num_states, std::vector<float>(qsteps)));

  std::vector<std::string> nklst;//human-readiable state labels (easy plotting)

  //pre-calculate the spherical Bessel function look-up table for efficiency
  std::vector< std::vector< std::vector<float> > > jLqr_f;
  AKF::sphericalBesselTable(jLqr_f,max_L,qgrid,wf.r);

  //Calculate the AK
  std::cout<<"\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)\n"
        ,demin/keV,demax/keV,demin,demax);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)\n"
        ,qmin/qMeV,qmax/qMeV,qmin,qmax);

  //Store state info (each orbital) [just useful for plotting!]
  for(size_t is=0; is<wf.nlist.size(); is++){
    int k = wf.kappa[is];
    int l = ATI::l_k(k);
    //if(l>max_l) continue;
    int twoj = ATI::twoj_k(k);
    int n=wf.nlist[is];
    std::string nk=std::to_string(n)+ATI::l_symbol(l)
        +"_{"+std::to_string(twoj)+"/2}";
    nklst.push_back(nk);
  }

  //Calculate K(q,E)
  std::cout<<"Running dE loops ("<<desteps<<").."<<std::flush;
  #pragma omp parallel for
  for(int ide=0; ide<desteps; ide++){
    double dE = Egrid.x(ide);
    //Loop over core (bound) states:
    // for(size_t is=0; is<wf.nlist.size(); is++){
    for(auto is : wf.stateIndexList){
      int k = wf.kappa[is];
      int l = ATI::l_k(k); //XXX check this!
      if(l>max_l) continue;
      if(plane_wave) AKF::calculateKpw_nk(wf,is,dE,jLqr_f[l],AK[ide][is]);
      else AKF::calculateK_nk(wf,is,max_L,dE,jLqr_f,AK[ide][is]);
    }// END loop over bound states
  }
  std::cout<<"..done :)\n";

  //Write out to text file (in gnuplot friendly form)
  if(text_out) AKF::writeToTextFile(fname,AK,nklst,qmin,qmax,demin,demax);
  // //Write out AK as binary file
  if(bin_out) AKF::akReadWrite(fname,true,AK,nklst,qmin,qmax,demin,demax);
  std::cout<<"Written to: "<<fname;
  if(text_out) std::cout<<".txt";
  if(text_out && bin_out) std::cout<<", and ";
  if(bin_out) std::cout<<".bin";
  std::cout<<"\n";

  std::cout<<"\n "<<sw.reading_str()<<"\n";
  return 0;
}
