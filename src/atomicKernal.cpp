#include "AKF_akFunctions.h"

//******************************************************************************
int main(void){

  struct timeval start, end;
  gettimeofday(&start, NULL);

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
  bool noExch = false;
  if(Gf<0){
    noExch = true;
    Gf = 0;
  }

  bool plane_wave = false;
  if(max_L<0) plane_wave = true;

  //default Hartree convergance goal:
  if(hart_del==0) hart_del=1.e-6;

  //allow for single-step in dE or q grid
  if(desteps==1) demax=demin;
  if(qsteps==1) qmax=qmin;

  // Fix maximum angular momentum values:
  if(max_l<0 || max_l>3) max_l=3; //default: all core states (no >f)
  if(plane_wave) max_L = max_l; //for spherical bessel.

  //alpha can't be zero:
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

  //Zeff: only for testing!
  bool use_Zeff = false; // XXX Make sure set to false!
  double en_zef = -43.;
  int n_zef = 3;
  int k_zef = -1;
  double Zeff = n_zef*sqrt(-2.*en_zef);

  //outut file name 9excluding extension):
  std::string fname = "ak-"+Z_str+"_"+label;

  //What to output?
  bool text_out = true;
  bool bin_out  = false;
  if(iout_format==1) text_out = false;
  if(iout_format>0)  bin_out  = true;

  //Print some info to screen:
  printf("\nRunning Atomic Kernal for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");
  if(Gf!=0) printf("Using Green potential: H=%.4f  d=%.4f\n",Gh,Gd);
  else if(noExch)printf("Using Hartree pot. (converge to %.0e)\n",hart_del);
  else printf("Using Hartree Fock (converge to %.0e)\n",hart_del);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  double h_target = (M_PI/20.)/sqrt(2.*demax);
  if(wf.h > h_target){
    int old_ngp = ngp;
    wf.DzubaRadialGrid(h_target,r0,rmax);
    ngp = wf.ngp;
    std::cout<<"\nWARNING 101: Grid not dense enough for contimuum state with "
         <<"ec="<<demax<<"au\n";
    std::cout<<"Updateing ngp: "<<old_ngp<<" --> "<<ngp<<"\n";
  }
  printf("Grid: pts=%i h=%6.4f r0=%.0e Rmax=%5.1f\n",wf.ngp,wf.h,wf.r[0]
         ,wf.r[wf.ngp-1]);

  //If non-zero A is given, use spherical nucleus.
  if(A>0) wf.sphericalNucleus();

  //Determine which states are in the core:
  int core_ok = wf.determineCore(str_core);
  if(core_ok==2){
    std::cout<<"Problem with core: "<<str_core<<"\n";
    return 1;
  }

  if(use_Zeff){
    //Use Zeff (single oribtal only, just for tests!)
    wf.solveZeff(n_zef,k_zef,Zeff,false);
  }else if(Gf==0){
    //use Hartree method:
    if(noExch) HF::hartreeCore(wf,hart_del);
    else       HF::hartreeFockCore(wf,hart_del);
  }else{
    //Use Green (local parametric) potential
    //Fill the electron part of the (local/direct) potential
    for(int i=0; i<wf.ngp; i++)
      wf.vdir.push_back(PRM::green(Z,wf.r[i],Gh,Gd));
    // Solve Dirac equation for each (bound) core state:
    wf.solveInitialCore();
  }


  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\n     n l_j    k Rinf its    eps     En (au)     En (/cm)    En (eV)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    int n=wf.nlist[i];
    int k=wf.kappa[i];
    int twoj = ATI::twoj_k(k);
    int l = ATI::l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    double x = wf.core_ocf[i];
    printf("%2i) %2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f   (%.2f)\n",
        i,n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*FPC::Hartree_invcm, eni*FPC::Hartree_eV,x);
  }

  //////////////////////////////////////////////////

  //Arrays to store results for outputting later:
  std::vector< std::vector< std::vector<float> > > AK; //float ok?

  int num_states = wf.nlist.size();
  AK.resize(desteps, std::vector< std::vector<float> >
    (num_states, std::vector<float>(qsteps)));

  std::vector<std::string> nklst;

  //pre-calculate the spherical Bessel function look-up table for efficiency
  std::vector< std::vector< std::vector<float> > > jLqr_f;
  AKF::sphericalBesselTable(jLqr_f,max_L,qmin,qmax,qsteps,wf.r);

  //Calculate the AK
  std::cout<<"\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)\n"
        ,demin/keV,demax/keV,demin,demax);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)\n"
        ,qmin/qMeV,qmax/qMeV,qmin,qmax);

  //Store state info (each orbital)
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
  #pragma omp parallel for if(desteps>3)
  for(int ide=0; ide<desteps; ide++){
    double y = double(ide)/(desteps-1);
    if(desteps==1) y=0;
    double dE = demin*pow(demax/demin,y);
    //Loop over core (bound) states:
    //std::vector< std::vector<float> > AK_nk;
    for(size_t is=0; is<wf.nlist.size(); is++){
      int k = wf.kappa[is];
      int l = ATI::l_k(k);
      if(l>max_l) continue;
      // if(plane_wave) AKF::calculateKpw_nk(wf,is,dE,jLqr_f[l],AK_nk);
      // else if(use_Zeff) AKF::calculateK_nk(wf,is,max_L,dE,jLqr_f,AK_nk,Zeff);
      // else AKF::calculateK_nk(wf,is,max_L,dE,jLqr_f,AK_nk);
      if(plane_wave) AKF::calculateKpw_nk(wf,is,dE,jLqr_f[l],AK[ide][is]);
      else if(use_Zeff)
        AKF::calculateK_nk(wf,is,max_L,dE,jLqr_f,AK[ide][is],Zeff);
      else AKF::calculateK_nk(wf,is,max_L,dE,jLqr_f,AK[ide][is]);
    }// END loop over bound states
    //AK[ide] = AK_nk;
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

  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.1f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.1f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  return 0;
}
