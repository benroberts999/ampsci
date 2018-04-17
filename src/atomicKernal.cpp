#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include "ContinuumOrbitals.h"
#include "HF_hartree.h"

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

  std::vector<std::string> str_core; //States for the core

  //Green potential parameters
  int Gf;
  double Gh,Gd;

  //q and dE grids:
  double qmin,qmax,demin,demax;
  int qsteps,desteps;

  std::string label; //label for output file

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("atomicKernal.in");
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
    ifs >> Gf >> Gh >> Gd;        getline(ifs,jnk);
    ifs >> hart_del;              getline(ifs,jnk);
    ifs >> varalpha;              getline(ifs,jnk);
    ifs >>demin>>demax>>desteps;  getline(ifs,jnk);
    ifs >>qmin>> qmax >> qsteps;  getline(ifs,jnk);
    ifs >> label;                 getline(ifs,jnk);
    ifs.close();
  }

  //default Hartree convergance goal:
  if(hart_del==0) hart_del=1.e-6;

  //allow for single-step in dE or q grid
  if(desteps==1) demax=demin;
  if(qsteps==1) qmax=qmin;

  //alpha can't be zero:
  if(varalpha==0) varalpha=1.e-25;

  //Convert units for input dE into atomic units
  double keV = (1.e3/HARTREE_EV);
  demin*=keV;
  demax*=keV;
  double qMeV = (1.e6/(HARTREE_EV*CLIGHT));
  //nb: I don't change the units for q (momentum transfer) for now

  //Look-up atomic number, Z, and also A
  int Z = ATI_get_z(Z_str);
  if(Z==0) return 2;
  if(A==-1) A=ATI_a[Z]; //if none given, get default A

  printf("\nRunning Atomic Kernal for %s, Z=%i A=%i\n",
    Z_str.c_str(),Z,A);
  printf("*************************************************\n");
  if(Gf!=0) printf("Using Green potential: H=%.4f  d=%.4f\n",Gh,Gd);
  else printf("Using Hartree potential (converge to %.0e)\n",hart_del);

  //Generate the orbitals object:
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  printf("Grid: pts=%i h=%6.4f r0=%.0e Rmax=%5.1f\n"
  , wf.ngp,wf.h,wf.r[0],wf.r[wf.ngp-1]);

  // Check if 'h' is small enough for oscillating region:
  double h_target = (M_PI/15)/sqrt(2.*demax);
  if(wf.h>2*h_target){
    std::cout<<"\nWARNING 101: Grid not dense enough for contimuum state with "
      <<"ec="<<demax<<" (h="<<wf.h<<", need h<"<<h_target<<")\n";
    std::cout<<"Program will continue, but continuum wfs may be bad.\n\n";
  }

  //If non-zero A is given, use spherical nucleus.
  if(A>0) wf.sphericalNucleus();

  //Determine which states are in the core:
  int core_ok = wf.determineCore(str_core);
  if(core_ok==2){
    std::cout<<"Problem with core: ";
    for(size_t i=0; i<str_core.size(); i++) std::cout<<str_core[i]<<" ";
    std::cout<<"\n";
    return 1;
  }

  if(Gf==0){
    HF_hartreeCore(wf,hart_del);
  }else{
    //Fill the electron part of the (local/direct) potential
    for(int i=0; i<wf.ngp; i++){
      wf.vdir.push_back(PRM_green(Z,wf.r[i],Gh,Gd));
    }
    // Solve Dirac equation for each (bound) core state:
    wf.solveInitialCore();
  }

  //make list of energy indices in sorted order:
  std::vector<int> sort_list;
  wf.sortedEnergyList(sort_list);

  //Output results:
  printf("\n n l_j    k Rinf its    eps     En (au)     En (/cm)    En (eV)\n");
  for(size_t m=0; m<sort_list.size(); m++){
    int i = sort_list[m];
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = ATI_twoj_k(k);
    int l = ATI_l_k(k);
    double rinf = wf.r[wf.pinflist[i]];
    double eni = wf.en[i];
    printf("%2i %s_%i/2 %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f\n",
        n,ATI_l(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        eni, eni*HARTREE_ICM, eni*HARTREE_EV);
  }

  //Continuum wavefunction object
  ContinuumOrbitals cntm(wf);

  //XXX
  /*
  ALSO:
    * Need angular factor
    * Different (l', L) states
    * -- check <||> is non-zero before calc!
  */
  std::vector< std::vector< std::vector<float> > > AK; //float ok?
  std::vector<float> qlst(qsteps);
  std::vector<float> dElst;
  std::vector<std::string> nklst;

  int ic=0; //cntm state!
  int max_l=0; //maximum bound-state l

  std::cout<<"\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)\n"
  ,   demin/keV,demax/keV,demin,demax);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)\n"
  ,   qmin/qMeV,qmax/qMeV,qmin,qmax);
  for(int ide=0; ide<desteps; ide++){
    int pc = int(100.*ide/desteps);
    std::cout<<" Running dE step "<<ide<<"/"<<desteps<<"  -  "<<pc<<"% done"
    <<"                                        \r";
    std::cout.flush();
    double y;
    if(desteps>1) y=ide/(desteps-1.);
    else y=0;
    double dE = demin*pow(demax/demin,y);

    std::vector< std::vector<float> > AK_nk;
    for(size_t is=0; is<wf.nlist.size(); is++){
      int k=wf.klist[is];
      int l = ATI_l_k(k);
      if(l>max_l) continue;
      int twoj = ATI_twoj_k(k);
      int n=wf.nlist[is];

      if(ide==0){
        std::string nk=std::to_string(n)+ATI_l(l)
          +"_{"+std::to_string(twoj)+"/2}";
        nklst.push_back(nk);
      }

      double ec = dE+wf.en[is];
      if(ec>0)cntm.solveLocalContinuum(ec,0);

      std::vector<float> AK_nk_q(qsteps);
      #pragma omp parallel for
      for(int iq=0; iq<qsteps; iq++){
        double x=iq/(qsteps-1.);
        double q = qmin*pow(qmax/qmin,x);

        double a=0;
        double jLqr;
        if(cntm.p.size()>0){
          int maxj = wf.pinflist[is]; //don't bother going further
          // int maxj = wf.ngp;
          for(int j=0; j<maxj; j++){
            jLqr = sin(q*wf.r[j]) / (q*wf.r[j]);
            a += (wf.p[is][j]*cntm.p[ic][j] + wf.q[is][j]*cntm.q[ic][j])
                 *jLqr*wf.drdt[j];// *h below!
          }
        }
        if(ide==0) qlst[iq]=q;
        double C_Lkk = 2.; //XXX only for all j=1/2 states!
        AK_nk_q[iq] = C_Lkk*pow(a*wf.h,2);
      }
      AK_nk.push_back(AK_nk_q);
    }
    dElst.push_back(dE);
    AK.push_back(AK_nk);
    cntm.clear(); //deletes cntm wfs for this energy
  }
  std::cout<<" Running dE step "<<desteps<<"/"<<desteps<<"  -  100% done  :)   "
  <<"                \n";// extra space to over-write any left-over junk.


  //ALSO: should write the AK array to a binary file here,
  //for ease of use in other applictaions.

  //Write out file (in gnuplot friendly form)
  std::ofstream ofile;
  std::string fname = "ak-test_"+label+".txt";
  ofile.open(fname);
  ofile<<"dE(au) q(au) ";
  for(size_t i=0; i<nklst.size(); i++) ofile<<nklst[i]<<" ";
  ofile<<"\n";
  for(size_t i=0; i<dElst.size(); i++) ofile<<dElst[i]<<" ";
  ofile<<"\n\n";
  for(size_t i=0; i<AK.size(); i++){
    for(size_t k=0; k<AK[0][0].size(); k++){
      ofile<<dElst[i]<<" "<<qlst[k]<<" ";
      for(size_t j=0; j<AK[0].size(); j++){
        ofile<<AK[i][j][k]<<" ";
      }
      ofile<<"\n";
    }
    ofile<<"\n";
  }
  ofile.close();

  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.1f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.1f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  return 0;
}
