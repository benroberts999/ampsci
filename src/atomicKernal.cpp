#include "ElectronOrbitals.h"
#include "INT_quadratureIntegration.h"
#include "PRM_parametricPotentials.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include "ContinuumOrbitals.h"
#include "HF_hartree.h"
#include "WIG_369j.h"
#include <gsl/gsl_sf_bessel.h>

//Declare angular coeficient. [See Phys.Rev.D 93, 115037 (2016).]
double CLkk(int L, int ka, int kb);

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

  // Max anglular momentums
  int max_l,max_lc,max_L;

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
    ifs >> max_l >> max_lc;       getline(ifs,jnk);
    ifs >> max_L;                 getline(ifs,jnk);
    ifs >> label;                 getline(ifs,jnk);
    ifs.close();
  }

  //default Hartree convergance goal:
  if(hart_del==0) hart_del=1.e-6;

  //allow for single-step in dE or q grid
  if(desteps==1) demax=demin;
  if(qsteps==1) qmax=qmin;

  // Fix maximum angular momentum values:
  if(max_lc<0) return 2;
  if(max_l<0 || max_l>3) max_l=3; //default: all core states (no >f)
  if(max_L>max_l+max_lc || max_L<0) max_L = max_l + max_lc; //triangle rule

  //alpha can't be zero:
  if(varalpha==0) varalpha=1.e-25;

  //Convert units for input q and dE range into atomic units
  double keV = (1.e3/HARTREE_EV);
  demin*=keV;
  demax*=keV;
  double qMeV = (1.e6/(HARTREE_EV*CLIGHT));
  qmin*=qMeV;
  qmax*=qMeV;


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

  //Arrays to store results for outputting later:
  std::vector< std::vector< std::vector<float> > > AK; //float ok?
  std::vector<float> qlst(qsteps);
  std::vector<float> dElst;
  std::vector<std::string> nklst;

  //pre-calculate the spherical Bessel function look-up table for efficiency
  std::cout<<std::endl;
  std::vector< std::vector< std::vector<float> > > jLqr_f;
  jLqr_f.resize(max_L+1, std::vector< std::vector<float> >
    (qsteps, std::vector<float>(wf.ngp)));
  for(int L=0; L<=max_L; L++){
    std::cout<<"\rCalculating spherical Bessel look-up table for L="
    <<L<<"/"<<max_L<<" .. "<<std::flush;
    #pragma omp parallel for
    for(int iq=0; iq<qsteps; iq++){
      double x=iq/(qsteps-1.);
      double q = qmin*pow(qmax/qmin,x);
      for(int ir=0; ir<wf.ngp; ir++){
        jLqr_f[L][iq][ir] = gsl_sf_bessel_jl(L, q*wf.r[ir]);
      }
    }
  }
  std::cout<<"done\n";

  //Calculate the AK
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

    //Loop over core (bound) states:
    std::vector< std::vector<float> > AK_nk;
    for(size_t is=0; is<wf.nlist.size(); is++){
      int k = wf.klist[is];
      int l = ATI_l_k(k);
      if(l>max_l) continue;
      int twoj = ATI_twoj_k(k);
      int n=wf.nlist[is];

      if(ide==0){
        std::string nk=std::to_string(n)+ATI_l(l)
          +"_{"+std::to_string(twoj)+"/2}";
        nklst.push_back(nk);
      }

      //Calculate continuum wavefunctions
      double ec = dE+wf.en[is];
      cntm.clear();
      if(ec>0)cntm.solveLocalContinuum(ec,max_lc);

      // Generate AK for each L, lc, and q
      //NB: L and lc summed, not stored indevidually
      std::vector<float> AK_nk_q(qsteps);
      for(int L=0; L<=max_L; L++){
        for(size_t ic=0; ic<cntm.klist.size(); ic++){
          int kc = cntm.klist[ic];
          int lc = ATI_l_k(kc);
          if(lc > max_lc) break;
          double dC_Lkk = CLkk(L,k,kc);
          if(dC_Lkk==0) continue;
          #pragma omp parallel for
          for(int iq=0; iq<qsteps; iq++){
            double x=iq/(qsteps-1.);
            double q = qmin*pow(qmax/qmin,x);
            double a = 0;
            double jLqr = 0;
            if(cntm.p.size()>0){
              if(ec<=0) std::cout<<"ERROR 244: !?!?\n";
              int maxj = wf.pinflist[is]; //don't bother going further
              //Do the radial integral:
              a=0;
              for(int j=0; j<maxj; j++){
                jLqr = jLqr_f[L][iq][j];
                a += (wf.p[is][j]*cntm.p[ic][j] + wf.q[is][j]*cntm.q[ic][j])
                     *jLqr*wf.drdt[j];// *h below!
              }
            }
            if(ide==0) qlst[iq]=q;
            AK_nk_q[iq] += dC_Lkk*pow(a*wf.h,2);
          }
        } // END loop over cntm states (ic)
      } // end L loop
      AK_nk.push_back(AK_nk_q);
      cntm.clear(); //deletes cntm wfs for this energy
    }// END loop over bound states
    dElst.push_back(dE);
    AK.push_back(AK_nk);
  }
  std::cout<<" Running dE step "<<desteps<<"/"<<desteps<<"  -  100% done  :)   "
  <<"                \n";// extra space to over-write any left-over junk.


  //ALSO: should write the AK array to a binary file here,
  //for ease of use in other applictaions.

  //Write out to text file (in gnuplot friendly form)
  std::ofstream ofile;
  std::string fname = "ak-"+Z_str+"_"+label+".txt";
  ofile.open(fname);
  ofile<<"dE(keV) q(MeV) ";
  for(size_t i=0; i<nklst.size(); i++) ofile<<nklst[i]<<" ";
  if(qsteps>1){
    ofile<<"\n";
    for(size_t i=0; i<dElst.size(); i++) ofile<<dElst[i]/keV<<" ";
  }
  ofile<<"\n\n";
  for(size_t i=0; i<AK.size(); i++){
    for(size_t k=0; k<AK[0][0].size(); k++){
      ofile<<dElst[i]/keV<<" "<<qlst[k]/qMeV<<" ";
      for(size_t j=0; j<AK[0].size(); j++){
        ofile<<AK[i][j][k]<<" ";
      }
      ofile<<"\n";
    }
    if(qsteps>1)ofile<<"\n";
  }
  ofile.close();

  // AK[i][j][k] // [dE][nk][q]
  //open file
  //write out 'header' info
  //write out array

  gettimeofday(&end, NULL);
  double total_time = (end.tv_sec-start.tv_sec)
  + (end.tv_usec - start.tv_usec)*1e-6;
  if(total_time<1) printf ("\nt=%.3f ms.\n",total_time*1000);
  else if(total_time<60) printf ("\nt=%.1f s.\n",total_time);
  else if(total_time<3600) printf ("\nt=%.1f mins.\n",total_time/60.);
  else printf ("\nt=%.1f hours.\n",total_time/3600.);

  return 0;
}



//******************************************************************************
double CLkk(int L, int ka, int kb)
/*
Calculates the angular coeficient (averaged over all m)
B. M. Roberts, V. A. Dzuba, V. V. Flambaum, M. Pospelov, and Y. V. Stadnik,
Phys. Rev. D 93, 115037 (2016). [arXiv:1604.04559]
*/
{
  int two_ja = ATI_twoj_k(ka);
  int two_jb = ATI_twoj_k(kb);
  double ja = 0.5*two_ja;
  double jb = 0.5*two_jb;
  int la = ATI_l_k(ka);
  int lb = ATI_l_k(kb);

  double tjB = WIG_3j(jb,L,ja,-0.5,0,0.5);
  if(fabs(tjB)==0) return 0;
  double B = 1./pow(tjB,2);

  //(-1)^(ja etc) -> calc sign
  int s1 = -1;
  if((two_ja+two_jb-2*(la+lb))%4==0) s1=1;

  double tj1 = WIG_3j(lb,la,L,0,0,0);
  double A = (1./4)*s1*(2*L+1)*pow(tj1,2);
  double X = s1*(two_ja+1)*(two_jb+1)*pow(tj1,2);
  double tj2 = WIG_3j(lb,la,L,-1,1,0);
  double Y = 8*sqrt(la*(la+1)*lb*(lb+1))*tj1*tj2;
  double Z = -4*(ka+1)*(kb+1)*pow(tj2,2);

  return (A*B)*(X+Y+Z);
}


//******************************************************************************
/*
Helper functions for the binary read in/out
*/
template<typename T>
int binary_write(std::fstream& stream, const T& value){
    stream.write(reinterpret_cast<const char*>(&value), sizeof(T));
    return 0;
}
template<typename T>
int binary_read(std::fstream& stream, T& value){
    stream.read(reinterpret_cast<char*>(&value), sizeof(T));
    return 0;
}
int binary_str_write(std::fstream& stream, std::string& value){
    binary_write(stream,value.length());
    stream.write(value.c_str(), value.length());
    return 0;
}
int binary_str_read(std::fstream& stream, std::string& ovalue){
    size_t temp_len;
    binary_read(stream,temp_len);
    char* value = new char[temp_len+1];
    stream.read(value, temp_len);
    value[temp_len] = '\0'; //null 'end of string' character
    ovalue = value;
    delete [] value;
    return 0;
}
template<typename T>
int binary_rw(std::fstream& stream, T& value, bool write){
    if(write) binary_write(stream, value);
    else binary_read(stream, value);
    return 0;
}
template<typename T>
int binary_str_rw(std::fstream& stream, T& value, bool write){
    if(write) binary_str_write(stream,value);
    else binary_str_read(stream,value);
    return 0;
}
//XXX remove first  ones???


//************
int akReadWrite(std::string fname, bool write,
  std::vector< std::vector< std::vector<double> > > &AK,
  std::vector<double> &dElst,
  std::vector<double> &nklst,
  std::vector<double> &qlst)
{
  std::fstream iof;
  if(write) iof.open(fname,std::ios_base::out|std::ios_base::binary);
  else      iof.open(fname,std::ios_base::in|std::ios_base::binary);

  if(write){
    binary_rw(iof,AK.size(),write);
    binary_rw(iof,AK[0].size(),write);
    binary_rw(iof,AK[0][0].size(),write);
  }else{
    int nq,ns,nde;
    binary_rw(iof,nde,write);
    binary_rw(iof,ns,write);
    binary_rw(iof,nq,write);
    AK.resize(nde,std::vector< std::vector<double> >(ns,std::vector<double>(nq)));
    dElst.resize(nde);
    nklst.resize(ns);
    qlst.resize(nq);
  }
  for(size_t ie=0; ie<AK.size(); ie++){
    binary_rw(iof,dElst[ie],write);
    for(size_t in=0; ie<AK[0].size(); ie++){
      if(ie==0) binary_rw(iof,nklst[in],write);
      for(size_t iq=0; ie<AK[0][0].size(); ie++){
        if(ie==0 && in==0) binary_rw(iof,qlst[iq],write);
        binary_rw(iof,AK[ie][in][iq],write);
      }
    }
  }


  return 0; //?
}
