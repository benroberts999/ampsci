#include "ElectronOrbitals.h"
#include "ATI_atomInfo.h"
#include "INT_quadratureIntegration.h"
#include "ChronoTimer.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "HartreeFockClass.h"


int main(void){

  ChronoTimer sw(true);//start timer

  int Z,A;
  int n_max,l_max;
  int ngp = 2000;
  double r0,rmax;
  double varalpha;
  int iextra;
  std::ifstream ifile;
  ifile.open("h-like.in");
  {
    std::string junk;
    ifile >> Z >> A;            getline(ifile,junk);
    ifile >> n_max >> l_max;    getline(ifile,junk);
    ifile >> r0 >> rmax >> ngp; getline(ifile,junk);
    ifile >> varalpha;          getline(ifile,junk);
    ifile >> iextra;            getline(ifile,junk);
  }
  ifile.close();
  // bool extra=false;
  // if(iextra==1) extra=true;

  printf("\nRunning SolveDBS for Local H-like potential, Z=%i\n",Z);
  printf("*************************************************\n");

  //Generate the orbitals object:
  Z=54;
  ElectronOrbitals wf(Z,A,ngp,r0,rmax,varalpha);
  if(A!=0) wf.sphericalNucleus();

  wf.determineCore("Xe");
  HartreeFock hf(wf);

  printf(" n l_j    k  R_inf its eps     En (au)            Error (au)\n");
  int num_states = wf.nlist.size();
  for(int i=0; i<num_states; i++){
    int n=wf.nlist[i];
    int k=wf.kappa[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double del = wf.en[i] - wf.diracen(wf.Z,n,k);
    double rinf = wf.r[wf.pinflist[i]];
    printf("%2i %s_%i/2 (%2i)  %3.0f %3i  %5.0e  %.15f  %7.0e\n",
        n,ATI::l_symbol(l).c_str(),twoj,k,rinf,wf.itslist[i],wf.epslist[i],
        wf.en[i],del);
  }



  std::cout<<"\n Total time: "<<sw.reading_str()<<"\n";

  return 0;
}
