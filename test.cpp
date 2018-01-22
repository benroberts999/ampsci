#include "ElectronOrbitals.h"
#include "atomInfo.h"
#include <iostream>


int main(void){

  clock_t t;
  t = clock();

  printf("\nRunning SolveDBS for Local potential\n");
  printf("*******************************\n");
  //printf("NGP=%i, Size of box: Rmax=%.1f a.u., h=%f\n\n",NGP,r(NGP-1),h);

  int Z=1;
  int A=1;
  ElectronOrbitals wf(Z,A,1000,1);

  wf.localBoundState(3);

  int num_states = wf.nlist.size();
  for(int i=0; i<num_states; i++){
    int n=wf.nlist[i];
    int k=wf.klist[i];
    int twoj = 2*abs(k)-1;
    int l = (abs(2*k+1)-1)/2;
    double del = wf.en[i] - wf.diracen(wf.Z,n,k);
    //std::cout<<wf.diracen(wf.Z,n,k)<<"\n";
    printf("%i %2i (%s %i/2) En=%.15f   %i %7.1e  %8.1e\n",
        n,k,atinfo_l(l).c_str(),twoj,wf.en[i],wf.itslist[i],wf.epslist[i],del);
  }


}
