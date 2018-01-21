#include "ElectronOrbitals.h"
//#include "physicalConstants.h"
#include "atomInfo.h"
//#include "adamsSolveLocalBS.h"
#include <iostream>


int main(void){

  std::cout<<"Hello world!\n";
  
  std::cout<<atinfo_a[15]<<"\n";

  ElectronOrbitals wf(1,1,1000);
  
  wf.localBoundState();

}
