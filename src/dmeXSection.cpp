#include "akFunctions.h"

//******************************************************************************
int main(void){


  //Arrays to store results for outputting later:
  std::vector< std::vector< std::vector<float> > > AK; //float ok?
  std::vector<float> qlst;
  std::vector<float> dElst;
  std::vector<std::string> nklst;

  //Read in AK file
  std::string akfn = "ak-Xe_test.bin"; //XXX
  akReadWrite(akfn,false,AK,dElst,nklst,qlst);


  return 0;
}
