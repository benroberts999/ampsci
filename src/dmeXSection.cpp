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

  int nqsteps = qlst.size();
  double qmin = qlst[0];
  double qmax = qlst[nqsteps-1];

  double dqonq = log(qmax/qmin)/(nqsteps-1); //need to multiply by q


  return 0;
}
