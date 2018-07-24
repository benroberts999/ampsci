#include "akFunctions.h"

//******************************************************************************
int main(void){


  std::string akfn;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("dmeXSection.in");
    std::string jnk;
    ifs >> akfn;         getline(ifs,jnk);
    ifs.close();
  }


  //Arrays to store results for outputting later:
  std::vector< std::vector< std::vector<float> > > AKenq;
  std::vector<std::string> nklst;
  double qmin,qmax,demin,demax;

  //Read in AK file
  //std::string akfn = "ak-Xe_L6"; //XXX
  std::cout<<"Opening file: "<<akfn<<".bin\n";
  akReadWrite(akfn,false,AKenq,nklst,qmin,qmax,demin,demax);

  int desteps = (int) AKenq.size();
  int num_states = (int) AKenq[0].size();
  int qsteps =  (int) AKenq[0][0].size();

  if(num_states != (int)nklst.size()) return 1;

  double dqonq = log(qmax/qmin)/(qsteps-1); //need to multiply by q
  //XXX check!

  std::cout<<qmin<<" "<<qmax<<" "<<demin<<" "<<demax<<"\n";
  std::cout<<qsteps<<" "<<num_states<<" "<<desteps<<"\n";


  //Min/max m_chi
  //min/max m_v (or mu) = if 0, just 1 step
  // Also: for super massive case: linear! 1 step, absorb into x-section

  //Integrate over q,v
  //Clever way to do v integral (or constant v?)
  //nb: there is a v_min!


  return 0;
}
