#include "AKF_akFunctions.h"

double fv(double v){
  if(v<0.1) return 0;
  else if(v>750.) return 0;
  else return 1./750.;
}



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
  AKF::akReadWrite(akfn,false,AKenq,nklst,qmin,qmax,demin,demax);

  int desteps = (int) AKenq.size();
  int num_states = (int) AKenq[0].size();
  int qsteps =  (int) AKenq[0][0].size();

  if(num_states != (int)nklst.size()) return 1;

  // Do q derivative on i grid:
  double dqonq = log(qmax/qmin)/(qsteps-1); //need to multiply by q

  std::cout<<qmin<<" "<<qmax<<" "<<demin<<" "<<demax<<"\n";
  std::cout<<qsteps<<" "<<num_states<<" "<<desteps<<"\n";

  double v=137*1.e-3; // typical v..integrate later...
  double m=1000; //XXX check units!!!
  double dE = 2000/27.2; //XXX check...later, function!
  double mv = 100.;//XXX check units!!!
  double qminus = m*v - sqrt(m*m*v*v-2*m*dE);
  double qplus  = m*v + sqrt(m*m*v*v-2*m*dE);

  double Aconst = 8*M_PI*FPC::c2;

  //find clostest dE to given "target"!
  //How to do q-grid integration?

  double a=0;
  for(int iq=0; iq<qsteps; iq++){
    double x = double(iq)/(qsteps-1);
    double q = qmin*pow(qmax/qmin,x);
    if(q<qminus || q>qplus) continue;
    a+=AKenq[1][1][iq]*pow(q,-3)*dqonq;
  }
  std::cout<<"a="<<a<<"\n";


  //Min/max m_chi
  //min/max m_v (or mu) = if 0, just 1 step
  // Also: for super massive case: linear! 1 step, absorb into x-section

  //Integrate over q,v
  //Clever way to do v integral (or constant v?)
  //nb: there is a v_min!


  return 0;
}
