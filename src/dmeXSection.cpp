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
  std::cout<<"\n\n dqonq="<<dqonq<<"\n\n";

  std::cout<<qmin<<" "<<qmax<<" "<<demin<<" "<<demax<<"\n";
  std::cout<<qsteps<<" "<<num_states<<" "<<desteps<<"\n";


  double m=1; //XXX in GeV
  //convert WIMP masses from GeV to a.u.
  m *= (1.e3/FPC::m_e_MeV);

  bool finite_med = true;
  //Vector mass:
  //NB: have both '0' case (easy), and 'inf' case! (changes coupling const!)
  double mv = 0.1;
  //convert from MeV to au:
  mv *= (1./FPC::m_e_MeV);


  double v=137*1.e-3; // typical v..integrate later...

  double dE = demin;

  double arg = pow(m*v,2)-2.*m*dE; //may be negative; skip!
  double qminus = m*v - sqrt(arg);
  double qplus  = m*v + sqrt(arg);

  //Numberical constant. Note: Inlcudes dqonq
  //NOTE: if alpgha_chi = alpha in the code.. kill c2 !!
  double Aconst = 8*M_PI*FPC::c2*dqonq;
  std::cout<<"\n\n Aconst="<<Aconst<<"\n\n";

  std::cout<<"\n\n demax="<<demax*FPC::Hartree_eV<<"\n\n";

  //find clostest dE to given "target"!
  //How to do q-grid integration?

  double a=0;
  for(int iq=0; iq<qsteps; iq++){
    double x = double(iq)/(qsteps-1);
    double q = qmin*pow(qmax/qmin,x);
    double dq_on_dqonq = q; //devide by dqonq - just a const. Include in Acont!
    double Fq = q*dq_on_dqonq;
    if(finite_med) Fq /= pow(q*q+mv*mv,2);
    if(q<qminus || q>qplus) continue;
    a+=(1./v)*Fq*AKenq[1][2][iq];// xxx choose state + dE!
  }
  std::cout<<"a="<<a*Aconst<<"\n";


  //Min/max m_chi
  //min/max m_v (or mu) = if 0, just 1 step
  // Also: for super massive case: linear! 1 step, absorb into x-section

  //Integrate over q,v
  //Clever way to do v integral (or constant v?)
  //nb: there is a v_min!

  //XXX also: integrate over dE (in bins?)


  return 0;
}
