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


  double m=1.; //XXX in GeV
  //convert WIMP masses from GeV to a.u.
  m *= (1.e3/FPC::m_e_MeV);

  bool finite_med = true;
  //Vector mass:
  //NB: have both '0' case (easy), and 'inf' case! (changes coupling const!)
  double mv = 10.;
  //convert from MeV to au:
  mv *= (1./FPC::m_e_MeV);


  double v=137*1.e-2; // typical v..integrate later...

  //double dE = demin;



  //Numberical constant. Note: Inlcudes dqonq
  //NOTE: if alpgha_chi = alpha in the code.. kill c2 !!
  double Aconst = 8*M_PI*FPC::c2*dqonq;
  std::cout<<"\nAconst="<<Aconst<<"\n";
  std::cout<<"\n\n demin="<<demin*FPC::Hartree_eV<<"\n";
  std::cout<<"demax="<<demax*FPC::Hartree_eV<<"\n\n";

  //find clostest dE to given "target"!
  //How to do q-grid integration?

  //XXX XXX XXX the 3s state showing zero... WHY?!?!?!

  //XXX define a v grid. Sort f(v) into Array.
  // This a) speeds it up
  // b) easier function-swapping
  // c) convert to atomic units

for(int ie=0; ie<desteps; ie++){
  double a=0;
  double xe = double(ie)/(desteps-1); //XXX allow for single dE !!
  double dE = demin*pow(demax/demin,xe);
  double arg = pow(m*v,2)-2.*m*dE; //may be negative; skip!
  if(arg<0) continue;
  double qminus = m*v - sqrt(arg);
  double qplus  = m*v + sqrt(arg);
  std::cout<<ie<<" "<<qminus<<" "<<qplus<<"\n";
  for(int ink=0; ink<num_states; ink++){
    double a_nk=0;
    for(int iq=0; iq<qsteps; iq++){
      double x = double(iq)/(qsteps-1);
      double q = qmin*pow(qmax/qmin,x);
      if(q<qminus || q>qplus) continue;
      double dq_on_dqonq = q; //devide by dqonq - just a const. Include in Acont!
      double Fq = q*dq_on_dqonq;
      if(finite_med) Fq /= pow(q*q+mv*mv,2);
      a_nk+=(1./v)*fv(v)*Fq*AKenq[ie][ink][iq];// xxx choose state + dE!
    }
    std::cout<<nklst[ink]<<"; a="<<a_nk*Aconst<<"\n";
    a+=a_nk;
  }
  std::cout<<" "<<ie<<"; a="<<a*Aconst<<"\n";
}

std::cout<<AKenq[0][2][100]<<"\n";


  //Min/max m_chi
  //min/max m_v (or mu) = if 0, just 1 step
  // Also: for super massive case: linear! 1 step, absorb into x-section

  //Integrate over q,v
  //Clever way to do v integral (or constant v?)
  //nb: there is a v_min!

  //XXX also: integrate over dE (in bins?)


  return 0;
}
