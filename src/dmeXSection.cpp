#include "AKF_akFunctions.h"
#include "SHM_standardHaloModel.h"

// double fv(double v){
//   if(v<0.1) return 0;
//   else if(v>750.) return 0;
//   else return 1./750.;
// }

// double fv(double v, double phi=0);
double fv_au(double v_au, double phi=0);

//******************************************************************************
double fv_au(double v_au, double phi){
  double v = v_au * (FPC::c_SI/FPC::c); //will be in m/s
  v/=1.e3; //convert from m/s -> km/s
  return SHM::fv(v,phi);
}


//******************************************************************************
int main(void){

  //Units Conversion Factors.
  //Conver FROM a.u. TO rel./other units
  double E_to_keV = FPC::Hartree_eV/1000.; //au -> keV (energy)
  double M_to_GeV = FPC::m_e_MeV/1000.;
  double M_to_MeV = FPC::m_e_MeV;
  double V_to_kms = (FPC::c_SI/FPC::c)/1000.;
  double Q_to_MeV = FPC::Hartree_eV*FPC::c/1.e6;

  //define input parameters
  std::string akfn; //name of K file to read in
  int vsteps;
  double mmin,mmax,mvmin,mvmax; //m_chi and m_v masses
  int n_m,i_mv,n_mv;
  double de_target;
  std::string label="testx";

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("dmeXSection.in");
    std::string jnk;
    ifs >> akfn;                    getline(ifs,jnk);
    ifs >> vsteps;                  getline(ifs,jnk);
    ifs >> mmin >> mmax >> n_m;     getline(ifs,jnk);
    ifs >> i_mv;                    getline(ifs,jnk);
    ifs >> mvmin >> mvmax >> n_mv;  getline(ifs,jnk);
    ifs >> de_target;               getline(ifs,jnk);
    ifs >> label;                   getline(ifs,jnk);
    ifs.close();
  }

  //DM mass:
  mmin /= M_to_GeV;
  mmax /= M_to_GeV;
  if(mmax<=mmin || n_m<=1){
    mmax = mmin;
    n_m = 1;
  }
  //Mediator mass: convert units + set-up finite/infinite case
  bool finite_med = true;
  if(i_mv==0){
    //massless case. Do sepperately
    mvmin = mvmax = 0;
    n_mv = 1;
  }else if(i_mv==2){
    //Heavy-mediator case (contact interaction)
    finite_med = false;
    mvmin = mvmax = 1./0.;
    n_mv = 1;
  }else if(i_mv==1){
    mvmin /= M_to_MeV;
    mvmax /= M_to_MeV;
  }else{
    std::cout<<"Wrong m_v input given. Try again. Or not, whatevs\n";
    return 1;
  }
  if(mvmax<=mvmin || n_mv<=1){
    mvmax = mvmin;
    n_mv = 1;
  }

  //convert E target to au
  std::cout<<de_target<<"\n";
  de_target /= E_to_keV;
  std::cout<<de_target<<"\n";




  //DM mass (in GeV)
  double m=1.;
  //convert WIMP masses from GeV to a.u.
  m /= M_to_GeV;

  double mv = mvmin;

  //Arrays/values to be filled from input AK file:
  std::vector< std::vector< std::vector<float> > > AKenq;
  std::vector<std::string> nklst;
  double qmin,qmax,demin,demax;

  //Read in AK file
  std::cout<<"Opening file: "<<akfn<<".bin\n";
  AKF::akReadWrite(akfn,false,AKenq,nklst,qmin,qmax,demin,demax);

  int desteps    = (int) AKenq.size();
  int num_states = (int) AKenq[0].size();
  int qsteps     = (int) AKenq[0][0].size();

  if(num_states != (int)nklst.size()) return 1; //just sanity check


  bool plotv = false;

  //If only doing a single dE:
  int i_et = -1;
  if(de_target>0 || desteps==1){
    plotv = true;
    double tmp = (desteps-1)*log(de_target/demin)/log(demax/demin);
    i_et = (int) ceil(tmp);
    double xe = double(i_et)/(desteps-1);
    de_target = demin*pow(demax/demin,xe);
    if(desteps==1){
      de_target = demin;
      i_et = 0;
    }
    demin = de_target;
    demin = de_target;
    std::cout<<de_target<<"\n";
    if(i_et >= desteps || i_et<0) std::cout<<"Bad E target\n";
    desteps = 1;
  }


  // Do q derivative on i grid:
  double dqonq = log(qmax/qmin)/(qsteps-1); //need to multiply by q

  //converts from MeV to au
  // double qMeV = 1./Q_to_MeV;

  //Grid of f_v(v). Can use to change vel profiles
  //int vsteps = 1000;
  std::vector<double> arr_fv(vsteps);
  double max_v = (SHM::max_v)/V_to_kms;
  double dv = max_v/vsteps;
  for(int i=0; i<vsteps; i++){
    double v = (i+1)*dv; //don't include zero
    arr_fv[i] = fv_au(v);
  }


  printf("q  grid: %6.2f -> %6.2f  MeV, %5i steps\n"
    ,qmin*Q_to_MeV,qmax*Q_to_MeV,qsteps);
  printf("E  grid: %6.2f -> %6.2f  keV, %5i steps\n"
    ,demin*E_to_keV,demax*E_to_keV,desteps);
  printf("v  grid: %6.2f -> %6.1f km/s, %5i steps\n"
    ,dv*V_to_kms,max_v*V_to_kms,vsteps);

  printf("Mx grid: %6.2f -> %6.1f  GeV, %5i steps\n"
    ,mmin*M_to_GeV,mmax*M_to_GeV,n_m);
  printf("Mv grid: %6.3f -> %6.3f  MeV, %5i steps\n"
    ,mvmin*M_to_MeV,mvmax*M_to_MeV,n_mv);


  printf("m_chi = %6.1f GeV\n",m*M_to_GeV);
  printf("m_v=mu= %6.3f MeV\n",mv*M_to_MeV);

  //Numberical constant. Note: Inlcudes dqonq
  //NOTE: if alpgha_chi = alpha in the code.. kill c2 !!
  //NB: includes dv - Not good if not integrating over v!?!
  double A = 8*M_PI*FPC::c2;

  std::cout<<"A = "<<A<<"\n";

  std::ofstream ofv;
  ofv.open("vplot-"+label+".txt");


  for(int ie=0; ie<desteps; ie++){
    //double a=0;
    double xe = double(ie)/(desteps-1);
    double dE = demin*pow(demax/demin,xe);
    if(desteps==1){
      dE = de_target;
      ie = i_et;
    }
    double vmin = sqrt(dE*2/m);
    double ds_de=0;
    for(int iv=0; iv<vsteps; iv++){
      double v = (iv+1)*dv;
      if(v<vmin) continue;
      double fvonv = arr_fv[iv]/v;
      //if(plotv) fvonv = 1./v; //don't mult by fv if plotting sig(v)
      double arg = pow(m*v,2)-2.*m*dE; //may be negative; skip!
      if(arg<0 || fvonv==0) continue;
      double qminus = m*v - sqrt(arg);
      double qplus  = m*v + sqrt(arg);
      double ds_de_dv=0;
      for(int ink=0; ink<num_states; ink++){
        for(int iq=0; iq<qsteps; iq++){
          double x = double(iq)/(qsteps-1);
          double q = qmin*pow(qmax/qmin,x);
          if(q<qminus || q>qplus) continue;
          double dq_on_dqonq = q; //devide by dqonq - just a const.
          //Include dqonq in Aconst!
          double Fq = q*dq_on_dqonq; //extra q factor from dqonq (Jacobian)
          if(finite_med) Fq /= pow(q*q+mv*mv,2);
          ds_de_dv += Fq*AKenq[ie][ink][iq];  //dq/q included in aAconst
        }
      }
      ds_de += ds_de_dv*fvonv;
      if(plotv) ofv<<v*V_to_kms<<" "<<(ds_de_dv/(v*v))*A*dqonq<<"\n";
      if(plotv) std::cout<<v*V_to_kms<<" "<<(ds_de_dv/(v*v))*A*dqonq<<"\n";
    }
     std::cout<<"E="<<dE*E_to_keV<<" keV"
     <<"; <sig.v>/dE = "<<ds_de*A*dv*dqonq<<"\n";
  }


  return 0;
}
