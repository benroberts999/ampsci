#include "AKF_akFunctions.h"
#include "SHM_standardHaloModel.h"
#include <iomanip>

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
template<typename T>
double dsdE_iEdEvum_qg(std::vector< std::vector< std::vector<T> > > &K_enq,
  int iE, double dE, double v, double mu, double mx,
  double qmin, double qmax)
{
  double A_au = 8*M_PI*FPC::c2;
  double v2 = pow(v,2);

  int desteps    = (int) K_enq.size();
  int num_states = (int) K_enq[0].size();
  int qsteps     = (int) K_enq[0][0].size();
  if(iE>=desteps) return 0; //or -1?

  bool finite_med = true;
  if(mu<0) finite_med = false;

  // Do q derivative on i grid:
  double dqonq = log(qmax/qmin)/(qsteps-1); //need to multiply by q for dq

  double arg = pow(mx*v,2)-2.*mx*dE;
  if(arg<0) return 0;

  double qminus = mx*v - sqrt(arg);
  double qplus  = mx*v + sqrt(arg);
  if(qminus>qmax || qplus<qmin) return 0; //or -1 ?
  double dsdE = 0;
  for(int ink=0; ink<num_states; ink++){
    #pragma omp parallel for
    for(int iq=0; iq<qsteps; iq++){
      double x = double(iq)/(qsteps-1);
      double q = qmin*pow(qmax/qmin,x);
      if(q<qminus || q>qplus) continue;
      double t_Fq = q*q; //(dq/q) is constant, multiply at end
      if(finite_med) t_Fq /= pow(q*q+mu*mu,2);
      #pragma omp critical
      {
        dsdE += t_Fq*K_enq[iE][ink][iq];  //dq/q included in aAconst
      }
    }
  }

  dsdE *= (A_au/v2)*dqonq;
  return dsdE;
}

//******************************************************************************
int main(void){

  //Units Conversion Factors.
  //Conver FROM a.u. TO rel./other units
  double E_to_keV = FPC::Hartree_eV/1000.; //au -> keV (energy)
  double M_to_GeV = FPC::m_e_MeV/1000.;
  double M_to_MeV = FPC::m_e_MeV;
  double V_to_kms = (FPC::c_SI/FPC::c)/1000.;
    double V_to_cms = (FPC::c_SI*FPC::alpha)*100.; // au -> cm/s
    double V_to_cmday = V_to_cms*(24*60*60); // cm/s -> cm/day
  double Q_to_MeV = FPC::Hartree_eV*FPC::c/1.e6;
  //cross section: ds/dE

  double dsdE_to_cm2keV = pow(FPC::aB_cm,2)/E_to_keV; // au -> cm^2/keV
  double dsvdE_to_cm3keVday = dsdE_to_cm2keV*V_to_cmday;

  //Numberical constant. For alpha_chi = 1 !
  // (if want alpgha_chi = alpha in the code, kill c2 )
  double A_au = 8*M_PI*FPC::c2;
  //Convert from au -> cm^2/keV
  double A_cmkeV = A_au*dsdE_to_cm2keV;

  //DM energy density (in GeV/cm^3)
  //note: will be divided by m_chi, also in GeV
  double rhoDM_GeVcm3 = 0.4; // GeV/cm^3

  //define input parameters
  std::string akfn; //name of K file to read in
  int vsteps;
  double mxmin,mxmax,mvmin,mvmax; //m_chi and m_v masses
  int n_mx,i_mv,n_mv;
  double de_target;
  std::string label="testx";

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("dmeXSection.in");
    std::string jnk;
    ifs >> akfn;                    getline(ifs,jnk);
    ifs >> vsteps;                  getline(ifs,jnk);
    ifs >> mxmin >> mxmax >> n_mx;   getline(ifs,jnk);
    ifs >> i_mv;                    getline(ifs,jnk);
    ifs >> mvmin >> mvmax >> n_mv;  getline(ifs,jnk);
    ifs >> de_target;               getline(ifs,jnk);
    ifs >> label;                   getline(ifs,jnk);
    ifs.close();
  }
  //DM mass:
  mxmin /= M_to_GeV;
  mxmax /= M_to_GeV;
  if(mxmax<=mxmin || n_mx<=1){
    mxmax = mxmin;
    n_mx = 1;
  }
  //Mediator mass: convert units + set-up finite/infinite case
  if(mvmax<=mvmin || n_mv<=1){
    mvmax = mvmin;
    n_mv = 1;
  }
  if(i_mv==0){
    //massless case. Do sepperately
    mvmin = mvmax = 0;
    n_mv = 1;
  }else if(i_mv==2){
    //Heavy-mediator case (contact interaction)
    mvmin = mvmax = -1; //1./0.;
    n_mv = 1;
  }else if(i_mv==1){
    mvmin /= M_to_MeV;
    mvmax /= M_to_MeV;
  }else{
    std::cout<<"Wrong m_v input given. Try again. Or not, whatevs\n";
    return 1;
  }

  //convert E target to au
  de_target /= E_to_keV;

  // //masses (loop over later!) XXX
  // double mx = mxmin;
  // double mv = mvmin;

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

  //If only doing a single dE:
  bool plotv = false;
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
  double dqonq = log(qmax/qmin)/(qsteps-1); //need to multiply by q for dq

  //Grid of f_v(v). Can use to change vel profiles
  std::vector<double> arr_fv(vsteps);
  double max_v = (SHM::max_v)/V_to_kms;
  double dv = max_v/vsteps;
  for(int i=0; i<vsteps; i++){
    double v = (i+1)*dv; //don't include zero
    arr_fv[i] = fv_au(v);
  }


  //Print the grid info to screen:
  printf("q  grid: %6.2f -> %6.2f  MeV, %5i steps\n"
    ,qmin*Q_to_MeV,qmax*Q_to_MeV,qsteps);
  printf("E  grid: %6.2f -> %6.2f  keV, %5i steps\n"
    ,demin*E_to_keV,demax*E_to_keV,desteps);
  printf("v  grid: %6.2f -> %6.1f km/s, %5i steps\n"
    ,dv*V_to_kms,max_v*V_to_kms,vsteps);
  printf("Mx grid: %6.2f -> %6.1f  GeV, %5i steps\n"
    ,mxmin*M_to_GeV,mxmax*M_to_GeV,n_mx);
  printf("Mv grid: %6.3f -> %6.3f  MeV, %5i steps\n"
    ,mvmin*M_to_MeV,mvmax*M_to_MeV,n_mv);
  printf("A = %.1f au = %.2e cm^2/keV\n",A_au,A_cmkeV);
  std::cout<<"ds/dE conversion factor:   "<<dsdE_to_cm2keV<<" cm^2/keV\n"
    <<"ds.v/dE conversion factor: "
    <<dsvdE_to_cm3keVday<<"   cm^3/keV/day\n";

  // return 1;

  std::ofstream ofv,ofE;
  ofv.open("vplot-"+label+".txt");
  ofE.open("Eplot-"+label+".txt");

  double dE = de_target;
  int ie = i_et;

  std::vector< std::vector< std::vector<float> > > ds_mv_mx_x;

  if(plotv) ds_mv_mx_x.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,std::vector<float>(vsteps)));
  else ds_mv_mx_x.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,std::vector<float>(desteps)));


  for(int imv=0; imv<n_mv; imv++){
    double xmv = double(imv)/(n_mv-1);
    double mv = mvmin*pow(mvmax/mvmin,xmv);
    if(n_mv==1) mv = mvmin;


    for(int imx=0; imx<n_mx; imx++){
      double xmx = double(imx)/(n_mx-1);
      double mx = mxmin*pow(mxmax/mxmin,xmx);
      if(n_mx==1) mx = mxmin;

      //#pragma omp parallel for
      for(int ie=0; ie<desteps; ie++){
        //double a=0;
        double xe = double(ie)/(desteps-1);
        double dE = demin*pow(demax/demin,xe);
        if(desteps==1){
          dE = de_target;
          ie = i_et;
        }

        double dsvdE = 0;
        double vmin = sqrt(dE*2/mx);
        for(int iv=0; iv<vsteps; iv++){
          double v = (iv+1)*dv;
          if(v<vmin) continue;
          double dsdE = dsdE_iEdEvum_qg(AKenq,ie,dE,v,mv,mx,qmin,qmax);
          if(plotv) ds_mv_mx_x[imv][imx][iv] = v*dsdE;
          dsvdE += arr_fv[iv]*v*dsdE;
        }//v
        dsvdE *= dv;
        //ofE<<dE*E_to_keV<<" "<<dsvdE*dsvdE_to_cm3keVday<<"\n";
        if(!plotv) ds_mv_mx_x[imv][imx][ie] = dsvdE;
      }//dE
    }//mx
  }//mv




  //Plot as function of v:

  ofv<<"# m_v blocks: ";
  for(int imv=0; imv<n_mv; imv++){
    double xmv = double(imv)/(n_mv-1);
    double mv = mvmin*pow(mvmax/mvmin,xmv);
    if(n_mv==1) mv = mvmin;
    ofv<<imv<<","<<mv*M_to_MeV<<" ";
  }
  ofv<<"\n";

  for(int imv=0; imv<n_mv; imv++){
    if(!plotv) break;
    double xmv = double(imv)/(n_mv-1);
    double mv = mvmin*pow(mvmax/mvmin,xmv);
    if(n_mv==1) mv = mvmin;

    ofv<<"\""<<std::fixed<<std::setprecision(2)<<mv*M_to_MeV<<" MeV\"   ";
    for(int imx=0; imx<n_mx; imx++){
      double xmx = double(imx)/(n_mx-1);
      double mx = mxmin*pow(mxmax/mxmin,xmx);
      if(n_mx==1) mx = mxmin;
      //ofv<<mx*M_to_GeV<<" ";
      ofv<<"\""<<std::setprecision(1)<<mx*M_to_GeV<<" GeV\"   ";
    }
    ofv<<"\n"<<std::scientific<<std::setprecision(6);

    for(int iv=0; iv<vsteps; iv++){
      double v = (iv+1)*dv;
      ofv<<v*V_to_kms<<" ";

      for(int imx=0; imx<n_mx; imx++){
        ofv<<ds_mv_mx_x[imv][imx][iv]*dsvdE_to_cm3keVday<<" ";
      }//mx
      ofv<<"\n";
    }//v
    ofv<<"\n";
  }//mv

  return 0;
}
