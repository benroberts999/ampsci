#include "AKF_akFunctions.h"
#include "SHM_standardHaloModel.h"
#include <iomanip>

////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////


//******************************************************************************
double g(double s, double x){
  double a = 0.398942/s;
  double y = (x/s)*(x/s);
  // if(y<0.2) return a*(1.-0.5*pow(x,2)+0.125*pow(x,4)-0.02083*pow(x,6)
  //  +0.0026042*pow(x,8)-0.00026042*pow(x,10));
  return a*exp(-0.5*y);
}


//******************************************************************************
double fv_au(double v_au, double cosphi=0, double dves=0, double dv0=0);
//******************************************************************************
double fv_au(double v_au, double cosphi, double dves, double dv0){
  double v = v_au * (FPC::c_SI/FPC::c); //will be in m/s
  v/=1.e3; //convert from m/s -> km/s
  return SHM::fv(v,cosphi,dves,dv0);
}


//******************************************************************************
template<typename T>
double dsdE_iEdEvum_qg(std::vector< std::vector< std::vector<T> > > &K_enq,
  int iE, double dE, double v, double mv, double mx,
  double qmin, double qmax)
/*
calcualtes cross-section ds/dE, for given E, v, mu, and mx.
Also needs min/max q (from grid) to re-construct q grids.
note: mu = mv c / hbar = mv/alpha!
*/
{
  //double A_au = 8*M_PI*FPC::c2; // alpha_chi = 1 in code.
  double A_au = 8*M_PI; //alpha_chi = eta*alpha; eta=1 in code
  double v2 = pow(v,2);

  int desteps    = (int) K_enq.size();
  int num_states = (int) K_enq[0].size();
  int qsteps     = (int) K_enq[0][0].size();
  if(iE>=desteps) return 0; //or -1?

  double mu = mv*FPC::c;

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
      #pragma omp critical (qint)
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

  //define input parameters
  std::string akfn; //name of K file to read in
  int vsteps;
  double mxmin,mxmax,mvmin,mvmax; //m_chi and m_v masses
  int n_mx,i_mv,n_mv;
  double de_target;
  std::string label="testx";
  bool plotv = false; //if single dE, plot fn of v!
  double Atot;
  double Ebi,Ebf,Ebw; // E bins: initial,final, width
  double cosp,dvesc,dv0;
  double dres;

  //Open and read the input file:
  {
    std::ifstream ifs;
    ifs.open("dmeXSection.in");
    std::string jnk;
    ifs >> akfn;                    getline(ifs,jnk);
    ifs >> vsteps;                  getline(ifs,jnk);
    ifs >> mxmin >> mxmax >> n_mx;  getline(ifs,jnk);
    ifs >> i_mv;                    getline(ifs,jnk);
    ifs >> mvmin >> mvmax >> n_mv;  getline(ifs,jnk);
    ifs >> cosp >> dvesc >> dv0;    getline(ifs,jnk);
    ifs >> dres;                    getline(ifs,jnk);
    ifs >> de_target;               getline(ifs,jnk);
    ifs >> Atot;                    getline(ifs,jnk);
    ifs >> Ebi>>Ebf>>Ebw;           getline(ifs,jnk);
    ifs >> label;                   getline(ifs,jnk);
    ifs.close();
  }

  if(cosp>1 || cosp<-1) return 1; //add message

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
    //std::cout<<mvmin<<" MeV = ";
    mvmin /= M_to_MeV;
    mvmax /= M_to_MeV;
    //std::cout<<mvmin<<" au ?\n ";
  }else{
    std::cout<<"Wrong m_v input given. Try again. Or not, whatevs\n";
    return 1;
  }
  //return 1;

  //convert E target to au
  de_target /= E_to_keV;

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

  //If only doing a single dE, find correct index for target dE:
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
    if(i_et >= desteps || i_et<0) std::cout<<"Bad E target\n";
    desteps = 1;
    std::cout<<"Signle dE: i="<<i_et<<" : "<<de_target<<"au = "
      <<de_target*E_to_keV<<"kev\n";
  }

  // Do q derivative on i grid:
  //double dqonq = log(qmax/qmin)/(qsteps-1); //need to multiply by q for dq

  //Grid of f_v(v). Can use to change vel profiles
  std::vector<double> arr_fv(vsteps);
  double max_v = (SHM::MAXV)/V_to_kms;
  double dv = max_v/vsteps;
  double fvNorm = SHM::normfv(cosp,dvesc,dv0);
  for(int i=0; i<vsteps; i++){
    double v = (i+1)*dv; //don't include zero
    arr_fv[i] = fvNorm*fv_au(v,cosp,dvesc,dv0);
  }


  //Print the grid info to screen:
  printf("\nq  grid: %6.2f -> %6.2f  MeV, %5i steps\n"
    ,qmin*Q_to_MeV,qmax*Q_to_MeV,qsteps);
  printf("E  grid: %6.2f -> %6.2f  keV, %5i steps\n"
    ,demin*E_to_keV,demax*E_to_keV,desteps);
  printf("v  grid: %6.2f -> %6.1f km/s, %5i steps\n"
    ,dv*V_to_kms,max_v*V_to_kms,vsteps);
  printf("Mx grid: %6.2f -> %6.1f  GeV, %5i steps\n"
    ,mxmin*M_to_GeV,mxmax*M_to_GeV,n_mx);
  if(mvmin<0) std::cout<<"Heavy meadiator\n";
  else printf("Mv grid: %6.3f -> %6.3f  MeV, %5i steps\n"
    ,mvmin*M_to_MeV,mvmax*M_to_MeV,n_mv);

  printf("\nA = %.1f au = %.2e cm^2/keV\n",A_au,A_cmkeV);
  std::cout<<"ds/dE conversion factor:   "<<dsdE_to_cm2keV<<" cm^2/keV\n"
    <<"ds.v/dE conversion factor: "
    <<dsvdE_to_cm3keVday<<"   cm^3/keV/day\n\n";

  // return 1;




  //Array to store cross-section
  // ds.v/dE (fun. of mv, mx, {v or E})
  std::vector< std::vector< std::vector<float> > > dsv_mv_mx_x;
  if(plotv) dsv_mv_mx_x.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,std::vector<float>(vsteps)));
  else dsv_mv_mx_x.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,std::vector<float>(desteps)));

  //Fill arrays:
  // If single E (plot v), calculate ds.v/dE for each v
  // If many E (no plot v), calc <ds.v>/dE for each E
  std::cout<<"Doing q integrations: ";
  printf("cos(phi)=%4.1f; dvesc=%4.1f; dv0=%4.1f.\n",cosp,dvesc,dv0);
  if(plotv) std::cout<<"Forming ds/dE as function of v:\n";
  else      std::cout<<"Forming <ds.v>/dE as function of dE:\n";
  std::cout.precision(1);
  std::cout<<std::fixed;
  for(int imv=0; imv<n_mv; imv++){
    double xmv = double(imv)/(n_mv-1);
    double mv = mvmin*pow(mvmax/mvmin,xmv);
    if(n_mv==1) mv = mvmin;
    for(int imx=0; imx<n_mx; imx++){
      double xmx = double(imx)/(n_mx-1);
      double mx = mxmin*pow(mxmax/mxmin,xmx);
      if(n_mx==1) mx = mxmin;
      std::cout<<"mv = "<<mv*M_to_MeV<<" MeV; mx = "<<mx*M_to_GeV<<" GeV\n";
      //#pragma omp parallel for
      for(int ie=0; ie<desteps; ie++){
        //double a=0;
        double xe = double(ie)/(desteps-1);
        double dE = demin*pow(demax/demin,xe);
        if(desteps==1){
          dE = de_target;
          ie = i_et;
        }else{
          //std::cout<<"\rE: "<<dE<<"/"<<demax<<" au       ";
          double pc = 100.*(ie+1.)/desteps;
          printf("\rE: %5.1f/%5.1f au  -  %5.1f%%        ",dE,demax,pc);
          std::cout<<std::flush;
        }
        double dsvdE = 0;
        double vmin = sqrt(dE*2/mx);
        for(int iv=0; iv<vsteps; iv++){
          double v = (iv+1)*dv;
          if(v<vmin) continue;
          double dsdE = dsdE_iEdEvum_qg(AKenq,ie,dE,v,mv,mx,qmin,qmax);
          if(plotv) dsv_mv_mx_x[imv][imx][iv] = dsdE; //nb: no v! XXX v^2?
          //#pragma omp critical (Eloop)
          {
            dsvdE += arr_fv[iv]*v*dsdE;
          }
        }//v
        dsvdE *= dv;
        if(!plotv) dsv_mv_mx_x[imv][imx][ie] = dsvdE;
      }//dE
      if(desteps!=1) std::cout<<"\n";
    }//mx
  }//mv
  std::cout<<"\n";

  if(Ebw==0){
    //Open relevant output files
    std::ofstream of;
    std::string fname,str_E;
    str_E = std::to_string(int(de_target*E_to_keV*1000));
    if(plotv) fname = "plot-sdE_v-"+str_E+"-"+label+".txt";
    else      fname = "plot-svdE_dE-"+label+".txt";
    of.open(fname.c_str());

    if(plotv) std::cout<<"Plotting ds/dE as function of v:\n";
    else      std::cout<<"Plotting <ds.v>/dE as function of dE:\n";
    std::cout<<"Filename: "<<fname<<"\n\n";

    //Output plots. Function of v, or dE
    of<<"# m_v blocks: ";
    for(int imv=0; imv<n_mv; imv++){
      double xmv = double(imv)/(n_mv-1);
      double mv = mvmin*pow(mvmax/mvmin,xmv);
      if(n_mv==1) mv = mvmin;
      of<<imv<<","<<mv*M_to_MeV<<" ";
    }
    of<<"\n";
    for(int imv=0; imv<n_mv; imv++){
      double xmv = double(imv)/(n_mv-1);
      double mv = mvmin*pow(mvmax/mvmin,xmv);
      if(n_mv==1) mv = mvmin;
      of<<"\""<<std::fixed<<std::setprecision(2)<<mv*M_to_MeV<<" MeV\"   ";
      for(int imx=0; imx<n_mx; imx++){
        double xmx = double(imx)/(n_mx-1);
        double mx = mxmin*pow(mxmax/mxmin,xmx);
        if(n_mx==1) mx = mxmin;
        of<<"\""<<std::setprecision(1)<<mx*M_to_GeV<<" GeV\"   ";
      }
      of<<"\n"<<std::scientific<<std::setprecision(6);
      for(int iv=0; iv<vsteps; iv++){
        if(!plotv) break;
        double v = (iv+1)*dv;
        of<<v*V_to_kms<<" ";
        for(int imx=0; imx<n_mx; imx++){
          //output ds/dE [not ds.v/dE]
          of<<(dsv_mv_mx_x[imv][imx][iv])*dsdE_to_cm2keV<<" ";
          //of<<(dsv_mv_mx_x[imv][imx][iv])*dsvdE_to_cm3keVday<<" ";
        }//mx
        of<<"\n";
      }//v
      for(int ie=0; ie<desteps; ie++){
        if(plotv) break;
        double xe = double(ie)/(desteps-1);
        double dE = demin*pow(demax/demin,xe);
        of<<dE*E_to_keV<<" ";
        for(int imx=0; imx<n_mx; imx++){
          of<<dsv_mv_mx_x[imv][imx][ie]*dsvdE_to_cm3keVday<<" ";
        }//mx
        of<<"\n";
      }//v
      of<<"\n";
    }//mv
    of.close();
  }


  if(desteps==1) return 0; //finished
  if(Ebw==0) return 0; //finished
  std::cout.precision(4);
  std::cout<<std::scientific;


  //XXX output dsv_mv_mx_x as binary
  //(only for all E - NOT v!!)
  //**************************************************************************


  int imv = 0;
  int imx = 0;
  double dEonE = log(demax/demin)/(desteps-1);

  double mx = mxmin*pow(mxmax/mxmin,double(imx)/(n_mx-1));

  //calculate Gaussian smearing (DAMA)
  //std::vector<double> y(desteps);
  std::vector<float> y = dsv_mv_mx_x[imv][imx];

  //calculate Gaussian smearing (DAMA)
  // Also: include detector efficiency!? Can ignore for dama..
  for(int i = 0; i<desteps; i++){
    double E = demin*pow(demax/demin,double(i)/(desteps-1));
    double y0 = 0;
    double alph = 0.45 + dres*0.04;
    double beta = 0.009 + dres*0.005;
    //triple check this! Super important!
    double s = (alph*sqrt(E*E_to_keV) + beta*(E*E_to_keV))/E_to_keV;
    //std::cout<<E*E_to_keV<<" "<<s*E_to_keV<<" "<<s<<"\n";
    for(int j = 0; j<desteps; j++){
      double Ep = demin*pow(demax/demin,double(j)/(desteps-1));
      y0 += g(s,E-Ep)*dsv_mv_mx_x[imv][imx][j]*Ep;
    }
    y[i] = y0*dEonE;
  }

  // //no smearing:
  // for(int i = 0; i<desteps; i++){
  //   y[i] = dsv_mv_mx_x[imv][imx][i];
  // }

  double MN = Atot*(FPC::u_NMU*FPC::m_e_kg); //Total atomic/mol. mass (in kg)
  double rho_on_mxc2 = rhoDM_GeVcm3/(mx*M_to_GeV);
  double rateFac = dsvdE_to_cm3keVday*rho_on_mxc2/MN;
  //in (cm^3/kev/day)/(cm^3 *kg) = 1/keV/kg/day

  // //double check units: flux of DM particles, per second
  // std::cout<<"DM flux through 1cm^2 per second: ~";
  // std::cout<<300.e5*rhoDM_GeVcm3/(mx*M_to_GeV)<<"\n\n";

  //Integrate smeared rate over energy bins
  std::cout<<"Input file: "<<akfn<<"\n";
  std::cout<<"Integrated energy bins:\n";
  printf("cos(phi)=%4.1f; dvesc=%4.1f; dv0=%4.1f.\n",cosp,dvesc,dv0);
  std::cout<<"E_a-E_b   E(kev)  R (cpd/kg/keV)\n";
  double Ea=Ebi/E_to_keV;
  while(Ea<(Ebf + Ebw)/E_to_keV){

    double Eb = Ea + Ebw/E_to_keV;
    double tmpa = (desteps-1)*log(Ea/demin)/log(demax/demin);
    double tmpb = (desteps-1)*log(Eb/demin)/log(demax/demin);
    int ieA = (int) ceil(tmpa);
    int ieB = (int) ceil(tmpb);
    if(ieB>=desteps) break;
    if(ieA<0)ieA=0;

    double Rate=0;
    for(int ie = ieA; ie<ieB; ie++){
      double xe = double(ie)/(desteps-1);
      double E = demin*pow(demax/demin,xe);
      Rate += y[ie]*E; //nb: E is from Jacobian; * dE/E below
    }
    Rate *= rateFac*dEonE/Ebw;
    // std::cout<<dEonE<<" "<<Rate<<"\n";
    printf("%3.1f-%3.1f: %6.3f   %.2e\n",
    Ea*E_to_keV,Eb*E_to_keV,0.5*(Ea+Eb)*E_to_keV,Rate);
    Ea += Ebw/E_to_keV;
  }


  return 0;
}
