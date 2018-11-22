#include "AKF_akFunctions.h"
#include "StandardHaloModel.h"
#include "SHM_standardHaloModel.h"
#include <iomanip>

//Convert FROM a.u. TO rel./other units
double E_to_keV = FPC::Hartree_eV/1000.; //au -> keV (energy)
double M_to_GeV = FPC::m_e_MeV/1000.;
double M_to_MeV = FPC::m_e_MeV;
double V_to_kms = (FPC::c_SI/FPC::c)/1000.;
double V_to_cms = (FPC::c_SI*FPC::alpha)*100.; // au -> cm/s
double V_to_cmday = V_to_cms*(24*60*60); // cm/s -> cm/day
double Q_to_MeV = FPC::Hartree_eV*FPC::c/1.e6;

//DM energy density (in GeV/cm^3)
double rhoDM_GeVcm3 = 0.4; // GeV/cm^3

double sbe_1e37_cm2 = 1.e-37;
//With with sig-bar_e = sbe_cm2 = 1.e-37cm^2;
double dsdE_to_cm2keV = sbe_1e37_cm2/E_to_keV; // au -> cm^2/keV
double dsvdE_to_cm3keVday = dsdE_to_cm2keV*V_to_cmday;


//******************************************************************************
struct ExpGrid{
  int N;
  double min;
  double max;
  double dxonx;
  double x(int i);
  double dxdi(int i);
  int findNextIndex(double x);
  ExpGrid(int in_N, double in_min, double in_max);
};

ExpGrid::ExpGrid(int in_N, double in_min, double in_max){
  N = in_N;
  min = in_min;
  max = in_max;
  if(min>max){
    min = in_max;
    max = in_min;
  }
  if(in_N<=1){
    max=min;
    N=1;
  }
  dxonx = log(max/min)/(N-1);
}

double ExpGrid::x(int i){
  double y = double(i)/(N-1);
  if(N==1) return min;
  return min*pow(max/min,y);
}
double ExpGrid::dxdi(int i){
  return dxonx*x(i);
}

int ExpGrid::findNextIndex(double x){
  double tmp = (N-1)*log(x/min)/log(max/min);
  int i = (int) ceil(tmp);
  return i;
}

//******************************************************************************
double g(double s, double x)
// Simple Gaussian
{
  double a = 0.398942/s;
  double y = (x/s)*(x/s);
  return a*exp(-0.5*y);
}

//******************************************************************************
double fv_au(double v_au, double cosphi, double dves, double dv0)
// SHM vel. distribution, v in atomic units.
{
  double v = v_au * (FPC::c_SI/FPC::c); //will be in m/s
  v/=1.e3; //convert from m/s -> km/s
  return SHM::fv(v,cosphi,dves,dv0);
}

//******************************************************************************
template<typename T>
double dsdE_Evmvmx(std::vector< std::vector<T> > &Ke_nq,
  double E, double v, double mv, double mx, ExpGrid &qgrid)
/*
Note: just takes in _PART_ of K (for given E)
calcualtes cross-section ds/dE, for given E, v, mu, and mx.
Does q integrations, and sums over states
note: mu = mv c / hbar = mv/alpha!
Output in units of sig-bar_e
*/
{
  double arg = pow(mx*v,2)-2.*mx*E;
  if(arg<0) return 0;

  int num_states = (int) Ke_nq.size();
  int qsteps     = qgrid.N;
  double dqonq = qgrid.dxonx;

  double mu = mv*FPC::c;

  bool finite_med = true;
  if(mu<0) finite_med = false;

  double qminus = mx*v - sqrt(arg);
  double qplus  = mx*v + sqrt(arg);
  if(qminus>qgrid.max || qplus<qgrid.min) return 0;

  double dsdE = 0;
  for(int ink=0; ink<num_states; ink++){
    for(int iq=0; iq<qsteps; iq++){
      double q = qgrid.x(iq);
      if(q<qminus || q>qplus) continue;
      double qdq_on_dqonq = q*q; //(dq/q) is constant, multiply at end
      double F_chi = 1.;
      if(finite_med) F_chi = 1./pow(q*q+mu*mu,2);
      dsdE += qdq_on_dqonq*F_chi*Ke_nq[ink][iq];  //dq/q included below
    }//q int
  }//states

  double A = 0.5; //in units of sig-bar_e !
  dsdE *= (A/pow(v,2))*dqonq;
  return dsdE;
}

//******************************************************************************
template<typename T>
double dsvdE_Evmvmx(std::vector< std::vector<T> > &Ke_nq,
  double E, double mv, double mx, ExpGrid &qgrid,
  std::vector<double> &arr_fv, double dv)
/*
Calculates <ds.v>/dE for given E, mv, mx
Does the v integration
Note: only takes _part_ of the K array! (for given E)
*/
{
  int vsteps = (int) arr_fv.size();
  double vmin = sqrt(2*E/mx);
  double dsvdE = 0;
  for(int iv=0; iv<vsteps; iv++){
    double v = (iv+1)*dv;
    if(v<vmin) continue;
    double dsdE = dsdE_Evmvmx(Ke_nq,E,v,mv,mx,qgrid);
    dsvdE += arr_fv[iv]*v*dsdE;
  }//v
  return dsvdE*dv;
}

//******************************************************************************
template<typename T>
void form_dsvdE(
  std::vector<float> &dsvde,
  std::vector< std::vector< std::vector<T> > > &K_enq,
  double mv, double mx,
  ExpGrid &Egrid, ExpGrid &qgrid,
  std::vector<double> &arr_fv, double dv)
/*
Forms <dsv>/de array (for given mx, mv)
Note: mv<0 means "heavy" mediator [Fx=1]
*/
{
  //Loop through E, create dsvde array
  int desteps = Egrid.N;
  #pragma omp parallel for
  for(int ie=0; ie<desteps; ie++){
    double E = Egrid.x(ie);
    //Do v (and q) integrations:
    double dsvdE = dsvdE_Evmvmx(K_enq[ie],E,mv,mx,qgrid,arr_fv,dv);
    dsvde[ie] = dsvdE;
  }//dE
}

//******************************************************************************
void writeForGnuplot_mvBlock(
  std::vector< std::vector< std::vector<float> > > &X_mv_mx_x,
  ExpGrid mvgrid, ExpGrid mxgrid, ExpGrid Egrid,
  std::string fname)
{

  std::ofstream of(fname.c_str());

  int n_mv = mvgrid.N;
  int n_mx = mxgrid.N;
  int desteps = Egrid.N;

  of<<"# m_v blocks: ";
  for(int imv=0; imv<n_mv; imv++){
    double mv = mvgrid.x(imv);
    of<<imv<<","<<mv*M_to_MeV<<" ";
  }
  of<<"\n";

  for(int imv=0; imv<n_mv; imv++){
    double mv = mvgrid.x(imv);
    of<<"\""<<std::fixed<<std::setprecision(2)<<mv*M_to_MeV<<" MeV\"   ";
    for(int imx=0; imx<n_mx; imx++){
      double mx = mxgrid.x(imx);
      of<<"\""<<std::setprecision(1)<<mx*M_to_GeV<<" GeV\"   ";
    }
    of<<"\n"<<std::scientific<<std::setprecision(6);
    for(int ie=0; ie<desteps; ie++){
      double E = Egrid.x(ie);
      of<<E*E_to_keV<<" ";
      for(int imx=0; imx<n_mx; imx++){
        of<<X_mv_mx_x[imv][imx][ie]*dsvdE_to_cm3keVday<<" ";
      }//mx
      of<<"\n";
    }//E
    of<<"\n";
  }//mv

  of.close();
}
//******************************************************************************
void writeForGnuplot_mxBlock(
  std::vector< std::vector< std::vector<float> > > &X_mv_mx_x,
  ExpGrid mvgrid, ExpGrid mxgrid, ExpGrid Egrid,
  std::string fname)
{

  std::ofstream of(fname.c_str());

  int n_mv = mvgrid.N;
  int n_mx = mxgrid.N;
  int desteps = Egrid.N;

  of<<"# m_x blocks: ";
  for(int imx=0; imx<n_mx; imx++){
    double mx = mxgrid.x(imx);
    of<<imx<<","<<mx*M_to_GeV<<" ";
  }
  of<<"\n";

  for(int imx=0; imx<n_mx; imx++){
    double mx = mxgrid.x(imx);
    of<<"\""<<std::fixed<<std::setprecision(2)<<mx*M_to_GeV<<" GeV\"   ";
    for(int imv=0; imv<n_mv; imv++){
      double mv = mvgrid.x(imv);
      of<<"\""<<std::setprecision(1)<<mv*M_to_MeV<<" MeV\"   ";
    }
    of<<"\n"<<std::scientific<<std::setprecision(6);
    for(int ie=0; ie<desteps; ie++){
      double E = Egrid.x(ie);
      of<<E*E_to_keV<<" ";
      for(int imv=0; imv<n_mv; imv++){
        of<<X_mv_mx_x[imv][imx][ie]*dsvdE_to_cm3keVday<<" ";
      }//mv
      of<<"\n";
    }//E
    of<<"\n";
  }//mx

  of.close();
}



//******************************************************************************
//******************************************************************************
//******************************************************************************
int main(void){

  //define input parameters
  std::string akfn; //name of K file to read in
  int vsteps;
  double mxmin,mxmax,mvmin,mvmax; //m_chi and m_v masses
  int n_mx,i_mv,n_mv;
  //double de_target;
  std::string label="testx";
  //bool plotv = false; //if single dE, plot fn of v!
  double Atot;
  double iEbin,fEbin,wEbin; // E bins: initial,final, width
  double cosp,dvesc,dv0;
  double dres,err_PEkeV; //detector resolution, PE-keV errors [-1]

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
    ifs >> dres >> err_PEkeV;       getline(ifs,jnk);
    ifs >> Atot;                    getline(ifs,jnk);
    ifs >> iEbin >> fEbin >> wEbin; getline(ifs,jnk);
    ifs >> label;                   getline(ifs,jnk);
    ifs.close();
  }

  if(cosp>1 || cosp<-1) return 1; //add message

  //DM mass: Convert from GeV to au:
  mxmin /= M_to_GeV;
  mxmax /= M_to_GeV;
  ExpGrid mxgrid(n_mx,mxmin,mxmax);

  //Mediator mass: convert units + set-up finite/infinite case
  if(i_mv==0){
    //massless case. Do sepperately
    mvmin = mvmax = 0;
    n_mv = 1;
  }else if(i_mv==2){
    //Heavy-mediator case (contact interaction)
    mvmin = mvmax = -1; //1./0.;
    n_mv = 1;
  }else{
    mvmin /= M_to_MeV;
    mvmax /= M_to_MeV;
  }
  ExpGrid mvgrid(n_mv,mvmin,mvmax);

  //Energy bins for integrating/averaging
  iEbin/=E_to_keV;
  fEbin/=E_to_keV;
  wEbin/=E_to_keV;

  //Arrays/values to be filled from input AK file:
  std::vector< std::vector< std::vector<float> > > Kenq;
  std::vector<std::string> nklst;
  double qmin,qmax,demin,demax;
  //Read in AK file
  std::cout<<"Opening file: "<<akfn<<".bin\n";
  AKF::akReadWrite(akfn,false,Kenq,nklst,qmin,qmax,demin,demax);
  int desteps    = (int) Kenq.size();
  int num_states = (int) Kenq[0].size();
  int qsteps     = (int) Kenq[0][0].size();
  if(num_states != (int)nklst.size()) return 1; //just sanity check
  //Create the E and q grids:
  ExpGrid Egrid(desteps,demin,demax);
  ExpGrid qgrid(qsteps,qmin,qmax);

  //Grid of f_v(v). Can use to change vel profiles
  std::vector<double> arr_fv(vsteps);
  double max_v = (SHM::MAXV)/V_to_kms;
  double dv = max_v/vsteps;
  double fvNorm = SHM::normfv(cosp,dvesc,dv0);
  for(int i=0; i<vsteps; i++){
    double v = (i+1)*dv; //don't include zero
    arr_fv[i] = fvNorm*fv_au(v,cosp,dvesc,dv0);
  }

  std::vector<double> arr_fv2(vsteps);
  StandardHaloModel shm(cosp,dvesc,dv0);
  for(int i=0; i<vsteps; i++){
    double v = (i+1)*dv; //don't include zero
    double vkms = v * (FPC::c_SI/FPC::c); //will be in m/s
    vkms/=1.e3; //convert from m/s -> km/s
    arr_fv2[i] = shm.fv(vkms);
  }

  double a=0, b=0;
  for(int i=0; i<vsteps; i++){
    double v = (i+1)*dv; //don't include zero
    std::cout<<v<<" "<<arr_fv[i]<<" "<<arr_fv2[i]<<" "<<arr_fv[i]-arr_fv2[i]<<"\n";
    a += arr_fv2[i]*dv*(FPC::c_SI/FPC::c)/1.e3;
    b += arr_fv[i]*dv;
  }
  std::cout<<a<<" "<<b<<"\n";

  std::cout<<"NOTE: Fchi NOT normalised correctly??\n";

  return 1;

  //Print the grid info to screen:
  printf("\nq :%6.2f -> %6.2f MeV, N=%4i\n",qmin*Q_to_MeV,qmax*Q_to_MeV,qsteps);
  printf("E :%6.2f -> %6.2f keV, N=%4i\n"
    ,demin*E_to_keV,demax*E_to_keV,desteps);
  printf("v :%6.2f -> %6.2fkm/s, N=%4i\n"
    ,dv*V_to_kms,max_v*V_to_kms,vsteps);
  printf("Mx:%6.2f -> %6.2f GeV, N=%4i\n"
    ,mxmin*M_to_GeV,mxmax*M_to_GeV,n_mx);
  if(mvmin<0) std::cout<<"Heavy meadiator\n";
  else printf("Mv:%6.2f -> %6.2f MeV, N=%4i\n"
    ,mvmin*M_to_MeV,mvmax*M_to_MeV,n_mv);

  // Units + conversions for dsvde..etc
  std::cout<<"Doing calculations with sig-bar_e =  "<<sbe_1e37_cm2<<"cm2\n"
    <<"ds/dE conversion factor:   "<<dsdE_to_cm2keV<<" cm^2/keV\n"
    <<"ds.v/dE conversion factor: "
    <<dsvdE_to_cm3keVday<<"   cm^3/keV/day\n\n";

  double al_x = (sqrt(sbe_1e37_cm2)/FPC::aB_cm)*(FPC::alpha/sqrt(16.*M_PI));
  double al_xH = al_x/(FPC::m_e_MeV*FPC::alpha2);
  std::cout<<"Note: sig-bar_e =  "<<sbe_1e37_cm2<<"cm2 corresponds to:\n"
    <<" Light mediator: al_x = "<<al_x<<"\n"
    <<" Heavy mediator: al_x = "<<al_xH<<"*(mv/MeV)^2\n\n";

  //Array to store cross-section
  // ds.v/dE (fun. of mv, mx, E)
  std::vector< std::vector< std::vector<float> > > dsv_mv_mx_E;
  dsv_mv_mx_E.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,
      std::vector<float>(desteps)
    )
  );

  //Fill arrays:
  // If single E (plot v), calculate ds.v/dE for each v
  // If many E (no plot v), calc <ds.v>/dE for each E
  std::cout<<"Doing q and v integrations: ";
  printf("cos(phi)=%4.1f; dvesc=%4.1f; dv0=%4.1f.\n",cosp,dvesc,dv0);
  for(int imv=0; imv<n_mv; imv++){
    double mv = mvgrid.x(imv);
    for(int imx=0; imx<n_mx; imx++){
      double mx = mxgrid.x(imx);
      printf("M_chi=%5.2f GeV",mx*M_to_GeV);
      if(mv>=0) printf(" ; M_v=%6.3f MeV",mv*M_to_MeV);
      std::cout<<" .. "<<std::flush;
      form_dsvdE(dsv_mv_mx_E[imv][imx],Kenq,mv,mx,Egrid,qgrid,arr_fv,dv);
      std::cout<<" .. Done\n";
    }//mx
  }//mv
  std::cout<<"\n";

  bool write_dsvde = true;

  //output <ds.v>/dE for gnuplot:
  if(write_dsvde){
    std::string fn_dsvde = "plot-dsvde_mx-"+label+".out";
    std::cout<<"Writing to file: "<<fn_dsvde<<"\n";
    writeForGnuplot_mvBlock(dsv_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dsvde);
    if(n_mv>1){
      fn_dsvde = "plot-dsvde_mv-"+label+".out";
      std::cout<<"Writing to file: "<<fn_dsvde<<"\n";
      writeForGnuplot_mvBlock(dsv_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dsvde);
    }
  }


  // *********************
  //    Here, is just for DAMA! Move into sepperate function!!
  // *********************


  //Array to store observable Rate, S
  std::vector< std::vector< std::vector<float> > > dSdE_mv_mx_E;
  dSdE_mv_mx_E.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,
      std::vector<float>(desteps)
    )
  );


  //Hardware threshold: should be between [-1,1]
  double PE_per_keV = 6.5 + err_PEkeV*1.;
  double E_thresh_HW = 1./PE_per_keV/E_to_keV;

  //calculate Gaussian smearing (DAMA)
  double alpha = 0.45 + dres*0.04;
  double beta = 0.009 + dres*0.005;

  double dEonE = Egrid.dxonx;
  double MN = Atot*(FPC::u_NMU*FPC::m_e_kg); //Total atomic/mol. mass (in kg)

  for(int imv=0; imv<n_mv; imv++){
    for(int imx=0; imx<n_mx; imx++){
      double mx = mxgrid.x(imx);
      double rho_on_mxc2 = rhoDM_GeVcm3/(mx*M_to_GeV);
      double rateFac = dsvdE_to_cm3keVday*rho_on_mxc2/MN;
      for(int i = 0; i<desteps; i++){
        double E = Egrid.x(i);
        if(E<E_thresh_HW){
          dSdE_mv_mx_E[imv][imx][i] = 0;
          continue;
        }
        //triple check this! Super important!
        double s = (alpha*sqrt(E*E_to_keV) + beta*(E*E_to_keV))/E_to_keV;
        double y0 = 0;
        //integrate over Eprime:
        for(int j = 0; j<desteps; j++){
          double Ep = Egrid.x(j);
          if(Ep<E_thresh_HW) continue; //hardware threshold
          y0 += g(s,E-Ep)*dsv_mv_mx_E[imv][imx][j]*Ep;
        }
        dSdE_mv_mx_E[imv][imx][i] = y0*dEonE*rateFac;
      }
    }
  }


  bool write_dS = true;
  //output dS/dE for gnuplot:
  if(write_dS){
    std::string fn_dSde = "plot-dSdE_mx-"+label+".out";
    std::cout<<"Writing to file: "<<fn_dSde<<"\n";
    writeForGnuplot_mvBlock(dSdE_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dSde);
    if(n_mv>1){
      fn_dSde = "plot-dSdE_mv-"+label+".out";
      std::cout<<"Writing to file: "<<fn_dSde<<"\n";
      writeForGnuplot_mvBlock(dSdE_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dSde);
    }
  }


  // Integrate/average into energy bins:
  int num_bins = ceil((fEbin-iEbin)/wEbin);

  //Array to store averaged Rate, S
  std::vector< std::vector< std::vector<float> > > S_mv_mx_E;
  S_mv_mx_E.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,
      std::vector<float>(num_bins)
    )
  );

  std::cout<<"\nIntegrated (averaged) each energy bin\n";
  std::cout<<"(Just outputting for first mx/mv: ";
    printf("M_chi=%5.2f GeV",mxmin*M_to_GeV);
    if(mvmin>=0) printf(" ; M_v=%6.3f MeV",mvmin*M_to_MeV);
  std::cout<<"\n";
  for(int imv=0; imv<n_mv; imv++){
    for(int imx=0; imx<n_mx; imx++){
      for(int i=0; i<num_bins; i++){
        int ieA = Egrid.findNextIndex(iEbin+i*wEbin);
        int ieB = Egrid.findNextIndex(iEbin+(i+1)*wEbin);
        if(ieA<0) ieA=0;

        double Rate = 0;
        for(int ie = ieA; ie<ieB; ie++){
          if(ieB>=desteps) break;
          double E = Egrid.x(ie);
          Rate += dSdE_mv_mx_E[imv][imx][ie]*E;
          //nb: E is from Jacobian; * dE/E below
        }
        Rate *= dEonE/wEbin;
        S_mv_mx_E[imv][imx][i] = Rate;
        double EaKev = Egrid.x(ieA)*E_to_keV;
        double EbKev = Egrid.x(ieB)*E_to_keV;
        if(imv==0 && imx==0){
          //only print first one to screen
          printf("%3.1f-%3.1f: %6.3f   %.2e     %i\n",
          EaKev,EbKev,0.5*(EaKev+EbKev),Rate,ieB-ieA);
        }
      }
    }
  }

  return 0;
}
