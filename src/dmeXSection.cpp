#include "AKF_akFunctions.h"
#include "StandardHaloModel.h"
#include "FPC_physicalConstants.h"
#include "ChronoTimer.h"
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
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


// *** === *** === *** === *** === *** === *** === *** === *** === *** === ***
// Some macro/typedef short cuts
using FloatVec3D = std::vector< std::vector< std::vector<float> > >;
using FloatVec2D = std::vector< std::vector<float> >;

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

// double dsdE = dsdE_Evmvmx(Ke_nq,E,v,mv,mx,qgrid);
//******************************************************************************
template<typename T>
double dsdE_Evmvmx(const std::vector< std::vector<T> > &Ke_nq,
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
  int num_states = (int) (Ke_nq.size());
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

// double dsvdE = dsvdE_Evmvmx(K_enq[ie],E,mv,mx,qgrid,arr_fv,dv);
//******************************************************************************
template<typename T>
double dsvdE_Evmvmx(
  const std::vector< std::vector<T> > &Ke_nq,
  double E, double mv, double mx,
  ExpGrid &qgrid,
  const std::vector<double> &arr_fv, double dv)
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
  const std::vector< std::vector< std::vector<T> > > &K_enq,
  double mv, double mx,
  ExpGrid &Egrid, ExpGrid &qgrid,
  const std::vector<double> &arr_fv, double dv)
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
  const FloatVec3D &X_mv_mx_x,
  ExpGrid mvgrid, ExpGrid mxgrid, ExpGrid Egrid,
  std::string fname, double y_unit_convert)
{
//dsvdE_to_cm3keVday
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
        of<<X_mv_mx_x[imv][imx][ie]*y_unit_convert<<" ";
      }//mx
      of<<"\n";
    }//E
    of<<"\n";
  }//mv

  of.close();
}
//******************************************************************************
void writeForGnuplot_mxBlock(
  const FloatVec3D &X_mv_mx_x,
  ExpGrid &mvgrid, ExpGrid &mxgrid, ExpGrid &Egrid,
  const std::string &fname)
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
template<typename T>
void calculate_dsvde_array(
  const std::vector< std::vector< std::vector<T> > > &Kenq,
  FloatVec3D &dsv_mv_mx_E,
  ExpGrid &mvgrid, ExpGrid &mxgrid, ExpGrid &Egrid, ExpGrid &qgrid,
  const std::vector< std::vector<double> > &arr_fv, double dv,
  bool do_anMod)
{

  int n_mv = mvgrid.N;
  int n_mx = mxgrid.N;
  int desteps = Egrid.N;

  dsv_mv_mx_E.resize(n_mv,std::vector< std::vector<float> >(n_mx,
      std::vector<float>(desteps)));

  FloatVec3D dsv_mv_mx_Emax;
  if(do_anMod){
    dsv_mv_mx_Emax.resize(n_mv,std::vector< std::vector<float> >(n_mx,
      std::vector<float>(desteps)));
  }

  // Calculate <ds.v>/dE for each E (for each mx, mv)
  std::cout<<"Calculating <ds.v>/dE (doing q and v integrations): \n";
  for(int imv=0; imv<n_mv; imv++){
    double mv = mvgrid.x(imv);
    #pragma omp parallel for if(n_mx>15)
    for(int imx=0; imx<n_mx; imx++){
      double mx = mxgrid.x(imx);
      //printf("M_chi=%5.2f GeV",mx*M_to_GeV);
      //if(mv>=0) printf(" ; M_v=%6.3f MeV",mv*M_to_MeV);
      std::cout<<"\r .. "<<std::flush;
      if(do_anMod){
        form_dsvdE(dsv_mv_mx_E[imv][imx],Kenq,mv,mx,Egrid,qgrid,arr_fv[1],dv);
        std::cout<<"\r  ... "<<std::flush;
        form_dsvdE(dsv_mv_mx_Emax[imv][imx],Kenq,mv,mx,Egrid,qgrid,arr_fv[2]
          ,dv);
      }else{
        form_dsvdE(dsv_mv_mx_E[imv][imx],Kenq,mv,mx,Egrid,qgrid,arr_fv[0],dv);
      }
      std::cout<<"\r .. .. "<<std::flush;
    }//mx
  }//mv
  std::cout<<"\n";

  //If doing an. mod, form "amplitude" for ds.v/dE, used below
  if(do_anMod){
    std::cout<<"\nAnnual modulation.\n";
    std::cout<<"Forming <ds.v>_mod = (<ds.v>_max - <ds.v>_min)/2 \n";
    for(int imv=0; imv<n_mv; imv++){
      for(int imx=0; imx<n_mx; imx++){
        for(int ie=0; ie<desteps; ie++){
          dsv_mv_mx_E[imv][imx][ie] =
           0.5*(dsv_mv_mx_Emax[imv][imx][ie]-dsv_mv_mx_E[imv][imx][ie]);
        }
      }
    }
  }

}



//******************************************************************************
void doDAMA(const FloatVec3D &dsv_mv_mx_E,
  ExpGrid &mvgrid, ExpGrid &mxgrid, ExpGrid &Egrid,
  bool do_anMod, bool write_dSdE, bool write_SofM,
  double dres, double err_PEkeV, double Atot,
  double iEbin, double fEbin, double wEbin,
  std::string label
)
{
  std::cout<<"\n*** Doing S1 for DAMA ***\n";

  int n_mv = mvgrid.N;
  int n_mx = mxgrid.N;
  int desteps = Egrid.N;

  //Array to store observable Rate, S
  FloatVec3D dSdE_mv_mx_E;
  dSdE_mv_mx_E.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,
      std::vector<float>(desteps)
    )
  );

  //Hardware threshold: should be between [-1,1]
  double PE_per_keV = 6.5 + err_PEkeV*1.;
  double E_thresh_HW = 1./PE_per_keV/E_to_keV;

  //DAMA parameters for Gaussian resolution (smearing)
  double alpha = 0.45 + dres*0.04;
  double beta = 0.009 + dres*0.005;

  double dEonE = Egrid.dxonx;
  double MN = Atot*(FPC::u_NMU*FPC::m_e_kg); //Total atomic/mol. mass (in kg)

  //Calculate _observable_ rate, S, for DAMA
  // inlcuding Gaussian resolution, +
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

  //(Only write dS/dE if not writing )
  //output dS/dE for gnuplot:
  if(write_dSdE){
    std::string spref = do_anMod ? "dSmdE" : "dSdE";
    std::string fn_dSde = spref+"_mx-"+label+".out";
    std::cout<<"Writing to file: "<<fn_dSde<<"\n";
    double u = 1; //already converted!
    if(n_mv==1 || n_mx>1){
      writeForGnuplot_mvBlock(dSdE_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dSde,u);
    }
    if(n_mv>1){
      fn_dSde = spref+"_mv-"+label+".out";
      std::cout<<"Writing to file: "<<fn_dSde<<"\n";
      writeForGnuplot_mvBlock(dSdE_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dSde,u);
    }
  }

  // Integrate/average into energy bins:
  int num_bins = ceil((fEbin-iEbin)/wEbin);

  //Array to store averaged Rate, S
  FloatVec3D S_mv_mx_E;
  S_mv_mx_E.resize(n_mv,
    std::vector< std::vector<float> >(n_mx,
      std::vector<float>(num_bins)
    )
  );

  //Creates S (or Sm) [actually S/E_w] - rate averaged over each energy bin
  //Also prints energy bins to screen (only for first mx, mv)
  std::cout<<"\nCalculating ";
  if(do_anMod)std::cout<<"Sm/Ew, annual modulation amplitude (/day/kg/keV)\n";
  else        std::cout<<"S/Ew, average event rate (/day/kg/keV)\n";
  std::cout<<"Integrated (averaged) each energy bin.\n";
  std::cout<<"(Just outputting for first mx/mv: ";
    printf("M_chi=%5.2f GeV",mxgrid.min*M_to_GeV);
    if(mvgrid.min>=0) printf(" ; M_v=%6.3f MeV",mvgrid.min*M_to_MeV);
  std::cout<<")\n";
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


  // For each energy bin, output as function of m_chi.
  // Note: if doing this, don't output above files! (too many m_chi's)
  // Each collumn a different energy bin
  //Each m_v
  if(write_SofM){
    std::string fn_S = do_anMod ? "Sm" : "S";
    fn_S = fn_S+"_mx-"+label+".out";
    std::cout<<"Writing to file: "<<fn_S<<"\n";
    std::ofstream of(fn_S);

    for(int imv=0; imv<n_mv; imv++){
      double mv = mvgrid.x(imv);
      of<<"\""<<std::fixed<<std::setprecision(2)<<mv*M_to_MeV<<" MeV\"   ";
      for(int i=0; i<num_bins; i++){
        double EaKev = (iEbin+i*wEbin)*E_to_keV;
        double EbKev = EaKev + wEbin*E_to_keV;
        of<<"\""
          <<std::fixed<<std::setprecision(1)<<EaKev<<"-"<<EbKev<<" keV\"   ";
      }
      of<<"\n"<<std::scientific<<std::setprecision(6);
      for(int imx=0; imx<n_mx; imx++){
        double mx = mxgrid.x(imx)*M_to_GeV;
        of<<mx<<" ";
        for(int ie=0; ie<num_bins; ie++){
          of<<S_mv_mx_E[imv][imx][ie]<<" ";
        }//mx
        of<<"\n";
      }//E
      of<<"\n";
    }//mv
    of.close();
  }

}


//******************************************************************************
//******************************************************************************
//******************************************************************************
int main(void){

  ChronoTimer sw(true);//start timer

  //define input parameters
  std::string akfn;
  int vsteps;
  double mxmin,mxmax,mvmin,mvmax; //m_chi and m_v masses
  int n_mx,i_mv,n_mv;
  std::string label;
  double Atot;
  double iEbin,fEbin,wEbin; // E bins: initial,final, width
  double dvesc,dv0;
  double dres,err_PEkeV; //detector resolution, PE-keV errors [-1]
  bool do_anMod;
  bool do_DAMA;
  bool write_dSdE, write_SofM;

  //Open and read the input file:
  {
    std::ifstream ifs("dmeXSection.in");
    int idodama, ianmod, iSM;
    std::string jnk;
    ifs >> akfn;                    getline(ifs,jnk);
    ifs >> label;                   getline(ifs,jnk);
    ifs >> vsteps;                  getline(ifs,jnk);
    ifs >> mxmin >> mxmax >> n_mx;  getline(ifs,jnk);
    ifs >> i_mv;                    getline(ifs,jnk);
    ifs >> mvmin >> mvmax >> n_mv;  getline(ifs,jnk);
    ifs >> ianmod >> dvesc >> dv0;  getline(ifs,jnk);
    ifs >> idodama;                 getline(ifs,jnk);
    ifs >> dres >> err_PEkeV;       getline(ifs,jnk);
    ifs >> Atot;                    getline(ifs,jnk);
    ifs >> iEbin >> fEbin >> wEbin; getline(ifs,jnk);
    ifs >> iSM;                     getline(ifs,jnk);
    label = label=="na" ? akfn : akfn+"-"+label;
    do_anMod = ianmod==1 ? true : false;
    do_DAMA  = idodama==1 ? true : false;
    write_SofM = iSM==0 ? false : true;
    write_dSdE = iSM==1 ? false : true;
  }

  if(!do_DAMA) do_anMod = false; //don't do anmod if outputting dsvde!
  if(do_anMod){
    std::cout<<"Running for annual modulation\n";
  }

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
    //Do for a range of m_v's (almost never do this..)
    mvmin /= M_to_MeV;
    mvmax /= M_to_MeV;
  }
  ExpGrid mvgrid(n_mv,mvmin,mvmax);

  //Energy bins for integrating/averaging
  iEbin/=E_to_keV;
  fEbin/=E_to_keV;
  wEbin/=E_to_keV;

  //Arrays/values to be filled from input AK file:
  FloatVec3D Kenq;
  std::vector<std::string> nklst;
  double qmin,qmax,demin,demax;
  //Read in AK file
  std::cout<<"Opening file: "<<akfn<<".bin\n";
  int iok=AKF::akReadWrite(akfn,false,Kenq,nklst,qmin,qmax,demin,demax);
  if(iok!=0) return 1;
  int desteps    = (int) Kenq.size();
  int num_states = (int) Kenq[0].size();
  int qsteps     = (int) Kenq[0][0].size();
  if(num_states != (int)nklst.size()) return 1; //just sanity check
  //Create the E and q grids:
  ExpGrid Egrid(desteps,demin,demax);
  ExpGrid qgrid(qsteps,qmin,qmax);

  // Grid of SHM vel. distro f_v(v). Can use to change vel profiles
  // Note: SHM is in km/s units, both for v and f!
  // f.dv = 1 => [f] = [1/v]
  double max_v = (SHMCONSTS::MAXV)/V_to_kms;
  double dv = max_v/vsteps;
  int num_cp = do_anMod ? 3 : 1; //just 1 v. dist? or 3 (for an mod.)?
  std::vector< std::vector<double> > arr_fv(num_cp,
      std::vector<double>(vsteps));
  for(int icp = 0; icp<num_cp; icp++){
    double cosp = icp==0 ? 0 : pow(-1,icp); //0, -1, 1
    StandardHaloModel shm(cosp,dvesc,dv0);
    for(int i=0; i<vsteps; i++){
      double v = (i+1)*dv; //don't include zero
      double vkms = v*(FPC::c_SI/FPC::c)/1.e3; //will be in km/s
      arr_fv[icp][i] = shm.fv(vkms)/(1./V_to_kms); //convert to a.u. [f]=[1/v]
    }
  }

  //Print the grid info to screen:
  std::cout<<"\nGrids:\n";
  printf("q :%6.2f -> %6.2f MeV, N=%4i\n",qmin*Q_to_MeV,qmax*Q_to_MeV,qsteps);
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
  std::cout<<"\nDoing calculations with sig-bar_e =  "<<sbe_1e37_cm2<<"cm2\n"
    <<"ds/dE conversion factor:   "<<dsdE_to_cm2keV<<" cm^2/keV\n"
    <<"ds.v/dE conversion factor: "
    <<dsvdE_to_cm3keVday<<"   cm^3/keV/day\n\n";

  double al_x = (sqrt(sbe_1e37_cm2)/FPC::aB_cm)*(FPC::alpha/sqrt(16.*M_PI));
  double al_xH = al_x/(FPC::m_e_MeV*FPC::alpha2);
  std::cout<<"Note: sig-bar_e =  "<<sbe_1e37_cm2<<"cm2 corresponds to:\n"
    <<" Light mediator: al_x = "<<al_x<<"\n"
    <<" Heavy mediator: al_x = "<<al_xH<<"*(mv/MeV)^2\n\n";


  // ********************************************************
  // Calculate <ds.v>/dE for each E (for each mx, mv)
  //ds.v/dE (fun. of mv, mx, E)
  FloatVec3D dsv_mv_mx_E; //Array to store cross-section
  sw.start();
  calculate_dsvde_array(Kenq,dsv_mv_mx_E,mvgrid,mxgrid,Egrid,qgrid,
    arr_fv,dv,do_anMod);
  std::cout<<"<ds.v>/dE: "<<sw.lap_reading_str()<<"\n";



  //Output <ds.v>/dE for gnuplot:
  bool write_dsvde = do_DAMA? false : true; //otherwise, get's printed too often
  if(write_dsvde){
    std::string fn_dsvde = "dsvde_mx-"+label+".out";
    std::cout<<"Writing to file: "<<fn_dsvde<<"\n";
    double u = dsvdE_to_cm3keVday; //convert units for <ds.v>/dE
    if(n_mv==1 || n_mx>1){
      writeForGnuplot_mvBlock(dsv_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dsvde,u);
    }
    if(n_mv>1){
      fn_dsvde = "dsvde_mv-"+label+".out";
      std::cout<<"Writing to file: "<<fn_dsvde<<"\n";
      writeForGnuplot_mvBlock(dsv_mv_mx_E,mvgrid,mxgrid,Egrid,fn_dsvde,u);
    }
  }

  // *********************
  //    Here, is just for DAMA! Move into sepperate function!!
  // *********************


  //Calculate and output obervable rate for DAMA
  if(do_DAMA)
    doDAMA(dsv_mv_mx_E,mvgrid, mxgrid, Egrid, do_anMod, write_dSdE, write_SofM,
    dres, err_PEkeV, Atot, iEbin, fEbin, wEbin, label);


  std::cout<<"\nTotal: "<<sw.reading_str()<<"\n";
  return 0;
}
