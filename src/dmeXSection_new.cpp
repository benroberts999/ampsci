#include "AKF_akFunctions.h"
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
};
void ExpGrid::ExpGrid(int in_N, double in_min, double in_max){
  N = in_N;
  min = in_min;
  max = in_max;
  dxonx = log(max/min)/(N-1);
}
double ExpGrid::x(int i){
  double y = double(i)/(N-1);
  return min*pow(max/min,y);
}
double ExpGrid::dxdi(int i){
  return dxonx*x(int i);
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
note: mu = mv c / hbar = mv/alpha!
Output in units of sig-bar_e
*/
{
  double arg = pow(mx*v,2)-2.*mx*E;
  if(arg<0) return 0;

  int num_states = (int) K_enq.size();
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
    #pragma omp parallel for
    for(int iq=0; iq<qsteps; iq++){
      double x = double(iq)/(qsteps-1);
      double q = qgrid.x(i);
      if(q<qminus || q>qplus) continue;
      double qdq_on_dqonq = q*q; //(dq/q) is constant, multiply at end
      double F_chi = 1.;
      if(finite_med) F_chi = 1./pow(q*q+mu*mu,2);
      #pragma omp critical (qint)
      {
        dsdE += qdq_on_dqonq*F_chi*K_enq[iE][ink][iq];  //dq/q included below
      }
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
Note: only takes _part_ of the K array! (for given E)
*/
{
  int vsteps = (int) arr_fv.size();
  double vmin = sqrt(2*E/mx);
  double dsvdE = 0;
  for(int iv=0; iv<vsteps; iv++){
    double v = (iv+1)*dv;
    if(v<vmin) continue;
    double dsdE = dsdE_iEdEvum_qg(Ke_nq,E,v,mv,mx,qmin,qgrid);
    dsvdE += arr_fv[iv]*v*dsdE;
  }//v
  return dsvdE*dv;
}
