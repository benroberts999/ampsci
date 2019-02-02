#include "AKF_akFunctions.h"
#include "ChronoTimer.h"
#include "ExponentialGrid.h"
#include "FPC_physicalConstants.h"
#include "FileIO_fileReadWrite.h"
#include "StandardHaloModel.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Convert FROM a.u. (hb=e=me=1, c=1/alpha) TO relativistic (hb=c=1)
double E_to_keV = FPC::Hartree_eV / 1000.;         // au -> keV (energy)
double Q_to_MeV = FPC::Hartree_eV * FPC::c / 1.e6; // momentum (q transfer)
double M_to_GeV = FPC::m_e_MeV / 1000.;
double M_to_MeV = FPC::m_e_MeV;
// Convert velocity FROM au to SI units
double V_to_kms = (FPC::c_SI / FPC::c) / 1000.;
double V_to_cms = (FPC::c_SI * FPC::alpha) * 100.; // au -> cm/s
double V_to_cmday = V_to_cms * (24 * 60 * 60);     // cm/s -> cm/day

// DM energy density (in GeV/cm^3)
double rhoDM_GeVcm3 = 0.4; // GeV/cm^3

// Default value used for bar-sigma_e + conversions from a.u.
double sbe_1e37_cm2 = 1.e-37;
double dsdE_to_cm2keV = sbe_1e37_cm2 / E_to_keV; // au -> cm^2/keV
double dsvdE_to_cm3keVday = dsdE_to_cm2keV * V_to_cmday;

// *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### *** ### ***
// Some macro/typedef short cuts
// Just to save time/space
using FloatVec3D = std::vector<std::vector<std::vector<float>>>;
using FloatVec2D = std::vector<std::vector<float>>;

//******************************************************************************
double g(double s, double x)
// Simple Gaussian
{
  double a = 0.398942 / s;
  double y = (x / s) * (x / s);
  return a * exp(-0.5 * y);
}

//******************************************************************************
/*
DM Form factors, simplified for ultra-light, and ultra-heavy cases
*/
double F_chi_2_heavy(double, double)
// Limit of heavy mediatior (mu >> q)
{
  return 1;
}
double F_chi_2_light(double, double q)
// Limit of light mediatior (mu << q)
{
  return 1. / (q * q * q * q);
}
double F_chi_2_intermediate(double mu, double q) {
  return pow((mu * mu + 1.) / (mu * mu + q * q), 2);
}

//******************************************************************************
void writeForGnuplot_mvBlock(const FloatVec3D &X_mv_mx_x, ExpGrid &mvgrid,
                             ExpGrid &mxgrid, ExpGrid &Egrid, std::string fname,
                             double y_unit_convert)
// Write to gnu-plot formatted text file.
// Each column is a m_chi (DM mass). Each block is different m_v
{
  std::ofstream of(fname.c_str());

  int n_mv = mvgrid.N();
  int n_mx = mxgrid.N();
  int desteps = Egrid.N();

  of << "# m_v blocks: ";
  for (int imv = 0; imv < n_mv; imv++) {
    double mv = mvgrid.x(imv);
    if (mv >= 0)
      of << imv << "," << mv * M_to_MeV << " ";
    else
      of << imv << ","
         << "Ultra-massive"
         << " ";
  }
  of << "\n";

  for (int imv = 0; imv < n_mv; imv++) {
    double mv = mvgrid.x(imv);
    if (mv >= 0)
      of << "\"" << std::fixed << std::setprecision(2) << mv * M_to_MeV
         << " MeV\"   ";
    else
      of << "\""
         << "infty"
         << " MeV\"   ";
    for (int imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.x(imx);
      of << "\"" << std::setprecision(1) << mx * M_to_GeV << " GeV\"   ";
    }
    of << "\n" << std::scientific << std::setprecision(6);
    for (int ie = 0; ie < desteps; ie++) {
      double E = Egrid.x(ie);
      of << E * E_to_keV << " ";
      for (int imx = 0; imx < n_mx; imx++) {
        of << double(X_mv_mx_x[imv][imx][ie]) * y_unit_convert << " ";
      } // mx
      of << "\n";
    } // E
    of << "\n";
  } // mv

  of.close();
}
//******************************************************************************
void writeForGnuplot_mxBlock(const FloatVec3D &X_mv_mx_x, ExpGrid &mvgrid,
                             ExpGrid &mxgrid, ExpGrid &Egrid,
                             const std::string &fname, double y_unit_convert)
// Write to gnu-plot formatted text file
// Each column is a m_v (mediator mass). Each block is different m_chi
{
  std::ofstream of(fname.c_str());

  int n_mv = mvgrid.N();
  int n_mx = mxgrid.N();
  int desteps = Egrid.N();

  of << "# m_x blocks: ";
  for (int imx = 0; imx < n_mx; imx++) {
    double mx = mxgrid.x(imx);
    of << imx << "," << mx * M_to_GeV << " ";
  }
  of << "\n";

  for (int imx = 0; imx < n_mx; imx++) {
    double mx = mxgrid.x(imx);
    of << "\"" << std::fixed << std::setprecision(2) << mx * M_to_GeV
       << " GeV\"   ";
    for (int imv = 0; imv < n_mv; imv++) {
      double mv = mvgrid.x(imv);
      of << "\"" << std::setprecision(2) << mv * M_to_MeV << " MeV\"   ";
    }
    of << "\n" << std::scientific << std::setprecision(6);
    for (int ie = 0; ie < desteps; ie++) {
      double E = Egrid.x(ie);
      of << E * E_to_keV << " ";
      for (int imv = 0; imv < n_mv; imv++) {
        of << double(X_mv_mx_x[imv][imx][ie]) * y_unit_convert << " ";
      } // mv
      of << "\n";
    } // E
    of << "\n";
  } // mx

  of.close();
}

//******************************************************************************
double dsdE_Evmvmx(const std::vector<std::vector<float>> &Ke_nq, double E,
                   double v, double mv, double mx, ExpGrid &qgrid,
                   double (*F_chi_2)(double, double))
/*
calcualtes cross-section ds/dE, for given E, v, mu, and mx.
Note: just takes in _PART_ of K (for given E)
Does q integrations, and sums over states
note: mu = mv c / hbar = mv/alpha!
If given negative value for mv, will use super-heavy mediator case
Output in units of sig-bar_e
Uses a function pointer for DM form factor. F_chi_2(mu,q) := |F_chi|^2
*/
{
  double arg = pow(mx * v, 2) - 2. * mx * E;
  if (arg < 0)
    return 0;
  int num_states = (int)(Ke_nq.size());
  int qsteps = qgrid.N();
  double dqonq = qgrid.dxonx();

  double mu = mv * FPC::c; // mu = m_v*c

  double qminus = mx * v - sqrt(arg);
  double qplus = mx * v + sqrt(arg);
  if (qminus > qgrid.max() || qplus < qgrid.min())
    return 0;

  double dsdE = 0;
  for (int ink = 0; ink < num_states; ink++) {
    for (int iq = 0; iq < qsteps; iq++) {
      double q = qgrid.x(iq);
      if (q < qminus || q > qplus)
        continue;
      double qdq_on_dqonq = q * q; //(dq/q) is constant, multiply at end
      double FX2 = F_chi_2(mu, q); // DM form factr (^2) [uses function pointer]
      dsdE += qdq_on_dqonq * FX2 * double(Ke_nq[ink][iq]);
      // dq/q included below
    } // q int
  }   // states

  double A = 0.5; // in units of sig-bar_e !
  dsdE *= (A / pow(v, 2)) * dqonq;
  return dsdE;
}

//******************************************************************************
double dsvdE_Evmvmx(const std::vector<std::vector<float>> &Ke_nq, double E,
                    double mv, double mx, ExpGrid &qgrid,
                    const std::vector<double> &arr_fv, double dv,
                    double (*F_chi_2)(double, double))
/*
Calculates <ds.v>/dE for given E, mv, mx
Does the v integration
Note: only takes _part_ of the K array! (for given E)
*/
{
  int vsteps = (int)arr_fv.size();
  double vmin = sqrt(2 * E / mx);
  double dsvdE = 0;
  for (int iv = 0; iv < vsteps; iv++) {
    double v = (iv + 1) * dv;
    if (v < vmin)
      continue;
    double dsdE = dsdE_Evmvmx(Ke_nq, E, v, mv, mx, qgrid, F_chi_2);
    dsvdE += arr_fv[iv] * v * dsdE;
  } // v
  return dsvdE * dv;
}

//******************************************************************************
void form_dsvdE(std::vector<float> &dsvde,
                const std::vector<std::vector<std::vector<float>>> &K_enq,
                double mv, double mx, ExpGrid &Egrid, ExpGrid &qgrid,
                const std::vector<double> &arr_fv, double dv,
                double (*F_chi_2)(double, double))
/*
Forms <ds.v>/dE array (for given mx, mv)
Note: mv<0 means "heavy" mediator [Fx=1]
*/
{
  // Loop through E, create dsvde array
  int desteps = Egrid.N();
  //#pragma omp parallel for
  for (int ie = 0; ie < desteps; ie++) {
    double E = Egrid.x(ie);
    // Do v (and q) integrations:
    double dsvdE =
        dsvdE_Evmvmx(K_enq[ie], E, mv, mx, qgrid, arr_fv, dv, F_chi_2);
    dsvde[ie] = (float)dsvdE;
  } // dE
}

//******************************************************************************
void calculate_dsvde_array(
    const std::vector<std::vector<std::vector<float>>> &Kenq,
    FloatVec3D &dsv_mv_mx_E, ExpGrid &mvgrid, ExpGrid &mxgrid, ExpGrid &Egrid,
    ExpGrid &qgrid, const std::vector<std::vector<double>> &arr_fv, double dv,
    bool do_anMod, double (*F_chi_2)(double, double))
/*
Fills the <ds.v>/dE array for each value of m_chi and m_v
Note: If doing annual modulation, then will calculate:
  <ds.v>_mod = (<ds.v>_max - <ds.v>_min)/2
instead
*/
{
  int n_mv = mvgrid.N();
  int n_mx = mxgrid.N();
  int desteps = Egrid.N();

  dsv_mv_mx_E.resize(
      n_mv, std::vector<std::vector<float>>(n_mx, std::vector<float>(desteps)));

  FloatVec3D dsv_mv_mx_Emax;
  if (do_anMod) {
    dsv_mv_mx_Emax.resize(n_mv, std::vector<std::vector<float>>(
                                    n_mx, std::vector<float>(desteps)));
  }

  // Calculate <ds.v>/dE for each E (for each mx, mv)
  std::cout << "Calculating <ds.v>/dE (doing q and v integrations): \n";
  for (int imv = 0; imv < n_mv; imv++) {
    double mv = mvgrid.x(imv);
#pragma omp parallel for
    for (int imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.x(imx);
      std::cout << "\r .. " << std::flush;
      if (do_anMod) {
        form_dsvdE(dsv_mv_mx_E[imv][imx], Kenq, mv, mx, Egrid, qgrid, arr_fv[1],
                   dv, F_chi_2);
        std::cout << "\r  ... " << std::flush;
        form_dsvdE(dsv_mv_mx_Emax[imv][imx], Kenq, mv, mx, Egrid, qgrid,
                   arr_fv[2], dv, F_chi_2);
      } else {
        form_dsvdE(dsv_mv_mx_E[imv][imx], Kenq, mv, mx, Egrid, qgrid, arr_fv[0],
                   dv, F_chi_2);
      }
      std::cout << "\r .. .. " << std::flush;
    } // mx
  }   // mv
  std::cout << "\n";

  // If doing an. mod, form "amplitude" for ds.v/dE, used below
  if (do_anMod) {
    std::cout << "\nAnnual modulation.\n";
    std::cout << "Forming <ds.v>_mod = (<ds.v>_max - <ds.v>_min)/2 \n";
    for (int imv = 0; imv < n_mv; imv++) {
      for (int imx = 0; imx < n_mx; imx++) {
        for (int ie = 0; ie < desteps; ie++) {
          dsv_mv_mx_E[imv][imx][ie] =
              0.5f * (dsv_mv_mx_Emax[imv][imx][ie] - dsv_mv_mx_E[imv][imx][ie]);
        }
      }
    }
  }
}

//******************************************************************************
void doDAMA(const FloatVec3D &dsv_mv_mx_E, ExpGrid &mvgrid, ExpGrid &mxgrid,
            ExpGrid &Egrid, bool do_anMod, bool write_dSdE, bool write_SofM,
            double dres, double err_PEkeV, double Atot, double iEbin,
            double fEbin, double wEbin, std::string label)
/*
Takes in <ds.v>/dE, and calculated dS/dE, for DAMA.
Depending on above, will be either S_0 or S_m (mod. amplitude)
Includes detector resolution etc.
Optionally further integrates into energy bins
*/
{
  std::cout << "\n*** Doing S1 for DAMA ***\n";
  int n_mv = mvgrid.N();
  int n_mx = mxgrid.N();
  int desteps = Egrid.N();

  // Array to store observable Rate, S
  FloatVec3D dSdE_mv_mx_E;
  dSdE_mv_mx_E.resize(
      n_mv, std::vector<std::vector<float>>(n_mx, std::vector<float>(desteps)));

  // Hardware threshold: should be between [-1,1]
  double PE_per_keV = 6.5 + err_PEkeV * 1.;
  double E_thresh_HW = 1. / PE_per_keV / E_to_keV;

  // DAMA parameters for Gaussian resolution (smearing)
  double alpha = 0.45 + dres * 0.04;
  double beta = 0.009 + dres * 0.005;

  double dEonE = Egrid.dxonx();
  double MN =
      Atot * (FPC::u_NMU * FPC::m_e_kg); // Total atomic/mol. mass (in kg)

  // Calculate _observable_ rate, S, for DAMA
  // inlcuding Gaussian resolution, +
  for (int imv = 0; imv < n_mv; imv++) {
    for (int imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.x(imx);
      double rho_on_mxc2 = rhoDM_GeVcm3 / (mx * M_to_GeV);
      double rateFac = dsvdE_to_cm3keVday * rho_on_mxc2 / MN;
      for (int i = 0; i < desteps; i++) {
        double E = Egrid.x(i);
        if (E < E_thresh_HW) {
          dSdE_mv_mx_E[imv][imx][i] = 0;
          continue;
        }
        // triple check this! Super important!
        double s =
            (alpha * sqrt(E * E_to_keV) + beta * (E * E_to_keV)) / E_to_keV;
        double y0 = 0;
        // integrate over Eprime:
        for (int j = 0; j < desteps; j++) {
          double Ep = Egrid.x(j);
          if (Ep < E_thresh_HW)
            continue; // hardware threshold
          y0 += g(s, E - Ep) * double(dsv_mv_mx_E[imv][imx][j]) * Ep;
        }
        dSdE_mv_mx_E[imv][imx][i] = (float)(y0 * dEonE * rateFac);
      }
    }
  }

  //(Only write dS/dE if not writing )
  // output dS/dE for gnuplot:
  if (write_dSdE) {
    std::string spref = do_anMod ? "dSmdE" : "dSdE";
    std::string fn_dSde = spref + "_mx-" + label + ".out";
    std::cout << "Writing to file: " << fn_dSde << "\n";
    double u = 1; // already converted!
    if (n_mv == 1 || n_mx > 1) {
      writeForGnuplot_mvBlock(dSdE_mv_mx_E, mvgrid, mxgrid, Egrid, fn_dSde, u);
    }
    if (n_mv > 1) {
      fn_dSde = spref + "_mv-" + label + ".out";
      std::cout << "Writing to file: " << fn_dSde << "\n";
      writeForGnuplot_mxBlock(dSdE_mv_mx_E, mvgrid, mxgrid, Egrid, fn_dSde, u);
    }
  }

  if (wEbin == 0 || fEbin == 0 || iEbin == fEbin)
    return;
  // Integrate/average into energy bins:
  int num_bins = (int)ceil((fEbin - iEbin) / wEbin);

  // Array to store averaged Rate, S
  FloatVec3D S_mv_mx_E;
  S_mv_mx_E.resize(n_mv, std::vector<std::vector<float>>(
                             n_mx, std::vector<float>(num_bins)));

  // Creates S (or Sm) [actually S/E_w] - rate averaged over each energy bin
  // Also prints energy bins to screen (only for first mx, mv)
  std::cout << "\nCalculating ";
  if (do_anMod)
    std::cout << "Sm/Ew, annual modulation amplitude (/day/kg/keV)\n";
  else
    std::cout << "S/Ew, average event rate (/day/kg/keV)\n";
  std::cout << "Integrated (averaged) each energy bin.\n";
  std::cout << "(Just outputting for first mx/mv: ";
  printf("M_chi=%5.2f GeV", mxgrid.min() * M_to_GeV);
  if (mvgrid.min() >= 0)
    printf(" ; M_v=%6.3f MeV", mvgrid.min() * M_to_MeV);
  std::cout << ")\n";
  for (int imv = 0; imv < n_mv; imv++) {
    for (int imx = 0; imx < n_mx; imx++) {
      for (int i = 0; i < num_bins; i++) {
        int ieA = Egrid.findNextIndex(iEbin + i * wEbin);
        int ieB = Egrid.findNextIndex(iEbin + (i + 1) * wEbin);
        if (ieA < 0)
          ieA = 0;

        double Rate = 0;
        for (int ie = ieA; ie < ieB; ie++) {
          if (ieB >= desteps)
            break;
          double E = Egrid.x(ie);
          Rate += double(dSdE_mv_mx_E[imv][imx][ie]) * E;
          // nb: E is from Jacobian; * dE/E below
        }
        Rate *= dEonE / wEbin;
        S_mv_mx_E[imv][imx][i] = (float)Rate;
        double EaKev = Egrid.x(ieA) * E_to_keV;
        double EbKev = Egrid.x(ieB) * E_to_keV;
        if (imv == 0 && imx == 0) {
          // only print first one to screen
          printf("%3.1f-%3.1f: %6.3f   %.2e     %i\n", EaKev, EbKev,
                 0.5 * (EaKev + EbKev), Rate, ieB - ieA);
        }
      }
    }
  }

  if (!write_SofM)
    return;
  // For each energy bin, output as function of m_chi.
  // Note: if doing this, don't output above files! (too many m_chi's)
  // Each collumn a different energy bin
  // Each m_v
  std::string fn_S = do_anMod ? "Sm" : "S";
  fn_S = fn_S + "_mx-" + label + ".out";
  std::ofstream of; //(fn_S);

  if (n_mv == 1 || n_mx > 1) {
    // Write Sm/[E] as function of m_chi, each m_v as block
    std::cout << "Writing to file: " << fn_S << "\n";
    of.open(fn_S);
    for (int imv = 0; imv < n_mv; imv++) {
      double mv = mvgrid.x(imv);
      if (mv >= 0)
        of << "\"" << std::fixed << std::setprecision(2) << mv * M_to_MeV
           << " MeV\"   ";
      else
        of << "\""
           << "infty"
           << " MeV\"   ";
      for (int i = 0; i < num_bins; i++) {
        double EaKev = (iEbin + i * wEbin) * E_to_keV;
        double EbKev = EaKev + wEbin * E_to_keV;
        of << "\"" << std::fixed << std::setprecision(1) << EaKev << "-"
           << EbKev << " keV\"   ";
      }
      of << "\n" << std::scientific << std::setprecision(6);
      for (int imx = 0; imx < n_mx; imx++) {
        double mx = mxgrid.x(imx) * M_to_GeV;
        of << mx << " ";
        for (int ie = 0; ie < num_bins; ie++) {
          of << S_mv_mx_E[imv][imx][ie] << " ";
        } // mx
        of << "\n";
      } // E
      of << "\n";
    } // mv
    of.close();
  }

  if (n_mv > 1) {
    // Write Sm/[E] as function of m_v, each m_chi as block
    fn_S = do_anMod ? "Sm" : "S";
    fn_S = fn_S + "_mv-" + label + ".out";
    std::cout << "Writing to file: " << fn_S << "\n";
    of.open(fn_S);
    for (int imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.x(imx);
      of << "\"" << std::fixed << std::setprecision(2) << mx * M_to_GeV
         << " GeV\"   ";
      for (int i = 0; i < num_bins; i++) {
        double EaKev = (iEbin + i * wEbin) * E_to_keV;
        double EbKev = EaKev + wEbin * E_to_keV;
        of << "\"" << std::fixed << std::setprecision(1) << EaKev << "-"
           << EbKev << " keV\"   ";
      }
      of << "\n" << std::scientific << std::setprecision(6);
      for (int imv = 0; imv < n_mv; imv++) {
        double mv = mvgrid.x(imv) * M_to_MeV;
        of << mv << " ";
        for (int ie = 0; ie < num_bins; ie++) {
          of << S_mv_mx_E[imv][imx][ie] << " ";
        } // mv
        of << "\n";
      } // E
      of << "\n";
    } // mx
    of.close();
  }
}

//******************************************************************************
//******************************************************************************
//******************************************************************************
int main() {
  ChronoTimer sw(true); // start timer

  // define input parameters
  std::string akfn;
  double mxmin, mxmax, mvmin, mvmax; // m_chi and m_v masses
  int n_mx, n_mv;                    // number of steps in m_chi, and m_v grids
  std::string label;                 // label to append to output file names
  double Atot;                       // Total nuclear mass number (Na+I)
  double iEbin, fEbin, wEbin;        // E bins: initial,final, width
  double dvesc, dv0;                 // error terms for SHM
  double dres, err_PEkeV;            // detector resolution, PE-keV errors [-1]
  // Bools to control which parts to calculate:
  bool do_anMod;
  bool do_DAMA;
  bool write_dsvde, write_dSdE, write_SofM;
  // Which type of mediator (for DM form factor)
  enum class Mediator { light, intermediate, heavy };
  Mediator mediator;

  // Open and read the input file:
  {
    int i_mv, iwr_dsvde, idodama, ianmod, iSM; // temp settings
    auto tp = std::forward_as_tuple(akfn, label, mxmin, mxmax, n_mx, i_mv,
                                    mvmin, mvmax, n_mv, dvesc, dv0, ianmod,
                                    iwr_dsvde, idodama, dres, err_PEkeV, Atot,
                                    iEbin, fEbin, wEbin, iSM);
    FileIO::setInputParameters("dmeXSection.in", tp);
    label = label == "na" ? akfn : akfn + "-" + label;
    // what to write to file:
    write_dsvde = iwr_dsvde == 1 ? true : false;
    write_SofM = iSM == 0 ? false : true;
    write_dSdE = iSM == 1 ? false : true;
    // What to calclate:
    do_anMod = ianmod == 1 ? true : false;
    do_DAMA = idodama == 1 ? true : false;
    switch (i_mv) {
    case 0:
      mediator = Mediator::light;
      break;
    case 1:
      mediator = Mediator::intermediate;
      break;
    case 2:
      mediator = Mediator::heavy;
      break;
    default:
      std::cout << "Invalid mediator\n";
      return 1;
    }
  }

  std::cout << "XXX"
            << "\n"
            << std::flush;

  // Number of points in the velocity integration.
  const int vsteps = 100; // no need to be input

  // DM mass: Convert from GeV to au + define grids
  mxmin /= M_to_GeV;
  mxmax /= M_to_GeV;
  if (n_mx == 1)
    mxmax = mxmin;
  ExpGrid mxgrid(n_mx, mxmin, mxmax);

  // Append 'l'('h') to label for light(heavy) mediator
  switch (mediator) {
  case Mediator::light:
    label = label + "_l";
    break;
  case Mediator::intermediate:
    label = label + "_m";
    break;
  case Mediator::heavy:
    label = label + "_h";
    break;
  }

  // Mediator mass: convert units + set-up finite/infinite case
  // Note: For super-heavy mediator (contact interaction)
  // will set m_v to be negative! (this is silly, but it works)
  if (mediator == Mediator::light) {
    // massless case. Do sepperately
    mvmin = mvmax = 0;
    n_mv = 1;
  } else if (mediator == Mediator::heavy) {
    // Heavy-mediator case (contact interaction)
    mvmin = mvmax = -1.; // switch..silly..
    n_mv = 1;
  } else {
    // Do for a range of m_v's
    if (mvmin == 0 || mvmax == 0) {
      std::cout << "cannot have 0 mv mass for interdiatie mediator\n";
      return 1;
    }
    mvmin /= M_to_MeV;
    mvmax /= M_to_MeV;
  }
  ExpGrid mvgrid(n_mv, mvmin, mvmax);

  // Energy bins for integrating/averaging (convert units)
  iEbin /= E_to_keV;
  fEbin /= E_to_keV;
  wEbin /= E_to_keV;

  // Arrays/values to be filled from input AK file:
  FloatVec3D Kenq;
  std::vector<std::string> nklst;
  double qmin, qmax, demin, demax;
  // Read in AK file
  std::cout << "Opening file: " << akfn << ".bin\n";
  int iok =
      AKF::akReadWrite(akfn, false, Kenq, nklst, qmin, qmax, demin, demax);
  if (iok != 0)
    return 1;
  int desteps = (int)Kenq.size();
  int num_states = (int)Kenq[0].size();
  int qsteps = (int)Kenq[0][0].size();
  if (num_states != (int)nklst.size())
    return 1; // just sanity check
  // Create the E and q grids:
  ExpGrid Egrid(desteps, demin, demax);
  ExpGrid qgrid(qsteps, qmin, qmax);

  // Grid of SHM vel. distro f_v(v). Can use to change vel profiles
  // Note: SHM is in km/s units, both for v and f!
  // f.dv = 1 => [f] = [1/v]
  double max_v = (SHMCONSTS::MAXV) / V_to_kms;
  double dv = max_v / vsteps;
  int num_cp = do_anMod ? 3 : 1; // just 1 v. dist? or 3 (for an mod.)?
  std::vector<std::vector<double>> arr_fv(num_cp, std::vector<double>(vsteps));
  for (int icp = 0; icp < num_cp; icp++) {
    double cosp = icp == 0 ? 0 : pow(-1, icp); // 0, -1, 1
    StandardHaloModel shm(cosp, dvesc, dv0);
    for (int i = 0; i < vsteps; i++) {
      double v = (i + 1) * dv;                       // don't include zero
      double vkms = v * (FPC::c_SI / FPC::c) / 1.e3; // will be in km/s
      arr_fv[icp][i] =
          shm.fv(vkms) / (1. / V_to_kms); // convert to a.u. [f]=[1/v]
    }
  }

  // Print the grid info to screen:
  if (do_anMod)
    std::cout << "Running for annual modulation amplitude\n";
  std::cout << "\nGrids:\n";
  printf("q :%6.2f -> %6.2f MeV, N=%4i\n", qmin * Q_to_MeV, qmax * Q_to_MeV,
         qsteps);
  printf("E :%6.2f -> %6.2f keV, N=%4i\n", demin * E_to_keV, demax * E_to_keV,
         desteps);
  printf("v :%6.2f -> %6.2fkm/s, N=%4i\n", dv * V_to_kms, max_v * V_to_kms,
         vsteps);
  printf("Mx:%6.2f -> %6.2f GeV, N=%4i\n", mxmin * M_to_GeV, mxmax * M_to_GeV,
         n_mx);
  switch (mediator) {
  case Mediator::light:
    std::cout << "Ultra-light meadiator (F_x = q^-2)\n";
    break;
  case Mediator::heavy:
    std::cout << "Heavy meadiator  (F_x = 1)\n";
    break;
  default:
    printf("Mv:%6.2f -> %6.2f MeV, N=%4i\n", mvmin * M_to_MeV, mvmax * M_to_MeV,
           n_mv);
  }

  // Units + conversions for dsvde..etc
  std::cout << "\nDoing calculations with sig-bar_e =  " << sbe_1e37_cm2
            << "cm2\n"
            << "ds/dE conversion factor:   " << dsdE_to_cm2keV << " cm^2/keV\n"
            << "ds.v/dE conversion factor: " << dsvdE_to_cm3keVday
            << "   cm^3/keV/day\n\n";

  double al_x =
      (sqrt(sbe_1e37_cm2) / FPC::aB_cm) * (FPC::alpha / sqrt(16. * M_PI));
  double al_xH = al_x / (FPC::m_e_MeV * FPC::alpha2);
  std::cout << "Note: sig-bar_e =  " << sbe_1e37_cm2 << "cm2 corresponds to:\n"
            << " Light mediator: al_x = " << al_x << "\n"
            << " Heavy mediator: al_x = " << al_xH << "*(mv/MeV)^2\n\n";

  // Define + set function pointer for DM form-factor:
  double (*F_chi_2)(double, double);
  switch (mediator) {
  case Mediator::light:
    F_chi_2 = &F_chi_2_light;
    break;
  case Mediator::heavy:
    F_chi_2 = &F_chi_2_heavy;
    break;
  case Mediator::intermediate:
    F_chi_2 = &F_chi_2_intermediate;
    break;
  }

  // ********************************************************
  // Calculate <ds.v>/dE for each E (for each mx, mv)
  // ds.v/dE (fun. of mv, mx, E)
  FloatVec3D dsv_mv_mx_E; // Array to store cross-section
  sw.start();
  calculate_dsvde_array(Kenq, dsv_mv_mx_E, mvgrid, mxgrid, Egrid, qgrid, arr_fv,
                        dv, do_anMod, F_chi_2);
  std::cout << "<ds.v>/dE: " << sw.lap_reading_str() << "\n";

  // Output <ds.v>/dE for gnuplot:
  if (write_dsvde) {
    std::string prefix = do_anMod ? "dsvde_mod" : "dsvde";
    std::string fn_dsvde = prefix + "_mx-" + label + ".out";
    std::cout << "Writing to file: " << fn_dsvde << "\n";
    double u = dsvdE_to_cm3keVday; // convert units for <ds.v>/dE
    if (n_mv == 1 || n_mx > 1) {
      writeForGnuplot_mvBlock(dsv_mv_mx_E, mvgrid, mxgrid, Egrid, fn_dsvde, u);
    }
    if (n_mv > 1) {
      fn_dsvde = prefix + "_mv-" + label + ".out";
      std::cout << "Writing to file: " << fn_dsvde << "\n";
      writeForGnuplot_mxBlock(dsv_mv_mx_E, mvgrid, mxgrid, Egrid, fn_dsvde, u);
    }
  }

  // ********************************************************
  // Calculate and output obervable rate for DAMA
  if (do_DAMA) {
    doDAMA(dsv_mv_mx_E, mvgrid, mxgrid, Egrid, do_anMod, write_dSdE, write_SofM,
           dres, err_PEkeV, Atot, iEbin, fEbin, wEbin, label);
  }

  std::cout << "\nTotal: " << sw.reading_str() << "\n";
  return 0;
}
