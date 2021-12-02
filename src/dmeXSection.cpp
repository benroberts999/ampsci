#include "DMionisation/AKF_akFunctions.hpp"
#include "DMionisation/StandardHaloModel.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#pragma GCC diagnostic ignored "-Wsign-conversion"

// Convert FROM a.u. (hb=e=me=1, c=1/alpha) TO relativistic (hb=c=1)
const double E_to_keV = PhysConst::Hartree_eV / 1000.; // au -> keV (energy)
const double Q_to_MeV =
    PhysConst::Hartree_eV * PhysConst::c / 1.e6; // momentum (q transfer)
const double M_to_GeV = PhysConst::m_e_MeV / 1000.;
const double M_to_MeV = PhysConst::m_e_MeV;
// Convert velocity FROM au to SI units
const double V_to_kms = (PhysConst::c_SI / PhysConst::c) / 1000.;
const double V_to_cms =
    (PhysConst::c_SI * PhysConst::alpha) * 100.;     // au -> cm/s
const double V_to_cmday = V_to_cms * (24 * 60 * 60); // cm/s -> cm/day

// DM energy density (in GeV/cm^3)
const double rhoDM_GeVcm3 = 0.4; // GeV/cm^3

// Default value used for bar-sigma_e + conversions from a.u.
const double sbe_1e37_cm2 = 1.e-37;
const double dsdE_to_cm2keV = sbe_1e37_cm2 / E_to_keV; // au -> cm^2/keV
const double dsvdE_to_cm3keVday = dsdE_to_cm2keV * V_to_cmday;

const double dsdE_to_cm2_perau = sbe_1e37_cm2; // au -> cm^2/keV
const double dsvdE_to_cm3_per_auday = dsdE_to_cm2_perau * V_to_cmday;

// Some typedef short cuts
// Just to save time/space
using FloatVec3D = std::vector<std::vector<std::vector<float>>>;
using DoubleVec3D = std::vector<std::vector<std::vector<double>>>;
using DoubleVec2D = std::vector<std::vector<double>>;

//******************************************************************************
double quickExp(double x) {
  // Good to parts in 1.e-5
  if (x < 0)
    return 1.0 / quickExp(-x);
  if (x < 0.05)
    return 1.0 + x + 0.5 * x * x;
  if (x < 0.25)
    return 1.0 + x + 0.5 * x * x + 0.16666667 * x * x * x +
           0.041666667 * x * x * x * x;
  return std::exp(x);
}

//******************************************************************************
double gaussian(double s, double x)
// Simple Gaussian
{
  double a = 0.398942 / s;
  double y = (x / s) * (x / s);
  return a * quickExp(-0.5 * y);
}

//******************************************************************************
struct SumLogk
/*
Calculates sum [log(n), n=1 -> k]
Uses a lookup-table. Every time given new k value, calcs up to that k.
If given a k value smaller than max already calc'd, just looks up answer.
Designed to be included as a static object inside some other function
*/
{

public:
  SumLogk(int k = 0) {
    sum_logk_array.reserve(30);
    if (k > 0)
      calculate_new_term(k);
  };

  double sum_logk(int k) {
    if (k <= 1)
      return 0;
    if (k <= (int)sum_logk_array.size())
      return sum_logk_array[k - 1];
    else
      return calculate_new_term(k);
  }

  double calculate_new_term(int k) {
    int max_k_sofar = (int)sum_logk_array.size();
    double sum = sum_logk_array.back();
    for (int ik = max_k_sofar + 1; ik <= k; ik++) {
      // std::cerr << "\n k=" << ik << "\n";
      sum += log(ik);
      sum_logk_array.push_back(sum);
    }
    return sum;
  }

private:
  std::vector<double> sum_logk_array;
};

//******************************************************************************
double Pois(int k, double lambda)
/*
Safely returns Poisson distribution, without danger of overflows.
Also, uses a static 'on-the-fly' lookup table for  Sum[log[n]]
P_k(lambda) = Exp[-lambda] * (lambda^k) / k!
log(P) = ln(p) = -lambda + k * log(lambda) - Sum[log[n],{n,1,k}]
k must be >= 0. No check is performed.
*/
{
  static SumLogk slk;
  double logP = -lambda + k * log(lambda) - slk.sum_logk(k);
  return quickExp(logP);
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
  // return 1. / std::pow(q, 4);
  return 1. / (q * q * q * q); // wow. MUCH faster? wtf compiler can't fix?
}
double F_chi_2_intermediate(double mu, double q) {
  double a = (mu * mu + 1.) / (mu * mu + q * q);
  return a * a;
}

//******************************************************************************
void writeForGnuplot_mvBlock(const DoubleVec3D &X_mv_mx_x, const Grid &mvgrid,
                             const Grid &mxgrid, const Grid &Egrid,
                             const std::string &fname, double y_unit_convert,
                             double x_unit_convert = E_to_keV)
// Write to gnu-plot formatted text file.
// Each column is a m_chi (DM mass). Each block is different m_v
{
  std::ofstream of(fname.c_str());
  of << "# sigma^bar_e = " << sbe_1e37_cm2 << " cm^2\n";

  auto n_mv = mvgrid.num_points();
  auto n_mx = mxgrid.num_points();
  auto desteps = Egrid.num_points(); //.num_points;

  of << "# m_v blocks: ";
  for (auto imv = 0ul; imv < n_mv; imv++) {
    double mv = mvgrid.r()[imv];
    if (mv >= 0)
      of << imv << "," << mv * M_to_MeV << " ";
    else
      of << imv << ","
         << "Ultra-massive"
         << " ";
  }
  of << "\n";

  for (std::size_t imv = 0; imv < n_mv; imv++) {
    double mv = mvgrid.r()[imv];
    if (mv >= 0)
      of << "\"" << std::fixed << std::setprecision(2) << mv * M_to_MeV
         << " MeV\"   ";
    else
      of << "\""
         << "infty"
         << " MeV\"   ";
    // for (int imx = 0; imx < n_mx; imx++) {
    //   double mx = mxgrid.r()[imx];
    //   of << "\"" << std::setprecision(1) << mx * M_to_GeV << " GeV\"   ";
    // }
    for (auto mx : mxgrid.r()) {
      of << "\"" << std::setprecision(1) << mx * M_to_GeV << " GeV\"   ";
    }
    of << "\n" << std::scientific << std::setprecision(6);
    for (std::size_t ie = 0; ie < desteps; ie++) {
      double E = Egrid.r()[ie]; // x(ie);
      of << E * x_unit_convert << " ";
      for (std::size_t imx = 0; imx < n_mx; imx++) {
        of << double(X_mv_mx_x[imv][imx][ie]) * y_unit_convert << " ";
      } // mx
      of << "\n";
    } // E
    of << "\n";
  } // mv

  of.close();
}
//******************************************************************************
void writeForGnuplot_mxBlock(const DoubleVec3D &X_mv_mx_x, const Grid &mvgrid,
                             const Grid &mxgrid, const Grid &Egrid,
                             const std::string &fname, double y_unit_convert,
                             double x_unit_convert = E_to_keV)
// Write to gnu-plot formatted text file
// Each column is a m_v (mediator mass). Each block is different m_chi
{
  std::ofstream of(fname.c_str());
  of << "# sigma^bar_e = " << sbe_1e37_cm2 << " cm^2\n";

  std::size_t n_mv = mvgrid.num_points();
  std::size_t n_mx = mxgrid.num_points();
  std::size_t desteps = Egrid.num_points(); // num_points;

  of << "# m_x blocks: ";
  for (std::size_t imx = 0; imx < n_mx; imx++) {
    double mx = mxgrid.r()[imx];
    of << imx << "," << mx * M_to_GeV << " ";
  }
  of << "\n";

  for (std::size_t imx = 0; imx < n_mx; imx++) {
    double mx = mxgrid.r()[imx];
    of << "\"" << std::fixed << std::setprecision(2) << mx * M_to_GeV
       << " GeV\"   ";
    for (std::size_t imv = 0; imv < n_mv; imv++) {
      double mv = mvgrid.r()[imv];
      of << "\"" << std::setprecision(2) << mv * M_to_MeV << " MeV\"   ";
    }
    of << "\n" << std::scientific << std::setprecision(6);
    for (std::size_t ie = 0; ie < desteps; ie++) {
      double E = Egrid.r()[ie]; // x(ie);
      of << E * x_unit_convert << " ";
      for (std::size_t imv = 0; imv < n_mv; imv++) {
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
                   double v, double mv, double mx, const Grid &qgrid,
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
  double arg = std::pow(mx * v, 2) - 2. * mx * E;
  if (arg < 0)
    return 0;
  std::size_t num_states = (Ke_nq.size());
  std::size_t qsteps = qgrid.num_points();

  double mu = mv * PhysConst::c;

  double qminus = mx * v - std::sqrt(arg);
  double qplus = mx * v + std::sqrt(arg);
  if (qminus > qgrid.rmax() || qplus < qgrid.r0())
    return 0;

  double dsdE = 0;
  for (std::size_t ink = 0; ink < num_states; ink++) {
    for (std::size_t iq = 0; iq < qsteps; iq++) {
      double q = qgrid.r()[iq];
      if (q < qminus || q > qplus)
        continue;
      double qdq_ondu = q * qgrid.drdu()[iq];
      double FX2 = F_chi_2(mu, q); // DM form factr (^2) [uses function pointer]
      dsdE += qdq_ondu * FX2 * double(Ke_nq[ink][iq]);
    } // q int
  }   // states

  double A = 0.5; // in units of sig-bar_e !
  dsdE *= (A / std::pow(v, 2)) * qgrid.du();
  return dsdE;
}

//******************************************************************************
double dsvdE_Evmvmx(const std::vector<std::vector<float>> &Ke_nq, double E,
                    double mv, double mx, const Grid &qgrid,
                    const std::vector<double> &arr_fv, double dv,
                    double (*F_chi_2)(double, double))
/*
Calculates <ds.v>/dE for given E, mv, mx
Does the v integration
Note: only takes _part_ of the K array! (for given E)
*/
{
  std::size_t vsteps = arr_fv.size();
  double vmin = std::sqrt(2 * E / mx);
  double dsvdE = 0;
  for (std::size_t iv = 0; iv < vsteps; iv++) {
    double v = double(iv + 1) * dv;
    if (v < vmin)
      continue;
    double dsdE = dsdE_Evmvmx(Ke_nq, E, v, mv, mx, qgrid, F_chi_2);
    dsvdE += arr_fv[iv] * v * dsdE;
  } // v
  return dsvdE * dv;
}

//******************************************************************************
void form_dsvdE(std::vector<double> &dsvde,
                const std::vector<std::vector<std::vector<float>>> &K_enq,
                double mv, double mx, const Grid &Egrid, const Grid &qgrid,
                const std::vector<double> &arr_fv, double dv,
                double (*F_chi_2)(double, double))
/*
Forms <ds.v>/dE array (for given mx, mv)
Note: mv<0 means "heavy" mediator [Fx=1]
*/
{
  // Loop through E, create dsvde array
  std::size_t desteps = Egrid.num_points();
  for (std::size_t ie = 0; ie < desteps; ie++) {
    double E = Egrid.r()[ie];
    // Do v (and q) integrations:
    double dsvdE =
        dsvdE_Evmvmx(K_enq[ie], E, mv, mx, qgrid, arr_fv, dv, F_chi_2);
    dsvde[ie] = dsvdE;
  } // dE
}

//******************************************************************************
void calculate_dsvde_array(
    const std::vector<std::vector<std::vector<float>>> &Kenq,
    DoubleVec3D &dsv_mv_mx_E, const Grid &mvgrid, const Grid &mxgrid,
    const Grid &Egrid, const Grid &qgrid,
    const std::vector<std::vector<double>> &arr_fv, double dv, bool do_anMod,
    double (*F_chi_2)(double, double))
/*
Fills the <ds.v>/dE array for each value of m_chi and m_v
Note: If doing annual modulation, then will calculate:
  <ds.v>_mod = (<ds.v>_max - <ds.v>_min)/2
instead
*/
{
  std::size_t n_mv = mvgrid.num_points();
  std::size_t n_mx = mxgrid.num_points();
  std::size_t desteps = Egrid.num_points();

  dsv_mv_mx_E.resize(n_mv, std::vector<std::vector<double>>(
                               n_mx, std::vector<double>(desteps)));

  DoubleVec3D dsv_mv_mx_Emax;
  if (do_anMod) {
    dsv_mv_mx_Emax.resize(n_mv, std::vector<std::vector<double>>(
                                    n_mx, std::vector<double>(desteps)));
  }

  // Calculate <ds.v>/dE for each E (for each mx, mv)
  std::cout << "Calculating <ds.v>/dE (doing q and v integrations): \n";
  for (std::size_t imv = 0; imv < n_mv; imv++) {
    double mv = mvgrid.r()[imv];
#pragma omp parallel for
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.r()[imx];
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
    for (std::size_t imv = 0; imv < n_mv; imv++) {
      for (std::size_t imx = 0; imx < n_mx; imx++) {
        for (std::size_t ie = 0; ie < desteps; ie++) {
          dsv_mv_mx_E[imv][imx][ie] =
              0.5 * (dsv_mv_mx_Emax[imv][imx][ie] - dsv_mv_mx_E[imv][imx][ie]);
        } // energy
      }   // mx
    }     // mv
  }       // end an. mod.
}

//******************************************************************************
std::vector<double>
convolvedRate(const std::vector<double> &in_rate, const Grid &in_grid,
              const std::vector<std::vector<double>> &f_conv,
              double convert_units) {

  auto Nout = f_conv.size();

  // std::size_t Nout = f_smear.size();
  std::vector<double> out_rate; //(Nout);
  out_rate.reserve(Nout);
  for (std::size_t i = 0; i < Nout; i++) {
    double g =
        // NumCalc::integrate({&in_rate, &f_conv[i], &in_grid.drdu()},
        // in_grid.du());
        NumCalc::integrate(in_grid.du(), 0, 0, in_rate, f_conv[i],
                           in_grid.drdu());
    out_rate.push_back(convert_units * g);
  }

  return out_rate;
}

//******************************************************************************
void doDAMA(const DoubleVec3D &dsv_mv_mx_E, const Grid &mvgrid,
            const Grid &mxgrid, const Grid &Egrid, bool do_anMod,
            bool write_dSdE, bool write_SofM, double dres, double err_PEkeV,
            double Atot, double iEbin, double fEbin, double wEbin,
            const std::string &label)
/*
Takes in <ds.v>/dE, and calculated dS/dE, for DAMA.
Depending on above, will be either S_0 or S_m (mod. amplitude)
Includes detector resolution etc.
Optionally further integrates into energy bins
*/
{
  std::cout << "\n*** Doing S1 for DAMA ***\n";
  std::size_t n_mv = mvgrid.num_points();
  std::size_t n_mx = mxgrid.num_points();
  std::size_t desteps = Egrid.num_points();

  // DAMA resolution/threshold parameters:
  // Hardware threshold: should be between [-1,1]
  double PE_per_keV = 6.5 + err_PEkeV * 1.;
  double E_thresh_HW = 1. / PE_per_keV / E_to_keV;
  // DAMA parameters for Gaussian resolution (smearing)
  double alpha = 0.45 + dres * 0.04;
  double beta = 0.009 + dres * 0.005;
  std::cout << "DAMA detector resolution:";
  if (dres != 0)
    std::cout << " with error term: " << dres << ".";
  std::cout << "\n ==> alpha = " << alpha << ", beta=" << beta << " (keV).\n";
  std::cout << "Hardware threshold (1 PE, from PE->keV conversion):";
  if (err_PEkeV != 0)
    std::cout << " with error term: " << err_PEkeV << ".";
  std::cout << "\n ==> PE per keV = " << PE_per_keV << " (PE).";
  std::cout << "\n ==> E_HW       = " << E_thresh_HW * E_to_keV << " (keV).\n";

  // Create the Gaussian-smearing array (includes HW threshold)
  std::vector<std::vector<double>> gausVec(desteps);
  for (std::size_t i = 0; i < desteps; i++) {
    double Eobs = Egrid.r(i);
    double sigma =
        (alpha * std::sqrt(Eobs * E_to_keV) + beta * (Eobs * E_to_keV)) /
        E_to_keV;
    for (auto Eer : Egrid.r()) {
      double g = gaussian(sigma, Eobs - Eer);
      if (Eer < E_thresh_HW || Eobs < E_thresh_HW)
        g = 0;
      gausVec[i].push_back(g);
    }
  }

  // Array to store observable Rate, S
  DoubleVec3D dSdE_mv_mx_E;
  dSdE_mv_mx_E.resize(n_mv, std::vector<std::vector<double>>(
                                n_mx, std::vector<double>(desteps)));

  // Total atomic/mol. mass (in kg) [for units]
  double MN = Atot * (PhysConst::u_NMU * PhysConst::m_e_kg);

  // Calculate _observable_ rate, S, for DAMA
  // inlcuding Gaussian resolution, + hard-ware threshold
  // Converts units to counts/day/kg/keV
  for (std::size_t imv = 0; imv < n_mv; imv++) {
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.r()[imx];
      double rho_on_mxc2 = rhoDM_GeVcm3 / (mx * M_to_GeV);
      double rate_units = dsvdE_to_cm3keVday * rho_on_mxc2 / MN;
      dSdE_mv_mx_E[imv][imx] =
          convolvedRate(dsv_mv_mx_E[imv][imx], Egrid, gausVec, rate_units);
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
  DoubleVec3D S_mv_mx_E;
  S_mv_mx_E.resize(n_mv, std::vector<std::vector<double>>(
                             n_mx, std::vector<double>(num_bins)));

  // Creates S (or Sm) [actually S/E_w] - rate averaged over each energy bin
  // Also prints energy bins to screen (only for first mx, mv)
  std::cout << "\nCalculating ";
  if (do_anMod)
    std::cout << "Sm/Ew, annual modulation amplitude (/day/kg/keV)\n";
  else
    std::cout << "S/Ew, average event rate (/day/kg/keV)\n";
  std::cout << "Integrated (averaged) each energy bin.\n";
  std::cout << "(Just outputting for first mx/mv: ";
  printf("M_chi=%5.2f GeV", mxgrid.r0() * M_to_GeV);
  if (mvgrid.r0() >= 0)
    printf(" ; M_v=%6.3f MeV", mvgrid.r0() * M_to_MeV);
  std::cout << ")\n";
  for (std::size_t imv = 0; imv < n_mv; imv++) {
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      for (int i = 0; i < num_bins; i++) {
        auto ieA = Egrid.getIndex(iEbin + i * wEbin);
        auto ieB = Egrid.getIndex(iEbin + (i + 1) * wEbin);
        // if (ieA < 0)
        //   ieA = 0;

        double Rate = 0;
        for (auto ie = ieA; ie < ieB; ie++) {
          // if (ieB >= desteps)
          //   break;
          double dEdu = Egrid.drdu()[ie]; // r[ie];
          Rate += double(dSdE_mv_mx_E[imv][imx][ie]) * dEdu;
          // nb: E is from Jacobian; * dE/E below
        }
        Rate *= Egrid.du() / wEbin;
        S_mv_mx_E[imv][imx][i] = Rate;
        double EaKev = Egrid.r()[ieA] * E_to_keV;
        double EbKev = Egrid.r()[ieB] * E_to_keV;
        if (imv == 0 && imx == 0) {
          // only print first one to screen
          printf("%3.1f-%3.1f: %6.3f   %.2e     %i\n", EaKev, EbKev,
                 0.5 * (EaKev + EbKev), Rate, (int)(ieB - ieA));
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
    of << "#Rate for DAMA: integrated (averaged) into bins. "
       << "Units: counts/day/kg/keV\n";
    of << "# sigma^bar_e = " << sbe_1e37_cm2 << " cm^2\n";
    for (std::size_t imv = 0; imv < n_mv; imv++) {
      double mv = mvgrid.r()[imv];
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
      for (std::size_t imx = 0; imx < n_mx; imx++) {
        double mx = mxgrid.r()[imx] * M_to_GeV;
        of << mx << " ";
        for (std::size_t ie = 0; ie < (std::size_t)num_bins; ie++) {
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
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.r()[imx];
      of << "\"" << std::fixed << std::setprecision(2) << mx * M_to_GeV
         << " GeV\"   ";
      for (int i = 0; i < num_bins; i++) {
        double EaKev = (iEbin + i * wEbin) * E_to_keV;
        double EbKev = EaKev + wEbin * E_to_keV;
        of << "\"" << std::fixed << std::setprecision(1) << EaKev << "-"
           << EbKev << " keV\"   ";
      }
      of << "\n" << std::scientific << std::setprecision(6);
      for (std::size_t imv = 0; imv < n_mv; imv++) {
        double mv = mvgrid.r()[imv] * M_to_MeV;
        of << mv << " ";
        for (std::size_t ie = 0; ie < (std::size_t)num_bins; ie++) {
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
void doXe100(const DoubleVec3D &dsv_mv_mx_E, const Grid &mvgrid,
             const Grid &mxgrid, const Grid &Egrid, double Atot, double N_err,
             double sPMT_err, bool write_dSds1, bool do_anMod,
             bool write_integrated, double s1_a, double s1_b,
             const std::string &label)
/*
Mostly, coming from:
 - The XENON100 Collaboration, Phys. Rev. D 90, 062009 (2014).
 - The XENON100 Collaboration, Astropart. Phys. 54, 11 (2014).
*/
{

  std::cout << "\n*** Doing S1 for Xe100 ***\n\n";
  std::size_t n_mv = mvgrid.num_points();
  std::size_t n_mx = mxgrid.num_points();

  // NofE: Conversion between energy recoil and PhotoElectrons
  // From: Fig. 2 of Phys. Rev. D 90, 062009 (2014).
  auto NofE = [](double E_au, double err) {
    double max = 1.23 * std::pow(E_au * E_to_keV, 1.43482);
    double min = 0.75 * std::pow(E_au * E_to_keV, 1.63265);
    double a = 0.5 * (err + 1.);
    double b = 1. - a;
    return a * max + b * min;
  };

  // Generate P_n(N(E)) on n, E grids
  std::size_t n_max = 20; // input?
  std::vector<std::vector<double>> P(n_max);
  for (std::size_t n = 0; n < n_max; n++) {
    // I don't use n=0; but P must be rectangular..
    std::vector<double> Pn;
    Pn.reserve(Egrid.num_points());
    for (auto E : Egrid.r()) {
      auto Ne = NofE(E, N_err);
      Pn.push_back(Pois((int)n, Ne));
    }
    P[n] = Pn;
  }
  std::cout << "Xe100: eER -> PE conversion N(e): ";
  if (N_err != 0)
    std::cout << "with error term: " << N_err << ".";
  std::cout << "\n ==> N(1keV) = " << NofE(1. / E_to_keV, N_err)
            << ", N(2keV) = " << NofE(2. / E_to_keV, N_err) << " PE\n";

  // Total atomic/mol. mass (in kg)
  double MN = Atot * (PhysConst::u_NMU * PhysConst::m_e_kg);

  // Array to store Poisson-smear S1 rate, F
  DoubleVec3D F_mv_mx_n;
  F_mv_mx_n.resize(
      n_mv, std::vector<std::vector<double>>(n_mx, std::vector<double>(n_max)));

  // Calculate Poiss-smeared rate, F [mv, mv, n]
  // F has units: counts/kg/day
  for (std::size_t imv = 0; imv < n_mv; imv++) {
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      double mx = mxgrid.r()[imx];
      double rho_on_mxc2 = rhoDM_GeVcm3 / (mx * M_to_GeV);
      double rate_units = dsvdE_to_cm3_per_auday * rho_on_mxc2 / MN;
      F_mv_mx_n[imv][imx] =
          convolvedRate(dsv_mv_mx_E[imv][imx], Egrid, P, rate_units);
    }
  }

  // Convert Poisson-smeared rate F, to Gaussian-smeared observable, S
  // S has units counts/day/s1 [s1=PE]

  // create s1 grid
  std::size_t num_s1 = Egrid.num_points(); // use same number of points!
  double s1min = 1.;                       // Threshold...  always 1 or above
  double s1max = 15.;                      // 6 kev gives less than 16 PE
  // Grid s1grid(s1min, s1max, num_s1, GridType::logarithmic);
  std::cout << "\nWARNING: Change Grid construction, not tested! XXX\n";
  // Grid s1grid(s1min, s1max, num_s1, GridType::logarithmic);
  Grid s1grid({num_s1, s1min, s1max, 0, GridType::logarithmic});

  // Acceptance
  // Fig. 1 of Phys. Rev. D 90, 062009 (2014).
  auto Aacc = [](double s1) { return 0.88 * (1. - std::exp(-0.33 * s1)); };

  // Gaussian s1 (PE) resolution
  double sigma_pmt = 0.5 * (1 + sPMT_err * 0.05);
  // PMT resolution: 0.5 PE, Astropart. Phys. 54, 11 (2014).
  // Note: No uncertainty in this is given. I take 5%. Only for OoM estimates!

  std::cout << "PMT resolution: ";
  if (sPMT_err != 0)
    std::cout << "with error term: " << sPMT_err << " (" << 5. * sPMT_err
              << "%)";
  std::cout << "\n ==> sigma_PMT = " << sigma_pmt << " (PE)\n";

  std::vector<std::vector<double>> GA_s1_n(num_s1);
  for (std::size_t is1 = 0; is1 < num_s1; is1++) {
    double s1 = s1grid.r()[is1];
    std::vector<double> G_n(n_max);
    double A_acc = Aacc(s1);
    G_n[0] = 0;
    for (std::size_t n = 1; n < n_max; n++) {
      double sigma = std::sqrt((double)n) * sigma_pmt;
      // GA_s1_n[is1]
      G_n[n] = gaussian(sigma, s1 - (double)n) * A_acc;
    }
    GA_s1_n[is1] = G_n;
  }

  DoubleVec3D dS_mv_mx_s1; // dS/ds1
  dS_mv_mx_s1.resize(n_mv, std::vector<std::vector<double>>(
                               n_mx, std::vector<double>(num_s1)));

  // Sum over n [F -> dS/ds1]
  // dS/ds1 units: counts/day/kg/PE
  for (std::size_t imv = 0; imv < n_mv; imv++) {
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      for (std::size_t is1 = 0; is1 < num_s1; is1++) {
        double sum_over_n = 0;
        for (std::size_t n = 1; n < n_max; n++) {
          sum_over_n += F_mv_mx_n[imv][imx][n] * GA_s1_n[is1][n];
        }
        dS_mv_mx_s1[imv][imx][is1] = sum_over_n;
      }
    }
  }

  // output dS/ds1 for gnuplot:
  if (write_dSds1) {
    std::string spref = do_anMod ? "dSmds1" : "dSds1";
    std::string fn_dSds1 = spref + "_mx-" + label + ".out";
    std::cout << "Writing to file: " << fn_dSds1 << "\n";
    double ux = 1; // already converted!
    double uy = 1; // already converted!
    if (n_mv == 1 || n_mx > 1) {
      writeForGnuplot_mvBlock(dS_mv_mx_s1, mvgrid, mxgrid, s1grid, fn_dSds1, uy,
                              ux);
    }
    if (n_mv > 1) {
      fn_dSds1 = spref + "_mv-" + label + ".out";
      std::cout << "Writing to file: " << fn_dSds1 << "\n";
      writeForGnuplot_mxBlock(dS_mv_mx_s1, mvgrid, mxgrid, s1grid, fn_dSds1, uy,
                              ux);
    }
  }

  // s1 integrations. Final rates (in counts/day/kg)
  auto is1_a = s1grid.getIndex(s1_a, true);
  auto is1_b = s1grid.getIndex(s1_b, true);
  std::vector<std::vector<double>> rate(n_mv, std::vector<double>(n_mx));
  for (std::size_t imv = 0; imv < n_mv; imv++) {
    for (std::size_t imx = 0; imx < n_mx; imx++) {
      rate[imv][imx] = NumCalc::integrate(s1grid.du(), is1_a, is1_b,
                                          dS_mv_mx_s1[imv][imx], s1grid.drdu());
    }
  }

  // Output some of the integrated s1 rates to screen.
  // Only outputs the first mv, and first few mxs
  std::cout << "\nIntegrate s1 rate, for: " << s1_a << " - " << s1_b << " PE\n";
  std::cout << "Mx (Gev)    Rate (counts/day/kg)\n";
  for (std::size_t imx = 0; imx < n_mx; imx += 2) {
    printf("  %6.2f      %.2e\n", mxgrid.r()[imx] * M_to_GeV, rate[0][imx]);
    if (imx > 8) {
      std::cout << "(only printed a few Mx's)\n\n";
      break;
    }
  }

  if (!write_integrated)
    return;

  // Write out integrated S1 rates (Xe100)
  std::string range =
      "_" + std::to_string((int)s1_a) + "-" + std::to_string((int)s1_b) + "PE";
  std::string fn = do_anMod ? "S1m" : "S1";
  fn = fn + range;

  // Write as function of Mx (for each Mv)
  if (mxgrid.num_points() > 1) {
    auto fn_S1 = fn + "_mx-mv_" + label + ".out";
    std::ofstream of(fn_S1);
    of << "# Integrated rate for Xe100 (" << s1_a << " - " << s1_b << " PE)."
       << " units: counts/day/kg\n";
    of << "# sigma^bar_e = " << sbe_1e37_cm2 << " cm^2\n";
    of << "Mx(GeV)\\Mv(MeV) ";
    for (std::size_t imv = 0; imv < mvgrid.num_points(); imv++) {
      double mv = mvgrid.r()[imv];
      if (mv >= 0)
        of << "\"" << std::fixed << std::setprecision(2) << mv * M_to_MeV
           << " MeV\"   ";
      else
        of << "\"infty MeV\"   ";
    }
    of << "\n";
    for (std::size_t imx = 0; imx < mxgrid.num_points(); imx++) {
      of << std::scientific << std::setprecision(4);
      of << mxgrid.r()[imx] * M_to_GeV << " ";
      for (std::size_t imv = 0; imv < mvgrid.num_points(); imv++) {
        of << std::scientific << std::setprecision(4) << rate[imv][imx] << " ";
      }
      of << "\n";
    }
    of.close();
    std::cout << "Written Integrated S1 rate for Xe100 to file: " << fn_S1
              << "\n";
  }

  // Write as function of Mv (for each Mx)
  if (mvgrid.num_points() > 1) {
    auto fn_S1 = fn + "_mv-mx_" + label + ".out";
    std::ofstream of(fn_S1);
    of << "# Integrated rate for Xe100 (" << s1_a << " - " << s1_b << " PE)."
       << " units: counts/day/kg\n";
    of << "# sigma^bar_e = " << sbe_1e37_cm2 << " cm^2\n";
    of << "Mv(MeV)\\Mx(GeV) ";
    for (std::size_t imx = 0; imx < mxgrid.num_points(); imx++) {
      double mx = mxgrid.r()[imx];
      of << "\"" << std::fixed << std::setprecision(2) << mx * M_to_GeV
         << " GeV\"   ";
    }
    of << "\n";
    for (std::size_t imv = 0; imv < mvgrid.num_points(); imv++) {
      of << std::scientific << std::setprecision(4);
      of << mvgrid.r()[imv] * M_to_MeV << " ";
      for (std::size_t imx = 0; imx < mxgrid.num_points(); imx++) {
        of << std::scientific << std::setprecision(4) << rate[imv][imx] << " ";
      }
      of << "\n";
    }
    of.close();
    std::cout << "Written Integrated S1 rate for Xe100 to file: " << fn_S1
              << "\n";
  }

  return;
}

//******************************************************************************
//******************************************************************************
//******************************************************************************
int main(int argc, char *argv[]) {
  IO::ChronoTimer sw; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "dmeXSection.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // define input parameters
  std::string akfn;
  double mxmin, mxmax, mvmin, mvmax; // m_chi and m_v masses
  int n_mx, n_mv;                    // number of steps in m_chi, and m_v grids
  std::string label;                 // label to append to output file names
  double Atot_DAMA;                  // Total nuclear mass number (Na+I)
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

  bool do_Xe100;
  double N_err, sPMT_err;
  bool write_dSds1;
  bool write_integrated;
  double Atot_XE;
  double s1_a;
  double s1_b;

  // Open and read the input file:
  {
    int i_mv, iwr_dsvde, idodama, ianmod, iSM; // temp settings
    int iXe100, iwdS1;
    auto tp = std::forward_as_tuple(
        akfn, label, mxmin, mxmax, n_mx, i_mv, mvmin, mvmax, n_mv, dvesc, dv0,
        ianmod, iwr_dsvde, idodama, dres, err_PEkeV, Atot_DAMA, iEbin, fEbin,
        wEbin, iSM, iXe100, N_err, sPMT_err, Atot_XE, s1_a, s1_b, iwdS1);
    IO::FRW::setInputParameters(input_file, tp);
    label = label == "na" ? akfn : akfn + "-" + label;
    // what to write to file:
    write_dsvde = iwr_dsvde == 1;
    write_SofM = iSM == 0 ? false : true;
    write_dSdE = iSM == 1 ? false : true;
    // What to calclate:
    do_anMod = ianmod == 1 ? true : false;
    do_DAMA = idodama == 1 ? true : false;
    do_Xe100 = iXe100 == 1 ? true : false;
    write_dSds1 = iwdS1 == 1 ? false : true;
    write_integrated = iwdS1 == 0 ? false : true;
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

  // Number of points in the velocity integration.
  const int vsteps = 100; // no need to be input

  // DM mass: Convert from GeV to au + define grids
  mxmin /= M_to_GeV;
  mxmax /= M_to_GeV;
  if (n_mx == 1)
    mxmax = mxmin;
  // Grid mxgrid(mxmin, mxmax, n_mx, GridType::logarithmic);
  Grid mxgrid({std::size_t(n_mx), mxmin, mxmax, 0, GridType::logarithmic});

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
  Grid mvgrid({std::size_t(n_mv), mvmin, mvmax, 0, GridType::logarithmic});

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
  auto desteps = Kenq.size();
  auto num_states = Kenq[0].size();
  auto qsteps = Kenq[0][0].size();
  if (num_states != nklst.size())
    return 1; // just sanity check
  // Create the E and q grids:
  // Grid Egrid(demin, demax, desteps, GridType::logarithmic);
  // Grid qgrid(qmin, qmax, qsteps, GridType::logarithmic);
  Grid Egrid({desteps, demin, demax, 0, GridType::logarithmic});
  Grid qgrid({qsteps, qmin, qmax, 0, GridType::logarithmic});

  // Grid of SHM vel. distro f_v(v). Can use to change vel profiles
  // Note: SHM is in km/s units, both for v and f!
  // f.dv = 1 => [f] = [1/v]

  double max_v = (AstroConsts::v_max) / V_to_kms;
  double dv = max_v / vsteps;
  int num_cp = do_anMod ? 3 : 1; // just 1 v. dist? or 3 (for an mod.)?
  std::vector<std::vector<double>> arr_fv(num_cp, std::vector<double>(vsteps));
  for (int icp = 0; icp < num_cp; icp++) {
    double cosp = icp == 0 ? 0 : std::pow(-1, icp); // 0, -1, 1
    StandardHaloModel shm(cosp, dvesc, dv0);
    for (int i = 0; i < vsteps; i++) {
      double v = (i + 1) * dv; // don't include zero
      double vkms =
          v * (PhysConst::c_SI / PhysConst::c) / 1.e3; // will be in km/s
      arr_fv[icp][i] =
          shm.fv(vkms) / (1. / V_to_kms); // convert to a.u. [f]=[1/v]
    }
  }
  std::cout << "\nUsing SHM velocity distro";
  if (dvesc != 0 || dv0 != 0)
    std::cout << " with error terms (dvesc, dv0)=" << dvesc << "," << dv0
              << ".";
  std::cout << "\n ==> v0=" << AstroConsts::v0 + dv0 * AstroConsts::delta_v0
            << ", "
            << "vesc =" << AstroConsts::v_esc + dvesc * AstroConsts::delta_v_esc
            << " km/s\n";

  // Print the grid info to screen:
  if (do_anMod)
    std::cout << "Running for annual modulation amplitude\n";
  std::cout << "\nGrids:\n";
  printf("q :%6.2f -> %6.2f MeV, N=%4i\n", qmin * Q_to_MeV, qmax * Q_to_MeV,
         (int)qsteps);
  printf("E :%6.2f -> %6.2f keV, N=%4i\n", demin * E_to_keV, demax * E_to_keV,
         (int)desteps);
  printf("v :%6.2f -> %6.2fkm/s, N=%4i\n", dv * V_to_kms, max_v * V_to_kms,
         (int)vsteps);
  printf("Mx:%6.2f -> %6.2f GeV, N=%4i\n", mxmin * M_to_GeV, mxmax * M_to_GeV,
         (int)n_mx);
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

  double al_x = (std::sqrt(sbe_1e37_cm2) / PhysConst::aB_cm) *
                (PhysConst::alpha / std::sqrt(16. * M_PI));
  double al_xH = al_x / (PhysConst::m_e_MeV * PhysConst::alpha2);
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
  DoubleVec3D dsv_mv_mx_E; // Array to store cross-section
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

  // Calculate and output obervable rates for XENON100
  if (do_Xe100)
    doXe100(dsv_mv_mx_E, mvgrid, mxgrid, Egrid, Atot_XE, N_err, sPMT_err,
            write_dSds1, do_anMod, write_integrated, s1_a, s1_b, label);

  // Calculate and output obervable rate for DAMA
  if (do_DAMA)
    doDAMA(dsv_mv_mx_E, mvgrid, mxgrid, Egrid, do_anMod, write_dSdE, write_SofM,
           dres, err_PEkeV, Atot_DAMA, iEbin, fEbin, wEbin, label);

  std::cout << "\nTotal: " << sw.reading_str() << "\n";
  return 0;
}
