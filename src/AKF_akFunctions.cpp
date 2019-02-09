#include "AKF_akFunctions.h"
#include "ATI_atomInfo.h"
#include "ContinuumOrbitals.h"
#include "ElectronOrbitals.h"
#include "ExponentialGrid.h"
#include "FPC_physicalConstants.h"
#include "FileIO_fileReadWrite.h"
#include "SBF_sphericalBesselFunctions.h"
#include "Wigner_369j.h"
#include <cmath>
#include <fstream>
#include <iostream>

namespace AKF {

//******************************************************************************
double CLkk(int L, int ka, int kb)
/*
Angular coeficient (nb: is already squared)
*/
{
  int la = ATI::l_k(ka);
  int lb = ATI::l_k(kb);
  int two_ja = ATI::twoj_k(ka);
  int two_jb = ATI::twoj_k(kb);
  double ja = 0.5 * two_ja;
  double jb = 0.5 * two_jb;

  if ((la + lb + L) % 2 != 0)
    return 0; // Parity rule
  if ((la + lb < L) || (abs(la - lb) > L))
    return 0; // triangle rule (l)

  double tjB = Wigner::threej(jb, ja, L, -0.5, 0.5, 0);
  return (2 * ja + 1) * (2 * jb + 1) * (2 * L + 1) * pow(tjB, 2);
}

//******************************************************************************
void writeToTextFile(std::string fname,
                     const std::vector<std::vector<std::vector<float>>> &AK,
                     const std::vector<std::string> &nklst, double qmin,
                     double qmax, double demin, double demax)
/*
Writes the K factor to a text-file, in GNU-plot readable format
*/
{
  int desteps = (int)AK.size();       // dE
  int num_states = (int)AK[0].size(); // nk
  int qsteps = (int)AK[0][0].size();  // q

  double qMeV = (1.e6 / (FPC::Hartree_eV * FPC::c));
  double keV = (1.e3 / FPC::Hartree_eV);

  std::ofstream ofile;
  ofile.open(fname + ".txt");
  ofile << "dE(keV) q(MeV) ";
  for (size_t i = 0; i < nklst.size(); i++)
    ofile << nklst[i] << " "; // xxx can range!
  ofile << "Sum\n\n";
  for (int i = 0; i < desteps; i++) {
    for (int k = 0; k < qsteps; k++) {
      double x = double(k) / (qsteps - 1);
      if (qsteps == 1)
        x = 0;
      double q = qmin * pow(qmax / qmin, x);
      double y = double(i) / (desteps - 1);
      if (desteps == 1)
        y = 0;
      double dE = demin * pow(demax / demin, y);
      ofile << dE / keV << " " << q / qMeV << " ";
      float sum = 0.f;
      for (int j = 0; j < num_states; j++) {
        sum += AK[i][j][k];
        ofile << AK[i][j][k] << " ";
      }
      ofile << sum << "\n";
    }
    if (qsteps > 1)
      ofile << "\n";
  }
  ofile.close();
}

//******************************************************************************
int akReadWrite(std::string fname, bool write,
                std::vector<std::vector<std::vector<float>>> &AK,
                std::vector<std::string> &nklst, double &qmin, double &qmax,
                double &dEmin, double &dEmax)
/*
Writes K function (+ all required size etc.) values to a binary file.
The binary file is read by other programs (e.g., dmeXSection)
Uses FileIO_fileReadWrite
*/
{
  FileIO::RoW row = write ? FileIO::write : FileIO::read;

  std::fstream iof;
  fname = fname + ".bin";
  FileIO::open_binary(iof, fname, row);

  if (iof.fail()) {
    std::cout << "Can't open " << fname << "\n";
    return 1;
  }

  if (write) {
    int nde = (int)AK.size();      // dE
    int ns = (int)AK[0].size();    // nk
    int nq = (int)AK[0][0].size(); // q
    FileIO::binary_rw(iof, nde, row);
    FileIO::binary_rw(iof, ns, row);
    FileIO::binary_rw(iof, nq, row);
  } else {
    int nq, ns, nde;
    FileIO::binary_rw(iof, nde, row);
    FileIO::binary_rw(iof, ns, row);
    FileIO::binary_rw(iof, nq, row);
    AK.resize(nde, std::vector<std::vector<float>>(ns, std::vector<float>(nq)));
    nklst.resize(ns);
  }
  FileIO::binary_rw(iof, qmin, row);
  FileIO::binary_rw(iof, qmax, row);
  FileIO::binary_rw(iof, dEmin, row);
  FileIO::binary_rw(iof, dEmax, row);
  for (size_t ie = 0; ie < AK.size(); ie++) {
    for (size_t in = 0; in < AK[0].size(); in++) {
      if (ie == 0)
        FileIO::binary_str_rw(iof, nklst[in], row);
      for (size_t iq = 0; iq < AK[0][0].size(); iq++) {
        FileIO::binary_rw(iof, AK[ie][in][iq], row);
      }
    }
  }

  return 0;
}

//******************************************************************************
int calculateK_nk(const ElectronOrbitals &wf, int is, int max_L, double dE,
                  std::vector<std::vector<std::vector<float>>> &jLqr_f,
                  std::vector<float> &AK_nk_q, double Zeff)
/*
Calculates the atomic factor for a given core state (is) and energy.
Note: dE = I + ec is depositied energy, not cntm energy
Zeff is '-1' by default. If Zeff > 0, will solve w/ Zeff model
Zeff no longer works at main() level.
*/
{
  ContinuumOrbitals cntm(wf); // create cntm object [survives locally only]

  int k = wf.ka(is);
  int l = wf.lorb(is);

  int qsteps = (int)jLqr_f[0].size();

  // Calculate continuum wavefunctions
  double ec = dE + wf.en[is];
  cntm.clear();
  int lc_max = l + max_L;
  int lc_min = l - max_L;
  if (lc_min < 0)
    lc_min = 0;
  if (ec > 0) {
    if (Zeff > 0)
      cntm.solveZeffContinuum(ec, Zeff, lc_min, lc_max); // Zeff version
    else
      cntm.solveLocalContinuum(ec, lc_min, lc_max);
  }

  double x_ocf = wf.occ_frac[is]; // occupancy fraction. Usually 1

  // Generate AK for each L, lc, and q
  // L and lc are summed, not stored indevidually
  for (int L = 0; L <= max_L; L++) {
    for (size_t ic = 0; ic < cntm.kappa.size(); ic++) {
      int kc = cntm.kappa[ic];
      double dC_Lkk = CLkk(L, k, kc);
      if (dC_Lkk == 0)
        continue;
      //#pragma omp parallel for
      for (int iq = 0; iq < qsteps; iq++) {
        double a = 0.;
        double jLqr = 0.;
        if (cntm.f.size() > 0) {
          int maxj = wf.pinflist[is]; // don't bother going further
          // Do the radial integral:
          a = 0;
          for (int j = 0; j < maxj; j++) {
            jLqr = (double)jLqr_f[L][iq][j];
            a += (wf.f[is][j] * cntm.f[ic][j] + wf.g[is][j] * cntm.g[ic][j]) *
                 jLqr * wf.drdt[j]; // *h below!
          }
        }
        AK_nk_q[iq] += (float)(dC_Lkk * pow(a * wf.h, 2) * x_ocf);
      } // q
    }   // END loop over cntm states (ic)
  }     // end L loop
  // cntm.clear(); //deletes cntm wfs for this energy
  return 0;
}

//******************************************************************************
int calculateKpw_nk(const ElectronOrbitals &wf, int nk, double dE,
                    std::vector<std::vector<float>> &jl_qr,
                    std::vector<float> &tmpK_q)
/*
For plane-wave final state.
Only has f-part....Can restore g-part, but need to be sure of plane-wave!
Chi(q) - Int[ f_nk*j_l(qr)*r , {r,0,inf}]
Should be called once per initial state

XXX Note sure if correct! esp, (q) angular part!? XXX

*/
{
  if (nk >= (int)wf.f.size())
    return 1; // should never occur

  int twoj = wf.twoj(nk);

  int qsteps = (int)jl_qr.size();

  double eps = dE - wf.en[nk];
  int maxir = wf.pinflist[nk]; // don't bother going further

  for (int iq = 0; iq < qsteps; iq++) {
    if (eps <= 0)
      break;
    double chi_q = 0.;
    for (int ir = 0; ir < maxir; ir++)
      chi_q += wf.f[nk][ir] * double(jl_qr[iq][ir]) * wf.r[ir] * wf.drdt[ir];
    chi_q *= wf.h;
    tmpK_q[iq] =
        (float)((2. / M_PI) * (twoj + 1) * pow(chi_q, 2) * sqrt(2. * eps));
    // tmpK_q[iq] = pow(4*3.14159,2)*pow(chi_q,2); //XXX XXX just cf KOPP
  }

  return 0;
}

//******************************************************************************
void sphericalBesselTable(std::vector<std::vector<std::vector<float>>> &jLqr_f,
                          int max_L, const ExpGrid &qgrid,
                          const std::vector<double> &r)
/*
Creates a look-up table w/ spherical Bessel functions. For speed.
Uses SBF_sphericalBesselFunctions
*/
{
  std::cout << std::endl;
  int ngp = (int)r.size();
  int qsteps = qgrid.N();

  jLqr_f.resize(max_L + 1, std::vector<std::vector<float>>(
                               qsteps, std::vector<float>(ngp)));
  for (int L = 0; L <= max_L; L++) {
    std::cout << "\rCalculating spherical Bessel look-up table for L=" << L
              << "/" << max_L << " .. " << std::flush;
#pragma omp parallel for
    for (int iq = 0; iq < qsteps; iq++) {
      double q = qgrid.x(iq);
      for (int ir = 0; ir < ngp; ir++) {
        double tmp = SBF::JL(L, q * r[ir]);
        // If q(dr) is too large, "missing" j_L oscillations
        //(overstepping them). This helps to fix that.
        // By averaging the J_L function. Note: only works if wf is smooth
        int num_extra = 0;
        if (ir < ngp - 1) {
          double qdrop = q * (r[ir + 1] - r[ir]) / M_PI;
          double min_qdrop = 0.01; // require 100 pts per half wavelength!
          if (qdrop > min_qdrop)
            num_extra = int(qdrop / min_qdrop) + 3;
        }
        { // Include 'extra' points into j_L (avg):
          for (int i = 0; i < num_extra; i++) {
            double b = (i + 1.) / (num_extra + 1.);
            double a = 1. - b;
            double qrtmp = q * (a * r[ir] + b * r[ir + 1]);
            tmp += SBF::JL(L, qrtmp);
          }
          tmp /= (num_extra + 1);
        }
        jLqr_f[L][iq][ir] = (float)tmp;
      }
    }
  }
  std::cout << "done\n";
}

} // namespace AKF
