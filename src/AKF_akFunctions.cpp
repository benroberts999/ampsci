#include "AKF_akFunctions.hpp"
#include "AtomInfo.hpp"
#include "ContinuumOrbitals.hpp"
#include "FileIO_fileReadWrite.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "PhysConst_constants.hpp"
#include "SBF_sphericalBessel.hpp"
#include "Wavefunction.hpp"
#include "Wigner_369j.hpp"
#include <fstream>
#include <iostream>

namespace AKF {

//******************************************************************************
double CLkk(int L, int ka, int kb)
// /*
// Angular coeficient (nb: is already squared)
// C_{k}^{k',L} = [j][j'][L] * (j,j',L, -1/,1/2,0)^2 * pi(l+l'+L)
// */
{
  int la = AtomInfo::l_k(ka);
  int lb = AtomInfo::l_k(kb);
  int two_ja = AtomInfo::twoj_k(ka);
  int two_jb = AtomInfo::twoj_k(kb);
  double ja = 0.5 * two_ja;
  double jb = 0.5 * two_jb;

  if ((la + lb + L) % 2 != 0)
    return 0; // Parity rule
  if ((la + lb < L) || (abs(la - lb) > L))
    return 0; // triangle rule (l)
  // Note: triangle rule included in 3j, so this is not needed (but faster)
  // But, parity rule not included in 3j, so must be checked!

  double tjB = Wigner::threej(jb, ja, L, -0.5, 0.5, 0);
  return (2 * ja + 1) * (2 * jb + 1) * (2 * L + 1) * pow(tjB, 2);
}

//******************************************************************************
void writeToTextFile(std::string fname,
                     const std::vector<std::vector<std::vector<float>>> &AK,
                     const std::vector<std::string> &nklst, double qmin,
                     double qmax, double demin, double demax)
// /*
// Writes the K factor to a text-file, in GNU-plot readable format
// XXX NOTE: Re-creates grids! Could use Grid class!
// XXX This mean we MUST use exponential Grid! Fix this! XXX
// */
{
  int desteps = (int)AK.size();       // dE
  int num_states = (int)AK[0].size(); // nk
  int qsteps = (int)AK[0][0].size();  // q

  double qMeV = (1.e6 / (PhysConst::Hartree_eV * PhysConst::c));
  double keV = (1.e3 / PhysConst::Hartree_eV);

  std::ofstream ofile;
  ofile.open(fname + ".txt");
  ofile << "dE(keV) q(MeV) ";
  for (auto nk : nklst) {
    ofile << nk << " ";
  }
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
// /*
// Writes K function (+ all required size etc.) values to a binary file.
// The binary file is read by other programs (e.g., dmeXSection)
// Uses FileIO_fileReadWrite
// XXX NOTE: Re-creates grids! Could use Grid class!
// XXX This mean we MUST use exponential Grid! Fix this! XXX
// */
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
  for (std::size_t ie = 0; ie < AK.size(); ie++) {
    for (std::size_t in = 0; in < AK[0].size(); in++) {
      if (ie == 0)
        FileIO::binary_str_rw(iof, nklst[in], row);
      for (std::size_t iq = 0; iq < AK[0][0].size(); iq++) {
        FileIO::binary_rw(iof, AK[ie][in][iq], row);
      }
    }
  }

  return 0;
}

//******************************************************************************
int calculateK_nk(const Wavefunction &wf, std::size_t is, int max_L, double dE,
                  std::vector<std::vector<std::vector<double>>> &jLqr_f,
                  std::vector<float> &AK_nk_q, double)
// Calculates the atomic factor for a given core state (is) and energy.
// Note: dE = I + ec is depositied energy, not cntm energy
// Zeff is '-1' by default. If Zeff > 0, will solve w/ Zeff model
// Zeff no longer works at main() level.
{
  ContinuumOrbitals cntm(wf); // create cntm object [survives locally only]
  auto &psi = wf.core_orbitals[is];

  int k = psi.k;   // wf.ka(is);
  int l = psi.l(); // wf.lorb(is);

  int qsteps = (int)jLqr_f[0].size();

  // Calculate continuum wavefunctions
  double ec = dE + wf.core_orbitals[is].en;
  cntm.clear();
  int lc_max = l + max_L;
  int lc_min = l - max_L;
  if (lc_min < 0)
    lc_min = 0;
  if (ec > 0) {
    cntm.solveLocalContinuum(ec, lc_min, lc_max);
  }

  double x_ocf = psi.occ_frac; // occupancy fraction. Usually 1

  // Generate AK for each L, lc, and q
  // L and lc are summed, not stored indevidually
  for (int L = 0; L <= max_L; L++) {
    for (const auto &phic : cntm.orbitals) {
      int kc = phic.k;
      double dC_Lkk = CLkk(L, k, kc);
      if (dC_Lkk == 0)
        continue;
      //#pragma omp parallel for
      for (int iq = 0; iq < qsteps; iq++) {
        double a = 0.;
        auto maxj = psi.pinf; // don't bother going further
        double af = NumCalc::integrate(psi.f, phic.f, jLqr_f[L][iq],
                                       wf.rgrid.drdu, 1., 0, maxj);
        double ag = NumCalc::integrate(psi.g, phic.g, jLqr_f[L][iq],
                                       wf.rgrid.drdu, 1., 0, maxj);
        a = af + ag;
        AK_nk_q[iq] += (float)(dC_Lkk * pow(a * wf.rgrid.du, 2) * x_ocf);
      } // q
    }   // END loop over cntm states (ic)
  }     // end L loop
  return 0;
}

//******************************************************************************
int calculateKpw_nk(const Wavefunction &wf, std::size_t nk, double dE,
                    std::vector<std::vector<double>> &jl_qr,
                    std::vector<float> &tmpK_q)
// /*
// For plane-wave final state.
// Only has f-part....Can restore g-part, but need to be sure of plane-wave!
// Chi(q) - Int[ f_nk*j_l(qr)*r , {r,0,inf}]
// Should be called once per initial state
{

  auto &psi = wf.core_orbitals[nk];

  int twoj = psi.twoj(); // wf.twoj(nk);

  auto qsteps = jl_qr.size();

  double eps = dE - psi.en;
  auto maxir = psi.pinf; // don't bother going further

  if (eps <= 0)
    return 0;

  for (auto iq = 0ul; iq < qsteps; iq++) {
    double chi_q = NumCalc::integrate(psi.f, jl_qr[iq], wf.rgrid.r,
                                      wf.rgrid.drdu, wf.rgrid.du, 0, maxir);
    tmpK_q[iq] =
        (float)((2. / M_PI) * (twoj + 1) * pow(chi_q, 2) * sqrt(2. * eps));
    // tmpK_q[iq] = pow(4*3.14159,2)*pow(chi_q,2); // just cf KOPP
  }

  return 0;
}

//******************************************************************************
void sphericalBesselTable(std::vector<std::vector<std::vector<double>>> &jLqr_f,
                          int max_L, const std::vector<double> &q_array,
                          const std::vector<double> &r)
// /*
// Creates a look-up table w/ spherical Bessel functions. For speed.
// Uses SBF_sphericalBessel
// */
{
  std::cout << std::endl;
  int ngp = (int)r.size();
  int qsteps = (int)q_array.size();

  jLqr_f.resize(max_L + 1, std::vector<std::vector<double>>(
                               qsteps, std::vector<double>(ngp)));
  for (int L = 0; L <= max_L; L++) {
    std::cout << "\rCalculating spherical Bessel look-up table for L=" << L
              << "/" << max_L << " .. " << std::flush;
#pragma omp parallel for
    for (int iq = 0; iq < qsteps; iq++) {
      double q = q_array[iq];
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
        jLqr_f[L][iq][ir] = tmp;
      }
    }
  }
  std::cout << "done\n";
}

} // namespace AKF
