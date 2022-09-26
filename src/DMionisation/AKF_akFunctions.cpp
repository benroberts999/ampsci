#include "DMionisation/AKF_akFunctions.hpp"
#include "Angular/Wigner369j.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Maths/SphericalBessel.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>

namespace AKF {

//******************************************************************************
double CLkk_DLkk(int L, int ka, int kb, std::string dmec)
// /*
// Selects angular coefficient based on DM-electron coupling
// Vector & scalar cases:
// C_{k}^{k',L} = [j][j'][L] * (j,j',L, -1/,1/2,0)^2 * pi(l+l'+L)
// Pseudovector & pseudoscalar cases:
// D_{k}^{k',L} related to C_{k}^{k',L} via transformation
// k -> k_tilde = -k, l -> l_tilde = |k_tilde + 1/2| - 1/2
// */
{
  int la;
  int lb = Angular::l_k(kb);
  int two_ja = Angular::twoj_k(ka);
  int two_jb = Angular::twoj_k(kb);

  if (dmec == "Pseudovector" || dmec == "Pseudoscalar") {
    int la_tilde = Angular::l_k(-ka);
    la = la_tilde;
  } else {
    la = Angular::l_k(ka);
  }

  if ((la + lb + L) % 2 != 0)
    return 0; // Parity rule

  const double tjs = Angular::threej_2(two_jb, two_ja, 2 * L, -1, 1, 0);
  return (two_ja + 1) * (two_jb + 1) * (2 * L + 1) * tjs * tjs;
}

//******************************************************************************
double Heaviside(double x)
// Heaviside step function: Theta(x) = {(1 for x>0) and (0 for x<=0)}
{
  return (x > 0.0) ? 1.0 : 0.0;
}

//==============================================================================
std::vector<float>
calculateK_nk(const Wavefunction &wf, const DiracSpinor &psi, int max_L,
              double dE,
              const std::vector<std::vector<std::vector<double>>> &jLqr_f,
              bool alt_akf, bool force_rescale, bool subtract_self,
              bool force_orthog, std::string dmec, bool zeff_cont)
// Calculates the atomic factor for a given core state (is) and energy.
// Note: dE = I + ec is depositied energy, not cntm energy
// Zeff is '-1' by default. If Zeff > 0, will solve w/ Zeff model
// Zeff no longer works at main() level.
{
  ContinuumOrbitals cntm(wf); // create cntm object [survives locally only]

  const auto qsteps = jLqr_f[0].size();

  std::vector<float> AK_nk_q;
  AK_nk_q.resize(qsteps);

  // Calculate continuum wavefunctions
  double ec = dE + psi.en();
  // cntm.clear();
  const int l = psi.l();
  int lc_max = l + max_L;
  int lc_min = std::max(l - max_L, 0);
  if (ec > 0) {
    if (zeff_cont) {
      cntm.solveContinuumZeff(ec, lc_min, lc_max, psi.en(), (double)psi.n(),
                              force_orthog, &psi);
    } else {
      cntm.solveContinuumHF(ec, lc_min, lc_max, force_rescale, subtract_self,
                            force_orthog, &psi);
    }
  }

  const double x_ocf = psi.occ_frac(); // occupancy fraction. Usually 1

  // Generate AK for each L, lc, and q
  // L and lc are summed, not stored individually
  for (std::size_t L = 0; L <= std::size_t(max_L); L++) {
    for (const auto &phic : cntm.orbitals) {
      // double dC_Lkk_sqrt = CLkk_sqrt((int)L, psi.k, phic.k);
      double dC_Lkk = CLkk_DLkk((int)L, psi.kappa(), phic.kappa(), dmec);
      if (dC_Lkk == 0)
        continue;
      //#pragma omp parallel for
      for (std::size_t iq = 0; iq < qsteps; iq++) {
        double a = 0.;
        // double a_noj = 0;
        auto maxj = psi.max_pt(); // don't bother going further
        // Need to change this so that aff, agg, afg, & agf are all
        // calculated everytime -> separate function or if checks?
        double aff = NumCalc::integrate(1.0, 0, maxj, psi.f(), phic.f(),
                                        jLqr_f[L][iq], wf.grid().drdu());
        double agg = NumCalc::integrate(1.0, 0, maxj, psi.g(), phic.g(),
                                        jLqr_f[L][iq], wf.grid().drdu());
        double afg = NumCalc::integrate(1.0, 0, maxj, psi.f(), phic.g(),
                                        jLqr_f[L][iq], wf.grid().drdu());
        double agf = NumCalc::integrate(1.0, 0, maxj, psi.g(), phic.f(),
                                        jLqr_f[L][iq], wf.grid().drdu());

        // Not working?
        // double aff_noj = NumCalc::integrate(1.0, 0, maxj, psi.f(), phic.f(),
        //                                     wf.grid().drdu());
        // double agg_noj = NumCalc::integrate(1.0, 0, maxj, psi.g(), phic.g(),
        //                                     wf.grid().drdu());
        // a = (aff + agg) * wf.rgrid->du();
        // a_noj = (aff_noj + agg_noj) * wf.rgrid->du();
        // if ((alt_akf) && (psi.k == phic.k)) {
        //   AK_nk_q[iq] += (float)((dC_Lkk * a * a + a_noj * a_noj -
        //                           2.0 * dC_Lkk_sqrt * a * a_noj) *
        //                          x_ocf);
        // } else {
        //   AK_nk_q[iq] += (float)(dC_Lkk * a * a * x_ocf);
        // }

        // if ((alt_akf) && (psi.k == phic.k)) {
        //   aff -= NumCalc::integrate(1.0, 0, maxj, psi.f(), phic.f(),
        //                             wf.rgrid->drdu());
        //   agg -= NumCalc::integrate(1.0, 0, maxj, psi.g(), phic.g(),
        //                             wf.rgrid->drdu());
        //   afg -= NumCalc::integrate(1.0, 0, maxj, psi.f(), phic.g(),
        //                             wf.rgrid->drdu());
        //   agf -= NumCalc::integrate(1.0, 0, maxj, psi.g(), phic.f(),
        //                             wf.rgrid->drdu());
        // }
        if (alt_akf) {
          // getting rid of compiler warning for now
        }

        if (dmec == "Vector") {
          a = (aff + agg);
        } else if (dmec == "Scalar") {
          a = aff - agg;
        } else if (dmec == "Pseudovector") {
          a = afg - agf;
        } else if (dmec == "Pseudoscalar") {
          a = afg + agf;
        }
        AK_nk_q[iq] +=
            (float)(dC_Lkk * std::pow(a * wf.grid().du(), 2) * x_ocf);
      } // q
    }   // END loop over cntm states (ic)
  }     // end L loop
  return AK_nk_q;
}

//******************************************************************************
std::vector<float> stepK_nk(const DiracSpinor &psi, double dE,
                            const std::vector<float> &AFBE_table)
/*
Approximates atomic factor using table generated by AFBindingEnergy
*/
{
  // q range needs to be the same for the table and AK
  const auto qsteps = AFBE_table.size();

  std::vector<float> AK_nk_step;
  AK_nk_step.resize(qsteps);

  for (std::size_t iq = 0; iq < qsteps; iq++) {
    AK_nk_step[iq] += (float)Heaviside(dE + psi.en()) * AFBE_table[iq];
  }
  // }
  return AK_nk_step;
}

//******************************************************************************
std::vector<float>
calculateKpw_nk(const Wavefunction &wf, const DiracSpinor &psi, double dE,
                const std::vector<std::vector<double>> &jl_qr)
// /*
// For plane-wave final state.
// Only has f-part....Can restore g-part, but need to be sure of plane-wave!
// Chi(q) - Int[ f_nk*j_l(qr)*r , {r,0,inf}]
// Should be called once per initial state
{

  const int twoj = psi.twoj(); // wf.twoj(nk);

  const auto qsteps = jl_qr.size();
  std::vector<float> tmpK_q(qsteps);

  const double eps = dE - psi.en();
  const auto maxir = psi.max_pt(); // don't bother going further

  if (eps <= 0)
    return tmpK_q;

  for (auto iq = 0ul; iq < qsteps; iq++) {
    const double chi_q =
        NumCalc::integrate(wf.grid().du(), 0, maxir, psi.f(), jl_qr[iq],
                           wf.grid().r(), wf.grid().drdu());
    tmpK_q[iq] = (float)((2. / M_PI) * (twoj + 1) * std::pow(chi_q, 2) *
                         std::sqrt(2. * eps));
    // tmpK_q[iq] = std::pow(4*3.14159,2)*std::pow(chi_q,2); // just cf KOPP
  }

  return tmpK_q;
}

//==============================================================================
void write_Knk_plaintext(const std::string &fname,
                         const std::vector<std::vector<std::vector<float>>> &AK,
                         const std::vector<std::string> &nklst,
                         const Grid &qgrid, const Grid &Egrid)
// /*
// Writes the K factor to a text-file, in GNU-plot readable format
// */
{
  const auto desteps = AK.size();       // dE
  const auto num_states = AK[0].size(); // nk
  const auto qsteps = AK[0][0].size();  // q

  const double qMeV = (1.e6 / (PhysConst::Hartree_eV * PhysConst::c));
  const double keV = (1.e3 / PhysConst::Hartree_eV);

  std::ofstream ofile;
  ofile.open(fname + ".txt");
  ofile << "dE(keV) q(MeV) ";
  for (const auto &nk : nklst) {
    ofile << nk << " ";
  }
  ofile << "Sum\n\n";
  for (auto idE = 0ul; idE < desteps; idE++) {
    for (auto iq = 0ul; iq < qsteps; iq++) {
      const auto q = qgrid.r(iq);
      const auto dE = Egrid.r(idE);
      ofile << dE / keV << " " << q / qMeV << " ";
      float sum = 0.0f;
      for (auto j = 0ul; j < num_states; j++) {
        sum += AK[idE][j][iq];
        ofile << AK[idE][j][iq] << " ";
      }
      ofile << sum << "\n";
    }
    if (qsteps > 1)
      ofile << "\n";
  }
  ofile.close();
}

//==============================================================================
void write_Ktot_plaintext(
    const std::string &fname,
    const std::vector<std::vector<std::vector<float>>> &AK, const Grid &qgrid,
    const Grid &Egrid)
// /*
// Writes the K factor to a text-file, in GNU-plot readable format
// XXX NOTE: Re-creates grids! Could use Grid class!
// XXX This mean we MUST use exponential Grid! Fix this! XXX
// */
{
  const auto desteps = AK.size();       // dE
  const auto num_states = AK[0].size(); // nk
  const auto qsteps = AK[0][0].size();  // q

  const double qMeV = (1.e6 / (PhysConst::Hartree_eV * PhysConst::c));
  const double keV = (1.e3 / PhysConst::Hartree_eV);

  std::ofstream ofile;
  ofile.open(fname + "_tot.txt");
  ofile
      << "# K_tot(dE,q): dE is rows, q is cols; first row/col is dE/q values\n";
  ofile << "q(MeV)/dE(keV)";
  for (auto idE = 0ul; idE < desteps; idE++) {
    const auto dE = Egrid.r(idE);
    ofile << " " << dE / keV;
  }
  ofile << "\n";

  for (auto iq = 0ul; iq < qsteps; iq++) {
    const auto q = qgrid.r(iq);
    ofile << q / qMeV;
    for (auto idE = 0ul; idE < desteps; idE++) {

      float AK_dEq = 0.0f;
      for (auto j = 0ul; j < num_states; j++) {
        AK_dEq += AK[idE][j][iq];
      }
      ofile << " " << AK_dEq;
    }
    ofile << "\n";
  }
  ofile.close();
}

//******************************************************************************
void writeToTextFile_AFBE(
    const std::string &fname,
    const std::vector<std::vector<std::vector<float>>> &AK,
    const std::vector<std::string> &nklst, const Grid &qgrid,
    const std::vector<double> deion)
// /*
// Writes the K factor to a text-file, in GNU-plot readable format
// XXX NOTE: Re-creates grids! Could use Grid class!
// XXX This mean we MUST use exponential Grid! Fix this! XXX
// */
{
  // int desteps = (int)AK.size();       // dE
  const auto num_states = AK[0].size(); // nk
  const auto qsteps = AK[0][0].size();  // q

  const double qMeV = (1.e6 / (PhysConst::Hartree_eV * PhysConst::c));
  const double keV = (1.e3 / PhysConst::Hartree_eV);

  std::ofstream ofile;
  ofile.open(fname + ".txt");
  ofile << "q(MeV) ";
  for (const auto &nk : nklst) {
    ofile << "dE(keV) " << nk << " ";
  }
  ofile << "Sum\n\n";
  for (auto iq = 0ul; iq < qsteps; iq++) {
    const auto q = qgrid.r(iq);
    ofile << q / qMeV << " ";
    float sum = 0.0f;
    for (auto j = 0ul; j < num_states; j++) {
      sum += AK[0][j][iq];
      ofile << deion[j] / keV << " " << AK[0][j][iq] << " ";
    }
    ofile << sum << "\n";
  }
  if (qsteps > 1)
    ofile << "\n";
  ofile.close();
}

//******************************************************************************
int akReadWrite(const std::string &fname, bool write,
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
  IO::FRW::RoW row = write ? IO::FRW::write : IO::FRW::read;

  std::fstream iof;
  IO::FRW::open_binary(iof, fname + ".bin", row);

  if (iof.fail()) {
    std::cout << "Can't open " << fname << ".bin\n";
    return 1;
  }

  if (write) {
    auto nde = AK.size();      // dE
    auto ns = AK[0].size();    // nk
    auto nq = AK[0][0].size(); // q
    IO::FRW::binary_rw(iof, nde, row);
    IO::FRW::binary_rw(iof, ns, row);
    IO::FRW::binary_rw(iof, nq, row);
  } else {
    std::size_t nq, ns, nde;
    IO::FRW::binary_rw(iof, nde, row);
    IO::FRW::binary_rw(iof, ns, row);
    IO::FRW::binary_rw(iof, nq, row);
    std::cout << "nde = " << nde << ", ns = " << ns << ", nq = " << nq
              << std::endl;
    AK.resize(nde, std::vector<std::vector<float>>(ns, std::vector<float>(nq)));
    nklst.resize(ns);
  }
  IO::FRW::binary_rw(iof, qmin, row);
  IO::FRW::binary_rw(iof, qmax, row);
  IO::FRW::binary_rw(iof, dEmin, row);
  IO::FRW::binary_rw(iof, dEmax, row);
  for (std::size_t ie = 0; ie < AK.size(); ie++) {
    for (std::size_t in = 0; in < AK[0].size(); in++) {
      if (ie == 0)
        IO::FRW::binary_str_rw(iof, nklst[in], row);
      for (std::size_t iq = 0; iq < AK[0][0].size(); iq++) {
        IO::FRW::binary_rw(iof, AK[ie][in][iq], row);
      }
    }
  }

  return 0;
}

//******************************************************************************
int akReadWrite_AFBE(const std::string &fname, bool write,
                     std::vector<std::vector<std::vector<float>>> &AK,
                     std::vector<std::string> &nklst, double &qmin,
                     double &qmax, std::vector<double> &deion)
// /*
// Writes K function (+ all required size etc.) values to a binary file.
// The binary file is read by other programs (e.g., dmeXSection)
// Uses FileIO_fileReadWrite
// XXX NOTE: Re-creates grids! Could use Grid class!
// XXX This mean we MUST use exponential Grid! Fix this! XXX
// */
{
  IO::FRW::RoW row = write ? IO::FRW::write : IO::FRW::read;

  std::fstream iof;
  IO::FRW::open_binary(iof, fname + ".bin", row);

  if (iof.fail()) {
    std::cout << "Can't open " << fname << ".bin\n";
    return 1;
  }

  if (write) {
    auto ns = AK[0].size();    // nk
    auto nq = AK[0][0].size(); // q
    IO::FRW::binary_rw(iof, ns, row);
    IO::FRW::binary_rw(iof, nq, row);
  } else {
    std::size_t nq, ns;
    IO::FRW::binary_rw(iof, ns, row);
    IO::FRW::binary_rw(iof, nq, row);
    AK.resize(1, std::vector<std::vector<float>>(ns, std::vector<float>(nq)));
    nklst.resize(ns);
  }
  IO::FRW::binary_rw(iof, qmin, row);
  IO::FRW::binary_rw(iof, qmax, row);
  for (std::size_t in = 0; in < AK[0].size(); in++) {
    IO::FRW::binary_str_rw(iof, nklst[in], row);
    IO::FRW::binary_rw(iof, deion[in], row);
    for (std::size_t iq = 0; iq < AK[0][0].size(); iq++) {
      IO::FRW::binary_rw(iof, AK[0][in][iq], row);
    }
  }

  return 0;
} // namespace AKF

/*
int calculateK_nk()
*/

//******************************************************************************
std::vector<std::vector<std::vector<double>>>
sphericalBesselTable(int max_L, const std::vector<double> &q_array,
                     const std::vector<double> &r)
// /*
// Creates a look-up table w/ spherical Bessel functions. For speed.
// Uses SphericalBessel
// */
{
  std::cout << "\n";
  const auto num_points = r.size();
  const auto qsteps = q_array.size();

  std::vector<std::vector<std::vector<double>>> jLqr_f;
  jLqr_f.resize(std::size_t(max_L) + 1,
                std::vector<std::vector<double>>(
                    qsteps, std::vector<double>(num_points)));
  for (auto L = 0u; L <= unsigned(max_L); L++) {
    std::cout << "\rCalculating spherical Bessel look-up table for L=" << L
              << "/" << max_L << " .. " << std::flush;
#pragma omp parallel for collapse(2)
    for (auto iq = 0ul; iq < qsteps; iq++) {
      for (auto ir = 0ul; ir < num_points; ir++) {
        const double q = q_array[iq];
        double tmp = SphericalBessel::JL(int(L), q *r[ir]);
        // If q(dr) is too large, "missing" j_L oscillations
        //(overstepping them). This helps to fix that.
        // By averaging the J_L function. Note: only works if wf is smooth
        int num_extra = 0;
        if (ir < num_points - 1) {
          const double qdrop = q * (r[ir + 1] - r[ir]) / M_PI;
          const double min_qdrop = 0.01; // require 100 pts per half wavelength!
          if (qdrop > min_qdrop)
            num_extra = int(qdrop / min_qdrop) + 3;
        }
        { // Include 'extra' points into j_L (avg):
          for (auto i = 0; i < num_extra; i++) {
            const double b = (i + 1.0) / (num_extra + 1.0);
            const double a = 1.0 - b;
            const double qrtmp = q * (a * r[ir] + b * r[ir + 1]);
            tmp += SphericalBessel::JL(int(L), qrtmp);
          }
          tmp /= (num_extra + 1);
        }
        jLqr_f[L][iq][ir] = tmp;
      }
    }
  }
  std::cout << "done\n" << std::flush;
  return jLqr_f;
}

} // namespace AKF
