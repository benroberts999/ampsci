#include "DMionisation/AKF_akFunctions.hpp"
#include "Angular/Wigner369j.hpp"
#include "DiracOperator/DiracOperator.hpp"
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

std::vector<double> K_q(double dE, const HF::HartreeFock *hf,
                        const DiracOperator::jL &jl) {
  const auto q_points = jl.q_grid().num_points();
  std::vector<double> kq_out(q_points);

  const int lc_max = DiracSpinor::max_l(hf->core()) + jl.max_L();

  ContinuumOrbitals cntm(hf); // create cntm object [survives locally only]
  for (auto &Fa : hf->core()) {
    const auto ec = dE + Fa.en();
    if (ec < 0.0)
      continue;
    cntm.clear();
    cntm.solveContinuumHF(ec, lc_max, &Fa);
    const auto x = Fa.occ_frac(); // normally 1
#pragma omp parallel for
    for (std::size_t iq = 0; iq < q_points; ++iq) {
      for (const auto &Fe : cntm.orbitals) {
        for (std::size_t L = 0; L <= jl.max_L(); ++L) {
          if (jl.is_zero(Fa, Fe, L))
            continue;
          const auto q = jl.q_grid().r(iq);
          const auto me = jl.rme(Fa, Fe, L, q);
          kq_out.at(iq) += double(2 * L + 1) * me * me * x;
        }
      }
    }
  }
  return kq_out;
}

//==============================================================================
double CLkk(int L, int ka, int kb)
// /*
// Angular coeficient (nb: is already squared)
// C_{k}^{k',L} = [j][j'][L] * (j,j',L, -1/,1/2,0)^2 * pi(l+l'+L)
// */
{
  const int la = Angular::l_k(ka);
  const int lb = Angular::l_k(kb);
  const int two_ja = Angular::twoj_k(ka);
  const int two_jb = Angular::twoj_k(kb);

  if ((la + lb + L) % 2 != 0)
    return 0; // Parity rule
  // Note: parity rule not included in 3j, so must be checked!

  const double tjs = Angular::threej_2(two_jb, two_ja, 2 * L, -1, 1, 0);
  return (two_ja + 1) * (two_jb + 1) * (2 * L + 1) * tjs * tjs;
}

//==============================================================================
std::vector<float>
calculateK_nk(const Wavefunction &wf, const DiracSpinor &psi, int max_L,
              double dE,
              const std::vector<std::vector<std::vector<double>>> &jLqr_f)
// Calculates the atomic factor for a given core state (is) and energy.
// Note: dE = I + ec is depositied energy, not cntm energy
{
  ContinuumOrbitals cntm(wf); // create cntm object [survives locally only]

  const auto qsteps = jLqr_f[0].size();

  std::vector<float> AK_nk_q;
  AK_nk_q.resize(qsteps); // use +=, so can't use .push_back()

  // Calculate continuum wavefunctions
  const double ec = dE + psi.en();
  // XXX Correct? Seems fine. Triangle rule is actually for j (not l)
  const int l = psi.l();
  const int lc_max = l + max_L; // +1?
  const int lc_min = std::max(l - max_L, 0);
  if (ec > 0) {
    cntm.solveContinuumHF(ec, lc_min, lc_max, &psi);
  }

  const double x_ocf = psi.occ_frac(); // occupancy fraction. Usually 1

  // Generate AK for each L, lc, and q
  // L and lc are summed, not stored indevidually
  for (std::size_t L = 0; L <= std::size_t(max_L); L++) {
    for (const auto &phic : cntm.orbitals) {
      const double dC_Lkk = CLkk(int(L), psi.kappa(), phic.kappa());
      if (dC_Lkk == 0)
        continue;
      for (std::size_t iq = 0; iq < qsteps; iq++) {
        const auto rad_int = psi * (jLqr_f[L][iq] * phic);
        AK_nk_q[iq] += (float)(dC_Lkk * rad_int * rad_int * x_ocf);
      } // q
    }   // END loop over cntm states (ic)
  }     // end L loop

  return AK_nk_q;
}

//==============================================================================
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
  ofile << "# K_tot(dE,q): dE is rows, q is cols; first row/col is dE/q "
           "values\n";
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

//==============================================================================
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

//==============================================================================
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
        double tmp = SphericalBessel::JL(int(L), q * r[ir]);
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
