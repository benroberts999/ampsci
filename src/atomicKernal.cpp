#include "AKF_akFunctions.hpp"
#include "AtomInfo.hpp"
#include "ChronoTimer.hpp"
#include "ContinuumOrbitals.hpp"
#include "FileIO_fileReadWrite.hpp"
#include "Grid.hpp"
#include "HartreeFockClass.hpp"
#include "Parametric_potentials.hpp"
#include "PhysConst_constants.hpp"
#include "Wavefunction.hpp"
#include <iostream>
#include <tuple>

//******************************************************************************
int main(int argc, char *argv[]) {
  ChronoTimer timer; // start the overall timer

  std::string input_file = (argc > 1) ? argv[1] : "atomicKernal.in";
  std::cout << "Reading input from: " << input_file << "\n";

  // Input options
  std::string Z_str;
  int A;
  double r0, rmax;
  int ngp;
  double varalpha;      // for non-relativistic approx
  double hart_del;      // HART convergance
  std::string str_core; // States for the core
  // Green potential parameters
  int Gf;
  double Gh, Gd;
  // q and dE grids:
  double qmin, qmax, demin, demax;
  int qsteps, desteps;
  // Max anglular momentums
  int max_l, max_L;
  int iout_format;
  std::string label; // label for output file

  {
    auto tp =
        std::forward_as_tuple(Z_str, A, str_core, r0, rmax, ngp, Gf, Gh, Gd,
                              hart_del, varalpha, demin, demax, desteps, qmin,
                              qmax, qsteps, max_l, max_L, iout_format, label);
    FileIO::setInputParameters(input_file, tp);
  }
  if (Gf < 0 || Gf > 2)
    Gf = 0;

  // If L<0, will use plane-waves (instead of cntm fns)
  bool plane_wave = (max_L < 0) ? true : false;

  // default Hartree convergance goal:
  if (hart_del == 0)
    hart_del = 1.e-6;

  // Fix maximum angular momentum values:
  if (max_l < 0 || max_l > 3)
    max_l = 3; // default: all core states (no >f)
  if (plane_wave)
    max_L = max_l; // for spherical bessel.

  // alpha can't be zero, just make v. small
  if (varalpha == 0)
    varalpha = 1.e-25;

  // allow for single-step in dE or q grid
  if (desteps == 1)
    demax = demin;
  if (qsteps == 1)
    qmax = qmin;

  // Convert units for input q and dE range into atomic units
  double keV = (1.e3 / PhysConst::Hartree_eV);
  demin *= keV;
  demax *= keV;
  double qMeV = (1.e6 / (PhysConst::Hartree_eV * PhysConst::c));
  qmin *= qMeV;
  qmax *= qMeV;

  // Set up the E and q grids
  Grid Egrid(demin, demax, desteps, GridType::logarithmic);
  Grid qgrid(qmin, qmax, qsteps, GridType::logarithmic);

  // Look-up atomic number, Z
  int Z = AtomInfo::get_z(Z_str);

  // Make sure h (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  double du_target = (M_PI / 20.) / sqrt(2. * demax);
  double du = Grid::calc_du_from_ngp(r0, rmax, ngp, GridType::loglinear, 3.5);
  if (du > du_target) {
    auto new_ngp =
        Grid::calc_ngp_from_du(r0, rmax, du_target, GridType::loglinear, 3.5);
    int old_ngp = ngp;
    ngp = (int)new_ngp;
    std::cout
        << "\nWARNING 101: Grid not dense enough for contimuum state with "
        << "ec=" << demax << "au\n";
    std::cout << "Updating ngp: " << old_ngp << " --> " << ngp << "\n";
  }

  // Generate the orbitals object:
  Wavefunction wf(Z, A, ngp, r0, rmax, varalpha);

  // outut file name (excluding extension):
  std::string fname = "ak-" + Z_str + "_" + label;

  // Write out as text and/or binary file
  bool text_out = (iout_format == 1) ? false : true;
  bool bin_out = (iout_format > 0) ? true : false;

  // Print some info to screen:
  std::cout << "\nRunning Atomic Kernal for " << Z_str << ", Z=" << Z
            << " A=" << wf.Anuc() << "\n";
  std::cout << "*************************************************\n";
  if (Gf == 1)
    printf("Using Green parametric potential: H=%.4f  d=%.4f\n", Gh, Gd);
  else if (Gf == 2)
    printf("Using Hartree [no exchange] (converge to %.0e)\n", hart_del);
  else
    printf("Using Hartree Fock (converge to %.0e)\n", hart_del);

  std::cout << "Radial " << wf.rgrid.gridParameters() << "\n\n";

  // Do Hartree-fock (or parametric potential) for Core
  timer.start();
  if (Gf != 1) {
    bool excludeExchange = (Gf == 2) ? true : false;
    HartreeFock hf(wf, str_core, hart_del, excludeExchange);
  } else {
    // Use Green (local parametric) potential
    // Fill the electron part of the (local/direct) potential
    wf.vdir.reserve(wf.rgrid.ngp);
    for (auto r : wf.rgrid.r)
      wf.vdir.push_back(Parametric::green(Z, r, Gh, Gd));
    wf.solveInitialCore(str_core); // solves w/ Green
  }
  std::cout << "Time for HF: " << timer.lap_reading_str() << "\n";

  // Output HF results:
  std::cout << "\n     state  k Rinf its    eps      En (au)     En (/cm)    "
            << "En (eV)   Oc.Frac.\n";
  int i = 0;
  for (auto &phi : wf.core_orbitals) {
    double rinf = wf.rinf(phi);
    printf("%2i)%7s %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f   (%.2f)\n",
           i++, phi.symbol().c_str(), phi.k, rinf, phi.its, phi.eps, phi.en,
           phi.en * PhysConst::Hartree_invcm, phi.en * PhysConst::Hartree_eV,
           phi.occ_frac);
  }
  //////////////////////////////////////////////////

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<float>>> AK; // float ok?
  int num_states = (int)wf.core_orbitals.size();
  AK.resize(desteps, std::vector<std::vector<float>>(
                         num_states, std::vector<float>(qsteps)));

  // Store state info (each orbital) [just useful for plotting!]
  std::vector<std::string> nklst; // human-readiable state labels (easy
                                  // plotting)
  nklst.reserve(wf.core_orbitals.size());
  // for (auto i : wf.stateIndexList)
  //   nklst.emplace_back(wf.seTermSymbol(i, true));
  for (auto &phi : wf.core_orbitals)
    nklst.emplace_back(phi.symbol(true));

  // pre-calculate the spherical Bessel function look-up table for efficiency
  timer.start();
  std::vector<std::vector<std::vector<double>>> jLqr_f;
  AKF::sphericalBesselTable(jLqr_f, max_L, qgrid.r, wf.rgrid.r);
  std::cout << "Time for SB table: " << timer.lap_reading_str() << "\n";

  // Calculate the AK (print to screen)
  std::cout << "\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)  [N=%i]\n", demin / keV,
         demax / keV, demin, demax, desteps);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)  [N=%i]\n", qmin / qMeV,
         qmax / qMeV, qmin, qmax, qsteps);

  // Calculate K(q,E)
  timer.start();
  std::cout << "Running dE loops (" << desteps << ").." << std::flush;
#pragma omp parallel for
  for (int ide = 0; ide < desteps; ide++) {
    double dE = Egrid.r[ide];
    // Loop over core (bound) states:
    // for (auto is : wf.stateIndexList) {
    for (std::size_t is = 0; is < wf.core_orbitals.size(); is++) {
      int l = wf.core_orbitals[is].l(); // lorb(is);
      if (l > max_l)
        continue;
      if (plane_wave)
        AKF::calculateKpw_nk(wf, is, dE, jLqr_f[l], AK[ide][is]);
      else
        AKF::calculateK_nk(wf, is, max_L, dE, jLqr_f, AK[ide][is]);
    } // END loop over bound states
  }
  std::cout << "..done :)\n";
  std::cout << "Time for AK: " << timer.lap_reading_str() << "\n";

  // Write out to text file (in gnuplot friendly form)
  if (text_out)
    AKF::writeToTextFile(fname, AK, nklst, qmin, qmax, demin, demax);
  // //Write out AK as binary file
  if (bin_out)
    AKF::akReadWrite(fname, true, AK, nklst, qmin, qmax, demin, demax);
  std::cout << "Written to: " << fname;
  if (text_out)
    std::cout << ".txt";
  if (text_out && bin_out)
    std::cout << ", and ";
  if (bin_out)
    std::cout << ".bin";
  std::cout << "\n";

  std::cout << "\n " << timer.reading_str() << "\n";
  return 0;
}
