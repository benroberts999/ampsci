#include "AKF_akFunctions.h"
#include "ATI_atomInfo.h"
#include "ChronoTimer.h"
#include "ContinuumOrbitals.h"
#include "ElectronOrbitals.h"
#include "ExponentialGrid.h"
#include "FPC_physicalConstants.h"
#include "FileIO_fileReadWrite.h"
#include "HartreeFockClass.h"
#include "PRM_parametricPotentials.h"
#include <cmath>
#include <iostream>
#include <tuple>

//******************************************************************************
int main(int argc, char *argv[]) {
  ChronoTimer sw; // start the overall timer

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
  if (Gf < 0)
    Gf = 0;

  // If L<0, will use plane-waves (instead of cntm fns)
  bool plane_wave = (max_L < 0) ? true : false;

  // default Hartree convergance goal:
  if (hart_del == 0)
    hart_del = 1.e-6;

  // allow for single-step in dE or q grid
  if (desteps == 1)
    demax = demin;
  if (qsteps == 1)
    qmax = qmin;
  // Set up the E and q grids
  ExpGrid Egrid(desteps, demin, demax);
  ExpGrid qgrid(qsteps, qmin, qmax);

  // Fix maximum angular momentum values:
  if (max_l < 0 || max_l > 3)
    max_l = 3; // default: all core states (no >f)
  if (plane_wave)
    max_L = max_l; // for spherical bessel.

  // alpha can't be zero, just make v. small
  if (varalpha == 0)
    varalpha = 1.e-25;

  // Convert units for input q and dE range into atomic units
  double keV = (1.e3 / FPC::Hartree_eV);
  demin *= keV;
  demax *= keV;
  double qMeV = (1.e6 / (FPC::Hartree_eV * FPC::c));
  qmin *= qMeV;
  qmax *= qMeV;

  // Look-up atomic number, Z, and also A
  int Z = ATI::get_z(Z_str);
  if (Z == 0)
    return 2;

  // Generate the orbitals object:
  ElectronOrbitals wf(Z, A, ngp, r0, rmax, varalpha);

  // outut file name (excluding extension):
  std::string fname = "ak-" + Z_str + "_" + label;

  // Write out as text and/or binary file
  bool text_out = (iout_format == 1) ? false : true;
  bool bin_out = (iout_format > 0) ? true : false;

  // Print some info to screen:
  printf("\nRunning Atomic Kernal for %s, Z=%i A=%i\n", Z_str.c_str(), Z,
         wf.Anuc());
  printf("*************************************************\n");
  if (Gf != 0)
    printf("Using Green potential: H=%.4f  d=%.4f\n", Gh, Gd);
  else
    printf("Using Hartree Fock (converge to %.0e)\n", hart_del);

  // Make sure h (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  // Also need to re-form nuclear potential.
  // Updated class method will avoid this!
  double h_target = (M_PI / 20.) / sqrt(2. * demax);
  if (wf.h > h_target) {
    int old_ngp = ngp;
    wf.logLinearRadialGrid(h_target, r0, rmax);
    if (wf.Anuc() > 0)
      wf.formNuclearPotential(NucleusType::Fermi);
    else
      wf.formNuclearPotential(NucleusType::zero);
    ngp = wf.ngp;
    std::cout
        << "\nWARNING 101: Grid not dense enough for contimuum state with "
        << "ec=" << demax << "au\n";
    std::cout << "Updating ngp: " << old_ngp << " --> " << ngp << "\n";
  }
  printf("Grid: pts=%i h=%6.4f r0=%.0e Rmax=%5.1f\n", wf.ngp, wf.h,
         wf.r.front(), wf.r.back());

  // Do Hartree-fock (or parametric potential) for Core
  if (Gf == 0) {
    HartreeFock hf(wf, str_core, hart_del);
  } else {
    // Use Green (local parametric) potential
    // Fill the electron part of the (local/direct) potential
    wf.vdir.reserve(wf.ngp);
    for (auto r : wf.r)
      wf.vdir.push_back(PRM::green(Z, r, Gh, Gd));
    std::cerr << "143\n";
    std::cin.get();
    wf.solveInitialCore(str_core); // solves w/ Green
  }

  // Output results:
  std::cout << "\n     state  k Rinf its    eps      En (au)     En (/cm)    "
            << "En (eV)   Oc.Frac.\n";
  for (int i : wf.stateIndexList) {
    auto nlj = wf.seTermSymbol(i);
    int k = wf.ka(i);
    double rinf = wf.rinf(i);
    double eni = wf.en[i];
    double x = wf.occ_frac[i];
    printf("%2i)%7s %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f   (%.2f)\n", i,
           nlj.c_str(), k, rinf, wf.itslist[i], wf.epslist[i], eni,
           eni * FPC::Hartree_invcm, eni * FPC::Hartree_eV, x);
  }

  //////////////////////////////////////////////////

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<float>>> AK; // float ok?
  int num_states = (int)wf.nlist.size();
  AK.resize(desteps, std::vector<std::vector<float>>(
                         num_states, std::vector<float>(qsteps)));

  // Store state info (each orbital) [just useful for plotting!]
  std::vector<std::string> nklst; // human-readiable state labels (easy
                                  // plotting)
  nklst.reserve(wf.nlist.size());
  for (auto i : wf.stateIndexList)
    nklst.emplace_back(wf.seTermSymbol(i, true));

  // pre-calculate the spherical Bessel function look-up table for efficiency
  std::vector<std::vector<std::vector<float>>> jLqr_f;
  AKF::sphericalBesselTable(jLqr_f, max_L, qgrid, wf.r);

  // Calculate the AK (print to screen)
  std::cout << "\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)\n", demin / keV,
         demax / keV, demin, demax);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)\n", qmin / qMeV,
         qmax / qMeV, qmin, qmax);

  // Calculate K(q,E)
  std::cout << "Running dE loops (" << desteps << ").." << std::flush;
#pragma omp parallel for
  for (int ide = 0; ide < desteps; ide++) {
    double dE = Egrid.x(ide);
    // Loop over core (bound) states:
    for (auto is : wf.stateIndexList) {
      int l = wf.lorb(is);
      if (l > max_l)
        continue;
      if (plane_wave)
        AKF::calculateKpw_nk(wf, is, dE, jLqr_f[l], AK[ide][is]);
      else
        AKF::calculateK_nk(wf, is, max_L, dE, jLqr_f, AK[ide][is]);
    } // END loop over bound states
  }
  std::cout << "..done :)\n";

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

  std::cout << "\n " << sw.reading_str() << "\n";
  return 0;
}
