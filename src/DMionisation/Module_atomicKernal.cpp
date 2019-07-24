#include "Module_atomicKernal.hpp"
#include "../ChronoTimer.hpp"
#include "../ContinuumOrbitals.hpp"
#include "../Grid.hpp"
#include "../HartreeFockClass.hpp"
#include "../Physics/PhysConst_constants.hpp"
#include "../UserInput.hpp"
#include "../Wavefunction.hpp"
#include "AKF_akFunctions.hpp"
#include <iostream>

namespace Module {

//******************************************************************************
void atomicKernal(const UserInputBlock &input, const Wavefunction &wf) {
  ChronoTimer timer; // start the overall timer

  auto demin = input.get<double>("Emin");
  auto demax = input.get<double>("Emax");
  auto desteps = input.get<int>("Esteps");

  auto qmin = input.get<double>("qmin");
  auto qmax = input.get<double>("qmax");
  auto qsteps = input.get<int>("qsteps");

  // allow for single-step in dE or q grid
  if (desteps == 1)
    demax = demin;
  if (qsteps == 1)
    qmax = qmin;

  // Convert units for input q and dE range into atomic units
  double keV = (1.0e3 / PhysConst::Hartree_eV);
  demin *= keV;
  demax *= keV;
  double qMeV = (1.0e6 / (PhysConst::Hartree_eV * PhysConst::c));
  qmin *= qMeV;
  qmax *= qMeV;

  // Set up the E and q grids
  Grid Egrid(demin, demax, desteps, GridType::logarithmic);
  Grid qgrid(qmin, qmax, qsteps, GridType::logarithmic);

  auto max_l_core = wf.maxCore_l();
  auto max_l = input.get<int>("max_l_bound", max_l_core);
  if (max_l < 0 || max_l > max_l_core)
    max_l = max_l_core;
  auto max_L = input.get<int>("max_L");

  bool plane_wave = input.get<bool>("use_plane_waves", false);
  if (plane_wave)
    max_L = max_l; // for spherical bessel.

  auto label = input.get<std::string>("label", "");

  // output format
  auto text_out = input.get<bool>("output_text", false);
  auto bin_out = input.get<bool>("output_binary", false);
  if (!text_out && !bin_out)
    bin_out = true; // print message?

  // Make sure h (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  double du_target = (M_PI / 20.) / sqrt(2. * demax);
  auto du = wf.rgrid.du;
  if (du > du_target) {
    auto new_ngp = Grid::calc_ngp_from_du(wf.rgrid.r0, wf.rgrid.rmax, du_target,
                                          GridType::loglinear, 3.5);
    auto old_ngp = wf.rgrid.ngp;
    // ngp = (int)new_ngp;
    std::cerr
        << "\nWARNING 118: Grid not dense enough for contimuum state with "
        << "ec=" << demax << "au\n";
    std::cerr << "You should update ngp from " << old_ngp << " --> " << new_ngp
              << "\n";
    std::cerr << "Program will continue, but may fail\n";
  }

  // outut file name (excluding extension):
  std::string fname = "ak-" + wf.atomicSymbol(); // + "_" + label;
  if (label != "")
    fname += "_" + label;

  // Print some info to screen:
  std::cout << "\nRunning Atomic Kernal for " << wf.atom() << "\n";
  std::cout << "*************************************************\n";
  // std::cout << "Radial " << wf.rgrid.gridParameters() << "\n\n";

  // Output HF results:
  std::cout << "  state   k     En (au)    En (eV)   Oc.Frac.\n";
  for (const auto &phi : wf.core_orbitals) {
    // double rinf = wf.rinf(phi);
    // printf("%2i)%7s %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f   (%.2f)\n",
    //        i++, phi.symbol().c_str(), phi.k, rinf, phi.its, phi.eps, phi.en,
    //        phi.en * PhysConst::Hartree_invcm, phi.en * PhysConst::Hartree_eV,
    //        phi.occ_frac);
    printf(" %7s %2i %11.5f %10.2f   [%3.2f]", phi.symbol().c_str(), phi.k,
           phi.en, phi.en * PhysConst::Hartree_eV, phi.occ_frac);
    if (phi.l() > max_l)
      std::cout << " (excluded from K)";
    std::cout << "\n";
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
  return;
}

} // namespace Module
