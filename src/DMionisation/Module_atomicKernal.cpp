#include "DMionisation/Module_atomicKernal.hpp"
#include "DMionisation/AKF_akFunctions.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>

namespace Module {

//******************************************************************************
void atomicKernal(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer;

  input.check(
      {{"Emin", "[keV] minimum energy transfer (dE) ~0.1"},
       {"Emax", "[keV] maximum dE"},
       {"Esteps", "numer of steps along dE grid (logarithmic grid)"},
       {"qmin", "[MeV] minimum momentum transfer (q) ~0.01"},
       {"qmax", "[MeV] maximum q"},
       {"qsteps", "numer of steps along q grid (logarithmic grid)"},
       {"max_l_bound",
        "Max. orbital ang. mom. for bound states to include (only for tests)"},
       {"max_L", "Maximum multipolarity used in exp(iqr) expansion"},
       {"use_plane_waves", "true/false"},
       {"label", "label for output files"},
       {"output_text", "true/false (output K(dE,q) to .txt)"},
       {"output_binary", "true/false (output K(dE,q) to .bin)"}});

  // Read input: q and dE grids:
  const auto demin_kev = input.get<double>("Emin", 0.1);
  const auto demax_kev = input.get<double>("Emax", demin_kev);
  auto desteps = input.get<std::size_t>("Esteps", 1);
  const auto qmin_Mev = input.get<double>("qmin", 0.01);
  const auto qmax_Mev = input.get<double>("qmax", qmin_Mev);
  auto qsteps = input.get<std::size_t>("qsteps", 1);

  // Convert units for input q and dE range into atomic units
  const double keV = (1.0e3 / PhysConst::Hartree_eV);
  auto demin = demin_kev * keV;
  auto demax = (desteps == 1) ? demin : demax_kev * keV;
  const double qMeV = (1.0e6 / (PhysConst::Hartree_eV * PhysConst::c));
  auto qmin = qmin_Mev * qMeV;
  auto qmax = (qsteps == 1) ? qmin : qmax_Mev * qMeV;

  // Set up the E and q grids
  const Grid Egrid({desteps, demin, demax, 0, GridType::logarithmic});
  const Grid qgrid({qsteps, qmin, qmax, 0, GridType::logarithmic});

  // read in max l and L
  const auto max_l_core = wf.maxCore_l();
  auto max_l = input.get<int>("max_l_bound", max_l_core);
  if (max_l < 0 || max_l > max_l_core)
    max_l = max_l_core;
  auto max_L = input.get<int>("max_L", 2 * max_l); // random default..

  const bool plane_wave = input.get<bool>("use_plane_waves", false);
  if (plane_wave)
    max_L = max_l; // for spherical bessel.

  const auto label = input.get<std::string>("label", "");

  // output format
  auto text_out = input.get<bool>("output_text", false);
  auto bin_out = input.get<bool>("output_binary", false);
  if (!text_out && !bin_out)
    bin_out = true; // print message?

  // Make sure h (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  const double du_target = (M_PI / 20.) / std::sqrt(2. * demax);
  const auto du = wf.rgrid->du();
  if (du > du_target) {
    const auto new_num_points = Grid::calc_num_points_from_du(
        wf.rgrid->r0(), wf.rgrid->rmax(), du_target, GridType::loglinear, 4.0);
    const auto old_num_points = wf.rgrid->num_points();
    // num_points = (int)new_num_points;
    std::cerr
        << "\nWARNING 118: Grid not dense enough for contimuum state with "
        << "ec=" << demax << "au\n";
    std::cerr << "You should update num_points from " << old_num_points
              << " --> " << new_num_points << "\n";
    std::cerr << "Program will continue, but may fail\n";
  }

  // outut file name (excluding extension):
  std::string fname = "ak-" + wf.atomicSymbol();
  if (label != "")
    fname += "_" + label;

  // Print some info to screen:
  std::cout << "\nRunning Atomic Kernal for " << wf.atom() << "\n";
  std::cout << "*************************************************\n";

  // Output HF results:
  std::cout << "  state   k     En (au)    En (eV)   Oc.Frac.\n";
  for (const auto &phi : wf.core) {
    printf(" %7s %2i %11.5f %10.2f   [%3.2f]", phi.symbol().c_str(), phi.k,
           phi.en(), phi.en() * PhysConst::Hartree_eV, phi.occ_frac());
    if (phi.l() > max_l)
      std::cout << " (excluded from K)";
    std::cout << "\n";
  }
  //////////////////////////////////////////////////

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<float>>> AK; // float ok?
  // AK[i_dE][n_core]
  const auto num_states = wf.core.size();
  AK.resize(desteps, std::vector<std::vector<float>>(num_states));

  // Store state info (each orbital) [just useful for plotting!]
  std::vector<std::string> nklst; // human-readiable state labels (easy
                                  // plotting)
  nklst.reserve(wf.core.size());
  for (auto &phi : wf.core)
    nklst.emplace_back(phi.symbol(true));

  // pre-calculate the spherical Bessel function look-up table for efficiency
  timer.start();
  const auto jLqr_f =
      AKF::sphericalBesselTable(max_L, qgrid.r(), wf.rgrid->r());
  std::cout << "Time for SB table: " << timer.lap_reading_str() << "\n";

  // Calculate the AK (print to screen)
  std::cout << "\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)  [N=%i]\n", demin / keV,
         demax / keV, demin, demax, (int)desteps);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)  [N=%i]\n", qmin / qMeV,
         qmax / qMeV, qmin, qmax, (int)qsteps);

  // Calculate K(q,E)
  timer.start();
  std::cout << "Running dE loops (" << desteps << ").." << std::flush;
#pragma omp parallel for
  for (std::size_t ide = 0; ide < desteps; ide++) {
    const double dE = Egrid.r(ide);
    // Loop over core (bound) states:
    for (std::size_t is = 0; is < wf.core.size(); is++) {
      const auto &phi_c = wf.core[is];
      const auto l = std::size_t(phi_c.l());
      if ((int)l > max_l)
        continue;
      if (plane_wave)
        AK[ide][is] = AKF::calculateKpw_nk(wf, phi_c, dE, jLqr_f[l]);
      else
        AK[ide][is] = AKF::calculateK_nk(wf, phi_c, max_L, dE, jLqr_f);
    } // END loop over bound states
  }
  std::cout << "..done :)\n";
  std::cout << "Time for AK: " << timer.lap_reading_str() << "\n";

  // Write out to text file (in gnuplot friendly form)
  if (text_out) {
    AKF::write_Knk_plaintext(fname, AK, nklst, qgrid, Egrid);
    AKF::write_Ktot_plaintext(fname, AK, qgrid, Egrid);
  }
  // //Write out AK as binary file
  if (bin_out) {
    AKF::akReadWrite(fname, true, AK, nklst, qmin, qmax, demin, demax);
  }
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
