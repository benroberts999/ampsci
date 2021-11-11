#include "DMionisation/AKF_akFunctions.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>

#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace Module {

//******************************************************************************
void AFBindingEnergy(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer; // start the overall timer

  input.checkBlock({"qmin", "qmax", "qsteps", "max_l_bound", "max_L",
                    "use_plane_waves", "output_text", "label", "use_alt_akf",
                    "force_rescale", "subtract_self", "force_orthog",
                    "dme_coupling", "Etune_mult", "Etune_add"});

  auto qmin = input.get<double>("qmin", 1.0);
  auto qmax = input.get<double>("qmax", 1.0);
  auto qsteps = input.get<int>("qsteps", 1.0);

  // allow for single-step in dE or q grid
  if (qsteps == 1)
    qmax = qmin;

  // Convert units for input q and dE range into atomic units
  // double keV = (1.0e3 / PhysConst::Hartree_eV);
  double qMeV = (1.0e6 / (PhysConst::Hartree_eV * PhysConst::c));
  qmin *= qMeV;
  qmax *= qMeV;

  // Set up the E and q grids
  // Grid Egrid(demin, demax, desteps, GridType::logarithmic);
  // Grid qgrid(qmin, qmax, qsteps, GridType::logarithmic);

  // GridParameters(std::size_t innum_points, double inr0, double inrmax,
  //                double inb = 4.0, GridType intype = GridType::loglinear,
  //                double indu = 0);

  // double demin, demax;
  int desteps = 1;

  Grid qgrid({(std::size_t)qsteps, qmin, qmax, 0, GridType::logarithmic});

  auto max_l_core = wf.maxCore_l();
  auto max_l = input.get<int>("max_l_bound", max_l_core);
  if (max_l < 0 || max_l > max_l_core)
    max_l = max_l_core;
  auto max_L = input.get<int>("max_L", 2 * max_l); // random default..

  bool plane_wave = input.get<bool>("use_plane_waves", false);
  if (plane_wave)
    max_L = max_l; // for spherical bessel.

  // output format
  // auto bin_out = input.get<bool>("output_binary", false);

  // if alt_akf then exp(iqr) -> exp(iqr) - 1 (i.e. j_L -> j_L - 1)
  auto alt_akf = input.get<bool>("use_alt_akf", false);
  // Options used by solveContinuumHF (called in AKF)
  auto force_rescale = input.get<bool>("force_rescale", false);
  auto subtract_self = input.get<bool>("subtract_self", false);
  auto force_orthog = input.get<bool>("force_orthog", false);

  auto label = input.get<std::string>("label", "");

  // DM-electron couplings
  std::vector<std::string> dmec_opt = {"Vector", "Scalar", "Pseudovector",
                                       "Pseudoscalar"};
  // std::string dme_coupling = input.get<std::string>("dme_coupling", dmec[0]);
  std::string dmec = input.get<std::string>("dme_coupling", "Vector");
  int dmec_check = 0;
  for (const auto &option : dmec_opt) {
    dmec_check += (dmec == option) ? 1 : 0;
  }
  if (dmec_check == 0) {
    std::cerr << "\nWARNING: dm-electron coupling '" << dmec
              << "' unknown, defaulting to Vector\n";
    dmec = "Vector";
  }

  auto Etune_mult = input.get<double>("Etune_mult", 1.1);
  auto Etune_add = input.get<double>("Etune_add", 0.01);

  // dmec = (check_dmec == true) ? ;

  // Make sure h (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  // ***should change this to work for each new energy
  // double du_target = (M_PI / 20.) / std::sqrt(2. * demax);
  // auto du = wf.rgrid->du();
  // if (du > du_target) {
  //   auto new_num_points = Grid::calc_num_points_from_du(
  //       wf.rgrid->r0(), wf.rgrid->rmax(), du_target,
  //       GridType::loglinear, 4.0);
  //   auto old_num_points = wf.rgrid->num_points();
  //   // num_points = (int)new_num_points;
  //   std::cerr
  //       << "\nWARNING 118: Grid not dense enough for contimuum state with "
  //       << "ec=" << demax << "au\n";
  //   std::cerr << "You should update num_points from " << old_num_points
  //             << " --> " << new_num_points << "\n";
  //   std::cerr << "Program will continue, but may fail\n";
  // }

  // outut file name (excluding extension):
  // std::string fname = "afbe-" + wf.atomicSymbol(); // + "_" + label;
  // if (label != "")
  //   fname += "_" + label;

  std::string fname_table = "afbe_table-" + wf.atomicSymbol();
  if (label != "")
    fname_table += "_" + label;

  // Print some info to screen:
  std::cout << "\nRunning Atomic Kernal for " << wf.atom()
            << " at E = 1.1*I_njl\n";
  std::cout << "*************************************************\n";
  // std::cout << "Radial " << wf.rgrid->gridParameters() << "\n\n";

  // Output HF results:
  std::cout << "  state   k     En (au)    En (eV)   Oc.Frac.\n";
  for (const auto &phi : wf.core) {
    // double rinf = wf.rinf(phi);
    // printf("%2i)%7s %2i  %3.0f %3i  %5.0e  %11.5f %12.0f %10.2f   (%.2f)\n",
    //        i++, phi.symbol().c_str(), phi.k, rinf, phi.its(), phi.eps(),
    //        phi.en(), phi.en() * PhysConst::Hartree_invcm, phi.en() *
    //        PhysConst::Hartree_eV, phi.occ_frac());
    printf(" %7s %2i %11.5f %10.2f   [%3.2f]", phi.symbol().c_str(), phi.k,
           phi.en(), phi.en() * PhysConst::Hartree_eV, phi.occ_frac());
    if (phi.l() > max_l)
      std::cout << " (excluded from K)";
    std::cout << "\n";
  }
  //////////////////////////////////////////////////

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<float>>> AK; // float ok?
  int num_states = (int)wf.core.size();
  AK.resize(desteps, std::vector<std::vector<float>>(
                         num_states, std::vector<float>(qsteps)));

  // Store state info (each orbital) [just useful for plotting!]
  std::vector<std::string> nklst; // human-readiable state labels (easy
                                  // plotting)
  nklst.reserve(wf.core.size());
  for (auto &phi : wf.core)
    nklst.emplace_back(phi.symbol(true));

  // pre-calculate the spherical Bessel function look-up table for efficiency
  timer.start();
  std::vector<std::vector<std::vector<double>>> jLqr_f;
  AKF::sphericalBesselTable(jLqr_f, max_L, qgrid.r(), wf.rgrid->r());
  std::cout << "Time for SB table: " << timer.lap_reading_str() << "\n";

  // Calculate the AK (print to screen)
  // std::cout << "\nCalculating atomic kernal AK(q,dE):\n";
  // printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)  [N=%i]\n", demin / keV,
  //        demax / keV, demin, demax, desteps);
  // printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)  [N=%i]\n", qmin / qMeV,
  //        qmax / qMeV, qmin, qmax, qsteps);

  // Calculate K(q,E)
  timer.start();
  std::cout << "Generating AF table...\," << std::flush;
  std::vector<double> eabove;
  for (std::size_t is = 0; is < wf.core.size(); is++) {
    // Storing the energies that are used
    double dE = -Etune_mult * wf.core[is].en() + Etune_add;
    eabove.push_back(dE);

    int l = wf.core[is].l(); // lorb(is);
    if (l > max_l)
      continue;
    if (plane_wave)
      AKF::calculateKpw_nk(wf, wf.core[is], dE, jLqr_f[l], AK[0][is]);
    else
      AKF::calculateK_nk(wf, wf.core[is], max_L, dE, jLqr_f, AK[0][is], alt_akf,
                         force_rescale, subtract_self, force_orthog, dmec);
  } // END loop over bound states
  // }
  std::cout << "..done :)\n";
  std::cout << "Time for AFBE: " << timer.lap_reading_str() << "\n";

  // Write out to text file (in gnuplot friendly form)
  AKF::writeToTextFile_AFBE(fname_table, AK, nklst, qmin, qmax, eabove);
  AKF::akReadWrite_AFBE(fname_table, true, AK, nklst, qmin, qmax, eabove);
  std::cout << "Written to: " << fname_table << ".txt, and .bin";

  std::cout << "\n " << timer.reading_str() << "\n";
  return;
}

} // namespace Module
