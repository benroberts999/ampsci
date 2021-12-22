#include "DMionisation/AKF_akFunctions.hpp"
#include "DMionisation/Module_atomicKernal.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>

#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace Module {

//******************************************************************************
void AFStepFunction(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer; // start the overall timer

  input.checkBlock(
      {{"Emin", "[keV] minimum energy transfer (dE) ~0.1"},
       {"Emax", "[keV] maximum dE"},
       {"Esteps", "numer of steps along dE grid (logarithmic grid)"},
       {"qmin", "[MeV] minimum momentum transfer (q) ~0.01"},
       {"qmax", "[MeV] maximum q"},
       {"qsteps", "number of steps along q grid (logarithmic grid)"},
       {"max_l_bound",
        "Max. orbital ang. mom. for bound states to include (only for tests)"},
       {"table_label", "label for atomic factor table"},
       {"output_label", "label for output files"},
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
  auto demin = demin_kev * UnitConv::Energy_keV_to_au;
  auto demax = (desteps == 1) ? demin : demax_kev * UnitConv::Energy_keV_to_au;
  auto qmin = qmin_Mev * UnitConv::Momentum_MeV_to_au;
  auto qmax = (qsteps == 1) ? qmin : qmax_Mev * UnitConv::Momentum_MeV_to_au;

  // Set up the E and q grids
  const Grid Egrid({desteps, demin, demax, 0, GridType::logarithmic});
  const Grid qgrid({qsteps, qmin, qmax, 0, GridType::logarithmic});

  const auto max_l_core = wf.maxCore_l();
  auto max_l = input.get<int>("max_l_bound", max_l_core);
  if (max_l < 0 || max_l > max_l_core)
    max_l = max_l_core;

  auto table_label = input.get<std::string>("table_label", "");
  auto output_label = input.get<std::string>("output_label", "");

  // output format
  auto text_out = input.get<bool>("output_text", false);
  auto bin_out = input.get<bool>("output_binary", false);
  if (!text_out && !bin_out)
    bin_out = true; // print message?

  // Make sure h (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  double du_target = (M_PI / 20.) / std::sqrt(2. * demax);
  auto du = wf.rgrid->du();
  if (du > du_target) {
    auto new_num_points = Grid::calc_num_points_from_du(
        wf.rgrid->r0(), wf.rgrid->rmax(), du_target, GridType::loglinear, 4.0);
    auto old_num_points = wf.rgrid->num_points();
    // num_points = (int)new_num_points;
    std::cerr
        << "\nWARNING 118: Grid not dense enough for contimuum state with "
        << "ec=" << demax << "au\n";
    std::cerr << "You should update num_points from " << old_num_points
              << " --> " << new_num_points << "\n";
    std::cerr << "Program will continue, but may fail\n";
  }

  // outut file name (excluding extension):
  std::string fname = "afsf-" + wf.atomicSymbol(); // + "_" + label;
  if (output_label != "")
    fname += "_" + output_label;

  // Print some info to screen:
  std::cout << "\nRunning Step Function Approximated Atomic Kernal for "
            << wf.atom() << "\n";
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

  //

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<float>>> AK; // float ok?
  const auto num_states = wf.core.size();
  AK.resize(desteps, std::vector<std::vector<float>>(num_states));

  // Arrays to store input K table
  std::vector<std::vector<std::vector<float>>> AFBE_table;
  AFBE_table.resize(1, std::vector<std::vector<float>>(
                           num_states, std::vector<float>(qsteps)));

  std::vector<std::string> nklst;
  nklst.reserve(wf.core.size());
  for (auto &phi : wf.core)
    nklst.emplace_back(phi.symbol(true));

  std::vector<double> eabove;
  eabove.reserve(wf.core.size());

  // Start timer
  timer.start();
  // std::cout << "Looking for atomic factor table...\n";
  std::string fname_table = "afbe_table-" + wf.atomicSymbol();
  if (table_label != "")
    fname_table += "_" + table_label;

  std::cout << "Reading atomic factor table...\n";
  AKF::akReadWrite_AFBE(fname_table, false, AFBE_table, nklst, qmin, qmax,
                        eabove);

  // Calculate the AK (print to screen)
  std::cout << "\nCalculating atomic kernal AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)  [N=%i]\n",
         demin * UnitConv::Energy_au_to_keV, demax * UnitConv::Energy_au_to_keV,
         demin, demax, (int)desteps);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)  [N=%i]\n",
         qmin * UnitConv::Momentum_au_to_MeV,
         qmax * UnitConv::Momentum_au_to_MeV, qmin, qmax, (int)qsteps);

  // Calculate K(q,E)
  timer.start();
  std::cout << "Running dE loops (" << desteps << ")..\n" << std::flush;
#pragma omp parallel for
  for (std::size_t ide = 0; ide < desteps; ide++) {
    double dE = Egrid.r()[ide];
    // Loop over core (bound) states:
    for (std::size_t is = 0; is < wf.core.size(); is++) {
      const auto &psi = wf.core[is];
      const auto l = std::size_t(psi.l());
      if ((int)l > max_l)
        continue;
      AK[ide][is] = AKF::stepK_nk(psi, dE, AFBE_table[0][is]);
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
