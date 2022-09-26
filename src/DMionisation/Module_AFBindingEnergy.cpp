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
void AFBindingEnergy(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer; // start the overall timer

  input.check(
      {{"qmin", "[MeV] minimum momentum transfer (q) ~0.01"},
       {"qmax", "[MeV] maximum q"},
       {"qsteps", "number of steps along q grid (logarithmic grid)"},
       {"max_l_bound",
        "Max. orbital ang. mom. for bound states to include (only for tests)"},
       {"max_L", "Maximum multipolarity used in exp(iqr) expansion"},
       {"use_plane_waves", "true/false"},
       {"label", "label for output files"},
       {"output_text", "true/false (output K(dE,q) to .txt)"},
       {"use_alt_akf", " "},
       {"force_rescale", " "},
       {"subtract_self", " "},
       {"force_orthog", " "},
       {"dme_coupling", " "},
       {"Etune_mult", " "},
       {"Etune_add", " "},
       {"use_Zeff_cont", "use Zeff model for continuum states"}});

  const auto qmin_Mev = input.get<double>("qmin", 0.01);
  const auto qmax_Mev = input.get<double>("qmax", qmin_Mev);
  auto qsteps = input.get<std::size_t>("qsteps", 1);

  // Convert units for input q and dE range into atomic units
  auto qmin = qmin_Mev * UnitConv::Momentum_MeV_to_au;
  auto qmax = (qsteps == 1) ? qmin : qmax_Mev * UnitConv::Momentum_MeV_to_au;

  std::size_t desteps = 1;

  const Grid qgrid({qsteps, qmin, qmax, 0, GridType::logarithmic});

  const auto max_l_core = DiracSpinor::max_l(wf.core());
  auto max_l = input.get<int>("max_l_bound", max_l_core);
  if (max_l < 0 || max_l > max_l_core)
    max_l = max_l_core;
  auto max_L = input.get<int>("max_L", 2 * max_l); // random default..

  const bool plane_wave = input.get<bool>("use_plane_waves", false);
  if (plane_wave)
    max_L = max_l; // for spherical bessel.

  // if alt_akf then exp(iqr) -> exp(iqr) - 1 (i.e. j_L -> j_L - 1)
  auto alt_akf = input.get<bool>("use_alt_akf", false);
  // Options used by solveContinuumHF (called in AKF)
  auto force_rescale = input.get<bool>("force_rescale", false);
  auto subtract_self = input.get<bool>("subtract_self", false);
  auto force_orthog = input.get<bool>("force_orthog", false);

  auto zeff_cont = input.get<bool>("use_Zeff_cont", false);

  auto label = input.get<std::string>("label", "");

  // DM-electron couplings
  std::vector<std::string> dmec_opt = {"Vector", "Scalar", "Pseudovector",
                                       "Pseudoscalar"};
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

  // Table output file name (excluding extension)
  std::string fname_table = "afbe_table-" + wf.atomicSymbol();
  if (label != "")
    fname_table += "_" + label;

  // Print some info to screen:
  std::cout << "\nRunning Atomic Kernal for " << wf.atom()
            << " at E = 1.1*I_njl\n";
  std::cout << "*************************************************\n";

  // Output HF results:
  std::cout << "  state   k     En (au)    En (eV)   Oc.Frac.\n";
  for (const auto &phi : wf.core()) {
    printf(" %7s %2i %11.5f %10.2f   [%3.2f]", phi.symbol().c_str(),
           phi.kappa(), phi.en(), phi.en() * PhysConst::Hartree_eV,
           phi.occ_frac());
    if (phi.l() > max_l)
      std::cout << " (excluded from K)";
    std::cout << "\n";
  }
  //////////////////////////////////////////////////

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<float>>> AK; // float ok?
  const auto num_states = wf.core().size();
  AK.resize(desteps, std::vector<std::vector<float>>(num_states));

  // Store state info (each orbital) [just useful for plotting!]
  std::vector<std::string> nklst; // human-readiable state labels (easy
                                  // plotting)
  nklst.reserve(wf.core().size());
  for (auto &phi : wf.core())
    nklst.emplace_back(phi.symbol(true));

  // pre-calculate the spherical Bessel function look-up table for efficiency
  timer.start();
  const auto jLqr_f =
      AKF::sphericalBesselTable(max_L, qgrid.r(), wf.grid().r());
  std::cout << "Time for SB table: " << timer.lap_reading_str() << "\n";

  // Calculate the AK (print to screen)
  std::cout << "\nCalculating atomic kernal AK(q):\n";
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)  [N=%i]\n",
         qmin * UnitConv::Momentum_au_to_MeV,
         qmax * UnitConv::Momentum_au_to_MeV, qmin, qmax, (int)qsteps);

  // Calculate K(q,E)
  timer.start();
  std::cout << "Generating AF table...\n" << std::flush;
  std::vector<double> eabove;
  for (std::size_t is = 0; is < wf.core().size(); is++) {
    // Storing the energies that are used
    double dE = -Etune_mult * wf.core()[is].en() + Etune_add;
    eabove.push_back(dE);

    int l = wf.core()[is].l(); // lorb(is);
    if (l > max_l)
      continue;
    if (plane_wave)
      AK[0][is] = AKF::calculateKpw_nk(wf, wf.core()[is], dE, jLqr_f[l]);
    else
      AK[0][is] = AKF::calculateK_nk(wf, wf.core()[is], max_L, dE, jLqr_f,
                                     alt_akf, force_rescale, subtract_self,
                                     force_orthog, dmec, zeff_cont);
  } // END loop over bound states
  std::cout << "..done :)\n";
  std::cout << "Time for AFBE: " << timer.lap_reading_str() << "\n";

  // Write out to text file (in gnuplot friendly form)
  AKF::writeToTextFile_AFBE(fname_table, AK, nklst, qgrid, eabove);
  AKF::akReadWrite_AFBE(fname_table, true, AK, nklst, qmin, qmax, eabove);
  std::cout << "Written to: " << fname_table << ".txt, and .bin";

  std::cout << "\n " << timer.reading_str() << "\n";
  return;
}

} // namespace Module
