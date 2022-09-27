#include "DMionisation/Module_atomicKernel.hpp"
#include "DMionisation/AKF_akFunctions.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/UnitConv_conversions.hpp"
#include "Wavefunction/ContinuumOrbitals.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>

namespace Module {

//==============================================================================
void atomicKernel(const IO::InputBlock &input, const Wavefunction &wf) {
  IO::ChronoTimer timer("atomicKernel");

  input.check(
      {{"Emin", "Minimum energy transfer (dE) in keV [0.1]"},
       {"Emax", "Maximum energy transfer (dE) in keV [Emin]"},
       {"Esteps", "Numer of steps along dE grid (logarithmic grid) [1]"},
       {"qmin", "Minimum momentum transfer (q) in MeV [0.01]"},
       {"qmax", "Maximum momentum transfer (q) in MeV [qmin]"},
       {"qsteps", "Number of steps along q grid (logarithmic grid) [1]"},
       {"max_l_bound",
        "Max. orbital ang. mom. for bound states to include (only for tests)"},
       {"max_L", "Maximum multipolarity used in exp(iqr) expansion [6]"},
       {"label", "label for output files"},
       {"output_text", "true/false (output K(dE,q) to .txt)"},
       {"output_binary", "true/false (output K(dE,q) to .bin)"},
       {"use_alt_akf", "Replace e^(iqr) -> e^(iqr)-1"},
       {"force_rescale", "Rescale V(r) when solving cntm orbitals [false]"},
       {"subtract_self", "Subtract Hartree-Fock self-interaction [true]"},
       {"force_orthog", "Force orthogonality of cntm orbitals [true]"},
       {"dme_coupling",
        "Vector, Scalar (g0), Pseudovector (g5), Pseudoscalar (g0g5) [Vector]"},
       {"use_Zeff_cont", "use Zeff model for continuum states [false]"}});
  if (input.has_option("help"))
    return;

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

  // read in max l and L
  const auto max_l_core = DiracSpinor::max_l(wf.core());
  auto max_l = input.get<int>("max_l_bound", max_l_core);
  if (max_l < 0 || max_l > max_l_core)
    max_l = max_l_core;
  auto max_L = input.get<int>("max_L", 6);

  const auto label = input.get<std::string>("label", "");

  // output format
  const auto text_out = input.get<bool>("output_text", true);
  const auto bin_out = input.get<bool>("output_binary", true);

  // if alt_akf then exp(iqr) -> exp(iqr) - 1 (i.e. j_L -> j_L - 1)
  const auto alt_akf = input.get<bool>("use_alt_akf", false);
  // Options used by solveContinuumHF (called in AKF)
  const auto force_rescale = input.get<bool>("force_rescale", false);
  const auto subtract_self = input.get<bool>("subtract_self", true);
  const auto force_orthog = input.get<bool>("force_orthog", true);
  const auto zeff_cont = input.get<bool>("use_Zeff_cont", false);

  // Make sure du (large-r step size) is small enough to
  // calculate (normalise) cntm functions with energy = demax
  const double du_target = (M_PI / 15.0) / std::sqrt(2.0 * demax);
  const auto du = wf.grid().du();
  if (du > du_target) {
    const auto new_num_points = Grid::calc_num_points_from_du(
        wf.grid().r0(), wf.grid().rmax(), du_target, GridType::loglinear,
        wf.grid().loglin_b());
    const auto old_num_points = wf.grid().num_points();
    // num_points = (int)new_num_points;
    std::cout
        << "\nWARNING 118: Grid not dense enough for contimuum state with "
        << "ec=" << demax << "au\n";
    std::cout << "You should update num_points from " << old_num_points
              << " --> " << new_num_points << "\n";
    std::cout << "Program will continue, but may fail\n";
  }

  //----------------------------------------------------------------------------

  // Print some info to screen:
  std::cout << "\nRunning Atomic Kernel for " << wf.atom() << "\n";
  std::cout << "=================================================\n";

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

  //----------------------------------------------------------------------------

  // DM-electron couplings
  const std::vector<std::string> dmec_opt = {"Vector", "Scalar", "Pseudovector",
                                             "Pseudoscalar"};
  const std::string dmec = input.get<std::string>("dme_coupling", "Vector");
  // check this option exists:
  const bool dmec_ok =
      std::find(dmec_opt.begin(), dmec_opt.end(), dmec) != dmec_opt.end();
  if (!dmec_ok) {
    std::cerr << "\nWARNING: dm-electron coupling '" << dmec
              << "' unknown, defaulting to Vector\n";
  } else {
    std::cout << "\nUsing " << dmec << " electron coupling:\n";
  }

  const DiracOperator::jL jl_v(wf.grid(), qgrid, std::size_t(max_L));
  const DiracOperator::g0jL jl_s(jl_v);
  const DiracOperator::ig5jL jl_pv(jl_v);
  const DiracOperator::ig0g5jL jl_ps(jl_v);
  // Use polymorphism to set correct operator
  const DiracOperator::jL &jl = dmec == "Vector" ?
                                    jl_v :
                                    dmec == "Scalar" ?
                                    jl_s :
                                    dmec == "Pseudovector" ?
                                    jl_pv :
                                    dmec == "Pseudoscalar" ? jl_ps : jl_v;
  std::cout << jl.name() << "\n";

  // Calculate the AK (print to screen)
  std::cout << "\nCalculating atomic kernel AK(q,dE):\n";
  printf(" dE: %5.2f -- %5.1f keV  (%.2f -- %.1f au)  [N=%i]\n",
         demin * UnitConv::Energy_au_to_keV, demax * UnitConv::Energy_au_to_keV,
         demin, demax, (int)desteps);
  printf("  q: %5.0e -- %5.1g MeV  (%.2f -- %.1f au)  [N=%i]\n",
         qmin * UnitConv::Momentum_au_to_MeV,
         qmax * UnitConv::Momentum_au_to_MeV, qmin, qmax, (int)qsteps);

  //----------------------------------------------------------------------------

  // Arrays to store results for outputting later:
  std::vector<std::vector<std::vector<double>>> K_E_n_q;
  const auto num_states = wf.core().size();
  K_E_n_q.resize(desteps, std::vector<std::vector<double>>(num_states));

  // Store state info (each orbital) [just useful for plotting!]
  // human-readiable state labels (easy plotting)
  std::vector<std::string> nklst;
  nklst.reserve(wf.core().size());
  for (auto &Fa : wf.core())
    nklst.emplace_back(Fa.symbol(true));

  // Calculate K_n(q,E) = K(q,n,E) [store core state seperately]
  std::cout << "Running dE loops (" << desteps << ").." << std::flush;
#pragma omp parallel for
  for (std::size_t ide = 0; ide < desteps; ide++) {
    const double dE = Egrid.r(ide);
    // Loop over core (bound) states:
    for (std::size_t ia = 0; ia < wf.core().size(); ia++) {
      const auto &Fa = wf.core().at(ia);
      if (Fa.l() > max_l)
        continue;
      K_E_n_q[ide][ia] =
          AKF::calculateK_nk(wf, Fa, max_L, dE, jl, alt_akf, force_rescale,
                             subtract_self, force_orthog, zeff_cont);
    } // END loop over bound states
  }
  std::cout << "..done :)\n";

  // Sum over core states to form total K(E,q);
  LinAlg::Matrix K_Eq(Egrid.num_points(), qgrid.num_points());
  for (std::size_t ie = 0; ie < Egrid.num_points(); ie++) {
    for (std::size_t iq = 0; iq < qgrid.num_points(); iq++) {
      for (std::size_t ia = 0; ia < wf.core().size(); ia++) {
        K_Eq(ie, iq) += K_E_n_q[ie][ia][iq];
      }
    }
  }

  //----------------------------------------------------------------------------

  // outut file name (excluding extension):
  std::string fname = "ak-" + wf.atomicSymbol();
  if (label != "")
    fname += "_" + label;

  // Write out to text file (in gnuplot friendly form)
  if (text_out) {
    AKF::write_Knk_plaintext(fname, K_E_n_q, nklst, qgrid, Egrid);
    AKF::write_Ktot_plaintext(fname, K_Eq, Egrid, qgrid);
  }
  // Write out AK as binary file (to be read in bu dmeXSection)
  if (bin_out) {
    AKF::akReadWrite(fname, true, K_E_n_q, nklst, qmin, qmax, demin, demax);
  }
  std::cout << "Written to: " << fname;
  if (text_out) {
    std::cout << ".txt, ";
    std::cout << "_tot.txt";
  }
  if (text_out && bin_out)
    std::cout << ", and ";
  if (bin_out)
    std::cout << ".bin";
  std::cout << "\n";
}

} // namespace Module
