#include "Modules/exampleModule.hpp"
#include "DiracOperator/Operators.hpp" //For E1 operator
#include "IO/UserInput.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
//
#include "IO/ChronoTimer.hpp"
#include "MBPT/StructureRad.hpp"

namespace Module {

void exampleModule(const IO::UserInputBlock &input, const Wavefunction &wf) {
  // This is an example module, designed to help you write new modules

  //****************************************************************************
  // Find core/valence energy: allows distingush core/valence states
  const auto ec_max =
      std::max_element(cbegin(wf.core), cend(wf.core), DiracSpinor::comp_en)
          ->en;
  const auto ev_min = std::min_element(cbegin(wf.valence), cend(wf.valence),
                                       DiracSpinor::comp_en)
                          ->en;
  const auto en_core = 0.5 * (ev_min + ec_max);

  const auto h = DiracOperator::E1(*wf.rgrid);

  MBPT::StructureRad sr(wf.basis, en_core);
  for (const auto &v : wf.valence) {
    for (const auto &w : wf.valence) {
      IO::ChronoTimer t("x");
      std::cout << w.symbol() << "-" << v.symbol() << " " << h.reducedME(w, v)
                << " + " << sr.srTB(&h, w, v, 0.0) << "\n";
    }
  }

  return;
  //****************************************************************************

  // In this example, we will solve a new wavefunction, assuming a different
  // nuclear charge distribution, and see the difference in the energies and E1
  // matrix elements this produces.

  // Read in some optional input options: A, rms, and type
  // "checkBlock" is optional, but helps prevent errors on input options:
  input.checkBlock({"A", "rrms", "type"});
  // If no A is specified, use 0 (i.e., poinlike)
  auto A = input.get("A", 0);
  // if rrms<0, the default value (for the given A) will be looked up
  auto rrms = input.get("rrms", -1.0);
  // Will assume a Fermi nucleus, unless A or r =0
  auto nuc_type = input.get<std::string>("type", "Fermi");
  // Set nuc. type explicitely to 'pointlike' if A=0, or r_rms = 0.0
  if (A == 0 || rrms == 0.0)
    nuc_type = "pointlike";
  if (nuc_type == "pointlike") {
    A = 0;
    rrms = 0.0;
  }

  // Use the same Grid and alpha, but different nuclear parameters (except Z)
  Wavefunction wfpt(wf.rgrid->params(),
                    {wf.Znuc(), A, nuc_type, rrms, Nuclear::default_t},
                    wf.alpha / PhysConst::alpha);

  std::cout << "\n";
  IO::print_line();
  std::cout << "Calculating finite-nuclear size corrections for \n"
            << wf.atom() << ", " << wf.nuclearParams() << "\nvs.\n"
            << wfpt.atom() << ", " << wfpt.nuclearParams() << "\n\n";

  // Note: Assume only Hartree-Fock approximation, no Breit, QED, Sigma etc.

  // Solve Hartree-Fock core for new wavefuinction:
  wfpt.hartreeFockCore("HartreeFock", 0.0, wf.coreConfiguration());
  // print the new core energies:
  wfpt.printCore();

  // Solve Hartree-Fock valence
  wfpt.hartreeFockValence(DiracSpinor::state_config(wf.valence));
  wfpt.printValence();

  // Calculate the energy shifts in atomic units, and GHz
  std::cout << "\nFinite nuclear charge energy shifts:\n";
  std::cout << "  state            au           GHz\n";
  for (auto i = 0ul; i < wf.valence.size(); ++i) {
    const auto del_e = wf.valence[i].en - wfpt.valence[i].en;
    printf("%7s  %12.5e  %12.5e\n", wf.valence[i].symbol().c_str(), del_e,
           del_e * PhysConst::Hartree_GHz);
  }

  // Now, we calculate E1 matrix elements.
  std::cout << "\nFinite nuclear charge shifts to E1 reduced MEs (no RPA):\n";

  std::cout << "                      A=" << wf.Anuc()
            << "         A=" << wfpt.Anuc() << "         Shift\n";

  // 1) Create E1 operator:
  const auto hE1 = DiracOperator::E1(*wf.rgrid);

  // 2) Loop through each pair of valence states, calc E1 matrix elements:
  for (auto a = 0ul; a < wf.valence.size(); ++a) {
    const auto &Fa = wf.valence[a];
    const auto &F0a = wfpt.valence[a]; // pointlike orbital
    for (auto b = 0ul; b < a; ++b) {
      const auto &Fb = wf.valence[b];
      const auto &F0b = wfpt.valence[b]; // pointlike orbital

      // Skip the MEs which are zero due to selection rules:
      if (hE1.isZero(Fa.k, Fb.k))
        continue;

      const auto e1 = hE1.reducedME(Fa, Fb);
      const auto e10 = hE1.reducedME(F0a, F0b);
      const auto delta = e1 - e10;

      std::cout << hE1.rme_symbol(Fa, Fb) << ": ";
      printf("%12.5e  %12.5e  %12.5e\n", e1, e10, delta);
    }
  }
}

} // namespace Module
