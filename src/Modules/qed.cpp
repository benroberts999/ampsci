#include "DiracOperator/include.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <numeric>
#include <string>
#include <vector>

namespace Module {

// Declare, register, then define below.
void QED(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_QED{"QED", "QED corrections to energies at Hartree-Fock level",
                     &QED};
} // namespace

void QED(const IO::InputBlock &input, const Wavefunction &wf) {

  // Check input options for spelling mistakes etc.:
  input.check(
    {{"rcut", "Maximum radius (au) to calculate Rad Pot for [5.0]"},
     {"scale_rN", "Scale factor for Nuclear size within QED potential. 0 for "
                  "pointlike, 1 for typical [1.0]"},
     {"scale_l", "List of doubles. Extra scaling factor for each l. e.g., "
                 "1,0,1 => include for s and d, but not for p. Normally should "
                 "be left blank. [1.0]"},
     {"units", "Energy units for output: au (default) or cm (cm^-1)"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // We should not run this module if outer wavefunction already includes QED
  if (wf.vrad() != nullptr) {
    fmt::print(
      "\nError: Do not run QED module if QED corrections are included at "
      "main level!\nModule not run\n\n");
    return;
  }

  // Just HF level: complicated to re-solve Sigma with QED - just use ampsci
  if (wf.Sigma()) {
    fmt::print(
      "\nError: QED module is not compatible with correlation potential "
      "(Brueckner orbitals). Run with HartreeFock (not Brueckner). To "
      "include QED with correlations, add QED directly to the HF "
      "potential and run manually.\n");
    return;
  }

  const auto units_in = input.get("units", std::string{"au"});
  const bool use_cm = (units_in == "cm");
  const double un = use_cm ? PhysConst::Hartree_invcm : 1.0;
  const std::string unit_label = use_cm ? "cm^-1" : "au";

  // This module calculates QED corections:
  // 1. Energy corrections, without relaxation
  // 2. Energy corrections, with relaxation

  // Read in options for QED:
  const auto rcut = input.get("rcut", 5.0);
  const auto scale_rn = input.get("scale_rN", 1.0);
  const auto r_N_au =
    std::sqrt(5.0 / 3.0) * wf.nucleus().r_rms() * scale_rn / PhysConst::aB_fm;
  const auto x_spd = input.get("scale_l", std::vector{1.0});

  fmt::print("rcut = {} au\n", rcut);
  if (scale_rn != 1.0)
    fmt::print("scale_rN = {}\n", scale_rn);
  if (x_spd != std::vector{1.0})
    fmt::print("scale_l = {}\n", fmt::join(x_spd, ", "));

  // Construct the radiative potentials (seperately, and a 'total' one)
  // nb: don't write to disk, allows easier update of, e.g., rcut
  std::cout << "\n";
  auto Ueh = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {1.0, 0.0, 0.0, 0.0, 0.0}, x_spd, true, false);
  auto SEh = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {0.0, 1.0, 0.0, 0.0, 0.0}, x_spd, true, false);
  auto SEl = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {0.0, 0.0, 1.0, 0.0, 0.0}, x_spd, true, false);
  auto SEm = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                         {0.0, 0.0, 0.0, 1.0, 0.0}, x_spd, true, false);
  auto WK = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                        {0.0, 0.0, 0.0, 0.0, 1.0}, x_spd, true, false);
  auto Vtot = QED::RadPot(wf.grid().r(), wf.Znuc(), r_N_au, rcut,
                          {1.0, 1.0, 1.0, 1.0, 1.0}, x_spd, false, false);

  //----------------------------------------------------------------------------

  // Create operators, used to find first-order energy corrections
  DiracOperator::Vrad Vu(Ueh);
  DiracOperator::Vrad Vh(SEh);
  DiracOperator::Vrad Vl(SEl);
  DiracOperator::Vrad Vm(SEm);
  DiracOperator::Vrad Vw(WK);
  DiracOperator::Vrad Vt(Vtot);

  // Include QED into Hartree-Fock:
  // We 'copy' everything over from outer wavefunction, then add QED potential.
  Wavefunction wf_u = wf;
  Wavefunction wf_h = wf;
  Wavefunction wf_l = wf;
  Wavefunction wf_m = wf;
  Wavefunction wf_w = wf;
  Wavefunction wf_t = wf;
  if (wf.vHF()) {
    wf_u.vHF()->set_Vrad(Ueh);
    wf_h.vHF()->set_Vrad(SEh);
    wf_l.vHF()->set_Vrad(SEl);
    wf_m.vHF()->set_Vrad(SEm);
    wf_w.vHF()->set_Vrad(WK);
    wf_t.vHF()->set_Vrad(Vtot);
  } else {
    fmt::print("Error: no HF potential found.\n");
    return;
  }

  fmt::print("\nRe-solving Hartree-Fock for core states, including:\n");
  fmt::print("Uehling potential:\n");
  wf_u.solve_core();
  fmt::print("Self-energy (high frequency):\n");
  wf_h.solve_core();
  fmt::print("Self-energy (low frequency):\n");
  wf_l.solve_core();
  fmt::print("Self-energy (magnetic):\n");
  wf_m.solve_core();
  fmt::print("Wichmann-Kroll (approximate):\n");
  wf_w.solve_core();
  fmt::print("Total radiative potential:\n");
  wf_t.solve_core();

  // We clear the valence wavefunctions, and re-solve them from scratch.
  wf_u.valence().clear();
  wf_h.valence().clear();
  wf_l.valence().clear();
  wf_m.valence().clear();
  wf_w.valence().clear();
  wf_t.valence().clear();

  // Solve for valence states with QED:
  // nb: hartreeFockBrueckner() has no effect if Sigma doesn't exist
  const auto valence_list = DiracSpinor::state_config(wf.valence());
  fmt::print("\nRe-solving Hartree-Fock for valence states {}, including:\n",
             valence_list);

  fmt::print("Uehling potential:\n");
  wf_u.solve_valence(valence_list);

  fmt::print("Self-energy (high frequency):\n");
  wf_h.solve_valence(valence_list);

  fmt::print("Self-energy (low frequency):\n");
  wf_l.solve_valence(valence_list);

  fmt::print("Self-energy (magnetic):\n");
  wf_m.solve_valence(valence_list);

  fmt::print("Wichmann-Kroll (approximate):\n");
  wf_w.solve_valence(valence_list);

  fmt::print("Total radiative potential:\n");
  wf_t.solve_valence(valence_list);

  fmt::print("\nCore energies (including total QED):\n");
  wf_t.printCore();
  fmt::print("\nValence energies (including total QED):\n");
  wf_t.printValence();

  //----------------------------------------------------------------------------
  // Output:

  // common table "header"
  const auto qed_header = [&]() {
    fmt::print("{:5s} {:>11s} {:>11s} {:>11s} {:>11s} {:>11s} {:>11s}  [{}]\n",
               "State", "Uehl", "SE(h)", "SE(l)", "SE(m)", "WK", "Total",
               unit_label);
  };

  // 1. First-order QED corrections to energies (no relaxation)
  const auto de1_row_printer = [&](const DiracSpinor &s) {
    fmt::print("{:5s} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n",
               s.shortSymbol(), Vu.radialIntegral(s, s) * un,
               Vh.radialIntegral(s, s) * un, Vl.radialIntegral(s, s) * un,
               Vm.radialIntegral(s, s) * un, Vw.radialIntegral(s, s) * un,
               Vt.radialIntegral(s, s) * un);
  };

  fmt::print("\nCore energy corrections (first-order, no relaxation)  [{}]:\n",
             unit_label);
  qed_header();
  for (const auto &a : wf.core()) {
    de1_row_printer(a);
  }

  if (!wf.valence().empty()) {
    fmt::print(
      "\nValence energy corrections (first-order, no relaxation)  [{}]:\n",
      unit_label);
    qed_header();
    for (const auto &v : wf.valence()) {
      de1_row_printer(v);
    }
  }

  // 2. Total QED corrections to energies (with relaxation)
  const auto de_relax_row = [&](const DiracSpinor &s0, const Wavefunction &wu,
                                const Wavefunction &wh, const Wavefunction &wl,
                                const Wavefunction &wm, const Wavefunction &ww,
                                const Wavefunction &wt) {
    const auto e0 = s0.en();
    fmt::print("{:5s} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e} {:11.4e}\n",
               s0.shortSymbol(),
               (wu.getState(s0.n(), s0.kappa())->en() - e0) * un,
               (wh.getState(s0.n(), s0.kappa())->en() - e0) * un,
               (wl.getState(s0.n(), s0.kappa())->en() - e0) * un,
               (wm.getState(s0.n(), s0.kappa())->en() - e0) * un,
               (ww.getState(s0.n(), s0.kappa())->en() - e0) * un,
               (wt.getState(s0.n(), s0.kappa())->en() - e0) * un);
  };

  fmt::print("\nCore energy corrections (with relaxation)  [{}]:\n",
             unit_label);
  qed_header();
  for (const auto &a : wf.core())
    de_relax_row(a, wf_u, wf_h, wf_l, wf_m, wf_w, wf_t);

  if (!wf.valence().empty()) {
    fmt::print("\nValence energy corrections (with relaxation)  [{}]:\n",
               unit_label);
    qed_header();
    for (const auto &v0 : wf.valence())
      de_relax_row(v0, wf_u, wf_h, wf_l, wf_m, wf_w, wf_t);
  }
}

} // namespace Module
