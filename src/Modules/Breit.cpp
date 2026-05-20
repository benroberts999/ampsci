#include "HF/Breit.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include <string>

namespace Module {

// Declare, register, then define below.
void Breit(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_Breit{
  "Breit", "Breit corrections to energies at the Hartree-Fock levels", &Breit};
} // namespace

void Breit(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"units", "Energy units for output: au (default) or cm (cm^-1)"},
     {"lambda", "Frequency scaling for f-dependent Breit (default 1.0)"},
     {"second_order",
      "Do second-order (Johnson) formula. Requires basis. [false]"}});

  if (input.has_option("help")) {
    return;
  }

  if (wf.vHF()->vBreit()) {
    fmt::print("Error: This module assumes Breit is not included in "
               "wavefunction (it is included within this module!).\n");
    return;
  }

  if (wf.Sigma()) {
    fmt::print("Error: Breit module is not compatible with correlation "
               "potential (Brueckner orbitals). Run with HartreeFock (not "
               "Brueckner). To include Breit with correlations, add Breit "
               "directly to the MBPT/Sigma module and run manually.\n");
    return;
  }

  const auto units_str = input.get("units", std::string{"au"});
  const bool use_cm = (units_str == "cm");
  const double un = use_cm ? PhysConst::Hartree_invcm : 1.0;
  const std::string unit_label = use_cm ? "cm^-1" : "au";

  const double lambda = input.get("lambda", 1.0);
  const bool second_order = input.get("second_order", false);

  std::cout << "Using scale: " << lambda
            << " for f-dependent Breit (scales the frequencies)\n";

  // Static Gaunt-only, retardation-only, full static, full f-dependent
  auto G = HF::Breit();
  auto R = HF::Breit();
  const auto Br = HF::Breit();
  auto Br_f = HF::Breit();
  G.update_scale(1.0, 1.0, 1.0, 0.0, 0.0);
  R.update_scale(1.0, 0.0, 0.0, 1.0, 1.0);
  Br_f.update_lambda_f(lambda);

  const auto de1_header = [&] {
    fmt::print("      {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}  [{}]\n", "E(HF)",
               "Gaunt", "Ret.", "Freq.", "Total", unit_label);
  };
  const auto de1_row = [&](const DiracSpinor &s,
                           const std::vector<DiracSpinor> &core_states) {
    const auto eg = s * G.VbrFa(s, core_states);
    const auto er = s * R.VbrFa(s, core_states);
    const auto ef = s * Br_f.VbrFa(s, core_states);
    fmt::print("{:4s}  {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}\n",
               s.shortSymbol(), s.en() * un, eg * un, er * un,
               (ef - eg - er) * un, ef * un);
  };

  fmt::print("\nCore energy corrections (de(B1,1), without relaxation):\n");
  de1_header();
  for (const auto &a : wf.core())
    de1_row(a, wf.core());

  if (!wf.valence().empty()) {
    fmt::print(
      "\nValence energy corrections (de(B1,1), without relaxation):\n");
    de1_header();
    for (const auto &v : wf.valence())
      de1_row(v, wf.core());
  }

  // Solve HF three times: Gaunt only, Gaunt+Ret., Gaunt+Ret.+Freq.
  // Ret. and Freq. corrections extracted as differences.
  const auto core_str = wf.coreConfiguration();
  const auto val_str = DiracSpinor::state_config(wf.valence());

  fmt::print("\nSolve HF (Gaunt only)\n");
  auto wf_G = wf;
  wf_G.valence().clear();
  wf_G.solve_core("HartreeFock",
                  HF::Breit::Params{1.0, 1.0, 1.0, 0.0, 0.0, 0.0}, core_str,
                  1.0e-13, true);
  wf_G.solve_valence(val_str, true);

  fmt::print("\nSolve HF (Gaunt + Ret., static)\n");
  auto wf_GR = wf;
  wf_GR.valence().clear();
  wf_GR.solve_core("HartreeFock",
                   HF::Breit::Params{1.0, 1.0, 1.0, 1.0, 1.0, 0.0}, core_str,
                   1.0e-13, true);
  wf_GR.solve_valence(val_str, true);

  fmt::print("\nSolve HF (Gaunt + Ret. + Freq., f-dependent)\n");
  auto wf_GRF = wf;
  wf_GRF.valence().clear();
  wf_GRF.solve_core("HartreeFock",
                    HF::Breit::Params{1.0, 1.0, 1.0, 1.0, 1.0, lambda},
                    core_str, 1.0e-13, true);
  wf_GRF.solve_valence(val_str, true);

  // Relaxed tables: differences between solved wavefunctions
  const auto print_relax_table = [&](const std::vector<DiracSpinor> &orbitals) {
    de1_header();
    for (const auto &s : orbitals) {
      const auto e0 = s.en();
      const auto eG = wf_G.getState(s.symbol())->en();
      const auto eGR = wf_GR.getState(s.symbol())->en();
      const auto eGRF = wf_GRF.getState(s.symbol())->en();
      fmt::print("{:4s}  {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}\n",
                 s.shortSymbol(), e0 * un, (eG - e0) * un, (eGR - eG) * un,
                 (eGRF - eGR) * un, (eGRF - e0) * un);
    }
  };

  fmt::print("\nCore energy corrections (de(B1,1), with relaxation):\n");
  print_relax_table(wf.core());

  if (!wf.valence().empty()) {
    fmt::print("\nValence energy corrections (de(B1,1), with relaxation):\n");
    print_relax_table(wf.valence());
  }

  // Johnson formula - not that useful (excludes relaxation, which matters a lot)
  if (second_order && !wf.valence().empty() && !wf.basis().empty()) {
    fmt::print(
      "\n"
      "Using the `second-order Breit` correction as in Johnson Sec 7.4.5\n"
      "Does not fully incorperate relaxation, so should be used only for "
      "comparisons.\n"
      "Requires basis\n"
      "\n"
      "de(B1,1) = <v|Vbr|v> - first-order one-body Breit correction.\n"
      "           Gaunt and Ret. columns are static; Freq. is the additional\n"
      "           frequency-dependent correction. Total includes all three.\n"
      "de(C,2)  = <v|Sigma_C|v> - second-order Coulomb correction (normal "
      "MBPT)\n"
      "de(B1,2) = HF (one-body) part of second-order Breit correction\n"
      "de(B2,2) = Sigma (two-body) part of second-order Breit correction\n"
      "Total    = de(B1,1) + de(B1,2) + de(B2,2)\n");

    const auto eFermi = wf.FermiLevel();
    const auto holes =
      qip::select_if(wf.basis(), [=](auto &a) { return a.en() < eFermi; });
    const auto excited =
      qip::select_if(wf.basis(), [=](auto &a) { return a.en() > eFermi; });

    fmt::print("\n2nd-order energy corrections (no relaxation)  [{}]:\n",
               unit_label);
    fmt::print("      {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s}\n",
               "E(HF)", "de(C,2)", "de(B1,1)", "de(B1,2)", "de(B2,2)", "Total");

    const MBPT::StructureRad sr(wf.basis(), eFermi, {0, 999}, "", 99, {}, {},
                                false);

    for (const auto &v : wf.valence()) {
      const auto e0 = v.en();
      const auto de1 = v * Br.VbrFa(v, wf.core());
      const auto de2_C = sr.Sigma_vw(v, v);
      const auto de2_Z1 = Br.de2_HF(v, holes, excited);
      const auto de2_Z2 = Br.de2(v, holes, excited);
      const auto de2 = de2_Z1 + de2_Z2;
      fmt::print("{:4s} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} "
                 "{:10.3e}\n",
                 v.shortSymbol(), e0 * un, de2_C * un, de1 * un, de2_Z1 * un,
                 de2_Z2 * un, (de1 + de2) * un);
    }
  }
}

} // namespace Module
