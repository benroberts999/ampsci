#include "Amplitudes/MatrixElements.hpp"
#include "CI/CI_Integrals.hpp"
#include "Coulomb/meTable.hpp"
#include "DiracOperator/include.hpp" //For E1 operator
#include "ExternalField/TDHF.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Modules/Modules.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/format.hpp"
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace Module {

// Declare, register, then define below.
void CI_Pol(const IO::InputBlock &input, const Wavefunction &wf);
namespace {
const Register r_CI_Pol{
  "CI_Pol", "Sum over states calculation of two valence polarisability",
  &CI_Pol};
} // namespace

void CI_Pol(const IO::InputBlock &input, const Wavefunction &wf) {
  // This module computes matrix elements and energies for two valence systems

  //these quantities are then used to compute second order amplitudes using the sum over states (SOS) method

  input.check(
    {{"",
      "Static scalar polarisability of a CI state, by sum-over-states:\n"
      "alpha = (2/3[J]) sum_n |<n||d||v>|^2/(E_n - E_v),\n"
      "where the sum runs over the CI solutions of each allowed intermediate "
      "(J,parity). Requires a CI{} block (which is where the CI options are "
      "set)."},
     {"state",
      "The CI state, as J{+,-}:index - e.g., '0+' is the lowest J=0 "
      "even-parity solution, and '1-:2' the third J=1 odd-parity one (the "
      "index counts from 0, in order of energy, and may be omitted). The "
      "form 'e0', 'o1:2' is also accepted"},
     {"J", "Angular momentum of wavefunction [deprecated: use state]"},
     {"parity", "Parity of state [deprecated: use state]"},
     {"state_number", "specify which state with angular momentum J and "
                      "parity PI is required [deprecated: use state]"},
     {"StructureRadiation{}",
      "Options for structure radiation and normalisation. If this block is "
      "included, SR+N is added to the single-particle matrix element used in "
      "the sum"}});

  // Check for Structure Radiation
  const auto t_SR_input = input.getBlock("StructureRadiation");
  auto SR_input =
    t_SR_input ? *t_SR_input : IO::InputBlock{"StructureRadiation"};
  if (input.has_option("help")) {
    SR_input.add("help;");
  }
  SR_input.check(
    {{"", "If this block is included, SR + Normalisation corrections will be "
          "included"},
     {"Qk_file", "true/false/filename - SR: filename for QkTable file. If "
                 "blank will not use QkTable; if exists, will read it in; if "
                 "doesn't exist, will create it and write to disk. If 'true' "
                 "will use default filename"},
     {"n_minmax", "list; min,max n for core/excited (internal): [1,inf]"},
     {"n_max_legs",
      "Max n of the CI basis states that SR+N is applied to. SR+N is only "
      "meaningful between physical states, and the high-n states of the CI "
      "basis are cavity states [default: max_n_core + 3]"},
     {"norm", "Include the normalisation of states? If false, only the "
              "structure radiation is included [true]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  fmt2::styled_print(fg(fmt::color::blue), "\n* Note: ");
  std::cout
    << "This module is deprecated and will be removed.\n"
    << "Instead, use the general `CI_secondOrder' module, with t=E1, s=E1. "
       "Then, set method=SOS to use the sum-over-states method.\n\n";

  using namespace std::string_literals;

  //external field operator-----------------------------------------------------
  const DiracOperator::E1 h1(wf.grid());
  const auto k1 = h1.rank();
  const auto p1 = h1.parity();
  //----------------------------------------------------------------------------

  //states
  // state = J{+,-}:index takes precedence; J/parity/state_number are deprecated
  const auto input_state = input.get("state", ""s);
  const auto level = CI::parse_level(input_state);
  if (!input_state.empty() && !level) {
    fmt2::error();
    fmt::print(": Could not parse '{}' as a CI level: give as J{{+,-}}:index, "
               "e.g., 0+ or 1-:2\n",
               input_state);
    return;
  }
  const auto twoJ = level ? level->twoJ : 2 * input.get("J", 0);
  const auto parity = level ? level->parity : input.get("parity", 1);
  const auto nv = level ? level->index : input.get("state_number", 0ul);
  //get CI wavefunction for "initial state"
  const auto wfV = wf.CIwf(twoJ, parity);
  if (wfV == nullptr || nv >= wfV->num_solutions()) {
    fmt2::error();
    std::cout << ": Could not find the requested CI state ("
              << CI::to_string({twoJ, parity, nv})
              << "). Requires a CI{} block.\n"
                 "Available (J,parity):";
    for (const auto &psi : wf.CIwfs()) {
      fmt::print(" {}{}", CI::J_string(psi.twoJ()),
                 psi.parity() == 1 ? '+' : '-');
    }
    std::cout << "\n";
    return;
  }

  // The matrix elements are only required between the CI basis states
  const auto &orbitals = wf.CI_integrals().ci_basis;

  // Structure radiation.
  // SR+N is only meaningful between physical states,
  // so it is applied only to the low-n part of the CI basis.
  const auto sr_n_max =
    SR_input.get("n_max_legs", DiracSpinor::max_n(wf.core()) + 3);
  const auto sr_norm = SR_input.get("norm", true);

  std::optional<MBPT::StructureRad> sr;
  if (t_SR_input) {
    const auto n_minmax = SR_input.get("n_minmax", std::vector{1});
    const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
    const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
    const auto Qk_file_t = SR_input.get("Qk_file", "false"s);
    const std::string Qk_file =
      Qk_file_t != "false" ?
        (Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t) :
        "";
    std::cout << "\nIncluding structure radiation";
    std::cout << (sr_norm ? " and normalisation:\n" : " (no normalisation):\n");
    fmt::print("Applied to CI basis states with n <= {}\n", sr_n_max);
    sr =
      MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file);
  }

  ExternalField::TDHF tdhf(&h1, wf.vHF());
  tdhf.solve_core(0.0);
  if (sr) {
    sr->solve_core(&h1, &tdhf);
  }

  // Store and fill table of single orbtial matrix elements
  // The normalisation of states is applied to the single-particle matrix
  // elements, as in the single-valence case
  const auto sTable = Amplitudes::me_table(
    orbitals, &h1, &tdhf, sr ? &*sr : nullptr, {}, sr_n_max, sr_norm);

  //Number of states for final allowed angular momentum and parity obtained when solving CI+MBPT

  // get allowed states for matrix elements of oprator of rank K

  //for polarisability
  double pol = 0.0;
  double pol_factor = 2.0 / (3.0 * (twoJ + 1.0));

  for (int twoJ2 = twoJ - 2 * k1; twoJ2 <= twoJ + 2 * k1; twoJ2 += 2) {
    if (twoJ2 < 0) {
      continue;
    }
    if (twoJ2 == 0 && twoJ == 0) {
      continue;
    }
    const auto wfn = wf.CIwf(twoJ2, -parity);
    if (wfn == nullptr) {
      continue;
    }
    const auto num = wfn->num_solutions();
    double conj_phase = Angular::neg1pow_2(twoJ - twoJ2);
    std::cout << "\n Matrix elements and energy denominator contributions to "
                 "polarisability \n\n\n";
    for (std::size_t i = 0; i < num; i++) {
      double CI_ME_1 = CI::ReducedME(*wfV, nv, *wfn, i, sTable, k1, p1);
      double CI_ME_2 = CI::ReducedME(*wfn, i, *wfV, nv, sTable, k1, p1);
      double deltaE = wfn->energy(i) - wfV->energy(nv);
      //get configurations
      const auto cf = (*wfV).info(nv).config;
      const auto Tf = (*wfV).info(nv).config;
      const auto cn = (*wfn).info(i).config;
      //get term symbols
      std::cout << i << " " << parity << " " << "Matrix element for states"
                << cf << "--->" << cn << " " << CI_ME_1 << "  " << CI_ME_2
                << " " << "Energy denominator " << deltaE << "\n";
      double pol_cont = CI_ME_1 * conj_phase * CI_ME_2 / deltaE;
      std::cout << "Contribution to polarisability " << pol_cont << " a_0^3\n";
      pol += pol_cont;
    }
  }
  //wfV.info(nv).config;
  const auto cf = (*wfV).info(nv).config;
  pol = pol_factor * pol;
  std::cout << "Polarisability of state with configuration " << cf
            << ", angular momentum " << CI::J_string(twoJ) << " and parity "
            << parity << " is " << pol << " a_0^3\n";
}

} // namespace Module
