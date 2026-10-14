#include "Amplitudes/MatrixElements.hpp"
#include "CI/include.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/CorePolarisation.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/StructureRad.hpp"
#include "Modules/Modules.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace Module {

/*!
  @brief Matrix elements of any operator for HF/Brueckner valence states.
  @details
  Computes matrix elements \f$ \redmatel{a}{h}{b} \f$ of any registered
  operator @p operator between valence states, with optional RPA corrections.
  Output can be reduced matrix elements, stretched-state matrix elements
  (with \f$ m = j \f$), or hyperfine constants \f$ A, B, \ldots \f$.

  Optionally:
  - Includes RPA corrections (TDHF, diagram, or basis method).
  - Includes core-state matrix elements.
  - Uses the spectrum instead of HF valence states.
  - Solves RPA at a fixed frequency, or at each transition frequency ('each').

  Frequency-dependent operators are evaluated at each pair's transition
  frequency (their physical frequency); omega_operator pins them at a fixed
  frequency instead (rarely meaningful).

  Calculations are performed by Amplitudes::matrix_elements; this module
  only parses input and prints.

  @note For hyperfine operators (hfs, MLVP), the default output is HFS
        constants rather than reduced matrix elements.
*/
void matrixElements(const IO::InputBlock &input, const Wavefunction &wf);

//------------------------------------------------------------------------------

/*!
  @brief Matrix elements of any operator between CI many-body states.
  @details
  Computes matrix elements between CI (configuration interaction) many-body
  wavefunctions. Loops over angular momentum \f$ J \f$ and parity, applying
  selection rules. Optionally includes:
  - RPA corrections (diagram method).
  - Structure radiation and normalisation corrections (via @p StructureRadiation{}).
  - Frequency-dependent operators, solved at each transition frequency or at a
    fixed value.

  @note Requires CI wavefunctions to be computed (via the CI{} module). The
        single-particle basis is taken from the stored CI integrals (see
        Wavefunction::CI_integrals), so it need not be re-listed here.
*/
void CI_matrixElements(const IO::InputBlock &input, const Wavefunction &wf);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

namespace {
const Register r_matrixElements{"matrixElements",
                                "Calculates matrix elements of any operator",
                                &matrixElements};

const Register r_CI_matrixElements{
  "CI_matrixElements",
  "Calculates matrix elements of any operator for CI wavefunctions",
  &CI_matrixElements};

} // namespace

//==============================================================================
// As Module::matrixElements, but all calculations done by
// Amplitudes::matrix_elements; this module only parses input and prints.
void matrixElements(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({
    {"", "Matrix elements of any operator for HF/Brueckner valence states. "
         "Supports RPA, diagonal and off-diagonal elements, core states, "
         "and optional use of the spectrum instead of valence states.\n"
         "Note: SR and Norm are not included here; calculated seperately in "
         "own modules"},
    {"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
    {"options{}", "options specific to operator (see ampsci -o 'operator')"},
    {"rpa",
     "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram [true]"},
    {"rpa_options{}", "Block: some further options for RPA"},
    {"omega",
     "Text or number. Frequency RPA is solved at. Put 'each' to solve at "
     "correct frequency for each transition. [0.0]"},
    {"omega_operator",
     "Frequency-dependent operators are evaluated at the transition "
     "frequency for each element, which is the physical frequency. Set this "
     "to pin the operator at a fixed frequency instead (rarely meaningful; "
     "mainly for comparison to older calculations)."},
    {"printBoth", "print <a|h|b> and <b|h|a> [false]"},
    {"include_core", "If true, includes core states in calculation. Will "
                     "use HF core, unless use_spectrum is true [false]"},
    {"use_spectrum",
     "If true (and spectrum available), will use spectrum for valence "
     "states [false]"},
    {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
    {"off-diagonal",
     "Calculate off-diagonal matrix elements (if non-zero) [true]"},
    {"what",
     "What to calculate? Options are: Reduced (reduced matric elements), "
     "Stetched (stretched states, with j=m= [j=min(ja,jb) for off-diagonal]), "
     "or HFConstant for (hyperfine A,B,etc. constants). Default is Reduced, "
     "except for hyperfine operator, for which it is HFConstant"},
  });

  const auto t_rpa_input = input.getBlock("rpa_options");
  auto rpa_input = t_rpa_input ? *t_rpa_input : IO::InputBlock{"rpa_options"};
  if (input.has_option("help")) {
    rpa_input.add("help;");
  }
  rpa_input.check({{"eps", "Convergance goal [1.0e-10]"},
                   {"eta", "Damping factor - be careful with this [0.35]"},
                   {"max_iterations", "Maximum number of iterations. 1 should "
                                      "correspond to first-order RPA [128]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timer("matrixElements");

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = DiracOperator::generate(oper, h_options, wf);

  // treat hyperfine operator differently: constants instead of RME
  const bool is_hyperfine =
    qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  std::cout << "\n"
            << "Matrix Elements - Operator: " << h->name() << "\n";

  // Determine "what" to calculate:
  const auto what_str =
    input.get<std::string>("what", is_hyperfine ? "HFConstant" : "Reduced");
  const auto matel_type = DiracOperator::parse_MatrixElementType(what_str);

  if (matel_type == DiracOperator::MatrixElementType::Reduced) {
    std::cout << "Reduced matrix elements\n";
  } else if (matel_type == DiracOperator::MatrixElementType::HFConstant) {
    const auto EM = h->rank() % 2 == 0 ? "E" : "M";
    const std::string str = "ABCDEFGHIJKLMNOP";
    const auto sk = std::size_t(h->rank());
    const auto sym = sk < str.size() ? str.at(sk - 1) : ' ';
    std::cout << "Hyperfine " << sym << " constants (" << EM << h->rank()
              << ")\n";
  } else if (matel_type == DiracOperator::MatrixElementType::Stretched) {
    std::cout << "Stretched states with m=J [J=min(j_a,j_b) for off-diag]\n";
  } else {
    fmt2::warning();
    std::cout << " - Unkown matrix element type?\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  // Determine which states to calculate for:
  const bool print_both = input.get("printBoth", false);
  const auto include_core = input.get("include_core", false);
  const auto use_spectrum =
    wf.spectrum().empty() ? false : input.get("use_spectrum", false);
  if (include_core) {
    std::cout << "Including core-state matrix elements\n";
  }
  if (use_spectrum) {
    std::cout << "Using Spectrum (instead of valence) for matrix elements\n";
  }

  // RPA:
  auto rpa_method_str = input.get("rpa", std::string("true"));
  if (wf.core().empty())
    rpa_method_str = "false";
  auto rpa = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                     wf.basis(), wf.identity());

  const auto rpa_eps = rpa_input.get("eps", 1.0e-16);
  const auto rpa_its = rpa_input.get("max_iterations", 128);
  const auto rpa_eta = rpa_input.get("eta", 0.35);
  if (rpa) {
    rpa->set_eta(rpa_eta);
    rpa->eps_target() = rpa_eps;
  }

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = qip::ci_compare(str_om, "each");
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  const auto omega_operator = input.get<double>("omega_operator");

  if (h->freqDependantQ()) {
    std::cout << "Frequency-dependent operator; at omega = ";
    if (omega_operator)
      std::cout << *omega_operator << " (pinned)\n";
    else
      std::cout << "each transition frequency\n";
  }
  if (rpa) {
    std::cout << "RPA solved at omega = ";
    if (eachFreqQ)
      std::cout << "each transition frequency\n";
    else
      std::cout << omega << "\n";
  }

  // ability to use spectrum instead of valence, and optionally include core
  std::vector<DiracSpinor> orbs;

  // Add core states to list (either from HF or from spectrum):
  if (include_core) {
    if (use_spectrum) {
      for (const auto &a : wf.core()) {
        const auto t =
          std::find(wf.spectrum().cbegin(), wf.spectrum().cend(), a);
        if (t != wf.spectrum().cend()) {
          orbs.push_back(*t);
        }
      }
    } else {
      orbs.insert(orbs.end(), wf.core().begin(), wf.core().end());
    }
  }
  // add valcence states:
  if (use_spectrum) {
    for (const auto &v : wf.valence()) {
      const auto t = std::find(wf.spectrum().cbegin(), wf.spectrum().cend(), v);
      if (t != wf.spectrum().cend()) {
        orbs.push_back(*t);
      }
    }
  } else {
    orbs.insert(orbs.end(), wf.valence().begin(), wf.valence().end());
  }

  // Frequency choices: the operator is pinned only if omega_operator was
  // given; the RPA is solved here unless it follows each transition
  using Amplitudes::Frequency;
  Amplitudes::MEoptions options{
    omega_operator ? Frequency::fixed : Frequency::transition,
    eachFreqQ ? Frequency::transition : Frequency::fixed};
  if (omega_operator) {
    h->updateFrequency(*omega_operator);
  }
  if (!eachFreqQ && rpa) {
    rpa->solve_core(omega, rpa_its, true);
  }
  options.diagonal = input.get("diagonal", true);
  options.off_diagonal = input.get("off-diagonal", true);
  options.calculate_both = print_both;
  options.type = matel_type;
  options.rpa_iterations = rpa_its;

  const auto mes =
    Amplitudes::matrix_elements(orbs, h.get(), rpa.get(), nullptr, options);

  // Print:
  std::cout << "\n" << h->name() << "\n";
  std::cout << "\n   a    b    w_ab        t0_ab";
  if (rpa)
    std::cout << "          +RPA ";
  std::cout << "\n";
  for (const auto &m : mes) {
    fmt::print(" {:4s} {:4s}  {:10.7f}  {:13.6e}", m.a, m.b, m.omega,
               m.value0());
    if (m.dv != 0.0) {
      fmt::print("  {:13.6e}", m.value());
    }
    fmt::print("\n");
  }
  std::cout << "\n";
}

//============================================================================
// Calculates Structure Radiation + Normalisation of States
//============================================================================
// Calculates matrix elements for CI wavefunctions
void CI_matrixElements(const IO::InputBlock &input, const Wavefunction &wf) {
  //
  input.check(
    {{"",
      "Matrix elements of any operator between CI many-body states. "
      "Loops over J and parity, applies selection rules, optionally includes "
      "RPA and structure radiation + normalisation corrections."},
     {"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
     {"options{}", "options specific to operator"},
     {"rpa", "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram"},
     {"omega",
      "Text or number. Freq. for RPA (and freq. dependent operators). Put "
      "'each' to solve at correct frequency for each transition. [0.0]"},
     {"J", "List of angular momentum Js to calculate matrix elements for. If "
           "blank, all available Js will be calculated. Must be integers "
           "(two-electron only)."},
     {"J+", "As above, but for EVEN CSFs only (takes precedence over J)."},
     {"J-", "As above, but for ODD CSFs (takes precedence over J)."},
     {"num_solutions", "Maximum solution number to calculate MEs for. If "
                       "blank, will calculate all."},
     {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
     {"off-diagonal",
      "Calculate off-diagonal matrix elements (if non-zero) [true]"},
     {"what",
      "What to calculate? Options are: Reduced (reduced matric elements), "
      "Stetched (stretched states, with j=m= [j=min(ja,jb) for "
      "off-diagonal]), "
      "or HFConstant for (hyperfine A,B,etc. constants). Default is Reduced, "
      "except for hyperfine operator, for which it is HFConstant"},
     {"StructureRadiation{}",
      "Options for Structure Radiation and normalisation (details below)"}});

  // Check for Struc Rad
  const auto t_SR_input = input.getBlock("StructureRadiation");
  auto SR_input =
    t_SR_input ? *t_SR_input : IO::InputBlock{"StructureRadiation"};
  if (input.has_option("help")) {
    SR_input.add("help;");
  }
  SR_input.check(
    {{"", "If this block is included, SR + Normalisation "
          "corrections will be included"},
     {"Qk_file",
      "true/false/filename - SR: filename for QkTable file. If blank will "
      "not use QkTable; if exists, will read it in; if doesn't exist, will "
      "create it and write to disk. If 'true' will use default filename. "
      "Save time (10x) at cost of memory. Note: Using QkTable "
      "implies splines used for diagram legs"},
     {"n_minmax", "list; min,max n for core/excited: [1,inf]"},
     {"n_max_legs",
      "Max n of the CI basis states that SR+N is applied to. SR+N is only "
      "meaningful between physical states, and the high-n states of the CI "
      "basis are cavity states [default: maximum n for core + 3, as for "
      "cis2_basis]"},
     {"norm", "Include the normalisation of states? If false, only the "
              "structure radiation is included [true]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer t("CI_matrixElements");

  // The CI basis, as used to construct the CI solutions.
  // Only matrix elements between these states are required
  const auto &ci_basis = wf.CI_integrals().ci_basis;
  if (wf.CIwfs().empty() || ci_basis.empty()) {
    fmt2::error();
    std::cout << ": Requires CI solutions and the single-particle basis used "
                 "to construct them. Include a CI{} block\n";
    return;
  }

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = DiracOperator::generate(oper, h_options, wf);

  const bool is_hyperfine = oper == "hfs";
  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = qip::ci_compare(str_om, "each");
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  if (h->freqDependantQ()) {
    std::cout << "Frequency-dependent operator; at omega = ";
    if (eachFreqQ)
      std::cout << "each transition frequency\n";
    else
      std::cout << omega << "\n";
  }

  // Determine "what" to calculate:
  const auto what_str =
    input.get<std::string>("what", is_hyperfine ? "HFConstant" : "Reduced");
  const auto matel_type = DiracOperator::parse_MatrixElementType(what_str);

  if (matel_type == DiracOperator::MatrixElementType::Reduced) {
    std::cout << "Reduced matrix elements\n";
  } else if (matel_type == DiracOperator::MatrixElementType::HFConstant) {
    const auto EM = h->rank() % 2 == 0 ? "E" : "M";
    std::cout << "Hyperfine constants " << EM << h->rank() << "\n";
  } else if (matel_type == DiracOperator::MatrixElementType::Stretched) {
    std::cout << "Stretched states with m=J [J=min(j_a,j_b) for off-diag]\n";
  } else {
    fmt2::warning();
    std::cout << " - Unkown matrix element type?\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  // RPA:
  auto rpa_method_str = input.get("rpa", std::string("false"));

  if (wf.core().empty())
    rpa_method_str = "false";
  auto rpa = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                     wf.basis(), wf.identity());

  // SR+N is only meaningful between physical states, so it is applied only to
  // the low-n part of the CI basis
  const auto sr_n_max =
    SR_input.get("n_max_legs", DiracSpinor::max_n(wf.core()) + 3);
  const auto sr_norm = SR_input.get("norm", true);

  std::optional<MBPT::StructureRad> sr;
  if (t_SR_input) {
    // min/max n (for core/excited basis)
    const auto n_minmax = SR_input.get("n_minmax", std::vector{1});
    const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
    const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
    const auto Qk_file_t = SR_input.get("Qk_file", std::string{"false"});
    std::string Qk_file =
      Qk_file_t != "false" ?
        Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t :
        "";

    std::cout << "\nIncluding Structure radiation";
    std::cout << (sr_norm ? " and normalisation of states:\n" :
                            " (no normalisation of states):\n");
    if (n_min > 1)
      std::cout << "Including from n = " << n_min << "\n";
    if (n_max < 999)
      std::cout << "Including to n = " << n_max << "\n";
    std::cout << "Applied to CI basis states with n <= " << sr_n_max << "\n";
    if (!Qk_file.empty()) {
      std::cout
        << "Will read/write Qk integrals to file: " << Qk_file
        << "\n  -- Note: means spline/basis states used for spline legs\n";
    } else {
      std::cout << "Will calculate Qk integrals on-the-fly\n";
    }
    std::cout << std::flush;

    sr =
      MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file);
  }
  MBPT::StructureRad *const p_sr = sr ? &*sr : nullptr;

  if (!eachFreqQ && h->freqDependantQ()) {
    std::cout << "Frequency-dependent operator at fixed frequency: w=" << omega
              << "\n";
    h->updateFrequency(omega);
  }
  if (eachFreqQ && h->freqDependantQ()) {
    std::cout
      << "Frequency-dependent operator at frequency of each transition\n";
  }

  if (eachFreqQ && p_sr) {
    std::cout << "Warning: SR+N will take a long time..\n";
  }

  if (!eachFreqQ && rpa) {
    std::cout << "Solving RPA at fixed frequency: w=" << omega << "\n";
    rpa->solve_core(omega, 300);
  }
  // The single-particle table follows each transition if anything in it
  // (operator, RPA, or the SR tables that contain them) changes with omega
  const bool table_each = eachFreqQ && (h->freqDependantQ() || rpa);

  // SR requires solve_core, whether or not there is anything frequency dependent
  if (sr && !table_each) {
    sr->solve_core(h.get(), rpa.get());
  }
  if (eachFreqQ && rpa) {
    std::cout << "Solving RPA at each frequency\n";
  }

  // The normalisation is applied to the single-particle matrix elements
  Coulomb::meTable<double> me_tab;
  if (!table_each) {
    std::cout << "Calculate matrix element table.." << std::flush;
    me_tab = Amplitudes::me_table(ci_basis, h.get(), rpa.get(), p_sr, omega,
                                  sr_n_max, sr_norm);
    std::cout << "..done\n" << std::flush;
  }

  const auto J_list =
    input.get("J", std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
  const auto J_even_list = input.get("J+", J_list);
  const auto J_odd_list = input.get("J-", J_list);
  const auto num_solutions = input.get("num_solutions", 5ul);

  // With omega = each, the operator, the RPA, and the SR tables follow each
  // transition (see Amplitudes::solve_at_frequency)
  using Amplitudes::Frequency;
  Amplitudes::MEoptions ci_options{Frequency::transition,
                                   Frequency::transition};
  ci_options.rpa_iterations = 100;
  ci_options.print = false;

  auto calc_me = [&](const CI::PsiJPi &wfA, std::size_t iA,
                     const CI::PsiJPi &wfB, std::size_t iB) {
    const auto t_omega = wfA.energy(iA) - wfB.energy(iB);

    if (table_each) {
      Amplitudes::solve_at_frequency(h.get(), rpa.get(), p_sr, t_omega,
                                     ci_options);
      std::cout << "Re-Calculate matrix element table.." << std::flush;
      me_tab = Amplitudes::me_table(ci_basis, h.get(), rpa.get(), p_sr, t_omega,
                                    sr_n_max, sr_norm);
      std::cout << "..done\n" << std::flush;
    }

    const auto factor = h->matel_factor(matel_type, wfA.twoJ(), wfB.twoJ());

    const auto me =
      factor * CI::ReducedME(wfA, iA, wfB, iB, me_tab, h->rank(), h->parity());

    auto p1 = wfA.parity() == 1 ? '+' : '-';
    auto p2 = wfB.parity() == 1 ? '+' : '-';

    if (eachFreqQ && rpa) {
      fmt::print(
        "{}{} {:2} {:5s} {:3s} - {}{} {:2} {:5s} {:3s}  {:2} {:.0e}  {: .5f} "
        "{: 12.5e}\n",
        wfA.twoJ() / 2, p1, iA, wfA.info(iA).config,
        CI::Term_Symbol((int)wfA.info(iA).L, (int)wfA.info(iA).twoS,
                        wfA.parity()),
        wfB.twoJ() / 2, p2, iB, wfB.info(iB).config,
        CI::Term_Symbol((int)wfB.info(iB).L, (int)wfB.info(iB).twoS,
                        wfB.parity()),
        rpa->last_its(), rpa->last_eps(), t_omega, me);
    } else {
      fmt::print(
        "{}{} {:2} {:5s} {:3s} - {}{} {:2} {:5s} {:3s}  {: .5f} {: 12.5e}\n",
        wfA.twoJ() / 2, p1, iA, wfA.info(iA).config,
        CI::Term_Symbol((int)wfA.info(iA).L, (int)wfA.info(iA).twoS,
                        wfA.parity()),
        wfB.twoJ() / 2, p2, iB, wfB.info(iB).config,
        CI::Term_Symbol((int)wfB.info(iB).L, (int)wfB.info(iB).twoS,
                        wfB.parity()),
        t_omega, me);
    }
  };

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------

  std::cout << "\n";

  auto me_calculator = [&](const auto &list1, int pi1, const auto &list2,
                           int pi2, bool diagonal) {
    for (auto ej : list1) {
      const auto wf_e = wf.CIwf(ej, pi1);
      if (wf_e == nullptr)
        continue;
      for (auto oj : list2) {
        const auto wf_o = wf.CIwf(oj, pi2);
        if (wf_o == nullptr)
          continue;

        // selection rules:
        if (!h->selectrion_rule(2 * ej, pi1, 2 * oj, pi2))
          continue;

        for (std::size_t i = 0; i < num_solutions && i < wf_e->num_solutions();
             ++i) {
          for (std::size_t j = 0;
               j < num_solutions && j < wf_o->num_solutions(); ++j) {

            if (diagonal && pi1 == pi2 && ej == oj && i == j) {
              calc_me(*wf_e, i, *wf_o, j);
            }
            if (!diagonal && (pi1 != pi2 || ej != oj || i != j)) {
              calc_me(*wf_e, i, *wf_o, j);
            }
          }
        }
      }
    }
  };

  std::cout << "<a||t||b> -- "
            << DiracOperator::parse_MatrixElementType(matel_type) << "\n";
  if (eachFreqQ && rpa) {
    std::cout << "Ja  # conf      - Jb  # conf      its eps    w_ab     t_ab\n";
  } else {
    std::cout << "Ja  # conf      - Jb  # conf       w_ab     t_ab\n";
  }

  const bool diagonal = input.get("diagonal", true);
  const bool off_diagonal = input.get("off-diagonal", true);

  // diagonal:
  if (diagonal) {
    me_calculator(J_even_list, 1, J_even_list, 1, true);
    me_calculator(J_odd_list, -1, J_odd_list, -1, true);
  }

  // off-diagonal:
  if (off_diagonal) {
    me_calculator(J_even_list, 1, J_even_list, 1, false);
    me_calculator(J_even_list, 1, J_odd_list, -1, false);
    me_calculator(J_odd_list, -1, J_even_list, 1, false);
    me_calculator(J_odd_list, -1, J_odd_list, -1, false);
  }

  std::cout << "\n";
}

} // namespace Module
