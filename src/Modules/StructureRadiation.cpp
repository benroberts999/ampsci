#include "Amplitudes/MatrixElements.hpp"
#include "DiracOperator/GenerateOperator.hpp"
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

void structureRadiation(const IO::InputBlock &input, const Wavefunction &wf);

namespace {
const Register r_structureRadiation{
  "StructureRadiation",
  "Structure radiation + normalisation corrections to matrix elements",
  &structureRadiation};
// Aliases for the same module:
const Register r_alias_StructureRad{
  "StructureRad", "Alias for StructureRadiation", &structureRadiation};
const Register r_alias_StrucRad{"StrucRad", "Alias for StructureRadiation",
                                &structureRadiation};
} // namespace

//==============================================================================
// The calculations are done by Amplitudes::matrix_element (with the
// frequencies set by Amplitudes::solve_at_frequency); this module only
// parses input and prints every part of each element.
void structureRadiation(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
    {{"", "Calculates structure radiation, normalisation of states, and "
          "Brueckner orbital corrections to matrix elements using "
          "perturbation theory"},
     {"operator", "e.g., E1, hfs"},
     {"options{}", "options specific to operator; blank by dflt"},
     {"rpa", "true(=TDHF), false, TDHF, basis, diagram, cntm [true]"},
     {"omega",
      "Text or number. Frequency for RPA and the SR tables/denominators. "
      "Put 'each' to solve at correct frequency for each transition. [0.0]"},
     {"omega_operator",
      "Frequency-dependent operators are evaluated at the transition "
      "frequency for each element, which is the physical frequency. Set this "
      "to pin the operator at a fixed frequency instead (rarely meaningful; "
      "mainly for comparison to older calculations)."},
     {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
     {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
     {"off-diagonal",
      "Calculate off-diagonal matrix elements (if non-zero) [true]"},
     {"what",
      "What to calculate? Options are: Reduced (reduced matric elements), "
      "Stetched (stretched states, with j=m= [j=min(ja,jb) for off-diagonal]), "
      "or HFConstant for (hyperfine A,B,etc. constants). Default is Reduced, "
      "except for hyperfine operator, for which it is HFConstant"},
     {"Qk_file",
      "true/false/filename - SR: filename for QkTable file. If blank will "
      "not use QkTable; if exists, will read it in; if doesn't exist, will "
      "create it and write to disk. If 'true' will use default filename. "
      "Save time (10x) at cost of memory. Note: Using QkTable implies "
      "legs=basis [true]"},
     {"n_minmax", "list; min,max n for core/excited: (1,inf)dflt"},
     {"k_cut",
      "Maximum multipolarity k to include in Coulomb Qk. Default: all"},
     {"include_core", "If true, includes core states in calculation. Will "
                      "use HF core, unless legs=spectrum [false]"},
     {"legs",
      "Which states to use for diagram legs? hf/basis/spectrum/brueckner "
      "[hf]"},
     {"fk", "list: Coulomb screening factors for each k"},
     {"etak", "list: Hole-particle factors for each k"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timerSR("StructureRadiation");

  // Get input options:
  const auto oper = input.get<std::string>("operator", "E1");
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  // treat hyperfine operator differently: constants instead of RME
  const bool is_hyperfine =
    qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  const auto h = DiracOperator::generate(oper, h_options, wf);

  const auto Qk_file_t = input.get("Qk_file", std::string{"true"});
  const std::string Qk_file =
    Qk_file_t != "false" ?
      Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t :
      "";

  const auto k_cut = input.get("k_cut", 99);
  if (k_cut < 99) {
    std::cout << "Cutting Qk integrals at k = " << k_cut << "\n";
  }

  const auto include_core = input.get("include_core", false);
  if (include_core) {
    std::cout << "Including core-state matrix elements\n";
  }

  const auto legs_str =
    Qk_file == "" ? input.get("legs", std::string{"hf"}) : "basis";

  // External-leg states: core (optionally), then valence, taken from the
  // requested source (hf/basis/spectrum/brueckner)
  std::vector<DiracSpinor> orbs;
  if (include_core) {
    if (qip::ci_wc_compare(legs_str, "spectrum")) {
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
  const auto &orbs2 = //
    qip::ci_wc_compare(legs_str, "basis")    ? wf.basis() :
    qip::ci_wc_compare(legs_str, "spectrum") ? wf.spectrum() :
    qip::ci_wc_compare(legs_str, "bru*")     ? wf.valence() :
                                               wf.hf_valence();
  for (const auto &v : wf.valence()) {
    const auto t = std::find(orbs2.cbegin(), orbs2.cend(), v);
    if (t != orbs2.cend()) {
      orbs.push_back(*t);
    }
  }

  // If the legs are already Brueckner orbitals, the BO correction is
  // already included in them
  const auto have_brueckner =
    wf.Sigma() && (qip::ci_wc_compare(legs_str, "spectrum") ||
                   qip::ci_wc_compare(legs_str, "bru*"));

  // effective screening factors (Coulomb)
  const auto fk = input.get("fk", std::vector<double>{});
  const auto etak = input.get("etak", std::vector<double>{});

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = qip::ci_compare(str_om, "each");
  const auto const_omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);
  const auto omega_operator = input.get<double>("omega_operator");

  if (h->freqDependantQ()) {
    std::cout << "Frequency-dependent operator; at omega = ";
    if (omega_operator)
      std::cout << *omega_operator << " (pinned)\n";
    else
      std::cout << "each transition frequency\n";
  }

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

  // min/max n (for core/excited basis)
  const auto n_minmax = input.get("n_minmax", std::vector{1, 999});
  const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;

  // RPA:
  const auto rpa_method_str = input.get("rpa", std::string("true"));
  auto dV = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                    wf.basis(), wf.identity());

  std::cout << "\nStructure radiation and normalisation of states:\n";
  std::cout << "h=" << h->name() << "\n";
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";
  if (!fk.empty()) {
    std::cout << "Including effective Coulomb screening, with: fk = ";
    for (const auto &f : fk) {
      std::cout << f << ", ";
    }
    std::cout << "\n";
  }
  if (!etak.empty()) {
    std::cout << "Including effective hole-particle interation, with: etak = ";
    for (const auto &f : etak) {
      std::cout << f << ", ";
    }
    std::cout << "\n";
  }
  std::cout << "Evaluated at ";
  if (eachFreqQ)
    std::cout << "each transition frequency\n";
  else
    std::cout << "constant frequency: w = " << const_omega << "\n";

  if (qip::ci_wc_compare(legs_str, "basis"))
    std::cout << "Using basis (splines) for diagram legs (external states)\n";
  else if (qip::ci_wc_compare(legs_str, "spectrum"))
    std::cout << "Using spectrum for diagram legs (external states)\n";
  else if (qip::ci_wc_compare(legs_str, "bru*"))
    std::cout
      << "Using HF/Brueckner states for diagram legs (external states)\n";
  else
    std::cout << "Using HF valence states for diagram legs (external states)\n";

  if (!Qk_file.empty()) {
    std::cout << "Will read/write Qk integrals to file: " << Qk_file << "\n";
  } else {
    std::cout << "Will calculate Qk integrals on-the-fly\n";
  }

  // Check orthogonality of splines to valence
  if (qip::ci_wc_compare(legs_str, "basis") ||
      qip::ci_wc_compare(legs_str, "spectrum")) {
    std::cout << "\nUsing basis/spectrum for diagram legs: Check "
                 "orthonormality:\n";
    const auto &tmp_val =
      qip::ci_wc_compare(legs_str, "basis") ? wf.hf_valence() : wf.valence();
    for (const auto &v : tmp_val) {
      const auto bp = std::find(cbegin(orbs), cend(orbs), v);
      if (bp == cend(orbs))
        continue;
      const auto eps1 = std::abs(v * *bp - 1.0);
      const auto eps2 = std::abs(v.en() - bp->en()) / std::abs(v.en());
      const auto eps = std::max(eps1, eps2);
      const std::string warn = eps > 1.0e-2 ? "***" :
                               eps > 1.0e-3 ? "**" :
                               eps > 1.0e-4 ? "*" :
                                              "";
      fmt::print("{:4s} : |<v|v'>-1| = {:.1e},  dE = {:.1e}  {}\n",
                 v.shortSymbol(), eps1, eps2, warn);
    }
  }

  if (wf.core().empty() || wf.valence().empty() || wf.basis().empty())
    return;

  // Construct SR object (holds basis, Qk integrals, screening):
  MBPT::StructureRad sr(wf.basis(), wf.FermiLevel(), {n_min, n_max}, Qk_file,
                        k_cut, fk, etak);
  std::cout << std::flush;

  // Frequency choices: the operator is pinned only if omega_operator was
  // given; the RPA, and the SR tables, are solved here unless they follow
  // each transition
  using Amplitudes::Frequency;
  Amplitudes::MEoptions options{
    omega_operator ? Frequency::fixed : Frequency::transition,
    eachFreqQ ? Frequency::transition : Frequency::fixed};
  if (omega_operator && h->freqDependantQ()) {
    h->updateFrequency(*omega_operator);
  }
  if (!eachFreqQ && dV) {
    dV->solve_core(const_omega, 100, true);
  }
  if (!eachFreqQ && (omega_operator || !h->freqDependantQ())) {
    sr.solve_core(h.get(), dV.get());
  }
  options.diagonal = input.get("diagonal", true);
  options.off_diagonal = input.get("off-diagonal", true);
  options.calculate_both = input.get("printBoth", false);
  options.type = matel_type;
  options.rpa_iterations = 100;

  // Each element, with every part printed as it is calculated. The Brueckner
  // orbital correction (with its frequency-derivative term) is kept alongside
  std::vector<Amplitudes::MEdata> mes;
  std::vector<double> bo_terms;
  std::optional<double> solved_at{};
  for (const auto [w, v] :
       Amplitudes::orbital_pairs(orbs, orbs, h.get(), options)) {
    IO::ChronoTimer timer("time");

    const auto w_ab = w->en() - v->en();
    fmt::print("\n{} -> {}: w = {:.8f}\n", v->shortSymbol(), w->shortSymbol(),
               w_ab);

    // Solved once per frequency (the diagonal elements all share w = 0)
    if (solved_at != w_ab) {
      Amplitudes::solve_at_frequency(h.get(), dV.get(), &sr, w_ab, options);
      solved_at = w_ab;
    }
    auto me = Amplitudes::matrix_element(
      *w, *v, h.get(), dV.get(), &sr, matel_type, w_ab,
      Amplitudes::sr_frequency(dV.get(), w_ab, options));
    const auto factor = me.factor;

    fmt::print("{:8s}  {:12.5e}\n", "t0", factor * me.t0);
    if (dV) {
      fmt::print("{:8s}  {:12.5e}\n", "dV", factor * me.dv);
      fmt::print("{:8s}  {:12.5e}\n", "RPA", me.value());
    }
    fmt::print("{:8s}  {:12.5e}\n", "SR", factor * me.sr);
    // Normalisation of states, of the RPA element, and of the RPA+SR one
    fmt::print("{:8s}  {:12.5e}\n", "Norm(0)", factor * me.norm);
    fmt::print("{:8s}  {:12.5e}\n", "Norm(SR)",
               factor * me.norm_factor * (me.t0 + me.dv + me.sr));

    // Brueckner orbital correction (only when the legs are not already
    // Brueckner orbitals), with the frequency-derivative term when the
    // operator is at its transition frequency
    double bo = 0.0;
    if (!have_brueckner) {
      bo = sr.BO(*w, *v);
      fmt::print("{:8s}  {:12.5e}\n", "BO", factor * bo);
      if (h->freqDependantQ() && !omega_operator && w != v) {
        const auto tderiv =
          Amplitudes::frequency_derivative_term(*w, *v, h.get(), &sr, w_ab);
        fmt::print("+ {:7s} {:12.5e}\n", "Tderiv", factor * tderiv);
        bo += tderiv;
      }
    }
    fmt::print("{:8s}  {:12.5e}\n", "Total", me.total() + factor * bo);
    std::cout << std::flush;

    mes.push_back(me);
    bo_terms.push_back(bo);
  }

  // Print summary table:
  std::cout << "\n"
            << "Structure Radiation + Normalisation of states: " << h->name()
            << "\n";
  std::cout << "Units: " << h->units() << "\n\n";
  std::cout << "           T1           SR           Norm         ";
  if (!have_brueckner)
    std::cout << "BO           ";
  std::cout << "Total\n";
  for (std::size_t i = 0; i < mes.size(); ++i) {
    const auto &m = mes[i];
    const auto bo = m.factor * bo_terms[i];
    fmt::print("{:4s} {:4s} {:+.4e}  {:+.4e}  {:+.4e}", m.a, m.b, m.value(),
               m.factor * m.sr, m.factor * m.norm);
    if (!have_brueckner) {
      fmt::print("  {:+.4e}", bo);
    }
    fmt::print("  {:+.4e}\n", m.total() + bo);
  }
  std::cout << std::endl;
}

} // namespace Module
