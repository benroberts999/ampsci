#include "Modules/matrixElements.hpp"
#include "CI/CI.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/Sigma2.hpp"
#include "MBPT/StructureRad.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/ostream.hpp"
#include "qip/String.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace Module {

//==============================================================================
void matrixElements(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check(
      {{"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"options{}", "options specific to operator (see ampsci -o 'operator')"},
       {"rpa",
        "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram [true]"},
       {"rpa_options{}", "Block: some further options for RPA"},
       {"omega",
        "Text or number. Freq. for RPA (and freq. dependent operators). Put "
        "'each' to solve at correct frequency for each transition. [0.0]"},
       {"printBoth", "print <a|h|b> and <b|h|a> [false]"},
       {"use_spectrum",
        "If true (and spectrum available), will use spectrum for valence "
        "states [false]"},
       {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
       {"off-diagonal",
        "Calculate off-diagonal matrix elements (if non-zero) [true]"},
       {"StructureRadiation{}",
        "Options for Structure Radiation and normalisation (details below)"}});

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
       {"n_minmax", "list; min,max n for core/excited: [1,inf]"}});

  const auto t_rpa_input = input.getBlock("rpa_options");
  auto rpa_input = t_rpa_input ? *t_rpa_input : IO::InputBlock{"rpa_options"};
  if (input.has_option("help")) {
    rpa_input.add("help;");
  }
  rpa_input.check({{"eps", "Convergance goal [1.0e-10]"},
                   {"eta", "Damping factor - be carful with this [0.4]"},
                   {"max_iterations", "Maximum number of iterations. 1 should "
                                      "correspond to first-order RPA [100]"}});

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

  // treat hyperfine operator differently: A constants instead of RME
  const bool hf_AB =
      qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  const bool diagonal = input.get("diagonal", true);
  const bool off_diagonal = input.get("off-diagonal", true);

  std::cout << "\n"
            << "Matrix Elements - Operator: " << h->name() << "\n";
  if (hf_AB && h->rank() % 2 != 0) {
    std::cout << "Hyperfine A constants (magnetic type), K=" << h->rank()
              << "\n";
  } else if (hf_AB && h->rank() % 2 == 0) {
    std::cout << "Hyperfine B constants (electric type), K=" << h->rank()
              << "\n";
  } else {
    std::cout << "Reduced matrix elements\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  const bool print_both = input.get("printBoth", false);
  const auto use_spectrum =
      wf.spectrum().empty() ? false : input.get("use_spectrum", false);
  if (use_spectrum) {
    std::cout << "Using Spectrum (instead of valence) for matrix elements\n";
  }

  // RPA:
  auto rpa_method_str = input.get("rpa", std::string("true"));
  if (wf.core().empty())
    rpa_method_str = "false";
  auto rpa = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                     wf.basis(), wf.identity());

  const auto rpa_eps = rpa_input.get("eps", 1.0e-10);
  const auto rpa_its = rpa_input.get("max_iterations", 100);
  const auto rpa_eta = rpa_input.get("eta", 0.4);
  if (rpa) {
    rpa->set_eta(rpa_eta);
    rpa->eps_target() = rpa_eps;
  }

  const auto rpaQ = rpa != nullptr;

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

  // For SR+N
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

    std::cout
        << "\nIncluding Structure radiation and normalisation of states:\n";
    if (n_min > 1)
      std::cout << "Including from n = " << n_min << "\n";
    if (n_max < 999)
      std::cout << "Including to n = " << n_max << "\n";
    if (!Qk_file.empty()) {
      std::cout
          << "Will read/write Qk integrals to file: " << Qk_file
          << "\n  -- Note: means spline/basis states used for spline legs\n";
    } else {
      std::cout << "Will calculate Qk integrals on-the-fly\n";
    }
    std::cout << std::flush;

    sr = MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max},
                            Qk_file);
  }

  const auto pi = h->parity();

  // ability to use spectrum instead of valence
  std::vector<DiracSpinor> t_orbs;
  if (use_spectrum) {
    for (const auto &v : wf.valence()) {
      const auto t = std::find(wf.spectrum().cbegin(), wf.spectrum().cend(), v);
      if (t != wf.spectrum().cend()) {
        t_orbs.push_back(*t);
      }
    }
  }
  const auto &orbs = use_spectrum ? t_orbs : wf.valence();

  if (h->freqDependantQ() && !eachFreqQ) {
    h->updateFrequency(omega);
  }

  if (rpa && !eachFreqQ) {
    rpa->solve_core(omega, rpa_its);
  }

  std::stringstream os;

  //----------------------------------------------------
  // First, diagonal:
  if (diagonal && h->parity() == 1) {

    if (eachFreqQ && h->freqDependantQ()) {
      h->updateFrequency(0.0);
    }
    if (eachFreqQ && rpa) {
      if (rpa_its == 1) {
        rpa->clear();
      }
      rpa->solve_core(0.0, rpa_its);
    }

    for (const auto &a : orbs) {

      if (h->isZero(a.kappa(), a.kappa()))
        continue;

      const auto factor = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                      h->rank(), a.kappa(), a.kappa()) :
                                  1.0;

      const auto hab = h->reducedME(a, a);
      const auto dv = rpa ? rpa->dV(a, a) : 0.0;
      const auto sub_tot = factor * (hab + dv);

      const auto ww = 0.0;

      fmt::print(os, " {:4s} {:4s}  {:10.7f}  {:13.6e}", a.shortSymbol(),
                 a.shortSymbol(), ww, factor * hab);
      if (dv != 0.0) {
        fmt::print(os, "  {:13.6e}", factor * (hab + dv));
      }
      if (sr) {
        fmt::print("\n{}", a.shortSymbol());
        fmt::print("  ME : {:15.8e}", factor * hab);
        if (rpa)
          fmt::print(" + {:15.8e} = {:15.8e}", factor * dv, sub_tot);
        fmt::print("\n");
        fmt::print("    SR0 : ");
        std::cout << std::flush;
        const auto [tb, dvtb] = sr->srTB(h.get(), a, a, 0.0, rpa.get());
        fmt::print("{:15.8e} + ", factor * tb);
        std::cout << std::flush;
        const auto [c, dvc] = sr->srC(h.get(), a, a, rpa.get());
        fmt::print("{:15.8e} + ", factor * c);
        std::cout << std::flush;
        const auto [n, dvn] = sr->norm(h.get(), a, a, rpa.get());
        const auto sr0 = (tb + c + n) * factor;
        fmt::print("{:15.8e} = {:15.8e}\n", factor * n, sr0);
        std::cout << std::flush;
        const auto sr_rpa = rpa ? (dvtb + dvc + dvn) * factor : sr0;
        if (rpa)
          fmt::print(" SR+RPA : {:15.8e} + {:15.8e} + {:15.8e} = {:15.8e}\n",
                     factor * dvtb, factor * dvc, factor * dvn, sr_rpa);
        fmt::print("  Total : {:15.8e}\n", sub_tot + sr_rpa);
        std::cout << std::flush;
        fmt::print(os, "  {:13.6e}", sub_tot + sr_rpa);
      }

      fmt::print(os, "\n");
    }
  }

  //----------------------------------------------------
  // Then, off-diagonal:
  if (off_diagonal) {
    if ((eachFreqQ && rpa) && !sr)
      std::cout << "\n";
    for (std::size_t ib = 0; ib < orbs.size(); ib++) {
      const auto &b = orbs.at(ib);
      for (std::size_t ia = 0; ia < orbs.size(); ia++) {
        const auto &a = orbs.at(ia);

        if (a == b)
          continue;
        if (h->isZero(a.kappa(), b.kappa()))
          continue;

        // Ensure even-parity state on right for odd-parity operators
        if (pi == -1) {
          if (!print_both && b.parity() == -1)
            continue;
        } else {
          if (!print_both && ib > ia)
            continue;
        }

        const auto ww = eachFreqQ ? std::abs(a.en() - b.en()) : 0.0;
        const auto ww_s = a.en() - b.en();

        if (sr)
          std::cout << "\n";
        if ((eachFreqQ && rpa) || sr)
          fmt::print("{} - {} : w = {:.8f}\n", a.shortSymbol(), b.shortSymbol(),
                     ww_s);

        if (eachFreqQ && h->freqDependantQ()) {
          h->updateFrequency(ww);
        }
        if (eachFreqQ && rpa) {
          if (rpa->last_eps() > 1.0e-5 || rpa_its == 1 ||
              std::isnan(rpa->last_eps()))
            rpa->clear();
          std::cout << " RPA(w) : ";
          rpa->solve_core(ww, rpa_its);
        }

        const auto factor = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                        h->rank(), a.kappa(), b.kappa()) :
                                    1.0;

        const auto hab = h->reducedME(a, b);
        const auto dv = rpa ? rpa->dV(a, b) : 0.0;
        const auto sub_tot = factor * (hab + dv);

        fmt::print(os, " {:4s} {:4s}  {:10.7f}  {:13.6e}", a.shortSymbol(),
                   b.shortSymbol(), ww_s, factor * hab);
        if (dv != 0.0) {
          fmt::print(os, "  {:13.6e}", factor * (hab + dv));
        }
        if (sr) {
          fmt::print("     ME : {:15.8e}", factor * hab);
          if (rpa)
            fmt::print(" + {:15.8e} = {:15.8e}", factor * dv, sub_tot);
          fmt::print("\n");
          fmt::print("    SR0 : ");
          std::cout << std::flush;
          const auto [tb, dvtb] = sr->srTB(h.get(), a, b, ww, rpa.get());
          fmt::print("{:15.8e} + ", factor * tb);
          std::cout << std::flush;
          const auto [c, dvc] = sr->srC(h.get(), a, b, rpa.get());
          fmt::print("{:15.8e} + ", factor * c);
          std::cout << std::flush;
          const auto [n, dvn] = sr->norm(h.get(), a, b, rpa.get());
          const auto sr0 = (tb + c + n) * factor;
          std::cout << std::flush;
          fmt::print("{:15.8e} = {:15.8e}\n", factor * n, sr0);
          const auto sr_rpa = rpa ? (dvtb + dvc + dvn) * factor : sr0;
          if (rpa)
            fmt::print(" SR+RPA : {:15.8e} + {:15.8e} + {:15.8e} = {:15.8e}\n",
                       factor * dvtb, factor * dvc, factor * dvn, sr_rpa);
          fmt::print("  Total : {:15.8e}\n", sub_tot + sr_rpa);
          std::cout << std::flush;
          fmt::print(os, "  {:13.6e}", sub_tot + sr_rpa);
        }
        fmt::print(os, "\n");
      }
    }
  }

  std::cout << "\n" << h->name() << "\n";
  std::cout << "\n   a    b    w_ab        t0_ab";
  if (rpaQ)
    std::cout << "          +RPA ";
  if (sr)
    std::cout << "          +SRN";
  std::cout << "\n";
  std::cout << os.str();
  std::cout << "\n";
}

//============================================================================
// Calculates Structure Radiation + Normalisation of States
void structureRad(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Calculates structure radiation, normalisation of states, and "
            "Brueckner orbital corrections to matrix elements using "
            "perturbation theory"},
       {"operator", "e.g., E1, hfs"},
       {"options{}", "options specific to operator; blank by dflt"},
       {"rpa", "true(=TDHF), false, TDHF, basis, diagram [true]"},
       {"omega", "freq. for RPA"},
       {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
       {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
       {"off-diagonal",
        "Calculate off-diagonal matrix elements (if non-zero) [true]"},
       {"Qk_file",
        "true/false/filename - SR: filename for QkTable file. If blank will "
        "not use QkTable; if exists, will read it in; if doesn't exist, will "
        "create it and write to disk. If 'true' will use default filename. "
        "Save time (10x) at cost of memory. Note: Using QkTable implies "
        "legs=basis"},
       {"n_minmax", "list; min,max n for core/excited: (1,inf)dflt"},
       {"legs",
        "Which states to use for diagram legs? hf/basis/spectrum/brueckner "
        "[hf]"},
       {"fk", "list: Coulomb screening factors for each k"},
       {"etak", "list: Hole-particle factors for each k"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timerSR("structureRad");

  // Get input options:
  const auto oper = input.get<std::string>("operator", "E1");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  // treat hyperfine operator differently: A constants instead of RME
  const bool hf_AB =
      qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  // const auto h = generateOperator(oper, h_options, wf, true);
  const auto h = DiracOperator::generate(oper, h_options, wf);

  const auto Qk_file_t = input.get("Qk_file", std::string{"false"});
  std::string Qk_file =
      Qk_file_t != "false" ?
          Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t :
          "";
  // note: Using QkFile is ~10x faster (not including time to construct QkTable)
  // - but requires large amount of memory. Trade-off

  const auto legs_str =
      Qk_file == "" ? input.get("legs", std::string{"hf"}) : "basis";

  const auto &valence = qip::ci_wc_compare(legs_str, "basis") ? wf.basis() :
                        qip::ci_wc_compare(legs_str, "spectrum") ?
                                                                wf.spectrum() :
                        qip::ci_wc_compare(legs_str, "bru*") ? wf.valence() :
                                                               wf.hf_valence();

  const auto have_brueckner =
      wf.Sigma() && (qip::ci_wc_compare(legs_str, "spectrum") ||
                     qip::ci_wc_compare(legs_str, "bru*"));

  // effective screening factors (Coulomb)
  const auto fk = input.get("fk", std::vector<double>{});
  const auto etak = input.get("etak", std::vector<double>{});

  const auto print_both = input.get("printBoth", false);

  const auto diagonal = input.get("diagonal", true);
  const auto off_diagonal = input.get("off-diagonal", true);

  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto const_omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  // min/max n (for core/excited basis)
  const auto n_minmax = input.get("n_minmax", std::vector{1, 999});
  const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;

  // For freq-dependent operators:
  if (h->freqDependantQ() && !eachFreqQ)
    h->updateFrequency(const_omega);

  // RPA:
  const auto rpa_method_str = input.get("rpa", std::string("true"));
  auto dV = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                    wf.basis(), wf.identity());

  const auto rpaQ = dV != nullptr;

  if (rpaQ && !eachFreqQ) {
    dV->solve_core(const_omega);
  }

  std::cout << "\nStructure radiation and normalisation of states:\n";
  std::cout << "h=" << h->name() << "\n";
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";

  if (!fk.empty()) {
    std::cout << "Including effective Coulomb screening, with: fk = ";
    for (auto &f : fk) {
      std::cout << f << ", ";
    }
    std::cout << "\n";
  }
  if (!etak.empty()) {
    std::cout << "Including effective hole-particle interation, with: etak = ";
    for (auto &f : etak) {
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
  std::vector<std::string> v_warnings;
  if (qip::ci_wc_compare(legs_str, "basis") ||
      qip::ci_wc_compare(legs_str, "spectrum")) {
    std::cout << "\nUsing basis/spectrum for diagram legs: Check "
                 "orthonormality:\n";

    const auto &tmp_val =
        qip::ci_wc_compare(legs_str, "basis") ? wf.hf_valence() : wf.valence();
    for (const auto &v : tmp_val) {
      const auto bp = std::find(cbegin(valence), cend(valence), v);
      if (bp == cend(valence))
        continue;
      const auto &b = *bp;
      const auto eps1 = std::abs(v * b - 1.0);
      const auto eps2 = std::abs(v.en() - b.en()) / std::abs(v.en());
      const auto eps = std::max(eps1, eps2);

      std::string warn = eps > 1.0e-2 ? "***" :
                         eps > 1.0e-3 ? "**" :
                         eps > 1.0e-4 ? "*" :
                                        "";
      fmt::print("{:4s} : |<v|v'>-1| = {:.1e},  dE = {:.1e}  {}\n",
                 v.shortSymbol(), eps1, eps2, warn);
      v_warnings.push_back(warn);
    }
  }

  // Lambda to format the output
  const auto printer = [](auto str, auto t, auto dv) {
    fmt::print("{:6s}  {:12.5e} + {:12.5e} = {:12.5e}\n", str, t, dv - t, dv);
    std::cout << std::flush;
    // nb: the 'flush' is to force a cout flush; particularly when piping
    // output to a file, this wasn't happening soon enough
  };

  if (wf.core().empty() || wf.valence().empty() || wf.basis().empty())
    return;

  // Find core/valence energy: allows distingush core/valence states
  const auto en_core = wf.FermiLevel();

  // ----------- ** Actual Calculations ** -----------

  // Construct SR object:
  MBPT::StructureRad sr(wf.basis(), en_core, {n_min, n_max}, Qk_file, fk, etak);
  std::cout << std::flush;

  const auto pi = h->parity();
  std::stringstream os;

  for (const auto diag : {true, false}) {

    if (!diagonal && diag)
      continue;
    if (!off_diagonal && !diag)
      continue;

    for (std::size_t ib = 0; ib < wf.valence().size(); ib++) {
      const auto &v = wf.valence().at(ib);
      for (std::size_t ia = 0; ia < wf.valence().size(); ia++) {
        const auto &w = wf.valence().at(ia);
        if (h->isZero(w.kappa(), v.kappa()))
          continue;
        if (pi == -1) {
          if (!print_both && v.parity() == -1)
            continue;
        } else {
          if (!print_both && ib > ia)
            continue;
        }
        if (diag && v != w)
          continue;
        if (!diag && v == w)
          continue;

        const auto factor = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                        h->rank(), v.kappa(), w.kappa()) :
                                    1.0;

        // Option to use splines (or valence states) to compute Struc Rad (use
        // splines for legs)
        const auto ws = std::find(cbegin(valence), cend(valence), w);
        const auto vs = std::find(cbegin(valence), cend(valence), v);

        if ((ws == cend(valence) || vs == cend(valence))) {
          std::cout << "Don't have requested spline for: " << w.symbol()
                    << " or " << v.symbol() << "\n";
          continue;
        }

        IO::ChronoTimer timer("time");

        std::cout << "\n" << w << " - " << v << ":\n";

        const auto ww = eachFreqQ ? std::abs(ws->en() - vs->en()) : const_omega;
        const auto ww_s = eachFreqQ ? ws->en() - vs->en() : const_omega;
        if (eachFreqQ && h->freqDependantQ()) {
          h->updateFrequency(ww);
        }
        if (eachFreqQ && rpaQ) {
          if (dV->last_eps() > 1.0e-3 || std::isnan(dV->last_eps()))
            dV->clear();
          dV->solve_core(ww);
        }

        // Zeroth-order MEs:
        const auto twvs = factor * h->reducedME(*ws, *vs); // splines here
        const auto dvs = twvs + (dV ? factor * dV->dV(*ws, *vs) : 0.0);
        printer("t0", twvs, dvs);

        // "Top" + "Bottom" SR terms:
        const auto [tb, tb_dv] = sr.srTB(h.get(), *ws, *vs, ww, dV.get());
        printer("SR(TB)", factor * tb, factor * tb_dv);
        // "Centre" SR term:
        const auto [c, c_dv] = sr.srC(h.get(), *ws, *vs, dV.get());
        printer("SR(C)", factor * c, factor * c_dv);

        // "Normalisation"
        const auto [n, n_dv] = sr.norm(h.get(), *ws, *vs, dV.get());
        printer("Norm", factor * n, factor * n_dv);

        double T_bo = 0.0, T_bo_dv = 0.0;
        if (!have_brueckner) {
          if (eachFreqQ && h->freqDependantQ() && *ws != *vs) {
            const auto de2_w = sr.Sigma_vw(*ws, *ws);
            const auto de2_v = sr.Sigma_vw(*vs, *vs);
            const auto dw = (de2_w - de2_v) * (ww_s / ww); //
            // Make sign match that for ww used in operator etc.
            const auto del = 0.05 * ww;
            h->updateFrequency(ww + del);
            // ignore frequency dependendence in RPA solution,
            // but not in RPA matrix element?
            const auto h_plus = h->reducedME(*ws, *vs);
            h->updateFrequency(ww - del);
            const auto h_minus = h->reducedME(*ws, *vs);
            const auto d_h = (h_plus - h_minus) / (2 * del);
            const auto T_deriv = d_h * dw;
            fmt::print(" {:6s} {:12.5e}\n", "d(ww)", dw);
            fmt::print(" {:6s} {:12.5e}  (*added to BO)\n", "Tderiv",
                       factor * T_deriv);
            h->updateFrequency(ww);
            T_bo_dv += T_deriv;
            T_bo += T_deriv;
          }

          auto [bo, dvbo] = sr.BO(h.get(), *ws, *vs, dV.get());
          printer("BO", factor * (bo), factor * (dvbo));
          T_bo_dv += dvbo;
          T_bo += bo;
        }

        const auto total = dvs + factor * (tb_dv + c_dv + n_dv + T_bo_dv);

        printer("Total", twvs + factor * (tb + c + n + T_bo), total);

        fmt::print(os,
                   "{:4s} {:4s} {:+.3e} {:+.3e} {:+.3e} {:+.3e} "
                   "{:+.3e}",
                   w.shortSymbol(), v.shortSymbol(), dvs, factor * (tb + c),
                   factor * (tb_dv + c_dv), factor * (n), factor * n_dv);
        if (!have_brueckner)
          fmt::print(os, " {:+.3e} {:+.3e}", factor * T_bo_dv, total);
        fmt::print(os, "\n");
      }
    }
  }

  std::cout << "\n";

  std::cout << "\n"
            << "Structure Radiation + Normalisation of states: " << h->name()
            << "\n";
  if (hf_AB && h->rank() % 2 != 0) {
    std::cout << "Hyperfine A constants (magnetic type), K=" << h->rank()
              << "\n";
  } else if (hf_AB && h->rank() % 2 == 0) {
    std::cout << "Hyperfine B constants (electric type), K=" << h->rank()
              << "\n";
  } else {
    std::cout << "Reduced matrix elements\n";
  }
  std::cout << "Units: " << h->units() << "\n\n";

  std::cout << "           t0+dV      SR0        SR+dV    "
               "  Norm0      Norm+dV    ";
  if (!have_brueckner)
    std::cout << "BO         Total";
  std::cout << "\n";
  std::cout << os.str() << std::endl;

  return;
}

//============================================================================
void normalisation(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check(
      {{"", "Calculates normalisation correction, via derivative of "
            "Correlation potential. Uses correlation potential from main "
            "wavefunction"},
       {"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"options{}", "options specific to operator (see ampsci -o 'operator')"},
       {"rpa",
        "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram [true]"},
       {"omega",
        "Text or number. Freq. for RPA (and freq. dependent operators). Put "
        "'each' to solve at correct frequency for each transition. [0.0]"},
       {"printBoth", "print <a|h|b> and <b|h|a> [false]"},
       {"use_basis",
        "If true, will use basis states for valence states [false]"},
       {"diagonal", "Calculate diagonal matrix elements (if non-zero) [true]"},
       {"off-diagonal",
        "Calculate off-diagonal matrix elements (if non-zero) [true]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timer("normalisation");

  const auto delta = input.get("delta", 1.0e-4);

  if (!wf.Sigma()) {
    std::cout << "Correlation potential required to run this module!\n";
    return;
  }

  std::cout << "Setup Correlation potentials to calculate derivative:\n";
  const auto &Sigma0 = *wf.Sigma();
  const auto &v0 = wf.hf_valence().front();

  auto Sigma1 = *wf.Sigma();
  Sigma1.clear();
  Sigma1.formSigma(v0.kappa(), v0.en() + delta, v0.n(), &v0);

  auto Sigma2 = Sigma1;
  Sigma2.clear();

  std::cout << "\nCalculate correlation potentials:\n";
  for (const auto &v : wf.hf_valence()) {
    Sigma1.formSigma(v.kappa(), v.en() + delta, v.n(), &v);
    Sigma2.formSigma(v.kappa(), v.en() - delta, v.n(), &v);
    const auto lambda = Sigma0.getLambda(v.kappa(), v.n());
    std::cout << "Lambda = " << lambda << "\n\n";
  }

  //-------------------------------------------

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = DiracOperator::generate(oper, h_options, wf);

  // treat hyperfine operator differently: A constants instead of RME
  const bool hf_AB =
      qip::ci_compare(oper, "hfs") || qip::ci_compare(oper, "MLVP");

  const bool diagonal = input.get("diagonal", true);
  const bool off_diagonal = input.get("off-diagonal", true);

  std::cout << "\n"
            << "Matrix Elements - Operator: " << h->name() << "\n";
  if (hf_AB && h->rank() % 2 != 0) {
    std::cout << "Hyperfine A constants (magnetic type), K=" << h->rank()
              << "\n";
  } else if (hf_AB && h->rank() % 2 == 0) {
    std::cout << "Hyperfine B constants (electric type), K=" << h->rank()
              << "\n";
  } else {
    std::cout << "Reduced matrix elements\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  const bool print_both = input.get("printBoth", false);
  const auto use_basis =
      wf.basis().empty() ? false : input.get("use_basis", false);
  if (use_basis) {
    std::cout << "Using basis (instead of valence) for matrix elements\n";
  }

  // RPA:
  auto rpa_method_str = input.get("rpa", std::string("true"));
  if (wf.core().empty())
    rpa_method_str = "false";
  auto rpa = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                     wf.basis(), wf.identity());

  const auto rpaQ = rpa != nullptr;

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

  //-------------------------------------------

  // ability to use spectrum instead of valence
  std::vector<DiracSpinor> t_orbs;
  if (use_basis) {
    for (const auto &v : wf.hf_valence()) {
      const auto t = std::find(wf.basis().cbegin(), wf.basis().cend(), v);
      if (t != wf.basis().cend()) {
        t_orbs.push_back(*t);
      }
    }
  }
  const auto &orbs = use_basis ? t_orbs : wf.hf_valence();

  if (h->freqDependantQ() && !eachFreqQ) {
    h->updateFrequency(omega);
  }

  if (rpa && !eachFreqQ) {
    rpa->solve_core(omega);
  }

  std::stringstream os;

  const auto sr = MBPT::StructureRad(wf.basis(), wf.FermiLevel());

  //----------------------------------------------------
  // First, diagonal:
  if (diagonal && h->parity() == 1) {

    if (eachFreqQ && h->freqDependantQ()) {
      h->updateFrequency(0.0);
    }
    if (eachFreqQ && rpa) {
      rpa->solve_core(0.0);
    }

    for (const auto &v : orbs) {

      if (h->isZero(v.kappa(), v.kappa()))
        continue;

      const auto factor = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                      h->rank(), v.kappa(), v.kappa()) :
                                  1.0;

      const auto lambda_v = Sigma0.getLambda(v.kappa(), v.n());

      const auto &Sv1 = *Sigma1.getSigma(v.kappa(), v.n());
      const auto &Sv2 = *Sigma2.getSigma(v.kappa(), v.n());
      const auto dSv = lambda_v * v * ((Sv1 - Sv2) * v) / (2 * delta);

      const auto h0 = factor * h->reducedME(v, v);
      const auto dV = rpaQ ? factor * rpa->dV(v, v) : 0.0;

      fmt::print(os, "{:4s} {:4s} {:+.5e} {:+.5e} {:+.5e} {:+.5e} {:+.5e}\n",
                 v.shortSymbol(), v.shortSymbol(), h0, dV, dSv, dSv,
                 (h0 + dV) * dSv);
    }
  }
  // Then, off-diagonal:
  if (off_diagonal) {

    for (std::size_t ib = 0; ib < orbs.size(); ib++) {
      const auto &w = orbs.at(ib);
      for (std::size_t ia = 0; ia < orbs.size(); ia++) {
        const auto &v = orbs.at(ia);

        if (v == w)
          continue;
        if (h->isZero(v.kappa(), w.kappa()))
          continue;
        // ensure even parity on left
        if (!print_both && w.parity() == -1)
          continue;

        const auto ww = eachFreqQ ? std::abs(v.en() - w.en()) : 0.0;

        if (eachFreqQ && h->freqDependantQ()) {
          h->updateFrequency(ww);
        }
        if (eachFreqQ && rpa) {
          if (rpa->last_eps() > 1.0e-5 || std::isnan(rpa->last_eps()))
            rpa->clear();
          std::cout << " RPA(w) : ";
          rpa->solve_core(ww);
        }

        const auto factor = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                        h->rank(), v.kappa(), w.kappa()) :
                                    1.0;

        const auto lambda_v = Sigma0.getLambda(v.kappa(), v.n());
        const auto lambda_w = Sigma0.getLambda(w.kappa(), w.n());

        const auto &Sv1 = *Sigma1.getSigma(v.kappa(), v.n());
        const auto &Sv2 = *Sigma2.getSigma(v.kappa(), v.n());
        const auto &Sw1 = *Sigma1.getSigma(w.kappa(), w.n());
        const auto &Sw2 = *Sigma2.getSigma(w.kappa(), w.n());
        const auto dSv = lambda_v * v * ((Sv1 - Sv2) * v) / (2 * delta);
        const auto dSw = lambda_w * w * ((Sw1 - Sw2) * w) / (2 * delta);

        const auto h0 = factor * h->reducedME(v, w);
        const auto dV = rpaQ ? factor * rpa->dV(v, w) : 0.0;

        fmt::print(os, "{:4s} {:4s} {:+.5e} {:+.5e} {:+.5e} {:+.5e} {:+.5e}\n",
                   v.shortSymbol(), w.shortSymbol(), h0, dV, dSv, dSw,
                   0.5 * (h0 + dV) * (dSv + dSw));
      }
    }
  }

  std::cout
      << "\na    b     t_ab         dv_ab        dΣ_a         dΣ_b         "
         "Norm\n";
  std::cout << os.str();
  std::cout << "\n";
}

//============================================================================
// Calculates matrix elements for CI wavefunctions
void CI_matrixElements(const IO::InputBlock &input, const Wavefunction &wf) {
  //
  input.check(
      {{"operator", "e.g., E1, hfs (see ampsci -o for available operators)"},
       {"options{}", "options specific to operator"},
       {"rpa", "Method used for RPA: true(=TDHF), false, TDHF, basis, diagram"},
       {"ci_basis",
        "Re-list the CI basis used in CI{} for more efficient ME calculations "
        "[20spdf]. Only matrix elements between these states included"},
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
       {"n_minmax", "list; min,max n for core/excited: [1,inf]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer t("CI_matrixElements");

  const auto oper = input.get<std::string>("operator", "");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  const auto h = DiracOperator::generate(oper, h_options, wf);

  const bool hf_AB = oper == "hfs";
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

  std::cout << "\n"
            << "CI Matrix Elements -  Operator: " << h->name() << "\n";
  if (hf_AB && h->rank() % 2 != 0) {
    std::cout << "Hyperfine A constants (magnetic type), K=" << h->rank()
              << "\n";
  } else if (hf_AB && h->rank() % 2 == 0) {
    std::cout << "Hyperfine B constants (electric type), K=" << h->rank()
              << "\n";
  } else {
    std::cout << "Reduced matrix elements\n";
  }
  std::cout << "Units: " << h->units() << "\n";

  // RPA:
  auto rpa_method_str = input.get("rpa", std::string("false"));

  if (wf.core().empty())
    rpa_method_str = "false";
  auto rpa = ExternalField::make_rpa(rpa_method_str, h.get(), wf.vHF(), true,
                                     wf.basis(), wf.identity());

  const auto basis_string =
      input.get("ci_basis", DiracSpinor::state_config(wf.basis()));
  const std::vector<DiracSpinor> ci_basis =
      CI::basis_subset(wf.basis(), basis_string, wf.coreConfiguration());

  std::optional<MBPT::StructureRad> sr;
  if (t_SR_input) {
    // min/max n (for core/excited basis)
    const auto n_minmax = SR_input.get("n_minmax", std::vector{1});
    const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
    const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
    const auto Qk_file_t = input.get("Qk_file", std::string{"false"});
    std::string Qk_file =
        Qk_file_t != "false" ?
            Qk_file_t == "true" ? wf.identity() + ".qk.abf" : Qk_file_t :
            "";

    std::cout
        << "\nIncluding Structure radiation and normalisation of states:\n";
    if (n_min > 1)
      std::cout << "Including from n = " << n_min << "\n";
    if (n_max < 999)
      std::cout << "Including to n = " << n_max << "\n";
    if (!Qk_file.empty()) {
      std::cout
          << "Will read/write Qk integrals to file: " << Qk_file
          << "\n  -- Note: means spline/basis states used for spline legs\n";
    } else {
      std::cout << "Will calculate Qk integrals on-the-fly\n";
    }
    std::cout << std::flush;

    sr = MBPT::StructureRad(wf.basis(), wf.FermiLevel(), {n_min, n_max},
                            Qk_file);
  }
  const MBPT::StructureRad *const p_sr = sr ? &*sr : nullptr;

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
  if (eachFreqQ && rpa) {
    std::cout << "Solving RPA at each frequency\n";
  }

  Coulomb::meTable<double> me_tab;
  if (!eachFreqQ || !h->freqDependantQ()) {
    std::cout << "Calculate matrix element table.." << std::flush;
    me_tab = ExternalField::me_table(ci_basis, h.get(), rpa.get(), p_sr, omega);
    std::cout << "..done\n" << std::flush;
  }

  const auto J_list =
      input.get("J", std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
  const auto J_even_list = input.get("J+", J_list);
  const auto J_odd_list = input.get("J-", J_list);
  const auto num_solutions = input.get("num_solutions", 5ul);

  auto calc_me = [&](const CI::PsiJPi &wfA, std::size_t iA,
                     const CI::PsiJPi &wfB, std::size_t iB) {
    const auto t_omega = std::abs(wfA.energy(iA) - wfB.energy(iB));

    if (eachFreqQ && h->freqDependantQ()) {
      h->updateFrequency(t_omega);
    }
    if (eachFreqQ && rpa) {
      rpa->solve_core(t_omega, 100, false);
    }
    if (eachFreqQ && h->freqDependantQ()) {
      std::cout << "Re-Calculate matrix element table.." << std::flush;
      me_tab =
          ExternalField::me_table(ci_basis, h.get(), rpa.get(), p_sr, t_omega);
      std::cout << "..done\n" << std::flush;
    }

    const auto factor = hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB_2J(
                                    h->rank(), wfA.twoJ(), wfB.twoJ()) :
                                1.0;

    const auto me = factor * CI::ReducedME(wfA, iA, wfB, iB, me_tab, h->rank(),
                                           h->parity());

    auto p1 = wfA.parity() == 1 ? '+' : '-';
    auto p2 = wfB.parity() == 1 ? '+' : '-';

    if (eachFreqQ && rpa) {
      fmt::print(
          "{}{} {:2} {:5s} {:3s} - {}{} {:2} {:5s} {:3s}  {:2} {:.0e}  {:.5f} "
          "{:12.5e}\n",
          wfA.twoJ() / 2, p1, iA, wfA.info(iA).config,
          CI::Term_Symbol((int)wfA.info(iA).L, (int)wfA.info(iA).twoS,
                          wfA.parity()),
          wfB.twoJ() / 2, p2, iB, wfB.info(iB).config,
          CI::Term_Symbol((int)wfB.info(iB).L, (int)wfB.info(iB).twoS,
                          wfB.parity()),
          rpa->last_its(), rpa->last_eps(), t_omega, me);
    } else {
      fmt::print(
          "{}{} {:2} {:5s} {:3s} - {}{} {:2} {:5s} {:3s}  {:.5f} {:12.5e}\n",
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
