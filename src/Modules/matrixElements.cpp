#include "Modules/matrixElements.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/DiagramRPA.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "ExternalField/calcMatrixElements.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/StructureRad.hpp"
#include "Maths/Interpolator.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
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
  input.check({{"operator", "e.g., E1, hfs"},
               {"options{}", "options specific to operator; blank by dflt"},
               {"rpa", "true(=TDHF), false, TDHF, basis, diagram"},
               {"omega", "Text or number. Freq. for RPA. Put 'each' to solve "
                         "at correct frequency for each transition. [0.0]"},
               {"what", "What to calculate: rme (reduced ME), A (hyperfine A/B "
                        "coeficient), radial_integral. [rme]"},
               {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
               {"diagonal", "only <a|h|a> (dflt false)"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
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

  const auto default_what = oper == "hfs" ? "A" : "rme";
  const auto what = input.get<std::string>("what", default_what);
  const bool radial_int = qip::ci_wc_compare(what, "radial*");
  const bool hf_AB =
      oper == "hfs" && h->rank() <= 2 &&
      (qip::ci_wc_compare(what, "A*") || qip::ci_wc_compare(what, "B*"));
  const bool diagonal_only = hf_AB ? true : input.get("diagonal", false);

  const auto which_str = radial_int              ? "(radial integral)." :
                         hf_AB && h->rank() == 1 ? "(HFS constant A)." :
                         hf_AB && h->rank() == 2 ? "(HFS constant B)." :
                                                   "(reduced).";

  std::cout << "\n"
            << "Matrix Elements - " << which_str << " Operator: " << h->name()
            << "\n";
  std::cout << "Units: " << h->units() << "\n";

  const bool print_both = input.get("printBoth", false);

  const auto rpa_method_str = input.get("rpa", std::string("TDHF"));
  auto rpa_method =
      (rpa_method_str == "true")  ? ExternalField::method::TDHF :
      (rpa_method_str == "false") ? ExternalField::method::none :
                                    ExternalField::ParseMethod(rpa_method_str);

  if (rpa_method == ExternalField::method::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 77: ");
    fmt::print(
        "RPA method {} not found - check spelling? Defaulting to NO rpa\n",
        rpa_method_str);
    rpa_method = ExternalField::method::none;
  }

  if (wf.core().empty())
    rpa_method = ExternalField::method::none;
  const auto rpaQ = rpa_method != ExternalField::method::none;

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

  if ((h->parity() == 1) && rpa_method == ExternalField::method::TDHF) {
    std::cout << "\n\n*CAUTION*:\n RPA (TDHF method) may not work for this "
                 "operator.\n Consider using diagram or basis method\n\n";
  }

  std::unique_ptr<ExternalField::CorePolarisation> rpa{nullptr};
  if (rpaQ)
    std::cout << "Including RPA: ";
  if (rpa_method == ExternalField::method::TDHF) {
    std::cout << "TDHF method\n";
    rpa = std::make_unique<ExternalField::TDHF>(h.get(), wf.vHF());
  } else if (rpa_method == ExternalField::method::basis) {
    std::cout << "TDHF/basis method\n";
    rpa = std::make_unique<ExternalField::TDHFbasis>(h.get(), wf.vHF(),
                                                     wf.basis());
  } else if (rpa_method == ExternalField::method::diagram) {
    std::cout << "diagram method\n";
    rpa = std::make_unique<ExternalField::DiagramRPA>(h.get(), wf.basis(),
                                                      wf.vHF(), wf.identity());
  }

  const auto mes =
      ExternalField::calcMatrixElements(wf.valence(), h.get(), rpa.get(), omega,
                                        eachFreqQ, diagonal_only, print_both);

  std::cout << (rpaQ ? ExternalField::MEdata::title() :
                       ExternalField::MEdata::title_noRPA())
            << "\n";
  for (const auto &[a, b, w_ab, hab, dv] : mes) {

    const auto Fa = *wf.getState(a);
    const auto Fb = *wf.getState(b);
    const auto factor = radial_int ? 1.0 / h->angularF(Fa.kappa(), Fb.kappa()) :
                        hf_AB ? DiracOperator::Hyperfine::convert_RME_to_AB(
                                    h->rank(), Fa.kappa(), Fb.kappa()) :
                                1.0;

    printf(" %4s %4s  %8.5f  %13.6e", a.c_str(), b.c_str(), w_ab, factor * hab);
    if (dv != 0.0) {
      printf("  %13.6e", factor * (hab + dv));
    }
    std::cout << "\n";
  }
}

//============================================================================
// Calculates Structure Radiation + Normalisation of States
void structureRad(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"operator", "e.g., E1, hfs"},
       {"options{}", "options specific to operator; blank by dflt"},
       {"rpa", "true(=TDHF), false, TDHF, basis, diagram [false]"},
       {"omega", "freq. for RPA"},
       {"printBoth", "print <a|h|b> and <b|h|a> (dflt false)"},
       {"onlyDiagonal", "only <a|h|a> (dflt false)"},
       {"Qk_file",
        "filename for QkTable file. If blank will not use QkTable; if exists, "
        "will read it in; if doesn't exist, will create it and write to disk. "
        "Save time (10x) at cost of memory. Note: Using QkTable implies "
        "splineLegs=true"},
       {"n_minmax", "list; min,max n for core/excited: (1,inf)dflt"},
       {"splineLegs", "Use splines for diagram legs (false dflt)"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  // Get input options:
  const auto oper = input.get<std::string>("operator", "E1");
  // Get optional 'options' for operator
  auto h_options = IO::InputBlock(oper, {});
  const auto tmp_opt = input.getBlock("options");
  if (tmp_opt) {
    h_options = *tmp_opt;
  }

  // const auto h = generateOperator(oper, h_options, wf, true);
  const auto h = DiracOperator::generate(oper, h_options, wf);

  const auto Qk_file = input.get("Qk_file", std::string{""});
  // note: Using QkFile is ~10x faster (not including time to construct QkTable)
  // - but requires large amount of memory. Trade-off

  // Use spline states as diagram legs?
  // nb: Using QkTable imples using splines for legs
  const auto spline_legs =
      Qk_file.empty() ? input.get("splineLegs", false) : true;

  const auto print_both = input.get("printBoth", false);
  const auto only_diagonal = input.get("onlyDiagonal", false);
  // const auto rpaQ = input.get("rpa", true);
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

  // Parse method for RPA:
  const auto rpa_method_str = input.get("rpa", std::string("false"));
  auto rpa_method =
      (rpa_method_str == "true")  ? ExternalField::method::TDHF :
      (rpa_method_str == "false") ? ExternalField::method::none :
                                    ExternalField::ParseMethod(rpa_method_str);

  if (rpa_method == ExternalField::method::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 77: ");
    fmt::print(
        "RPA method {} not found - check spelling? Defaulting to NO rpa\n",
        rpa_method_str);
    rpa_method = ExternalField::method::none;
  }
  const auto rpaQ = rpa_method != ExternalField::method::none;

  // do RPA:
  std::unique_ptr<ExternalField::CorePolarisation> dV{nullptr};
  if (rpaQ)
    std::cout << "Including RPA: ";
  if (rpa_method == ExternalField::method::TDHF) {
    std::cout << "TDHF method\n";
    dV = std::make_unique<ExternalField::TDHF>(h.get(), wf.vHF());
  } else if (rpa_method == ExternalField::method::basis) {
    std::cout << "TDHF/basis method\n";
    dV = std::make_unique<ExternalField::TDHFbasis>(h.get(), wf.vHF(),
                                                    wf.basis());
  } else if (rpa_method == ExternalField::method::diagram) {
    std::cout << "diagram method\n";
    dV = std::make_unique<ExternalField::DiagramRPA>(h.get(), wf.basis(),
                                                     wf.vHF(), wf.identity());
  }

  if (rpaQ && !eachFreqQ) {
    dV->solve_core(const_omega);
  }

  std::cout << "\nStructure radiation and normalisation of states:\n";
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";
  std::cout << "h=" << h->name() << " (reduced ME)\n";
  std::cout << "Evaluated at ";
  if (eachFreqQ)
    std::cout << "each transition frequency\n";
  else
    std::cout << "constant frequency: w = " << const_omega << "\n";
  if (spline_legs)
    std::cout << "Using splines for diagram legs (external states)\n";
  else
    std::cout << "Using valence states for diagram legs (external states)\n";
  if (!Qk_file.empty()) {
    std::cout << "Will read/write Qk integrals to file: " << Qk_file << "\n";
  } else {
    std::cout << "Will calculate Qk integrals on-the-fly\n";
  }

  // Lambda to format the output
  const auto printer = [rpaQ](auto str, auto t, auto dv) {
    if (rpaQ)
      printf(" %6s: %12.5e  %12.5e\n", str, t, dv);
    else
      printf(" %6s: %12.5e\n", str, t);
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
  MBPT::StructureRad sr(wf.basis(), en_core, {n_min, n_max}, Qk_file);
  std::cout << std::flush;

  struct Output {
    std::string ab{};
    double t0{0.0};
    double t0dv{0.0};
    double sr{0.0};
    double n{0.0};
    double srdv{0.0};
    double ndv{0.0};
  };
  std::vector<Output> out;

  const auto pi = h->parity();

  for (std::size_t ib = 0; ib < wf.valence().size(); ib++) {
    const auto &v = wf.valence().at(ib);
    for (std::size_t ia = 0; ia < wf.valence().size(); ia++) {
      const auto &w = wf.valence().at(ia);
      if (h->isZero(w.kappa(), v.kappa()))
        continue;
      if (only_diagonal && v != w)
        continue;
      if (pi == -1) {
        if (!print_both && v.parity() == -1)
          continue;
      } else {
        if (!print_both && ib > ia)
          continue;
      }

      Output t_out;

      // Option to use splines (or valence states) to compute Struc Rad (use
      // splines for legs)
      const auto ws = std::find(cbegin(wf.basis()), cend(wf.basis()), w);
      const auto vs = std::find(cbegin(wf.basis()), cend(wf.basis()), v);
      if (spline_legs && (ws == cend(wf.basis()) || vs == cend(wf.basis()))) {
        std::cout << "Don't have requested spline for: " << w.symbol() << "-"
                  << v.symbol() << "\n";
        continue;
      }
      const auto *vp = spline_legs ? &*vs : &v;
      const auto *wp = spline_legs ? &*ws : &w;

      IO::ChronoTimer timer("time");

      std::cout << "\n" << h->rme_symbol(w, v) << ":\n";
      t_out.ab = w.shortSymbol() + " " + v.shortSymbol();

      // const auto ww = eachFreqQ ? std::abs(wp->en() - vp->en()) : const_omega;
      const auto ww = eachFreqQ ? wp->en() - vp->en() : const_omega;
      if (eachFreqQ && h->freqDependantQ()) {
        h->updateFrequency(ww);
      }
      if (eachFreqQ && rpaQ) {
        if (dV->get_eps() > 1.0e-3)
          dV->clear();
        dV->solve_core(std::abs(ww));
      }

      // Zeroth-order MEs:
      const auto twvs = h->reducedME(*ws, *vs); // splines here
      const auto twv = h->reducedME(w, v);
      const auto dvs = dV ? twvs + dV->dV(*wp, *vp) : 0.0;
      const auto dv = dV ? twv + dV->dV(w, v) : 0.0;
      printer("t(spl)", twvs, dvs);
      printer("t(val)", twv, dv);

      t_out.t0 = spline_legs ? twvs : twv;
      t_out.t0dv = spline_legs ? dvs : dv;

      // "Top" + "Bottom" SR terms:
      const auto [tb, tb_dv] = sr.srTB(h.get(), *wp, *vp, ww, dV.get());
      printer("SR(TB)", tb, tb_dv);
      // "Centre" SR term:
      const auto [c, c_dv] = sr.srC(h.get(), *wp, *vp, dV.get());
      printer("SR(C)", c, c_dv);

      std::cout << "========\n";
      printer("SR", tb + c, tb_dv + c_dv);

      t_out.sr = tb + c;
      t_out.srdv = tb_dv + c_dv;

      // "Normalisation"
      const auto [n, n_dv] = sr.norm(h.get(), *wp, *vp, dV.get());
      printer("Norm", n, n_dv);

      t_out.n = n;
      t_out.ndv = n_dv;

      printer("Total", tb + c + n, tb_dv + c_dv + n_dv);
      printer("as %", 100.0 * (tb + c + n) / twvs,
              100.0 * (tb_dv + c_dv + n_dv) / dvs);
      out.push_back(t_out);
    }
  }
  std::cout << "\n";
  std::cout << "Structure Radiation + Normalisation of states.\n"
               "Reduced matrix elements (au)\n"
            << "            t0         SR         Norm       SR+N     ";
  if (rpaQ)
    std::cout << " |  t0+dV      SR+dV      Norm+dV";
  std::cout << "\n";
  for (auto &[ab, t0, t0dv, sr0, n, srdv, ndv] : out) {
    printf(" %9s %10.3e %10.3e %10.3e %10.3e", ab.c_str(), t0, sr0, n, sr0 + n);
    if (rpaQ) {
      printf(" | %10.3e %10.3e %10.3e", t0dv, srdv, ndv);
    }
    std::cout << "\n";
  }

  return;
}

} // namespace Module
