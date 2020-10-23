#include "MBPT/StructureRad.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "Modules/matrixElements.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace Module {

// Calculates Structure Radiation + Normalisation of States
void structureRad(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"operator", "options", "rpa", "printBoth", "onlyDiagonal",
                    "omega", "n_minmax", "splineLegs"});

  // Get input options:
  const auto oper = input.get<std::string>("operator", "E1");
  const auto h_options =
      IO::UserInputBlock(oper, input.get<std::string>("options", ""));
  const auto h = generateOperator(oper, h_options, wf, true);

  // Use spline states as diagram legs?
  const auto spline_legs = input.get("splineLegs", false);

  const auto print_both = input.get("printBoth", false);
  const auto only_diagonal = input.get("onlyDiagonal", false);
  const auto rpaQ = input.get("rpa", true);
  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto const_omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);

  // min/max n (for core/excited basis)
  const auto n_minmax = input.get_list("n_minmax", std::vector{1, 999});
  const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;

  // For freq-dependent operators:
  if (h->freqDependantQ && !eachFreqQ)
    h->updateFrequency(const_omega);

  // do RPA:
  std::unique_ptr<HF::ExternalField> dV{nullptr};
  if (rpaQ) {
    dV = std::make_unique<HF::ExternalField>(h.get(), wf.getHF());
    if (!eachFreqQ)
      dV->solve_TDHFcore(const_omega);
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

  if (wf.core.empty() || wf.valence.empty() || wf.basis.empty())
    return;

  // Find core/valence energy: allows distingush core/valence states
  const auto en_core = wf.en_coreval();

  // ----------- ** Actual Calculations ** -----------

  // Construct SR object:
  MBPT::StructureRad sr(wf.basis, en_core, {n_min, n_max});
  std::cout << std::flush;

  // Loop through all valence states, calc SR+NS
  for (const auto &v : wf.valence) {
    for (const auto &w : wf.valence) {
      if (h->isZero(w.k, v.k))
        continue;

      if (only_diagonal && w != v)
        continue;
      if (!print_both && v > w)
        continue;

      // Option to use splines (or valence states) to compute Struc Rad (use
      // splines for legs)
      const auto ws = std::find(cbegin(wf.basis), cend(wf.basis), w);
      const auto vs = std::find(cbegin(wf.basis), cend(wf.basis), v);
      if (spline_legs && (ws == cend(wf.basis) || vs == cend(wf.basis))) {
        std::cout << "Don't have requested spline for: " << w.symbol() << "-"
                  << v.symbol() << "\n";
        continue;
      }
      const auto *vp = spline_legs ? &*vs : &v;
      const auto *wp = spline_legs ? &*ws : &w;

      IO::ChronoTimer timer("time");

      std::cout << "\n" << h->rme_symbol(w, v) << ":\n";

      const auto ww = eachFreqQ ? std::abs(wp->en - vp->en) : const_omega;
      if (eachFreqQ && h->freqDependantQ) {
        h->updateFrequency(ww);
      }
      if (eachFreqQ && rpaQ) {
        if (dV->get_eps() > 1.0e-3)
          dV->clear_dPsi();
        dV->solve_TDHFcore(ww);
      }

      // Zeroth-order MEs:
      const auto twvs = h->reducedME(*ws, *vs); // splines here
      const auto twv = h->reducedME(w, v);
      const auto dvs = dV ? twvs + dV->dV(*wp, *vp) : 0.0;
      const auto dv = dV ? twv + dV->dV(w, v) : 0.0;
      printer("t(spl)", twvs, dvs);
      printer("t(val)", twv, dv);

      // "Top" + "Bottom" SR terms:
      const auto [tb, tb_dv] = sr.srTB(h.get(), *wp, *vp, ww, dV.get());
      printer("SR(TB)", tb, tb_dv);
      // "Centre" SR term:
      const auto [c, c_dv] = sr.srC(h.get(), *wp, *vp, dV.get());
      printer("SR(C)", c, c_dv);

      std::cout << "========\n";
      printer("SR", tb + c, tb_dv + c_dv);

      // "Normalisation"
      const auto [n, n_dv] = sr.norm(h.get(), *wp, *vp, dV.get());
      printer("Norm", n, n_dv);

      printer("Total", tb + c + n, tb_dv + c_dv + n_dv);
      printer("as %", 100.0 * (tb + c + n) / twvs,
              100.0 * (tb_dv + c_dv + n_dv) / dvs);
    }
  }
  std::cout << "\n";

  return;
  //****************************************************************************
}

} // namespace Module
