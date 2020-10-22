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

//! Calculates Structure Radiation + Normalisation of States
/*!
Note: Most input options are similar to MatrixElements module:

  Module::structureRad{ operator; options; rpa; printBoth; onlyDiagonal; omega;
n_minmax;  }

n_minmax: is input as list of ints:
 * n_minmax = min,max;
 * min: minimum n for core states kept in summations
 * max: maximum n for excited states kept in summations

For explanation of the rest, see MatrixElements module.

*/
void structureRad(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"operator", "options", "rpa", "printBoth", "onlyDiagonal",
                    "omega", "n_minmax"});

  // Get input options:
  const auto oper = input.get<std::string>("operator", "E1");
  const auto h_options =
      IO::UserInputBlock(oper, input.get<std::string>("options", ""));
  const auto h = generateOperator(oper, h_options, wf, true);

  const auto rpaQ = input.get("rpa", true);
  const auto str_om = input.get<std::string>("omega", "_");
  const bool eachFreqQ = str_om == "each" || str_om == "Each";
  const auto omega = eachFreqQ ? 0.0 : input.get("omega", 0.0);
  const auto print_both = input.get("printBoth", false);
  const auto only_diagonal = input.get("onlyDiagonal", false);

  // min/max n (for core/excited basis)
  const auto n_minmax = input.get_list("n_minmax", std::vector{1, 999});
  const auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  const auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";

  // For freq-dependent operators:
  if (h->freqDependantQ && !eachFreqQ)
    h->updateFrequency(omega);

  // do RPA:
  std::unique_ptr<HF::ExternalField> dV{nullptr};
  if (rpaQ) {
    dV = std::make_unique<HF::ExternalField>(h.get(), wf.getHF());
    if (!eachFreqQ)
      dV->solve_TDHFcore(omega);
  }

  std::cout << "\nStructure radiation and normalisation of states:\n";
  std::cout << "h=" << h->name() << " (reduced ME)\n";
  if (rpaQ) {
    std::cout << "RPA (TDHF) at";
    if (eachFreqQ)
      std::cout << " each frequency\n";
    else
      std::cout << "w=" << omega << "\n";
  }

  // Lambda to format the output
  const auto printer = [rpaQ](auto str, auto t, auto dv) {
    if (rpaQ)
      printf(" %6s: %12.5e  %12.5e\n", str, t, dv);
    else
      printf(" %6s: %12.5e\n", str, t);
  };

  if (wf.core.empty() || wf.valence.empty() || wf.basis.empty())
    return;

  // Find core/valence energy: allows distingush core/valence states
  const auto en_core = wf.en_coreval();

  // ----------- ** Actual Calculations ** -----------

  // Construct SR object:
  MBPT::StructureRad sr(wf.basis, en_core, {n_min, n_max});

  // Loop through all valence states, calc SR+NS
  for (const auto &v : wf.valence) {
    for (const auto &w : wf.valence) {
      if (h->isZero(w.k, v.k))
        continue;

      if (only_diagonal && w != v)
        continue;
      if (!print_both && v > w)
        continue;

      // nb: need splines to compute Struc Rad (use splines for legs)
      const auto ws = std::find(cbegin(wf.basis), cend(wf.basis), w);
      const auto vs = std::find(cbegin(wf.basis), cend(wf.basis), v);

      if (ws == cend(wf.basis) || vs == cend(wf.basis)) {
        std::cout << "Don't have requested spline for: " << w.symbol() << "-"
                  << v.symbol() << "\n";
        continue;
      }

      IO::ChronoTimer timer("time");

      std::cout << "\n" << h->rme_symbol(w, v) << ":\n";

      const auto ww = std::abs(ws->en - vs->en);
      if (eachFreqQ && h->freqDependantQ) {
        h->updateFrequency(ww);
      }
      if (eachFreqQ && rpaQ) {
        if (dV->get_eps() > 1.0e-3)
          dV->clear_dPsi();
        dV->solve_TDHFcore(ww);
      }

      // Zeroth-order MEs:
      const auto twvs = h->reducedME(*ws, *vs);
      const auto twv = h->reducedME(w, v);
      const auto dvs = dV ? twvs + dV->dV(*ws, *vs) : 0.0;
      const auto dv = dV ? twv + dV->dV(w, v) : 0.0;
      printer("t(spl)", twvs, dvs);
      printer("t(val)", twv, dv);

      // "Top" + "Bottom" SR terms:
      const auto [tb, tb_dv] = sr.srTB(h.get(), *ws, *vs, omega, dV.get());
      printer("SR(TB)", tb, tb_dv);
      // "Centre" SR term:
      const auto [c, c_dv] = sr.srC(h.get(), *ws, *vs, dV.get());
      printer("SR(C)", c, c_dv);

      std::cout << "========\n";
      printer("SR", tb + c, tb_dv + c_dv);

      // "Normalisation"
      const auto [n, n_dv] = sr.norm(h.get(), *ws, *vs, dV.get());
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
