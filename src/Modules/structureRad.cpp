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

void structureRad(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"operator", "options", "rpa", "printBoth", "onlyDiagonal",
                    "omega", "n_minmax"});

  // Find core/valence energy: allows distingush core/valence states
  const auto ec_max =
      std::max_element(cbegin(wf.core), cend(wf.core), DiracSpinor::comp_en)
          ->en;
  const auto ev_min = std::min_element(cbegin(wf.valence), cend(wf.valence),
                                       DiracSpinor::comp_en)
                          ->en;
  const auto en_core = 0.5 * (ev_min + ec_max);

  const auto oper = input.get<std::string>("operator", "E1");
  const auto h_options =
      IO::UserInputBlock(oper, input.get<std::string>("options", ""));
  const auto h = generateOperator(oper, h_options, wf, true);

  const auto rpaQ = input.get("rpa", true);
  const auto omega = input.get("omega", 0.0);

  const auto print_both = input.get("printBoth", false);
  const auto only_diagonal = input.get("onlyDiagonal", false);

  const auto n_minmax = input.get_list("n_minmax", std::vector{1, 999});
  auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";

  std::unique_ptr<HF::ExternalField> dV{nullptr};
  if (rpaQ) {
    dV = std::make_unique<HF::ExternalField>(h.get(), wf.getHF());
    dV->solve_TDHFcore(omega);
  }

  auto printer = [rpaQ](auto str, auto t, auto dv) {
    if (rpaQ)
      printf(" %6s: %12.5e  %12.5e\n", str, t, dv);
    else
      printf(" %6s: %12.5e\n", str, t);
  };

  std::cout << "\nStructure radiation and normalisation of states:\n";
  std::cout << "h=" << h->name() << " (reduced ME)\n";

  MBPT::StructureRad sr(wf.basis, en_core, {n_min, n_max});

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

      IO::ChronoTimer timer();

      std::cout << "\n" << h->rme_symbol(w, v) << ":\n";

      const auto twvs = h->reducedME(*ws, *vs);
      const auto twv = h->reducedME(w, v);
      const auto dvs = dV ? twvs + dV->dV(*ws, *vs) : 0.0;
      const auto dv = dV ? twv + dV->dV(w, v) : 0.0;

      printer("t(spl)", twvs, dvs);
      printer("t(val)", twv, dv);

      const auto [tb, tb_dv] = sr.srTB(h.get(), *ws, *vs, omega, dV.get());
      printer("SR(TB)", tb, tb_dv);

      const auto [c, c_dv] = sr.srC(h.get(), *ws, *vs, dV.get());
      printer("SR(C)", c, c_dv);

      std::cout << "========\n";
      printer("SR", tb + c, tb_dv + c_dv);
      const auto [n, n_dv] = sr.norm(h.get(), *ws, *vs, dV.get());
      printer("Norm", n, n_dv);
    }
  }
  std::cout << "\n";

  return;
  //****************************************************************************
}

} // namespace Module
