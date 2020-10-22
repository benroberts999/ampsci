#include "MBPT/StructureRad.hpp"
#include "DiracOperator/Operators.hpp"
#include "HF/ExternalField.hpp"
#include "IO/ChronoTimer.hpp"
#include "IO/UserInput.hpp"
#include "Modules/matrixElements.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace Module {

void structureRad(const IO::UserInputBlock &input, const Wavefunction &wf) {

  input.checkBlock({"operator", "options", "printBoth", "onlyDiagonal", "omega",
                    "n_minmax"});

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

  const auto print_both = input.get("printBoth", false);
  const auto only_diagonal = input.get("onlyDiagonal", false);

  const auto n_minmax = input.get_list("n_minmax", std::vector{1, 999});
  auto n_min = n_minmax.size() > 0 ? n_minmax[0] : 1;
  auto n_max = n_minmax.size() > 1 ? n_minmax[1] : 999;
  if (n_min > 1)
    std::cout << "Including from n = " << n_min << "\n";
  if (n_max < 999)
    std::cout << "Including to n = " << n_max << "\n";

  std::cout
      << "\nStructure radiation and normalisation of states (reduced ME):\n";
  MBPT::StructureRad sr(wf.basis, en_core, {1, 60});
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

      IO::ChronoTimer timer("t:");

      std::cout << w.symbol() << "-" << v.symbol() << ":\n"
                << h->reducedME(*ws, *vs) << " " << h->reducedME(w, v)
                << std::endl;

      const auto n = sr.norm(h.get(), *ws, *vs);
      std::cout << " Norm : " << n << "\n" << std::flush;

      const auto [t, b] = sr.srTB(h.get(), *ws, *vs, 0.0);
      std::cout << " SR(T): " << t << "\n" << std::flush;
      std::cout << " SR(B): " << b << "\n" << std::flush;

      const auto c = sr.srC(h.get(), *ws, *vs);
      std::cout << " SR(C): " << c << "\n" << std::flush;
      std::cout << " SR   : " << t + b + c << "\n" << std::flush;
      std::cout << "\n";
    }
  }

  return;
  //****************************************************************************
}

} // namespace Module
