#include "Modules/isotopeShift.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

#include <gsl/gsl_fit.h>

namespace Module {

void fieldShift(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Calculates field shift: F = d(E)/d(<r^2>)"},
       {"print", "Print each step? [true]"},
       {"min_pc", "Minimum percentage shift in r [1.0e-3]"},
       {"max_pc", "Maximum percentage shift in r [1.0e-1]"},
       {"num_steps", "Number of steps for derivative (for each sign)? [3]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto print = input.get("print", true);

  const auto min_pc = input.get("min_pc", 1.0e-3);
  const auto max_pc = input.get("max_pc", 1.0e-1);
  const auto num_steps = input.get<unsigned long>("num_steps", 3);

  Wavefunction wfB(wf.grid_sptr(), wf.nucleus(), wf.alpha() / PhysConst::alpha);

  std::cout << "Calculating field shift corrections for \n"
            << wf.atom() << ", " << wf.nucleus() << "\n"
            << "By fitting de = F<dr^2> for small delta r\n";

  wfB.copySigma(wf.Sigma());
  const auto core_string = wf.coreConfiguration();
  const auto val_string = DiracSpinor::state_config(wf.valence());
  const auto r0 = wf.get_rrms();

  std::vector<std::vector<std::pair<double, double>>> data(wf.valence().size());
  const auto delta_grid = Grid(r0 * min_pc / 100.0, r0 * max_pc / 100.0,
                               num_steps, GridType::logarithmic);

  if (print) {
    std::cout << "\n   r_rms (fm)    del(r)     del(r^2)     dE (GHz)   F "
                 "(GHz/fm^2)\n";
  } else {
    std::cout << "\nRunning...\n";
  }
  for (const auto pm : {-1, 1}) {
    for (const auto del : delta_grid.r()) {
      const auto rB = r0 + pm * del;
      const auto dr2 = rB * rB - r0 * r0;

      auto nuc_b = wf.nucleus();
      nuc_b.set_rrms(rB);

      wfB.update_Vnuc(Nuclear::formPotential(nuc_b, wf.grid().r()));

      wfB.solve_core("HartreeFock", 0.0, core_string, 0.0, false);
      wfB.solve_valence(val_string, false);
      wfB.hartreeFockBrueckner(false);

      for (auto i = 0ul; i < wfB.valence().size(); ++i) {
        const auto &Fv = wfB.valence()[i];
        const auto &Fv0 = *wf.getState(Fv.n(), Fv.kappa());
        const auto dE = (Fv.en() - Fv0.en()) * PhysConst::Hartree_GHz;
        const auto tF = dE / dr2;
        if (print)
          printf("%4s  %7.5f  %+7.5f  %11.4e  %11.4e  %10.3e\n",
                 Fv.shortSymbol().c_str(), rB, rB - r0, dr2, dE, tF);
        auto &data_v = data[i];
        data_v.emplace_back(dr2, dE);
      }
      if (print)
        std::cout << "\n";
    }
  }

  std::cout << "\n";

  auto sorter = [](auto p1, auto p2) { return p1.first < p2.first; };

  for (auto i = 0ul; i < wfB.valence().size(); ++i) {
    const auto &Fv = wfB.valence()[i];
    auto &data_v = data[i];
    std::sort(begin(data_v), end(data_v), sorter);

    [[maybe_unused]] double c0, c1, cov00, cov01, cov11, sumsq;

    // Fit, without c0
    gsl_fit_mul(&data_v[0].first, 2, &data_v[0].second, 2, data_v.size(), &c1,
                &cov11, &sumsq);
    std::cout << Fv.symbol() << " "
              << "F = " << c1 << " GHz/fm^2, sd = " << std::sqrt(sumsq) << "\n";
  }
}

} // namespace Module
