#include "Modules/isotopeShift.hpp"
#include "DiracOperator/DiracOperator.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

#include <gsl/gsl_fit.h>

namespace Module {

void fieldShift(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check(
      {{"", "Calculates field shift: F = d(E)/d(<r^2>) by direct calculation"},
       {"core_relaxation", "Include Core relaxation (equiv to RPA)? [true]"},
       {"print", "Print each step to screen? [true]"},
       {"write", "Write dE(r^2) to file? [false]"},
       {"minmax_delta", "Minimum relative shift in r [1.0e-5]"},
       {"num_steps", "Number of steps for fit (for each sign)? [5]"},
       {"grid", "Logarithmic or linear grid for dr2 [logarithmic]"}});
  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto core_relax = input.get("core_relax*", true);
  const auto print = input.get("print", true);
  const auto write = input.get("write", false);

  const auto [min_d, max_d] =
      input.get("minmax_delta", std::array{1.0e-5, 1.0e-2});

  const auto num_steps = input.get<unsigned long>("num_steps", 5);

  const auto grid_type = input.get("grid", std::string{"logarithmic"});

  Wavefunction wfB(wf.grid_sptr(), wf.nucleus(), wf.alpha() / PhysConst::alpha);

  std::cout << "Calculating field shift corrections for \n"
            << wf.atom() << ", " << wf.nucleus() << "\n"
            << "By fitting de = F<dr^2> for small delta r\n";

  wfB.copySigma(wf.Sigma());
  const auto core_string = wf.coreConfiguration();
  const auto val_string = DiracSpinor::state_config(wf.valence());
  const auto r0 = wf.get_rrms();

  std::cout << "\n";
  wfB.solve_core("HartreeFock", 0.0, core_string, 0.0, true);
  if (!core_relax) {
    std::cout << "Not including Core relaxtion\n";
  } else {
    std::cout << "Including Core relaxtion\n";
  }

  std::ofstream of;
  if (write) {
    of.open(wf.identity() + "_FS.txt");
    of << "dr dr2 ";
    for (auto &v : wf.valence()) {
      of << " " << v;
    }
    of << "\n";
  }

  const auto gtype = qip::ci_wc_compare(grid_type, "log*") ?
                         GridType::logarithmic :
                         GridType::linear;

  std::vector<std::vector<std::pair<double, double>>> data(wf.valence().size());
  const auto delta_grid = Grid(r0 * min_d, r0 * max_d, num_steps, gtype);

  if (print) {
    std::cout << "\n   r_rms (fm)    del(r)     del(r^2)     dE (GHz)   F "
                 "(GHz/fm^2)\n";
  } else {
    std::cout << "\nRunning...\n";
  }
  for (const auto pm : {-1, 1}) {
    const auto dr_list =
        pm == -1 ? qip::reverse(delta_grid.r()) : delta_grid.r();
    for (const auto del : dr_list) {
      const auto rB = r0 + pm * del;
      const auto dr2 = rB * rB - r0 * r0;

      auto nuc_b = wf.nucleus();
      nuc_b.set_rrms(rB);

      wfB.update_Vnuc(Nuclear::formPotential(nuc_b, wf.grid().r()));

      if (write) {
        of << rB - r0 << " " << " " << dr2;
      }

      if (core_relax)
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

        if (write) {
          of << " " << dE;
        }
      }
      if (write) {
        of << "\n";
      }

      if (print)
        std::cout << "\n";
    }
  }

  std::cout << "\n";

  // Fit straight line to data
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
