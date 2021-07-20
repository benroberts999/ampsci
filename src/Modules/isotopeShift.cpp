#include "Modules/isotopeShift.hpp"
#include "DiracOperator/Operators.hpp" //For E1 operator
#include "IO/InputBlock.hpp"
#include "Physics/PhysConst_constants.hpp" // For GHz unit conversion
#include "Wavefunction/Wavefunction.hpp"

#include <gsl/gsl_fit.h>

namespace Module {

void fieldShift(const IO::InputBlock &, const Wavefunction &wfA) {

  Wavefunction wfB(wfA.rgrid->params(), wfA.get_nuclearParameters(),
                   wfA.alpha / PhysConst::alpha);

  std::cout << "\n";
  IO::print_line();
  std::cout << "Calculating field shift corrections for \n"
            << wfA.atom() << ", " << wfA.nuclearParams() << "\nvs.\n"
            << wfB.atom() << ", " << wfB.nuclearParams() << "\n\n";

  wfB.copySigma(wfA.getSigma());
  wfB.solve_core("HartreeFock", 0.0, wfA.coreConfiguration());
  wfB.solve_valence(DiracSpinor::state_config(wfA.valence));
  wfB.hartreeFockBrueckner();
  wfB.printValence();

  const auto vA = wfA.vnuc;
  const auto vB = wfB.vnuc;

  const auto r0B = wfB.get_rrms();

  // const auto enB0 = wfB.valence[0].en();
  std::vector<std::vector<std::pair<double, double>>> data(wfB.valence.size());
  // data.emplace_back(0.0, 0.0);

  const auto min_pc = 0.001;
  const auto max_pc = 1.0;
  const auto num_steps = 10;

  const auto delta_grid = Grid(r0B * min_pc / 100.0, r0B * max_pc / 100.0,
                               num_steps, GridType::logarithmic);
  const auto t = 2.3;

  std::cout << "\n   r_rms (fm),   del(r),    del(r^2),    dE (GHz)\n";
  for (const auto pm : {-1, 1}) {
    for (const auto del : delta_grid.r()) {
      const auto rB = r0B + pm * del;
      const auto dr2 = r0B * r0B - rB * rB;

      wfB.vnuc = Nuclear::fermiNuclearPotential(
          wfB.Znuc(), t, Nuclear::c_hdr_formula_rrms_t(rB, t), wfA.rgrid->r());

      wfB.solve_core("", 0.0, "", 0.0, false);
      wfB.solve_valence("", false);
      wfB.hartreeFockBrueckner(false);

      for (auto i = 0ul; i < wfB.valence.size(); ++i) {
        const auto &Fv = wfB.valence[i];
        const auto &Fv0 = *wfA.getState(Fv.n, Fv.k);
        const auto dE = -(Fv.en() - Fv0.en()) * PhysConst::Hartree_GHz;

        printf("%4s, %7.5f, %+7.5f, %11.4e, %11.4e\n", Fv.shortSymbol().c_str(),
               rB, rB - r0B, dr2, dE);
        auto &data_v = data[i];
        data_v.emplace_back(dr2, dE);
      }
    }
  }

  std::cout << "\n";
  // for (const auto &[dr2, de] : data) {
  //   std::cout << dr2 << " " << de << "\n";
  // }
  // std::cout << "\n";

  auto sorter = [](auto p1, auto p2) { return p1.first < p2.first; };

  for (auto i = 0ul; i < wfB.valence.size(); ++i) {
    const auto &Fv = wfB.valence[i];
    auto &data_v = data[i];
    std::sort(begin(data_v), end(data_v), sorter);

    [[maybe_unused]] double c0, c1, cov00, cov01, cov11, sumsq;

    /*
        // nb: abusing fact that vector and pair store doubles in contiguous
       memory gsl_fit_linear(&data_v[0].first, 2, &data_v[0].second, 2,
       data_v.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq); printf("Eb -
       Eb0 = %.5e d<r2> + %2e\n", c1, c0); std::cout << Fv.symbol() << " "
                  << "F = " << c1 << " GHz/fm^2, sd = " << std::sqrt(sumsq) <<
       "\n";
    */

    // Fit, without c0
    gsl_fit_mul(&data_v[0].first, 2, &data_v[0].second, 2, data_v.size(), &c1,
                &cov11, &sumsq);
    std::cout << Fv.symbol() << " "
              << "F = " << c1 << " GHz/fm^2, sd = " << std::sqrt(sumsq) << "\n";
  }
}

} // namespace Module
