#pragma once
#include "DiracOperator/Operators.hpp"
#include "MixedStates.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include <iomanip>
#include <string>

namespace UnitTest {

namespace helper {
struct Case {
  double eps;
  std::string name = "";
  double w = 0.0;
};
} // namespace helper

//******************************************************************************
namespace helper {
std::pair<Case, Case> MS_loops(const Wavefunction &wf,
                               const DiracOperator::TensorOperator *h) {
  // compare the best and worst! (of each operator)
  Case best{999.0};
  Case worst{0.0};

  const auto omega_mults = std::vector{0.0, 0.25};

// for (const auto &Fv : wf.valence) {
#pragma omp parallel for
  for (auto i = 0ul; i < wf.valence.size(); ++i) {
    const auto &Fv = wf.valence[i];
    const auto vl = wf.get_Vlocal(Fv.l());
    for (const auto &Fm : wf.valence) {
      if (Fm == Fv)
        continue; // gives 0/0
      if (h->isZero(Fv.k, Fm.k))
        continue;

      if (Fm.k != Fv.k && std::abs(Fm.n - Fv.n) != 0)
        continue;

      // std::cout << h->name() << " " << Fv.symbol() << " " << Fm.symbol()
      //           << "\n";

      const auto h_mv = h->reducedME(Fm, Fv);
      // Don't check in cases where ME is extremely small. NOTE: This requires
      // a good choice of units, since ME may be small due to small
      // coeficient, which should not count as 'small'
      if (std::abs(h_mv) < 1.0e-8) {
        // std::cout << h_mv << "\n";
        continue;
      }

      auto hFv = h->reduced_rhs(Fm.k, Fv);
      if (Fm.k == Fv.k) {
        auto de = h->reducedME(Fv, Fv);
        hFv -= de * Fv;
        // std::cout << "de=" << de << "\n";
      }

      for (const auto &w_mult : omega_mults) {
        const auto w = std::abs(Fv.en * w_mult);

        const auto dFv =
            HF::solveMixedState(Fm.k, Fv, w, vl, wf.alpha, wf.core, hFv);

        const auto lhs = Fm * dFv;
        const auto rhs = h_mv / (Fv.en - Fm.en + w);
        const auto eps = (lhs - rhs) / lhs;

// find worst-case:
#pragma omp critical(find_max)
        {
          if (std::abs(eps) > std::abs(worst.eps)) {
            worst.eps = eps;
            worst.name = Fm.shortSymbol() + "|" + Fv.shortSymbol();
            worst.w = w;
          }
          if (std::abs(eps) < std::abs(best.eps)) {
            best.eps = eps;
            best.name = Fm.shortSymbol() + "|" + Fv.shortSymbol();
            best.w = w;
          }
        }
      }
    }
  }
  return {best, worst};
}
} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests Mixed States (TDHF method, solving non-local DE)
bool MixedStates(std::ostream &obuff) {
  bool passQ = true;

  // Create wavefunction object
  Wavefunction wf({10000, 1.0e-7, 200.0, 0.33 * 200.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.hartreeFockCore("HartreeFock", 0.0, "[Xe]");
  wf.hartreeFockValence("7sp5d4f");

  auto hE1 = DiracOperator::E1(*wf.rgrid);
  auto hE2 = DiracOperator::Ek(*wf.rgrid, 2);
  auto c = Nuclear::c_hdr_formula_rrms_t(wf.get_rrms());
  auto hpnc = DiracOperator::PNCnsi(c, Nuclear::default_t, *wf.rgrid);
  auto hM1 = DiracOperator::M1(*wf.rgrid, wf.alpha, 0.0);
  // spherical ball"
  auto hhfs = DiracOperator::Hyperfine(
      1.0, 1.0, std::sqrt(5.0 / 3) * wf.get_rrms(), *wf.rgrid,
      DiracOperator::Hyperfine::sphericalBall_F());

  std::vector<DiracOperator::TensorOperator *> hs{&hE1, &hE2, &hpnc, &hhfs,
                                                  &hM1};

  for (const auto h : hs) {
    // nb: could do each sepperately ..

    const auto [best, worst] = helper::MS_loops(wf, h);

    std::stringstream wbest, wworst;
    wbest << std::setprecision(2) << best.w;
    wworst << std::setprecision(2) << worst.w;

    passQ &= qip::check_value(&obuff,
                              "MS:" + h->name() + " worst (" + worst.name +
                                  " w=" + wworst.str() + ")",
                              worst.eps, 0.0, 4e-5);
    passQ &= qip::check_value(&obuff,
                              "MS:" + h->name() + " best (" + best.name +
                                  " w=" + wbest.str() + ")",
                              best.eps, 0.0, 6e-7);
  }

  return passQ;
}

} // namespace UnitTest
