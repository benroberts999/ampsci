#pragma once
#include "DiracOperator/DiracOperator.hpp"
#include "MixedStates.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Check.hpp"
#include "qip/String.hpp"
#include <iomanip>
#include <string>
//
#include "ExternalField/TDHF.hpp"

namespace UnitTest {

//******************************************************************************
namespace helper {

struct Case {
  double eps;
  std::string name = "";
  double w = 0.0;
};

// Compares: <B||dA> to <B||h||A> / (e_A - e_B + w)
// Where dA is solution to: (H - e_A - w)|dA> = -(h - de_A)|A>
// Calculates above comparison for each pair of valence states (A,B), for which
// the matrix element (and denominator) is non-zero. Also loops over different w
// (frequency) values.
// Returns the best and worst cases
inline std::pair<Case, Case> MS_loops(const Wavefunction &wf,
                                      const DiracOperator::TensorOperator *h);

} // namespace helper

//******************************************************************************
//******************************************************************************
//! Unit tests Mixed States (TDHF method, solving non-local DE)
bool MixedStates(std::ostream &obuff) {
  bool passQ = true;

  // Create wavefunction object
  Wavefunction wf({6000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.solve_valence("7sp5d4f");

  // Define E1, E2, pnc, M1, and hyperfine operators:
  const auto hE1 = DiracOperator::E1(*wf.rgrid);
  const auto hE2 = DiracOperator::Ek(*wf.rgrid, 2);
  const auto c = Nuclear::c_hdr_formula_rrms_t(wf.get_rrms());
  const auto hpnc = DiracOperator::PNCnsi(c, Nuclear::default_t, *wf.rgrid);
  // M1 little problematic
  // const auto hM1 = DiracOperator::M1(*wf.rgrid, wf.alpha, 0.0);
  // Use "spherical ball" model for hyperfine (Works best.)
  // Fails for some d states with pointlike (?)
  const auto hhfs = DiracOperator::HyperfineA(
      1.0, 1.0, std::sqrt(5.0 / 3) * wf.get_rrms(), *wf.rgrid,
      DiracOperator::Hyperfine::sphericalBall_F());

  const std::vector<const DiracOperator::TensorOperator *> hs{&hE1, &hE2, &hpnc,
                                                              &hhfs /*, &hM1*/};

  // For each operator
  for (const auto h : hs) {

    const auto [best, worst] = helper::MS_loops(wf, h);

    const auto wbest = qip::fstring("%.2g", best.w);
    const auto wworst = qip::fstring("%.2g", worst.w);

    passQ &= qip::check_value(&obuff,
                              "MS:" + h->name() + " worst (" + worst.name +
                                  " w=" + wworst + ")",
                              worst.eps, 0.0, 4e-5);
    passQ &= qip::check_value(
        &obuff, "MS:" + h->name() + " best (" + best.name + " w=" + wbest + ")",
        best.eps, 0.0, 6e-7);
  }

  // Since we have trouble with TDHF and HFS, do it again more thoroughly here.
  // Note: This makes it seem the problem is in solve_core, NOT in solve dPsi
  // This should be easier to fix!
  // NOTE: This loop is very much the same as the other one... could combine

  // Hyperfine
  std::cout << "\nTest hyperfine (again, but more)\n";
  auto &h = hhfs;

  for (int max_its : {0, 1, 99}) {
    double worst_eps = 0.0;
    std::string worst_set{};
    double best_eps = 1.0;
    std::string best_set{};
    ExternalField::TDHF dv(&h, wf.getHF());
    dv.solve_core(0.0, max_its);

    // to test the .get() X,Y's
    ExternalField::TDHF dPsi(&h, wf.getHF());
    dPsi.solve_core(0.0, 1); // 1 it; no dV, but solve for dPsi

    int count = 0;
    for (const auto &Fv : wf.valence) {
      for (const auto &Fm : wf.valence) {
        if (Fm == Fv || h.isZero(Fm.k, Fv.k))
          continue;

        const auto Xb =
            dv.solve_dPsi(Fv, 0.0, ExternalField::dPsiType::X, Fm.k);
        const auto Yb = dv.solve_dPsi(Fv, 0.0, ExternalField::dPsiType::Y, Fm.k,
                                      nullptr, ExternalField::StateType::bra);
        const auto h_mv = h.reducedME(Fm, Fv) + dv.dV(Fm, Fv);
        const auto lhs = Fm * Xb;
        const auto rhs = h_mv / (Fv.en() - Fm.en());
        const auto eps = (lhs - rhs) / (lhs + rhs);

        const auto h_vm = h.reducedME(Fv, Fm) + dv.dV(Fv, Fm);
        const auto lhsY = Yb * Fm;
        const auto rhsY = h_vm / (Fv.en() - Fm.en());
        const auto epsY = (lhsY - rhsY) / (lhsY + rhsY);

        if (count % 5 == 0 || count == 0) {
          // only print a subset (too many)
          std::cout << "<" << Fv.shortSymbol() << "|h|" << Fm.shortSymbol()
                    << "> , <" << Fv.shortSymbol() << "|X_" << Xb.shortSymbol()
                    << "> : "; // << rhs << " " << lhs << "\n";
          printf("%10.7g, %10.7g  %6.0e\n", rhs, lhs, eps);
          std::cout << "<" << Fm.shortSymbol() << "|h|" << Fv.shortSymbol()
                    << "> , <Y_" << Yb.shortSymbol() << "|" << Fv.shortSymbol()
                    << "> : "; // << rhsY << " " << lhsY << "\n";

          printf("%10.7g, %10.7g  %6.0e", rhsY, lhsY, epsY);
          if (std::abs(epsY + eps) > 1.0e-3)
            std::cout << "  ***";
          std::cout << "\n";
        }
        ++count;

        if (std::abs(eps) > std::abs(worst_eps)) {
          worst_eps = eps;
          worst_set = "<" + Fv.shortSymbol() + "|h|" + Fm.shortSymbol() +
                      ">/<" + Fv.shortSymbol() + "|X_" + Xb.shortSymbol() + ">";
        }
        if (std::abs(epsY) > std::abs(worst_eps)) {
          worst_eps = epsY;
          worst_set = "<" + Fm.shortSymbol() + "|h|" + Fv.shortSymbol() +
                      ">/<Y_" + Xb.shortSymbol() + "|" + Fv.shortSymbol() + ">";
        }
        if (std::abs(eps) < std::abs(best_eps)) {
          best_eps = eps;
          best_set = "<" + Fv.shortSymbol() + "|h|" + Fm.shortSymbol() + ">/<" +
                     Fv.shortSymbol() + "|X_" + Xb.shortSymbol() + ">";
        }
        if (std::abs(epsY) < std::abs(best_eps)) {
          best_eps = epsY;
          best_set = "<" + Fm.shortSymbol() + "|h|" + Fv.shortSymbol() +
                     ">/<Y_" + Xb.shortSymbol() + "|" + Fv.shortSymbol() + ">";
        }
      }
    }
    std::cout << worst_set << " " << worst_eps << "\n";
    // the "best" are all ~1.e-y
    passQ &= qip::check_value(&obuff, "MS: hfs(B) " + best_set, best_eps, 0.0,
                              1.0e-7);
    passQ &= qip::check_value(&obuff, "MS: hfs(W) " + worst_set, worst_eps, 0.0,
                              1.0e-4);
  }

  return passQ;
}

} // namespace UnitTest

//******************************************************************************
inline std::pair<UnitTest::helper::Case, UnitTest::helper::Case>
UnitTest::helper::MS_loops(const Wavefunction &wf,
                           const DiracOperator::TensorOperator *h) {
  // compare the best and worst! (of each operator)
  Case best{999.0};
  Case worst{0.0};

  const auto omega_mults = std::vector{0.0, 0.25};

#pragma omp parallel for
  for (auto i = 0ul; i < wf.valence.size(); ++i) {
    const auto &Fv = wf.valence[i];
    const auto vl = wf.get_Vlocal(Fv.l());
    for (const auto &Fm : wf.valence) {

      // Only do for cases that make sense:
      if (Fm == Fv)
        continue; // gives 0/0
      if (h->isZero(Fv.k, Fm.k))
        continue;
      if (Fm.k != Fv.k && std::abs(Fm.n - Fv.n) != 0)
        continue;

      const auto h_mv = h->reducedME(Fm, Fv);

      // Don't do in cases where ME is extremely small. NOTE: This requires
      // a good choice of units, since ME may be small due to small
      // coeficient, which should not count as 'small'
      if (std::abs(h_mv) < 1.0e-8) {
        continue;
      }

      // Form 'rhs': (h - de_A)|A>
      auto hFv = h->reduced_rhs(Fm.k, Fv);
      if (Fm.k == Fv.k) {
        auto de = h->reducedME(Fv, Fv);
        hFv -= de * Fv;
      }

      // loop over few frequencies:
      for (const auto &w_mult : omega_mults) {
        const auto w = std::abs(Fv.en() * w_mult);

        const auto dFv = ExternalField::solveMixedState(Fm.k, Fv, w, vl,
                                                        wf.alpha, wf.core, hFv);

        const auto lhs = Fm * dFv;
        const auto rhs = h_mv / (Fv.en() - Fm.en() + w);
        const auto eps = (lhs - rhs) / (lhs + rhs);

// find the best and worst case:
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
