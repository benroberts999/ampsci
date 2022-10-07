#include "MixedStates.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/String.hpp"
#include <iomanip>
#include <string>

//==============================================================================
TEST_CASE("External Field: Mixed-states",
          "[ExternalField][MixedStates][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: Mixed-states (unit)\n";

  // Create wavefunction object
  Wavefunction wf({1200, 1.0e-5, 80.0, 20.0, "loglinear", -1.0},
                  {"Rb", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Kr]");

  std::cout << "\nTest Mixed-States (TDHF) equation (no dV): \n"
            << "(H-e)Xc = -(h-de)Fc\n"
            << "Compare <a|Xc>/(ea-ec) to <a||h||c>\n";

  std::stringstream summary;

  // test E1 (odd+vector+real), pnc (odd+scalar+imag), hfs (even+vector+real)
  for (std::string_view oper : {"E1", "hfs", "pnc"}) {
    std::cout << "\n" << oper << "\n";
    auto h = DiracOperator::generate(oper, {}, wf);

    std::string worst{""};
    double weps{0.0};
    std::string best{""};
    double beps{10.0};

    for (auto omega : {0.0, 0.05, -0.1}) {
      if (oper != "E1" && omega != 0.0)
        continue;
      for (auto &Fc : wf.core()) {
        for (const auto &Fm : wf.core()) {
          if (Fc <= Fm || h->isZero(Fm, Fc) || std::abs(Fm.n() - Fc.n()) > 1)
            continue;
          if (omega == 0.0 && Fc == Fm)
            continue;

          auto hmc1 = h->reducedME(Fm, Fc);
          if (std::abs(hmc1) < 1.0e-5)
            continue;

          const auto hFc = h->reduced_rhs(Fm.kappa(), Fc);
          auto imag = h->imaginaryQ();
          auto de = Fc.kappa() == Fm.kappa() && !imag ? Fc * hFc : 0.0;
          const auto &dF = ExternalField::solveMixedState(
              Fc, omega, hFc - de * Fc, wf.vHF());

          auto hmc2 = (Fm * dF) * (Fc.en() - Fm.en() + omega);
          auto eps = std::abs((hmc1 - hmc2) / hmc1);
          auto del = std::abs((hmc1 - hmc2));
          eps = std::min(eps, del);
          std::cout << h->name() << " " << Fm << "|" << Fc << " : ";
          printf("%+.4e [%+.4e] %.0e w=%.2f\n", hmc2, hmc1, eps, omega);
          if (eps > weps) {
            weps = eps;
            worst = Fm.shortSymbol() + "|" + Fc.shortSymbol();
          }
          if (eps < beps) {
            beps = eps;
            best = Fm.shortSymbol() + "|" + Fc.shortSymbol();
          }
        }
      }
    }

    REQUIRE(weps < 1.0e-3);
    REQUIRE(beps < 1.0e-5);
    summary << h->name() << " worst: " << worst << " "
            << qip::fstring(".1e", weps) << ", best: " << best << " "
            << qip::fstring(".1e", beps) << "\n";
  }
  std::cout << "\nMixed-States summary:\n" << summary.str() << "\n";
}

//==============================================================================
TEST_CASE("External Field: Mixed-states2",
          "[ExternalField][MixedStates][integration][!mayfail]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: Mixed-states2\n";

  // This test illucidates the issue with TDHF for hyperfine:

  // Create wavefunction object
  Wavefunction wf({2000, 1.0e-5, 80.0, 20.0, "loglinear", -1.0},
                  {"Rb", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Kr]");
  wf.solve_valence("5sp");

  SplineBasis::Parameters bspl_param;
  {
    bspl_param.states = "80spdfg";
    bspl_param.n = 90;
    bspl_param.k = 9;
    bspl_param.r0 = 1.0e-6;
    bspl_param.reps = 0.0;
    bspl_param.rmax = 60.0;
    bspl_param.positronQ = false;
  }

  wf.formBasis(bspl_param);

  for (std::string_view oper : {"E1", "hfs"}) {
    std::cout << "\n" << oper << "\n";
    auto h = DiracOperator::generate(oper, {}, wf);

    ExternalField::TDHF tdhf(h.get(), wf.vHF());
    ExternalField::TDHFbasis basis(h.get(), wf.vHF(), wf.basis());

    const auto r = wf.grid().r();

    for (auto omega : {0.0, 0.05}) {

      tdhf.clear();
      basis.clear();
      tdhf.solve_core(omega, 1, false);
      basis.solve_core(omega, 1, false);

      if (oper != "E1" && omega != 0.0)
        continue;
      for (auto &Fc : wf.core()) {
        for (const auto &Fm : wf.core()) {
          if (Fc <= Fm || h->isZero(Fm, Fc) || std::abs(Fm.n() - Fc.n()) > 1)
            continue;
          if (omega == 0.0 && Fc == Fm)
            continue;

          auto hmc1 = h->reducedME(Fm, Fc);
          if (std::abs(hmc1) < 1.0e-5)
            continue;

          const auto hFc = h->reduced_rhs(Fm.kappa(), Fc);
          auto imag = h->imaginaryQ();
          auto de = Fc.kappa() == Fm.kappa() && !imag ? Fc * hFc : 0.0;
          const auto &dF = ExternalField::solveMixedState(
              Fc, omega, hFc - de * Fc, wf.vHF());

          auto dF2 =
              tdhf.get_dPsi_x(Fc, ExternalField::dPsiType::X, hFc.kappa());
          auto dF3 =
              basis.get_dPsi_x(Fc, ExternalField::dPsiType::X, hFc.kappa());

          auto hmc2 = (Fm * dF) * (Fc.en() - Fm.en() + omega);
          auto eps = std::abs((hmc1 - hmc2) / hmc1);
          auto del = std::abs((hmc1 - hmc2));
          eps = std::min(eps, del);
          std::cout << Fm << "|" << Fc << " ";
          printf("w=%.2f %+.3e %+.3e ", omega, hmc1, hmc2);

          auto hmc3 = (Fm * dF2) * (Fc.en() - Fm.en() + omega);
          auto hmc4 = (Fm * dF3) * (Fc.en() - Fm.en() + omega);

          printf("%+.3e %+.3e  %.1e", hmc3, hmc4, eps);
          if (eps > 1.0e-2)
            std::cout << "  ****";
          std::cout << "\n";

          auto n1 = dF2.norm();
          auto n2 = dF3.norm();
          auto d3 = (dF2 - dF3).norm() / dF2.norm();
          printf("               %+.3e %+.3e  %.1e", n1, n2, d3);
          if (d3 > 1.0e-2)
            std::cout << "  ******";
          std::cout << "\n";

          auto dv1 = tdhf.dV(Fc, Fm);
          auto dv2 = basis.dV(Fc, Fm);
          // auto dv3 = basis.dV1(Fc, Fm);
          auto eps3 = std::abs((dv1 - dv2) / dv1);
          printf("               %+.3e %+.3e  %.1e", dv1, dv2, eps3);
          if (eps3 > 1.0e-2)
            std::cout << "  ******";
          std::cout << "\n";
        }
      }

      double worst = 0.0;
      for (auto &Fc : wf.core()) {
        const auto &dF_t = tdhf.get_dPsis(Fc, ExternalField::dPsiType::X);
        const auto &dF_b = basis.get_dPsis(Fc, ExternalField::dPsiType::X);
        const auto &dY_t = tdhf.get_dPsis(Fc, ExternalField::dPsiType::Y);
        const auto &dY_b = basis.get_dPsis(Fc, ExternalField::dPsiType::Y);
        assert(dF_t.size() == dF_b.size());
        for (std::size_t ix = 0; ix < dF_t.size(); ix++) {
          const auto &dF_tx = dF_t[ix];
          const auto &dF_bx = dF_b[ix];
          const auto &dY_tx = dY_t[ix];
          const auto &dY_bx = dY_b[ix];
          assert(dF_tx.kappa() == dF_bx.kappa());
          std::cout << Fc << "/" << dF_tx << ": " << dF_tx.norm() << " "
                    << dF_bx.norm() << " | ";
          std::cout << dY_tx.norm() << " " << dY_bx.norm() << "\n";
          auto eps = std::abs(dF_tx.norm() / dF_bx.norm() - 1.0);
          if (eps > worst) {
            worst = eps;
          }
        }
      }
      REQUIRE(worst < 1.0e-2);
    }
  }
}

//==============================================================================
//==============================================================================
// //==============================================================================
// namespace helper {

// struct Case {
//   double eps;
//   std::string name = "";
//   double w = 0.0;
// };

// // Compares: <B||dA> to <B||h||A> / (e_A - e_B + w)
// // Where dA is solution to: (H - e_A - w)|dA> = -(h - de_A)|A>
// // Calculates above comparison for each pair of valence states (A,B), for which
// // the matrix element (and denominator) is non-zero. Also loops over different w
// // (frequency) values.
// // Returns the best and worst cases
// inline std::pair<Case, Case> MS_loops(const Wavefunction &wf,
//                                       const DiracOperator::TensorOperator *h);

// } // namespace helper

// //==============================================================================
// //==============================================================================
// //! Unit tests Mixed States (TDHF method, solving non-local DE)
// TEST_CASE("External Field: Mixed-states (old)",
//           "[ExternalField][MixedStates][integration]") {
//   std::cout << "\n----------------------------------------\n";
//   std::cout << "External Field: Mixed-states\n";

//   // Create wavefunction object
//   Wavefunction wf({6000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
//                   {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
//   wf.solve_core("HartreeFock", 0.0, "[Xe]");
//   wf.solve_valence("7sp5d4f");

//   // Define E1, E2, pnc, M1, and hyperfine operators:
//   const auto hE1 = DiracOperator::E1(wf.grid());
//   const auto hE2 = DiracOperator::Ek(wf.grid(), 2);
//   const auto c = Nuclear::c_hdr_formula_rrms_t(wf.get_rrms());
//   const auto hpnc = DiracOperator::PNCnsi(c, Nuclear::default_t, wf.grid());
//   // M1 little problematic
//   // const auto hM1 = DiracOperator::M1(wf.grid(), wf.alpha(), 0.0);
//   // Use "spherical ball" model for hyperfine (Works best.)
//   // Fails for some d states with pointlike (?)
//   const auto hhfs = DiracOperator::HyperfineA(
//       1.0, 1.0, std::sqrt(5.0 / 3) * wf.get_rrms(), wf.grid(),
//       DiracOperator::Hyperfine::sphericalBall_F());

//   const std::vector<const DiracOperator::TensorOperator *> hs{&hE1, &hE2, &hpnc,
//                                                               &hhfs /*, &hM1*/};

//   // For each operator
//   for (const auto h : hs) {

//     const auto [best, worst] = helper::MS_loops(wf, h);

//     const auto wbest = qip::fstring("%.2g", best.w);
//     const auto wworst = qip::fstring("%.2g", worst.w);

//     std::cout << "Mixed states: " + h->name() + " worst (" + worst.name +
//                      " w=" + wworst + ") "
//               << worst.eps << "\n";
//     std::cout << "Mixed states: " + h->name() + " best (" + best.name +
//                      " w=" + wbest + ") "
//               << worst.eps << "\n";
//     REQUIRE(std::abs(worst.eps) < 4.0e-5);
//     REQUIRE(std::abs(best.eps) < 6.0e-7);
//   }

//   // Since we have trouble with TDHF and HFS, do it again more thoroughly here.
//   // Note: This makes it seem the problem is in solve_core, NOT in solve dPsi
//   // This should be easier to fix!
//   // NOTE: This loop is very much the same as the other one... could combine

//   // Hyperfine
//   std::cout << "\nTest hyperfine (again, but more)\n";
//   auto &h = hhfs;

//   for (int max_its : {0, 1, 99}) {
//     double worst_eps = 0.0;
//     std::string worst_set{};
//     double best_eps = 1.0;
//     std::string best_set{};
//     ExternalField::TDHF dv(&h, wf.vHF());
//     dv.solve_core(0.0, max_its);

//     // to test the .get() X,Y's
//     ExternalField::TDHF dPsi(&h, wf.vHF());
//     dPsi.solve_core(0.0, 1); // 1 it; no dV, but solve for dPsi

//     int count = 0;
//     for (const auto &Fv : wf.valence()) {
//       for (const auto &Fm : wf.valence()) {
//         if (Fm == Fv || h.isZero(Fm.kappa(), Fv.kappa()))
//           continue;

//         const auto Xb =
//             dv.solve_dPsi(Fv, 0.0, ExternalField::dPsiType::X, Fm.kappa());
//         const auto Yb =
//             dv.solve_dPsi(Fv, 0.0, ExternalField::dPsiType::Y, Fm.kappa(),
//                           nullptr, ExternalField::StateType::bra);
//         const auto h_mv = h.reducedME(Fm, Fv) + dv.dV(Fm, Fv);
//         const auto lhs = Fm * Xb;
//         const auto rhs = h_mv / (Fv.en() - Fm.en());
//         const auto eps = (lhs - rhs) / (lhs + rhs);

//         const auto h_vm = h.reducedME(Fv, Fm) + dv.dV(Fv, Fm);
//         const auto lhsY = Yb * Fm;
//         const auto rhsY = h_vm / (Fv.en() - Fm.en());
//         const auto epsY = (lhsY - rhsY) / (lhsY + rhsY);

//         if (count % 5 == 0 || count == 0) {
//           // only print a subset (too many)
//           std::cout << "<" << Fv.shortSymbol() << "|h|" << Fm.shortSymbol()
//                     << "> , <" << Fv.shortSymbol() << "|X_" << Xb.shortSymbol()
//                     << "> : "; // << rhs << " " << lhs << "\n";
//           printf("%10.7g, %10.7g  %6.0e\n", rhs, lhs, eps);
//           std::cout << "<" << Fm.shortSymbol() << "|h|" << Fv.shortSymbol()
//                     << "> , <Y_" << Yb.shortSymbol() << "|" << Fv.shortSymbol()
//                     << "> : "; // << rhsY << " " << lhsY << "\n";

//           printf("%10.7g, %10.7g  %6.0e", rhsY, lhsY, epsY);
//           if (std::abs(epsY + eps) > 1.0e-3)
//             std::cout << "  ***";
//           std::cout << "\n";
//         }
//         ++count;

//         if (std::abs(eps) > std::abs(worst_eps)) {
//           worst_eps = eps;
//           worst_set = "<" + Fv.shortSymbol() + "|h|" + Fm.shortSymbol() +
//                       ">/<" + Fv.shortSymbol() + "|X_" + Xb.shortSymbol() + ">";
//         }
//         if (std::abs(epsY) > std::abs(worst_eps)) {
//           worst_eps = epsY;
//           worst_set = "<" + Fm.shortSymbol() + "|h|" + Fv.shortSymbol() +
//                       ">/<Y_" + Xb.shortSymbol() + "|" + Fv.shortSymbol() + ">";
//         }
//         if (std::abs(eps) < std::abs(best_eps)) {
//           best_eps = eps;
//           best_set = "<" + Fv.shortSymbol() + "|h|" + Fm.shortSymbol() + ">/<" +
//                      Fv.shortSymbol() + "|X_" + Xb.shortSymbol() + ">";
//         }
//         if (std::abs(epsY) < std::abs(best_eps)) {
//           best_eps = epsY;
//           best_set = "<" + Fm.shortSymbol() + "|h|" + Fv.shortSymbol() +
//                      ">/<Y_" + Xb.shortSymbol() + "|" + Fv.shortSymbol() + ">";
//         }
//       }
//     }
//     std::cout << worst_set << " " << worst_eps << "\n";
//     REQUIRE(std::abs(best_eps) < 1.0e-7);
//     REQUIRE(std::abs(worst_eps) < 1.0e-4);
//   }
// }

// //==============================================================================
// inline std::pair<helper::Case, helper::Case>
// helper::MS_loops(const Wavefunction &wf,
//                  const DiracOperator::TensorOperator *h) {
//   // compare the best and worst! (of each operator)
//   Case best{999.0};
//   Case worst{0.0};

//   const auto omega_mults = std::vector{0.0, 0.25};

// #pragma omp parallel for
//   for (auto i = 0ul; i < wf.valence().size(); ++i) {
//     const auto &Fv = wf.valence()[i];
//     const auto vl = wf.vlocal(Fv.l());
//     for (const auto &Fm : wf.valence()) {

//       // Only do for cases that make sense:
//       if (Fm == Fv)
//         continue; // gives 0/0
//       if (h->isZero(Fv.kappa(), Fm.kappa()))
//         continue;
//       if (Fm.kappa() != Fv.kappa() && std::abs(Fm.n() - Fv.n()) != 0)
//         continue;

//       const auto h_mv = h->reducedME(Fm, Fv);

//       // Don't do in cases where ME is extremely small. NOTE: This requires
//       // a good choice of units, since ME may be small due to small
//       // coeficient, which should not count as 'small'
//       if (std::abs(h_mv) < 1.0e-8) {
//         continue;
//       }

//       // Form 'rhs': (h - de_A)|A>
//       auto hFv = h->reduced_rhs(Fm.kappa(), Fv);
//       if (Fm.kappa() == Fv.kappa()) {
//         auto de = h->reducedME(Fv, Fv);
//         hFv -= de * Fv;
//       }

//       // loop over few frequencies:
//       for (const auto &w_mult : omega_mults) {
//         const auto w = std::abs(Fv.en() * w_mult);

//         const auto dFv = ExternalField::solveMixedState(Fv, w, vl, wf.alpha(),
//                                                         wf.core(), hFv);

//         const auto lhs = Fm * dFv;
//         const auto rhs = h_mv / (Fv.en() - Fm.en() + w);
//         const auto eps = (lhs - rhs) / (lhs + rhs);

// // find the best and worst case:
// #pragma omp critical(find_max)
//         {
//           if (std::abs(eps) > std::abs(worst.eps)) {
//             worst.eps = eps;
//             worst.name = Fm.shortSymbol() + "|" + Fv.shortSymbol();
//             worst.w = w;
//           }
//           if (std::abs(eps) < std::abs(best.eps)) {
//             best.eps = eps;
//             best.name = Fm.shortSymbol() + "|" + Fv.shortSymbol();
//             best.w = w;
//           }
//         }
//       }
//     }
//   }
//   return {best, worst};
// }
