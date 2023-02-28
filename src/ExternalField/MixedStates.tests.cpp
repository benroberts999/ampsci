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
