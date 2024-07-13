#include "MixedStates.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHF.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "IO/ChronoTimer.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/String.hpp"
#include <iomanip>
#include <string>

//==============================================================================
TEST_CASE("External Field: Mixed-states (unit)",
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
            << fmt::format("{:.1e}", weps) << ", best: " << best << " "
            << fmt::format("{:.1e}", beps) << "\n";
  }
  std::cout << "\nMixed-States summary:\n" << summary.str() << "\n";
}

//==============================================================================
TEST_CASE("External Field: Mixed-states (full)",
          "[ExternalField][MixedStates][integration]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: Mixed-states - full\n";

  // Create wavefunction object
  Wavefunction wf({4000, 1.0e-6, 90.0, 20.0, "loglinear", -1.0},
                  {"Rb", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Kr]");
  wf.solve_valence("5sp");

  SplineBasis::Parameters bspl_param;
  {
    bspl_param.states = "85spdfg";
    bspl_param.n = 90;
    bspl_param.k = 9;
    bspl_param.r0 = 1.0e-5;
    bspl_param.reps = 0.0;
    bspl_param.rmax = 75.0;
    bspl_param.positronQ = false;
  }

  wf.formBasis(bspl_param);

  for (std::string_view oper : {"E1", "hfs"}) {
    std::cout << "\n" << oper << "\n";
    auto h = DiracOperator::generate(oper, {}, wf);

    // const auto r = wf.grid().r();

    double omega = 0.00;

    // <a|Xn>
    fmt::print("\n{:3s} {:3s} {:>11s} {:>11s} {:>11s} | {:>12s} {:>12s}  eps\n",
               "n", "a", "<a|h|n>", "TDHF", "basis", "<a|Xn>_TDHF",
               "<a|Xn>_basis");

    const auto states = qip::merge(wf.core(), wf.valence());

    // use basis states,

    for (auto Fm : states) {
      for (auto Fc : wf.core()) {

        // Use basis states, instead of HF, for "exact" orthogonality
        // Fm = *DiracSpinor::find(Fm.n(), Fm.kappa(), wf.basis());
        // Fc = *DiracSpinor::find(Fc.n(), Fc.kappa(), wf.basis());

        if (Fc > Fm || h->isZero(Fm, Fc))
          continue;

        const auto hFc = h->reduced_rhs(Fm.kappa(), Fc);
        auto imag = h->imaginaryQ();
        auto de = Fc.kappa() == Fm.kappa() && !imag ? Fc * hFc : 0.0;

        // Solves: (H - e - de)*dF = -h * F
        const auto &dF_t =
            ExternalField::solveMixedState(Fc, omega, hFc - de * Fc, wf.vHF());

        // Directly finds: dF = \sum_n |n><n||h|c> / (ec - en + w)
        const auto &dF_b =
            ExternalField::solveMixedState_basis(Fc, hFc, omega, wf.basis());

        auto h0 = Fm * hFc;
        auto h1 = (Fm * dF_t) * (Fc.en() - Fm.en() + omega);
        auto h2 = (Fm * dF_b) * (Fc.en() - Fm.en() + omega);

        if (Fm == Fc) {
          // Fm*dF is not defined in this case!
          // Should be zero exactly from sum-over-state method - isn't, because
          // states aren't exactly orthogonal. Should diverge in TDHF case?
          fmt::print(
              "{:3s} {:3s} {:11.4e} {:^11s} {:^11s} |  {:11.4e}  {:11.4e}\n",
              Fc.shortSymbol(), Fm.shortSymbol(), h0, "-----", "-----",
              Fm * dF_t, Fm * dF_b);
        } else {
          const auto eps =
              0.5 * (std::abs(h1 / h0 - 1.0) + std::abs(h2 / h0 - 1.0));
          fmt::print("{:3s} {:3s} {:11.4e} {:11.4e} {:11.4e} |  {:11.4e}  "
                     "{:11.4e}  {:.0e}\n",
                     Fc.shortSymbol(), Fm.shortSymbol(), h0, h1, h2, Fm * dF_t,
                     Fm * dF_b, eps);

          // base target on overlap of wavefunctions:
          const auto dn = std::abs(Fc.n() - Fm.n());
          const auto target = dn == 0 ? 1.0e-4 :
                              dn == 1 ? 1.0e-3 :
                              dn == 2 ? 1.0e-3 :
                                        1.0e-2;
          REQUIRE(eps < target);
        }
      }
      std::cout << "\n";
    }
  }
}