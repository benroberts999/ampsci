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
