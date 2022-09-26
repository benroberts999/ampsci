#include "DiracOperator/DiracOperator.hpp"
#include "TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "calcMatrixElements.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <vector>

TEST_CASE("External Field: calcMatrixElements", "[ExternalField][unit]") {
  std::cout << "\n----------------------------------------\n";
  std::cout << "External Field: TDHF - calcMatrixElements\n";

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Li", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[He]");
  wf.solve_valence("2sp");

  auto dE1 = DiracOperator::E1(wf.grid());
  auto rpa = ExternalField::TDHF(&dE1, wf.vHF());
  // rpa.solve_core(0.0);

  auto mes =
      ExternalField::calcMatrixElements(wf.valence(), &dE1, &rpa, 0.0, false);

  for (auto &me : mes) {
    std::cout << me << "\n";
    auto [a, b, h, dv1, dv] = me;

    auto Fa = wf.getState(a);
    auto Fb = wf.getState(b);
    REQUIRE(Fa != nullptr);
    REQUIRE(Fb != nullptr);
    auto h0 = dE1.reducedME(*Fa, *Fb);
    REQUIRE(h == Approx(h0));

    auto hdv10 = rpa.dV1(*Fa, *Fb);
    REQUIRE(dv1 == Approx(hdv10));

    auto hdv0 = rpa.dV(*Fa, *Fb);
    REQUIRE(dv == Approx(hdv0));
  }
}