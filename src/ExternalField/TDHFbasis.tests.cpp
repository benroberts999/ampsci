
#include "TDHFbasis.hpp"
#include "DiagramRPA.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "calcMatrixElements.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
TEST_CASE("External Field: TDHFbasis - basic unit tests",
          "[ExternalField][TDHFbasis][RPA][unit]") {

  Wavefunction wf({500, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ne]", 1.0e-10, false);
  wf.solve_valence("3sp", false);
  // nb: use very small basis.
  // Don't care about numerical results, just that eveything is working correctly.
  SplineBasis::Parameters params{"10sp8d", 20, 7, 1.0e-3, 1.0e-3, 30.0};
  params.verbose = false;
  wf.formBasis(params);

  auto dE1 = DiracOperator::E1(wf.grid());

  auto rpa = ExternalField::TDHFbasis(&dE1, wf.vHF(), wf.basis());
  auto rpa2 =
    ExternalField::make_rpa("basis", &dE1, wf.vHF(), true, wf.basis(), "");

  rpa.solve_core(0.0, 15);
  rpa2->solve_core(0.0, 15);

  for (const auto &v : wf.valence()) {
    for (const auto &w : wf.valence()) {

      const auto dv_0 = rpa.dV(v, w);
      const auto dv_1 = rpa.dV(w, v) * dE1.symm_sign(v, w);

      const auto dv_2 = rpa2->dV(v, w);

      if (!dE1.isZero(v, w))
        REQUIRE(dv_0 != 0.0);
      REQUIRE(dv_0 == Approx(dv_1));
      REQUIRE(dv_0 == Approx(dv_2));
    }
  }

  for (const auto &Fc : wf.core()) {
    auto dPsis1 =
      rpa.form_dPsis(Fc, 0.0, ExternalField::dPsiType::X, wf.basis());
    auto dPsis2 = rpa.get_dPsis(Fc, ExternalField::dPsiType::X);
    REQUIRE(dPsis1.size() == dPsis2.size());
    REQUIRE(dPsis1.size() != 0);
    for (std::size_t i = 0; i < dPsis1.size(); ++i) {
      REQUIRE(dPsis1.at(i).norm2() == Approx(dPsis2.at(i).norm2()));
    }
  }

  auto F6s = wf.getState("3s");
  auto F6p = wf.getState("3p-");
  REQUIRE(F6s != nullptr);
  REQUIRE(F6p != nullptr);

  rpa.clear();
  const auto dv_00 = rpa.dV(*F6s, *F6p);
  REQUIRE(dv_00 == 0.0);
}

//==============================================================================
//==============================================================================

TEST_CASE("External Field: TDHFbasis (RPA)",
          "[ExternalField][TDHFbasis][DiagramRPA][RPA][integration]") {

  // Create wavefunction object, solve HF for core+valence
  Wavefunction wf({1600, 1.0e-6, 130.0, 10.0, "loglinear", -1.0},
                  {"K", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Ar]", 1.0e-13, false);
  wf.solve_valence("4sp3d4f", false);

  // Use reasonably small basis here
  // Use Diagram method as benchmark
  // Should be OK - since Diagram method is well tested
  SplineBasis::Parameters params{"25spdfg", 40, 7, 1.0e-4, 1.0e-5, 40.0};
  params.verbose = false;
  wf.formBasis(params);

  // To compare, should only include excited in the basis
  // (Makes no noticable difference)
  const auto [core, excited] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());

  auto h1 = DiracOperator::E1(wf.grid());
  auto h2 = DiracOperator::Ek(wf.grid(), 2);
  auto h3 = DiracOperator::hfs(1, 1.0, 1.0, wf.grid());

  std::cout << "\nCompare Basis to Diagram methods\n";

  for (const DiracOperator::TensorOperator *h :
       std::array<const DiracOperator::TensorOperator *, 3>{&h1, &h2, &h3}) {

    auto Basis = ExternalField::TDHFbasis(h, wf.vHF(), excited);
    auto Diagram =
      ExternalField::DiagramRPA(h, wf.basis(), wf.vHF(), "", false);

    std::cout << "\n" << h->name() << ":\n";
    double omega = h->name() == "E1" ? 0.1 : 0.0;
    Basis.solve_core(omega);
    Diagram.solve_core(omega);

    REQUIRE(Basis.last_eps() < 1.0e-10);
    REQUIRE(Diagram.last_eps() < 1.0e-10);
    REQUIRE(Basis.last_omega() == Approx(omega));

    fmt::print("{:3} {:3}  {:11}  [{:11}]   ({:5})\n", "v", "w", "Basis",
               "Diagram", "eps");
    for (const auto &v : wf.valence()) {
      for (const auto &w : wf.valence()) {
        if (h->isZero(v, w))
          continue;
        // False: prevent calculateing dV^dag for v>w ?
        // Because Basis inherits from TDHF, which distingushes dV and dV^dag
        const auto dv_basis = Basis.dV(v, w, false);
        const auto dv_diagram = Diagram.dV(v, w);
        const auto eps = std::abs((dv_basis - dv_diagram) / dv_diagram);
        fmt::print("{:3} {:3}  {:+.4e}  [{:+.4e}]   ({:.0e})\n",
                   v.shortSymbol(), w.shortSymbol(), dv_basis, dv_diagram, eps);

        const auto eps_target = h->name() == "hfs1" ? 1.0e-3 : 1.0e-4;
        CHECK(dv_basis == Approx(dv_diagram).epsilon(eps_target));
      }
    }
  }
}
