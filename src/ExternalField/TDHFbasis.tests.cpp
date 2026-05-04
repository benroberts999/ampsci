
#include "TDHFbasis.hpp"
#include "DiagramRPA.hpp"
#include "DiracOperator/include.hpp"
#include "TDHFbasis.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "calcMatrixElements.hpp"
#include "catch2/catch.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
TEST_CASE("External Field: TDHFbasis - basic unit tests",
          "[ExternalField][TDHFbasis][RPA][unit]") {

  Wavefunction wf({350, 1.0e-4, 80.0, 20.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", std::nullopt, "[Ne]", 1.0e-10, false);
  wf.solve_valence("3sp", false);
  // nb: use very small basis.
  // Don't care about numerical results, just that eveything is working correctly.
  SplineBasis::Parameters params{"9sp8d", 20, 7, 1.0e-3, 1.0e-3, 30.0};
  params.verbose = false;
  wf.formBasis(params);

  auto dE1 = DiracOperator::E1(wf.grid());

  auto rpa = ExternalField::TDHFbasis(&dE1, wf.vHF(), wf.basis());
  auto rpa2 =
    ExternalField::make_rpa("basis", &dE1, wf.vHF(), true, wf.basis(), "");

  rpa.solve_core(0.0, 18);
  rpa2->solve_core(0.0, 18);

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

  // dV_rhs(Fa) should return a DiracSpinor with the requested kappa
  for (int kappa : {-1, 1, -2, 2, -3}) {
    REQUIRE(rpa.dV_rhs(kappa, *F6s).kappa() == kappa);
    REQUIRE(rpa.dV_rhs(kappa, *F6p).kappa() == kappa);
  }

  // Fb * dV_rhs(kappa_b, Fa) should equal dV(Fb, Fa)
  {
    const auto dv_ps = rpa.dV(*F6p, *F6s);
    const auto dv_sp = rpa.dV(*F6s, *F6p);
    const auto rhs_ps = *F6p * rpa.dV_rhs(F6p->kappa(), *F6s);
    const auto rhs_sp = *F6s * rpa.dV_rhs(F6s->kappa(), *F6p);
    REQUIRE(rhs_ps == Approx(dv_ps));
    REQUIRE(rhs_sp == Approx(dv_sp));
  }

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
  wf.solve_core("HartreeFock", std::nullopt, "[Ar]", 1.0e-13, false);
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
        const auto dv_basis = Basis.dV(v, w);
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

//==============================================================================
// Breit-specific tests:
//==============================================================================
namespace helper {
std::pair<double, std::string>
do_dV_breit_basis(const Wavefunction &wf, const Wavefunction &wfB,
                  const DiracOperator::TensorOperator *const h,
                  double omega = 0.0, int max_l = 99);
} // namespace helper

//==============================================================================
//! Compared Breit contribution to dV using TDHF and TDHFbasis methods
TEST_CASE("External Field: TDHF (basis) Breit",
          "[ExternalField][TDHF][TDHFbasis][RPA][Breit][slow][integration]") {

  std::cout << "\nCompare Breit contribution to dV using TDHF and TDHFbasis "
               "methods\n"
               "Test of internal consistancy, not validation of Breit "
               "correction\n\n";
  // Test of internal consistancy, not validation of Breit correction

  // nb: don't use many grid points, since we don't need overall accuracy.
  // Need fair splines though, since splineTDHF compared to TDHF

  // Compares the Breit correction to dV using TDHF and TDHF basis methods.
  // NB: also a test of BASIS+Breit

  std::string valence = "6s5p4d";
  SplineBasis::Parameters bspl_param;
  {
    bspl_param.states = "40spd30f";
    bspl_param.n = 50;
    bspl_param.k = 7;
    bspl_param.r0 = 1.0e-3;
    bspl_param.reps = 0.0;
    bspl_param.rmax = 40.0;
  }
  bspl_param.verbose = false;

  // Create wavefunction object, solve HF for core+valence
  Wavefunction wf({1000, 1.0e-6, 120.0, 50.0, "loglinear", -1.0},
                  {"Rb", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", std::nullopt, "[Kr]", 1.0e-12, false); // no breit
  wf.solve_valence(valence, false);
  wf.formBasis(bspl_param);
  // wf.printValence();

  // Again, including Breit:
  Wavefunction wfB(wf.grid_sptr(), wf.nucleus(), 1.0);
  wfB.solve_core("HartreeFock", HF::Breit::Params{1.0}, "[Kr]", 1.0e-12, false); // w/Breit
  wfB.solve_valence(valence, false);
  wfB.formBasis(bspl_param);

  const auto hE1 = DiracOperator::E1(wf.grid());
  const auto hPNC =
    DiracOperator::PNCnsi(wf.nucleus().c(), wf.nucleus().t(), wf.grid());

  using namespace helper;
  const auto [eE1, sE1] = do_dV_breit_basis(wf, wfB, &hE1, 0.0);
  const auto [eE1w, sE1w] = do_dV_breit_basis(wf, wfB, &hE1, 0.05);
  const auto max_l = 1; // don't include d-states for pnc?
  const auto [ePNC, sPNC] = do_dV_breit_basis(wf, wfB, &hPNC, 0.0, max_l);

  std::cout << "Breit TDHF(bs) E1 w=0 " << sE1 << " " << eE1 << "\n";
  std::cout << "Breit TDHF(bs) E1 w=0.05 " << sE1w << " " << eE1w << "\n";
  std::cout << "Breit TDHF(bs) PNC " << sPNC << " " << ePNC << "\n";
  REQUIRE(std::abs(eE1) < 0.1);
  REQUIRE(std::abs(eE1w) < 0.1);
  // XXX Not sure why this one isn't great?
  REQUIRE(std::abs(ePNC) < 2.0);
}

//==============================================================================
std::pair<double, std::string>
helper::do_dV_breit_basis(const Wavefunction &wf, const Wavefunction &wfB,
                          const DiracOperator::TensorOperator *const h,
                          double omega, int max_l) {

  auto dV_tdhf = ExternalField::TDHF(h, wf.vHF());
  auto dV_basis = ExternalField::TDHFbasis(h, wf.vHF(), wf.basis());
  // w/ Breit:
  auto dVB_tdhf = ExternalField::TDHF(h, wfB.vHF());
  auto dVB_basis = ExternalField::TDHFbasis(h, wfB.vHF(), wfB.basis());

  // double omega = 0.0;
  const auto max_iterations = 50;
  const auto print_details = true;
  std::cout << "No Breit:\n";
  dV_tdhf.solve_core(omega, max_iterations, print_details);
  std::cout << "With Breit:\n";
  dVB_tdhf.solve_core(omega, max_iterations, print_details);
  std::cout << "No Breit:\n";
  dV_basis.solve_core(omega, max_iterations, print_details);
  std::cout << "With Breit:\n";
  dVB_basis.solve_core(omega, max_iterations, print_details);

  auto eps = 0.0;
  std::string worst = "";

  std::cout << "\n(Breit correction (%) to the dV correction)\n";
  std::cout << "dBr(%)  TDHF       TDHF(basis) " << h->name() << "\n";
  for (const auto &Fv : wf.valence()) {
    for (const auto &Fw : wf.valence()) {
      if (Fw > Fv || h->isZero(Fv.kappa(), Fw.kappa()) || Fv.l() > max_l ||
          Fw.l() > max_l)
        continue;

      // lookup states in Breit wf (do this way, so no depend on order)
      const auto &FvB = *wfB.getState(Fv.n(), Fv.kappa());
      const auto &FwB = *wfB.getState(Fw.n(), Fw.kappa());

      const auto dv = dV_tdhf.dV(Fv, Fw);
      const auto dvB = dVB_tdhf.dV(FvB, FwB);
      const auto dv_basis = dV_basis.dV(Fv, Fw);
      const auto dvB_basis = dVB_basis.dV(FvB, FwB);

      const auto dBr = (dvB - dv) / dv;
      const auto dBr_basis = (dvB_basis - dv_basis) / dv_basis;
      const auto str = Fv.shortSymbol() + "-" + Fw.shortSymbol();
      printf("%7s %+9.4e %+9.4e\n", str.c_str(), dBr * 100.0,
             dBr_basis * 100.0);
      const auto eps_this =
        2.0 * std::abs((dBr - dBr_basis) / (dBr + dBr_basis));
      if (eps_this > eps) {
        eps = eps_this;
        worst = str;
      }
    }
  }
  std::cout << "\n";
  return {eps, worst};
}
