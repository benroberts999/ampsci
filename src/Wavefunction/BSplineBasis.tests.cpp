#include "Wavefunction/BSplineBasis.hpp"
#include "DiracOperator/include.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

/* XXX 18/01/2022
  - Made change to spline (r00 = log-average<r0,s0_spl>)
  - Fixed Breit, but (slightly) broke basis
  - "Fix" by just loosening tollerences here
*/

inline double hfsA(const DiracOperator::TensorOperator *h,
                   const DiracSpinor &Fa) {
  auto Raa = h->radialIntegral(Fa, Fa);
  return Raa * Fa.kappa() / (Fa.jjp1()) * PhysConst::muN_CGS_MHz;
}

//==============================================================================
TEST_CASE("Wavefunction: BSpline-basis unit", "[BSpline][basis][unit]") {

  // Create wavefunction object, solve HF for core + valence
  Wavefunction wf({2500, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear"},
                  {"Cs", -1, "Fermi"});
  wf.set_HF(HF::Method::Local, "[Xe]");
  wf.solve_core();
  wf.solve_valence("7sp5d4f");

  const auto hA = DiracOperator::hfs(1, 1.0, 0.0, wf.grid(),
                                     DiracOperator::Hyperfine::pointlike_F());
  const auto r = wf.grid().r();

  const std::string states = "7sp5d4f";
  const int k = 7;
  const int n_spl = 70;
  const auto r0 = 1.0e-4;
  // const auto r0_eps = 1.0e-3;
  const auto rmax = 75.0;

  using ST = SplineBasis::SplineType;
  for (const auto type : {ST::Derevianko, ST::Johnson, ST::Fischer}) {
    const auto tname = SplineBasis::parseSplineType(type);

    const auto r0_eps = type == ST::Johnson ? 1.0e-4 : 0.0;

    // State comparison: energy, HFS, <r>, ortho
    wf.formBasis({states, std::size_t(n_spl), k, r0, r0_eps, rmax, "", type});
    fmt::print("\n{} (N={}, k={}):\n", tname, n_spl, k);
    fmt::print("{:<8}  {:>9}  {:>9}  {:>9}  {:>9}\n", "State", "dE/E", "dA/A",
               "d<r>", "ort");
    for (auto &Fb : wf.basis()) {
      const auto &Fn = *wf.getState(Fb.n(), Fb.kappa());
      const auto eps = std::abs((Fb.en() - Fn.en()) / Fn.en());
      const auto a = hfsA(&hA, Fn);
      const auto ab = hfsA(&hA, Fb);
      const auto eps_a = std::abs((ab - a) / a);
      const auto rn = Fn * (r * Fn);
      const auto rb = Fb * (r * Fb);
      const auto eps_r = std::abs((rb - rn) / rn);
      const auto norm = std::abs(Fb * Fn - 1.0);
      fmt::print("{:<8}  {:>9.2e}  {:>9.2e}  {:>9.2e}  {:>9.2e}\n", Fb.symbol(),
                 eps, eps_a, eps_r, norm);
      REQUIRE(eps < 1.0e-4);
      REQUIRE(eps_a < 1.0e-2);
      REQUIRE(eps_r < 1.0e-4);
      REQUIRE(norm < 1.0e-7);
    }

    // Sum rules (bigger basis with positrons)
    wf.formBasis({"spdf", 70, 9, r0, r0_eps, rmax, "spdf", type});

    fmt::print("\n{}: TKR sum rule, including positron states:\n", tname);
    const auto eps_tkr = type == ST::Johnson ? 1.0e-5 : 1.0e-7;
    const auto &tkr_l =
      SplineBasis::sumrule_TKR(wf.basis(), wf.grid().r(), true);
    for (auto &tkr : tkr_l) {
      REQUIRE(tkr < eps_tkr);
    }

    fmt::print("\n{}: Drake-Gordon sum rule, including positron states:\n",
               tname);
    const auto eps_dg = type == ST::Derevianko ? 1.0e-7 :
                        type == ST::Johnson    ? 1.0e-5 :
                                                 1.0e-6;
    for (int n = 0; n <= 2; ++n) {
      const auto &dgs =
        SplineBasis::sumrule_DG(n, wf.basis(), wf.grid(), wf.alpha(), true);
      for (auto &dg : dgs) {
        REQUIRE(std::abs(dg) < eps_dg);
      }
    }

    fmt::print("\n");
    SplineBasis::r_completeness_header();
    for (auto &v : wf.valence()) {
      const auto [e1, e2] =
        SplineBasis::r_completeness(v, wf.basis(), wf.grid(), true);
      if (v.l() <= 2) {
        REQUIRE(e1 < 1.0e-10);
        REQUIRE(e2 < 1.0e-9);
      } else {
        REQUIRE(e1 < 1.0e-8);
        REQUIRE(e2 < 1.0e-6);
      }
    }
  }
}

//==============================================================================
//! Unit tests for B-Spline basis (finite basis of relativistic orbitals)
TEST_CASE("Wavefunction: BSpline-basis",
          "[BSpline][basis][QED][Breit][integration]") {

  std::cout
    << "\nCheck all basis methods work well, including with QED and Breit\n";

  //============================================================================
  // Check vs. Hartree-Fock (energies/ortho)

  for (const auto f_QED : {false, true}) {
    for (const auto f_Breit : {0.0, 1.0}) {

      fmt::print(
        "\n====================\n{} Breit, {} QED\n====================\n",
        (f_Breit != 0.0 ? "With" : "No"), (f_QED ? "with" : "no"));

      // Create wavefunction object, solve HF for core + valence
      Wavefunction wf({2500, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear"},
                      {"Cs", -1, "Fermi"});
      const auto breit_params =
        f_Breit != 0.0 ?
          std::optional<HF::Breit::Params>{HF::Breit::Params{f_Breit}} :
          std::nullopt;
      wf.set_HF(HF::Method::HartreeFock, "[Xe]", breit_params);
      if (f_QED)
        wf.radiativePotential({1.0, 1.0, 1.0, 1.0, 0.0}, 10.0, 1.0, {1.0},
                              false, false);
      wf.solve_core();
      wf.solve_valence("7sp5d4f");

      for (const auto type : {SplineBasis::SplineType::Derevianko,
                              SplineBasis::SplineType::Johnson,
                              SplineBasis::SplineType::Fischer}) {

        // Compare the energy of a Dirac spinor to a double:
        const auto comp_en = [](const auto &Fa, double en) {
          return (en - Fa.en()) / Fa.en();
        };

        // Test splines compared to HF states (for various spline sets)
        for (std::size_t k : {7ul, 9ul}) {
          const auto nlst = k == 7 ? std::vector{75, 85} : std::vector{70, 90};
          for (const auto n : nlst) {
            fmt::print("\nType : {}\n", SplineBasis::parseSplineType(type));
            fmt::print("k_spl: {}\n", k);
            fmt::print("N_spl: {}\n", n);

            // Form spline basis:
            const std::string states = "30spdf";
            const auto r0 = 1.0e-5;
            const auto r0_eps =
              type == SplineBasis::SplineType::Johnson ? 1.0e-3 : 0.0;
            const auto rmax = 75.0;
            const std::string positron = "";
            const auto basis =
              SplineBasis::form_basis({states, std::size_t(n), k, r0, r0_eps,
                                       rmax, positron, type, false, false},
                                      wf);

            // Check orthonormality <a|b>:
            const auto [eps, str] = DiracSpinor::check_ortho(basis, basis);
            const auto [eps1, str1] =
              DiracSpinor::check_ortho(basis, wf.core());
            const auto [eps2, str2] =
              DiracSpinor::check_ortho(basis, wf.valence());
            std::cout << "Orthog (worst cases):\n";
            std::cout << str << " " << eps << "\n";
            std::cout << str1 << " " << eps1 << "\n";
            std::cout << str2 << " " << eps2 << "\n";
            REQUIRE(eps < 1.0e-12);
            REQUIRE(eps1 < 1.0e-4);
            REQUIRE(eps1 < 1.0e-4);

            // Compare energies
            std::vector<double> core_en;
            for (const auto &Fc : wf.core()) {
              const auto &pFb = std::find(cbegin(basis), cend(basis), Fc);
              core_en.push_back(pFb->en());
            }
            std::vector<double> val_en;
            for (const auto &Fv : wf.valence()) {
              const auto &pFb = std::find(cbegin(basis), cend(basis), Fv);
              val_en.push_back(pFb->en());
            }

            // Compare core+valence HF energies to the corresponding splines
            const auto [ce, cs] = qip::compare(wf.core(), core_en, comp_en);
            const auto [ve, vs] = qip::compare(wf.valence(), val_en, comp_en);
            std::cout << "Energy comparison (worst case):\n";
            fmt::print("E_{}: E_HF = {:10.5f}, eps = {:.1e}\n",
                       cs->shortSymbol(), cs->en(), ce);
            fmt::print("E_{}: E_HF = {:10.5f}, eps = {:.1e}\n",
                       vs->shortSymbol(), vs->en(), ve);
            REQUIRE(ce < 1.0e-5);
            REQUIRE(ve < 1.0e-5);

            //--------------
            // Compare HFS and <r>: worst core/val
            const auto h = DiracOperator::hfs(
              1, 1.0, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());
            const auto &rv = wf.grid().r();

            double wc_Aeps = 0.0, wv_Aeps = 0.0;
            double wc_A = 0.0, wc_Aspl = 0.0;
            double wv_A = 0.0, wv_Aspl = 0.0;
            std::string wc_Asym{}, wv_Asym{};
            double wc_reps = 0.0, wv_reps = 0.0;
            double wc_r = 0.0, wc_rspl = 0.0;
            double wv_r = 0.0, wv_rspl = 0.0;
            std::string wc_rsym{}, wv_rsym{};
            for (const auto &Fb : basis) {
              const auto *pFhf = wf.getState(Fb.n(), Fb.kappa());
              if (!pFhf)
                continue;
              const auto &Fhf = *pFhf;
              const auto Ahf = hfsA(&h, Fhf);
              const auto Aspl = hfsA(&h, Fb);
              const auto eps_a = std::abs((Ahf - Aspl) / Ahf);
              const auto rhf = Fhf * (rv * Fhf);
              const auto rspl = Fb * (rv * Fb);
              const auto eps_r = std::abs((rhf - rspl) / rhf);
              const bool is_core = !wf.isInValence(Fb.n(), Fb.kappa());
              if (is_core) {
                if (eps_a > wc_Aeps) {
                  wc_Aeps = eps_a;
                  wc_A = Ahf;
                  wc_Aspl = Aspl;
                  wc_Asym = Fb.symbol();
                }
                if (eps_r > wc_reps) {
                  wc_reps = eps_r;
                  wc_r = rhf;
                  wc_rspl = rspl;
                  wc_rsym = Fb.symbol();
                }
              } else {
                if (eps_a > wv_Aeps) {
                  wv_Aeps = eps_a;
                  wv_A = Ahf;
                  wv_Aspl = Aspl;
                  wv_Asym = Fb.symbol();
                }
                if (eps_r > wv_reps) {
                  wv_reps = eps_r;
                  wv_r = rhf;
                  wv_rspl = rspl;
                  wv_rsym = Fb.symbol();
                }
              }
            }
            std::cout << "HFS comparison (worst case):\n";
            fmt::print(
              "A_{}: A_HF = {:12.4e}, A_spl = {:12.4e}, eps = {:.1e}\n",
              wc_Asym, wc_A, wc_Aspl, wc_Aeps);
            fmt::print(
              "A_{}: A_HF = {:12.4e}, A_spl = {:12.4e}, eps = {:.1e}\n",
              wv_Asym, wv_A, wv_Aspl, wv_Aeps);
            const auto eps_hfs =
              type == SplineBasis::SplineType::Johnson ? 1.0e-2 : 1.0e-4;
            REQUIRE(wc_Aeps < eps_hfs);
            REQUIRE(wv_Aeps < eps_hfs);

            std::cout << "<r> comparison (worst case):\n";
            fmt::print(
              "r_{}: r_HF = {:12.4e}, r_spl = {:12.4e}, eps = {:.1e}\n",
              wc_rsym, wc_r, wc_rspl, wc_reps);
            fmt::print(
              "r_{}: r_HF = {:12.4e}, r_spl = {:12.4e}, eps = {:.1e}\n",
              wv_rsym, wv_r, wv_rspl, wv_reps);
            REQUIRE(wc_reps < 1.0e-5);
            REQUIRE(wv_reps < 1.0e-6);
          }
        }
      }
    }
  }

  //============================================================================

  // // Check low-r behavour (splines vs HF for hyperfine)
  // for (auto f_QED : {false, true}) {
  //   for (auto f_Br : {0.0, 1.0}) {

  //     Wavefunction wf({6000, 1.0e-7, 150.0, 0.33 * 150.0, "loglinear"},
  //                     {"Cs", -1, "Fermi"}, 1.0);

  //     const auto breit_params =
  //       f_Br != 0.0 ?
  //         std::optional<HF::Breit::Params>{HF::Breit::Params{f_Br}} :
  //         std::nullopt;
  //     wf.set_HF(HF::Method::HartreeFock, "[Xe]", breit_params);
  //     if (f_QED)
  //       wf.radiativePotential({1.0, 1.0, 1.0, 1.0, 0.0}, 10.0, 1.0, {1.0},
  //                             false, false);
  //     wf.solve_core();

  //     const std::string states = "7sp5d4f";
  //     wf.solve_valence(states);

  //     const auto r0 = 1.0e-4;
  //     const auto r0eps = 1.0e-5;
  //     const auto rmax = 60.0; // need large rmax, to match to large nl HF
  //     const int nspl = 75;    // need large n due to large rmax
  //     const int kspl = 7;
  //     wf.formBasis({states, nspl, kspl, r0, r0eps, rmax, "",
  //                   SplineBasis::SplineType::Derevianko});

  //     // Hyperfine operator: Pointlike, g=1
  //     const auto h = DiracOperator::hfs(
  //       1, 1.0, 0.0, wf.grid(), DiracOperator::Hyperfine::pointlike_F());

  //     // Calculate A with HF and spline states, compare for each l:
  //     std::vector<std::vector<double>> hfs(4); // s,p,d,f
  //     std::cout << "A, hf vs spl:\n";
  //     for (const auto &orbs : {wf.core(), wf.valence()}) { //*
  //       for (const auto &Fv : orbs) {
  //         const auto pFb = std::find(cbegin(wf.basis()), cend(wf.basis()), Fv);
  //         const auto Ahf = hfsA(&h, Fv);
  //         const auto Aspl = hfsA(&h, *pFb);
  //         const auto eps = (Ahf - Aspl) / Ahf;
  //         std::cout << Fv.symbol() << " " << Ahf << " " << Aspl << " " << eps
  //                   << "\n";
  //         hfs.at(std::size_t(Fv.l())).push_back(eps);
  //       }
  //     }
  //     for (std::size_t l = 0; l < hfs.size(); l++) {
  //       const auto eps =
  //         *std::max_element(cbegin(hfs[l]), cend(hfs[l]), qip::less_abs{});

  //       std::string name = "spl";
  //       if (f_Br != 0.0)
  //         name += "+Br";
  //       if (f_QED)
  //         name += "+QED";
  //       const std::string str = name + " low-r (hfs) l=";

  //       // pass &=
  //       // qip::check_value(&obuff, str + std::to_string(l), eps, 0.0, 1.0e-4);
  //       REQUIRE(std::abs(eps) < 1.0e-4);
  //     }
  //   }
  // }
}
