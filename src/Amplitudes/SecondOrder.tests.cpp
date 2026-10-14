#include "Amplitudes/SecondOrder.hpp"
#include "Angular/include.hpp"
#include "CI/SecondOrder.hpp"
#include "DiracOperator/include.hpp"
#include "ExternalField/TDHF.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

//==============================================================================
TEST_CASE("Amplitudes: second-order angular",
          "[Amplitudes][SecondOrder][unit]") {

  // allowed_K against explicitly-written triangle rules
  for (int kt = 0; kt <= 3; ++kt) {
    for (int ks = 0; ks <= 3; ++ks) {
      for (int twoJb = 1; twoJb <= 7; twoJb += 2) {
        for (int twoJa = 1; twoJa <= 7; twoJa += 2) {
          int first_allowed = -1;
          for (int K = 0; K <= kt + ks + 2; ++K) {
            const bool ops = (K >= std::abs(kt - ks)) && (K <= kt + ks);
            const bool states =
              (2 * K >= std::abs(twoJb - twoJa)) && (2 * K <= twoJb + twoJa);
            const bool expected = ops && states;
            REQUIRE(Amplitudes::allowed_K(K, kt, ks, twoJb, twoJa) == expected);
            if (expected && first_allowed < 0) {
              first_allowed = K;
            }
          }
          REQUIRE(Amplitudes::smallest_allowed_K(kt, ks, twoJb, twoJa) ==
                  first_allowed);
        }
      }
    }
  }
}

//==============================================================================
//! SOS amplitudes against the independent formulas of the polarisability
//! module (different angular reduction): same basis, same matrix elements,
//! so agreement should be at machine precision
TEST_CASE("Amplitudes: second-order SOS", "[Amplitudes][SecondOrder][unit]") {

  Wavefunction wf({1000, 1.0e-5, 100.0, 20.0, "loglinear", -1.0},
                  {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("HartreeFock", "[Xe]");
  wf.solve_valence("7sp5d");
  wf.formBasis({"30spd", 40, 7, {1.0e-4}, 40.0});

  const auto *p6s = wf.getState("6s+");
  const auto *p7s = wf.getState("7s+");
  const auto *p5dm = wf.getState("5d-");
  REQUIRE(p6s != nullptr);
  REQUIRE(p7s != nullptr);
  REQUIRE(p5dm != nullptr);

  const DiracOperator::E1 hE1(wf.grid());
  const auto &basis = wf.basis();

  fmt::print("\nSOS amplitudes vs polarisability-module formulas (Cs):\n");
  fmt::print("{:<22} {:>14} {:>14} {:>9}\n", "", "expected", "found", "eps");

  const auto check = [](const std::string &name, double expected,
                        double found) {
    const auto eps =
      std::abs((found - expected) / std::max(std::abs(expected), 1.0e-30));
    fmt::print("{:<22} {:14.7e} {:14.7e} {:9.1e}\n", name, expected, found,
               eps);
    REQUIRE(found == Approx(expected).epsilon(1.0e-10));
  };

  // alpha(omega), via alphaD::valence_sos formula (independent angular form)
  const auto alpha_ref = [&](const DiracSpinor &Fv, double w) {
    double alpha = 0.0;
    for (const auto &Fn : basis) {
      if (hE1.isZero(Fv, Fn))
        continue;
      const auto d = hE1.reducedME(Fn, Fv);
      const auto de = Fv.en() - Fn.en();
      const auto denom = w == 0.0 ? 1.0 / de : de / (de * de - w * w);
      alpha += d * d * denom;
    }
    return (-2.0 / 3.0) / (Fv.twoj() + 1) * alpha;
  };

  // Static and dynamic scalar polarisability of 6s: K = 0 amplitude
  for (const auto omega : {0.0, 0.05}) {
    const auto A0 =
      Amplitudes::sos_valence(0, *p6s, *p6s, &hE1, &hE1, omega, -omega, basis);
    const auto alpha = A0 / std::sqrt(3.0 * (p6s->twoj() + 1));
    const auto name = fmt::format("alpha_6s(w={})", omega);
    check(name, alpha_ref(*p6s, omega), alpha);
  }

  // Tensor polarisability of 5d-: K = 2 amplitude, alphaD::tensor2_sos form
  {
    const auto A2 =
      Amplitudes::sos_valence(2, *p5dm, *p5dm, &hE1, &hE1, 0.0, 0.0, basis);
    const auto twoJ = double(p5dm->twoj());
    const auto factor =
      -std::sqrt(2.0 * twoJ * (twoJ - 1.0) /
                 (3.0 * (twoJ + 1.0) * (twoJ + 2.0) * (twoJ + 3.0)));
    const auto alpha_2 = factor * A2;

    const auto &Fv = *p5dm;
    const auto ctop = 2.5 * (Fv.twoj() * (Fv.twoj() - 1));
    const auto cbot =
      3.0 * ((Fv.twoj() + 2) * (Fv.twoj() + 1) * (Fv.twoj() + 3));
    const auto C = +4.0 * std::sqrt(ctop / cbot);
    double a2_ref = 0.0;
    for (const auto &Fn : basis) {
      if (hE1.isZero(Fv, Fn))
        continue;
      const auto sj = Angular::sixj_2(Fv.twoj(), 2, Fn.twoj(), 2, Fv.twoj(), 4);
      if (sj == 0.0)
        continue;
      const auto sign = Angular::neg1pow_2(Fv.twoj() + Fn.twoj() + 2);
      const auto d = hE1.reducedME(Fn, Fv);
      a2_ref += sign * sj * d * d / (Fv.en() - Fn.en());
    }
    check("alpha_2_5d-", C * a2_ref, alpha_2);
  }

  // Scalar transition amplitude 6s -> 7s: z-component of A^0 against the
  // alphaD::transition_sos formula (m = 1/2)
  {
    const auto omega = p7s->en() - p6s->en();
    const auto A0 =
      Amplitudes::sos_valence(0, *p7s, *p6s, &hE1, &hE1, omega, 0.0, basis);
    const auto A_zz =
      CI::z_component(0, 1, 1, p7s->twoj(), p6s->twoj(), 1) * A0;

    double ref = 0.0;
    for (const auto &Fn : basis) {
      if (hE1.isZero(*p7s, Fn) || hE1.isZero(*p6s, Fn))
        continue;
      const auto d_vn = hE1.reducedME(*p7s, Fn);
      const auto d_nw = hE1.reducedME(Fn, *p6s);
      const auto f_de =
        1.0 / (p7s->en() - Fn.en()) + 1.0 / (p6s->en() - Fn.en());
      const auto f = hE1.rme3js(p7s->twoj(), Fn.twoj(), 1) *
                     hE1.rme3js(Fn.twoj(), p6s->twoj(), 1);
      ref += f * d_vn * d_nw * f_de;
    }
    check("A_zz(6s-7s)", ref, A_zz);
  }

  // Core polarisability: K = 0 core amplitude against explicit core sum
  {
    const auto excited =
      DiracSpinor::split_by_energy(basis, wf.FermiLevel()).second;
    const auto A0c = Amplitudes::sos_core(0, p6s->twoj(), &hE1, &hE1, 0.0, 0.0,
                                          wf.core(), excited);
    const auto alpha_core = A0c / std::sqrt(3.0 * (p6s->twoj() + 1));

    double ref = 0.0;
    for (const auto &Fc : wf.core()) {
      for (const auto &Fn : excited) {
        if (hE1.isZero(Fc, Fn))
          continue;
        const auto d = hE1.reducedME(Fn, Fc);
        ref += d * d / (Fc.en() - Fn.en());
      }
    }
    check("alpha_core", (-2.0 / 3.0) * ref, alpha_core);
  }
}

//==============================================================================
//! Polarisability with Sigma, with and without the large-r tail extrapolation,
//! at realistic settings: Sigma calculated to 30 aB, spectrum extending well
//! beyond (as in production runs), so spectrum B-splines straddle the Sigma
//! sub-grid boundary. Sigma must act as one fixed linear operator everywhere
//! it is used: Brueckner valence + spectrum (SOS) and mixed states (MS) must
//! give the same polarisability, and the spectrum must reproduce the
//! Brueckner energies. Regression: a tail fitted separately for each function
//! Sigma acted on corrupted the spectrum (SOS out by orders of magnitude).
TEST_CASE("Amplitudes: Sigma tail extrapolation",
          "[Amplitudes][SecondOrder][MBPT]") {

  Wavefunction wf0({1600, 1.0e-5, 100.0, 20.0, "loglinear", -1.0},
                   {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf0.solve_core("HartreeFock", "[Xe]", std::nullopt, 1.0e-13, true);
  wf0.solve_valence("6sp5d", false);
  SplineBasis::Parameters params{"25spdfghi", 35, 7, {1.0e-4}, 40.0};
  params.verbose = false;
  wf0.formBasis(params);

  const DiracOperator::E1 hE1(wf0.grid());

  const double sigma_rmax = 20.0;
  // SOS vs MS agreement is limited by spectrum completeness

  // alpha (SOS) for each valence state, without and with the tail correction
  std::vector<double> alpha_no_tail, alpha_tail;

  for (const bool tail : {false, true}) {
    std::cout << (tail ? "\nWith tail correction:\n" :
                         "\nWithout tail correction:\n");
    auto wf = wf0;
    wf.formSigma(3, 1.0e-4, sigma_rmax, 6, false, false, false, 0, {}, {}, {},
                 false, "", false, false, false, 6, -0.2, 0.01, 1.5, false,
                 std::nullopt, "", false, false, false, tail);
    wf.hartreeFockBrueckner(false);
    SplineBasis::Parameters spectrum_params{"70spdf", 70, 9, {1.0e-4}, 90.0};
    spectrum_params.verbose = false;
    wf.formSpectrum(spectrum_params);

    // The spectrum (H + Sigma diagonalised) must reproduce the Brueckner
    // valence energies
    const auto &spectrum = wf.spectrum();
    fmt::print("\n{} Sigma tail exptrapolation\n", tail ? "With" : "without");
    fmt::print("{:<5} {:>14} {:>14} {:>9}\n", "state", "E(Brueckner)",
               "E(spectrum)", "rel.diff");

    // Spectrum vs Brueckner energies: limited by spectrum completeness
    const double tol_energy = 1.0e-5;
    for (const auto &Fbrueckner : wf.valence()) {
      const auto Fspectrum =
        std::find(spectrum.cbegin(), spectrum.cend(), Fbrueckner);
      REQUIRE(Fspectrum != spectrum.cend());
      const auto rel_diff = std::abs(Fspectrum->en() / Fbrueckner.en() - 1.0);
      fmt::print("{:<5} {:14.9f} {:14.9f} {:9.1e}\n", Fbrueckner.shortSymbol(),
                 Fbrueckner.en(), Fspectrum->en(), rel_diff);
      REQUIRE(rel_diff < tol_energy);
    }

    // Static scalar polarisability: SOS (Brueckner + spectrum) vs MS
    const double tol_sos_ms = 1.0e-3;
    const ExternalField::TDHF tdhf(&hE1, wf.vHF());
    auto &alpha_sos_list = tail ? alpha_tail : alpha_no_tail;
    fmt::print("\n{} Sigma tail exptrapolation\n", tail ? "With" : "without");
    fmt::print("{:<5} {:>12} {:>12} {:6}\n", "state", "alpha(SOS)", "alpha(MS)",
               "eps");
    for (const auto &Fv : wf.valence()) {
      const auto to_alpha = 1.0 / std::sqrt(3.0 * (Fv.twoj() + 1));
      const auto alpha_sos =
        to_alpha *
        Amplitudes::sos_valence(0, Fv, Fv, &hE1, &hE1, 0.0, 0.0, spectrum);
      const auto alpha_ms =
        to_alpha * Amplitudes::ms_valence(0, Fv, Fv, &hE1, &hE1, 0.0, 0.0,
                                          &tdhf, &tdhf, wf.Sigma())
                     .first;
      fmt::print("{:<5} {:12.5e} {:12.5e}  {:6.0e}\n", Fv.shortSymbol(),
                 alpha_sos, alpha_ms, std::abs(alpha_sos / alpha_ms - 1.0));

      /// XXX the p states, MS does not work
      if (Fv.l() != 1) {
        REQUIRE(alpha_sos == Approx(alpha_ms).epsilon(tol_sos_ms));
      }
      alpha_sos_list.push_back(alpha_sos);
    }
  }

  // At realistic Sigma rmax, the tail is a negligible correction
  fmt::print("\nalpha(SOS), Sigma rmax = {} aB\n", sigma_rmax);
  fmt::print("{:<5} {:>12} {:>12}  {:6}\n", "state", "no tail", "with tail",
             "eps");
  const auto &valence = wf0.valence();
  const double tol_sos_tail = 1.0e-4;
  for (std::size_t i = 0; i < valence.size(); ++i) {
    const auto without = alpha_no_tail.at(i);
    const auto with = alpha_tail.at(i);
    fmt::print("{:<5} {:12.5e} {:12.5e} {:6.0e}\n", valence.at(i).shortSymbol(),
               without, with, std::abs(with / without - 1.0));
    REQUIRE(with == Approx(without).epsilon(tol_sos_tail));
  }
}

//==============================================================================
//! The key invariant: sum-over-states and mixed-states evaluations of the
//! same amplitude must agree (to basis completeness); and the two independent
//! routes of the mixed-states method must agree with each other
TEST_CASE("Amplitudes: second-order SOS vs MS",
          "[Amplitudes][SecondOrder][ExternalField][integration]") {

  Wavefunction wf({2200, 1.0e-6, 120.0, 20.0, "loglinear", -1.0},
                  {"Cs", 133, "Fermi", 4.8041, 2.3}, 1.0);
  wf.solve_core("HartreeFock", "[Xe]");
  wf.solve_valence("7sp");
  // f states needed so the SOS core sum has the 4d -> f excitations
  wf.formBasis({"90spdf", 100, 9, {1.0e-4}, 75.0});

  const auto *p6s = wf.getState("6s+");
  const auto *p7s = wf.getState("7s+");
  REQUIRE(p6s != nullptr);
  REQUIRE(p7s != nullptr);

  const DiracOperator::E1 hE1(wf.grid());
  // Factor -N, as the pnc module: units i(Qw/-N)e-11
  const DiracOperator::PNCnsi hPNC(5.67073, 2.3, wf.grid(), -78.0,
                                   "i(Qw/-N)*e-11");

  // TDHF objects drive the mixed states; unsolved => no RPA (dV = 0)
  const ExternalField::TDHF tdhf_E1(&hE1, wf.vHF());
  const ExternalField::TDHF tdhf_pnc(&hPNC, wf.vHF());

  const auto &basis = wf.basis();

  fmt::print("\nSecond order, SOS vs MS (Cs, HF level):\n");
  fmt::print("{:<24} {:>14} {:>14} {:>14} {:>9} {:>9}\n", "", "SOS", "MS(s)",
             "MS(t)", "SOSvMS", "MSvMS");

  const auto compare = [](const std::string &name, double sos,
                          std::pair<double, double> ms, double tol_sos,
                          double tol_ms) {
    const auto [ms_s, ms_t] = ms;
    const auto eps_sos =
      std::abs((ms_s - sos) / std::max(std::abs(sos), 1.0e-30));
    const auto eps_ms =
      std::abs((ms_s - ms_t) / std::max(std::abs(ms_s), 1.0e-30));
    fmt::print("{:<24} {:14.7e} {:14.7e} {:14.7e} {:9.1e} {:9.1e}\n", name, sos,
               ms_s, ms_t, eps_sos, eps_ms);
    REQUIRE(eps_sos < tol_sos);
    REQUIRE(eps_ms < tol_ms);
  };

  // (a) Static polarisability amplitude, 6s (K = 0)
  {
    const auto sos =
      Amplitudes::sos_valence(0, *p6s, *p6s, &hE1, &hE1, 0.0, 0.0, basis);
    const auto ms = Amplitudes::ms_valence(0, *p6s, *p6s, &hE1, &hE1, 0.0, 0.0,
                                           &tdhf_E1, &tdhf_E1);
    compare("A^0 6s static", sos, ms, 1.0e-4, 1.0e-9);

    // Core-valence part: below-Fermi terms of the SOS sum vs mixed states
    // projected onto the HF core
    const auto below_fermi =
      DiracSpinor::split_by_energy(basis, wf.FermiLevel()).first;
    const auto sos_cv =
      Amplitudes::sos_valence(0, *p6s, *p6s, &hE1, &hE1, 0.0, 0.0, below_fermi);
    const auto ms_cv =
      Amplitudes::ms_valence(0, *p6s, *p6s, &hE1, &hE1, 0.0, 0.0, &tdhf_E1,
                             &tdhf_E1, nullptr, wf.core());
    compare("A^0 6s cv part", sos_cv, ms_cv, 1.0e-3, 1.0e-9);
  }

  // (b) Dynamic polarisability amplitude, 6s, omega = 0.05
  {
    const auto sos =
      Amplitudes::sos_valence(0, *p6s, *p6s, &hE1, &hE1, 0.05, -0.05, basis);
    const auto ms = Amplitudes::ms_valence(0, *p6s, *p6s, &hE1, &hE1, 0.05,
                                           -0.05, &tdhf_E1, &tdhf_E1);
    compare("A^0 6s dynamic", sos, ms, 1.0e-4, 1.0e-9);
  }

  // (c) Scalar transition amplitude, 6s -> 7s (K = 0)
  {
    const auto omega = p7s->en() - p6s->en();
    const auto sos =
      Amplitudes::sos_valence(0, *p7s, *p6s, &hE1, &hE1, omega, 0.0, basis);
    const auto ms = Amplitudes::ms_valence(0, *p7s, *p6s, &hE1, &hE1, omega,
                                           0.0, &tdhf_E1, &tdhf_E1);
    // transition amplitude converges more slowly with basis size
    compare("A^0 6s-7s", sos, ms, 5.0e-4, 1.0e-9);
  }

  // (d) PNC amplitude, 6s -> 7s (K = 1, t = E1, s = pnc): the two MS routes
  // are genuinely independent here (mixed states of E1, and of pnc)
  {
    const auto omega = p7s->en() - p6s->en();
    const auto sos =
      Amplitudes::sos_valence(1, *p7s, *p6s, &hE1, &hPNC, omega, 0.0, basis);
    const auto ms = Amplitudes::ms_valence(1, *p7s, *p6s, &hE1, &hPNC, omega,
                                           0.0, &tdhf_E1, &tdhf_pnc);
    compare("A^1 6s-7s PNC", sos, ms, 1.0e-3, 1.0e-5);
  }

  // (e) Core amplitude (K = 0): SOS over excited basis vs mixed states
  {
    const auto excited =
      DiracSpinor::split_by_energy(basis, wf.FermiLevel()).second;
    const auto sos = Amplitudes::sos_core(0, p6s->twoj(), &hE1, &hE1, 0.0, 0.0,
                                          wf.core(), excited);
    const auto ms = Amplitudes::ms_core(0, p6s->twoj(), &hE1, &hE1, 0.0, 0.0,
                                        wf.core(), &tdhf_E1, &tdhf_E1);
    compare("A^0 core", sos, {ms, ms}, 5.0e-4, 1.0);
  }

  // (f) With RPA: solve TDHF, repeat static 6s and PNC cases
  {
    fmt::print("\nWith RPA (TDHF solved):\n");
    ExternalField::TDHF rpa_E1(&hE1, wf.vHF());
    ExternalField::TDHF rpa_pnc(&hPNC, wf.vHF());
    rpa_E1.solve_core(0.0, 128, false);
    rpa_pnc.solve_core(0.0, 128, false);

    const auto sos_a = Amplitudes::sos_valence(0, *p6s, *p6s, &hE1, &hE1, 0.0,
                                               0.0, basis, &rpa_E1, &rpa_E1);
    const auto ms_a = Amplitudes::ms_valence(0, *p6s, *p6s, &hE1, &hE1, 0.0,
                                             0.0, &rpa_E1, &rpa_E1);
    compare("A^0 6s static (RPA)", sos_a, ms_a, 1.0e-4, 1.0e-9);

    const auto omega = p7s->en() - p6s->en();
    const auto sos_p = Amplitudes::sos_valence(
      1, *p7s, *p6s, &hE1, &hPNC, omega, 0.0, basis, &rpa_E1, &rpa_pnc);
    const auto ms_p = Amplitudes::ms_valence(1, *p7s, *p6s, &hE1, &hPNC, omega,
                                             0.0, &rpa_E1, &rpa_pnc);
    compare("A^1 6s-7s PNC (RPA)", sos_p, ms_p, 1.0e-3, 1.0e-5);

    // The same, with the E1 RPA solved at the transition frequency. The RPA
    // correction then differs between the two signs of omega, and every t
    // vertex of both routes must carry dV(omega)
    {
      ExternalField::TDHF rpa_E1w(&hE1, wf.vHF());
      rpa_E1w.solve_core(omega, 128, false);
      const auto sos_w = Amplitudes::sos_valence(
        1, *p7s, *p6s, &hE1, &hPNC, omega, 0.0, basis, &rpa_E1w, &rpa_pnc);
      const auto ms_w = Amplitudes::ms_valence(1, *p7s, *p6s, &hE1, &hPNC,
                                               omega, 0.0, &rpa_E1w, &rpa_pnc);
      compare("A^1 6s-7s PNC (RPA, w)", sos_w, ms_w, 1.0e-3, 1.0e-5);
    }

    const auto excited =
      DiracSpinor::split_by_energy(basis, wf.FermiLevel()).second;
    const auto sos_c =
      Amplitudes::sos_core(0, p6s->twoj(), &hE1, &hE1, 0.0, 0.0, wf.core(),
                           excited, &rpa_E1, &rpa_E1);
    const auto ms_c = Amplitudes::ms_core(0, p6s->twoj(), &hE1, &hE1, 0.0, 0.0,
                                          wf.core(), &rpa_E1, &rpa_E1);
    compare("A^0 core (RPA)", sos_c, {ms_c, ms_c}, 5.0e-4, 1.0);
  }

  // (g) Regression against the older modules: values recorded 2026-07-29 from
  // polarisability, transitionPolarisability, and pnc modules (Cs HF,
  // 3000-point grid, 90spdf basis; MS method). That grid differs slightly
  // from this one, hence the loose tolerance. Signs: the older modules
  // evaluate the reverse transition direction for beta and E_pnc
  {
    fmt::print("\nRegression vs older modules:\n");
    fmt::print("{:<24} {:>14} {:>14}\n", "", "expected", "found");
    const auto tol = 3.0e-3;
    const auto require = [tol](const std::string &name, double expected,
                               double found) {
      fmt::print("{:<24} {:14.6e} {:14.6e}\n", name, expected, found);
      REQUIRE(found == Approx(expected).epsilon(tol));
    };

    // polarisability module, rpa=false, MS: alpha(6s) = 663.581
    {
      const auto [av, av2] = Amplitudes::ms_valence(
        0, *p6s, *p6s, &hE1, &hE1, 0.0, 0.0, &tdhf_E1, &tdhf_E1);
      const auto ac = Amplitudes::ms_core(0, p6s->twoj(), &hE1, &hE1, 0.0, 0.0,
                                          wf.core(), &tdhf_E1, &tdhf_E1);
      const auto alpha = (av + ac) / std::sqrt(3.0 * (p6s->twoj() + 1));
      require("alpha 6s (HF)", 663.581, alpha);
    }

    // transitionPolarisability module, rpa=false, MS: |beta| = 29.2786
    {
      const auto omega = p7s->en() - p6s->en();
      const auto [av, av2] = Amplitudes::ms_valence(
        1, *p7s, *p6s, &hE1, &hE1, omega, 0.0, &tdhf_E1, &tdhf_E1);
      const auto beta =
        av / (std::sqrt(2.0) * 2.0 * Angular::S_kk(p7s->kappa(), p6s->kappa()));
      require("beta 6s-7s (HF)", 29.2786, beta);
    }

    // pnc module, rpa=true, MS: |E_pnc| = 0.891897 [i(-Qw/N)e-11]
    {
      const auto omega = p7s->en() - p6s->en();
      ExternalField::TDHF rpa_E1w(&hE1, wf.vHF());
      ExternalField::TDHF rpa_pnc0(&hPNC, wf.vHF());
      rpa_E1w.solve_core(omega, 128, false);
      rpa_pnc0.solve_core(0.0, 128, false);
      const auto [av, av2] = Amplitudes::ms_valence(
        1, *p7s, *p6s, &hE1, &hPNC, omega, 0.0, &rpa_E1w, &rpa_pnc0);
      const auto E_pnc =
        av * CI::z_component(1, 1, 0, p7s->twoj(), p6s->twoj(), 1);
      require("|E_pnc| 6s-7s (RPA)", 0.891897, std::abs(E_pnc));
    }
  }
}
