#include "Angular/Wigner369j.hpp"
#include "CI_Integrals.hpp"
#include "ConfigurationInteraction.hpp"
#include "Coulomb/include.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/InputBlock.hpp"
#include "LinAlg/include.hpp"
#include "MBPT/Sigma2.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "qip/Random.hpp"
#include <algorithm>
#include <array>
#include <cstdio>
#include <vector>

//==============================================================================
TEST_CASE("CI: Configuration Interaction unit tests", "[CI][unit]") {

  std::cout << "CI, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"He", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5);
  wf.formBasis(SplineBasis::Parameters("6spd", 20, 6, {1.0e-2}, 20.0));

  std::string qk_filename = "deleteme_" + qip::random_string(3) + ".qk.abf";
  std::string ci_filename = "deleteme_" + qip::random_string(3) + ".ci.abf";

  const IO::InputBlock input{"CI", "ci_basis = 6spd;"
                                   "n_min_core = 1;"
                                   "J + = 0, 1;"
                                   "J - = 0, 1;"
                                   "qk_file = " +
                                     qk_filename +
                                     ";"
                                     "num_solutions = 2;"
                                     "ci_file = " +
                                     ci_filename + ";"};

  // Initial CI calculation (writes to file):
  const auto CIWFs = CI::configuration_interaction(input, wf).levels;

  // Expected values: (not accurate, just simple unit test)
  std::vector J{0, 1, 0, 1};
  std::vector parity{1, 1, -1, -1};
  std::vector energy{std::vector{-2.84834207, -2.13740814},
                     {-2.17164001, -2.0657259},
                     {-2.13100411, -2.05532175},
                     {-2.13100096, -2.12092003}};
  std::vector gj{
    std::vector{0.0, 0.0}, {1.9999, 1.9999}, {0.0, 0.0}, {1.5, 1.0}};
  // Terms: {1S0, 1S0}, {3S1, 3S1}, {3P0, 3P0}, {3P1, 1P1}
  std::vector L2{std::vector{0.0, 0.0}, {0.0, 0.0}, {2.0, 2.0}, {2.0, 2.0}};
  std::vector S2{std::vector{0.0, 0.0}, {2.0, 2.0}, {2.0, 2.0}, {2.0, 0.0}};

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &ci_wf = CIWFs.at(i);
    REQUIRE(ci_wf.twoJ() == 2 * J[i]);
    REQUIRE(ci_wf.parity() == parity[i]);
    for (std::size_t j = 0ul; j < ci_wf.num_solutions(); ++j) {
      REQUIRE(ci_wf.energy(j) == Approx(energy[i][j]).epsilon(1.0e-2));
      REQUIRE(ci_wf.info(j).gJ == Approx(gj[i][j]).epsilon(1.0e-2));
      REQUIRE(ci_wf.info(j).L2 == Approx(L2[i][j]).margin(0.05));
      REQUIRE(ci_wf.info(j).S2 == Approx(S2[i][j]).margin(0.05));
    }
  }

  //-----------------------------------------------------------------------
  std::cout << "\nRe-run CI, reading in from file:\n";

  // Re-"run" identical CI calculation
  // This time: should read in prev solution from disk
  // Read from ci file; compare energies, gJ, and coefficients
  auto input2 = input;
  input2.add("read_only = true;");
  const auto solution_2 = CI::configuration_interaction(input2, wf);
  const auto CIWFs_2 = solution_2.levels;

  // read_only must still read in the integrals: the modules need them to
  // re-construct the CI Hamiltonian of any other J/parity
  REQUIRE(solution_2.integrals.availableQ());
  REQUIRE(solution_2.integrals.qk.count() ==
          CI::configuration_interaction(input, wf).integrals.qk.count());

  REQUIRE(CIWFs.size() == CIWFs_2.size());

  for (std::size_t i = 0; i < CIWFs.size(); ++i) {
    const auto &w = CIWFs.at(i);
    const auto &r = CIWFs_2.at(i);
    REQUIRE(w.twoJ() == r.twoJ());
    REQUIRE(w.parity() == r.parity());
    REQUIRE(w.num_solutions() == r.num_solutions());
    REQUIRE(w.CSFs().size() == r.CSFs().size());
    for (std::size_t s = 0; s < w.num_solutions(); ++s) {
      REQUIRE(w.energy(s) == Approx(r.energy(s)));
      REQUIRE(w.info(s).gJ == Approx(r.info(s).gJ));
      const auto wc = w.coefs(s);
      const auto rc = r.coefs(s);
      REQUIRE(wc.size() == rc.size());
      REQUIRE(wc.size() == w.CSFs().size());
      for (std::size_t k = 0; k < wc.size(); ++k) {
        REQUIRE(wc[k] == Approx(rc[k]));
      }
    }
  }

  //-----------------------------------------------------------------------
  std::cout << "\nRe-run CI, reading in from file: but ask for more states\n";
  // Now, request 3 solutions:
  // ci file has 2, so must re-solve and overwrite
  const IO::InputBlock input3{"CI", "ci_basis = 6spd;"
                                    "n_min_core = 1;"
                                    "J + = 0, 1;"
                                    "J - = 0, 1;"
                                    "qk_file = " +
                                      qk_filename +
                                      ";"
                                      "num_solutions = 3;"
                                      "print_details = false;"
                                      "ci_file = " +
                                      ci_filename + ";"};
  // don't read in
  const IO::InputBlock input4{"CI", "ci_basis = 6spd;"
                                    "n_min_core = 1;"
                                    "J + = 0, 1;"
                                    "J - = 0, 1;"
                                    "qk_file = false;"
                                    "num_solutions = 3;"
                                    "print_details = false;"
                                    "ci_file = " +
                                      ci_filename + ";"};

  // Reads from file (2 stored, 3 required, re-calculates)
  const auto CIWFs_3a = CI::configuration_interaction(input3, wf).levels;
  // calculate 3 from scratch
  const auto CIWFs_3b = CI::configuration_interaction(input3, wf).levels;

  REQUIRE(CIWFs_3a.size() == CIWFs.size());
  REQUIRE(CIWFs_3b.size() == CIWFs.size());
  for (std::size_t i = 0; i < CIWFs_3a.size(); ++i) {
    const auto &a = CIWFs_3a.at(i);
    const auto &b = CIWFs_3b.at(i);
    REQUIRE(a.num_solutions() == 3);
    REQUIRE(b.num_solutions() == 3);
    // First 2 solutions match CIWFs
    for (std::size_t s = 0; s < 2; ++s) {
      REQUIRE(a.energy(s) == Approx(CIWFs.at(i).energy(s)));
      REQUIRE(a.info(s).gJ == Approx(CIWFs.at(i).info(s).gJ));
    }
    // All 3 identical between solve and read
    for (std::size_t s = 0; s < 3; ++s) {
      REQUIRE(a.energy(s) == Approx(b.energy(s)));
      REQUIRE(a.info(s).gJ == Approx(b.info(s).gJ));
      const auto ac = a.coefs(s);
      const auto bc = b.coefs(s);
      REQUIRE(ac.size() == bc.size());
      for (std::size_t k = 0; k < ac.size(); ++k) {
        REQUIRE(ac[k] == Approx(bc[k]));
      }
    }
  }

  //-----------------------------------------------------------------------
  // Sigma_2 extrapolation

  {
    const auto &ints = solution_2.integrals;

    // Table with S^k = Q^k for a small subset: ratio S/Q = 1 exactly
    Coulomb::LkTable Sk_test;
    const auto S_eq_Q = [&ints](int k, const DiracSpinor &a,
                                const DiracSpinor &b, const DiracSpinor &c,
                                const DiracSpinor &d) {
      return ints.qk.Q(k, a, b, c, d);
    };
    const auto Qk_SR = [](int k, const DiracSpinor &a, const DiracSpinor &b,
                          const DiracSpinor &c, const DiracSpinor &d) {
      return Coulomb::Qk_abcd_SR(k, a, b, c, d);
    };
    const auto sub_basis = CI::basis_subset(ints.ci_basis, "3sp");
    Sk_test.fill(sub_basis, S_eq_Q, Qk_SR, 8, false);

    // Average correction ratios: exactly 1 for every k with data
    const auto hk = MBPT::average_hk(Sk_test, ints.qk, sub_basis, 8);
    bool any_data = false;
    for (const auto &h : hk) {
      if (h != 0.0) {
        any_data = true;
        REQUIRE(h == Approx(1.0));
      }
    }
    REQUIRE(any_data);

    // Extrapolate with hk = 1 for all k: then S^k = Q^k everywhere, and
    // Sigma2_AB must equal the Coulomb part for every CSF pair.
    // nb: Sk_test only covers sub_basis, so all other pairs are extrapolated
    const std::vector<double> hk_1(10, 1.0);
    for (const auto &psi0 : CIWFs) {
      const auto twoJ = psi0.twoJ();
      const auto n_csf = std::min(psi0.CSFs().size(), std::size_t(6));
      for (std::size_t i = 0; i < n_csf; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
          const auto &A = psi0.CSF(i);
          const auto &B = psi0.CSF(j);
          const auto s2_x = CI::Sigma2_AB(A, B, twoJ, Sk_test, &ints.qk, hk_1);
          const auto [v, w] = A.states;
          const auto [x, y] = B.states;
          const auto coulomb = CI::CSF2_Coulomb(ints.qk, v, w, x, y, twoJ);
          REQUIRE(s2_x == Approx(coulomb).margin(1.0e-12));
        }
      }
    }
  }

  //-----------------------------------------------------------------------
  // Multiple (J, pi) blocks in one ci solutions file

  {
    const auto fname = "deleteme_" + qip::random_string(3) + ".ci.abf";
    const auto t_basis = CI::basis_subset(wf.basis(), "4sp");
    CI::PsiJPi psi_a{0, +1, t_basis};
    psi_a.set_solution(-1.0, LinAlg::Vector<double>(psi_a.CSFs().size()));
    REQUIRE(psi_a.read_write(fname, IO::FRW::write, std::cout));
    // second block is appended; the first must survive
    CI::PsiJPi psi_b{2, +1, t_basis};
    psi_b.set_solution(-0.5, LinAlg::Vector<double>(psi_b.CSFs().size()));
    REQUIRE(psi_b.read_write(fname, IO::FRW::write, std::cout));
    CI::PsiJPi read_a{0, +1, t_basis};
    REQUIRE(read_a.read_write(fname, IO::FRW::read, std::cout));
    REQUIRE(read_a.energy(0) == Approx(-1.0));
    CI::PsiJPi read_b{2, +1, t_basis};
    REQUIRE(read_b.read_write(fname, IO::FRW::read, std::cout));
    REQUIRE(read_b.energy(0) == Approx(-0.5));
    // a block that was never written:
    CI::PsiJPi read_c{0, -1, t_basis};
    REQUIRE(!read_c.read_write(fname, IO::FRW::read, std::cout));
  }

  //-----------------------------------------------------------------------
  // basic/Misc tests

  // Term(int two_J, int L, int two_S, int parity)

  REQUIRE(CI::Term_Symbol(2, 3, 2, +1) == "3^F_1");
  REQUIRE(CI::Term_Symbol(2, 3, 2, -1) == "3^F°_1");
  REQUIRE(CI::Term_Symbol(1, 3, 2, -1) == "3^F°_1/2");
  REQUIRE(CI::Term_Symbol(1, 2, 1, -1) == "2^D°_1/2");
  REQUIRE(CI::Term_Symbol(6, 0, 0, 1) == "1^S_3");
}

//==============================================================================
TEST_CASE("CI: derivative correction", "[CI][unit]") {

  std::cout << "CI, derivative (dSigma/dE) correction for Sigma_1\n";

  //-----------------------------------------------------------------------
  // corrected_Sigma: formula and guard

  // Sigma = 0: no correction
  REQUIRE(CI::corrected_Sigma(0.0, 0.1, 0.5) == 0.0);

  // dE = 0: Sigma unchanged (Sigma^2/Sigma = Sigma)
  REQUIRE(CI::corrected_Sigma(-0.1, 0.05, 0.0) == Approx(-0.1).margin(1.0e-14));

  // Exact resummed value: Sigma*Sigma/(Sigma - dE*dSigma)
  REQUIRE(CI::corrected_Sigma(-0.1, -0.02, -0.5) ==
          Approx(0.01 / -0.11).margin(1.0e-14));

  // Small dE: matches linear expansion, Sigma + dE*dSigma
  {
    const auto dE = 1.0e-4;
    const auto lin = -0.1 + dE * 0.05;
    REQUIRE(CI::corrected_Sigma(-0.1, 0.05, dE) == Approx(lin).margin(1.0e-8));
  }

  // Guard: corrected value exceeds |Sigma| (near divergence): uncorrected
  // -0.1 - (-0.5)*(0.02) = -0.09, so |Sig^2/denom| = 0.111 > 0.1
  REQUIRE(CI::corrected_Sigma(-0.1, 0.02, -0.5) ==
          Approx(-0.1).margin(1.0e-14));

  //-----------------------------------------------------------------------
  // corrected_Sk: formula and pole guard

  // Sk = 0: no shift
  REQUIRE(CI::corrected_Sk(0.0, 0.1, 0.5) == 0.0);

  // dE0 = 0: Sk unchanged (Sk^2/Sk = Sk)
  REQUIRE(CI::corrected_Sk(-0.1, 0.05, 0.0) == Approx(-0.1).margin(1.0e-14));

  // Exact resummed value: Sk*Sk/(Sk - dE0*dSk)
  REQUIRE(CI::corrected_Sk(-0.1, -0.02, -0.5) ==
          Approx(0.01 / -0.11).margin(1.0e-14));

  // Unlike corrected_Sigma, enhancement is allowed:
  // 0.1 - 0.5*0.05 = 0.075 (same sign as Sk): |corrected| = 0.133 > 0.1
  REQUIRE(CI::corrected_Sk(0.1, 0.05, 0.5) ==
          Approx(0.01 / 0.075).margin(1.0e-14));

  // Small dE0: matches linear expansion, Sk + dE0*dSk
  {
    const auto dE0 = 1.0e-4;
    const auto lin = -0.1 + dE0 * 0.05;
    REQUIRE(CI::corrected_Sk(-0.1, 0.05, dE0) == Approx(lin).margin(1.0e-8));
  }

  // Exact for a single denominator: Sk = N/D, dSk = -N/D^2
  // => corrected = N/(D + dE0)
  {
    const auto N = 0.35, D = -1.2, dE0 = 0.4;
    const auto Sk = N / D;
    const auto dSk = -N / (D * D);
    REQUIRE(CI::corrected_Sk(Sk, dSk, dE0) ==
            Approx(N / (D + dE0)).margin(1.0e-14));
  }

  // Pole guard: resummed denominator changes sign (dE0 past the pole at
  // Sk/dSk = 2.0): return unshifted
  REQUIRE(CI::corrected_Sk(0.1, 0.05, 3.0) == Approx(0.1).margin(1.0e-14));
  // Exactly at the pole:
  REQUIRE(CI::corrected_Sk(0.1, 0.05, 2.0) == Approx(0.1).margin(1.0e-14));

  //-----------------------------------------------------------------------
  // Sigma1Correction::delta

  const auto ia = static_cast<DiracSpinor::Index>(Angular::nk_to_index(1, -1));
  const auto ib = static_cast<DiracSpinor::Index>(Angular::nk_to_index(2, -1));

  CI::Sigma1Correction corr;
  corr.E0 = -2.0;
  corr.e_sigma[-1] = -1.0;
  corr.en[ia] = -1.0;
  corr.en[ib] = -0.5;
  corr.S1.add(ia, ia, -0.1);
  corr.S1.add(ia, ib, -0.05);
  corr.S1.add(ib, ia, -0.05);
  corr.S1.add(ib, ib, -0.08);
  corr.dS1.add(ia, ia, -0.02);
  corr.dS1.add(ia, ib, -0.02);
  corr.dS1.add(ib, ia, -0.02);
  corr.dS1.add(ib, ib, -0.02);

  REQUIRE(!corr.empty());

  // delta(a,a; spectator b): dE = E0 - e_sigma - en(b) = -2 + 1 + 0.5 = -0.5
  // corrected = 0.01/(-0.1 - 0.01) = -0.0909..., delta = +0.00909...
  REQUIRE(corr.delta_h1(ia, ia, ib) ==
          Approx(0.01 / -0.11 + 0.1).margin(1.0e-14));

  // Symmetric in (a, b):
  REQUIRE(corr.delta_h1(ia, ib, ia) == Approx(corr.delta_h1(ib, ia, ia)));

  // Missing orbital/spectator: no correction
  const auto ix = static_cast<DiracSpinor::Index>(Angular::nk_to_index(3, 1));
  REQUIRE(corr.delta_h1(ix, ix, ia) == 0.0);
  REQUIRE(corr.delta_h1(ia, ia, ix) == 0.0);

  // Zero derivative: no correction (any dE)
  CI::Sigma1Correction corr0 = corr;
  corr0.dS1 = {};
  corr0.dS1.add(ia, ia, 0.0);
  REQUIRE(corr0.delta_h1(ia, ia, ib) == Approx(0.0).margin(1.0e-14));

  // Default (empty) correction:
  REQUIRE(CI::Sigma1Correction{}.empty());

  //-----------------------------------------------------------------------
  // read_write round-trip

  const auto ds1_fname = "deleteme_" + qip::random_string(3) + ".ds1.abf";
  REQUIRE(corr.read_write(ds1_fname, IO::FRW::write));

  CI::Sigma1Correction corr2;
  corr2.E0 = corr.E0; // E0 is not stored (per J/pi; set at solve time)
  REQUIRE(corr2.read_write(ds1_fname, IO::FRW::read));

  REQUIRE(corr2.S1.getv(ia, ib) == corr.S1.getv(ia, ib));
  REQUIRE(corr2.dS1.getv(ib, ib) == corr.dS1.getv(ib, ib));
  REQUIRE(corr2.e_sigma == corr.e_sigma);
  REQUIRE(corr2.en == corr.en);
  REQUIRE(corr2.delta_h1(ia, ia, ib) == corr.delta_h1(ia, ia, ib));

  // Missing file: read fails
  REQUIRE(!CI::Sigma1Correction{}.read_write("deleteme_no_such_file.abf",
                                             IO::FRW::read));
}

//==============================================================================
// Tests the Sigma_1 correction tables built from a correlation potential
// holding dSigma/dE, against the second-order (qk-route) tables
TEST_CASE("CI: dSdE from correlation potential", "[CI][integration]") {

  std::cout << "CI, dSigma/dE correction tables from correlation potential\n";

  Wavefunction wf({800, 1.0e-5, 40.0, 10.0, "loglinear"},
                  {"Na", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-7);
  wf.solve_valence("3sp");
  // wf.formBasis(SplineBasis::Parameters("12spd", 25, 7, {1.1e-3, 3.2e-2}, 30.0));
  wf.formBasis(SplineBasis::Parameters("12spd", 25, 7, {1.0e-3}, 30.0));

  // Goldstone correlation potential, with dSigma/dE (derivative = true);
  // written to file for the round-trip check below
  const auto sig_fname = "deleteme_" + qip::random_string(3) + ".s2.abf";
  wf.formSigma(1, 1.0e-4, 30.0, 2, false, false, false, 0, {}, {}, {}, true,
               sig_fname, false, false, false, 6, -0.2, 0.01, 1.5, false,
               std::nullopt, "", true);
  REQUIRE(wf.Sigma() != nullptr);
  REQUIRE(wf.Sigma()->has_derivative());

  // Tables from the correlation potential:
  const auto ci_basis =
    CI::basis_subset(wf.basis(), "5sp", wf.coreConfiguration());
  const auto corr_CP = CI::calculate_dSdE_correction(ci_basis, *wf.Sigma());
  REQUIRE(!corr_CP.empty());

  // Same tables via the second-order qk route:
  Coulomb::QkTable qk;
  const auto yk = Coulomb::YkTable(wf.basis());
  qk.fill(wf.basis(), yk, 8, false);
  const auto [core_s1, excited_s1] =
    DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel(), 1);
  const auto corr_G =
    CI::calculate_dSdE_correction(ci_basis, core_s1, excited_s1, qk);

  // Same physics (both 2nd-order Goldstone), different machinery (sub-grid
  // matrix application vs direct basis sums): agree to sub-grid accuracy
  std::cout << "     S1(CP)      S1(qk)      dS1(CP)     dS1(qk)\n";
  for (const auto &v : ci_basis) {
    const auto S_cp = corr_CP.S1.getv(v, v);
    const auto S_g = corr_G.S1.getv(v, v);
    const auto dS_cp = corr_CP.dS1.getv(v, v);
    const auto dS_g = corr_G.dS1.getv(v, v);
    fmt::print("{:4s} {:+.8f} {:+.8f} {:+.7f} {:+.7f}\n", v.shortSymbol(), S_cp,
               S_g, dS_cp, dS_g);
    REQUIRE(S_cp == Approx(S_g).epsilon(0.05).margin(1.0e-6));
    REQUIRE(dS_cp == Approx(dS_g).epsilon(0.05).margin(1.0e-5));
  }

  // e_sigma: energy Sigma was formed at (HF valence) vs first ci-basis state
  // (basis approximates HF): close, not identical
  for (const auto &[kappa, e_cp] : corr_CP.e_sigma) {
    REQUIRE(e_cp == Approx(corr_G.e_sigma.at(kappa)).epsilon(1.0e-3));
  }

  //-----------------------------------------------------------------------
  // Round-trip: dSigma/dE read back from the sigma file

  const MBPT::CorrelationPotential Sig2(sig_fname, wf.vHF(), wf.basis(), 1.0e-4,
                                        30.0, 2, 1,
                                        MBPT::SigmaMethod::Goldstone);
  REQUIRE(Sig2.has_derivative());
  for (const auto &v : ci_basis) {
    REQUIRE(v * Sig2.dSigmaFv(v) ==
            Approx(v * wf.Sigma()->dSigmaFv(v)).margin(1.0e-12));
  }

  //-----------------------------------------------------------------------
  // fk stored alongside Sigma: storage, average_fk, and file round-trip

  {
    const auto fk_test = std::vector{0.7, 0.8, 0.9};
    const auto fk_fname = "deleteme_" + qip::random_string(3) + ".s2.abf";
    MBPT::CorrelationPotential Sig3(fk_fname, wf.vHF(), wf.basis(), 1.0e-4,
                                    30.0, 4, 1, MBPT::SigmaMethod::Goldstone,
                                    false, false, 0, {}, false, fk_test);
    Sig3.formSigma(-1, -0.19, 3);
    Sig3.formSigma(1, -0.11, 3);
    REQUIRE(Sig3.get(-1) != nullptr);
    REQUIRE(Sig3.get(-1)->fk == fk_test);
    // both stored sets identical: average is exact
    REQUIRE(Sig3.average_fk() == fk_test);
    // l_max = 0: only the s set
    REQUIRE(Sig3.average_fk(0) == fk_test);
    Sig3.write(fk_fname);

    const MBPT::CorrelationPotential Sig4(fk_fname, wf.vHF(), wf.basis(),
                                          1.0e-4, 30.0, 4, 1,
                                          MBPT::SigmaMethod::Goldstone);
    REQUIRE(Sig4.get(-1) != nullptr);
    REQUIRE(Sig4.get(-1)->fk == fk_test);
    REQUIRE(Sig4.average_fk() == fk_test);
  }

  // No stored fk: average is empty
  REQUIRE(wf.Sigma()->average_fk().empty());
}

//==============================================================================
// Three-electron CI Hamiltonian (CSF3), against independent results:
//  1. Closed s1/2 pair plus one electron: the block is the energy of the pair,
//     plus h1, plus the Hartree-Fock potential of the pair (direct and
//     exchange), which is computed independently by the HF routines.
//  2. A third electron that does not interact with the other two: the
//     spectrum is the two-electron CI spectrum plus its orbital energy. Covers
//     Coulomb, Breit, Sigma_2 (incl. extrapolation and the Brillouin-Wigner
//     pair-energy shift), and the dSigma/dE correction.
//  3. Centre-of-gravity sum rule: sum_J [J] Tr(H_J) over all CSFs of one
//     configuration equals the trace over its determinants, which is formed
//     from the two-electron (pair) traces.
TEST_CASE("CI: three-electron Hamiltonian", "[CI][3VCI][unit]") {

  std::cout << "CI three-electron Hamiltonian\n";

  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[He]", std::nullopt, 1.0e-5, false);
  auto b_params = SplineBasis::Parameters("5spd", 20, 6, {1.0e-2}, 20.0);
  b_params.verbose = false;
  wf.formBasis(b_params);

  // Closed 1s pair, and excited basis states
  const auto &a = wf.core().front();
  REQUIRE(a.n() == 1);
  std::vector<DiracSpinor> orbitals{a};
  for (const auto &Fn : wf.basis()) {
    if (Fn.n() > 1 && Fn.l() <= 2) {
      orbitals.push_back(Fn);
    }
  }

  Coulomb::QkTable qk;
  qk.fill(orbitals, Coulomb::YkTable(orbitals), -1, false);

  // One-body table: orbital energies, plus (arbitrary) off-diagonal elements
  Coulomb::meTable<double> h1;
  for (const auto &v : orbitals) {
    for (const auto &w : orbitals) {
      if (v.kappa() == w.kappa()) {
        h1.add(v, w, v == w ? v.en() : 0.01 / (v.n() + w.n()));
      }
    }
  }

  //----------------------------------------------------------------------------
  // 1. Closed pair a^2 plus one electron z: J = j_z, one CSF per z
  double worst_closed = 0.0;
  {
    const auto ia = a.nk_index();
    const auto E_pair =
      2.0 * h1.getv(ia, ia) + CI::CSF2_Coulomb(qk, ia, ia, ia, ia, 0);

    // Direct potential of the closed pair
    auto v_direct = Coulomb::yk_ab(0, a, a);
    for (auto &v : v_direct) {
      v *= a.twojp1();
    }

    for (const auto &z : orbitals) {
      if (z == a)
        continue;
      const auto csf_z = CI::form_CSF3s({ia, ia, z.nk_index()}, z.twoj());
      REQUIRE(csf_z.size() == 1);
      for (const auto &zp : orbitals) {
        if (zp == a || zp.kappa() != z.kappa())
          continue;
        const auto csf_zp = CI::form_CSF3s({ia, ia, zp.nk_index()}, zp.twoj());
        const auto H = CI::Hab(csf_z[0], csf_zp[0], z.twoj(), h1, qk);
        const auto V_HF =
          z * (v_direct * zp) + z * HF::vexFa(zp, std::vector{a});
        const auto expected = (z == zp ? E_pair : 0.0) + h1.getv(z, zp) + V_HF;
        worst_closed = std::max(worst_closed, std::abs(H - expected));
      }
    }
  }

  //----------------------------------------------------------------------------
  // 2. Third electron z (3d3/2) does not interact with the pair space S:
  // no integrals involving z, and no one-body coupling to S
  double worst_spectator = 0.0;
  {
    std::vector<DiracSpinor> S;
    std::vector<DiracSpinor> S_small; // Sk tabulated here only
    for (const auto &Fn : wf.basis()) {
      if (Fn.n() > 1 && Fn.l() <= 1) {
        S.push_back(Fn);
        if (Fn.n() <= 3) {
          S_small.push_back(Fn);
        }
      }
    }
    const auto &z =
      *std::find_if(wf.basis().cbegin(), wf.basis().cend(), [](const auto &Fn) {
        return Fn.n() == 3 && Fn.kappa() == 2;
      });

    Coulomb::QkTable qk_S;
    qk_S.fill(S, Coulomb::YkTable(S), -1, false);

    // Breit and Sigma_2 tables: arbitrary multiples of Q^k
    const auto all_k = [](int, const DiracSpinor &, const DiracSpinor &,
                          const DiracSpinor &,
                          const DiracSpinor &) { return true; };
    Coulomb::WkTable Bk;
    Bk.fill(
      S,
      [&qk_S](int k, const DiracSpinor &v, const DiracSpinor &w,
              const DiracSpinor &x,
              const DiracSpinor &y) { return 0.3 * qk_S.Q(k, v, w, x, y); },
      all_k, -1, false);
    Coulomb::LkTable Sk, dSk;
    Sk.fill(
      S_small,
      [&qk_S](int k, const DiracSpinor &v, const DiracSpinor &w,
              const DiracSpinor &x,
              const DiracSpinor &y) { return 0.2 * qk_S.Q(k, v, w, x, y); },
      all_k, -1, false);
    dSk.fill(
      S_small,
      [&qk_S](int k, const DiracSpinor &v, const DiracSpinor &w,
              const DiracSpinor &x,
              const DiracSpinor &y) { return 0.05 * qk_S.Q(k, v, w, x, y); },
      all_k, -1, false);
    // Sigma_2 extrapolation ratios, beyond S_small
    const std::vector<double> hk{0.1, 0.12, 0.14, 0.16};

    // dSigma/dE correction: no entries for z
    CI::Sigma1Correction s1c_2e;
    for (const auto &v : S) {
      s1c_2e.en[v.nk_index()] = v.en();
      s1c_2e.e_sigma.insert({v.kappa(), v.en()});
      for (const auto &w : S) {
        if (v.kappa() == w.kappa()) {
          s1c_2e.S1.add(v, w, -0.02 / (v.n() * w.n()));
          s1c_2e.dS1.add(v, w, -0.005 / (v.n() + w.n()));
        }
      }
    }
    s1c_2e.en[z.nk_index()] = z.en();
    s1c_2e.E0 = -0.6;
    const auto dE0_2e = 0.03;

    // Three-electron: the total energy includes the third electron
    auto s1c_3e = s1c_2e;
    s1c_3e.E0 = s1c_2e.E0 + z.en();
    const auto dE0_3e = dE0_2e + z.en();

    for (const int parity : {1, -1}) {
      for (const int twoJ : {1, 3, 5}) {

        std::vector<CI::CSF3> csfs_3e;
        for (std::size_t i = 0; i < S.size(); ++i) {
          for (std::size_t j = i; j < S.size(); ++j) {
            for (auto &csf : CI::form_CSF3s(
                   {S[i].nk_index(), S[j].nk_index(), z.nk_index()}, twoJ)) {
              if (csf.parity() == parity) {
                csfs_3e.push_back(std::move(csf));
              }
            }
          }
        }
        const auto H_3e = CI::construct_Hci(csfs_3e, twoJ, h1, qk_S, &Bk, &Sk,
                                            &s1c_3e, hk, &dSk, dE0_3e);
        const auto E_3e = LinAlg::symmhEigensystem(H_3e).first;

        // Two-electron spectra of each J2 that couples with j_z to J
        std::vector<double> expected;
        for (int twoJ2 = std::abs(twoJ - z.twoj()); twoJ2 <= twoJ + z.twoj();
             twoJ2 += 2) {
          const CI::PsiJPi psi_2e(twoJ2, parity * z.parity(), S);
          if (psi_2e.CSFs().empty())
            continue;
          const auto H_2e = CI::construct_Hci(psi_2e, h1, qk_S, &Bk, &Sk,
                                              &s1c_2e, hk, &dSk, dE0_2e);
          for (const auto e : LinAlg::symmhEigensystem(H_2e).first) {
            expected.push_back(e + z.en());
          }
        }
        std::sort(expected.begin(), expected.end());

        REQUIRE(E_3e.size() == expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
          worst_spectator =
            std::max(worst_spectator, std::abs(E_3e[i] - expected[i]));
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  // 3. Centre-of-gravity sum rule, for each type of configuration
  double worst_trace = 0.0;
  {
    const auto index = [&wf](int n, int kappa) {
      return std::find_if(wf.basis().cbegin(), wf.basis().cend(),
                          [n, kappa](const auto &Fn) {
                            return Fn.n() == n && Fn.kappa() == kappa;
                          })
        ->nk_index();
    };

    // Trace of g over the determinants of orbitals {x, y}
    const auto pair_trace = [&qk](DiracSpinor::Index x, DiracSpinor::Index y) {
      const auto tjx = Angular::nkindex_to_twoj(x);
      const auto tjy = Angular::nkindex_to_twoj(y);
      double trace = 0.0;
      for (int twoJ2 = std::abs(tjx - tjy); twoJ2 <= tjx + tjy; twoJ2 += 2) {
        if (x == y && twoJ2 % 4 != 0)
          continue;
        trace += (twoJ2 + 1) * CI::CSF2_Coulomb(qk, x, y, x, y, twoJ2);
      }
      return trace;
    };

    const std::vector<std::array<DiracSpinor::Index, 3>> configs{
      {index(2, -1), index(2, 1), index(3, 2)},   // 2s 2p- 3d-
      {index(2, 1), index(3, 1), index(4, 1)},    // 2p- 3p- 4p-
      {index(2, -2), index(2, -2), index(3, -1)}, // 2p+^2 3s
      {index(2, -1), index(2, -2), index(2, -2)}, // 2s 2p+^2
      {index(2, -2), index(2, -2), index(2, -2)}, // 2p+^3
      {index(3, -3), index(3, -3), index(3, -3)}, // 3d+^3
      {index(3, -3), index(3, -3), index(4, -3)} // 3d+^2 4d+
    };

    for (auto config : configs) {
      std::sort(config.begin(), config.end());
      const auto [x, y, z] = config;
      const auto nx = Angular::nkindex_to_twoj(x) + 1;
      const auto ny = Angular::nkindex_to_twoj(y) + 1;
      const auto nz = Angular::nkindex_to_twoj(z) + 1;

      // Number of determinants, and two-body trace over them
      double num_det = 0.0;
      double trace_g = 0.0;
      if (x != y && y != z) {
        num_det = nx * ny * nz;
        trace_g =
          nz * pair_trace(x, y) + ny * pair_trace(x, z) + nx * pair_trace(y, z);
      } else if (x == y && y != z) {
        num_det = nx * (nx - 1) / 2.0 * nz;
        trace_g = nz * pair_trace(x, x) + (nx - 1) * pair_trace(x, z);
      } else if (x != y && y == z) {
        num_det = nx * ny * (ny - 1) / 2.0;
        trace_g = nx * pair_trace(y, y) + (ny - 1) * pair_trace(x, y);
      } else {
        num_det = nx * (nx - 1) * (nx - 2) / 6.0;
        trace_g = (nx - 2) * pair_trace(x, x);
      }
      const auto trace_h1 =
        num_det * (h1.getv(x, x) + h1.getv(y, y) + h1.getv(z, z));
      const auto expected = trace_h1 + trace_g;

      double trace = 0.0;
      for (int twoJ = 1; twoJ <= nx + ny + nz - 3; twoJ += 2) {
        for (const auto &csf : CI::form_CSF3s(config, twoJ)) {
          trace += (twoJ + 1) * CI::Hab(csf, csf, twoJ, h1, qk);
        }
      }
      worst_trace =
        std::max(worst_trace, std::abs(trace - expected) / std::abs(expected));
    }
  }

  fmt::print("<H> vs <h> + <V_HF>:    {:.1e}\n", worst_closed);
  fmt::print("E_3e vs E_2e + e_z:     {:.1e}\n", worst_spectator);
  fmt::print("Tr(H) vs Tr(det), rel.: {:.1e}\n", worst_trace);
  REQUIRE(worst_closed < 1.0e-10);
  REQUIRE(worst_spectator < 1.0e-10);
  REQUIRE(worst_trace < 1.0e-12);
}

//==============================================================================
// Three-electron CI, run through configuration_interaction: input parsing
// (fractions, twoJ), the odd-2J blocks, the ci file, and the plumbing
// (block energies reproduce the directly constructed Hamiltonian)
TEST_CASE("CI: three-electron run", "[CI][3VCI][unit]") {

  std::cout << "CI three-electron run\n";

  // Li, all three electrons in the CI space (bare-nucleus basis)
  Wavefunction wf({400, 1.0e-4, 20.0, 0.33 * 20.0, "loglinear"},
                  {"Li", -1, "pointlike"}, 1.0);
  wf.solve_core("HartreeFock", "[]", std::nullopt, 1.0e-5, false);
  wf.formBasis(SplineBasis::Parameters("3sp", 20, 6, {1.0e-2}, 20.0));

  const std::string qk_filename =
    "deleteme_" + qip::random_string(3) + ".qk.abf";
  const std::string ci_filename =
    "deleteme_" + qip::random_string(3) + ".ci.abf";
  const auto files = "qk_file = " + qk_filename + "; ci_file = " + ci_filename +
                     "; num_solutions = 3; print_details = false;";

  const IO::InputBlock input{"CI", "ci_basis = 3sp; num_electrons = 3;"
                                   "J = 1/2, 3/2;" +
                                     files};
  const auto solution = CI::configuration_interaction(input, wf);
  const auto &levels = solution.levels;
  REQUIRE(solution.integrals.num_electrons == 3);
  REQUIRE(levels.size() == 4);

  const auto &basis = solution.integrals.ci_basis;
  for (const auto &psi : levels) {
    REQUIRE(psi.num_electrons() == 3);
    REQUIRE(psi.twoJ() % 2 == 1);
    REQUIRE(psi.num_CSFs() ==
            CI::form_CSF3s(psi.twoJ(), psi.parity(), basis).size());
    REQUIRE(psi.num_solutions() == 3);
    // Energies reproduce the directly constructed Hamiltonian
    const auto Hci = CI::construct_Hci(psi, solution.integrals);
    const auto E = LinAlg::symmhEigensystem(Hci).first;
    for (std::size_t i = 0; i < psi.num_solutions(); ++i) {
      REQUIRE(psi.energy(i) == Approx(E[i]).epsilon(1.0e-10));
    }
    // Every CSF has three electrons, and the config strings reflect that
    REQUIRE(psi.csf(0).num_electrons() == 3);
    REQUIRE(psi.csf(0).parity() == psi.parity());
  }

  // Ground state: 1s^2 ns, J=1/2 even (in the bare-nucleus basis, the 2s
  // and 3s spline orbitals carry similar weight). Bounds: above the exact
  // energy (-7.478 au, non-rel), and below the single-configuration limit
  const auto *gs = &levels.front();
  REQUIRE(gs->twoJ() == 1);
  REQUIRE(gs->parity() == 1);
  REQUIRE(gs->info(0).config.rfind("1s^2 ", 0) == 0);
  REQUIRE(gs->info(0).config.back() == 's');
  REQUIRE(gs->energy(0) > -7.48);
  REQUIRE(gs->energy(0) < -7.3);
  fmt::print("Li 1s^2 2s 2S1/2: E = {:.6f} au ({:.0f}% {})\n", gs->energy(0),
             100.0 * gs->info(0).ci2, gs->info(0).config);

  // Same blocks via twoJ (which supersedes J), read from the ci file
  // (read_only: a failed read would leave the blocks unsolved)
  const IO::InputBlock input_twoJ{"CI",
                                  "ci_basis = 3sp; num_electrons = 3;"
                                  "J = 2.5; twoJ = 1, 3; read_only = true;" +
                                    files};
  const auto levels_twoJ = CI::configuration_interaction(input_twoJ, wf).levels;
  REQUIRE(levels_twoJ.size() == levels.size());
  for (std::size_t i = 0; i < levels.size(); ++i) {
    REQUIRE(levels_twoJ[i].twoJ() == levels[i].twoJ());
    REQUIRE(levels_twoJ[i].parity() == levels[i].parity());
    for (std::size_t s = 0; s < levels[i].num_solutions(); ++s) {
      REQUIRE(levels_twoJ[i].energy(s) == Approx(levels[i].energy(s)));
      for (std::size_t k = 0; k < levels[i].num_CSFs(); ++k) {
        REQUIRE(levels_twoJ[i].coef(s, k) == Approx(levels[i].coef(s, k)));
      }
    }
  }

  // Integer J is not possible for three electrons: no blocks
  std::cout << "\nExpect one warning below (deliberate): J = 1 is not possible "
               "for 3 electrons\n";
  const IO::InputBlock input_bad{"CI", "ci_basis = 3sp; num_electrons = 3;"
                                       "J = 1;" +
                                         files};
  REQUIRE(CI::configuration_interaction(input_bad, wf).levels.empty());

  // Parity-specific list supersedes the plain one, for its parity only:
  // odd parity takes twoJ- (rejected), even parity takes J
  std::cout << "\nExpect one warning below (deliberate, from twoJ- = 2): J = 1 "
               "is not possible for 3 electrons\n";
  const IO::InputBlock input_pi{"CI", "ci_basis = 3sp; num_electrons = 3;"
                                      "J = 1/2; twoJ- = 2; read_only = true;" +
                                        files};
  const auto levels_pi = CI::configuration_interaction(input_pi, wf).levels;
  REQUIRE(levels_pi.size() == 1);
  REQUIRE(levels_pi[0].twoJ() == 1);
  REQUIRE(levels_pi[0].parity() == 1);
  REQUIRE(levels_pi[0].num_solutions() == 3);
}
