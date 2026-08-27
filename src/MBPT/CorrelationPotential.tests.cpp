#include "CorrelationPotential.hpp"
#include "DiracODE/include.hpp"
#include "Feynman.hpp"
#include "Goldstone.hpp"
#include "Maths/Grid.hpp"
#include "Sigma2.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include "fmt/format.hpp"
#include "qip/Maths.hpp"
#include "qip/Random.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <string>

//==============================================================================
TEST_CASE("MBPT: Goldstone, unit tests", "[MBPT][Goldstone][unit]") {

  std::cout << "Goldstone diagram, unit tests (not meant to be accurate)\n";

  Wavefunction wf({400, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-5);
  wf.solve_valence("3sp");
  wf.formBasis(SplineBasis::Parameters("20spd", 20, 6, 1.0e-4, 1.0e-4, 30.0));

  // These parameters are not meant to be accurate
  const double r0{1.0e-2};
  const double rmax{20.0};
  const std::size_t stride = 6;
  const int n_min_core = 2;

  const auto i0 = wf.grid().getIndex(r0);
  const auto size = (wf.grid().getIndex(rmax) - i0) / stride + 1;

  // Construct Goldtone diagrams (second-order only)
  MBPT::Goldstone Gs =
    MBPT::Goldstone(wf.basis(), wf.core(), i0, stride, size, n_min_core, false);

  // Test the "parameter" getters
  REQUIRE(Gs.stride() == stride);
  REQUIRE(Gs.n_min() == n_min_core);
  REQUIRE(Gs.lmax() == 2);

  const auto &Yeh = Gs.Yeh();
  const auto &[holes, excited] = Gs.basis();

  for (auto &v : wf.valence()) {

    auto Sigma =
      Gs.Sigma_direct(v.kappa(), v.en()) + Gs.Sigma_exchange(v.kappa(), v.en());

    auto de0 = MBPT::Sigma_vw(v, v, Yeh, holes, excited);

    auto de1 = v * (Sigma * v);

    // Not an accuracy test
    REQUIRE(de1 == Approx(de0).epsilon(1.0e-3));
  }
}

//==============================================================================
TEST_CASE("MBPT: Feynman unit tests", "[MBPT][Feynman][unit]") {

  std::cout << "Feynman diagram, unit tests (not meant to be accurate)\n";

  Wavefunction wf({1000, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-5);
  wf.solve_valence("3sp");

  // These parameters are not meant to be accurate
  const double r0{1.0e-3};
  const double rmax{30.0};
  const std::size_t stride = 8;
  const auto omre = -0.33 * wf.energy_gap();
  const int lmax = 2;
  const double w0{0.1};
  const double wratio{3.0};
  const int n_min_core = 2;

  const auto i0 = wf.grid().getIndex(r0);
  const auto size = (wf.grid().getIndex(rmax) - i0) / stride + 1;

  // Construct Feynman diagrams (second-order only)
  MBPT::Feynman Fy(wf.vHF(), i0, stride, size,
                   {MBPT::Screening::exclude, MBPT::HoleParticle::exclude, lmax,
                    omre, w0, wratio},
                   n_min_core, true);

  // Test the "parameter" getters
  REQUIRE(Fy.screening() == false);
  REQUIRE(Fy.hole_particle() == false);
  REQUIRE(Fy.stride() == stride);
  REQUIRE(Fy.n_min() == n_min_core);
  REQUIRE(Fy.lmax() == lmax);
  REQUIRE(Fy.w0() == Approx(w0));
  REQUIRE(Fy.wratio() == Approx(wratio));
  REQUIRE(Fy.omre() == Approx(omre));

  // Test exchange potential: just that it works, not an accuracy test
  for (auto &v : wf.valence()) {
    const auto vx = Fy.get_Vx_kappa(v.kappa());
    const auto dv = vx * v;
    const auto dex = qip::inner_product(v.f(), dv.f()) / double(Fy.stride());
    const auto dex0 = v * (wf.vHF()->vexFa(v));
    const auto eps = std::abs(dex / dex0 - 1.0);
    REQUIRE(eps < 0.1);
  }

  // with screening
  MBPT::Feynman FyS(wf.vHF(), i0, stride, size,
                    {MBPT::Screening::include, MBPT::HoleParticle::exclude,
                     lmax, omre, w0, wratio},
                    n_min_core, false);

  // with screening + hole-particle (all-order)
  MBPT::Feynman FyAO(wf.vHF(), i0, stride, size,
                     {MBPT::Screening::include, MBPT::HoleParticle::include,
                      lmax, omre, w0, wratio},
                     n_min_core, false);

  // expected data
  const std::vector expected_de{-0.0049, -0.0015, -0.0015};
  const std::vector expected_sc_ratio{0.85, 0.90, 0.90};
  const std::vector expected_hp_ratio{1.22, 1.23, 1.23};
  const double epsilon = 0.25; // just test to 25% (not anccuracy test)

  for (std::size_t i = 0; i < wf.valence().size(); ++i) {
    const auto &v = wf.valence().at(i);

    // Check second-order Feynman against expected
    const auto Sd = Fy.Sigma_direct(v.kappa(), v.en());
    const auto de0 = v * (Sd * v);
    // only require to ~1%, since not an accuracy test
    REQUIRE(de0 == Approx(expected_de[i]).epsilon(epsilon));

    // Test screening and hole-particle:
    const auto Sd_S = FyS.Sigma_direct(v.kappa(), v.en());
    const auto Sd_ao = FyAO.Sigma_direct(v.kappa(), v.en());
    const auto des = v * (Sd_S * v);
    const auto deao = v * (Sd_ao * v);

    // Test screening: ratio of screened to unscreened:
    REQUIRE(des / de0 == Approx(expected_sc_ratio[i]).epsilon(epsilon));

    // Test hole-particle: ratio of all-orders to screening
    REQUIRE(deao / des == Approx(expected_hp_ratio[i]).epsilon(epsilon));
  }

  // Effective screening factors of the dressed Coulomb line: per k, the
  // ratio of the direct diagram to the same with the screening removed
  // (same hp). Compare against the ratio formed from separate objects
  // (screening included/excluded; same include_G as FyS/FyAO), which tests
  // the unscreened Q*Pi*Q formed by the screened object, and the per-state
  // store. References: no screening, without/with hp
  MBPT::Feynman Fy0(wf.vHF(), i0, stride, size,
                    {MBPT::Screening::exclude, MBPT::HoleParticle::exclude,
                     lmax, omre, w0, wratio},
                    n_min_core, false, false);
  MBPT::Feynman FyH(wf.vHF(), i0, stride, size,
                    {MBPT::Screening::exclude, MBPT::HoleParticle::include,
                     lmax, omre, w0, wratio},
                    n_min_core, false, false);
  for (const auto &v : wf.valence()) {
    // Without screening, all factors are 1
    for (const auto f : Fy0.screening_fk(v)) {
      REQUIRE(f == 1.0);
    }
    const auto Sd_k = Fy0.Sigma_direct_each_k(v.kappa(), v.en());
    const auto SdS_k = FyS.Sigma_direct_each_k(v.kappa(), v.en());
    const auto SdH_k = FyH.Sigma_direct_each_k(v.kappa(), v.en());
    const auto SdAO_k = FyAO.Sigma_direct_each_k(v.kappa(), v.en());
    const auto fk_S = FyS.screening_fk(v);
    const auto fk_AO = FyAO.screening_fk(v);
    REQUIRE(fk_S.size() == Sd_k.size());
    REQUIRE(fk_AO.size() == Sd_k.size());
    bool tested_nontrivial = false;
    for (auto k = 0ul; k < Sd_k.size(); ++k) {
      const auto de0 = v * (Sd_k[k] * v);
      const auto deH = v * (SdH_k[k] * v);
      if (de0 == 0.0 || deH == 0.0)
        continue;
      const auto expected_S = (v * (SdS_k[k] * v)) / de0;
      const auto expected_AO = (v * (SdAO_k[k] * v)) / deH;
      REQUIRE(fk_S[k] == Approx(expected_S).epsilon(1.0e-6));
      REQUIRE(fk_AO[k] == Approx(expected_AO).epsilon(1.0e-6));
      if (std::abs(expected_S - 1.0) > 0.01) {
        tested_nontrivial = true;
      }
    }
    REQUIRE(tested_nontrivial);
    // Stored: a second call (at any energy) returns the same factors
    REQUIRE(FyS.screening_fk(v, v.en() + 0.01) == fk_S);
  }
  std::cout << "\n";
}

//==============================================================================
TEST_CASE("MBPT: Feynman complex Green",
          "[MBPT][Feynman][ComplexGreen][unit]") {

  // Tests the two methods for the Green's function at complex energy:
  //  (1) Dyson method: solve at Re(en), correct to complex via resolvent
  //  (2) Direct method: solve the Dirac equation at complex energy
  // Each is checked against the exact <a|G|a> = 1/(en - e_a), and against the
  // symmetries G(en*) = G(en)* and G(r1,r2) = G(r2,r1).

  fmt::print("Green's function at complex energy: Dyson (1) vs direct (2), Na\n"
             "  cc     : conjugate symmetry, g(en*) vs g(en)* (direct method)\n"
             "  sym    : complex symmetry, g(r1,r2) vs g(r2,r1)\n"
             "  core1,2: <a|G|a> vs. <a|a><a|a>*(en-e_a) - 1, core states a\n"
             "  val1,2 : same, for the n=5 state (valence region)\n");

  Wavefunction wf({1000, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]");
  wf.formBasis(SplineBasis::Parameters("30spdf", 40, 7, 1.0e-4, 0.0, 40.0, "",
                                       SplineBasis::SplineType::Derevianko,
                                       false, false));

  const double r0{1.0e-3};
  const double rmax{30.0};
  const std::size_t stride = 4;
  const auto omre = -0.33 * wf.energy_gap();
  const int lmax = 4;
  const double w0{1.0};
  const double wratio{100.0};
  const int n_min_core = 2;

  const auto i0 = wf.grid().getIndex(r0);
  const auto size = (wf.grid().getIndex(rmax) - i0) / stride + 1;

  // Green's functions only (Q*Pi*Q is never formed: only calculated on use)
  const MBPT::Feynman Fy(wf.vHF(), i0, stride, size,
                         {MBPT::Screening::exclude, MBPT::HoleParticle::exclude,
                          lmax, omre, w0, wratio},
                         n_min_core, true, false, "");

  // Test at typical energies: en_v + omre + iw, and e_core +/- (omre + iw)
  const std::vector<std::complex<double>> energies{{-0.5, 0.1},  {-0.5, 2.0},
                                                   {-0.5, 30.0}, {-2.0, 0.5},
                                                   {-3.0, 1.0},  {-1.5, 50.0}};

  // Only kappas that occur in the core: these have core states (for the pole
  // test) and are spanned by the basis (for the basis comparison)
  const auto max_ki = Angular::l_to_max_kindex(DiracSpinor::max_l(wf.core()));

  fmt::print("\n{:>5} {:>8} {:>7} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}\n",
             "kappa", "Re(en)", "Im(en)", "cc", "sym", "core1", "core2", "val1",
             "val2");

  int num_pole_tests{0};

  for (auto ki = 0ul; ki <= max_ki; ++ki) {
    const auto kappa = Angular::kindex_to_kappa(ki);
    for (const auto &en : energies) {

      const auto g1 = Fy.green_dyson(kappa, en);
      const auto g2 = Fy.green_complex_dirac(kappa, en);

      // Conjugate symmetry for the direct-complex method
      const auto g2cc = Fy.green_complex_dirac(kappa, std::conj(en));
      const auto eps_cc = MBPT::max_delta(g2cc, g2.conj()) / max_element(g2);

      // G at complex energy is complex-SYMMETRIC (transpose), not Hermitian:
      // ff(1,2) = ff(2,1), gg(1,2) = gg(2,1), and fg(1,2) = gf(2,1)
      const auto asymmetry = [](const MBPT::ComplexGMatrix &G) {
        double dev{0.0};
        for (auto i = 0ul; i < G.size(); ++i) {
          for (auto j = 0ul; j < G.size(); ++j) {
            dev = std::max(dev, std::abs(G.ff(i, j) - G.ff(j, i)));
            dev = std::max(dev, std::abs(G.gg(i, j) - G.gg(j, i)));
            dev = std::max(dev, std::abs(G.fg(i, j) - G.gf(j, i)));
          }
        }
        return dev;
      };
      const auto eps_sym = std::max(asymmetry(g1) / MBPT::max_element(g1),
                                    asymmetry(g2) / MBPT::max_element(g2));

      // <a|G|a>, with the sub-grid integration measure: the same contraction
      // the correlation potential is evaluated with, v * (Sigma * v)
      const auto braket = [](const DiracSpinor &Fa,
                             const MBPT::ComplexGMatrix &G) {
        const auto Gdr = G.drj();
        return std::complex<double>{Fa * (Gdr.real() * Fa),
                                    Fa * (Gdr.imag() * Fa)};
      };

      // Sharp analytic test: <a|G(en)|a> = 1/(en - e_a) for any HF eigenstate.
      // Taken as a ratio to the same braket of the exact single-pole Green's
      // function, |a><a|/(en - e_a), which cancels the sub-grid quadrature of
      // the orbital and leaves the error in G itself.
      const auto pole_test = [&](const DiracSpinor &Fa,
                                 const MBPT::ComplexGMatrix &G) {
        const auto exact = Fy.green_basis(Fa.kappa(), en, {Fa});
        return std::abs(braket(Fa, G) / braket(Fa, exact) - 1.0);
      };

      // Core states: sample G at small r
      double core1{0.0}, core2{0.0};
      for (const auto &Fa : wf.core()) {
        if (Fa.kappa() != kappa)
          continue;
        core1 = std::max(core1, pole_test(Fa, g1));
        core2 = std::max(core2, pole_test(Fa, g2));
        ++num_pole_tests;
      }

      // Excited state: samples G over the valence region, where Sigma is built
      const auto *Fv = DiracSpinor::find(3, kappa, wf.basis());
      REQUIRE(Fv != nullptr);
      const auto val1 = pole_test(*Fv, g1);
      const auto val2 = pole_test(*Fv, g2);

      fmt::print("{:5} {:8.1f} {:7.1f} {:9.1e} {:9.1e} {:9.1e} {:9.1e} "
                 "{:9.1e} {:9.1e}\n",
                 kappa, en.real(), en.imag(), eps_cc, eps_sym, core1, core2,
                 val1, val2);

      CHECK(eps_cc < 1.0e-10);
      CHECK(eps_sym < 1.0e-5);

      // The Dyson G satisfies the discrete resolvent algebra, so its bound-pole
      // content is accurate at every Im(en); its error appears instead in the
      // diffuse probe, which weights the continuum part of G
      CHECK(core1 < 0.005);
      CHECK(val1 < 0.05);

      // The direct method samples the true kernel: accurate while the sub-grid
      // resolves the near-diagonal ridge (width ~ 1/sqrt(2*Im(en))), aliasing
      // above that. This is why green() only uses it below Im(en) ~ 10
      if (std::abs(en.imag()) < 10.0) {
        CHECK(core2 < 0.01);
        CHECK(val2 < 0.001);
      } else {
        CHECK(core2 < 0.05);
        CHECK(val2 < 0.2);
      }
    }
  }
  // guard against the pole tests passing vacuously (no core states tested)
  REQUIRE(num_pole_tests > 0);
  std::cout << "\n";
}

//==============================================================================
TEST_CASE("MBPT: Feynman omre stability", "[MBPT][Feynman][omre]") {

  // Sigma is (analytically) independent of the contour position Re(w) = omre;
  // the numerical spread against omre measures the total error
  // (frequency quadrature + Green's function accuracy).
  // Compares the two complex-Green methods: Dyson vs direct-complex.

  fmt::print("Sigma_direct(3s, Na) vs contour position, omre.\n"
             "Exact Sigma is omre-independent: spread = numerical error.\n");

  Wavefunction wf({1000, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-6);
  wf.solve_valence("3s");
  wf.formBasis(SplineBasis::Parameters("30spdf", 35, 7, 1.0e-4, 1.0e-4, 40.0,
                                       "", SplineBasis::SplineType::Derevianko,
                                       false, false));

  const double r0{1.0e-3};
  const double rmax{30.0};
  const std::size_t stride = 4;
  const int lmax = DiracSpinor::max_l(wf.basis());
  const int n_min_core = 2;

  const auto i0 = wf.grid().getIndex(r0);
  const auto size = (wf.grid().getIndex(rmax) - i0) / stride + 1;

  const auto &v = wf.valence().front();

  const auto wr_best = MBPT::best_omre(wf.core(), wf.valence(), true);

  // Goldstone (basis) value: the 2nd-order direct reference
  // (differs from Feynman by basis truncation + negative-energy states)
  const MBPT::Goldstone Gs(wf.basis(), wf.core(), i0, stride, size, n_min_core,
                           true);
  const auto de_G = v * (Gs.Sigma_direct(v.kappa(), v.en()) * v);
  fmt::print("Goldstone reference: de = {:.6f} au\n", de_G);

  std::vector<double> omres{-0.15, -0.30, -0.60};
  omres.push_back(wr_best);

  // nb: w0 = 0.01 (not larger): the [0,w0] trapezoid panel requires the
  // integrand to be flat across the panel (w0 small vs pole distances)
  const double w0{0.01};
  const double wratio{1.5};

  for (const bool complex_green : {false, true}) {
    std::cout << (complex_green ?
                    "\nDirectly solving Green's function at complex energy" :
                    "\nUsing Dyson equation to shift Green's to complex energy")
              << "\n";
    std::vector<double> des;
    for (const auto omre : omres) {
      MBPT::Feynman Fy(wf.vHF(), i0, stride, size,
                       {MBPT::Screening::exclude, MBPT::HoleParticle::exclude,
                        lmax, omre, w0, wratio, complex_green},
                       n_min_core, true, false);
      const auto Sd = Fy.Sigma_direct(v.kappa(), v.en());
      des.push_back(v * (Sd * v));
      fmt::print("  omre = {:6.3f} : de = {:.6f}\n", omre, des.back());
    }

    // Compare at the recommended omre (the last entry of omres)
    const auto de_best = des.back();
    const auto [min, max] = std::minmax_element(des.cbegin(), des.cend());
    const auto spread = std::abs((*max - *min) / de_best);
    fmt::print("Spread = {:.1e}, \nvs Goldstone = {:.1e}\n", spread,
               (de_best / de_G - 1.0));

    // nb: tolerances are empirical for THIS config (regression guards, not
    // universal statements: e.g., at the Cs production config the Dyson
    // method agrees with Goldstone to ~1%)
    if (complex_green) {
      // Direct-complex: omre-stable, and agrees with Goldstone limit here
      CHECK(spread < 0.0001);
      CHECK(std::abs(de_best / de_G - 1.0) < 0.01);
    } else {
      CHECK(spread < 0.001);
      CHECK(std::abs(de_best / de_G - 1.0) < 0.01);
    }
  }
  std::cout << "\n";
}

//==============================================================================
TEST_CASE("MBPT: CorrelationPotential", "[MBPT][CorrelationPotential][unit]") {

  Wavefunction wf({400, 1.0e-4, 50.0, 0.33 * 100.0, "loglinear"},
                  {"Na", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-5);
  wf.solve_valence("3sp");
  wf.formBasis(SplineBasis::Parameters("20spd", 20, 6, 1.0e-4, 1.0e-4, 30.0));

  // These parameters are not meant to be accurate
  const double r0{1.0e-2};
  const double rmax{20.0};
  const std::size_t stride = 6;
  const int n_min_core = 2;

  auto SigmaG = MBPT::CorrelationPotential(
    "", wf.vHF(), wf.basis(), r0, rmax, stride, n_min_core,
    MBPT::SigmaMethod::Goldstone, false, false);

  REQUIRE(SigmaG.empty());

  for (auto &v : wf.valence()) {
    SigmaG.formSigma(v.kappa(), v.en(), v.n(), &v);
  }

  REQUIRE(!SigmaG.empty());

  std::string file_name = "deleteme_" + qip::random_string(6) + ".sig.abf";

  // test write and read:
  SigmaG.write(file_name);

  // read in:
  auto SigmaG2 = MBPT::CorrelationPotential(
    file_name, wf.vHF(), wf.basis(), r0, rmax, stride, n_min_core,
    MBPT::SigmaMethod::Goldstone, false, false);

  SigmaG2.print_subGrid();

  std::cout << "\n";

  for (auto &v : wf.valence()) {

    const auto de1 = v * SigmaG2(v);
    const auto de2 = v * SigmaG(v);

    REQUIRE(de1 == Approx(de2));

    const auto S = SigmaG2.getSigma(v.kappa(), v.n());
    REQUIRE(S != nullptr);
    const auto de3 = v * (*S * v);

    REQUIRE(de1 == Approx(de3));
  }

  //------------------------------------------------------

  // test scaling:
  SigmaG2.scale_Sigma({2.0, 2.0, 2.0});
  SigmaG2.print_scaling();
  for (auto &v : wf.valence()) {

    const auto de1 = v * SigmaG2(v);
    const auto de2 = 2.0 * v * SigmaG(v);

    REQUIRE(de1 == Approx(de2));
  }

  for (auto &v : wf.valence()) {
    SigmaG2.scale_Sigma(v.n() * 1.5, v.kappa(), v.n());
    const auto de1 = v * SigmaG2(v);
    const auto de2 = v.n() * 1.5 * v * SigmaG(v);
    REQUIRE(de1 == Approx(de2));
  }
}

//==============================================================================
//==============================================================================

//==============================================================================
//! Unit tests for second-order MBPT energy correction
TEST_CASE("MBPT: 2nd Order de", "[MBPT][integration]") {

  { // Compare with  K. Beloy and A. Derevianko,
    // Comput. Phys. Commun. 179, 310 (2008).
    Wavefunction wf({4000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Xe]");
    wf.solve_valence("6s");
    const auto &Fv = wf.valence().front();

    {
      // K. Beloy and A. Derevianko, Comput. Phys. Commun. 179, 310 (2008).
      const auto partial_KBAD =
        std::vector{-0.0000130, -0.0020027, -0.0105623, -0.0039347,
                    -0.0007563, -0.0002737, -0.0001182};
      // const auto total_KBAD = -0.0176609;
      const auto error = 0.0000002; // allow ony difference of 1.5 in last digit

      double prev = 0.0;
      std::vector<double> vals;
      wf.formBasis({"100spdfghi", 101, 11, 0.0, 1.0e-6, 50.0});
      // wf.formSigma(1, false);
      // Coulomb::YkTable yk(wf.basis()); // this is a *huge* basis!
      const auto [holes, excited] =
        DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
      Coulomb::YkTable yk(holes, excited);
      // const auto Sigma = wf.Sigma();
      std::cout << "cf Table 2 from Beloy, Derevianko, Comput.Phys.Commun. "
                   "179, 310 (2008):\n";
      for (int l = 0; l <= 6; ++l) {

        const auto de = MBPT::Sigma_vw(Fv, Fv, yk, holes, excited, l);
        // const auto de_2 = Sigma->Sigma_vw(Fv, Fv, l);
        // REQUIRE(de == Approx(de_2).epsilon(1.0e-9));
        vals.push_back(de - prev);
        printf("%i %10.7f %10.7f  [%10.7f]\n", l, de, de - prev,
               partial_KBAD[std::size_t(l)]);
        prev = de;
      }
      for (auto l = 0ul; l <= 6; ++l) {
        auto del = vals[l] - partial_KBAD[l];
        REQUIRE(std::abs(del) < error);
      }
    }

    { // "smaller" basis set (not exactly same as Derev)
      wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0});
      const auto [holes, excited] =
        DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
      Coulomb::YkTable yk(holes, excited);
      // wf.formSigma(1, false);
      // const auto Sigma = wf.Sigma();
      // const auto de = Sigma->Sigma_vw(Fv, Fv);
      const auto de = MBPT::Sigma_vw(Fv, Fv, yk, holes, excited);
      auto ok = de >= -0.01767 && de <= -0.01748 ? 1 : 0;
      // pass &= qip::check_value(&obuff, "MBPT(2) 'small' Cs 6s", ok, 1, 0);
      REQUIRE(ok);
    }
  }
}

//==============================================================================
//! Tests for second-order correlation potential
TEST_CASE("MBPT: Correlation Potential: Sigma2",
          "[MBPT][Sigma2][slow][integration]") {

  //----------------------------------------------------------------------------
  // Test Sigma:
  // Cs:

  const auto fname = std::string{"deleteme_"} + qip::random_string(5);

  std::vector<double> first_run;
  { // Compare Dzuba, using up to l=6 for splines
    std::cout << "Test Sigma(2) Brueckner, for Cs:\n";
    auto dzuba_i = std::vector{
      -0.02013813, -0.00410942, -0.00792483, -0.00702407, -0.00220878,
      -0.00199737, -0.01551449, -0.01466935, -0.00035253, -0.00035234};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Xe]");
    wf.solve_valence("7sp5d4f");
    wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0});
    wf.formSigma(3, 1.0e-4, 30.0, 14 /*stride*/, false, false, false, 0, {}, {},
                 {}, true, fname + "Cs");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();
    wf.printValence();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    std::cout << "delta Sigma(2) Bruckner, cf Dzuba:\n";
    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      const auto eps = std::abs((de[i] - dzuba_i[i]) / dzuba_i[i]);
      std::cout << de[i] << " [" << dzuba_i[i] << "] " << eps << "\n";
    }
    first_run = de; // copy this data, test the next run against:

    const auto [eps, at] = qip::compare_eps(dzuba_i, de);
    // pass &= qip::check_value(&obuff, "Sigma2 Cs", eps, 0.0, 0.01);
    REQUIRE(std::abs(eps) < 0.01);
  }

  std::cout << "\n";

  //----------------------------------------------------------------------------
  // Test reading in Sigma:
  // Cs:
  {
    std::cout << "Test reading in Sigma(2) Brueckner file, for Cs:\n";
    std::sort(begin(first_run), end(first_run)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Xe]");
    wf.solve_valence("7sp5d4f");
    // Don't calculate Sigma, read it in from above example:
    wf.formSigma(1, 0.0, 0.0, 1 /*stride*/, false, false, false, 0, {}, {}, {},
                 true, fname + "Cs");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();
    wf.printValence();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    for (auto i = 0ul; i < first_run.size(); ++i) {
      const auto eps = std::abs((de[i] - first_run[i]) / first_run[i]);
      std::cout << wf.valence().at(i) << " " << de[i] << " [" << first_run[i]
                << "] " << eps << "\n";
    }

    const auto [eps, at] = qip::compare_eps(first_run, de);
    // pass &= qip::check_value(&obuff, "Sigma2 Cs (read)", eps, 0.0, 1.0e-16);
    REQUIRE(std::abs(eps) < 1.0e-16);
  }

  std::cout << "\n";

  //----------------------------------------------------------------------------
  // Fr:
  { // Compare Dzuba, using up to l=6 for splines
    std::cout << "Test Sigma(2) Brueckner, for Fr:\n";
    auto dzuba_i =
      std::vector{-0.0245075, -0.0098094, -0.0069442, -0.0153430, -0.0133382};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({2000, 1.0e-6, 150.0, 0.33 * 150.0, "loglinear", -1.0},
                    {"Fr", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Rn]");
    wf.solve_valence("7sp6d");
    wf.formBasis({"30spdfghi", 40, 7, 0.0, 1.0e-6, 40.0});
    wf.formSigma(4, 1.0e-4, 30.0, 12 /*stride*/, false, false, false, 0, {}, {},
                 {}, true, fname + "Fr");

    std::vector<double> hf, br2;
    for (const auto &Fv : wf.valence()) {
      hf.push_back(Fv.en());
    }

    wf.hartreeFockBrueckner();

    for (const auto &Fv : wf.valence()) {
      br2.push_back(Fv.en());
    }

    auto de = qip::compose([](auto a, auto b) { return a - b; }, br2, hf);
    std::sort(begin(de), end(de)); // sort: don't depend on order
    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      const auto eps = std::abs((de[i] - dzuba_i[i]) / dzuba_i[i]);
      std::cout << de[i] << " [" << dzuba_i[i] << "] " << eps << "\n";
    }

    const auto [eps, at] = qip::compare_eps(dzuba_i, de);
    REQUIRE(std::abs(eps) < 0.02);
  }
}

//==============================================================================
//! Unit tests for all-orders correlation potential
TEST_CASE("MBPT: Correlation Potential: SigmaAO",
          "[MBPT][SigmaAO][slow][integration]") {

  { // Compare Dzuba, All-order sigma
    auto dzuba_i =
      std::vector{-0.14332871, -0.05844404, -0.09244689, -0.08985968,
                  -0.04392404, -0.04309476, -0.07812666, -0.07759564};
    std::sort(begin(dzuba_i), end(dzuba_i)); // sort: don't depend on order

    Wavefunction wf({4000, 1.0e-6, 120.0, 0.33 * 120.0, "loglinear", -1.0},
                    {"Cs", -1, "Fermi", -1.0, -1.0}, 1.0);
    wf.solve_core("HartreeFock", "[Xe]");
    wf.solve_valence("7sp5d");
    wf.formBasis({"35spdfghi", 40, 9, 0.0, 1.0e-6, 40.0});

    const auto n_min_core = 3;
    const auto rmin = 1.0e-4;
    const auto rmax = 30.0;
    const auto stride =
      int(wf.grid().getIndex(30.0) - wf.grid().getIndex(1.0e-4)) / 150;

    const auto omre = -std::abs(0.33 * wf.energy_gap());
    const double w0 = 0.01;
    const double wratio = 1.5;
    const auto lmax = 6;

    const std::vector fk{0.71, 0.589, 0.84, 0.885, 0.95, 0.976, 0.991};
    // const std::vector fk{0.72, 0.62, 0.83, 0.89, 0.94, 1.0};
    // wf.formSigma(3, true, 1.0e-4, 30.0, 14 /*stride*/);
    wf.formSigma(n_min_core, rmin, rmax, stride, false, false, false, 0, {}, fk,
                 {}, false, "", true, true, true, lmax, omre, w0, wratio);

    wf.hartreeFockBrueckner();

    std::vector<double> br;
    for (const auto &Fv : wf.valence()) {
      br.push_back(Fv.en());
    }
    std::sort(begin(br), end(br)); // sort: don't depend on order

    for (auto i = 0ul; i < dzuba_i.size(); ++i) {
      std::cout << wf.valence().at(i) << " " << br[i] << " [" << dzuba_i[i]
                << "]\n";
    }

    auto [eps, at] = qip::compare_eps(dzuba_i, br);
    // pass &= qip::check_value(&obuff, "Sigma all-orders Cs", eps, 0.0, 5e-04);
    // Used to be 5e-4..?
    REQUIRE(std::abs(eps) < 1.0e-2);
  }
}

//==============================================================================
TEST_CASE("MBPT: Sigma2", "[MBPT][Sigma2][CI][unit]") {

  // note: does not test formulas: just checks class is working correctly.
  // Not meant to be accurate!
  Wavefunction wf({400, 1.0e-3, 30.0, 10.0, "loglinear", -1.0},
                  {"Na", -1, "Fermi", -1.0, -1.0}, 1.0);
  wf.solve_core("Local", "[Ne]", std::nullopt, 1.0e-4);
  wf.solve_valence("3spd");
  wf.formBasis({"5spdf", 20, 5, 1.0e-2, 1.0e-2, 20.0});

  const auto &[core, excited] =
    MBPT::split_basis(wf.basis(), wf.FermiLevel(), 2);
  const auto mbpt_basis = qip::merge(core, excited);

  int kmax = 4;
  Coulomb::QkTable qk;
  Coulomb::YkTable yk(mbpt_basis);
  const Angular::SixJTable &SixJ = yk.SixJ();
  qk.fill(mbpt_basis, yk, kmax);

  for (auto &v : wf.valence()) {
    for (auto &w : wf.valence()) {
      for (auto &x : wf.valence()) {
        for (auto &y : wf.valence()) {
          for (int k = 0; k <= kmax; ++k) {

            //  Lk symmetry:
            //  {abcd} = badc

            const double sk1 = MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited,
                                             SixJ, MBPT::Denominators::Fermi0);

            const double sk2 = MBPT::Sk_vwxy(k, w, v, y, x, qk, core, excited,
                                             SixJ, MBPT::Denominators::Fermi0);

            const double sk3 =
              MBPT::Sigma2::S_Sigma2_ab(k, v, w, x, y, qk, core, excited, SixJ,
                                        MBPT::Denominators::Fermi0) +
              MBPT::Sigma2::S_Sigma2_c1(k, v, w, x, y, qk, core, excited, SixJ,
                                        MBPT::Denominators::Fermi0) +
              MBPT::Sigma2::S_Sigma2_c2(k, v, w, x, y, qk, core, excited, SixJ,
                                        MBPT::Denominators::Fermi0) +
              MBPT::Sigma2::S_Sigma2_d(k, v, w, x, y, qk, core, excited, SixJ,
                                       MBPT::Denominators::Fermi0);

            // tests symmetry:
            REQUIRE(sk2 == Approx(sk1));

            // tests formula (SR is checked in Sk_vwxy, not in internal):
            REQUIRE(sk3 == Approx(sk1));

            // Check selectrion rules
            if (MBPT::Sk_vwxy_SR(k, v, w, x, y)) {
              // This isn't always true... ?? Means SRs should be improves!
              // REQUIRE(sk3 != 0.0);
            } else {
              REQUIRE(sk3 == 0.0);
            }

            // DFK denominators: Lk and bra-ket (Hermitian) symmetries
            const double dfk1 = MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited,
                                              SixJ, MBPT::Denominators::DFK);
            const double dfk2 = MBPT::Sk_vwxy(k, w, v, y, x, qk, core, excited,
                                              SixJ, MBPT::Denominators::DFK);
            const double dfk3 = MBPT::Sk_vwxy(k, x, y, v, w, qk, core, excited,
                                              SixJ, MBPT::Denominators::DFK);
            REQUIRE(dfk2 == Approx(dfk1));
            REQUIRE(dfk3 == Approx(dfk1));

            // BW denominators: same symmetries, for an arbitrary E0
            const auto E0_test = -0.5;
            const double bw1 =
              MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited, SixJ,
                            MBPT::Denominators::BW, {}, E0_test);
            const double bw2 =
              MBPT::Sk_vwxy(k, w, v, y, x, qk, core, excited, SixJ,
                            MBPT::Denominators::BW, {}, E0_test);
            const double bw3 =
              MBPT::Sk_vwxy(k, x, y, v, w, qk, core, excited, SixJ,
                            MBPT::Denominators::BW, {}, E0_test);
            REQUIRE(bw2 == Approx(bw1));
            REQUIRE(bw3 == Approx(bw1));
          }
        }
      }
    }
  }

  // For a diagonal element, S^k_vwvw, setting E0 to the energy of that pair
  // makes every BW denominator equal to the RS one: BW is
  // E0 - (intermediate valence energies), and RS is what that becomes when
  // E0 = e_v + e_w. (Non-trivial for diagrams c1, c2 and d, whose external
  // parts do not vanish for a diagonal element.)
  for (const auto &v : wf.valence()) {
    for (const auto &w : wf.valence()) {
      const auto E0 = v.en() + w.en();
      for (int k = 0; k <= kmax; ++k) {
        const double s_bw = MBPT::Sk_vwxy(k, v, w, v, w, qk, core, excited,
                                          SixJ, MBPT::Denominators::BW, {}, E0);
        const double s_rs = MBPT::Sk_vwxy(k, v, w, v, w, qk, core, excited,
                                          SixJ, MBPT::Denominators::RS);
        REQUIRE(s_bw == Approx(s_rs));
      }
    }
  }

  // When every external leg is the lowest excited state of its kappa,
  // e_bar equals the actual energy, so DFK, RS, and Fermi all coincide
  std::vector<DiracSpinor> lowest;
  for (const auto &n : excited) {
    const auto same_kappa = [&n](const auto &m) {
      return m.kappa() == n.kappa();
    };
    if (n.l() <= 2 && std::find_if(lowest.begin(), lowest.end(), same_kappa) ==
                        lowest.end()) {
      lowest.push_back(n);
    }
  }
  for (const auto &v : lowest) {
    for (const auto &w : lowest) {
      for (const auto &x : lowest) {
        for (const auto &y : lowest) {
          for (int k = 0; k <= kmax; ++k) {
            const double s_dfk = MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited,
                                               SixJ, MBPT::Denominators::DFK);
            const double s_rs = MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited,
                                              SixJ, MBPT::Denominators::RS);
            const double s_f = MBPT::Sk_vwxy(k, v, w, x, y, qk, core, excited,
                                             SixJ, MBPT::Denominators::Fermi);
            REQUIRE(s_dfk == Approx(s_rs));
            REQUIRE(s_dfk == Approx(s_f));
          }
        }
      }
    }
  }
}
