#include "Amplitudes/Normalisation.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <iostream>

//==============================================================================
TEST_CASE("Amplitudes: dSigma_dE", "[Amplitudes][MBPT][unit]") {

  Wavefunction wf({400, 1.0e-4, 50.0, 33.0, "loglinear"}, {"Na", -1, "Fermi"},
                  1.0);
  wf.solve_core("HartreeFock", "[Ne]", std::nullopt, 1.0e-5);
  wf.solve_valence("3sp");
  // wf.formBasis({"20spd", 20, 6, {3.4e-4, 1.7e-2}, 30.0});
  wf.formBasis({"20spd", 20, 6, {1.0e-4}, 30.0});

  // Small Sigma; parameters not meant to be accurate
  const double r0{1.0e-2};
  const double rmax{20.0};
  const std::size_t stride = 6;
  const int n_min_core = 2;
  const auto delta = 1.0e-4;

  auto Sigma0 = MBPT::CorrelationPotential(
    "", wf.vHF(), wf.basis(), r0, rmax, stride, n_min_core,
    MBPT::SigmaMethod::Goldstone, false, false);
  auto Sigma1 = Sigma0;
  auto Sigma2 = Sigma0;
  for (const auto &v : wf.valence()) {
    Sigma0.formSigma(v.kappa(), v.en(), v.n(), &v);
    Sigma1.formSigma(v.kappa(), v.en() + delta, v.n(), &v);
    Sigma2.formSigma(v.kappa(), v.en() - delta, v.n(), &v);
  }

  for (const auto &v : wf.valence()) {
    const auto dS = Amplitudes::dSigma_dE(v, Sigma0, Sigma1, Sigma2, delta);

    // Equals the central difference of the correlation energy de(e) = <v|S|v>
    const auto *S1 = Sigma1.getSigma(v.kappa(), v.n());
    const auto *S2 = Sigma2.getSigma(v.kappa(), v.n());
    REQUIRE(S1 != nullptr);
    REQUIRE(S2 != nullptr);
    const auto lambda = Sigma0.getLambda(v.kappa(), v.n());
    const auto manual =
      lambda * (v * (*S1 * v) - v * (*S2 * v)) / (2.0 * delta);
    REQUIRE(dS == Approx(manual).epsilon(1.0e-10));

    // dSigma/dE is negative (Sigma deepens with energy) and small
    REQUIRE(dS < 0.0);
    REQUIRE(std::abs(dS) < 0.1);
  }
}
