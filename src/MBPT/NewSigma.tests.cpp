#include "MBPT/NewSigma.hpp"
#include "DiracODE/DiracODE.hpp"
#include "MBPT/Feynman.hpp"
#include "MBPT/Goldstone.hpp"
#include "MBPT/NewSigma.hpp" // just compile yests
#include "MBPT/Sigma2.hpp"
#include "Maths/Grid.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "catch2/catch.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>

//==============================================================================
TEST_CASE("MBPT: NewSigma", "[MBPT][NewSigma]") { //

  Wavefunction wf({2000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
                  {"Cs", -1, "Fermi"}, 1.0);
  wf.solve_core("HartreeFock", 0.0, "[Xe]");
  wf.formBasis({"35spdfg", 40, 7, 1.0e-3, 1.0e-4, 40.0, false});
  wf.solve_valence("7sp5d");
  wf.printCore();

  // MBPT::Goldstone Gold(wf.basis(), wf.core(),
  //                      MBPT::rgrid_params{1.0e-4, 30.0, 6}, 1, true);

  // const auto [holes, excited] =
  //     DiracSpinor::split_by_energy(wf.basis(), wf.FermiLevel());
  // Coulomb::YkTable yk(holes, excited);

  // std::cout << "\nSigma(2):\n";
  // for (auto &v : wf.valence()) {
  //   const auto Sd = Gold.Sigma_direct(v.kappa(), v.en());
  //   const auto Sx = Gold.Sigma_exchange(v.kappa(), v.en());

  //   const auto de1 = v * ((Sd + Sx) * v);

  //   const auto de0 = MBPT::Sigma_vw(v, v, yk, holes, excited);

  //   std::cout << v << " " << v.en() << " " << de0 << " " << de1 << "\n";
  // }

  // auto Sigma = MBPT::NewSigma("test.sig", wf.vHF(), wf.basis(), 1.0e-4, 30.0, 8,
  //                             3, MBPT::SigmaMethod::Goldstone, true);

  auto SigmaF =
      MBPT::NewSigma("test.sig", wf.vHF(), wf.basis(), 1.0e-4, 30.0, 8, 3,
                     MBPT::SigmaMethod::Feynman, true,
                     {MBPT::Screening::include, MBPT::HoleParticle::include, 4,
                      -0.33 * wf.energy_gap(), 0.01, 1.5});

  // for (auto &v : wf.valence()) {
  //   Sigma.formSigma(v.kappa(), v.en(), v.n(), &v);
  //   // SigmaF.formSigma(v.kappa(), v.en(), v.n(), &v);
  // }

  for (auto &v : wf.valence()) {
    // Sigma.formSigma(v.kappa(), v.en(), v.n(), &v);
    SigmaF.formSigma(v.kappa(), v.en(), v.n(), &v);
  }

  // std::cout << "\n";

  // for (auto &v : wf.valence()) {

  //   const auto de1 = v * Sigma(v);
  //   const auto de2 = v * SigmaF(v);

  //   // const auto de0 = MBPT::Sigma_vw(v, v, yk, holes, excited);

  //   std::cout << v << " " << v.en() << " " << de1 << " " << de2 << "\n";
  // }

  SigmaF.write("out.sig");

  auto SigmaF2 =
      MBPT::NewSigma("out.sig", wf.vHF(), wf.basis(), 1.0e-4, 30.0, 8, 3,
                     MBPT::SigmaMethod::Feynman, true,
                     {MBPT::Screening::include, MBPT::HoleParticle::include, 4,
                      -0.33 * wf.energy_gap(), 0.01, 1.5});

  for (auto &v : wf.valence()) {
    // Sigma.formSigma(v.kappa(), v.en(), v.n(), &v);
    SigmaF2.formSigma(v.kappa(), v.en(), v.n(), &v);
  }

  std::cout << "\n";

  for (auto &v : wf.valence()) {

    const auto de1 = v * SigmaF2(v);
    const auto de2 = v * SigmaF(v);

    // const auto de0 = MBPT::Sigma_vw(v, v, yk, holes, excited);

    std::cout << v << " " << v.en() << " " << de1 - de2 << "\n";
  }
}
