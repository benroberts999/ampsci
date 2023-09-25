#include "MBPT/Goldstone.hpp"
#include "DiracODE/DiracODE.hpp"
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
TEST_CASE("MBPT: Goldstone", "[MBPT][Goldstone]") { //

  // Wavefunction wf({1000, 1.0e-6, 100.0, 0.33 * 100.0, "loglinear"},
  //                 {"Na", -1, "Fermi"}, 1.0);
  // wf.solve_core("HartreeFock", 0.0, "[Ne]");
  // wf.formBasis({"35spdfg", 40, 7, 1.0e-3, 1.0e-4, 40.0, false});
  // wf.solve_valence("4sp3d");

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
}
