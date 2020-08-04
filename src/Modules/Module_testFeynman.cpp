#include "Modules/Module_testFeynman.hpp"
#include "IO/UserInput.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <vector>

/*

struct Sigma_params {
  Method method;
  int min_n_core;
  // Following only for Feynman method
  int max_l_excited;
  bool GreenBasis;
  bool PolBasis;
  double real_omega;
  bool screenCoulomb;
};

struct rgrid_params {
  double r0;
  double rmax;
  std::size_t stride;
};

struct wgrid_params {
  // Only used in Feynman method:
  double w0;
  double wmax;
  double ratio;
};

*/

namespace Module {
void testFeynman(const IO::UserInputBlock &input, const Wavefunction &wf) {
  std::cout << "\ntestFeynman:\n";

  input.checkBlock({"real_omega", "max_l", "pol_basis", "screening", "rmin",
                    "rmax", "stride"});

  const double omre = input.get("real_omega", -0.2);
  const auto method = MBPT::Method::Feynman;
  const auto min_n_core = 1;                       // not used?
  const int max_l_excited = input.get("max_l", 4); // up to g
  const bool GreenBasis = false; // input.get("Green_basis", false);
  const bool PolBasis = input.get("pol_basis", true); //?
  const bool screenCoulomb = input.get("screening", false);
  const MBPT::Sigma_params sigp{method,   min_n_core, max_l_excited, GreenBasis,
                                PolBasis, omre,       screenCoulomb};

  const auto rmin = input.get("rmin", 1.0e-4);
  const auto rmax = input.get("rmax", 30.0);
  const auto stride = input.get("stride", std::size_t(4));
  const MBPT::rgrid_params gridp{rmin, rmax, stride};

  const MBPT::FeynmanSigma Sigma(wf.getHF(), wf.basis, sigp, gridp, {}, "NA");

  for (auto kappa : {-1, 1, -2, 2, -3}) {
    std::cout << "\nkappa: " << kappa << "\n";
    auto Gr1 = Sigma.Green_hf(kappa, omre);
    auto Gr2 = Sigma.Green_hf_basis(kappa, omre, 0.0, false);

    std::cout << "      Gr          Basis       Expected\n";
    for (const auto orbs : {&wf.core, &wf.valence}) {
      for (const auto &Fv : *orbs) {
        if (Fv.k != kappa)
          continue;
        auto vGv1 = Fv * Sigma.Sigma_G_Fv(Gr1.get_real(), Fv);
        auto vGv2 = Fv * Sigma.Sigma_G_Fv(Gr2.get_real(), Fv);
        auto expect = 1.0 / (omre - Fv.en);

        auto error1 = std::abs((vGv1 - expect) / expect);
        auto error2 = std::abs((vGv2 - expect) / expect);
        printf("%4s %11.4e %11.4e %11.4e | %.1e %.1e\n",
               Fv.shortSymbol().c_str(), vGv1, vGv2, expect, error1, error2);
      }
    }
  }
} // namespace Module
} // namespace Module
