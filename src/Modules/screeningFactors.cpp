#include "Coulomb/CoulombIntegrals.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Modules/testFeynman.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <iostream>
#include <vector>

namespace Module {

void screeningFactors(const IO::InputBlock &input, const Wavefunction &wf) {

  std::cout << "\nCalculate effective screening parameters:\n";

  input.checkBlock_old({"n_min_core"});

  const auto omre = -std::abs(0.33 * wf.energy_gap());

  const auto stride =
      (wf.rgrid->getIndex(30.0) - wf.rgrid->getIndex(1.0e-4)) / 150;

  double w0 = 0.01;
  double wratio = 1.5;
  const auto n_min_core = input.get("n_min_core", 3);

  const int max_l_excited = 6;

  const auto include_G = false;

  const MBPT::Sigma_params sigp_X{MBPT::Method::Feynman,
                                  n_min_core,
                                  include_G,
                                  max_l_excited,
                                  false,
                                  false,
                                  omre,
                                  w0,
                                  wratio,
                                  true,
                                  false,
                                  {}};
  const MBPT::Sigma_params sigp_0{MBPT::Method::Feynman,
                                  n_min_core,
                                  include_G,
                                  max_l_excited,
                                  false,
                                  false,
                                  omre,
                                  w0,
                                  wratio,
                                  false,
                                  false,
                                  {}};

  const MBPT::rgrid_params gridp{1.0e-4, 30.0, stride};

  const MBPT::FeynmanSigma Sigma0(wf.getHF(), wf.basis, sigp_0, gridp, "");
  const MBPT::FeynmanSigma SigmaX(wf.getHF(), wf.basis, sigp_X, gridp, "");

  std::cout << "\n";

  const auto unit = PhysConst::Hartree_invcm;

  std::vector<std::vector<double>> fk;
  std::vector<int> kappa;

  std::cout << "k     de(0)     de(X)     fk\n";
  for (const auto &Fv : wf.valence) {
    std::cout << Fv.symbol() << "\n";
    auto &fk_v = fk.emplace_back();
    kappa.push_back(Fv.k);
    for (int k = 0; k < 9; ++k) {
      auto Sd0 = Sigma0.FeynmanDirect(Fv.k, Fv.en(), k);
      auto SdX = SigmaX.FeynmanDirect(Fv.k, Fv.en(), k);

      const auto de0 = Fv * Sigma0.act_G_Fv(Sd0, Fv);
      const auto deX = Fv * SigmaX.act_G_Fv(SdX, Fv);
      if (std::abs(de0 * unit) < 0.01)
        break;
      auto fkk = deX / de0;

      fk_v.push_back(fkk);

      printf("%i %9.3f %9.3f %6.3f\n", k, de0 * unit, deX * unit, fkk);
    }
  }

  std::vector<double> fk_av(10);
  std::vector<int> fk_w(10);
  for (auto i = 0ul; i < fk.size(); ++i) {
    std::cout << "kappa = " << kappa[i] << "\n";
    std::cout << "fk= ";
    for (auto j = 0ul; j < fk[i].size(); ++j) {
      printf("%.3f, ", fk[i][j]);
      fk_av[j] += fk[i][j];
      fk_w[j] += 1;
    }
    std::cout << "\n";
  }

  std::cout << "\nAverage:\n";
  std::cout << "fk = ";
  for (auto k = 0ul; k < fk_av.size(); ++k) {
    if (fk_w[k] == 0)
      break;
    printf("%.3f, ", fk_av[k] / fk_w[k]);
  }
  std::cout << "\n";
}

} // namespace Module
