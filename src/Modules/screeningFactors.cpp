#include "Coulomb/CoulombIntegrals.hpp"
#include "IO/InputBlock.hpp"
#include "MBPT/FeynmanSigma.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Modules/testFeynman.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <iostream>
#include <vector>

namespace Module {

void screeningFactors(const IO::InputBlock &input, const Wavefunction &wf) {

  std::cout << "\nCalculate effective screening parameters:\n";

  input.check({{"n_min_core", "min n to include in core sum"}});

  const auto omre = -std::abs(0.33 * wf.energy_gap());

  const auto stride =
      (wf.rgrid->getIndex(30.0) - wf.rgrid->getIndex(1.0e-3)) / 100;

  double w0 = 0.01;
  double wratio = 1.5;
  const auto n_min_core = input.get("n_min_core", 3);

  const int max_l_excited = 6;

  const auto include_G = false;

  const MBPT::Sigma_params sigp_scr{MBPT::Method::Feynman,
                                    n_min_core,
                                    include_G,
                                    max_l_excited,
                                    false,
                                    false,
                                    omre,
                                    w0,
                                    wratio,
                                    true,  // screening
                                    false, // hole-particle
                                    {}};

  const MBPT::Sigma_params sigp_hp{MBPT::Method::Feynman,
                                   n_min_core,
                                   include_G,
                                   max_l_excited,
                                   false,
                                   false,
                                   omre,
                                   w0,
                                   wratio,
                                   false, // screening
                                   true,  // hole-particle
                                   {}};

  const MBPT::Sigma_params sigp_AO{MBPT::Method::Feynman,
                                   n_min_core,
                                   include_G,
                                   max_l_excited,
                                   false,
                                   false,
                                   omre,
                                   w0,
                                   wratio,
                                   true, // screening
                                   true, // hole-particle
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
                                  false, // screening
                                  false, // hole-particle
                                  {}};

  const MBPT::rgrid_params gridp{1.0e-3, 30.0, stride};

  const MBPT::FeynmanSigma Sigma0(wf.getHF(), wf.basis, sigp_0, gridp, "");
  const MBPT::FeynmanSigma SigmaScr(wf.getHF(), wf.basis, sigp_scr, gridp, "");
  const MBPT::FeynmanSigma Sigmahp(wf.getHF(), wf.basis, sigp_hp, gridp, "");
  const MBPT::FeynmanSigma SigmaAO(wf.getHF(), wf.basis, sigp_AO, gridp, "");

  std::cout << "\n";

  const auto unit = PhysConst::Hartree_invcm;

  double weight = 0.0;
  std::vector<double> fk_avg;
  std::vector<double> fk2_avg;
  std::vector<double> eta_avg;
  // std::vector<int> kappa;

  std::cout << "k   de(0)   de(X)  de(hp)  de(AO) |     fk   fkhp   eta\n";
  for (const auto &Fv : wf.valence) {

    std::vector<double> fk;
    std::vector<double> fk_hp;
    std::vector<double> eta_hp;

    // XXX Two issues:
    // 1) For fk including hp, would be better to include holeparticle into
    // Coulomb only, but not into loop, then divide by 2nd order.
    // 2) For eta (hole particle correction): would be better to do for each
    // core state, not each k!
    // Both of these require changing FeynmanSigma class...

    std::cout << Fv.symbol() << "\n";
    // auto &fk_v = fk.emplace_back();
    // kappa.push_back(Fv.kappa());
    for (int k = 0; k < 8; ++k) {
      auto Sd0 = Sigma0.FeynmanDirect(Fv.kappa(), Fv.en(), k);
      auto SdX = SigmaScr.FeynmanDirect(Fv.kappa(), Fv.en(), k);
      auto Sdhp = Sigmahp.FeynmanDirect(Fv.kappa(), Fv.en(), k);
      auto SdAO = SigmaAO.FeynmanDirect(Fv.kappa(), Fv.en(), k);

      const auto de0 = Fv * Sigma0.act_G_Fv(Sd0, Fv);
      const auto deX = Fv * SigmaScr.act_G_Fv(SdX, Fv);
      const auto dehp = Fv * Sigmahp.act_G_Fv(Sdhp, Fv);
      const auto deAO = Fv * SigmaAO.act_G_Fv(SdAO, Fv);

      if (std::abs(de0 * unit) < 0.01)
        break;

      auto t_fk = deX / de0;
      auto t_fk_hp = deAO / dehp;
      auto t_eta_hp = dehp / de0; // nb: do each core state, not each k!

      fk.push_back(t_fk);
      fk_hp.push_back(t_fk_hp);
      eta_hp.push_back(t_eta_hp);

      printf("%i %7.1f %7.1f %7.1f %7.1f | %6.3f %6.3f %6.3f\n", k, de0 * unit,
             deX * unit, dehp * unit, deAO * unit, t_fk, t_fk_hp, t_eta_hp);
    }

    std::cout << Fv.symbol() << "\n";

    const auto printer = [](auto x) { printf("%.3f, ", x); };
    std::cout << "fk = ";
    std::for_each(fk.begin(), fk.end(), printer);
    std::cout << "// no hp\n";
    std::cout << "fk = ";
    std::for_each(fk_hp.begin(), fk_hp.end(), printer);
    std::cout << "// with hp\n";
    std::cout << "eta= ";
    std::for_each(eta_hp.begin(), eta_hp.end(), printer);
    std::cout << "// hp only (no scr.)\n";
    std::cout << "\n";

    using namespace qip::overloads;
    fk_avg += Fv.en() * fk;
    fk2_avg += Fv.en() * fk_hp;
    eta_avg += Fv.en() * eta_hp;
    weight += Fv.en();
  }

  std::cout << "\nWeighted average:\n";

  using namespace qip::overloads;
  fk_avg /= weight;
  fk2_avg /= weight;
  eta_avg /= weight;

  const auto printer = [](auto x) { printf("%.3f, ", x); };
  std::cout << "fk = ";
  std::for_each(fk_avg.begin(), fk_avg.end(), printer);
  std::cout << " // no hp\n";
  std::cout << "fk = ";
  std::for_each(fk2_avg.begin(), fk2_avg.end(), printer);
  std::cout << "// with hp\n";
  std::cout << "eta= ";
  std::for_each(eta_avg.begin(), eta_avg.end(), printer);
  std::cout << "// hp only (no scr.)\n";
}

} // namespace Module
