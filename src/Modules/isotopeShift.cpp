#include "Modules/isotopeShift.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "IO/InputBlock.hpp"
#include "fmt/ostream.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"

#include "gsl/gsl_statistics_double.h"

#include <iostream>

namespace Module {

void isotopeShift(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check(
    {{"num_iter", "Number of iterations to compute changes for [10]"},
    {"QFS", "Include quadratic field shift [false]"},
    {"use_spectrum",
      "If true (and spectrum available), will use spectrum for QFS "
      "states [false]"},
    {"dr", "Change in rms charge radius [0.001]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  IO::ChronoTimer timer("isotopeShift");

  const auto num_iter = input.get<std::size_t>("num_iter", 10);
  const auto include_qfs = input.get<bool>("QFS", false);
  const auto dr = input.get<double>("dr", 0.001);
  const auto use_spectrum =
      wf.spectrum().empty() ? false : input.get("use_spectrum", false);

  using namespace qip::overloads;

  /* 
   * <r^4> = \int_0^\infty [\rho(r) r^4] 4 \pi r^2 dr / \int_0^\infty [\rho(r) r^2] 4 \pi r^2 dr,
   * and \int_0^\infty [\rho(r) r^2] 4 \pi r^2 dr = Z by normalisation nuc. charge density.
   * Require four factors of aB to convert output to femtometres.
   */
  constexpr auto unit_conv = 4.0 * M_PI * PhysConst::aB_fm * PhysConst::aB_fm * PhysConst::aB_fm * PhysConst::aB_fm;
  const auto rfour_factor = unit_conv / wf.nucleus().z();

  // Used for computing fourth radial moment integral.
  const auto rpow6 = wf.grid().rpow(6);

  auto nuc = wf.nucleus();
  auto r_rms = nuc.r_rms();

  const auto r_rms0 = r_rms;
  const auto rsq0 = r_rms0 * r_rms0;
  const auto vnuc0 = wf.vnuc();

  const auto rho0 = Nuclear::fermiNuclearDensity_tcN(
                      Nuclear::deformation_effective_t(nuc.c(), nuc.t(), nuc.beta()),
                      nuc.c(), nuc.z(), wf.grid());
  
  const auto rfour0 = NumCalc::integrate(wf.grid().du(), 0.0, 0.0,
                                    rpow6, rho0, wf.grid().drdu())
                                    * rfour_factor;

  std::cout << "Calculating first-order shift in valence state energies using TDHF(basis)\n"
            << "Field isotopic shift: dE = F d<r^2> + G^4 d<r^4>\n";

  std::cout << "\nInitial nuclear parameters:\nr_rms0 = " << r_rms0 << "fm, "
            << "<r^2>0 = " << rsq0 << "fm^2, "
            << "<r^4>0 = " << rfour0 << "fm^4\n";

  // G^2 factor is second-order PT correction and should be approx. constant regardless of \delta<r^2>.
  auto G2_datatotal = include_qfs ?
      std::vector<std::vector<double>>(wf.valence().size(), std::vector<double>(num_iter)) : std::vector<std::vector<double>>{};

  //std::stringstream os;
  std::cout << "\nr_rms (fm)   dr2 (fm2)   dr4 (fm4)   ";
  for(const auto& Fv : wf.valence()){
    fmt::print(stdout, " {:15s}", Fv.shortSymbol()+std::string(" dE (GHz)"));
  }
  std::cout << std::endl;

  for(std::size_t i = 0ul; i < num_iter; ++i) {
    nuc.set_rrms(r_rms + dr);

    r_rms = nuc.r_rms();

    auto vnucpt = Nuclear::formPotential(nuc, wf.grid().r());

    const auto rhopt = Nuclear::fermiNuclearDensity_tcN(
                        Nuclear::deformation_effective_t(nuc.c(), nuc.t(), nuc.beta()),
                        nuc.c(), nuc.z(), wf.grid());

    auto rfour = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rhopt, wf.grid().drdu()) * rfour_factor;

    const auto delta_r2 = r_rms * r_rms - rsq0;
    const auto delta_r4 = rfour - rfour0;
    const auto recip_delta_r2sq = include_qfs ? 1.0/(delta_r2 * delta_r2) : 0.0;

    const auto delta_vnuc = vnucpt - vnuc0;

    DiracOperator::fieldshift fis(delta_vnuc);

    auto dVfis = ExternalField::TDHFbasis(&fis, wf.vHF(), wf.basis());

    dVfis.solve_core(0.0, 100, false);

    fmt::print(stdout, "{:10.4f} {:11.7f} {:11.7f}", r_rms, delta_r2, delta_r4);

    for(std::size_t j = 0ul; j < wf.valence().size(); ++j) {
      const auto& Fv = wf.valence().at(j);

      auto redME = fis.reducedME(Fv, Fv);
      auto dv = dVfis.dV(Fv, Fv);

      fmt::print(stdout, " {:15.6e}",  (redME + dv)*PhysConst::Hartree_GHz / sqrt(Fv.twojp1()));

      if(include_qfs) {
        const auto& dPsi_basis = use_spectrum ? wf.spectrum() : wf.basis();

        auto lhs = dVfis.form_dPsi(Fv, 0.0, ExternalField::dPsiType::X, Fv.kappa(), dPsi_basis, ExternalField::StateType::bra);
        auto rhs = fis.reduced_rhs(Fv.kappa(), Fv) + dVfis.dV_rhs(Fv.kappa(), Fv);

        // Both lhs and rhs are "reduced"---need to divide by 3j symbol twice. 
        // For physics see, e.g., Eq. 8 Phys. Rev. A 103, (2021)
        G2_datatotal.at(j).at(i) = (lhs * rhs) * recip_delta_r2sq * PhysConst::Hartree_GHz / Fv.twojp1();
      }
    }

    std::cout << std::endl;
  }

  std::cout << "\nNB. dr2 := <r^2> - <r^2>0, dr4 := <r^4> - <r^4>0\n";

  if(include_qfs) {
    std::cout << "\n            G^2 (GHz)  Std. err. (GHz)\n";
    for(std::size_t i = 0ul; i < wf.valence().size(); ++i) {
        auto G2_data = G2_datatotal.at(i).data();
        auto mean = gsl_stats_mean(G2_data, 1, num_iter);
        auto SEM = gsl_stats_sd(G2_data, 1, num_iter) / sqrt(num_iter);
        fmt::print(stdout, "{:7s} {:13.6e}    {:13.6e}\n",
                    wf.valence().at(i).symbol(), mean, SEM);
    }
    std::cout << "\n";
  }
}

} // namespace Module::isotopeShift