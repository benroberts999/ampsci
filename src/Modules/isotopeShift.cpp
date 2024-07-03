#include "Modules/isotopeShift.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "ExternalField/TDHFbasis.hpp"
#include "IO/InputBlock.hpp"
#include "fmt/ostream.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"

#include <iostream>

namespace Module {

void IsotopeShift(const IO::InputBlock &input, const Wavefunction &wf) {
  input.check({{"num_iter", "Number of iterations to compute changes for [10]"},
                {"include_QFS", "Include quadratic field shift (requires spectrum!) [false]"},
                {"dr", "Change in rms charge radius [0.001]"}});

  // If we are just requesting 'help', don't run module:
  if (input.has_option("help")) {
    return;
  }

  const auto num_iter = input.get<int>("num_iter", 10);
  const auto include_qfs = input.get<bool>("include_QFS", false);
  const auto dr = input.get<double>("dr", 0.001);

  using namespace qip::overloads;

  constexpr auto abfm4factor = 4 * M_PI * PhysConst::aB_fm * PhysConst::aB_fm * PhysConst::aB_fm * PhysConst::aB_fm;

  const auto vnuc0 = wf.vnuc();

  // Used for computing fourth radial moment integral.
  auto rpow6 = wf.grid().rpow(6);

  auto nucpt = wf.nucleus();
  const auto r_rms0 = nucpt.r_rms();
  const auto rho0 = Nuclear::fermiNuclearDensity_tcN(
                      Nuclear::deformation_effective_t(nucpt.c(), nucpt.t(), nucpt.beta()),
                      nucpt.c(), nucpt.z(), wf.grid());
  
  const auto r4   = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rho0,  wf.grid().drdu()) * abfm4factor / nucpt.z();

  std::cout << "Calculating change in valence state energies using TDHF(basis).\n"
            << "Field isotopic shift: dE = F d<r^2> + G^2 [d<r^2>]^2 + G^4 d<r^4>\n";

  if(include_qfs && wf.spectrum().empty())
    std::cout << "Attempted to calculate second-order corrections (QFS) w/o spectrum. Skipping\n";

  std::stringstream os;

  for(int i = 1; i < num_iter+1; ++i) {

    nucpt.set_rrms(r_rms0 + i * dr);
    auto vnucpt = Nuclear::formPotential(nucpt, wf.grid().r());

    const auto rhopt = Nuclear::fermiNuclearDensity_tcN(
                        Nuclear::deformation_effective_t(nucpt.c(), nucpt.t(), nucpt.beta()),
                        nucpt.c(), nucpt.z(), wf.grid());

    auto r4pt = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rhopt, wf.grid().drdu()) * abfm4factor/ nucpt.z();

    const auto delta_r2 = nucpt.r_rms() * nucpt.r_rms() - r_rms0 * r_rms0;
    const auto delta_r4 = r4pt - r4;
    const auto delta_vnuc = vnucpt - vnuc0;

    DiracOperator::fieldshift fis(delta_vnuc, delta_r2, delta_r4);

    auto dVfis = ExternalField::TDHFbasis(&fis, wf.vHF(), wf.basis());

    dVfis.solve_core(0.0, 100, false);

    fmt::print(os, "{:9.7f} {:9.7f}", delta_r2, delta_r4);

    for(const auto& Fv : wf.valence()) {
      auto redME = fis.reducedME(Fv, Fv);
      auto dv = dVfis.dV(Fv, Fv);

      fmt::print(os, " {:13.6e}",  (redME + dv)*PhysConst::Hartree_GHz / sqrt(Fv.twojp1()));

      if(include_qfs && !wf.spectrum().empty()) {
        // Both lhs and rhs are "reduced"---need to divide by 3j symbol twice. 
        auto lhs = dVfis.form_dPsi(Fv, 0.0, ExternalField::dPsiType::X, Fv.kappa(), wf.spectrum(), ExternalField::StateType::bra);
        auto rhs = fis.reduced_rhs(Fv.kappa(), Fv) + dVfis.dV_rhs(Fv.kappa(), Fv);
        
        // See, e.g., Eq. 8 Phys. Rev. A 103, (2021)
        auto red_q2 = (lhs * rhs) / (delta_r2 * delta_r2);
        fmt::print(os, " {:13.6e}",  red_q2*PhysConst::Hartree_GHz / Fv.twojp1());
      }
    }
    os << "\n";
  }
  std::cout << "\ndr2 (fm2) dr4 (fm4) ";
  for(auto Fv : wf.valence()){
    fmt::print(stdout, " {:13s}", Fv.shortSymbol()+std::string(" dE (GHz)"));
    if (include_qfs && !wf.spectrum().empty())
      fmt::print(stdout, " {:13s}", "G^2(GHz/fm4)");
  }
  std::cout << "\n";
  std::cout << os.str();
}

} // namespace Module::IsotopeShift