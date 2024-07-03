#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Physics/NuclearPotentials.hpp"

#include "qip/Vector.hpp" // qip::overloads

namespace DiracOperator {

class fieldshift : public ScalarOperator {
public:
  fieldshift(const std::vector<double>& vec, double delta_r2, double delta_r4)
      : ScalarOperator(vec), m_delta_r2(delta_r2), m_delta_r4(delta_r4) {};

  std::string name() const override { return std::string("Field shift"); }
  std::string units() const override { return "au"; }

  double m_delta_r2;
  double m_delta_r4;
};

inline std::unique_ptr<DiracOperator::TensorOperator>
generate_fieldshift(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{"delta_rrms", "Change in nuclear rms (fm) [0.01]"},
              {"scale_factor", "Scale factor [1.0]"}});

  if (input.has_option("help")) {
    return nullptr;
  }

  const auto dr_rms = input.get("delta_rrms", 0.01);
  const auto scale = input.get("scale_factor", 1.0);

  using namespace qip::overloads;

  auto vnuc = wf.vnuc();
  auto nucpt = wf.nucleus();
  const auto r_rms0 = nucpt.r_rms();

  const auto rho0 = Nuclear::fermiNuclearDensity_tcN(Nuclear::deformation_effective_t(nucpt.c(), nucpt.t(), nucpt.beta()), nucpt.c(), nucpt.z(), wf.grid());

  nucpt.set_rrms(r_rms0 + dr_rms);
  auto vnucpt = Nuclear::formPotential(nucpt, wf.grid().r());

  const auto rhopt = Nuclear::fermiNuclearDensity_tcN(Nuclear::deformation_effective_t(nucpt.c(), nucpt.t(), nucpt.beta()), nucpt.c(), nucpt.z(), wf.grid());

  constexpr auto abfm4factor = 4 * M_PI * PhysConst::aB_fm * PhysConst::aB_fm * PhysConst::aB_fm * PhysConst::aB_fm;

  auto rpow6 = wf.grid().rpow(6);

  auto r4   = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rho0,  wf.grid().drdu()) * abfm4factor / nucpt.z();
  auto r4pt = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rhopt, wf.grid().drdu()) * abfm4factor/ nucpt.z();

  const auto delta_r2 = nucpt.r_rms() * nucpt.r_rms() - r_rms0 * r_rms0;
  const auto delta_r4 = r4pt - r4;
  const auto delta_vnuc = vnucpt - vnuc;

  std::cout << "Field shift:\ndelta<r^2> = " << delta_r2 << " fm^2\ndelta<r^4> = " << delta_r4 <<" fm^4\n";

  return std::make_unique<fieldshift>(scale * delta_vnuc, delta_r2, delta_r4);
}

}