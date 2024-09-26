#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/NuclearPotentials.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Wavefunction/Wavefunction.hpp"

#include "qip/Vector.hpp" // qip::overloads

namespace DiracOperator {

//! @brief Field shift operator, (e.g.) dV = V(r+dr) - V(r)
/*! @details 
*/
class fieldshift final : public ScalarOperator {
public:
  fieldshift(const std::vector<double> &vec) : ScalarOperator(vec) {};

  std::string name() const override { return std::string("Field shift"); }
  std::string units() const override { return "au"; }
};

inline std::unique_ptr<DiracOperator::TensorOperator>
generate_fieldshift(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"", "Field shift operator."},
       {"drrms", "Change in nuclear rms charge radius (fm) [0.01]"},
       {"scale_factor", "Scale factor [1.0]"},
       {"", "The following are normally not set:"},
       {"dt", {"Change in skin-thickness t (fm) [0.0]"}},
       {"dbeta", {"Change in quadrupole deformation parameter, beta [0.0]"}}});
  if (input.has_option("help")) {
    return nullptr;
  }

  const auto drrms = input.get("drrms", 0.01);
  const auto scale = input.get("scale_factor", 1.0);
  const auto dt = input.get("dt", 0.0);
  const auto dbeta = input.get("dbeta", 0.0);

  using namespace qip::overloads;

  const auto nuc0 = wf.nucleus();
  const auto vnuc0 = wf.vnuc();

  auto nuc1 = nuc0;
  nuc1.set_rrms(nuc1.r_rms() + drrms);
  nuc1.t() = nuc0.t() + dt;
  nuc1.beta() = nuc0.beta() + dbeta;
  const auto vnuc1 = Nuclear::formPotential(nuc1, wf.grid().r());

  const auto rho0 = Nuclear::fermiNuclearDensity_tcN(
      Nuclear::deformation_effective_t(nuc0.c(), nuc0.t(), nuc0.beta()),
      nuc0.c(), nuc0.z(), wf.grid());
  const auto rho1 = Nuclear::fermiNuclearDensity_tcN(
      Nuclear::deformation_effective_t(nuc1.c(), nuc1.t(), nuc1.beta()),
      nuc1.c(), nuc1.z(), wf.grid());

  constexpr auto abfm4factor = 4 * M_PI * PhysConst::aB_fm * PhysConst::aB_fm *
                               PhysConst::aB_fm * PhysConst::aB_fm;

  /* 
   * <r^4> = \int_0^\infty [\rho(r) r^4] 4 \pi r^2 dr / \int_0^\infty [\rho(r) r^2] 4 \pi r^2 dr,
   * and \int_0^\infty [\rho(r) r^2] 4 \pi r^2 dr = Z by normalisation nuc. charge density.
   * Require four factors of aB to convert output to femtometres.
   */
  const auto rpow6 = wf.grid().rpow(6);

  auto r4_0 = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rho0,
                                 wf.grid().drdu()) *
              abfm4factor / nuc0.z();
  auto r4_1 = NumCalc::integrate(wf.grid().du(), 0.0, 0.0, rpow6, rho1,
                                 wf.grid().drdu()) *
              abfm4factor / nuc0.z();

  const auto delta_r2 =
      nuc1.r_rms() * nuc1.r_rms() - nuc0.r_rms() * nuc0.r_rms();
  const auto delta_r4 = r4_1 - r4_0;
  const auto delta_vnuc = vnuc1 - vnuc0;

  std::cout << "Field shift:\ndelta<r^2> = " << delta_r2
            << " fm^2\ndelta<r^4> = " << delta_r4 << " fm^4\n";

  return std::make_unique<fieldshift>(scale * delta_vnuc);
}

} // namespace DiracOperator
