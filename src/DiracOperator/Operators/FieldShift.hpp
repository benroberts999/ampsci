#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/PhysConst_constants.hpp"
#include "Potentials/NuclearPotentials.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp" // qip::overloads

namespace DiracOperator {

//! @brief Field shift operator, (e.g.) dV = V(r+dr) - V(r)
/*! @details 
*/
class fieldshift final : public ScalarOperator {

private:
  double m_dr2{-1};
  double m_dr4{-1};

public:
  fieldshift(const Grid &rgrid, const Nuclear::Nucleus &nuc1,
             const Nuclear::Nucleus &nuc2, double scale = 1.0)
      : ScalarOperator(
            Parity::even, scale,
            qip::overloads::operator-(Nuclear::formPotential(nuc2, rgrid.r()),
                                      Nuclear::formPotential(nuc1, rgrid.r()))),
        m_dr2(nuc2.r_rms() * nuc2.r_rms() - nuc1.r_rms() * nuc1.r_rms()),
        m_dr4(r4(rgrid, nuc2) - r4(rgrid, nuc1)) {}

  std::string name() const override { return std::string("Field shift"); }
  std::string units() const override { return "au"; }

  //! \delta <r^2> (in fm^2)
  double dr2() const { return m_dr2; }
  //! \delta <r^4> (in fm^2)
  double dr4() const { return m_dr4; }

  //! Helper function: calculates <r^4> (nuclear) in fm^4
  static double r4(const Grid &grid, const Nuclear::Nucleus &nuc) {

    //  <r^4> = \int_0^\infty [\rho(r) r^4] 4 \pi r^2 dr / Z
    //  Z is from normalisation of rho:
    //  \int_0^\infty [\rho(r) r^2] 4 \pi r^2 dr = Z
    //  by normalisation of nuclear charge density.
    //  Require four factors of aB to convert output to femtometres.

    constexpr auto abfm4factor = 4.0 * M_PI * qip::pow<4>(PhysConst::aB_fm);
    const auto rho = Nuclear::fermiNuclearDensity_tcN(
        Nuclear::deformation_effective_t(nuc.c(), nuc.t(), nuc.beta()), nuc.c(),
        nuc.z(), grid);
    const auto rpow6 = grid.rpow(6);

    return NumCalc::integrate(grid.du(), 0.0, 0.0, rpow6, rho, grid.drdu()) *
           abfm4factor / nuc.z();
  }

  //! Helper function: calculates <r^2> (nuclear) in fm^2
  static double r2(const Grid &grid, const Nuclear::Nucleus &nuc) {

    constexpr auto abfm2factor = 4.0 * M_PI * qip::pow<2>(PhysConst::aB_fm);
    const auto rho = Nuclear::fermiNuclearDensity_tcN(
        Nuclear::deformation_effective_t(nuc.c(), nuc.t(), nuc.beta()), nuc.c(),
        nuc.z(), grid);
    const auto rpow4 = grid.rpow(4);

    return NumCalc::integrate(grid.du(), 0.0, 0.0, rpow4, rho, grid.drdu()) *
           abfm2factor / nuc.z();
  }
};

inline std::unique_ptr<DiracOperator::TensorOperator>
generate_fieldshift(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check(
      {{"", "Field shift operator."},
       {"drrms",
        "Change in nuclear rms charge radius (fm); must not be zero [0.01]"},
       {"scale_factor", "Scale factor [1.0]"},
       {"", "The following are normally not set:"},
       {"dt", "Change in skin-thickness t (fm) [0.0]"},
       {"dbeta", "Change in quadrupole deformation parameter, beta [0.0]"},
       {"print", "Print details? [true]"}});
  if (input.has_option("help")) {
    return nullptr;
  }

  const auto drrms = input.get("drrms", 0.01);
  const auto scale = input.get("scale_factor", 1.0);
  const auto dt = input.get("dt", 0.0);
  const auto dbeta = input.get("dbeta", 0.0);
  const auto print = input.get("print", true);

  const auto &nuc0 = wf.nucleus();

  // Update parameters for "2nd" nucleus:
  auto nuc1 = nuc0;
  nuc1.set_rrms(nuc1.r_rms() + drrms);
  nuc1.t() = nuc0.t() + dt;
  nuc1.beta() = nuc0.beta() + dbeta;
  const auto vnuc1 = Nuclear::formPotential(nuc1, wf.grid().r());

  auto h = std::make_unique<fieldshift>(wf.grid(), nuc0, nuc1, scale);

  if (print) {
    std::cout << "Field shift:\n"
              << "Nuc1: " << nuc0 << "\n"
              << "Nuc2: " << nuc1 << "\n"
              << "delta<r^2> = " << h->dr2() << " fm^2\n"
              << "delta<r^4> = " << h->dr4() << " fm^4\n";
  }

  return h;
}

} // namespace DiracOperator
