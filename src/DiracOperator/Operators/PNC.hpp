#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"

namespace DiracOperator {

//******************************************************************************
//! @brief Nuclear-spin independent PNC operator (Qw)
/*! @details
\f[
h_{PNCnsi} = -\frac{G_F \, Q_W}{2\sqrt{2}} \, \rho(r) \, \gamma^5.
\f]
\f[
h_{PNCnsi} = \frac{G_F \, Q_W}{2\sqrt{2}} \, \rho(r) \, \gamma_5.
\f]
  - Ouput is in units of (Qw * 1.e-11.) by default. To get (Qw/-N), multiply
  by (-N) [can go into optional 'factor']

XXX Note - check sign? Ambiguous compared to many sources, probably stemming
from writing gamma_5 but meaning gamma^5 ??

Scalar, ME = Radial integral (F'|h|F)
\f[
(F'|h|F) = -\frac{G_F \, Q_W}{2\sqrt{2}}\int (f'g - g'f) rho(r) dr
\f]

Generates rho(r) according to fermi distribution, given c and t [c and t in
FERMI / femptometers].
*/
class PNCnsi final : public ScalarOperator {
public:
  PNCnsi(double c, double t, const Grid &rgrid, double factor = 1.0,
         const std::string &in_units = "iQw*e-11")
      : ScalarOperator(Parity::odd, factor * PhysConst::GFe11 / std::sqrt(8.0),
                       Nuclear::fermiNuclearDensity_tcN(t, c, 1.0, rgrid),
                       {0, -1, +1, 0}, 0, Realness::imaginary),
        m_unit(in_units) {}
  std::string name() const override final { return "pnc-nsi"; }
  std::string units() const override final { return m_unit; }

private:
  const std::string m_unit{"iQw*e-11"};
};

//******************************************************************************
inline std::unique_ptr<DiracOperator::TensorOperator>
generate_pnc(const IO::InputBlock &input, const Wavefunction &wf) {
  using namespace DiracOperator;
  input.check({{{"c", "Half-density radius for Ferni rho(r). Will use default "
                      "for given nucleus"},
                {"t", "skin thickness [2.3]"},
                {"print", "Write details to screen [true]"}}});
  const auto r_rms = Nuclear::find_rrms(wf.Znuc(), wf.Anuc());
  const auto c = input.get("c", Nuclear::c_hdr_formula_rrms_t(r_rms));
  const auto t = input.get("t", Nuclear::default_t);
  if (input.get("print", true))
    std::cout << "pnc: with c=" << c << ", t=" << t << "\n";
  return std::make_unique<PNCnsi>(c, t, *(wf.rgrid), 1.0, "iQwe-11");
}

} // namespace DiracOperator
