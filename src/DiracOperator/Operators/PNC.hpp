#pragma once
#include "DiracOperator/TensorOperator.hpp"

namespace DiracOperator {

//******************************************************************************
//! @brief Nuclear-spin independent PNC operator (Qw)
/*! @details
\f[
h_{PNCnsi} = \frac{G_F \, Q_W}{2\sqrt{2}} \, \rho(r) \, \gamma_5.
\f]
  - Check -ve sign?
  - Ouput is in units of (Qw * 1.e-11.) by default. To get (Qw/-N), multiply
  by (-N) [can go into optional 'factor']

Generates rho(r) according to fermi distribution, given c and t [c and t in
FERMI / femptometers].
*/
class PNCnsi final : public ScalarOperator {
public:
  PNCnsi(double c, double t, const Grid &rgrid, double factor = 1,
         const std::string &in_units = "iQw*e-11")
      : ScalarOperator(Parity::odd, factor * PhysConst::GFe11 / std::sqrt(8.0),
                       Nuclear::fermiNuclearDensity_tcN(t, c, 1, rgrid),
                       {0, 1, -1, 0}, 0, Realness::imaginary),
        m_unit(in_units) {}
  std::string name() const override final { return "pnc-nsi"; }
  std::string units() const override final { return m_unit; }

private:
  const std::string m_unit{"iQw*e-11"};
};

} // namespace DiracOperator
