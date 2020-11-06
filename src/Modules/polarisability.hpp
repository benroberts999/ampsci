#pragma once
#include <vector>
class Wavefunction;
namespace IO {
class UserInputBlock;
class UserInput;
} // namespace IO
namespace DiracOperator {
class E1;
}
namespace ExternalField {
class TDHF;
}
namespace MBPT {
class CorrelationPotential;
}
class DiracSpinor;

namespace Module {

//! Calculate dipole polarisabilitities (static, dynamic, alpha, vector beta)
/*! @details

Uses both 'solving equations' (TDHF) and sum-over-states methods.

Allowed inputs:

  Module::polarisability{ rpa; omega; transition; omega_max; omega_steps;  }

 - rpa: true/false. Include RPA or not (TDHF ,method)
 - omega: frequency used for alpha_0 (dipole polarisability). default is 0.
 - transition: For scalar/vector a->b transition polarisability.
   - in form "a,b", e.g., "6s+,7s+" for 6s_1/2 -> 7s_1/2
 - omega_max: maximum frequency for dynamic polarisability. Default is 0.
   - nb: only runs dynamic pol. if omega_max>0
 - omega_steps: Number of steps used for dynamic. default = 30. (linear scale)

Note: transition polarisabilities written for s-states only.
They might be correct for other states too, but NOT checked.
Especially for beta, pretty sure it's wrong for non-s states.

*/
void polarisability(const IO::UserInputBlock &input, const Wavefunction &wf);

//------------------------------------------------------------------------------
// "private" helper functions:
namespace Polarisability {
double alpha_core_tdhf(const std::vector<DiracSpinor> &core,
                       const DiracOperator::E1 &he1, ExternalField::TDHF &dVE1,
                       double omega);

double
alpha_valence_tdhf(const DiracSpinor &Fa, const DiracSpinor &Fb,
                   const DiracOperator::E1 &he1, double omega,
                   ExternalField::TDHF &dVE1,
                   const MBPT::CorrelationPotential *const Sigma = nullptr);

double alpha_core_sos(const std::vector<DiracSpinor> &core,
                      const std::vector<DiracSpinor> &basis,
                      const DiracOperator::E1 &he1, ExternalField::TDHF &dVE1,
                      double omega);

double alpha_valence_sos(const DiracSpinor &Fv,
                         const std::vector<DiracSpinor> &basis,
                         const DiracOperator::E1 &he1,
                         ExternalField::TDHF &dVE1, double omega);

double alpha_valence_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                         const std::vector<DiracSpinor> &basis,
                         const DiracOperator::E1 &he1,
                         ExternalField::TDHF &dVE1, double omega = 0.0);

double alpha_v_SRN(const DiracSpinor &Fv,
                   const std::vector<DiracSpinor> &spectrum, int n_max_sum,
                   const std::vector<DiracSpinor> &hf_basis,
                   const double en_core, const DiracOperator::E1 &he1,
                   ExternalField::TDHF &dVE1, double omega);

std::pair<double, double> beta_sos(const DiracSpinor &Fv, const DiracSpinor &Fw,
                                   const std::vector<DiracSpinor> &basis,
                                   const DiracOperator::E1 &he1,
                                   ExternalField::TDHF &dVE1,
                                   double omega = 0.0);
} // namespace Polarisability

} // namespace Module
