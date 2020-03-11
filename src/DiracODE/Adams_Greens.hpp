#pragma once
class DiracSpinor;

namespace DiracODE {
namespace Adams {

void GreenSolution(DiracSpinor &phi, const DiracSpinor &phiI,
                   const DiracSpinor &phi0, const double alpha,
                   const DiracSpinor &Sr);

} // namespace Adams
} // namespace DiracODE
