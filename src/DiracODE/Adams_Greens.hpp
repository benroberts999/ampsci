#pragma once
class DiracSpinor;

namespace DiracODE {
namespace Adams {

void GreenSolution(DiracSpinor &Fa, const DiracSpinor &Finf,
                   const DiracSpinor &Fzero, const double alpha,
                   const DiracSpinor &Sr);

} // namespace Adams
} // namespace DiracODE
