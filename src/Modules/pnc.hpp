#pragma once
#include <utility>
#include <vector>
class Wavefunction;
namespace IO {
class InputBlock;
class InputBlock;
} // namespace IO
namespace DiracOperator {
class TensorOperator;
}
namespace ExternalField {
class TDHF;
}
class DiracSpinor;
namespace MBPT {
// class CorrelationPotential;
class CorrelationPotential;
} // namespace MBPT
namespace Module {

//! Calculates E1 PNC amplitude
/*! @details

Uses both 'solving equations' (TDHF) and sum-over-states methods.
For solving equations, calculates both:
  - <yA_w|d| B> + <A |d|xB_w>
  - <A |w|XB_d> + <YA_d|w| B>
  - Does not (yet) include DCP

Allowed inputs:

  Module::pnc{ c; t; transition; rpa; omega; nmain; }

 - c, t: half-density radius and skin-thickness (in fm) for rho(r). Will look up
default values by default.
 - transition: For E1_PNC a->b transition.
   - in form "a,b", uses the 'short' notation:
   - e.g., "6s+,7s+" for 6s_1/2 -> 7s_1/2
   - e.g., "6s+,5d-" for 6s_1/2 -> 5d_3/2
 - rpa: true/false. Include RPA or not (TDHF ,method)
 - omega: frequency used for RPA (default is transition frequency of valence).
 - nmain: highest n (prin. q. number) considered as part of 'main'.
   - If not given, will be max(n_core)+4
   - (Calculation broken into core, main, tail)

*/
void calculatePNC(const IO::InputBlock &input, const Wavefunction &wf);

namespace Pnc {

// A->B PNC amplitude, using the sum-over-states method. Note: can swap the
// E1/PNC operator order (for sos this is trivial). main_n is last n considered
// for 'main' term, and all states with en < en_core are assumed to be in the
// core.
std::pair<double, double> pnc_sos(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                  const DiracOperator::TensorOperator *hpnc,
                                  const ExternalField::TDHF *dVpnc,
                                  const DiracOperator::TensorOperator *he1,
                                  const ExternalField::TDHF *dVE1,
                                  const std::vector<DiracSpinor> &spectrum,
                                  int main_n, double en_core, bool print);

// A->B PNC amplitude, using the mixed-states/tdhf method. Note: can swap the
// E1/PNC operator order: this is a good test for this method.
// Note: "spectrum" only used to separate core/main: Final calculation does not
// use any basis states
std::pair<double, double> pnc_tdhf(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                   const DiracOperator::TensorOperator *hpnc,
                                   const ExternalField::TDHF *dVpnc,
                                   const DiracOperator::TensorOperator *he1,
                                   const ExternalField::TDHF *dVE1,
                                   const MBPT::CorrelationPotential *Sigma,
                                   const std::vector<DiracSpinor> &spectrum,
                                   int main_n, double en_core, bool print);

// Force state dF to be orthogonal to core
[[nodiscard]] DiracSpinor
orthog_to_core(DiracSpinor dF, const std::vector<DiracSpinor> &in_orbs,
               double en_core);

// Force state dF to be orthogonal to core+main
[[nodiscard]] DiracSpinor
orthog_to_coremain(DiracSpinor dF, const std::vector<DiracSpinor> &in_orbs,
                   double en_core, int n_main);
} // namespace Pnc

} // namespace Module
