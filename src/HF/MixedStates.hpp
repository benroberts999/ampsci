#pragma once
#include <vector>
class Wavefunction;
class DiracSpinor;
class Grid;
namespace MBPT {
class CorrelationPotential;
}

namespace HF {
class Breit;

//! @brief Solves Mixed States (Dalgarno-Lewis) equation, inhomogenous
//! equation, with Hartree-Fock hamiltonian, including exchange
/*! @details
Solves
\f[ (H_{\rm HF} - \epsilon - \omega)\delta\phi = -\hat h \phi \f]
for \f$\delta\phi\f$ (dF).
Requires kappa angular momentum number of solution (dF), unperturbed orbital
Fa, a local potential (vl, typically vnuc + vdir), set of core electrons (for
exchange). Note sign on hFa (this is \f$\hat h \phi\f$, not \f$-\hat h
\phi\f$). eps_target is convergance goal for soling the inhomogenous dif.
equation.
*/
DiracSpinor
solveMixedState(const int k, const DiracSpinor &Fa, const double omega,
                const std::vector<double> &vl, const double alpha,
                const std::vector<DiracSpinor> &core, const DiracSpinor &hFa,
                const double eps_target = 1.0e-9,
                const MBPT::CorrelationPotential *const Sigma = nullptr,
                const Breit *const VBr = nullptr,
                const std::vector<double> &H_mag = {});
//! @brief Solves Mixed States (Dalgarno-Lewis equation)
/*! @details
As above, but starts with existing solution dF (may be 'zero'). If existing
solution is already axproximate solution, this allows equation to be solved
much quicker.
*/
void solveMixedState(DiracSpinor &dF, const DiracSpinor &Fa, const double omega,
                     const std::vector<double> &vl, const double alpha,
                     const std::vector<DiracSpinor> &core,
                     const DiracSpinor &hFa, const double eps_target = 1.0e-9,
                     const MBPT::CorrelationPotential *const Sigma = nullptr,
                     const Breit *const VBr = nullptr,
                     const std::vector<double> &H_mag = {});

} // namespace HF
