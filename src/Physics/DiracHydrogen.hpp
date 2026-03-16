#pragma once
#include "qip/StrongType.hpp"

//! @brief Exact relativistic hydrogen-like (Coulomb) wavefuntions
/*!
@details in form phi = (1/r)(f O, ig O')

From: H. A. Bethe and E. E. Salpeter, Quantum Mechanics of One-and Two-Electron
Atoms (Plenum, New York, 1977).

Note: Uses some numerically unstable functions, including Gamma functions and
confluent hypergeometric functions. So, for some inputs, may be numerically
unstable. For reasonable inputs (i.e., Zeff=5, up to n=~10), good to at least
parts in 10^12

*/
namespace DiracHydrogen {

//------------------------------------------------------------------------------
// Uses Strong Types:
enum class DiracTypes { DiracQN, AlphaFS, Zeff, PrincipalQN, RaB };

using DiracQN = qip::StrongType<DiracTypes::DiracQN, int>;
using AlphaFS = qip::StrongType<DiracTypes::AlphaFS, double>;
using Zeff = qip::StrongType<DiracTypes::Zeff, double>;
// double (allow eff):
using PrincipalQN = qip::StrongType<DiracTypes::PrincipalQN, double>;
using RaB = qip::StrongType<DiracTypes::RaB, double>;

//------------------------------------------------------------------------------

//! Relativistic factor gamma = Sqrt[k^2 - (aZ)^2]
double gamma(DiracQN k, Zeff z, AlphaFS a);

//! Energy, without rest mass
double enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a);

//! Enk = enk + c^2
double Enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a);

//! Upper radial component
double f(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a);

//! Lower (small) radial component
double g(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a);

} // namespace DiracHydrogen
