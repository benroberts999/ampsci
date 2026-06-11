#pragma once
#include "qip/StrongType.hpp"

/*!
  @brief Exact relativistic hydrogen-like (Coulomb) wavefunctions
  @details In the form
  \f[ \psi_{n\kappa m}(\vb{r}) = \frac{1}{r}
      \begin{pmatrix} f_{n\kappa}(r)\,\Omega_{\kappa m}(\hat n) \\
      i\,g_{n\kappa}(r)\,\Omega_{-\kappa,m}(\hat n) \end{pmatrix}. \f]

  From: H. A. Bethe and E. E. Salpeter, Quantum Mechanics of One-and Two-Electron
  Atoms (Plenum, New York, 1977).

  Note: Uses some numerically unstable functions, including Gamma functions and
  confluent hypergeometric functions. So, for some inputs, may be numerically
  unstable. For reasonable inputs (i.e., Zeff=5, up to n=~10), good to at least
  parts in 10^12

  The optional electron mass parameter @p m defaults to 1 (atomic units).
  The full relativistic energy is E = m*c^2 + enk = m/alpha^2 + enk.
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

//! Energy, without rest mass. @p m is the electron mass (default 1 a.u.)
double enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m = 1.0);

//! Full energy: Enk = m/alpha^2 + enk. @p m is the electron mass (default 1 a.u.)
double Enk(PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m = 1.0);

//! Upper (large) radial component. @p m is the electron mass (default 1 a.u.)
double f(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m = 1.0);

//! Lower (small) radial component. @p m is the electron mass (default 1 a.u.)
double g(RaB r, PrincipalQN n, DiracQN k, Zeff z, AlphaFS a, double m = 1.0);

//! Ratio g/f at r, for given energy @p e (without rest mass) and mass @p m.
double gfratio(double r, int k, double z, double a, double e, double m = 1.0);

} // namespace DiracHydrogen
