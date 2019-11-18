# Code documentation:

(very unfinished..)

# Dirac

## class DiracSpinor

 * psi1*psi2 = RadialIntegral(psi1,psi2)

---

# Adams

## DiracODE::

 * All free functions in namespace _DiracODE::_
 * Declared in DiracODE.hpp (function definitions in Adams_bound.cpp, Adams_continuum.cpp, Adams_Greens.hpp); some constant coefficients defined in Adams_coefs.hpp

```cpp
void boundState(DiracSpinor &psi, const double en0,
     const std::vector<double> &v, const double alpha,
     int log_dele = 14);
```
 * Solves bound state (H_0 + v - en)psi = 0, for psi and en, with local potential v=v(r) [en stored in psi]
 * en0: initial energy guess (must be close enough to true value, or may fail)
 * v: local potential (electron + nuclear)
 * alpha: alpha = lambda*alpha_0 - effective fine-structure-constant
 * log_eps: log10(eps). eps is convergence target for energy (|en_i - en_{i-1}|/en).


```cpp
void solveContinuum(DiracSpinor &phi, const double en,
     const std::vector<double> &v, const Grid &ext_grid,
     const double r_asym0, const double alpha);
```
 * Solves continuum state (regular at origin) for constant energy en>0, and local potential v=v(r) [electron+nuclear]
 * ext_grid: Grid, must extend past r_asym0
 * r_asym0: initial guess for maximum asymptotic region (real asymptotic r must be _smaller_ than r_asym0)
 * nb: phi is stores on a smaller grid than ext_grid


```cpp
DiracSpinor solve_inhomog(const int kappa, const double en,
     const std::vector<double> &v, const double alpha,
     const DiracSpinor &source);

void solve_inhomog(DiracSpinor &phi, const double en,
     const std::vector<double> &v, const double alpha,
     const DiracSpinor &source);

void solve_inhomog(DiracSpinor &phi, DiracSpinor &phi0, 
     DiracSpinor &phiI, const double en,
     const std::vector<double> &v, const double alpha,
     const DiracSpinor &source);
```

* Functions to solve in-homogeneous Dirac equation: (H-en)Phi = Source, where H (and v) contain only local terms [see "Method" documentation], for constant energy en
* Three overloads do the exact same things, but some are optimised to avoid re-allocating DiracSpinors if the already exist
* These routines solve also for phi0, phiI, which are solutions to homogeneous equation  (H-en)Phi = 0 [reg @ origin, and infinity, respectively]. 
  * The first two throw these solutions away, the third keeps them (in some cases they can be re-used)
  * These Spinors are solved internally and over-written, they don't need to be solved first (i.e., they are out parameters, not in/out parameters)
---

