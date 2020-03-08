# Code documentation:

(very unfinished..)

------------

# Dirac

## DiracSpinor (class)

 * Class to store radial Dirac spinor: F_nk = (f, g)

```cpp
DiracSpinor(int in_n, int in_k, const Grid &rgrid);
```
 * Constructor. Takes in constant n and k=kappa values + grid
 * A pointer to the Grid is stored (to avoid many copies of Grid)
   * This means Grid must survive [ensured to be true, see WF class]
   * (Should probably move this to be a std::shared_prt ?)

### Public data:
```cpp
const Grid *const p_rgrid; //pointer to radial Grid
const int n;               // principal Q. number
const int k;               // kappa (Dirac Q. number)
double en = 0;             // Single-particle energy (W = en + c^2)
std::vector<double> f;     // Upper radial component
std::vector<double> g;     // Lower (small) radial component
std::size_t pinf;          // f(r) = 0 above r[pinf]. Save time in integrals
int its;                   // # iterations used to converge this orbital
double eps;                // # delta(en)/en: convergence for orbital
double occ_frac;           // Fractional occupation number (=1 for closed shell)
```
 * Note: some important data stored as mutable public data. Little dangerous, but makes interfacing simple+fast. en, f, g etc: would be better if private w/ getters..

### Public const functions:
```cpp
int l() const;       // l
double jjp1() const; // j(j+1)
int twoj() const;
int twojp1() const;
int parity() const;
int k_index() const; // kappa index (-1,1,-2,2,..)=(0,1,2,3,..)
std::string symbol(bool gnuplot = false) const;
std::string shortSymbol() const;
double norm() const;  // returns <Fa|Fa>
std::pair<double, double> r0pinfratio() const;
```
 * Most self explanatory.
 * symbol(): gnuplot=true, puts braces [e.g., 2p_{1/2}]
 * shortSymbol(): e.g., 2p_1/2 -> 2p-, 2p_3/2 -> 2p+
 * r0pinfratio: .first = f(r0)/fmax,  .second = f(pinf-1)/fmax

### Public non-const functions:
```cpp
void scale(const double factor);
void scale(const std::vector<double> &v);
void normalise(double norm_to = 1.0);
```
 * These do what you would expect.
 * By default, normalise normalises to 1, but can normalise to other number.

### operator overloads:
 * v * psi_a, where v is a double or vector works the obvious way
 * psi_a * psi_b = <Fa|Fb>
 * psi_a == psi_b returns true if {na,ka}=={nb,kb}
 * psia > psib : first compares n, and then kappa (via kappa_index)
 * psi_a + psi_b -- **Note**: careful w/ pinf. At the moment, uses pinf of LHS. This might not be what should happen, but needs some work to fix this. (Used when damping orbitals, for example)

---

# Adams

 * NB: H_mag is off-diagonal addition to Dirac Hamiltonian (to include magnetic part of QED corrections.) Typically, included in functions right after v [documentation not updated for this]

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

------------
