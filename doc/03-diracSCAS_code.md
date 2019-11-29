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

# Angular

## Angular::Ck_ab (class)

 * Lookup table for C^k and 3j symbols (special m=1/2 case)
 * Builds 3j symbol lookup table for given maximum k and maximum j (2j)
 * 3j symbols, special case: (ja jb k \\ -1/2, 1/2, 0)
 * Ck_ab      := <ka||C^k||kb>        [symmetric up to +/- sign]
 * TildeCk_ab := (-1)^{ja+1/2} * ckab [symmetric]
 * Slightly faster than calculating on-the-fly

```cpp
Ck_ab(const int in_max_K = 0, const int in_max_twoj = 0);
void fill_maxK_twojmax(const int in_max_K, const int in_max_twoj);
```
 * Constructor: fills 3j table with all non-zero 3j-symbols up to and including the maximum specified k and j values (entered as integer 2j)
 * 'fill' function will extend the lookup table to new maximum values

```cpp
double get_tildeCkab_mutable(int k, int ka, int kb);
double get_Ckab_mutable(int k, int ka, int kb);
double get_3jkab_mutable(int k, int ka, int kb);
double get_tildeCkab(int k, int ka, int kb) const;
double get_Ckab(int k, int ka, int kb) const;
double get_3jkab(int k, int ka, int kb) const;
```
 * The _getters_ -- lookup values from table + return them
 * Note: input is _kappa_ (not j or 2j) [Ck includes parity]
 * The _mutable_ versions will calculate (+store) the required 3j symbols if user asks for one that doesn't exist (up to new max_2j and max_K)
   * This is 'safer' (in that it won't segfault), but is slower and not thread-safe [calls std::vector::push_back]
 * not-mutable versions are faster, and thread-safe (read-only) - but perform no bounds checking, so will seg-fault if you ask for a symbol that isn't calculated.

## Angular::SixJ (class)

 * Lookup table for 6j symbols: { ja, jb, k \\ jc, jd, l}
 * j's half-integer (called using integer 2j). Integer k and l
 * Much faster than calculating on the fly.
 * Stores 6j symbols up to + including given max_K and max_2j (stores all allowed l)

```cpp
SixJ(int in_max_k, int in_max_twoj);
void fill(int in_max_k, int in_max_twoj);
```
 * Constructor: fills 6j table with all non-zero 6j-symbols up to and including the maximum specified k and j values (entered as integer 2j)
 * 'fill' function will extend the lookup table to new maximum values

```cpp
double get_6j_mutable(int tja, int tjb, int tjc, int tjd, int k, int l);
double get_6j(int tja, int tjb, int tjc, int tjd, int k, int l) const;
```
 * The _getters_ -- lookup values from table + return them
 * Note: input is integer _2j_ (not j or kappa)
 * The _mutable_ versions will calculate (+store) the required 6j symbols if user asks for one that doesn't exist (up to new max_2j and max_K)
   * This is 'safer' (in that it won't segfault), but is slower and not thread-safe [calls std::vector::push_back]
 * not-mutable versions are faster, and thread-safe (read-only) - but perform no bounds checking, so will seg-fault if you ask for a symbol that isn't calculated.

------------
