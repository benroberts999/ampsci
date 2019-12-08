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
