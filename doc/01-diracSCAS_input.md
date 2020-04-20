# Input options for: diracSCAS

* The **diracSCAS** program should be run as:
  * _./diracSCAS inputFile.in_
  * ``inputFile.in'' is a plain-text input file, that contains all input options (if no input file is given, program looks for the default one, named "diracSCAS.in")
* First, the program reads in the main input options from the four main input "blocks" (Atom, Nucleus, HartreeFock, and Grid). It will use these to generate wavefunction/Hartree-Fock etc. Then, any number of optional "modules" are run using the above-calculated wavefunctions (e.g., these might calculate matrix elements, run tests etc.). The input blocks and options can be given in any order
* In general, the input file will have the following format:

```cpp
Atom { <input_options> }
Nucleus { <input_options> }
Grid { <input_options> }
HartreeFock { <input_options> }
//Optional modules:
Module::firstModule { <input_options> }
Module::secondModule { <input_options> }
```

* Most options have a default, and can be left blank, explicitly set to 'default', or removed entirely.
* The curly-braces denote the start/end of each block. *Don't* use any other curly-braces (nested braces are not allowed)
* Uses c++ style line comments (not block comments). Any commented-out line will not be read. White-space is ignored.
* For example, the following four inputs are all equivalent

```cpp
Atom {
  Z = Cs;
  A;
}
Atom {
  Z = Cs;
//  A;
}
Atom { Z = Cs; }
Atom{Z=Cs;A=default;}
```

* All available inputs for each input block are listed below
  * Inputs are taken in as either text, boolean (true/false), integers, or real numbers.
  * These will be denoted by [t], [b], [i], [r]

********************************************************************************
# Each input block:

## Atom

```cpp
Atom {
  Z = Cs;    //[t/i] required
  A;         //[i] Will look-up default value if default
  varAlpha2; //[r] default = 1
}
```
* Z: Which atom. Can enter as symbol (e.g., Z = Cs;) or integer (e.g., Z = 55;). Required (no default)
* A: nuclear mass number. Leave blank to look up default value
* varAlpha2: Scaling factor for alpha^2 (c = 1/alpha in atomic units): alpha^2 = varAlpha2 * alpha_real^2.
  * default=1. put a very small number to achieve non-relativistic limit (useful for tests)


## HartreeFock
```cpp
HartreeFock {
  core = [Xe]; //[t] required
  valence;     //[t] default = none
  sortOutput;  //[b] default = true
  method;      //[t] default = HartreeFock
  convergence; //[r] default = 1.0e-12
  orthonormaliseValence; //[b] default = false
}
```
* core: Core configuration. Required (no default)
  * Format: [Noble gas],extra (comma separated, no spaces)
  * Can also enter entire configuration, e.g., 1s2,2s2,2p6,.... (As well as Noble gas, can use Zn,Cd,Hg,Cn)
  * Can also add negative values for occupations
  * E.g. :
    * Cs (V^N-1): '[Xe]'
    * Au (V^N-1): '[Xe],4f14,5d10' or '[Hg],6s-2'
    * Tl (V^N-1): '[Xe],4f14,5d10,6s2' or '[Hg]'
    * I (V^N): '[Cd],5p5' or '[Xe],5p-1'
    * H-like: enter as: []  (or 1s0) -- no electrons in core
* valence: which valence states to calculate
  * e.g., "7sp5df" will do s and p states up to n=7, and d and f up to n=5
* sortOutput: true or false. Sort output by energy.
* method: which method to use. can be:
  * HartreeFock(default), ApproxHF, Hartree
* convergence: level we try to converge to.
* orthonormaliseValence: true/false. Orthogonalise valence states? false by default. Only really for testing


## Nucleus
```cpp
Nucleus {
  type;    //[t] default = Fermi
  rrms;    //[r] will loop-up default value based on Z,A
  skin_t;  //[r] default = 2.3
}
```
* rrms: nuclear root-mean-square charge radius (in femptometres = 10^-15m)
* type: Which distribution to use for nucleus? Options are: Fermi (default), spherical, point
* skin_t: skin thickness [only used by Fermi distro]


## Grid
```cpp
Grid {
  r0;         //[r] default = 1.0e-6
  rmax;       //[r] default = 120.0
  num_points; //[i] default = 1600
  type;       //[t] default = loglinear
  b;          //[r] default = rmax/3
  fixed_du;   //[r] default = -1. fixed_du>0: calculate + override num_points
}
```
* r0: grid starting point (in atomic units)
* rmax: Final grid point (in atomic units)
* num_points: number of points in the grid
* type: options are: loglinear (default), logarithmic, linear
  * Note: 'linear' grid requires a _very_ large number of points to work, and should essentially never be used.
* b: only used for loglinear grid; the grid is roughly logarithmic below this value, and linear after it. Default is 4.0 (atomic units). If b<0 or b>rmax, will revert to using a logarithmic grid


## dVpol (effective polarisation potential)
```cpp
dVpol {
  a_eff;  //[r] default = 0.0, typical ~1
  r_cut;  //[r] default = 1.0
}
//nb: all of these are optional, hence entire block can be omitted
```
* Effective polarisation potential:
* dV = -0.5 * a_eff / (r^4 + r_cut^4)
* nb: Added to direct potential _after_ HF for core, but _before_ HF for valence


## ExtraPotential (read from text file)
```cpp
ExtraPotential {
  filename; //[t] default = ""
  factor;   //[r] default = 0.0
  beforeHF; //[b] default = false
}
//nb: all of these are optional, hence entire block can be omitted
```
* Reads in extra potential from text file (space separated: 'x y' format):
* Interpolated these points onto the grid (but does NOT extrapolate,
  potential is assumed to be zero outside the given range)
* Potential is multiplied by 'factor'
* May be added before or after HF (if before: added to vnuc, if after: to vdir)


## RadPot (Ginges/Flambaum QED Radiative Potential)
```cpp
RadPot {
  RadPot;   //[b] default = false, to include QED, set to true
  Simple;   //[r] default = 0.0 // Vrad = -Z^2 * alpha * exp(-r/alpha)
  Ueh;      //[r] default = 0.0 // Uehling (vac pol)
  SE_h;     //[r] default = 0.0 // high-f SE
  SE_l;     //[r] default = 0.0 // low-f SE
  SE_m;     //[r] default = 0.0 // Magnetic SE
  rcut;     //[r] default = 1.0
  scale_rN; //[r] default = 1.0
  scale_l;  //[r,r...] (List) default = 1.0
  core_qed; //[b] default = true
}
```
* Adds QED radiative potential to Hamiltonian.
* RadPot must be set to true to include QED
* Each factor is a scale; 0 means don't include. 1 means include full potential. Any positive number is valid.
* rcut: Only calculates potential for r < rcut [for speed; rcut in au]
* scale_rN: finite nucleus effects: rN = rN * scale_rN (for testing only)
* scale_l: Optional input: Scaling factors for the V_rad for each l state; for higher states, uses the last given input. Input as a list of real numbers. Best explained with examples:
    * scale_l = 1; // include QED for all states
    * scale_l = 0,1,0; //include QED for p states only
    * scale_l = 0,1; //inlcude QED for p,d,f.. but not s states.
    * don't need to be 1 or 0, can be any real number.
* core_qed: if true, will include QED effects into core in Hartree-Fock (relaxation). If false, will include QED only for valence states

## Basis (B-spline basis for MBPT)
* The 'basis' is used for summing over states in MBPT. (A second 'basis', called spectrum, may be used for summation over states in other problems)
```cpp
Basis {
  number;   //[i] default = 0
  order;    //[i] default = 0
  r0;       //[r] default = 0
  r0_eps;   //[r] default = 0
  rmax;     //[r] default = 0
  print;    //[b] default = false
  positron; //[b] default = false
  states;   //[t] default = ""
}
```
* Constructs basis using _number_ splines of order _order_
* on sub-grid (r0,rmax) [if zero, will use full grid]
* r0_eps: Only calculate splines for r where relative core density is larger than r0_eps (updates r0 for each l). Typically ~1.0e-8. Set to zero to use r0.
* If print = true, will print basis energies
* positron: include negative energy states into basis
* states: which basis states to store
  * e.g., "7sp5df" will store s and p states up to n=7, and d and f up to n=5
  * spd will store _all_ (number) states for l<=2


## Correlations (Correlation potential, Sigma)
* For including correlations. 'basis' must exist to calculate Sigma, but not to read Sigma in from file.
```cpp
Correlations {
  io_file;        //[t] default = ""
  stride;         //[i] default = 4
  n_min_core;     //[i] default = 1
  energyShifts;   //[b] default = false
  Brueckner;      //[b] default = false
  lambda_k;       //[r,r...] (list) default is blank.
  fitTo_cm;       //[r,r...] (list) default is blank.
}
```
* Includes correlation corrections. note: splines must exist already
* io_file: Filename (NOT including extension) to read in Correlation potential from. If file doesn't exist, will calculate CP and write to this file. Leave blank to calculate and not write to file.
  * If reading Sigma in from file, basis doesn't need to exist
* stride: Only calculates Sigma every nth point (Sigma is NxN matrix, so stride=4 leads to ~16x speed-up vs 1)
* n_min_core: minimum core n included in the Sigma calculation; lowest states often contribute little, so this speeds up the calculations
* energyShifts: If true, will calculate the second-order energy shifts (from scratch, according to MBPT) - compares to <v|Sigma|v> if it exists
  * Note: Uses basis. If reading Sigma from disk, and no basis given, energy shifts will all be 0.0
* Brueckner: Construct Brueckner valence orbitals using correlation potential method (i.e., include correlations into wavefunctions and energies for valence states)
* lambda_k: Rescale Sigma -> lambda*Sigma. One lambda for each kappa. If not given, assumed to be 1.
  * Note: Lambda's are not written/read to file, so these must be given (if required) even when reading Sigma from disk
* fitTo_cm: Provide list of energies (lowest valence states for each kappa); Sigma for each kappa will be automatically re-scaled to exactly reproduce these. Give as binding energies in inverse cm! It will print the lambda_k's that it calculated
  * e.g., fitTo_cm = -31406.46773, -20228.200, -19674.161, -16907.211, -16809.625; will fit for the lowest s,p,d states for Cs
  * Will over-write lambda_k


## Spectrum (B-spline basis for MBPT)
* The 'Sprectrum' is similar to basis, but also includes correlation corrections (if Sigma exists)
* Useful, since we often need a small basis to compute MBPT terms, but a large basis to complete other sum-over-states calculations.
```cpp
Spectrum {
  // exact same inputs as Basis{}
}
```


## Modules and MatrixElements

Modules and MatrixElements work in essentially the same way. Each MatrixElements/Modules block will be run in order. You can comment-out just the block name, and the block will be skipped.

MatrixElements blocks calculate reduced matrix elements of given operator, Modules can do anything.

For MatrixElements, there are some options that apply for any operator; and then there are some
options specific to each operator

```cpp
MatrixElements::ExampleOperator { //this is not a real operator..
  // Options that apply to all operators:
  printBoth;      //[t] default = false
  onlyDiagonal;   //[t] default = false
  radialIntegral; //[b] default = false
  rpa;            //[b] default = true
  omega;          //[r] default = 0.0, or [t] ('each')
}
```
* printBoth: Print <a|h|b> and <b|h|a> ? false by default. (For _some_ operators, e.g., involving derivatives, this is a good test of numerical error. For most operators, values will be trivially the same; reduced matrix elements, sign may be different.)
* onlyDiagonal: If true, will only print diagonal MEs <a|h|a>
* radialIntegral: calculate radial integral, or reduced matrix elements (difference depends on definition of operator in the code)
* rpa: Include RPA (core polarisation) corrections to MEs, using TDHF method (note: mostly works, but not 100% yet)
* omega: frequency for solving TDHF/RPA equations, should be positive. Put "omega=each;" to solve at the frequency for each transition (i.e., re-solve TDHF for each transition).


### Available operators:

```cpp
MatrixElements::E1 { //Electric dipole operator:
  gauge; //[t] lform, vform. default = lform
}
```

```cpp
MatrixElements::Ek { //Electric multipole operator:
  k; //[i] default = 1
}
```
* k=1 => E1, dipole. k=2 => E2, quadrupole etc.

```cpp
MatrixElements::r { //scalar r
  power; //[r] default = 1. Will calc <|r^n|>.
}
```

```cpp
MatrixElements::pnc {// spin-independent (Qw) PNC operator.
  // Output given in units of i(-Q/N)e-11
  c; //[r] half-density radius. By default, uses rrms from Z,A [see nucleus]
  t; //[r] skin thickness. default = 2.3
}
```
```cpp
MatrixElements::hfs { // Magnetic dipole hyperfine structure constant A
  mu;     //[r] Nuc. mag. moment. Will be looked up by default
  I;      //[r] Nuc. spin. Will be looked up by default
  rrms;   //[r] Nuc. rms radius. Will be looked up by default
  F(r);   //[t] Bohr-Weisskopf function. ball, shell, pointlike, VolotkaBW, doublyOddBW
  printF; //[b] default = false. Will write F(r)/r^2 to text file
  // -----
  // the following are only used for "VolotkaBW/doublyOddBW"
  // both are optional (will be deduced otherwise)
  parity; //[i] parity of unpaired valence nucleon (+1/-1)
  gl;     //[i] =1 for valence proton, =0 for valence neutron
  // -----
  // The following are only read in if F(r) = doublyOddBW,
  // but are _required_ in that case (current values are for 212-Fr)
  mu1 = 4.00; //see paper for explanation
  I1 = 4.5;
  l1 = 5.;
  gl1 = 1;
  I2 = 0.5;
  l2 = 1.;
}
```

### Modules:

-------------------
```cpp
Module::Tests {
  orthonormal;     //[b] Prints worst <a|b>'s. default = true
  orthonormal_all; //[b] Print all <a|b>'s. default = false
  Hamiltonian;     //[b] check eigenvalues of Hamiltonian. default = false
  boundaries;      //[b] check f(rmax)/fmax. default = false
  basisTests;        //[b] Tests basis by evaluating sum rules + HFS/Energies
}
```
* Various tests of numerical errors:

-------------------
```cpp
Module::pnc {
    transition = na, ka, nb, ka; //[i,i,i,i] - required
    rpa;      //[b] default = true
    omega;    //[r]
}
```
* Calculates pnc amplitude {na,ka}->{nb,kb}
* Uses Solving-equations and sum-over-states (spectrum)
* omega: frequency to solve RPA/TDHF equations at. By default is E_a-Eb (from orbitals), but can be anything.

-------------------
```cpp
Module::lifetimes{
  E1;   //[b] Include E1 transitions. default = true
  E1;   //[b] Include E2 transitions. default = false
}
```
Calculates lifetimes of valence states. Note: uses HF energies (prints all data to screen)

-------------------
```cpp
Module::BohrWeisskopf { //Calculates BW effect for Ball/Single-particle
  // Takes same input at MatrixElements::hfs
  // Except for F(r), since it runs for each F(r)
}
```

-------------------
```cpp
Module::WriteOrbitals { //writes orbitals to textfile:
    label = outputLabel; //[t] Optional. blank by default
}
```
Writes the core + valence orbitals (+ the total electron density) to a file, in GNUplot friendly format.
The (optional) label will be appended to the output file name. Plot file using GNUPLOT. For example, try this:
* _plot for [i=2:20] "file.txt" u 1:i every :::0::0  w l t columnheader(i)_
* _plot for [i=2:20] "file.txt" u 1:i every :::1::1  w l t columnheader(i)_
* _plot "file.txt" u 1:2 every :::2::2  w l t "Core Density"_

-------------------
```cpp
Module::FitParametric {
  method = Green;     //[t] Green, Tietz
  statesToFit = core; //[t] core, valence, both
  fitWorst;           //[b] false (default), true;
}
```
Performs a 2D fit to determine the best-fit values for the given two-parameter parametric potential (Green, or Tietz potentials), returns H/g d/t parameters for the best-fit. Does fit to Hartree-Fock energies. Will either do for core or valence states, or both (works best for one or the other). fitWorst: if true, will optimise fit for the worst state. If false, uses least squares for the fit. False is default

-------------------

```cpp
Module::AtomicKernal {
    // Some typical inputs. All are required.
  Emin = 0.01; // in keV
  Emax = 4.0;
  Esteps = 25;
  qmin = 0.001; // in MeV
  qmax = 4.0;
  qsteps = 100;
  max_l_bound = 1; // l for bound states
  max_L = 2;       // L is multipolarity
  output_text = true;
  output_binary = true;
  label = test_new;
  use_plane_waves = false;
}
```
Calculates the "Atomic Kernal" (for scattering/ionisation) for each core
orbital, as a function of momentum transfer (q), and energy deposition (dE).
Writes result to human-readable (and gnuplot-friendly) file, and/or binary.
* For definitions/details, see:
  * B.M. Roberts, V.V. Flambaum [Phys.Rev.D 100, 063017 (2019)](https://link.aps.org/doi/10.1103/PhysRevD.100.063017 "pay-walled"); [arXiv:1904.07127](https://arxiv.org/abs/1904.07127 "free download").
  * B.M.Roberts, V.A.Dzuba, V.V.Flambaum, M.Pospelov, Y.V.Stadnik, [Phys.Rev.D 93, 115037 (2016)](https://link.aps.org/doi/10.1103/PhysRevD.93.115037 "pay-walled"); [arXiv:1604.04559](https://arxiv.org/abs/1604.04559 "free download").
* Note: need quite a dense grid [large number of points] for
  * a) highly oscillating J_L function at low r, and
  * b) to solve equation for high-energy continuum states.
* Sums over 'all' continuum angular momentum states (and multipolarities)
  * Maximum values for l are input parameters
Binary output from this program is read in by dmeXSection program
Note: tested only for neutral atoms (V^N potential).
Also: tested mainly for high values of q
