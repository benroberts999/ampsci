# Input options for: ampsci

This outlines/describes the input options/usage for ampsci. For a description of the physics, see: [ampsci.pdf](https://benroberts999.github.io/ampsci/ampsci.pdf)

* The **ampsci** program should be run as:
  * _./ampsci inputFile.in_
  * "inputFile.in" is a plain-text input file, that contains all input options * If no input file is given, program looks for the default one, named "ampsci.in"
* Can also be run simply by giving an atomic symbol (or Z) as command-line option, which will run a simple Hartree-Fock calculation, e.g.,: _./ampsci Cs_
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

* Most options have a default; these may be left blank, explicitly set to 'default', or removed entirely.
* The curly-braces denote the start/end of each block. Nested blocks are allowed
* Uses c++ style comments. Any commented-out line will not be read. White-space is ignored.
* For example, the following inputs are all equivalent

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
Atom { Z = 55; A = 133; }
Atom{Z=Cs;A=default;}
```

* All available inputs for each input block are listed below
  * Inputs are taken in as either text, boolean (true/false), integers, or real numbers, or a "sub-block"
  * These will be denoted by [t], [b], [i], [r], [sub-block]
  * A sub-block is a bracketed list of sub-options, e.g.
    * options{a=1; b=2;}
* Program will _usually_ warn you if an incorrect option is given, and print all the available options to the screen -- you can use this fact to get the code to print all the available options for each block.

********************************************************************************
## Auto-documentation

 * This document may go out-of-sync with the code
 * The best fail-safe way to check the available options for a given input block is to use the code itself.
 * Do this by adding a blank 'help' option to the input file:
 * The code will then print a list of all available options, and (usually) an explanation for them

blockname {
  help;
}



********************************************************************************
All available ampsci options/blocks are:

```cpp
ampsci{
  Atom;  // Which atom to run for
  Grid;  // Set radial grid parameters
  HartreeFock;  // Expl
  Nucleus;  // Set nuclear parameters
  RadPot;  // Inlcude QED radiative potential
  Basis;  // Basis used for MBPT
  Spectrum;  // Like basis; used for sum-over-states
  Correlations;  // Options for correlations
  ExtraPotential;  // Include an extra potential
  dVpol;  // Approximate correlation (polarisation) potential
  Module::*;  // Run any number of modules (* -> module name)
}
```


********************************************************************************
# Some details for each input block:

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
  core;        //[t] default = [] (no core)
  valence;     //[t] default = none
  sortOutput;  //[b] default = false
  method;      //[t] default = HartreeFock
  Breit;       //[r] default = 0.0
  convergence; //[r] default = 1.0e-12
}
```
* core: Core configuration. Format: "[Atom],extra"
  * 'Extra' is comma-separated list of 'nLm' terms (eg: '1s2,2s2,2p6')
  * '[Atom]' is optional, can just list all configurations
  * '[Atom]' is usually a Noble gas, but can be any atom
  * Can also add negative values for occupations (m)
  * E.g. :
    * Cs (V^N-1): '[Xe]'
    * Cs (V^N-1): '[Cs],6s-1'  (equivalent to [Xe])
    * Au (V^N-1): '[Xe],4f14,5d10' or '[Hg],6s-2'
    * Tl (V^N-1): '[Xe],4f14,5d10,6s2' or '[Hg]'
    * I (V^N): '[Cd],5p5' or '[Xe],5p-1', or [I]
    * H-like: enter as: []  (or 1s0) -- no electrons in core
* valence: which valence states to calculate
  * e.g., "7sp5df" will do s and p states up to n=7, and d and f up to n=5
* sortOutput: true or false. Sort output by energy.
* method: which method to use. can be:
  * HartreeFock(default), ApproxHF, Hartree, KohnSham
  * Note: for KohnSham: Should be V^N - i.e., include lowest valence state into core.
    * e.g., for Cs: core=[Xe],6s1 (or 'core=[Cs]'); Then list each valence state in valence:
    * e.g., to solve valence 6s, 7s, and 6p, write valence="6s7s6p";
* Breit: Include Breit into HF with given scale (0 means don't include)
  * Note: Will go into spline basis, and RPA equations automatically
* convergence: level we try to converge to.


## Nucleus
```cpp
Nucleus {
  type;    //[t] default = Fermi
  rrms;    //[r] will loop-up default value based on Z,A
  c;       //[r] will loop-up default value based on Z,A
  t;       //[r] default = 2.3
}
```
* rrms: nuclear root-mean-square charge radius (in femptometres = 10^-15m)
* type: Which distribution to use for nucleus? Options are: Fermi (default), spherical, point
* t: skin thickness [only used by Fermi distro]
* c: half-density radius [only used by Fermi distro]
  * nb: if rrms and c are given, c takes priority, and rrms is calculated from c (and t)


## Grid
```cpp
Grid {
  r0;         //[r] default = 1.0e-6
  rmax;       //[r] default = 120.0
  num_points; //[i] default = 1600
  type;       //[t] default = loglinear
  b;          //[r] default = rmax/3
  du;         //[r] default = blank.
}
```
* r0: grid starting point (in atomic units)
* rmax: Final grid point (in atomic units)
* num_points: number of points in the grid
* type: options are: loglinear (default), logarithmic, linear
  * Note: 'linear' grid requires a _very_ large number of points to work, and should essentially never be used.
* b: only used for loglinear grid; the grid is roughly logarithmic below this value, and linear after it. Default is 4.0 (atomic units). If b<0 or b>rmax, will revert to using a logarithmic grid
* du: if du>0.0, it will calculate num_points to fix du (step-size in uniform 'u' grid); will over-ride 'num_points' option.


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
* Interpolates these points onto the grid (but does NOT extrapolate,
  potential is assumed to be zero outside the given range)
* Potential is multiplied by 'factor'
* May be added before or after HF (if before: added to vnuc, if after: to vdir)


## RadPot (Ginges/Flambaum QED Radiative Potential)
```cpp
RadPot {
  Ueh;      //[r] default = 1.0 // Uehling (vac pol)
  SE_h;     //[r] default = 1.0 // high-f SE
  SE_l;     //[r] default = 1.0 // low-f SE
  SE_m;     //[r] default = 1.0 // Magnetic SE
  WK;       //[r] default = 0.0 // Wickman-Kroll
  rcut;     //[r] default = 5.0
  scale_rN; //[r] default = 1.0
  scale_l;  //[r,r...] (List) default = 1.0
  core_qed; //[b] default = true
}
```
* Adds QED radiative potential to Hamiltonian.
* QED will be included if this block is present; else not
* Will read from file if it exists (e.g., Z_uhlmw.qed)
* Each factor (Ueh, SE_h,..) is a scale; 0 means don't include. 1 means include full potential. Any positive number is valid.
* rcut: Only calculates potential for r < rcut [for speed; rcut in au]
* scale_rN: finite nucleus effects: rN = rN * scale_rN (=0 means pointlike)
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
  read;           //[t] default = ""
  write;          //[t] default = ""
  n_min_core;     //[i] default = 1
  energyShifts;   //[b] default = false
  Brueckner;      //[b] default = false
  lambda_kappa;   //[r,r...] (list) default is blank.
  fk;             //[r,r...] (list) default is blank.
  fitTo_cm;       //[r,r...] (list) default is blank.
  // Following are "sub-grid" options:
  stride;         //[i] default chosen so there's ~150 pts in region [e-4,30]
  rmin;           //[i] 1.0e-4
  rmax;           //[i] 30.0
}
```
* Includes correlation corrections. note: splines must exist already
* read/write: Read/write from/to file. Set to 'false' to calculate from scratch (and not write to file). By default, the file name is: "Atom".sig.
  * Alternatively, put any text here to be a custom filename (e.g., read/write="Cs_new"; will read/write from/to Cs_new.sig). Don't include the '.sig' extension (uses sigf for Feynman method, sig2 for Goldstone). Grids must match exactly when reading in from a file.
  * If reading Sigma in from file, basis doesn't need to exist
* n_min_core: minimum core n included in the Sigma calculation; lowest states often contribute little, so this speeds up the calculations
* energyShifts: If true, will calculate the second-order energy shifts (from scratch, according to MBPT) - compares to <v|Sigma|v> if it exists
  * Note: Uses basis. If reading Sigma from disk, and no basis given, energy shifts will all be 0.0
* Brueckner: Construct Brueckner valence orbitals using correlation potential method (i.e., include correlations into wavefunctions and energies for valence states)
* stride: Only calculates Sigma every nth point (Sigma is NxN matrix, so stride=4 leads to ~16x speed-up vs 1)
* rmin/rmax: min/max points along radial Grid Sigma is calculated+stored.
* lambda_kappa: Rescale Sigma -> lambda*Sigma. One lambda for each kappa. If not given, assumed to be 1.
  * Note: Lambda's are not written/read to file, so these must be given (if required) even when reading Sigma from disk
* fk: Effective screening factors; only used for 2nd-order Goldstone method
  * Note: Included directly into Sigma
  * e.g., for Cs: fk = 0.72, 0.62, 0.83, 0.89, 0.94, 1.0;
* fitTo_cm: Provide list of energies (lowest valence states for each kappa); Sigma for each kappa will be automatically re-scaled to exactly reproduce these. Give as binding energies in inverse cm! It will print the lambda_kappa's that it calculated
  * e.g., fitTo_cm = -31406.5, -20228.2, -19674.1; will fit for the lowest s & p states for Cs
  * Will over-write lambda_kappa
  * -43487.11, -28583.45, -28583.11, -12204.03, -12203.99, -6856.91, -6856.91; // Li
  * -41449.45, -24493.28, -24476.08, -12276.56, -12276.61, -6862.53, -6862.53; // Na
  * -35009.81, -22024.63, -21966.92, -13472.83, -13475.13, -6881.96, -6881.96; // K
  * -33690.81, -21111.86, -20874.265, -14335.161, -14335.607, -6898.692, -6898.718; // Rb
  * -31406.468, -20228.200, -19674.161, -16907.211, -16809.625, -6934.241, -6934.422; // Cs
  * -80686.30, -60424.74, -58733.90, -75812.45, -75011.49, -32427.68, -32202.97; // Ba+
  * -32848.87, -20611.46, -18924.87, -16619.00, -16419.23; // Fr
  * -81842.5 -60491.2, -55633.6, -69758.2, -68099.5, -32854.6, -32570.4; // Ra+


## Spectrum (B-spline basis for MBPT)
* The 'Spectrum' is similar to basis, but also includes correlation corrections (if Sigma exists)
* Useful, since we often need a small basis to compute MBPT terms, but a large basis to complete other sum-over-states calculations.
```cpp
Spectrum {
  // exact same inputs as Basis{}
}
```


## Modules

Each Modules block will be run in order.
You can comment-out just the block name, and the block will be skipped.

### Module::matrixElements

Module to calculate Matrix Elements.
For matrixElements, there are some options that apply for any operator; and then there are some options specific to each operator; these operator-specific options are given as a [] bracketed list of options ("sub-block")

```cpp
Module::matrixElements {
  operator;       //[t] default = ""
  options;        //[sub-block], default = ""
  printBoth;      //[t] default = false
  onlyDiagonal;   //[t] default = false
  radialIntegral; //[b] default = false
  rpa;            //[t] default = "TDHF"
  omega;          //[r] default = 0.0;  or [t] ('each')
  A_vertex;       //[r] default = 0.0 - for QED vertex
  b_vertex;       //[r] default = 1.0
}
```
* operator: name of operator; see list below
* options: list any operator-specific options (most will be blank)
  * _see below_
* printBoth: Print <a|h|b> and <b|h|a> ? false by default. (For _some_ operators, e.g., involving derivatives, this is a good test of numerical error. For most operators, values will be trivially the same; reduced matrix elements, sign may be different.)
* onlyDiagonal: If true, will only print diagonal MEs <a|h|a>
* radialIntegral: if true, calculates the radial integral (definition depends on specific operator)
* rpa: Include RPA (core polarisation) corrections to MEs:
  * rpa = TDHF; (default) => uses TDHF method
  * rpa = basis; => uses TDHF/basis method (uses basis)
  * rpa = diagram; => uses diagram (Goldstone) method (uses basis)
  * rpa = false; no RPA included
* omega: frequency for solving TDHF/RPA equations, should be positive. Put "omega=each;" to solve at the frequency for each transition (i.e., re-solve TDHF for each transition).
* A_vertex, b_vertex: Effective QED vertex func: A * alpha * exp(-b * r / alpha)


### Available operators:

Here I list the available operators, and their possible options.
Remember; the program will print out the full list of available options if you ask it to.

```cpp
operator = E1; //Electric dipole operator:
options{
  gauge; //[t] lform, vform. default = lform
}
```

```cpp
operator = Ek; //Electric multipole operator:
options{
  k; //[i] default = 1
}
```
* k=1 => E1, dipole. k=2 => E2, quadrupole etc.

```cpp
operator = r; //scalar r
options{
  power; //[r] default = 1. Will calc <|r^n|>.
}
```

```cpp
operator = pnc; // spin-independent (Qw) PNC operator.
options{
  c; //[r] half-density radius. By default, uses rrms from Z,A [see nucleus]
  t; //[r] skin thickness. default = 2.3
}
```

```cpp
operator = Hrad; // QED radiative potential
```
 * Takes similar arguments as RadPot (except for scale_l and core_qed)
   * Simple, Ueh, SE_h, SE_l, SE_mag, rcut, scale_rN
 * Including RPA should be equivalent to including QED into core HF equations
   * Note: only diagram/basis method seems to work
 * Typically used with radialIntegral=true (get energy shifts)

```cpp
operator = hfs; // Magnetic dipole hyperfine structure constant A
options{
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

### Other Modules:

 * Note: 'Modules' documentation also available in the doxygen (html) documentation on github; that is likely more up-to-date

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
  Module::pnc{ c; t; transition; rpa; omega; nmain; }
```
Uses both 'solving equations' (TDHF) and sum-over-states methods.
For solving equations, calculates both:
  - <yA_w|d| B> + <A |d|xB_w>
  - <A |w|XB_d> + <YA_d|w| B>
  - Does not (yet) include DCP

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

-------------------
```cpp
Module::polarisability{ rpa; omega; transition; omega_max; omega_steps;  }
```
* Calculate dipole polarisabilitities (static, dynamic, alpha, vector beta)
* Uses both 'solving equations' (TDHF) and sum-over-states methods.

 - rpa: true/false. Include RPA or not (TDHF ,method)
 - omega: frequency used for alpha_0 (dipole polarisability). default is 0.
 - transition: For scalar/vector a->b transition polarisability.
   - in form "a,b", e.g., "6s+,7s+" for 6s_1/2 -> 7s_1/2
 - omega_max: maximum frequency for dynamic polarisability. Default is 0.
   - nb: only runs dynamic pol. if omega_max>0
 - omega_steps: Number of steps used for dynamic. default = 30. (linear scale)

Note: transition polarisabilities written for s-states only.
They might be correct for other states too, but NOT checked.
Especially for beta, pretty sure it's wrong for non-s states.

-------------------
```cpp
Module::structureRad{
  operator; options; rpa; printBoth; onlyDiagonal; omega; n_minmax;  
}
```
 * Calculates Structure Radiation + Normalisation of States
 * Note: Most input options are similar to matrixElements module:
 * n_minmax: is input as list of ints:
   * n_minmax = min,max;
   * min: minimum n for core states kept in summations
   * max: maximum n for excited states kept in summations
 * For explanation of the rest, see matrixElements module.

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
Module::HFAnomaly {
  Alist;  //[i,i,...] // Which A's to calculate for (blank for all)
 // ~ most inputs same as operator = hfs
}
```
Calculates the hyperfine anomaly (and BW effect) for all available odd
isotopes of given atom, relative to the 'main' isotope (see Atom).
 * Note: Only runs for odd isotopes; to get anomaly for even isotopes, use am even isotope as reference (For doubly-odd Bohr-Weisskopf input options, see operator = hfs)
 * Takes same input as operator = hfs, except for F(r), since it runs for each F(r)

-------------------
```cpp
Module::BohrWeisskopf { //Calculates BW effect for Ball/Single-particle
  // Takes same input at operator = hfs
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
