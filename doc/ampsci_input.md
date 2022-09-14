# ampsci: input options

\brief ampsci input: descriptions of file format and available options

[[Home](/README.md)]

## basic usage

The program is run with input options from the command line.

### Main method: input options from a text file

* `./ampsci filename`
  * Runs ampsci with input option specified in file "filename"
  * See below for full description of input format,
and a detailed list of input options + descriptions.
  * run `./ampsci -h` to get breif instructions for input options
  * Several example input files are given in: _doc/examples/_, along with their expected output; use these to test if everything is working.

The Output is printed to screen. It's recommended to forward this to a text file.
The input options and the ampsci version details are also printed, so that the
program output contains all required info to exactly reproduce it. e.g.,

* `./ampsci input |tee -a outout`
  * Runs ampsci using input options in file "input".
  * Output will both be written to screen, and appended to
    file "output".

### quick method (simple calculations)

For very simple (Hartree-Fock only) calculations, you can run ampsci directly from the command line:

* `./ampsci <At> <Core> <Valence>`
  * `./ampsci Cs`
    * Runs ampsci for Cs using Hartree Fock (V^N) approximation
  * `./ampsci Cs [Xe] 6sd5d`
    * Runs ampsci for Cs using Hartree Fock with Xe-like core and valence
      states up to n=6 for s,p-states and n=5 for d-states
  * `./ampsci Cs`
    * Runs ampsci for Cs using Hartree Fock (V^N) approximation

### Other command-line options

* `./ampsci -v`
  * Prints version info (same as --version)
* `./ampsci -h`
  * Print help info, including input options (same as --help, -?)
* `./ampsci -m  <ModuleName>`
  * Prints list of available Modules (same as --modules)
  * ModuleName is optional. If given, will list avaiable options for that Module
* `./ampsci -o <OperatorName>`
  * Prints list of available operators (same as --operators)
  * OperatorName is optional. If given, will list avaiable options for Operator
* `./ampsci -a <BlockName>`
  * Prints list of available top-level ampsci options (same as --ampsci)
  * BlockName is optional; if given will print options for given ampsci Block
  * e.g., `./ampsci -a Basis` will print all available 'Basis' options
* `./ampsci -p <At> <Isotope>`
  * Prints periodic table with electronic+nuclear info (same as --periodicTable)
  * At and Isotope are optional. If given, will print info for given isotope
  * e.g., `./ampsci -p Cs`, `./ampsci -p Cs 133`, `./ampsci -p Cs all`
* `./ampsci -c`
  * Prints some handy physical constants (same as --constants)

--------------------------------------------------------------------------------

## Input file format

Input is a plain text file that consists of sets of 'Blocks' and 'Options'.

* Blocks are followed by curly-braces: BlockName{}
* Options are followed by a semi-colon: OptionName = option_value;
* Generally, each Block will have a set of Options that may be set
* Nearly all are optional - leave them blank and a default value will be used
* Blocks may be nested inside other Blocks
* White-space is ignored, as are ' and " characters
* You may use C++-style line '//' and block '/**/' comments

The code is "self-documenting". At any level (i.e., in any Block or at 'global'
level outside of any Block), set the option 'help;', and the code will print:

* a list of all available Blocks and Options at that level
* a description of what they are for, and
* the default value if they are left unset.

For example, setting 'help' at the top-level will print a list of all available
top-level Blocks:

```cpp
  Atom{}         // InputBlock. Which atom to run for
  Grid{}         // InputBlock. Set radial grid parameters
  HartreeFock{}  // InputBlock. Options for Solving atomic system
  Nucleus{}      // InputBlock. Set nuclear parameters
  RadPot{}       // InputBlock. Inlcude QED radiative potential
  Basis{}        // InputBlock. Basis used for MBPT
  Spectrum{}     // InputBlock. Like basis; used for sum-over-states
  Correlations{} // InputBlock. Options for correlations
  ExtraPotential{} // InputBlock. Include an extra potential
  dVpol{}        // InputBlock. Approximate correlation (polarisation) potential
  Module::*{}    // InputBlock. Run any number of modules (* -> module name)
```

You can get the same output by running `./ampsci -a`

Set 'help;' inside any of these to get full set of options of each of these,
and so on. Full descriptions of each Block/Option are given in doc/ - but the
self-documentation of the code will always be more up-to-date.

You can get the same output by running `./ampsci -a BlockName`.
For example, `./ampsci -a Basis` will print all available 'Basis' options

The general usage of the code is to first use the main blocks to construct the
atomic wavefunction and basis states, then to add as many 'Module::' blocks as
required. Each module is a seperate routine that will take the calculated
wavefunction and compute any desired property (e.g., matrix elements). The code
is designed such that anyone can write a new Module (See [doc/writing_modules.md](/doc/writing_modules.md))

e.g., To calculate Cs wavefunctions at HF level with 6s, 6p, and 5d valence
states, and then calculate E1 matrix elements including core polarisation (RPA):

```cpp
  Atom {
    Z = Cs;
    A = 133;
  }
  Grid { } // Leave all default; can also just drop entire Block
  Nucleus { } // Default values set according to isotope
  HartreeFock {
    core = [Xe];
    valence = 6sp5d;
  }
  Module::matrixElements {
    operator = E1;
    rpa = true;
  }
```

--------------------------------------------------------------------------------

## Auto-documentation

* This document may go out-of-sync with the code
* The best fail-safe way to check the available options for a given input block is to use the code itself.
* Do this by adding a blank 'help' option to the input file:
* The code will then print a list of all available options, and (usually) an explanation for them

blockname {
  help;
}

You can also access most of the self-documenation directly from the command-line:

* `./ampsci -h`
  * Print help info, including input options (same as --help, -?)
* `./ampsci -m  <ModuleName>`
  * Prints list of available Modules (same as --modules)
  * ModuleName is optional. If given, will list avaiable options for that Module
* `./ampsci -o <OperatorName>`
  * Prints list of available operators (same as --operators)
  * OperatorName is optional. If given, will list avaiable options for Operator
* `./ampsci -a <BlockName>`
  * Prints list of available top-level ampsci options (same as --ampsci)
  * BlockName is optional; if given will print options for given ampsci Block
  * e.g., `./ampsci -a Basis` will print all available 'Basis' options

--------------------------------------------------------------------------------

## Modules

* The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
* Any number of _modules_ can be run by adding a `Module::moduleName{}' block.
* Get a list of available modules: `./ampsci -m`
* `./ampsci -m  <ModuleName>`
  * Prints list of available Modules (same as --modules)
  * ModuleName is optional. If given, will list avaiable options for that Module
* See [doc/modules.md](/doc/modules.md) for full details
* The code is designed so that you can easily create your own modules. See [doc/writing_modules.md](/doc/writing_modules.md) for details

--------------------------------------------------------------------------------

## Details for each input block

* Basic details for each input block are given here.
* It's generally better to get this info from the code (by setting the `help;` option in any given block), since that will always be up-to-date, while this document may fall out-of-date
* You can also get this directly from the command-line:
* `./ampsci -a <BlockName>`
  * Prints list of available top-level ampsci options (same as --ampsci)
  * BlockName is optional; if given will print options for given ampsci Block
  * e.g., `./ampsci -a Basis` will print all available 'Basis' options

### Atom

Available Atom options/blocks are:

```cpp
Atom{
  Z; // string or int (e.g., Cs equivilant to 55). Atomic number [default H]
  A; // int. Atomic mass number (set A=0 to use pointlike nucleus) [default based on Z]
  varAlpha2; // Fractional variation of the fine-structure constant, alpha^2: d(a^2)/a_0^2. Use to enforce the non-relativistic limit (c->infinity => alpha->0), or calculate sensitivity to variation of alpha. [1.0]
}
```

### HartreeFock

```cpp
HartreeFock{
  /*
  Options for solving lowest-order atomic wavefunction
  */
  core; // Core configuration. Either list entire core, or use [At] short-hand. e.g., [He] equivilant to 1s2; [Xe],6s1 equivilant to [Cs] and to 1s2,2s2,...,5p6,6s1. [blank by default]
  valence; // e.g., 7sp5d will include valence states up to n=7 for s and p, but n=5 for d states. Automatically excludes states in the core. [blank by default]
  eps; // HF convergance goal [1.0e-13]
  method; // HartreeFock, Hartree, KohnSham, Local [HartreeFock]
  Breit; // Scale for factor for Breit Hamiltonian. Usially 0.0 (no Breit) or 1.0 (full Breit), but can take any value. [0.0]
  sortOutput; // Sort energy tables by energy? [false]
}
```

### Nucleus

```cpp
Nucleus{
  /*
  Options for nuclear potential (finite nuclear size). All are optional. Default is a Fermi-like nucleus, with parameters chosen according to isotope (see Atom{A;})
  */
  rrms; // Root-mean-square charge radius, in fm [default depends on Z and A]
  c; // Half-density radius, in fm (will over-ride rms) [default depends on Z and A]
  t; // Nuclear skin thickness, in fm [2.3]
  type; // Fermi, spherical, pointlike, Gaussian [Fermi]
}
```

### Grid

```cpp
Grid{
  /*
  Options for radial grid (lattice) used for integrations, solving equations and storing oritals. All relevant quantities are in units of Bohr radius (aB).
  */
  r0; // Initial grid point, in aB [1.0e-6]
  rmax; // Finial grid point [120.0]
  num_points; // Number of grid points [2000]
  type; // Type of grid: loglinear, logarithmic, linear [loglinear]
  b; // Only used for loglinear: grid is ~ logarithmic for r<b, linear for r>b [rmax/3]
  du; // du is uniform grid step size; set this instead of num_points - will override num_points [default set by num_points]. Rarely used.
}
```

### dVpol (effective polarisation potential)

```cpp
Available dVpol options/blocks are:
dVpol{
  a_eff; // scale factor for effective pol. potential [1]
  r_cut; // cut-off parameter [=1]
}
```

* Effective polarisation potential:

* dV = -0.5 * a_eff / (r^4 + r_cut^4)
* nb: Added to direct potential _after_ HF for core, but _before_ HF for valence

### ExtraPotential (read from text file)

```cpp
ExtraPotential{
  /*
  Option to add an extra potential (to Vnuc), before HF solved.
  */
  filename; // Read potential from file (r v(r)) - will be interpolated [blank]
  factor; // potential is scaled by this value [default=1]
  beforeHF; // include before HF (into core states). default=false
}
```

* Reads in extra potential from text file (space separated: 'x y' format):

* Interpolates these points onto the grid (but does NOT extrapolate,
  potential is assumed to be zero outside the given range)
* Potential is multiplied by 'factor'
* May be added before or after HF (if before: added to vnuc, if after: to vdir)

### RadPot (Ginges/Flambaum QED Radiative Potential)

```cpp
RadPot{
  /*
  QED Radiative potential will be included if this block is present
  */
  /*
  The following 5 are all doubles. Scale to include * potential; usually either 0.0 or 1.0, but can take any value:
  */
  Ueh; //   Uehling (vacuum pol). [1.0]
  SE_h; //   self-energy high-freq electric. [1.0]
  SE_l; //   self-energy low-freq electric. [1.0]
  SE_m; //   self-energy magnetic. [1.0]
  WK; //   Wickman-Kroll. [0.0]
  rcut; // Maximum radius (au) to calculate Rad Pot for [5.0]
  scale_rN; // Scale factor for Nuclear size. 0 for pointlike, 1 for typical [1.0]
  scale_l; // List of doubles. Extra scaling factor for each l e.g., 1,0,1 => include for s and d, but not for p [1.0]
  core_qed; // Include rad pot into Hartree-Fock core (relaxation) [true]
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

### Basis (B-spline basis for MBPT)

* The 'basis' is used for summing over states in MBPT. (A second 'basis', called spectrum, may be used for summation over states in other problems)

```cpp
Basis{
  number; // Number of splines used in expansion [0]
  order; // order of splines ~7-9 [7]
  r0; // minimum cavity radius (first internal knot) [1.0e-4]
  r0_eps; // Select cavity radius r0 for each l by position where |psi(r0)/psi_max| falls below r0_eps [1.0e-3]
  rmax; // maximum cavity radius [Grid{rmax}]
  states; // states to keep (e.g., 30spdf20ghi)
  print; // Print all spline energies (for testing) [false]
  positron; // Include -ve energy states [false]]
  type; // Derevianko (DKB) or Johnson [Derevianko]
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

### Correlations (Correlation potential, Sigma)

* For including correlations. 'basis' must exist to calculate Sigma, but not to read Sigma in from file.

```cpp
Correlations{
  /*
  Options for inclusion of correlations (correlation potential method). It's become a bit of a mess, and will be refactored ~soon~
  */
  Brueckner; // Form Brueckner orbitals [false]
  energyShifts; // Calculate MBPT2 shift [false]
  n_min_core; // Minimum core n to polarise [1]
  fitTo_cm; // List of binding energies (in cm^-1) to scale Sigma for. Must be in same order as valence states
  lambda_kappa; // Scaling factors for Sigma. Must be in same order as valence states
  read; // Filename to read in Sigma [false=don't read]
  write; // Filename to write Sigma to [false=don't write]
  rmin; // minimum radius to calculate sigma for [1.0e-4]
  rmax; // maximum radius to calculate sigma for [30.0]
  stride; // Only calculate Sigma every <stride> points
  each_valence; // Different Sigma for each valence states? [false]
  ek; // Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127, 7s+=-0.552;}. Blank => HF energies
  Feynman; // Use Feynman method [false]
  fk; // List of doubles. Screening factors for effective all-order exchange. In Feynman method, used in exchange+ladder only; Goldstone, used direct also. If blank, will calculate them from scratch. []
  eta; // List of doubles. Hole-Particle factors. In Feynman method, used in ladder only; Goldstone, used direct also. []
  screening; // bool. Include Screening [false]
  holeParticle; // Include hole-particle interaction [false]
  ladder; // Experimental feature. Filename for ladder diagram file (generated in the ladder Module). If blank, ladder not included. Only in Feynman. []
  lmax; // Maximum l used for Feynman method [6]
  basis_for_Green; // Use basis for Feynman Greens function [false]
  basis_for_pol; // Use basis for Feynman polarisation op [false]
  real_omega; // [worked out by default]
  imag_omega; // w0, wratio for Im(w) grid [0.01, 1.5]
  include_G; // Inlcude lower g-part into Sigma [false]
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
  * If blank, will calculate these from scratch for each state (better, slower)
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

### Spectrum (B-spline basis for MBPT)

* The 'Spectrum' is similar to basis, but also includes correlation corrections (if Sigma exists)

* Useful, since we often need a small basis to compute MBPT terms, but a large basis to complete other sum-over-states calculations.

```cpp
Spectrum{
  /*
  Options for 'spectrum', Spectrum is the same as 'Basis', but includes correlations. Spectrum is used for sum-over-states (while basis is used for MBPT).
  */
  number; // Number of splines used in expansion
  order; // order of splines ~7-9
  r0; // minimum cavity radius
  r0_eps; // Select cavity radius r0 for each l by position where |psi(r0)/psi_max| falls below r0_eps
  rmax; // maximum cavity radius
  states; // states to keep (e.g., 30spdf20ghi)
  print; // Print all spline energies (for testing)
  positron; // Include -ve energy states (true/false)
  type; // Derevianko (DKB) or Johnson [Derevianko]
}
```
