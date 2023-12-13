# ampsci: input options

\brief ampsci input: descriptions of file format and available options

[[Home](/README.md)]

## basic usage

The program is run with input options from the command line.

### Main method: input options from a text file

* `./ampsci filename`
  * Runs ampsci with input option specified in file "filename".
  * See below for full description of input format,
and a detailed list of input options + descriptions.
  * Run `./ampsci -h` to get brief instructions for input options
  * Several example input files are given in: _doc/examples/_, along with their expected output;
    use these to test if everything is working.

The Output is printed to screen. It's recommended to forward this to a text file.
The input options and the ampsci version details are also printed, so that the
program output contains all required info to exactly reproduce it. e.g.,

* `./ampsci input |tee -a output`
  * Runs ampsci using input options in file "input".
  * Output will both be written to screen, and appended to
    file "output".

### Command-line options

* run `./ampsci -h` to get brief instructions for input options

```text
    <filename>
        Runs ampsci taking options specified in file "filename" (eg, ./ampsci filename). 
        See documentation (or option -a) for input file format options.
        Example:
        ./ampsci input.in
            -Runs ampsci taking input options from file 'input.in'
    
    <At> <Core> <Valence>
        For quick and simple HF calculation. 
        If core is not given, guesses core configuration and runs using V^N approximation.
        Examples:
        ./ampsci Cs
            - Runs ampsci for Cs using Hartree-Fock (V^N) approximation
        ./ampsci Cs [Xe] 6sd5d
            - Runs ampsci for Cs using Hartree-Fock with Xe-like core and 
              valence states up to n=6 for s,p-states and n=5 for d-states

    -v (--version)
        Prints ampsci version (and git commit) details

    -h (--help, -?)
        Print help info, including some detail on input options

    -a <BlockName> (--ampsci)
        Prints list of available top-level ampsci options. BlockName is optional; 
        if given it will print options for given ampsci Block. You may list any number of blocks (space separated)
        Example:
        ./ampsci -a Atom HartreeFock

    -m <ModuleName> (--modules)
        Prints list of available Modules. ModuleName is optional; if given, will list available options for that Module
        Example:
        ./ampsci -m MatrixElements

    -o <OperatorName> (--operators)
        Prints list of available operators. OperatorName is optional; if given, 
        will list available options for that operator (most operators take no options).
        Example:
        ./ampsci -o E1

    -p <Atom> <Isotope> (--periodicTable)
        Prints textual periodic table with electronic + nuclear information.
        Atom and Isotope are optional; if given, will print info for that isotope. 
        Atom should be atomic symbol (eg Cs), or Z (55).
        If Isotope is blank, will print for 'default' isotope. 
        Can also list 'all' known isotope info
        Examples:
        ./ampsci -p Cs
        ./ampsci -p Cs 131
        ./ampsci -p Cs all
    
    -c (--constants)
        Prints some handy physical constants
```

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

```java
// Available ampsci options/blocks
ampsci{
  // These are the top-level ampsci input blocks and options. Default values are
  // given in square brackets following the description: [default_value]. Blocks
  // end with '{}', options end with ';'. run `ampsci -a BlockName` (for any of
  // the following blocks) to see all the options available for that block.

  // Which atom to run for
  Atom{}
  // Set nuclear parameters
  Nucleus{}
  // Set radial grid/lattice parameters
  Grid{}
  // Options for solving atomic system
  HartreeFock{}
  // Include QED radiative potential
  RadPot{}
  // Include an extra effective potential. Rarely used.
  ExtraPotential{}
  // Basis of HF eigenstates used for MBPT
  Basis{}
  // Options for MBPT and correlation corrections
  Correlations{}
  // Like basis, but includes correlations. Used for sum-over-states
  Spectrum{}
  // Configuration Interaction
  CI{}
  // Run any number of modules (* -> module name). `ampsci -m` to see available
  // modules
  Module::*{}
}
```

You can get the same output by running `./ampsci -a`

Set `help;` inside any of these to get full set of options of each of these,
and so on. Full descriptions of each Block/Option are given in doc/ - but the
self-documentation of the code will always be more up-to-date.

You can get the same output by running `./ampsci -a BlockName`.
For example, `./ampsci -a Basis` will print all available 'Basis' options

The general usage of the code is to first use the main blocks to construct the
atomic wavefunction and basis states, then to add as many 'Module::' blocks as
required. Each module is a separate routine that will take the calculated
wavefunction and compute any desired property (e.g., matrix elements). The code
is designed such that anyone can write a new Module (See [doc/writing_modules.md](/doc/writing_modules.md))

e.g., To calculate Cs wavefunctions at HF level with 6s, 6p, and 5d valence
states, and then calculate E1 matrix elements including core polarisation (RPA):

```java
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

You can also access most of the self-documentation directly from the command-line:

* `./ampsci -h`
  * Print help info, including input options (same as --help, -?)
* `./ampsci -m  <ModuleName>`
  * Prints list of available Modules (same as --modules)
  * ModuleName is optional. If given, will list available options for that Module
* `./ampsci -o <OperatorName>`
  * Prints list of available operators (same as --operators)
  * OperatorName is optional. If given, will list available options for Operator
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
  * ModuleName is optional. If given, will list available options for that Module
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

```java
// Available Atom options/blocks
Atom{
  // Atomic number or symbol (e.g., 55 or Cs). [H]
  Z;
  // Atomic mass number, for nuclear parameters including finite nuclear size.
  // Default based on Z.
  A;
  // Fractional variation of the fine-structure constant, alpha^2: (a/a0)^2. Use
  // to enforce the non-relativistic limit (c->infinity => alpha->0), or
  // calculate sensitivity to variation of alpha. [1.0]
  varAlpha2;
  // Optional label for output identity - for distinguishing outputs with
  // different parameters
  run_label;
}
```

### HartreeFock

```java
// Available HartreeFock options/blocks
HartreeFock{
  // Options for solving lowest-order atomic wavefunction

  // Core configuration. Either list entire core, or use [At] short-hand. e.g.,
  // [He] equivalent to 1s2; [Xe],6s1 equivalent to [Cs] and to
  // 1s2,2s2,...,5p6,6s1. Instead of one of the commas, you may use a ':' -
  // states above this are included into the core, but not excluded from the
  // valence list. Use this method for KohnSham, for example. [blank by default]
  core;
  // Valence configuration in `basis string' format. e.g., 7sp5df will include
  // valence states up to  n=7 for s and p, and up to n=5 for d and f states.
  // Automatically excludes states in the core (except those above the optional
  // ':'). [blank by default]
  valence;
  // HF convergence goal [1.0e-13]
  eps;
  // Method for mean-field approximation: HartreeFock, Hartree, KohnSham, Local
  // [HartreeFock]
  method;
  // Scale for factor for Breit Hamiltonian. Usually 0.0 (no Breit) or 1.0 (full
  // Breit), but can take any value. [0.0]
  Breit;
}
```

### Nucleus

```java
// Available Nucleus options/blocks
Nucleus{
  // Options for nuclear potential (finite nuclear size). All are optional.
  // Default is a Fermi-like nucleus, with parameters chosen according to
  // isotope (see Atom{A;})

  // Root-mean-square charge radius, in fm [default depends on Z and A]
  rrms;
  // Half-density radius, in fm (will over-ride rms) [default depends on Z and
  // A]
  c;
  // Nuclear skin thickness, in fm [2.3]
  t;
  // Fermi, spherical, point-like, Gaussian [Fermi]
  type;
}
```

### Grid

```java
Grid{
  // Options for radial grid (lattice) used for integrations, solving equations
  // and storing orbitals. All relevant quantities are in units of Bohr radius
  // (aB).

  // Initial grid point, in aB [1.0e-6]
  r0;
  // Finial grid point [120.0]
  rmax;
  // Number of grid points [2000]
  num_points;
  // Type of grid: loglinear, logarithmic, linear [loglinear]
  type;
  // Only used for loglinear: grid is ~ logarithmic for r<b, linear for r>b
  // [rmax/3]
  b;
  // du is uniform grid step size; set this instead of num_points - will
  // override num_points [default set by num_points]. Rarely used.
  du;
}
```

### ExtraPotential (read from text file)

```java
// Available ExtraPotential options/blocks
ExtraPotential{
  // Adds an extra potential (to Vnuc), before HF solved. Either effective
  // polarisation potential, V(r) = -0.5/(r^4 + r_cut^4), or read in from a
  // file.

  // Read potential from file (r v(r)) - will be interpolated. If not given,
  // will use pol. potential
  filename;
  // Radial cut-off parameter for effective pol. potential [=1]
  r_cut;
  // Overall scaling factor for potential is scaled by this value [1]
  scale;
}
```

* Reads in extra potential from text file (space separated: 'x y' format):

* Interpolates these points onto the grid (but does NOT extrapolate,
  potential is assumed to be zero outside the given range)
* Potential is multiplied by 'factor'
* May be added before or after HF (if before: added to Vnuc, if after: to Vdir)

### RadPot (Ginges/Flambaum QED Radiative Potential)

```java
RadPot{
  // QED Radiative potential will be included if this block is present

  // The following 5 are all doubles. Scale to include * potential; usually
  // either 0.0 or 1.0, but can take any value:

  //   Uehling (vacuum pol). [1.0]
  Ueh;
  //   self-energy high-freq electric. [1.0]
  SE_h;
  //   self-energy low-freq electric. [1.0]
  SE_l;
  //   self-energy magnetic. [1.0]
  SE_m;
  //   Wickman-Kroll. [0.0]
  WK;
  // Maximum radius (au) to calculate Rad Pot for [5.0]
  rcut;
  // Scale factor for Nuclear size. 0 for point-like, 1 for typical [1.0]
  scale_rN;
  // List of doubles. Extra scaling factor for each l e.g., 1,0,1 => include for
  // s and d, but not for p [1.0]
  scale_l;
}
```

* Adds QED radiative potential to Hamiltonian.

* QED will be included if this block is present; else not
* Will read from file if it exists (e.g., Z_uhlmw.qed)
* Each factor (Ueh, SE_h,..) is a scale; 0 means don't include. 1 means include full potential. Any positive number is valid.
* rcut: Only calculates potential for r < rcut [for speed; rcut in au]
* scale_rN: finite nucleus effects: rN = rN * scale_rN (=0 means point-like)
* scale_l: Optional input: Scaling factors for the V_rad for each l state; for higher states, uses the last given input. Input as a list of real numbers. Best explained with examples:
  * scale_l = 1; // include QED for all states
  * scale_l = 0,1,0; //include QED for p states only
  * scale_l = 0,1; //include QED for p,d,f.. but not s states.
  * don't need to be 1 or 0, can be any real number.
* core_qed: if true, will include QED effects into core in Hartree-Fock (relaxation). If false, will include QED only for valence states

### Basis (B-spline basis for MBPT)

* The 'basis' is used for summing over states in MBPT. (A second 'basis', called spectrum, may be used for summation over states in other problems)

```java
// Available Basis options/blocks
Basis{
  // Number of splines used in expansion [0]
  number;
  // order of splines ~7-9 [7]
  order;
  // minimum cavity radius (first internal knot) [1.0e-4]
  r0;
  // Select cavity radius r0 for each l by position where |psi(r0)/psi_max|
  // falls below r0_eps [1.0e-3]
  r0_eps;
  // maximum cavity radius [Grid{rmax}]
  rmax;
  // states to keep (e.g., 30spdf20ghi)
  states;
  // Force orthogonal to core [false]
  orthogonalise;
  // Print all spline energies (for testing) [false]
  print;
  // Include -ve energy states [false]]
  positron;
  // Derevianko (DKB) or Johnson [Derevianko]
  type;
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

```java
// Available Correlations options/blocks
Correlations{
  // Options for inclusion of correlations (correlation potential method).

  // Minimum core n to polarise [1]
  n_min_core;
  // Construct separate Sigma for each valence state? [false]
  each_valence;
  // List of binding energies (in cm^-1) to scale Sigma for. Must be in same
  // order as valence states
  fitTo_cm;
  // Scaling factors for Sigma. Must be in same order as valence states
  lambda_kappa;
  // Filename to read in Sigma. Set read=false to not read in. By default, will
  // be, e.g., CsI.sig2 or CsI.sigf.
  read;
  // Filename to write Sigma to. Set write=false to not write. By default same
  // as read
  write;
  // minimum radius to calculate sigma for [1.0e-4]
  rmin;
  // maximum radius to calculate sigma for [30.0]
  rmax;
  // Only calculate Sigma every <stride> points. Default such that there are 150
  // points between (1e-4, 30)
  stride;
  // Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127;
  // 7s+=-0.552;}. Blank => HF energies. Takes precedence over each_valence.
  // [blank]
  ek{}
  // Use all-orders method (implies Feynman=true; screening=true;
  // holeParticle=true;) [false]
  allOrder;
  // Use Feynman method [false]
  Feynman;
  // List of doubles. Screening factors for effective all-order exchange. In
  // Feynman method, used in exchange only (and G-part); Goldstone, used direct
  // also. If blank, will calculate them from scratch. []
  fk;
  // List of doubles. Hole-Particle factors. In Feynman method, used only for G
  // part; Goldstone, used in direct also. []
  eta;
  // Include all-orders screening. Only applicable for Feynman method [false]
  screening;
  // Include all-orders hole-particle interaction. Only applicable for Feynman
  // method [false]
  holeParticle;
  // Maximum l used for internal lines in Feynman method [6]
  lmax;
  // Real part of frequency used in contour integral. By Default, ~1/3 of the
  // core/valence energy gap
  real_omega;
  // Pair of comma-separated doubles: w0, wratio. Initial point, and ratio, for
  // logarithmic Im(w) grid [0.01, 1.5]
  imag_omega;
  // Include lower g-part into Sigma [false]
  include_G;
}
```

* Includes correlation corrections. note: splines must exist already

* read/write: Read/write from/to file. Set to 'false' to calculate from scratch (and not write to file). By default, the file name is: "Atom".sig.
  * Alternatively, put any text here to be a custom filename (e.g., read/write="Cs_new"; will read/write from/to Cs_new.sig). Don't include the '.sig' extension (uses sigf for Feynman method, sig2 for Goldstone). Grids must match exactly when reading in from a file.
  * If reading Sigma in from file, basis doesn't need to exist
* n_min_core: minimum core n included in the Sigma calculation; lowest states often contribute little, so this speeds up the calculations
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

```java
// Available Spectrum options/blocks
Spectrum{
  // Options for 'spectrum', Spectrum is the same as 'Basis', but includes
  // correlations. Spectrum is used for sum-over-states (while basis is used for
  // MBPT).

  // Number of splines used in expansion
  number;
  // order of splines ~7-9
  order;
  // minimum cavity radius
  r0;
  // Select cavity radius r0 for each l by position where |psi(r0)/psi_max|
  // falls below r0_eps
  r0_eps;
  // maximum cavity radius
  rmax;
  // states to keep (e.g., 30spdf20ghi)
  states;
  // Force orthogonal to valence [false]
  orthogonalise;
  // Print all spline energies (for testing)
  print;
  // Include -ve energy states (true/false)
  positron;
  // Derevianko (DKB) or Johnson [Derevianko]
  type;
}
```

### CI

```java
// Available CI options/blocks
CI{
  // Basis used for CI expansion; must be a sub-set of full ampsci basis
  // [default: 10spdf]
  ci_basis;
  // List of total angular momentum J for CI solutions (comma separated). Must
  // be integers (two-electron only). []
  J;
  // As above, but for EVEN CSFs only (takes precedence over J).
  J+;
  // As above, but for ODD CSFs (takes precedence over J).
  J-;
  // Number of CI solutions to find (for each J/pi) [5]
  num_solutions;
  // Include one-body MBPT correlations? [false]
  sigma1;
  // Include two-body MBPT correlations? [false]
  sigma2;
  // The subset of ci_basis for which the two-body MBPT corrections are
  // calculated. Must be a subset of ci_basis. If existing sk file has more
  // integrals, they will be used. [default: Nspdf, where N is maximum n for
  // core + 3]
  cis2_basis;
  // Basis used for the one-body MBPT diagrams (Sigma^1). These are the most
  // important, so in general the default (all basis states) should be used.
  // Must be a subset of full ampsci basis. [default: full basis]
  //  - Note: if CorrelationPotential is available, it will be used instead of
  // calculating the Sigma_1 integrals
  s1_basis;
  // Basis used for internal lines of the two-body MBPT diagrams (Sigma^2). Must
  // be a subset of s1_basis. [default: s1_basis]
  s2_basis;
  // Minimum n for core to be included in MBPT [1]
  n_min_core;
  // Maximum k (multipolarity) to include when calculating new Coulomb
  // integrals. Higher k often contribute negligably. Note: if qk file already
  // has higher-k terms, they will be included. Set negative (or very large) to
  // include all k. [8]
  max_k;
  // Filename for storing two-body Coulomb integrals. By default, is At.qk,
  // where At is atomic symbol.
  qk_file;
  // Filename for storing two-body Sigma_2 integrals. By default, is
  // At_n_b_k.sk, where At is atomic symbol, n is n_min_core, b is cis2_basis, k
  // is max_k.
  sk_file;
  // Excludes the Sigma_2 box corrections that have 'wrong' parity when
  // calculating Sigma2 matrix elements. Note: If existing sk file already has
  // these, they will be included [false]
  exclude_wrong_parity_box;
  // Sort output by energy? Default is to sort by J and Pi first. [false]
  sort_output;
  // Run CI in parallel (solve each J/Pi in parallel). Faster, uses slightly
  // more memory [true]
  parallel_ci;
}
```
