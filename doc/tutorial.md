# Basic Tutorial

\brief Basic tutorial for using ampsci, including examples.

[[Home](/README.md)]

This assumes you already have ampsci compiled.

* See [doc/compilation.md](/doc/compilation.md) for compilation instructions
* See [doc/tutorial_advanced.md](/doc/tutorial_advanced.md) for a more advanced tutorial, including correlation corrections
* and [doc/tutorial_CI.md](/doc/tutorial_CI.md) for tutorial on Configuration Interaction (+MBPT) calculations for two-valence systems

## Contents

* [Getting started](#start)
* [Running basic calculations](#basic)
* [Setting up the input file](#input)
  * [Example: running calculation using input file](#example)
* [Output format](#output)
* [Modules: using the wavefunctions](#modules)

## Getting started <a name="start"></a>

Open up a terminal and navigate to the directory where ampsci was compiled.

* If you are unsure using the terminal, see [this tutorial (ubuntu)](https://ubuntu.com/tutorials/command-line-for-beginners) (It's specifically for ubuntu (linux), but is mostly the same on mac)

First, just run the executable:

```bash
./ampsci
```

No calculation was performed, however, some instructions should be printed to the screen. These should serve as a reminder if you forget any of the commands we will look at below. You can also run:

```bash
./ampsci -h
```

which will print the same, along with some more detailed 'help' information.
You can ignore this output for now, as we'll work through the basic examples.

Now, try the following, which should print the current _version_ info (including the git commit hash, if you're using git). You can use this to check which exact version of the code you are running. This is also automatically printed when you run any calculation.

```bash
./ampsci -v
```

## Running basic calculations <a name="basic"></a>

To actually perform an atomic calculation, we need to tell ampsci what we want to calculate.

Typically, this is done using an input file where we specify all the options.
For very simple calculations, however, we can do this directly from the command line using the form:

```bash
./ampsci <Atom> <Core> <Valence>
```

(where `<Core>` and `<Valence>` are both optional).
For example, the following should all produce the same result. They will calculate Hartree-Fock for Cs, including all the electrons in the "core" (so-called $V^N$ approximation).

```bash
./ampsci Cs
./ampsci Cs [Cs]
./ampsci Cs [Xe],6s1
```

* The format for the "core" configuration copies that used by NIST in their periodic table. A nice copy can be downloaded from [here](https://www.nist.gov/image/periodic-table).
* Usually the term in the brackets `[]` is a noble gas, but in ampsci you can put any atom (the core configuration for non-noble-gas atoms is guessed based on simple filling rules, and is sometimes not correct, so always check the output!).

The next few cases will calculate Hartree-Fock for Cs, using the $V^{N-1}$ approximation, including valence states up to $n=7$ for $s$ and $p$ states, and $n=5$ for $d$ states:

```bash
./ampsci Cs [Xe] 7sp5d
```

## Setting up the input file <a name="input"></a>

For anything more complicated, we must use an input file.
The input file is a plain-text file; a full description of the format is given [elsewhere](/doc/ampsci_input.md) - here we will run through a basic example.

* See [doc/ampsci_input.md](/doc/ampsci_input.md) for detail input descriptions

Firstly, we can use the code to tell us which input options are available using the `-a` (or `--ampsci`) command-line option:

```bash
./ampsci -a
```

This will print the list of available top-level "Input Blocks".
The output will look something like this:

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

This will be in the correct format for the input file, so we can start making our input file by copy+pasting these.
We don't need all the options, so let's start with just the basics: `Atom{}`, `Grid{}`, and `HartreeFock{}`.
In each of these blocks, we can specify certain options; most have default values, and can be thus left blank.
We can also ask the code to tell us which options are available for each block, for example:

```bash
./ampsci -a Atom
```

The output will list all options for the `Atom{}` block:

```text
Available Atom options/blocks are:
Atom{
  Z;    // string or int (e.g., Cs equivalent to 55). Atomic number [default H]
  A;    // int. Atomic mass number (set A=0 to use pointlike nucleus) [default based on Z]
  varAlpha2;    // Fractional variation of the fine-structure constant, alpha^2: d(a^2)/a_0^2. Use to enforce the non-relativistic limit (c->infinity => alpha->0), or calculate sensitivity to variation of alpha. [1.0]
}
```

This is also in the correct format, so we can copy+paste into our input file.

* Anything after `//` is a comment and will be ignored by the program.
* Don't forget to close the curly braces `}` at the end of each input block, and the semi-colon `;` after each option!

We can list more than one block name for this. For example, try:

```bash
./ampsci -a Grid HartreeFock
```

which will give something like:

```java
// Available Grid options/blocks
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


// Available HartreeFock options/blocks
HartreeFock{
  // Options for solving lowest-order atomic wavefunction

  // Core configuration. Either list entire core, or use [At] short-hand. e.g.,
  // [He] equivilant to 1s2; [Xe],6s1 equivilant to [Cs] and to
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

### Example: running calculation using input file <a name="example"></a>

We will use the above to create our first simple
We will call this text file `example.in`, but any filename is OK (the `.in` file extension is just by convention; you may use `.txt` or anything else).
We may then set up the input file like:

```java
// example.in
Atom{
  Z = Cs; // Or Z = 55;
  A = 133; // may leave blank for default isotope
}
Grid{
  r0 = 1.0e-6;
  rmax = 150.0;
  num_points = 3000;
}
HartreeFock{
  core = [Xe];
  valence = 7sp5d;
}
```

From the `Atom{}` block, we have chosen to run calculations for Cs-133. Since we did not set any parameters in the `Nucleus{}` block, the code will run using the default nuclear parameters (which are looked up based on the isotope).

In the `Grid{}` block, we specify the discrete radial grid to have `3000` grid points dispersed between the range `[1.0e-6, 150.0]` atomic units (atomic unit of length is the Bohr radius: $1\,{\rm au}=1\,a_0 = \tfrac{4\pi\epsilon_0\hbar^2}{m_e e^2}\approx 5.3\times10^{-11}\,{\rm m}$).
Since we did not specify a grid `type` the default (log-linear) will be used (see [ampsci.pdf](https://ampsci.dev/ampsci.pdf)).

Finally, in the `HartreeFock{}` block, we specify that we wish to solve Hartree-Fock equations for Cs with a Xe-like core, and consider valence states up to $n=5$ for $s,p$-states, and $n=5$ for $d$-states.
Since a Xenon-like core has 1 fewer electrons than neutral Cs, this is a $V^{N-1}$ calculation.
Note that the format for the core configuration copies that used by NIST in their periodic table. A nice copy can be downloaded from [here](https://www.nist.gov/image/periodic-table).

You may also use the built-in periodic table to see the possible ground-state electron configuration:

```bash
./ampsci -p
or
./ampsci -p Cs
```

Caution: this "guesses" the configuration based on simple filling rules, so may not be the true ground-state configuration.

## Output format <a name="output"></a>

If we run the code with the above input file:

```bash
./ampsci example.in
```

we will get something like the following output:

```text
********************************************************************************
AMPSCI v: 0.0 [dev/f06534f0]
Compiled: g++-11 [Ubuntu 11.1.0-1ubuntu1~18.04.1] 11.1.0 2022-09-14 12:42 AEST
Run time: 2022-09-14 17:01:02

********************************************************************************
Atom { 
  Z = Cs;
  A = 133;
}
Grid { 
  r0 = 1.0e-6;
  rmax = 150.0;
  num_points = 3000;
}
HartreeFock { 
  core = [Xe];
  valence = 7sp5d;
}

Running for Cs, Z=55 A=133
Fermi nucleus;  r_rms = 4.8041, c_hdr = 5.67073, t = 2.3
Log-linear (b=50) grid: 1e-06->150, N=3000, du=0.36389
========================================================
Hartree-Fock
Core   :  it: 28 eps=8.2e-14 for 5p+
Val    :  it: 38 eps=0.0e+00 for 6s+ [ 38 eps=0e+00 for 6s+]

CsI-133
Core: [Xe] (V^N-1)
     state  k   Rinf its   eps         En (au)        En (/cm)
0   1s_1/2 -1    0.7  2  3e-26 -1330.118948369  -291927365.862
1   2s_1/2 -1    1.7  2  2e-23  -212.564531430   -46652522.176
2   2p_1/2  1    1.7  2  9e-24  -199.429544824   -43769725.833
3   2p_3/2 -2    1.8  2  1e-23  -186.436652692   -40918115.622
4   3s_1/2 -1    3.6  2  1e-21   -45.969754343   -10089194.888
5   3p_1/2  1    3.8  2  7e-22   -40.448315684    -8877379.174
6   3p_3/2 -2    4.0  2  7e-22   -37.894321341    -8316842.207
7   3d_3/2  2    4.6  2  4e-22   -28.309520222    -6213221.515
8   3d_5/2 -3    4.6  2  5e-22   -27.775176728    -6095946.673
9   4s_1/2 -1    7.9  2  1e-20    -9.512819127    -2087822.471
10  4p_1/2  1    9.0  2  8e-21    -7.446283820    -1634270.397
11  4p_3/2 -2    9.3  2  9e-21    -6.920999547    -1518983.824
12  4d_3/2  2   13.2  2  6e-21    -3.485619409     -765005.035
13  4d_5/2 -3   13.3  2  6e-21    -3.396901915     -745533.796
14  5s_1/2 -1   20.3  2  7e-21    -1.489803443     -326974.061
15  5p_1/2  1   26.3  2  3e-21    -0.907896946     -199260.348
16  5p_3/2 -2   27.3  2  3e-21    -0.840338411     -184432.963
E_c = -7786.643737
Val: state  k   Rinf its   eps         En (au)        En (/cm)   En (/cm)
0   6s_1/2 -1   70.8  1  0e+00    -0.127368053      -27954.056       0.00
1   7s_1/2 -1  110.4  1  0e+00    -0.055187351      -12112.224   15841.83
2   6p_1/2  1   87.2  1  0e+00    -0.085615846      -18790.506    9163.55
3   7p_1/2  1  128.1  1  0e+00    -0.042021373       -9222.625   18731.43
4   6p_3/2 -2   88.2  1  0e+00    -0.083785436      -18388.778    9565.28
5   7p_3/2 -2  129.1  1  0e+00    -0.041368028       -9079.233   18874.82
6   5d_3/2  2  101.5  1  0e+00    -0.064419644      -14138.478   13815.58
7   5d_5/2 -3  101.5  1  0e+00    -0.064529777      -14162.649   13791.41

ampsci: T = 522.26 ms
```

For calculations that matter, the _entire_ output should be saved.
The best way to do this is to forward the output directly to a text file when running the code.
The input options and the ampsci version details are also printed, so that the
program output contains all required info to exactly reproduce it.
A nice way to do this is by `piping` the output through the unix `tee` program:
e.g.,

```bash
./ampsci example.in |tee -a example.out
```

* Runs ampsci using input options in file "example.in".
* Output will both be written to screen, and _appended_ (due to the `-a`) to file "example.out".

It's important to understand the meaning of the program output.
The first part of the output just prints the version information for ampsci.
This is important for reproducing old results, which may depend on the ampsci version.

Then, the entire set of input options is printed.
This is useful, since it tells you explicitly which calculations were run, and is very useful for catching mistakes and reproducing old calculations.
Importantly, the format is exactly what is required on input, so to re-run the calculation, this can simply be copy+pasted into a new input file.

```text
Running for Cs, Z=55 A=133
Fermi nucleus;  r_rms = 4.8041, c_hdr = 5.67073, t = 2.3
Log-linear (b=50) grid: 1e-06->150, N=3000, du=0.36389
```

* This tells you atom and isotope for which the calculations were run, the exact nuclear parameters used (these are for the finite-nuclear size effect), and the exact grid parameters used.

```text
Hartree-Fock
Core   :  it: 28 eps=8.2e-14 for 5p+
Val    :  it: 38 eps=0.0e+00 for 6s+ [ 38 eps=0e+00 for 6s+]
```

* This is one of the most important outputs to check.
* This tells us that Hartree-Fock equations for the core converged to parts in $10^{14}$ in 28 iterations, and the worst core state convergence was for $5p_{3/2}$.
* The worst valence state converged to `0` (i.e., floating point underflowed) in 38 iterations
* It's important to check that none of these 'epsilon' values are large (i.e., $\epsilon\lesssim10^{-6}$), otherwise it means Hartree-Fock didn't converge properly, and the calculations will be unreliable. This is rarely an issue.

```text
CsI-133
Core: [Xe] (V^N-1)
     state  k   Rinf its   eps         En (au)        En (/cm)
0   1s_1/2 -1    0.7  2  3e-26 -1330.118948369  -291927365.862
1   2s_1/2 -1    1.7  2  2e-23  -212.564531430   -46652522.176
2   2p_1/2  1    1.7  2  9e-24  -199.429544824   -43769725.833
3   2p_3/2 -2    1.8  2  1e-23  -186.436652692   -40918115.622
4   3s_1/2 -1    3.6  2  1e-21   -45.969754343   -10089194.888
5   3p_1/2  1    3.8  2  7e-22   -40.448315684    -8877379.174
6   3p_3/2 -2    4.0  2  7e-22   -37.894321341    -8316842.207
7   3d_3/2  2    4.6  2  4e-22   -28.309520222    -6213221.515
8   3d_5/2 -3    4.6  2  5e-22   -27.775176728    -6095946.673
9   4s_1/2 -1    7.9  2  1e-20    -9.512819127    -2087822.471
10  4p_1/2  1    9.0  2  8e-21    -7.446283820    -1634270.397
11  4p_3/2 -2    9.3  2  9e-21    -6.920999547    -1518983.824
12  4d_3/2  2   13.2  2  6e-21    -3.485619409     -765005.035
13  4d_5/2 -3   13.3  2  6e-21    -3.396901915     -745533.796
14  5s_1/2 -1   20.3  2  7e-21    -1.489803443     -326974.061
15  5p_1/2  1   26.3  2  3e-21    -0.907896946     -199260.348
16  5p_3/2 -2   27.3  2  3e-21    -0.840338411     -184432.963
E_c = -7786.643737
```

* This is the output for the core states.
* The columns are:
  * `state`: which core state
  * `k`: value of $\kappa$ (Dirac quantum number) for state
  * `Rinf`: Radius of the 'practical infinity', i.e., point on the radial grid where that orbital can be safely said to be zero. It's crucial to check all of these values lie well below `rmax` (set in `Grid{}`). If these are close to `rmax`, `rmax` should be increased.
  * `its`: Number of iterations it took to solve the single-particle Dirac equation upon the final Hartree-Fock iteration
  * `eps`: Extent to which the single-particle energy converged when solving the Dirac equation upon the last Hartree-Fock iteration. These numbers should be extremely small; if any are larger, it means there is possibly a serious numerical error. Normally, these will always be fine, so long as Hartree-Fock converged
  * `En (au)` and `En (/cm)` are the single-particle (binding) energies in atomic units and ${\rm cm}^{-1}$, respectively. Atomic unit for energy is the Hartree:
  $E_h = 2 Ry = \frac{\hbar^2}{m_e a_0^2} \approx 27.2\,{\rm eV}$.
  * The final number `E_c` is the total Hartree-Fock energy of the core, in atomic units

```text
Val: state  k   Rinf its   eps         En (au)        En (/cm)   En (/cm)
0   6s_1/2 -1   70.8  1  0e+00    -0.127368053      -27954.056       0.00
1   7s_1/2 -1  110.4  1  0e+00    -0.055187351      -12112.224   15841.83
2   6p_1/2  1   87.2  1  0e+00    -0.085615846      -18790.506    9163.55
3   7p_1/2  1  128.1  1  0e+00    -0.042021373       -9222.625   18731.43
4   6p_3/2 -2   88.2  1  0e+00    -0.083785436      -18388.778    9565.28
5   7p_3/2 -2  129.1  1  0e+00    -0.041368028       -9079.233   18874.82
6   5d_3/2  2  101.5  1  0e+00    -0.064419644      -14138.478   13815.58
7   5d_5/2 -3  101.5  1  0e+00    -0.064529777      -14162.649   13791.41
```

* Finally, we have the output for the valence states
* Most of these are the same as above, except there is an extra column, which gives the excitation energies (in ${\rm cm}^{-1}$) relative to the lowest valence state (useful for comparing with [NIST database](https://physics.nist.gov/PhysRefData/ASD/levels_form.html))

To do more complicated calculations, including constructing complete set of basis orbitals, and including core-valence correlations, see:

* [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full physics description of what the code can do
* [doc/ampsci_input.md](/doc/ampsci_input.md) for detail on all input options

## Modules: using the wavefunctions <a name="modules"></a>

Above, we ran ampsci, which calculated the atomic wavefunctions and printed their energies to screen.
If we want to actually _do_ anything with the wavefunctions, we have to run one or more **modules**.
We do this by adding a module block to the input file, which has the form `Module::ModuleName{}`

Here, we will just consider a simple example. For further detail:

* See [doc/modules.md](/doc/modules.md) for details of currently avalable modules
* The code is designed so that you can easily create your own modules. See [doc/writing_modules.md](/doc/writing_modules.md) for details

We can see a list of all available modules with the `-m` command-line option

```bash
./ampsci -m
```

Many are listed. Here we will consider the simple/common example of `matrixElements` module, which (unsurprisingly) calculates matrix elements.

To see the available options for this block, list the block name after `-m` on the command-line:

```bash
./ampsci -m matrixElements
```

which prints:

```java
// Available Module::matrixElements options/blocks
Module::matrixElements{
  // e.g., E1, hfs (see ampsci -o for available operators)
  operator;
  // options specific to operator (see ampsci -o 'operator')
  options{}
  // Method used for RPA: true(=TDHF), false, TDHF, basis, diagram [true]
  rpa;
  // Text or number. Freq. for RPA (and freq. dependent operators). Put 'each'
  // to solve at correct frequency for each transition. [0.0]
  omega;
  // What to calculate: rme (reduced ME), A (hyperfine A/B coeficient). Default
  // is rme, except when operator=hfs, in which case default is A
  what;
  // print <a|h|b> and <b|h|a> [false]
  printBoth;
  // If true (and spectrum available), will use spectrum for valence states AND
  // for RPA (if diagram method), AND for SR+N [false]
  use_spectrum;
  // Calculate diagonal matrix elements (if non-zero) [true]
  diagonal;
  // Calculate off-diagonal matrix elements (if non-zero) [true]
  off-diagonal;
  // Options for Structure Radiation and normalisation (details below)
  StructureRadiation{}
}


// Available StructureRadiation options/blocks
StructureRadiation{
  // If this block is included, SR + Normalisation corrections will be included

  // filename for QkTable file. If blank will not use QkTable; if exists, will
  // read it in; if doesn't exist, will create it and write to disk. Save time
  // (10x) at cost of memory. Note: Using QkTable implies splines used for
  // diagram legs
  Qk_file;
  // list; min,max n for core/excited: [1,inf]
  n_minmax;
}
```

The first option, `operator` is which operator we want to calculate matrix elements of.
You can see a list of operators with the `-o` command-line option:

```bash
./ampsci -o
```

The second option, which is a sub-input-block, `options` is the set of options for the given operator. Most operators have no extra options, so this can be left blank. Some (e.g., hyperfine operator `hfs` have many available options). Like above, you can see the available options by further passing the operator name after `-o`. For example:

```bash
./ampsci -o hfs
```

Here we will consider the simpler `E1` operator.
To our above `example.in` file, we can add the following block (note we may add as many Module:: blocks as we like, they will all be run one-by-one in order):

```java
// example.in
// ... above input options ...
Module::Matrixelements{
  operator = E1;
  rpa = true;
  omega = 0.0;
}
```

The `omega = 0.0;` option means we have solved RPA equations for each transition at zero frequency. It is also possible to automatically solve RPA for each transition frequency by setting: `omega = each;`.
See [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for description of RPA.

The output format will depend on the specific module.
In this case, we are shown the extent to which the RPA (/TDHF) equations converged, and then the valued of the reduced matrix elements are printed (at the Hartree-Fock, first-order core-polarisation, and all-orders RPA levels)
