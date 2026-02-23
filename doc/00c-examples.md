\page examples Basic Usage and Example
\ingroup getting_started
\brief A few very basic examples

Here is a very breif introduction to running ampsci with an example input file.
See the [detailed tutorials](\ref tutorials) for more.

## Input options, self-deocumentation

Typically, ampsci will be run taking input from a plain text input file (in this case, called `example.in`):

```bash
./ampsci example.in
```

There are many input options available.
We can use the code to tell us which input options are available using the `-i` (or `--ampsci`) command-line option. For example, try:

```bash
./ampsci -i
./ampsci -i HartreeFock
```

which should output:

```java
// Available ampsci options/blocks
ampsci{
  // These are the top-level ampsci input blocks and options. Default values are
  // given in square brackets following the description: [default_value]. Blocks
  // end with '{}', options end with ';'. run `ampsci -i BlockName` (for any of
  // the following blocks) to see all the options available for that block.

  // Which atom to run for
  Atom{}
  // Set nuclear parameters
  Nucleus{}
  // Set radial grid/lattice parameters
  Grid{}
  // Options for solving atomic system
  HartreeFock{}
  // Options for the QED radiative potential (usually defaults suffice)
  RadPot{}
  // Include an extra effective potential. Rarely used.
  ExtraPotential{}
  // Basis of HF eigenstates used for MBPT
  Basis{}
  // Options for MBPT and correlation corrections
  Correlations{}
  // Like basis, but includes correlations. Used for sum-over-states
  Spectrum{}
  // Option for including `exotic` (e.g., muonic) atom states.
  Exotic{}
  // Configuration Interaction
  CI{}
  // Run any number of modules (* -> module name). `ampsci -m` to see available
  // modules
  Module::*{}
}
```

and

```java
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
  // HF convergance goal [1.0e-13]
  eps;
  // Method for mean-field approximation: HartreeFock, Hartree, KohnSham, Local
  // [HartreeFock]
  method;
  // Include Breit into HF? true/false, or scale factor. Scale factor for Breit
  // Hamiltonian is usially 0.0 (no Breit) or 1.0 (full Breit), but can take any
  // value. [false]
  Breit;
  // Include QED? Three options: true, false, valence. If 'valence, will include
  // QED only into valence states, but not the core. Detailed QED options are
  // set within the RadPot{} block - if that block is not set, defaults will be
  // used. By default, this option is false, unless the RadPot{} block exists,
  // in which case it is true
  QED;
}
```

respectively

## Basic example

* Here's one basic example
* Several other are provided with the code in `doc/examples/`

```java
Atom {
  Z = Cs;
  A = 133;
}

HartreeFock {
  core = [Xe];
  valence = 7sp;
}

Grid {
  r0 = 1e-6;
  rmax = 120.0;
  num_points = 2000;
}

// E1 matrix elements, including RPA (TDHF method, solved at w=0)
Module::matrixElements {
  operator = E1;
  rpa = TDHF;
  omega = 0.0;
}

// Calculate hyperfine matrix elements
Module::matrixElements {
  operator = hfs;
  off-diagonal = false; // only the diagonal matrix elements
  options{
      F = pointlike;
  }
}
```

The output should look something like:

```text
Running for Cs, Z=55 A=133
Fermi nucleus;  r_rms = 4.8041, c_hdr = 5.67073, t = 2.3
Loglinear (b=40) grid: 1e-06 -> 120.0, N=2000, du=0.432
========================================================
Hartree-Fock
Core   :  it: 28 eps=8.8e-14 for 5p-
Val    :  it: 38 eps=0.0e+00 for 6s+ [ 38 eps=0e+00 for 6s+]

Cs-133
Core: [Xe] V^N-1
#  nk   r_rms  Rinf    eps          En (au)         En (/cm)
0  1s+   0.03   0.6  7e-25  -1330.119222417   -291927426.008
1  2s+   0.12   1.6  4e-22   -212.564636008    -46652545.129
2  2p-   0.10   1.7  2e-22   -199.429666779    -43769752.599
3  2p+   0.11   1.7  3e-22   -186.436769887    -40918141.343
4  3s+   0.32   3.5  2e-20    -45.969772036    -10089198.771
5  3p-   0.31   3.7  2e-20    -40.448337479     -8877383.958
6  3p+   0.32   3.8  2e-20    -37.894342931     -8316846.946
7  3d-   0.29   4.4  1e-20    -28.309546081     -6213227.190
8  3d+   0.30   4.5  1e-20    -27.775202576     -6095952.346
9  4s+   0.74   7.7  3e-19     -9.512815190     -2087821.607
10 4p-   0.77   8.7  2e-19     -7.446281887     -1634269.972
11 4p+   0.80   9.0  2e-19     -6.920996585     -1518983.174
12 4d-   0.90  12.8  1e-19     -3.485619405      -765005.034
13 4d+   0.91  13.0  1e-19     -3.396901664      -745533.740
14 5s+   1.88  19.9  1e-19     -1.489800625      -326973.443
15 5p-   2.15  25.7  6e-20     -0.907895439      -199260.017
16 5p+   2.25  26.7  6e-20     -0.840336747      -184432.598
E_c = -7786.640563

Valence: CsI
#  nk    r_rms   Rinf    eps        En (au)       En (/cm)       En (/cm)
0  6s+    6.52   69.6  0e+00   -0.127368022     -27954.050          0.000
1  7s+   14.58  109.1  0e+00   -0.055187338     -12112.221      15841.829
2  6p-    8.65   86.1  0e+00   -0.085615790     -18790.494       9163.556
3  7p-   18.16  120.0  0e+00   -0.042021354      -9222.621      18731.429
4  6p+    8.86   87.0  0e+00   -0.083785371     -18388.763       9565.286
5  7p+   18.47  120.0  0e+00   -0.041368005      -9079.228      18874.822

--------------------------------------------------------------------------------
Module: Module::matrixElements

Matrix Elements - Operator: E1
Reduced matrix elements
Units: |e|aB
Including RPA: TDHF method
TDHF E1 (w=0.0000): 18 4.8e-11 [3p+,d-]

E1

   a    b    w_ab        t0_ab          +RPA 
 6p-  6s+    0.0417522  -5.277688e+00  -4.974363e+00
 7p-  6s+    0.0853467  -3.717428e-01  -2.387087e-01
 6p+  6s+    0.0435827   7.426436e+00   7.013022e+00
 7p+  6s+    0.0860000   6.947454e-01   5.087232e-01
 6p-  7s+   -0.0304285   4.413150e+00   4.449380e+00
 7p-  7s+    0.0131660  -1.100887e+01  -1.092105e+01
 6p+  7s+   -0.0285980  -6.671033e+00  -6.712242e+00
 7p+  7s+    0.0138193   1.534479e+01   1.522742e+01

matrixElements: T = 4.30 s

--------------------------------------------------------------------------------
Module: Module::matrixElements

Hyperfine structure: Cs, Z=55 A=133
K=1 (magnetic dipole)
Using pointlike nuclear distro for F(r)
w/ r_N = 0fm = 0au  (r_rms=0fm)
Points inside nucleus: 0
mu = 2.5778, I = 3.5, g = 0.736514

Matrix Elements - Operator: hfs1
Hyperfine constants (magnetic type), K=1
Units: MHz
Including RPA: TDHF method
TDHF hfs1 (w=0.0000): 28 3.3e-09 [4d-,s+]

hfs1

   a    b    w_ab        t0_ab          +RPA 
 6s+  6s+    0.0000000   1.431372e+03   1.725280e+03
 7s+  7s+    0.0000000   3.933064e+02   4.732411e+02
 6p-  6p-    0.0000000   1.607591e+02   2.012689e+02
 7p-  7p-    0.0000000   5.755954e+01   7.153063e+01
 6p+  6p+    0.0000000   2.387758e+01   4.276946e+01
 7p+  7p+    0.0000000   8.625730e+00   1.533316e+01

matrixElements: T = 4.88 s

ampsci: T = 9.40 s
```

See the [detailed tutorials](\ref tutorials) for more.
