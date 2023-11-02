# Advanced Tutorial: CI+MBPT

\brief Adanced ampsci tutorial: CI+MBPT for two-valence atoms

[[Home](/README.md)]

This assumes you already have ampsci compiled and have a basic understanding of how to run and use it.

* See [doc/tutorial.md](/doc/tutorial.md) for the basic tutorial
* See [ampsci.dev/ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full descrition of the physics

## Contents

* [CI+MBPT overview](#cimbpt)
* [Basic CI example](#basic)

## CI+MBPT overview: Configuration Interaction with Many-Body Perturbation Theory <a name="cimbpt"></a>

For an $M$-valence atomic system, the effective Hamiltonian is
$\begin{equation}
    H_{\rm CI} = H_{\rm CI}^1 + H_{\rm CI}^2 =
    \sum_i^M \left(h^{\rm HF}(r_i) + \Sigma^1(r_i)\right)
    +\sum_{i<j}\left(r^{-1}_{ij}
    +\Sigma^2(r_i,r_j)\right),
\end{equation}$
where $h^{\rm HF}$ is one-particle the Hartree-Fock Hamiltonian (with HF potential due to the $N-M$ core electrons), $\Sigma^1$ accounts for the core-valence correlations, and $\Sigma^2$ accounts for the screening of the valence-valence Coulomb interaction by the core electrons.

The CI routines in ampsci are only for two-valence systems: $M=2$.

In the CI method, approximate valence-space wavefunctions, $\Psi$, are expanded over $M$-particle wavefunctions called Configuration-State Functions (CSFs), $\psi_I$:
$\begin{equation}
    \left|{\Psi,J^\pi J_z}\right\rangle = \sum_I c_I  \left|{I,J^\pi J_z}\right\rangle.
\end{equation}$
The CSFs are combinations of Slater-determinants formed from single-particle eigenfunctions.
The CSFs are eigenfunctions of $J^2$, $J_z$, and parity ($\pi$).
For each $J^\pi$ symmetry, the energies and wavefunctions (expansion coefficients) are found by solving the Schr\"odinger equation, which for a finite set of $N_{CSF}$ CSFs, is cast to an $N_{CSF}^2$ eigenvalue problem:
$\begin{equation}
    \sum_J c_J\left\langle{I}\right|H_{\rm eff}\left|{J}\right\rangle = E c_I,
\end{equation}$
For the single-particle basis, we use eigenfuctions of the same $h^{\rm HF}$ Hamiltonian from the CI Hamiltonian. This is known as the $V^{N-M}$ approximation, with simplifies the MBPT part of the calculation, and is very accurate for two-valence systems.

-----------

## Basic CI example <a name="basic"></a>

As always, we check the available input options: `./ampsci -a CI`

`./ampsci -a CI`

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

For a simple example, we'll consider neutral Mg.

The first part of the imput file will be familiar from previous examples.
We don't need valence states in Hartree-Fock, so, this can be left blank.
We don't include the two valence states in the core, hence use $V^{N-2}$ approximation.

```java
Atom {
  Z = Mg;
}

HartreeFock {
  core = [Ne];
}

Basis {
  number = 30;
  order = 7;
  r0 = 1.0e-4; 
  rmax = 50.0;
  states = 30spdfghi;
}
```

To begin the CI calculation

```java
CI{
  ci_basis = 20spdf;
  J+ = 0,1,2; 
  J- = 0,1,2; 
  num_solutions = 3;
}
```

The `ci_basis` options means we will use two-particle CSFs formed from combinations of single-particle basis states up to $n=20$ for $s$, $p$, $d$, and $f$ states.
We will find the lowest 3 solutions for each $J^\pi$ symmetry up to $J=2$.

The first part of the output:

```text
  Using 20spdf = 124 orbitals in CI expansion
  Calculate two-body Coulomb integrals: Q^k_abcd
  For: 20spdf
  fill: T = 4.09 s
  Summary: 
  k=0: 674541 [712697]
  k=1: 3221991 [3439651]
  k=2: 5519503 [5967347]
  k=3: 3943836 [4026031]
  k=4: 2017036 [2144977]
  k=5: 421821 [444487]
  k=6: 97903 [99733]
  Total: 15896631 non-zero integrals
  Calculated 15896631 new Coulomb integrals
  Writing 15896631 integrals to file: Mg2.qk..
```

The code pre-computes all the two-body Coulomb integrals $Q^k_{ijkl}$.
This uses a significant memory footprint, but makes the subsequent calculations much faster.
It only calculates the required integrals. Even though our total basis was up to `30spdfghi`, we only used `20spdf` in the CI expansion (the remaining basis states will be used for MBPT in the next example).

Then, the code runs CI routine for each symmetry.
The output for each symmetry will look something like this:

```text
  Run CI for J=2, even parity
  Total CSFs: 3322
  Find first 3 solutions
  Eigenvalues: T = 6.35 s
  2 + 0   -0.61245597 au   -134418.55 cm^-1         0.00 cm^-1
    3s+3d- 24.315%
    3s+4d- 6.641%
    3s+3d+ 36.577%
    3s+4d+ 9.997%
    3p-3p+ 13.492%
    3p+^2  6.672%
    --------------
    gJ = 0.999992
    3s3d   1^D_2

  2 + 1   -0.60402693 au   -132568.59 cm^-1      1849.96 cm^-1
    3s+3d- 36.479%
    3s+4d- 22.333%
    3s+3d+ 24.246%
    3s+4d+ 14.846%
    --------------
    gJ = 1.16666
    3s3d   3^D_2

  2 + 2   -0.58093715 au   -127500.97 cm^-1      6917.58 cm^-1
    3s+4d- 17.604%
    3s+5d- 16.042%
    3s+3d+ 1.168%
    3s+4d+ 26.463%
    3s+5d+ 24.122%
    3p-3p+ 8.217%
    3p+^2  3.988%
    --------------
    gJ = 0.999999
    3s4d   1^D_2
```

After each line, the code will list all CSFs that contribute to the solutiom at above the 1% level (i.e., with $c_I^2>0.01$).
The code also calculates the $g$-factor (without RPAd), which are useful for level identification.
The term symbol (${}^{2S+1}L_J$ e.g., 1^D_2) and leading configuration are given -- the term symbol is 'guessed' based on the g-factor. It is not well defined relativistically, so should be seen as indicative only.

The final output shows a summary of the calculation:

```text
  J  pi  #  conf.  %   Term   Energy(au)   Energy(/cm)   Level(/cm)  gJ
  0  +1  0  3s^2   89  1S    -0.81804187    -179539.44         0.00
  0  +1  1  3s4s   66  1S    -0.62364477    -136874.21     42665.23
  0  +1  2  3s6s   61  1S    -0.58273524    -127895.60     51643.83
  1  +1  0  3s4s   83  3S    -0.63395079    -139136.12     40403.32  2.0000
  1  +1  1  3s3d   61  3D    -0.60402720    -132568.65     46970.79  0.5000
  1  +1  2  3s5s   47  3S    -0.58579051    -128566.16     50973.28  2.0000
  2  +1  0  3s3d   61  1D    -0.61245597    -134418.55     45120.89  1.0000
  2  +1  1  3s3d   61  3D    -0.60402693    -132568.59     46970.85  1.1667
  2  +1  2  3s4d   44  1D    -0.58093715    -127500.97     52038.47  1.0000
  0  -1  0  3s3p   94  3Po   -0.72277949    -158631.76     20907.67
  0  -1  1  3s4p   65  3Po   -0.60427214    -132622.40     46917.03
  0  -1  2  3s6p   64  3Po   -0.57513793    -126228.19     53311.25
  1  -1  0  3s3p   94  3Po   -0.72268645    -158611.34     20928.09  1.5000
  1  -1  1  3s3p   71  1Po   -0.66089031    -145048.66     34490.78  1.0000
  1  -1  2  3s4p   65  3Po   -0.60425669    -132619.01     46920.42  1.5000
  2  -1  0  3s3p   94  3Po   -0.72249992    -158570.40     20969.03  1.5000
  2  -1  1  3s4p   65  3Po   -0.60422554    -132612.18     46927.26  1.5000
  2  -1  2  3s6p   64  3Po   -0.57511959    -126224.16     53315.28  1.5000
```

Where `conf` is the leading configuration (in non-relativistic notation), and the `%` column shows the combined contrutions from each relativistic configuration with the same non-relativistic configuration
(e.g., `3s3d` may have contributions from `3s3d-` and `3s3d+`).

| Level  |             | AMPSCI  | Exp.    | $\Delta$ |
|--------|-------------|---------|---------|----------|
| $3s^2$ | ${}^1S_0$   | -179539 | -182939 | -1.9%   |
| $3s4s$ | ${}^1S_0$   | 42665   | 43503   | -1.9%   |
| $3s6s$ | ${}^1S_0$   | 51644   | 52556   | -1.7%   |
| $3s4s$ | ${}^3S_1$   | 40403   | 41197   | -1.9%   |
| $3s3d$ | ${}^3D_1$   | 46971   | 47957   | -2.1%   |
| $3s5s$ | ${}^3S_1$   | 50973   | 51873   | -1.7%   |
| $3s3d$ | ${}^1D_2$   | 45121   | 46403   | -2.8%   |
| $3s3d$ | ${}^3D_2$   | 46971   | 47957   | -2.1%   |
| $3s4d$ | ${}^1D_2$   | 52038   | 53135   | -2.1%   |
| $3s3p$ | ${}^3P^o_0$ | 20908   | 21850   | -4.3%   |
| $3s4p$ | ${}^3P^o_0$ | 46917   | 47841   | -1.9%   |
| $3s6p$ | ${}^3P^o_0$ | 53311   | 54249   | -1.7%   |
| $3s3p$ | ${}^3P^o_1$ | 20928   | 21870   | -4.3%   |
| $3s3p$ | ${}^1P^o_1$ | 34491   | 35051   | -1.6%   |
| $3s4p$ | ${}^3P^o_1$ | 46920   | 47844   | -1.9%   |
| $3s3p$ | ${}^3P^o_2$ | 20969   | 21911   | -4.3%   |
| $3s4p$ | ${}^3P^o_2$ | 46927   | 47851   | -1.9%   |
| $3s6p$ | ${}^3P^o_2$ | 53315   | 54253   | -1.7%   |

The table shows the results of the calculation, and comparison to experimental excitation energies, in units of ${\rm cm}^{-1}$ (for the ground state, the ionisation potential is instead shown).
The agreement is at the ~few % level.

-----------

## CI+MBPT example <a name="basic"></a>

Here, we will improve the accuracy of the calculation by including the MBPT corrections.
We do this by setting `sigma1` and `sigma2` options to true.

```java
CI{
  ci_basis = 20spdf;
  J+ = 0,1,2;
  J- = 0,1,2;
  num_solutions = 3;
  sigma1 = true;
  sigma2 = true;
}
```

This will use the full basis (from `Basis{}`) for the MBPT part. In this case, it was `30spdfghi`.
(Different subsets of the full basis can be used for the internal lines of the $\Sigma^1$ and $\Sigma^2$ diagrams using the `s1_basis` and `s2_basis` options, respectively).

The code will now calculate more $Q^k$ Coulomb integrals, since they are required in the MBPT calculations. It will read in the existing file (if there is one), so it doesn't need to start from scratch, and only calculate the missing integrals:

```text
  Calculate two-body Coulomb integrals: Q^k_abcd
  Read 15896631 integrals from file: Mg2.qk
  For: 20spdf
  fill: T = 947.11 ms
  and: 30spdfghi
  fill: T = 1.07 mins
  Summary: 
  k=0: 1235215 [1447153]
  k=1: 6337133 [6984629]
  k=2: 10584749 [12117689]
  k=3: 9629901 [12117689]
  k=4: 7888486 [12117689]
  k=5: 5720196 [11200489]
  k=6: 3267217 [6456007]
  k=7: 997236 [1056323]
  k=8: 0 [2]
  Total: 45660133 non-zero integrals
  Calculated 29763502 new Coulomb integrals
  Writing 45660133 integrals to file: Mg2.qk..
```

This uses close to 5Gb of memory, so may already be difficuly on some laptops.
Then, is will calculate the matrix elements of the two-body $\Sigma^2$ operator.
This is the slow part of the calculation.
Fortunately, the $\Sigma^2$ correction is rather small, and good accuracy can be obtained by only including matrix elements between the lowest few valence-space basis states.
By default, it will include up to $n=n_{\rm core}+3$, where $n_{\rm core}$ is the largest $n$ in the core, for each $l$ in the `ci_basis`.
This can be controlled manually with the `cis2_basis` option.

```text
  Calculate two-body MBPT integrals: Sigma^k_abcd
  For: 5spdf, using 30spdfghi
  Count non-zero: 0.89 ms
  Reserve: 0.17 ms
  Fill w/ zeros: 4.61 ms
  Fill w/ values: 8.49 s
  Summary: 
  k=0: 5151 [5503]
  k=1: 32131 [33223]
  k=2: 45451 [45481]
  k=3: 32131 [33223]
  k=4: 12403 [12983]
  k=5: 2701 [2753]
  k=6: 300 [313]
  k=7: 10 [11]
  Total: 130278 non-zero integrals
  fill: T = 8.50 s
  Calculated 130278 new MBPT integrals
  Writing 130278 integrals to file: Mg2_1_30spdfghi_8.sk..
```

The `For: 5spdf, using 30spdfghi` means the external lines in $\Sigma_{ijkl}$ include up to `5spdf`, while the internal lines use the full `30spdfghi` basis.

This leads to a significant improvement in the accuracy:

| Level  |             | AMPSCI  | Exp.    | $\Delta$ |
|--------|-------------|---------|---------|----------|
| $3s^2$ | ${}^1S_0$   | -182804 | -182939 | -0.07%   |
| $3s4s$ | ${}^1S_0$   | 43490   | 43503   | -0.03%   |
| $3s6s$ | ${}^1S_0$   | 52526   | 52556   | -0.06%   |
| $3s4s$ | ${}^3S_1$   | 41208   | 41197   | 0.03%    |
| $3s3d$ | ${}^3D_1$   | 47957   | 47957   | 0.00%    |
| $3s5s$ | ${}^3S_1$   | 51875   | 51873   | 0.00%    |
| $3s3d$ | ${}^1D_2$   | 46385   | 46403   | -0.04%   |
| $3s3d$ | ${}^3D_2$   | 47960   | 47957   | 0.01%    |
| $3s4d$ | ${}^1D_2$   | 53113   | 53135   | -0.04%   |
| $3s3p$ | ${}^3P^o_0$ | 21822   | 21850   | -0.13%   |
| $3s4p$ | ${}^3P^o_0$ | 47837   | 47841   | -0.01%   |
| $3s6p$ | ${}^3P^o_0$ | 54271   | 54249   | 0.04%    |
| $3s3p$ | ${}^3P^o_1$ | 21844   | 21870   | -0.12%   |
| $3s3p$ | ${}^1P^o_1$ | 35090   | 35051   | 0.11%    |
| $3s4p$ | ${}^3P^o_1$ | 47840   | 47844   | -0.01%   |
| $3s3p$ | ${}^3P^o_2$ | 21876   | 21911   | -0.16%   |
| $3s4p$ | ${}^3P^o_2$ | 47845   | 47851   | -0.01%   |
| $3s6p$ | ${}^3P^o_2$ | 54272   | 54253   | 0.03%    |

The discrepancies are now at the level of 0.1% or below.

On my pc, this entire calculation took less than two minutes.

-----------

## Matrix elements <a name="matrix"></a>

This takes very similar options to the regular `MatrixElements{}` module:

```java
// Available Module::CI_matrixElements options/blocks
Module::CI_matrixElements{
  // e.g., E1, hfs (see ampsci -o for available operators)
  operator;
  // options specific to operator
  options{}
  // Method used for RPA: true(=TDHF), false, TDHF, basis, diagram
  rpa;
  // Text or number. Freq. for RPA (and freq. dependent operators). Put 'each'
  // to solve at correct frequency for each transition. [0.0]
  omega;
  // List of angular momentum Js to calculate matrix elements for. If blank, all
  // available Js will be calculated. Must be integers (two-electron only).
  J;
  // As above, but for EVEN CSFs only (takes precedence over J).
  J+;
  // As above, but for ODD CSFs (takes precedence over J).
  J-;
  // Maximum solution number to calculate MEs for. If blank, will calculate all.
  num_solutions;
}
```

It's typically not recommended to use `omega=each`, as this solves the RPA equations for each transition, and the frequency-dependence is usually small.

```java
Module::CI_MatrixElements{
  operator = E1;
  rpa = true;
  omega = 0.15970531;
}
```

The output is something like:

```text
  CI Matrix Elements -  Operator: E1
  Units: |e|aB
  Including RPA: TDHF method
  Solving RPA at fixed frequency: w=0.159705
  TDHF E1 (w=0.1597): 13 7.5e-09 [1s+,p-]

  Ja  # conf     - Jb  # conf      w_ab     t_ab
  0+  0 3s^2 1S  - 1-  0 3s3p 3Po  0.09949 -1.23847e-03
  0+  0 3s^2 1S  - 1-  1 3s3p 1Po  0.15989 -4.02693e+00
  0+  1 3s4s 1S  - 1-  0 3s3p 3Po  0.09865 -1.68797e-03
  0+  1 3s4s 1S  - 1-  1 3s3p 1Po  0.03825 -4.24587e+00
  1-  0 3s3p 3Po - 0+  0 3s^2 1S   0.09949  1.23847e-03
  1-  0 3s3p 3Po - 0+  1 3s4s 1S   0.09865  1.68797e-03
  1-  1 3s3p 1Po - 0+  0 3s^2 1S   0.15989  4.02693e+00
  1-  1 3s3p 1Po - 0+  1 3s4s 1S   0.03825  4.24587e+00
```

(this is for a reduced set of levels).

This implies the lifetime of the first ${}^1P^o_1$ level is 2.12 ns, in excellent agreement with experiment.
