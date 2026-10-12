\page tutorial_ci Two-valence atoms: CI+MBPT

\brief CI+MBPT calculations for two-valence atoms

This assumes you already have ampsci compiled and have a basic understanding of how to run and use it.

* See [Compilation](\ref compilation) for compilation instructions
  * For CI calculations, linking against a threaded BLAS (OpenBLAS or MKL) is recommended for faster matrix diagonalisation -- see the compilation page for details
* See [Basic Tutorial](\ref tutorial_basic) for getting started if unfamiliar
* See [ampsci.dev/ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full description of the physics

## Contents

* [CI+MBPT overview](#cimbpt)
* [Basic CI example](#basic)

## CI+MBPT overview: Configuration Interaction with Many-Body Perturbation Theory <a name="cimbpt"></a>

For an \f$M\f$-valence atomic system, the effective Hamiltonian is

\f[
    H_{\rm CI} = H_{\rm CI}^1 + H_{\rm CI}^2 =
    \sum_i^M \left(h^{\rm HF}(r_i) + \Sigma^1(r_i)\right)
    +\sum_{i < j}\left(r^{-1}_{ij}
    +\Sigma^2(r_i,r_j)\right),
\f]

where \f$h^{\rm HF}\f$ is one-particle the Hartree-Fock Hamiltonian (with HF potential due to the \f$N-M\f$ core electrons), \f$\Sigma^1\f$ accounts for the core-valence correlations, and \f$\Sigma^2\f$ accounts for the screening of the valence-valence Coulomb interaction by the core electrons.

The CI routines in ampsci are only for two-valence systems: \f$M=2\f$.

In the CI method, approximate valence-space wavefunctions, \f$\Psi\f$, are expanded over \f$M\f$-particle wavefunctions called Configuration-State Functions (CSFs), \f$\psi_I\f$:

\f[
    \left|{\Psi,J^\pi J_z}\right\rangle = \sum_I c_I  \left|{I,J^\pi J_z}\right\rangle.
\f]

The CSFs are combinations of Slater-determinants formed from single-particle eigenfunctions.
The CSFs are eigenfunctions of \f$J^2\f$, \f$J_z\f$, and parity (\f$\pi\f$).
For each \f$J^\pi\f$ symmetry, the energies and wavefunctions (expansion coefficients) are found by solving the Schr\"odinger equation, which for a finite set of \f$N_{CSF}\f$ CSFs, is cast to an \f$N_{CSF}^2\f$ eigenvalue problem:

\f[
    \sum_J c_J\left\langle{I}\right|H_{\rm eff}\left|{J}\right\rangle = E c_I,
\f]

For the single-particle basis, we use eigenfunctions of the same \f$h^{\rm HF}\f$ Hamiltonian from the CI Hamiltonian. This is known as the \f$V^{N-M}\f$ approximation, with simplifies the MBPT part of the calculation, and is very accurate for two-valence systems.

-----

## Basic CI example <a name="basic"></a>

As always, we check the available input options: `./ampsci -i CI`

`./ampsci -i CI`

```java
// Available CI options/blocks
CI{
  ci_basis;
    // Basis used for CI expansion; must be a sub-set of
    // full ampsci basis [default: 20spdf]
  J;
    // List of total angular momentum J for CI solutions
    // (comma separated). Must be integers (two-electron
    // only). []
  J+;
    // As above, but for EVEN CSFs only (takes precedence
    // over J).
  J-;
    // As above, but for ODD CSFs (takes precedence over J).
  num_solutions;
    // Number of CI solutions to find (for each J/pi) [5]
  all_below_cm;
    // Find all CI solutions for energies below this
    // threshold, in inverse cm. Note that this is the total
    // energy, not the excitation energy. If set,
    // num_solutions is ignored.
  sigma1;
    // Include one-body MBPT correlations? [false]
  sigma2;
    // Include two-body MBPT correlations? [false]
  Brueckner;
    // Use Brueckner (spectrum) states for CI basis? Must
    // have Correlations, Spectrum, and sigma1. [false]
  cis2_basis;
    // The subset of ci_basis for which the two-body MBPT
    // corrections are calculated. Must be a subset of
    // ci_basis. If existing sk file has more integrals,
    // they will be used. [default: Nspdf, where N is
    // maximum n for core + 3]
  Breit2;
    // Include two-body Breit? Default is true if Breit
    // included in HF. Ignored if Breit not included in HF.
    // [true]
  Breit_basis;
    // Subset of ci_basis used to include two-body Breit
    // corrections into CI matrix. Large basis is slow, uses
    // huge memory, and makes small contribution. [default:
    // Nspdf, where N is maximum n for core + 6]
  s1_basis;
    // Usually should be left as default. Basis used for the
    // one-body MBPT diagrams (Sigma^1) internal lines.
    // These are the most important, so in general the
    // default (all basis states) should be used. Must be a
    // subset of full ampsci basis. [default: full basis]
    //  - Note: if CorrelationPotential is available, it
    // will be used instead of calculating the Sigma_1
    // integrals
  s2_basis;
    // Usually should be left blank. Basis used for internal
    // lines of the two-body MBPT diagrams (Sigma^2)
    // internal lines. Must be a subset of s1_basis.
    // [default: s1_basis]
  n_min_core;
    // Minimum n for core to be included in MBPT [1]
  max_k;
    // Maximum k (multipolarity) to include when calculating
    // new Coulomb integrals. Higher k often contribute
    // negligably. Note: if qk file already has higher-k
    // terms, they will be included. Set negative (or very
    // large) to include all k. [8]
  denominators;
    // 'DFK', 'BW', 'RS', 'Fermi', 'Fermi0'. Denominators
    // used in Sigma2 matrix elements. DFK
    // (Dzuba-Flambaum-Kozlov): target-state legs use the
    // lowest excited state of their kappa,
    // intermediate-state legs use actual energies. BW
    // (Brillouin-Wigner): the denominator is E0 minus the
    // energy of the intermediate state, where E0 is the
    // total valence energy of the target level; DFK is this
    // with E0 taken from the leading configuration. BW
    // requires the E0 of each (J, parity) block, so it
    // enables iterative_correction (requires sigma1): S^k
    // and dS^k/dE0 are tabulated at an internal reference
    // and each S^k is resummed to the E0 of the block being
    // solved. This doubles the Sigma_2 calculation time and
    // stores a second sk file. RS uses actual energies for
    // all external legs, Fermi uses the lowest excited
    // state for each kappa (both legs), Fermi0 uses lowest
    // excited state for all kappas (and thus cancels in all
    // except diagram 'd'). Applies to Sigma_2 only. [DFK]
  iterative_correction;
    // Include the derivative (dSigma/dE) correction to the
    // Sigma_1 matrix elements. Restores the configuration
    // dependence lost by evaluating Sigma_1 at a fixed
    // energy. The reference energy E0 is the lowest energy
    // of each J/pi, and is found iteratively. With
    // denominators = BW, the same E0 also shifts the
    // Sigma_2 matrix elements (see denominators). Prefer
    // calculating dSigma/dE in the Correlations block
    // (option derivative) which is much faster, and can
    // include all-orders [false]
  qk_file;
    // Filename for storing two-body Coulomb integrals. By
    // default, is At.qk.abf, where At is atomic symbol +
    // 'identity'. Set to 'false' to disable read/write.
  sk_file;
    // Filename for storing two-body Sigma_2 integrals. By
    // default, is At_b_hash.sk.abf, where At is atomic
    // symbol + 'identity', b is the s2 (internal) basis,
    // and hash encodes other settings that changes the
    // integrals. Set to 'false' to disable read/write.
  dsk_file;
    // As sk_file, but for the second Sigma_2 table used
    // when denominators = BW. Default is At_b_hash.dsk.abf.
    // nb: it holds S^k at the shifted E0; the derivative is
    // formed from it after reading.
  ds1_file;
    // Filename for storing the Sigma_1
    // derivative-correction tables (iterative_correction).
    // Default: At_b_hash.ds1.abf, where b is the s1
    // (internal) basis, and hash encodes other settings
    // that change the tables. Not used if the Correlations
    // block provides dSigma/dE (option derivative). Set to
    // 'false' to disable read/write.
  bk_file;
    // Filename for storing two-body Breit integrals. By
    // default, is At.bk.abf, where At is atomic symbol +
    // 'identity'. Set to 'false' to disable read/write.
  ci_file;
    // Filename for storing CI solutions. Default:
    // At_b_hash.ci.abf, where At is atomic symbol +
    // 'identity', b is ci_basis, and hash encodes other
    // setting that  changes the CI Hamiltonian. Set to
    // 'false' to disable read/write.
  read_only;
    // If true, will *only* read in the existing CI
    // solutions, and will not calculate anything new, even
    // if new states are requested. You must set the
    // ci_basis (not read in) and the J/pi solutions you
    // want; only these will be read in. Use if you know CI
    // has already been completed. [false]
  no_new_integrals;
    // Usually false. If set to true, ampsci will not
    // calculate any new Coulomb or Sigma_2 integrals, even
    // if they are implied by the above settings. This saves
    // time when we know all required integrals already
    // exist, since the code doesn't need to check. [false]
  sort_output;
    // Sort output by energy? Default is to sort by J and Pi
    // first. [false]
  print_details;
    // Condition to print details of each CI solution
    // (otherwise just prints summary) [true]
  fk;
    // Effective QPQ ~ fk Q screening factors for the
    // Sigma_2 Coulomb lines. Give an explicit list to use
    // those values; set 'false' for no screening (or 'fk =
    // 1;'). If blank (default), takes from correlation
    // potential if it has them.
  extrapolate_sigma2;
    // Extrapolate Sigma_2 to diagrams outside cis2_basis,
    // using average correction ratios: S^k ~ h_k*Q^k, where
    // h_k = <S^k/Q^k> is averaged over the calculated
    // Sigma_2 integrals. Note: These are stored in the Sk
    // table, but NOT written to disk [true]
  exclude_wrong_parity_box;
    // Excludes the Sigma_2 box corrections that have
    // 'wrong' parity when calculating Sigma2 matrix
    // elements (i.e., that don't match the Q^k selection
    // rules). Note: If existing sk file already has these,
    // they will be included [false]
}
```

For a simple example, we'll consider neutral Mg.

The first part of the input file will be familiar from previous examples.
We don't need valence states in Hartree-Fock, so, this can be left blank.
We don't include the two valence states in the core, hence use \f$V^{N-2}\f$ approximation.

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

The `ci_basis` options means we will use two-particle CSFs formed from combinations of single-particle basis states up to \f$n=20\f$ for \f$s\f$, \f$p\f$, \f$d\f$, and \f$f\f$ states.
We will find the lowest 3 solutions for each \f$J^\pi\f$ symmetry up to \f$J=2\f$.

The first part of the output:

```text
  Using 20spdf = 124 orbitals in CI expansion
  Calculate two-body Coulomb integrals: Q^k_abcd
  For: 20spdf (for CI)
  Stored non-zero two-body integrals (Qk symmetry):
  k =  0: 674541
  k =  1: 3221991
  k =  2: 5519503
  k =  3: 3943836
  k =  4: 2017036
  k =  5: 421821
  k =  6: 97903
  Total: 15896631
  Calculated 15896631 new Coulomb integrals
  Writing 15896631 integrals to file: Mg2.qk.abf..
```

The code pre-computes all the two-body Coulomb integrals \f$Q^k_{ijkl}\f$.
This uses a significant memory footprint, but makes the subsequent calculations much faster.
It only calculates the required integrals. Even though our total basis was up to `30spdfghi`, we only used `20spdf` in the CI expansion (the remaining basis states will be used for MBPT in the next example).

Then, the code runs CI routine for each symmetry.
The output for each symmetry will look something like this:

```text
  CI: J=2, even parity
  Total CSFs: 3322
  Find first 3 solutions
  3 eigenvalues: T = 472.56 ms
  472.57 ms
  Writing J=2 pi=+ (3 solutions) to Mg2_20spdf_xovgxe.ci.abf

  2 + 0   -0.61243636 au   -134414.24 cm^-1         0.00 cm^-1
     3s+3d- 24.314%
     3s+4d- 6.633%
     3s+3d+ 36.576%
     3s+4d+ 9.986%
     3p-3p+ 13.506%
     3p+^2  6.679%
     --------------
     gJ = 0.999992
     L_eff = 2.0000, S_eff = 0.0000
     3s3d   1^D_2

  2 + 1   -0.60399682 au   -132561.98 cm^-1      1852.26 cm^-1
     3s+3d- 36.482%
     3s+4d- 22.331%
     3s+3d+ 24.247%
     3s+4d+ 14.844%
     --------------
     gJ = 1.16666
     L_eff = 2.0000, S_eff = 1.0000
     3s3d   3^D_2

  2 + 2   -0.58091120 au   -127495.27 cm^-1      6918.97 cm^-1
     3s+4d- 17.618%
     3s+5d- 16.030%
     3s+4d+ 26.483%
     3s+5d+ 24.106%
     3p-3p+ 8.221%
     --------------
     gJ = 0.999999
     L_eff = 2.0000, S_eff = 0.0000
     3s4d   1^D_2
```

After each line, the code will list all CSFs that contribute to the solution at above the 1% level (i.e., with \f$c_I^2>0.01\f$).
The code also calculates the \f$g\f$-factor (without RPAd), which are useful for level identification.
The term symbol (\f${}^{2S+1}L_J\f$ e.g., 1^D_2) and leading configuration are given -- the term symbol is determined from the expectation values \f$\langle L^2\rangle\f$ and \f$\langle S^2\rangle\f$ (evaluated in the non-relativistic limit, and printed as \f$L_{\rm eff}\f$ and \f$S_{\rm eff}\f$, where \f$L_{\rm eff}(L_{\rm eff}+1)=\langle L^2\rangle\f$). It is not well defined relativistically, so should be seen as indicative only.

The final output shows a summary of the calculation:

```text
  J  pi  #  conf.  %   Term   Energy(au)   Energy(/cm)   Level(/cm)  gJ
  0  +1  0  3s^2   89  1S    -0.81801329    -179533.17         0.00
  0  +1  1  3s4s   66  1S    -0.62361589    -136867.87     42665.30
  0  +1  2  3s6s   61  1S    -0.58270456    -127888.87     51644.30
  1  +1  0  3s4s   83  3S    -0.63391442    -139128.13     40405.03  2.0000
  1  +1  1  3s3d   61  3D    -0.60399697    -132562.01     46971.15  0.5000
  1  +1  2  3s5s   47  3S    -0.58575697    -128558.80     50974.37  2.0000
  2  +1  0  3s3d   61  1D    -0.61243636    -134414.24     45118.92  1.0000
  2  +1  1  3s3d   61  3D    -0.60399682    -132561.98     46971.19  1.1667
  2  +1  2  3s4d   44  1D    -0.58091120    -127495.27     52037.90  1.0000
  0  -1  0  3s3p   94  3Po   -0.72276324    -158628.20     20904.97
  0  -1  1  3s4p   65  3Po   -0.60424210    -132615.81     46917.35
  0  -1  2  3s6p   64  3Po   -0.57510650    -126221.29     53311.88
  1  -1  0  3s3p   94  3Po   -0.72267118    -158607.99     20925.18  1.5000
  1  -1  1  3s3p   71  1Po   -0.66085861    -145041.70     34491.47  1.0000
  1  -1  2  3s4p   65  3Po   -0.60422689    -132612.47     46920.69  1.5000
  2  -1  0  3s3p   94  3Po   -0.72248658    -158567.48     20965.69  1.5000
  2  -1  1  3s4p   65  3Po   -0.60419623    -132605.74     46927.42  1.5000
  2  -1  2  3s6p   64  3Po   -0.57508849    -126217.33     53315.83  1.5000
```

Where `conf` is the leading configuration (in non-relativistic notation), and the `%` column shows the combined constructions from each relativistic configuration with the same non-relativistic configuration
(e.g., `3s3d` may have contributions from `3s3d-` and `3s3d+`).

| Level  |             | AMPSCI  | Exp.    | \f$\Delta\f$ |
|-----|-----|-----|-----|-----|
| \f$3s^2\f$ | \f${}^1S_0\f$   | -179533 | -182939 | -1.9%   |
| \f$3s4s\f$ | \f${}^1S_0\f$   | 42665   | 43503   | -1.9%   |
| \f$3s6s\f$ | \f${}^1S_0\f$   | 51644   | 52556   | -1.7%   |
| \f$3s4s\f$ | \f${}^3S_1\f$   | 40405   | 41197   | -1.9%   |
| \f$3s3d\f$ | \f${}^3D_1\f$   | 46971   | 47957   | -2.1%   |
| \f$3s5s\f$ | \f${}^3S_1\f$   | 50974   | 51873   | -1.7%   |
| \f$3s3d\f$ | \f${}^1D_2\f$   | 45119   | 46403   | -2.8%   |
| \f$3s3d\f$ | \f${}^3D_2\f$   | 46971   | 47957   | -2.1%   |
| \f$3s4d\f$ | \f${}^1D_2\f$   | 52038   | 53135   | -2.1%   |
| \f$3s3p\f$ | \f${}^3P^o_0\f$ | 20905   | 21850   | -4.3%   |
| \f$3s4p\f$ | \f${}^3P^o_0\f$ | 46917   | 47841   | -1.9%   |
| \f$3s6p\f$ | \f${}^3P^o_0\f$ | 53312   | 54249   | -1.7%   |
| \f$3s3p\f$ | \f${}^3P^o_1\f$ | 20925   | 21870   | -4.3%   |
| \f$3s3p\f$ | \f${}^1P^o_1\f$ | 34491   | 35051   | -1.6%   |
| \f$3s4p\f$ | \f${}^3P^o_1\f$ | 46921   | 47844   | -1.9%   |
| \f$3s3p\f$ | \f${}^3P^o_2\f$ | 20966   | 21911   | -4.3%   |
| \f$3s4p\f$ | \f${}^3P^o_2\f$ | 46927   | 47851   | -1.9%   |
| \f$3s6p\f$ | \f${}^3P^o_2\f$ | 53316   | 54253   | -1.7%   |

The table shows the results of the calculation, and comparison to experimental excitation energies, in units of \f${\rm cm}^{-1}\f$ (for the ground state, the ionisation potential is instead shown).
The agreement is at the ~few % level.

-----

## Storing and reusing CI solutions <a name="cifile"></a>

CI solutions (energies and eigenvectors) can be saved to disk and read back on subsequent runs, which avoids re-diagonalising the Hamiltonian.
By default, a filename is generated automatically based on the atomic identity, CI basis, and a hash of the other settings that change the CI Hamiltonian (e.g., `Mg2_20spdf_u05ndb.ci.abf`).
This can be overridden with the `ci_file` option; set it to `false` to disable read/write entirely.

On the first run the solutions are computed and written to file.
On subsequent runs ampsci reads them back and skips diagonalisation.
If more solutions are requested than are stored, ampsci re-solves (from scratch) and overwrites those objects in the file.

If the CI solutions already exist on disk and you do not need to recalculate (e.g., you are only computing matrix elements), set `read_only = true`.
This skips all integral calculations and reads directly from the file:

```java
CI{
  ci_basis = 20spdf;
  J+ = 0,1,2;
  J- = 0,1,2;
  read_only = true;
}
```

* The ci_basis _must_ be set, and must be the same.
* Code will throw an error/warning if the ci_basis does not match that expected
* The required J must be set: only the requested J's will be read in.
  * They can be a subset of previous run

### Multiple runs

The same CI output file stores solutions for each J/parity independently.
For example, if you required 20 `J=1-` levels, but only a few from the other states, you could run twice like:

```java
CI{
  ci_basis = 20spdf;
  J- = 1;
  num_solutions = 20;
}
```

```java
CI{
  ci_basis = 20spdf;
  J+ = 0,1,2;
  J- = 0,1,2;
  num_solutions = 3;
}
```

* In the second run, 20 previous `J=1-` solutions will be read in; there are no existing solutions for the other states, so the requested 3 will be calculated.

* These have to be completely seperate ampsci runs - cannot have two CI blocks in the same ampsci input file.

-----

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
(Different subsets of the full basis can be used for the internal lines of the \f$\Sigma^1\f$ and \f$\Sigma^2\f$ diagrams using the `s1_basis` and `s2_basis` options, respectively).

The code will now calculate more \f$Q^k\f$ Coulomb integrals, since they are required in the MBPT calculations. It will read in the existing file (if there is one), so it doesn't need to start from scratch, and only calculate the missing integrals:

```text
  Calculate two-body Coulomb integrals: Q^k_abcd
  Reading: Mg2.qk.abf
  Read 15896631 integrals from file: Mg2.qk.abf
  For: 20spdf (for CI)
  and: 30spdfghi (for MBPT)
  Stored non-zero two-body integrals (Qk symmetry):
  k =  0: 1235215
  k =  1: 6337133
  k =  2: 10584749
  k =  3: 9629901
  k =  4: 7888486
  k =  5: 5720196
  k =  6: 3267217
  k =  7: 997236
  k =  8: 0
  Total: 45660133
  Calculated 29763502 new Coulomb integrals
  Writing 45660133 integrals to file: Mg2.qk.abf..
```

This uses close to 5Gb of memory, so may already be difficulty on some laptops.
Then, is will calculate the matrix elements of the two-body \f$\Sigma^2\f$ operator.
This is the slow part of the calculation.
Fortunately, the \f$\Sigma^2\f$ correction is rather small, and good accuracy can be obtained by only including matrix elements between the lowest few valence-space basis states.
By default, it will include up to \f$n=n_{\rm core}+3\f$, where \f$n_{\rm core}\f$ is the largest \f$n\f$ in the core, for each \f$l\f$ in the `ci_basis`.
This can be controlled manually with the `cis2_basis` option.

```text
  Calculate two-body MBPT integrals: Sigma^k_abcd
  For: 5spdf, using 30spdfghi
  Calculating 130278 integrals
  0%, 10%, 21%, 31%, 36%, 47%, 57%, 63%, 73%, 84%, 94%, 100%
  3.59 s
  Calculated 130278 new MBPT integrals
  Writing 130278 integrals to file: Mg2_30spdfghi_gqlh2r.sk.abf..
```

The `For: 5spdf, using 30spdfghi` means the external lines in \f$\Sigma_{ijkl}\f$ include up to `5spdf`, while the internal lines use the full `30spdfghi` basis. This is usually the bottle-neck; these integrals take a very long time to compute, though typically make only a modest difference to the calculations.

The combined MBPT corrections lead to a significant improvement in the accuracy:

| Level      |                 | AMPSCI  | Exp.    | \f$\Delta\f$ |
|------------|-----------------|---------|---------|--------------|
| \f$3s^2\f$ | \f${}^1S_0\f$   | -182703 | -182939 | -0.13%       |
| \f$3s4s\f$ | \f${}^1S_0\f$   | 43440   | 43503   | -0.14%       |
| \f$3s6s\f$ | \f${}^1S_0\f$   | 52473   | 52556   | -0.16%       |
| \f$3s4s\f$ | \f${}^3S_1\f$   | 41146   | 41197   | -0.12%       |
| \f$3s3d\f$ | \f${}^3D_1\f$   | 47895   | 47957   | -0.13%       |
| \f$3s5s\f$ | \f${}^3S_1\f$   | 51810   | 51873   | -0.12%       |
| \f$3s3d\f$ | \f${}^1D_2\f$   | 46326   | 46403   | -0.17%       |
| \f$3s3d\f$ | \f${}^3D_2\f$   | 47895   | 47957   | -0.13%       |
| \f$3s4d\f$ | \f${}^1D_2\f$   | 53052   | 53135   | -0.16%       |
| \f$3s3p\f$ | \f${}^3P^o_0\f$ | 21787   | 21850   | -0.29%       |
| \f$3s4p\f$ | \f${}^3P^o_0\f$ | 47777   | 47841   | -0.13%       |
| \f$3s6p\f$ | \f${}^3P^o_0\f$ | 54204   | 54249   | -0.08%       |
| \f$3s3p\f$ | \f${}^3P^o_1\f$ | 21809   | 21870   | -0.28%       |
| \f$3s3p\f$ | \f${}^1P^o_1\f$ | 35041   | 35051   | -0.03%       |
| \f$3s4p\f$ | \f${}^3P^o_1\f$ | 47781   | 47844   | -0.13%       |
| \f$3s3p\f$ | \f${}^3P^o_2\f$ | 21852   | 21911   | -0.27%       |
| \f$3s4p\f$ | \f${}^3P^o_2\f$ | 47788   | 47851   | -0.13%       |
| \f$3s6p\f$ | \f${}^3P^o_2\f$ | 54208   | 54253   | -0.08%       |

The discrepancies are now at the level of ~0.1%.

On my PC, this entire calculation took less than two minutes.

### Memory usage

For two-valence CI, the memory bottleneck is calculating the two-body Coulomb (\f$Q^k\f$) integrals; the time bottleneck is usually calculating the two-body Correlation.

You can estimate the memory required using the `ampsci -z` command-line option.
In the above example, we used:

* `20spdf` for CI
* `30spdfghi` for MBPT, with a multipolarity cut-off of 7 (default)

The CI calculation requires only excited-excited Coulomb integrals, while the MBPT only requires Coulomb integrals involving 1 or 2 core states.

We can therefore check the memory requirements by running:
<div class="shell-block">
```bash
./ampsci -z 20spdf [Ne]
./ampsci -z 30spdfghi [Ne] 7
```
</div>

Should return something like:

```text
Estimating memory usage (nb: may be very rough)
For basis         : 20spdf
Total orbitals    : 128
Orbitals in core  : 4
Counting integrals...
With maximum k    : 6
Total integrals   : 17488761
Excited-only      : 15896631
Core-Excited      : 1590740
Estimated Qk size :
All possible Qk   : 1 Gb  [SR]
Excited-only      : 1 Gb  [CI]
Core-Excited      : 0 Gb  (1 or 2 core states allowed)  [MBPT]
```

and

```text
Estimating memory usage (nb: may be very rough)
For basis         : 30spdfghi
Total orbitals    : 348
Orbitals in core  : 4
Counting integrals...
With maximum k    : 7
Total integrals   : 1234351841
Excited-only      : 1204586189
Core-Excited      : 29763502
Estimated Qk size :
All possible Qk   : 45 Gb  [SR]
Excited-only      : 44 Gb  [CI]
Core-Excited      : 1 Gb  (1 or 2 core states allowed)  [MBPT]
```

Combining, we see that we should require 15896631 for the CI part, and 29763502 for MBPT.
This exactly matched the calculated Qk integral summaries we saw above.
In this case, the total memory requirement is very small (~2 Gb).
Notice, however, if we had tried to run CI with the full MBPT basis, we would have required nearly 50Gb!

-----

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

This implies the lifetime of the first \f${}^1P^o_1\f$ level is 2.12 ns, in excellent agreement with experiment.

## Second-order amplitudes: CI_secondOrder

Second-order (in the external field) amplitudes between CI states --
polarisabilities (scalar, tensor, vector/transition) and PNC amplitudes --
are calculated with the `CI_secondOrder` module.
The sums over intermediate states are evaluated with CI mixed states, so
they are complete (no sum over individual CI solutions, and no truncation
of the spectrum); each sum is formed two independent ways as an internal
check.

For two one-body operators \f$ t \f$ (rank \f$ k_t \f$, frequency
\f$ \omega \f$) and \f$ s \f$ (rank \f$ k_s \f$, frequency \f$ \omega_s \f$),
the module reports the reduced amplitude of rank \f$ K \f$,

\f[
  A^K = \sum_n \left[
    c_1 \frac{\langle B||t||n\rangle\langle n||s||A\rangle}{E_A+\omega_s-E_n}
  + c_2 \frac{\langle B||s||n\rangle\langle n||t||A\rangle}{E_A+\omega-E_n}
  \right],
\f]

and the specific cases (\f$ [J]\equiv 2J+1 \f$):

\f[
  \alpha_0 = \frac{A^0}{\sqrt{3[J_A]}},
  \qquad
  \alpha_2 = -\sqrt{\frac{2J(2J-1)}{3(J+1)(2J+1)(2J+3)}}\;A^2,
  \qquad
  \beta = \frac{A^1}{\sqrt{2}\,\langle B||\boldsymbol\sigma||A\rangle},
\f]

for \f$ t=s=E1 \f$ and \f$ K = 0, 2, 1 \f$, and
\f$ E_{\rm PNC} = A^1_0 \f$ for \f$ t=E1 \f$, \f$ s=h_W \f$.
The definitions of \f$ A^{k_tk_s}_{q_1q_2} \f$, \f$ A^K_Q \f$ and
\f$ A^K \f$, and the sign convention of the coupling, are on the
\ref module_secondorder page (identical, with \f$ j \to J \f$) and in
\ref CI::A_K_coefs "CI::A_K_coefs()".

For example, the static scalar polarisability of the ground state:

```java
Module::CI_secondOrder{
  A = 0+;
  K = 0;
}
```

States are given as `J{+,-}:index` (so `0+` is the lowest even J=0 level,
and `1-:1` the second-lowest odd J=1 level). Run `./ampsci -m
CI_secondOrder` for the full list of options; see \ref CI::A_K "CI::A_K()"
for the formulas and conventions. The single-valence analogue is the
\ref module_secondorder module.

-----

* See physics documentation: [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full physics description of the employed methods
