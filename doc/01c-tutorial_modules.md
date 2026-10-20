\page tutorial_modules Matrix Elements

\brief Calculate matrix elements of various operators

Here, we'll see how to use an ampsci `Module` to calculate matrix elements of various opertors.

This assumes you already have ampsci compiled, and are familiar with running basic calculations.

* See [Getting Started](\ref getting_started) for basic compilation instructions
* See [basic Tutorial](\ref tutorial_basic) for running basic calculations
* See [MBPT and Correlations](\ref tutorial_mbpt) for running more advanced calculations.
* See [Modules](\ref modules) and [Operators](\ref using_operators) for detailed documentation of the most common modules and operators.

## Modules: using the wavefunctions <a name="modules"></a>

Above, we ran ampsci, which calculated the atomic wavefunctions and printed their energies to screen.
If we want to actually _do_ anything with the wavefunctions, we have to run one or more **modules**.
We do this by adding a module block to the input file, which has the form `Module::ModuleName{}`.
We can see a list of all available modules with the `-m` command-line option:

<div class="shell-block">
```bash
./ampsci -m
```
</div>

The general usage of the code is to first use the main blocks to construct the
atomic wavefunction and basis states, then to add as many `Module::` blocks as
required.
Each module is a separate routine that will take the calculated wavefunction and compute any desired property (e.g., matrix elements).
They are independent and do not talk to each other, though may write output to the disk.
There are several available modules, here we will just focus on calculating matrix elements.
The code is designed so that anyone can write a new Module to calculate anything else they may desire.

* See [Modules](\ref modules) for more detail on general modules
  * Use `ampsci -m` to see a list of available modules
* And [Writing Custom Modules and Operators](\ref modules_custom) for information on writing your own modules and operators

## Calculating matrix elements

When we ran `./ampsci -m`, we would have seen a large list of available modules.
One of them would have been called `matrixElements`.
We tell ampsci to run this module by adding it as an input block to the input file.
Let's walk through a complete example for Cs: compute the wavefunctions at
Hartree-Fock level, then the E1 matrix elements (with core polarisation, RPA)
between the valence states:

```java
Atom { Z = Cs; }
HartreeFock { core = [Xe]; valence = 6sp; }
Grid { r0 = 1e-6; rmax = 120.0; num_points = 3000; }

Module::matrixElements{
  operator = E1;
  rpa = true;
}
```

Run this in the usual way (`./ampsci example.in`). After the Hartree-Fock
output, the module runs and prints:

```text
Matrix Elements - Operator: E1
Reduced matrix elements
Units: |e|aB
Including RPA: TDHF method
RPA solved at omega = 0
TDHF E1 (w=0.0000): 19 9.7e-11 [3p+,d-]

E1

   a    b    w_ab        t0_ab          +RPA
 6p-  6s+    0.0417522  -5.277687e+00  -4.974398e+00
 6p+  6s+    0.0435826   7.426435e+00   7.013071e+00
```

The first lines show what is being calculated and how well the RPA (TDHF)
equations converged. Then, for each pair of valence states allowed by the
selection rules: the transition frequency `w_ab` (in au), the lowest-order
(Hartree-Fock) reduced matrix element `t0_ab`, and the value including the
RPA correction.

Most modules take input options; list them by giving the module name after
`-m` on the command-line (the output is in input-file format, so you can
copy+paste it):

<div class="shell-block">
```bash
./ampsci -m matrixElements
```
</div>

See [Matrix Elements](\ref module_matrixelements) for the full set of
options (RPA methods, frequencies, core states, etc.).

## Choosing the operator

The `operator` option names the operator; see a list with the `-o`
command-line option:

<div class="shell-block">
```bash
./ampsci -o
```
</div>

The `options{}` sub-block holds operator-specific options. Most operators
have none (like E1 above); some, e.g. the hyperfine operator `hfs`, have
many -- list them by passing the operator name: `./ampsci -o hfs`.
For example, hyperfine A constants with a "ball" nuclear magnetisation
distribution, for the diagonal matrix elements only:

```java
Module::matrixElements{
  operator = hfs;
  options { F = ball; }
  off-diagonal = false;
}
```

which outputs (note the units: for hyperfine, constants in MHz are printed
rather than reduced matrix elements):

```text
Hyperfine structure: Cs, Z=55 A=133
K=1 (magnetic dipole)
Using ball nuclear distro for F(r)
w/ r_N = 6.20207fm = 0.000117202au  (r_rms=4.8041fm)
mu = 2.5778, I = 3.5, g = 0.736514

Matrix Elements - Operator: hfs1
Hyperfine A constants (M1)
Units: MHz
Including RPA: TDHF method

   a    b    w_ab        t0_ab          +RPA
 6s+  6s+    0.0000000   1.421144e+03   1.712908e+03
 6p-  6p-    0.0000000   1.606247e+02   2.011687e+02
 6p+  6p+    0.0000000   2.387730e+01   4.268347e+01
```

Note we may add as many `Module::` blocks as we like; they are run
one-by-one in order. See [Operators](\ref using_operators) for more detail
and examples (E1, hyperfine, PNC).

-----

## MBPT for matrix elements: Core Polarisation/RPA <a name="rpa"></a>

Core polarisation (RPA) is included in the matrix elements with the `rpa`
option, as above.
The best method to use is TDHF (the default), which is numerically stable,
and includes contributions from negative energy states automatically; the
`diagram` and `basis` methods are also available (these require a basis).
The `omega` option sets the frequency the RPA is solved at
(`omega = each;` re-solves at each transition frequency; the frequency
dependence is usually small).
See [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for a description of RPA.

-----

## MBPT for matrix elements: Structure Radiation <a name="sructrad"></a>

To improve the accuracy of matrix elements, structure radiation and
normalisation corrections should be included.
These are calculated with the separate
[Structure Radiation](\ref module_structurerad) module
(`Module::StructureRadiation{}`), which also gives the Brueckner orbital
corrections.

If a `Qk_file` filename is given, the program will first calculate all required Q^k Coulomb integrals before calculating structure radiation. This speeds up the calculation, at a great memory cost.

Note: this can take a significant amount of memory.
You can estimate the memory required, without running it, using the `ampsci -z` command-line option with the input file:

<div class="shell-block">
```bash
./ampsci -z example.in
```
</div>

With a structure radiation module in the input, this includes the full table of \f$Q^k\f$ integrals over the basis (`Q^k (all; SR)`).

-----

## Where to next

* [Matrix Elements](\ref module_matrixelements) - full options for the matrixElements module
* [Structure Radiation](\ref module_structurerad) - SR + normalisation corrections
* [Second-Order Amplitudes](\ref module_secondorder) - polarisabilities, PNC amplitudes
* [Dynamic Polarisability](\ref module_dynamicpol) - alpha(omega) over a frequency range
* [Operators](\ref using_operators) - available operators and their options
