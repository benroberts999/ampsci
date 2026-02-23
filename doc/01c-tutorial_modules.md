\page tutorial_modules Matrix Elements

\brief Calculate matrix elements of various operators

Here, we'll see how to use an ampsci `Module` to calculate matrix elements of various opertors.

This assumes you already have ampsci compiled, and are familiar with running basic calculations.

* See [Getting Started](\ref getting_started) for basic compilation instructions
* See [basic Tutorial](\ref tutorial_basic) for running basic calculations
* See [MBPT and Correlations](\ref tutorial_mbpt) for running more advanced calculations.

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
* And [Writing Custom Modules](\ref modules_custom) for information on writing your own custom module

## Calculating matrix elements

When we ran `./ampsci -m`, we would have seen a large list of available modules.
One of them would have been called `matrixElements`.
We tell amplsci to run this module by adding it as an input block to the input file:

```java
Module::matrixElements{}
```

Most modules will take input options.
To see the available options for this module, list the block name after `-m` on the command-line:

<div class="shell-block">
```bash
./ampsci -m matrixElements
```
</div>

(or `./ampsci -i Module::matrixElements`), which prints:

```java
// Available Module::matrixElements options/blocks
Module::matrixElements{
  operator;
    // e.g., E1, hfs (see ampsci -o for available operators)
  options{}
    // options specific to operator (see ampsci -o
    // 'operator')
  rpa;
    // Method used for RPA: true(=TDHF), false, TDHF, basis,
    // diagram [true]
  rpa_options{}
    // Block: some further options for RPA
  omega;
    // Text or number. Freq. for RPA (and freq. dependent
    // operators). Put 'each' to solve at correct frequency
    // for each transition. [0.0]
  printBoth;
    // print <a|h|b> and <b|h|a> [false]
  include_core;
    // If true, includes core states in calculation. Will
    // use HF core, unless use_spectrum is true [false]
  use_spectrum;
    // If true (and spectrum available), will use spectrum
    // for valence states [false]
  diagonal;
    // Calculate diagonal matrix elements (if non-zero)
    // [true]
  off-diagonal;
    // Calculate off-diagonal matrix elements (if non-zero)
    // [true]
  moments;
    // true/false - calculate moments, or reduced matrix
    // elements. Default is false, except for hyperfine
    // operators
  StructureRadiation{}
    // Options for Structure Radiation and normalisation
    // (details below)
}

// Available StructureRadiation options/blocks
StructureRadiation{
  // If this block is included, SR + Normalisation
  // corrections will be included
  Qk_file;
    // true/false/filename - SR: filename for QkTable file.
    // If blank will not use QkTable; if exists, will read
    // it in; if doesn't exist, will create it and write to
    // disk. If 'true' will use default filename. Save time
    // (10x) at cost of memory. Note: Using QkTable implies
    // splines used for diagram legs
  n_minmax;
    // list; min,max n for core/excited: [1,inf]
}

// Available rpa_options options/blocks
rpa_options{
  eps;
    // Convergance goal [1.0e-10]
  eta;
    // Damping factor - be carful with this [0.4]
  max_iterations;
    // Maximum number of iterations. 1 should correspond to
    // first-order RPA [100]
}
```

The first option, `operator` is which operator we want to calculate matrix elements of.
You can see a list of operators with the `-o` command-line option:

<div class="shell-block">
```bash
./ampsci -o
```
</div>

The second option, which is a sub-input-block, `options` is the set of options for the given operator. Most operators have no extra options, so this can be left blank. Some (e.g., hyperfine operator `hfs` have many available options). Like above, you can see the available options by further passing the operator name after `-o`. For example:

<div class="shell-block">
```bash
./ampsci -o hfs
```
</div>

Here we will consider the simpler `E1` operator.
To our above `example.in` file, we can add the following block (note we may add as many Module:: blocks as we like, they will all be run one-by-one in order):

```java
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

-----

## MBPT for matrix elements: Core Polarisation/RPA <a name="rpa"></a>

Core polarisation (RPA) is included in the matrix elements.

The best method to use is TDHF, which is numerically stable, and includes contribution from negative energy states automatically.

-----

## MBPT for matrix elements: Structure Radiation <a name="sructrad"></a>

To improve tha accuracy of matrix elements, structure radiation and normalisation corrections should be included.
There is an option to do this in the MatrixElements module.
There is also a `Module::StructureRad` module, which gives some finer control.

```java
// Available Module::matrixElements options/blocks
Module::matrixElements{
  //.... same as before

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

If a Qk filename is given, program will first calculate all required Q^k Coulomb integrals before calculating structure radiation. This speeds up the calculation, at a great memory cost.
