# ampsci

## Atomic Many-body Perturbation theory in the Screened Coulomb Interaction

[_ampsci_](https://ampsci.dev/)
is a c++ program for high-precision atomic structure calculations of single-valence systems, 
developed and maintained by [Benjamin M. Roberts](https://broberts.io/), University of Queensland, Australia

It solves the correlated Dirac equation using the Hartree-Fock + correlation potential method (based on Dzuba-Flambaum-Sushkov method) to produce a set of atomic wavefunctions and energies.
The method is fully relativistic, includes electron correlations, all-orders screening and hole-particle interaction, finite-nuclear size, Breit interaction, radiative QED effects, RPA for matrix elements, and structure radiation/renormalisation.
QED is included via the Flambaum-Ginges radiative potential method.
Can solve for continuum states with high energy, and calculate ionisation cross sections with large energy/momentum transfer.

Designed to be fast, accurate, and easy to use.
The "modules" system (see [doc/modules.md](doc/modules.md)) makes it simple to add your own routines to use the atomic wavefunctions to calculate whatever properties you may be interested in.

* The code is on GitHub: [github.com/benroberts999/ampsci](https://github.com/benroberts999/ampsci)
* See [ampsci.dev/](https://ampsci.dev/) for full documentation
* A full description of the physics methods and approximations, including references,
is given in the physics documentation: [ampsci.pdf][man-url].
* **Important:** this is a _pre-release_ version of the code: not fully tested or documented, and should not be used for publishable calculations (without consultation)

[![github][github-badge]](https://github.com/benroberts999/ampsci)
[![doxygen][doxygen-badge]][docs-url]
[![manual][manual-badge]][man-url]

[![tests][tests-badge]][tests-url]
[![build][build-badge]][build-url]
[![macOS][macOS-badge]][macOS-url]
[![cov][cov-badge]][cov-url]

--------------------------------------------------------------------------------

## Compilation and usage

Full documentation available online: [ampsci.dev/](https://ampsci.dev/).
Divided into sections:

 1. Usage intructions, input options
    * Compilation instructions (for linux/mac/windows): [doc/compilation.md](doc/compilation.md)
    * Detailed info on all input options: [doc/ampsci_input.md](doc/ampsci_input.md)
    * A basic step-by-step tutorial: [doc/tutorial.md](doc/tutorial.md)
    * A more advanced tutorial: [doc/tutorial_advanced.md](doc/tutorial_advanced.md)
      * See also: _doc/examples/ampsci.in_ -- an example/template input file
      * In _doc/examples/_ there are several example input files, with the expected output; use these to test if everything is working

 2. Physics documentation: [ampsci.dev/ampsci.pdf](https://ampsci.dev/ampsci.pdf)
    * Description of physics/methods used in the code
    * Includes many references to the works where the methods implemented here were developed.

 3. Modules
    * The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated. See [doc/modules.md](doc/modules.md) for description
    * The code is designed so that you can easily create your own modules. See [doc/writing_modules.md](doc/writing_modules.md) for details

 4. Code documentation -- details on classes/functions in the code
    * Available online: [ampsci.dev/](https://ampsci.dev/)
    * This should only be required if you plan to edit the code or add new modules

--------------------------------------------------------------------------------

[tests-badge]: https://github.com/benroberts999/ampsci/actions/workflows/tests.yml/badge.svg
[tests-url]: https://github.com/benroberts999/ampsci/actions/workflows/tests.yml
[build-badge]: https://github.com/benroberts999/ampsci/actions/workflows/build.yml/badge.svg
[build-url]: https://github.com/benroberts999/ampsci/actions/workflows/build.yml
[macOS-badge]: https://github.com/benroberts999/ampsci/actions/workflows/macOS.yml/badge.svg
[macOS-url]: https://github.com/benroberts999/ampsci/actions/workflows/macOS.yml
[doxygen-badge]: https://img.shields.io/badge/documentation-ampsci.dev/-blue
[docs-url]: https://ampsci.dev/
[manual-badge]: https://img.shields.io/badge/documentation-physics%20(pdf)-blue
[man-url]: https://ampsci.dev/ampsci.pdf
[cov-badge]: https://codecov.io/gh/benroberts999/ampsci/branch/main/graph/badge.svg?token=3M5MH5QXLL
[cov-url]: https://codecov.io/gh/benroberts999/ampsci
[c++-badge]: https://img.shields.io/badge/c++-17-blue
[github-badge]: https://img.shields.io/badge/Code%20available:-GitHub-blueviolet?style=flat&logo=github&logoColor=white

[tests-badge-v2]: tests-badge.svg
[build-badge-v2]: build-badge.svg
[macOS-badge-v2]: macOS-badge.svg
[cov-badge-v2]: cov-badge.svg
