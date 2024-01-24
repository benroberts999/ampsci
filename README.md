# ampsci

## Atomic Many-body Perturbation theory in the Screened Coulomb Interaction

[_ampsci_](https://ampsci.dev/)
is a c++ program for high-precision atomic structure calculations of one and two valence atomic systems, developed and maintained by [Benjamin M. Roberts](https://broberts.io/), University of Queensland, Australia

It solves the correlated Dirac equation using the Hartree-Fock + correlation potential method (based on Dzuba-Flambaum-Sushkov method) to produce a set of atomic wavefunctions and energies.
For two-valence systems, uses the CI+MBPT (Configuration Interaction with many-body perturbation theory) method.
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

### Quick start

* The `setup.sh` bash script should compile and build ampsci (uses Make)
  * It uses only defaults, and may not work on all systems.
    See documentation for full guide to compilation.
  * It assumes dependencies have already been installed. If not, see next:
* The `install-dependencies.sh` bash script should install all required dependencies.
  * It uses only defaults, and may not work on all systems.
* If there are issue with compilation
  * See full compilation instructions: [doc/compilation.md](doc/compilation.md)
* Check out the example input files to get running:
  * [_doc/examples/ampsci.in_](doc/examples/ampsci.in) -- an example/template input file
  * More: in [_doc/examples/_](doc/examples/) there are several example input files, along with the expected output; use these to test if everything is working

### Tutorials/examples

* The fastest way to get familiar with ampsci is to follow the tutorials
* A basic step-by-step tutorial: [doc/tutorial.md](doc/tutorial.md)
* More advanced tutorials follow: [doc/tutorial_advanced.md](doc/tutorial_advanced.md), [doc/tutorial_CI.md](doc/tutorial_CI.md)
* See which ampsci input options are available: `./ampsci -a`
  * See available input options for each input block by following with its name
  * e.g., `./ampsci -a HartreeFock`
* Check which Modules are aviable: `./ampsci -m`
  * See available input options for each module by following with its name
  * e.g., `./ampsci -m MatrixElements`
* Check which operators are aviable: `./ampsci -o`
  * See available input options for each operator by following with its name
  * e.g., `./ampsci -o hfs`

### Looking for atomic ionisation form-factors for dark-matter-electron scattering?

* See `Kionisation` module (`./ampsci -m Kionisation`)

### Documentation

Full documentation available online: [ampsci.dev/](https://ampsci.dev/).
Divided into sections:

 1. Usage intructions, input options
    * Compilation instructions (for linux/mac/windows): [doc/compilation.md](doc/compilation.md)
    * Detailed info on all input options: [doc/ampsci_input.md](doc/ampsci_input.md)

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

#### A selection of publications resulting from ampsci

* _Accurate electron-recoil ionization factors for dark matter direct detection in xenon, krypton and argon_, A. R. Caddell, V. V. Flambaum, B. M. Roberts, [Physical Review D **108**, 083030 (2023)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.083030), [[arXiv:2305.05125](https://arxiv.org/abs/2305.05125)]
* _The neutrino fog for dark matter-electron scattering experiments_, B. Carew, A. R. Caddell, T. N. Maity, C. A. J. O'Hare, [arXiv:2312.04303](https://arxiv.org/abs/2312.04303)
* _Electric dipole transition amplitudes for atoms and ions with one valence electron_, B. M. Roberts, C. J. Fairhall, J. S. M. Ginges, [Physical Review A **107**, 052812 (2023)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.052812), [[arXiv:2211.11134](https://arxiv.org/abs/2211.11134)]
* _Experimental and Theoretical Study of Dynamic Polarizabilities in the_ $5S - 5D_{5/2}$ _Clock Transition in Rubidium-87 and Determination of Electric Dipole Matrix Elements_, R. Hamilton _et al._, [Physical Review Applied **19**, 054059 (2023)](https://link.aps.org/doi/10.1103/PhysRevApplied.19.054059), [[arXiv:2212.10743](https://arxiv.org/abs/2212.10743)]
* _QED radiative corrections to electric dipole amplitudes in heavy atoms_, C. J. Fairhall, B. M. Roberts, J. S. M. Ginges, [Physical Review A **107**, 022813 (2023)](https://link.aps.org/doi/10.1103/PhysRevA.107.022813), [[arXiv:2212.11490](https://arxiv.org/abs/2212.11490)]

--------------------------------------------------------------------------------

[tests-badge]: https://github.com/benroberts999/ampsci/actions/workflows/tests.yml/badge.svg
[tests-url]: https://github.com/benroberts999/ampsci/actions/workflows/tests.yml
[build-badge]: https://github.com/benroberts999/ampsci/actions/workflows/build.yml/badge.svg
[build-url]: https://github.com/benroberts999/ampsci/actions/workflows/build.yml
[macOS-badge]: https://github.com/benroberts999/ampsci/actions/workflows/macOS.yml/badge.svg
[macOS-url]: https://github.com/benroberts999/ampsci/actions/workflows/macOS.yml
[doxygen-badge]: https://img.shields.io/badge/code%20docs-ampsci.dev/-blue
[docs-url]: https://ampsci.dev/
[manual-badge]: https://img.shields.io/badge/physics%20docs-ampsci.pdf-blue
[man-url]: https://ampsci.dev/ampsci.pdf
[cov-badge]: https://codecov.io/gh/benroberts999/ampsci/branch/main/graph/badge.svg?token=3M5MH5QXLL
[cov-url]: https://codecov.io/gh/benroberts999/ampsci
[c++-badge]: https://img.shields.io/badge/c++-17-blue
[github-badge]: https://img.shields.io/badge/Code%20available:-GitHub-blueviolet?style=flat&logo=github&logoColor=white

[tests-badge-v2]: tests-badge.svg
[build-badge-v2]: build-badge.svg
[macOS-badge-v2]: macOS-badge.svg
[cov-badge-v2]: cov-badge.svg
