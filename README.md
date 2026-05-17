# ampsci

## Atomic Many-body Perturbation theory in the Screened Coulomb Interaction

[_ampsci_](https://ampsci.dev/)
is a c++ program for high-precision atomic structure calculations of one and two valence atomic systems, developed and maintained by [Benjamin M. Roberts](https://broberts.io/), University of Queensland, Australia

Designed to be fast, accurate, and easy to use.

* Solves the correlated Dirac equation using the Hartree-Fock + correlation potential method (based on Dzuba-Flambaum-Sushkov method) to produce a set of atomic wavefunctions and energies.
* For two-valence systems, uses the CI+MBPT (Configuration Interaction with many-body perturbation theory) method.
* Fully relativistic, includes electron correlations, all-orders screening and hole-particle interaction, finite-nuclear size, Breit interaction, radiative QED effects, RPA for matrix elements, and structure radiation/renormalisation.
* QED is included via the Flambaum-Ginges radiative potential method.
* Can solve for continuum states with high energy, and calculate ionisation cross sections with large energy/momentum transfer.
* Can solve for exotic atoms (e.g., muonic atoms), including electron screening.
* The "modules" system see [ampsci.dev/modules](https://ampsci.dev/modules.html)) makes it relatively simple to add your own routines to use the atomic wavefunctions to calculate whatever properties you may be interested in.
* You can also write wavefunction output to disk in a JSON format, which can be read in to python, for example

The code is on GitHub: [github.com/benroberts999/ampsci](https://github.com/benroberts999/ampsci)

* See [ampsci.dev/](https://ampsci.dev/) for full documentation
  * Includes basic step-by-step introduction and tutorials
  * [ampsci.dev/tutorials](https://ampsci.dev/tutorials.html)
* A full description of the physics methods and approximations, including references,
is given in the physics documentation: [ampsci.pdf][man-url].

**Important:** this is a _pre-release_ version of the code: not fully tested or documented, and should not be used for publishable calculations (without consultation)

[![github][github-badge]](https://github.com/benroberts999/ampsci)
[![doxygen][doxygen-badge]][docs-url]
[![manual][manual-badge]][man-url]

[![tests][tests-badge]][tests-url]
[![build][build-badge]][build-url]
[![macOS][macOS-badge]][macOS-url]
[![cov][cov-badge]][cov-url]

--------------------------------------------------------------------------------

## Compilation and usage

* Full documentation: [ampsci.dev/](https://ampsci.dev/)

### Quick start

```shell
./install-dependencies.sh
./configure.sh -y
make
```

* The `configure.sh` bash script should automatically set up a Makefile, and then `make` will compile ampsci
  * This should work on nearly all systems, but may not work on all
* It assumes dependencies have already been installed. **Requires:**
  * C++ compiler (e.g., g++)
  * lapack, blas, and GSL libraries
  * make (to compile)
  * [optionally] OpenMP
* The `install-dependencies.sh` bash script will check if you have all the required dependencies, and then will ask to install anything that's missing
  * It should work on most systems, but may not work on all
* If there are any issues, see full compilation instructions at [ampsci.dev/compilation.html](https://ampsci.dev/compilation.html)

### Tutorials/examples

* The fastest way to get familiar with ampsci is to follow the tutorials at [ampsci.dev/tutorials.html](https://ampsci.dev/tutorials.html)

Much of ampsci documentation can be seen from the command line:

* See which ampsci input options are available: `./ampsci -i`
  * See available input options for each input block by following with its name
  * e.g., `./ampsci -i HartreeFock`
* Check which Modules are aviable: `./ampsci -m`
  * See available input options for each module by following with its name
  * e.g., `./ampsci -m MatrixElements`
* Check which operators are aviable: `./ampsci -o`
  * See available input options for each operator by following with its name
  * e.g., `./ampsci -o hfs`

Also check out the example input files to get running:

* `doc/examples/ampsci.in` -- an example/template input file
* More: in `doc/examples/`
  * There are several other example input files, along with the expected output; use these to test if everything is working

### Looking for atomic ionisation form-factors for dark-matter-electron scattering?

* See `Kionisation` module (`./ampsci -m Kionisation`)

--------------------------------------------------------------------------------

## Publications

The ampsci code and methods have been described in the following papers

* A. R. Caddell, V. V. Flambaum, B. M. Roberts, _Accurate electron-recoil ionization factors for dark matter direct detection in xenon, krypton and argon_, [Physical Review D **108**, 083030 (2023)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.108.083030), [[arXiv:2305.05125](https://arxiv.org/abs/2305.05125)]
* B. M. Roberts, C. J. Fairhall, J. S. M. Ginges, _Electric dipole transition amplitudes for atoms and ions with one valence electron_, [Physical Review A **107**, 052812 (2023)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.052812), [[arXiv:2211.11134](https://arxiv.org/abs/2211.11134)]

### A selection of publications resulting from ampsci

* R. B. Cserveny, B. M. Roberts, _Theoretical characterization of the barium II and radium II ions_, [Physical Review A 112, 032816 (2025)](https://link.aps.org/doi/10.1103/1rtb-8ymc), [[arXiv:2505.05230](https://arxiv.org/abs/2505.05230)]
* B. Carew, A. R. Caddell, T. N. Maity, C. A. J. O'Hare, _The neutrino fog for dark matter-electron scattering experiments_, [arXiv:2312.04303](https://arxiv.org/abs/2312.04303)
* R. Hamilton _et al._, _Experimental and Theoretical Study of Dynamic Polarizabilities in the_ \f$ 5S - 5D_{5/2} \f$ _Clock Transition in Rubidium-87 and Determination of Electric Dipole Matrix Elements_,  [Physical Review Applied **19**, 054059 (2023)](https://link.aps.org/doi/10.1103/PhysRevApplied.19.054059), [[arXiv:2212.10743](https://arxiv.org/abs/2212.10743)]
* C. J. Fairhall, B. M. Roberts, J. S. M. Ginges, _QED radiative corrections to electric dipole amplitudes in heavy atoms_, [Physical Review A **107**, 022813 (2023)](https://link.aps.org/doi/10.1103/PhysRevA.107.022813), [[arXiv:2212.11490](https://arxiv.org/abs/2212.11490)]
* The COSINE-100 Collaboration, _Constraints on sub-GeV dark matter scattering on electrons with COSINE-100_, [Physical Review D 113, 072012 (2026)](https://journals.aps.org/prd/abstract/10.1103/kszv-g7tl).
* The NEON Collaboration, _First Direct Search for Light Dark Matter Using the NEON Experiment at a Nuclear Reactor_, [Physical Review Letters 134, 021802 (2025)](https://link.aps.org/doi/10.1103/PhysRevLett.134.021802)
* J. M. Cline, M. Puel, and T. Toma, _Boosted dark matter from a phantom fluid_, [Physics Letters B 848, 138377 (2024)](https://www.sciencedirect.com/science/article/pii/S0370269323007116)
* The XENON Collaboration, _Search for Electronic Recoil Event Rate Modulation with 4 Years of XENON100 Data_, [Physical Review Letters 118, 101101 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.118.101101)

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
[cov-badge]: https://img.shields.io/codecov/c/github/benroberts999/ampsci?token=3M5MH5QXLL
[cov-url]: https://codecov.io/gh/benroberts999/ampsci
[c++-badge]: https://img.shields.io/badge/c++-17-blue
[github-badge]: https://img.shields.io/badge/Code%20available:-GitHub-blueviolet?style=flat&logo=github&logoColor=white

[tests-badge-v2]: tests-badge.svg
[build-badge-v2]: build-badge.svg
[macOS-badge-v2]: macOS-badge.svg
[cov-badge-v2]: cov-badge.svg
