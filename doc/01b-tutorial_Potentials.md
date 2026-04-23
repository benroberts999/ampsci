\page tutorial_potentials Atomic Potentials

\brief Atomic Potentials: Hartree-Fock, Nuclear, Breit and QED

This assumes you already have ampsci compiled and have a basic understanding of how to run and use it.

* See [Compilation](\ref compilation) for compilation instructions
* See [Basic Tutorial](\ref tutorial_basic) for getting started if unfamiliar
* See [ampsci.dev/ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full description of the physics

## Contents

* [Starting approximation:](#begin)
* [Hartree-Fock](#hf)
* [Nuclear Potential](#nuclear)
* [Breit](#breit)
* [Radiative QED](#qed)

## Starting approximation <a name="begin"></a>

The full atomic Hamiltonian for an \f$N\f$-electron atom:
\f[
  H = \sum_i^N h_0(r_i) + \sum_{i < j}\frac{1}{r_i - r_j}
\f]

\f[
  h_0(\vec{r_i}) = c \vec\alpha_i\cdot\vec{p_i} + c^2 (\beta_i-1) + V_{\rm nuc}.
\f]

The starting approximation is the Hartree-Fock method:

\f[
  H \approx \sum_i^N [h_0(\vec r_i) + v^{\rm HF} ].
\f]

## Hartree-Fock <a name="hf"></a>

We focus on case of single-valence systems, and start with so-called \f$V^{N-1}\f$ approximation, in which Hartree-Fock potential is due to the \f$N-1\f$ core electrons:

\f[
  \hat v^{\rm HF}\phi_a(\vec{r_1}) = \sum_{i\neq a}^{N_c}\Bigg(
  \int \frac{\phi_i^\dagger(\vec{r_2})\phi_i(\vec{r_2})}{|r_{12}|}d^3\vec{r_2}\,\phi_a(\vec{r_1})
  -\int \frac{\phi_i^\dagger(\vec{r_2})\phi_a(\vec{r_2})}{|r_{12}|}d^3\vec{r_2}\,\phi_i(\vec{r_1})
  \Bigg),
\f]

First, Hartree-Fock equations are solved self-consistently for all core electrons, then the valence states are found in the Frozen Hartree-Fock potential due to the core.

-----

As always, begin by checking the available options:

```bash
 ./ampsci -i HartreeFock
```

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
  // Scale for factor for Breit Hamiltonian. Usially 0.0 (no Breit) or 1.0 (full
  // Breit), but can take any value. [0.0]
  Breit;
  // Include QED? Three options: true, false, valence. If 'valencel, will
  // include QED only into valence states, but not the core. Detailed QED
  // options are set within the RadPot{} block - if that block is not set,
  // defaults will be used. By default, this option is false, unless the
  // RadPot{} block exists, in which case it is true
  QED;
}
```

For example, to run Hartree-Fock for Cs with a \f$V^{N-1}\f$ potential (i.e. Xe-like core),

```java
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
}
```

If you're unsure of core/valence states, use ampsci's in-built periodic table

```ampsci -p```
or
```ampsci -p Cs```
which will give:

```text
Cs,  cesium.
Z = 55;  A = 133 (default)

Electron config: [Xe],6s1   (guess)
 = 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 | 6s1 

Isotpe data:
Cs-133 (Z=55, A=133)
r_rms = 4.8041, c = 5.67073, mu = 2.5778, I = 3.5, parity = 1
```

* **Note**: The core configuration is a _guess_ based on basic electron filling rules. It's usually a good configuration for the ground-state, but may not be exact

-----

Besides Hartree-Fock, other methods are available:

* Hartree: exchange term excluded
* KohnSham: Uses Kohn-Sham density functional, which includes a localised exchange potential and Latter correction (nb: lowest valence electron should usually be included into the 'core' here, i.e. `core=[Xe]:6s1')
  * The colon tells ampsci to include 6s1 into the core potential, but above the 'Fermi leve; (i.e., also allow it in the valence).
* ApproxHF: Approximate (localised) Hartree-Fock potential (only for tests)
* Local: Uses a local parametric potential (Green's potential)
  * Parameters are tuned to roughly match Hartree-Fock
  * Only should be used as a numerical check/test

e.g.,

```java
HartreeFock{
  core = [Xe]:6s1;
  valence = 6sp5d;
  method = KohnSham;
}
```

-----

## Nuclear Potential <a name="nuclear"></a>

For a point-like nucleus, \f$ V_{\rm nuc} = -Z/r \f$ ;
in reality, the nuclear charge is distributed across the finite-size nucleus.
By default, to form \f$ V_{\rm nuc}, \f$ we assume the nuclear charge follows a Fermi distribution,
\f[
    \rho(r)
        = \frac{\rho_0}{1+\exp[(r-c)/a]},
\f]
where \f$\rho_0\f$ is a normalisation factor (\f$\int\rho\,{\rm d} V = Z\f$),
_c_ is the half-density radius, and
_a_ is defined via the 90 - 10% density fall-off \f${t\equiv 4a\ln3}\f$ (known as the ``skin thickness''),
which we take to be \f$t=2.3\,{\rm fm}\f$ by default for all heavy isotopes.

If no inputs are given, the code will assume a Fermi-like distribution for the nuclear charge, and look up default values for the nuclear parameters (charge radius and skin thickness). The default values are chosen from the specified isotope in the `Atom{}` block. These may also be specified specifically:

```java
// Available Nucleus options/blocks
Nucleus{
  // Options for nuclear potential (finite nuclear size). All are optional.
  // Default is a Fermi-like nucleus, with parameters chosen according to
  // isotope (see Atom{A;})

  // Root-mean-square charge radius, in fm. If not given, will use default for
  // given Z and A.
  rrms;
  // Half-density radius, in fm. ( If c and rrms are defined, c will over-ride
  // rms). Default depends on Z and A.
  c;
  // Fermi, spherical, point-like, Gaussian [Fermi], custom
  type;
  // Nuclear skin thickness, in fm [2.3]
  t;
  // Intrinsic nuclear quadrupole deformation parameter beta (Currently only for
  // Fermi potential) [0.0]
  beta;
  // Read potential from given file of the form (r v(r)) (space delimetered, in
  // atomic units) - will be interpolated onto atomic grid. Note: potential will
  // be extrapolated to larger R, assuming ~1/R, but will NOT be extrapolated to
  // smaller R. So, ensure smallest required {r V(r)} is given
  input_file;
  // List of comma separated real numbers. Not currently unused, but may be used
  // in future for more complicated nuclear distros.
  parameters;
}
```

* rrms: nuclear root-mean-square charge radius, in fm
* type: Nuclear charge distribution. The available nuclear types are: `Fermi`, `spherical`, `point-like`, and `Gaussian`, with `Fermi` being the default.
  * c: nuclear half-density radius (only meaningful for Fermi distribution). Usually this is not set; if it is set, it will _override_ rrms: \f$3c^2=5r_{\rm rms}^2-7\pi^2a^2\f$, and \f$t=4a~{\rm ln}3\f$
  * t: nuclear ``skin thickness'' (90-10% fall-off radius), in fm. Only for Fermi distribution.
  * \f$3c^2=5r_{\rm rms}^2-7\pi^2a^2\f$, and \f$t=4a~{\rm ln}3\f$
  * beta is deformation paramater; only used with Fermi model, to approximately model nuclear deformation (see [Flambaum, Dzuba, Phys. Rev. A **100**, 032511 (2019)](https://doi.org/10.1103/PhysRevA.100.032511))
* **note** Code almost exclusively uses atomic units, however, as the convention is so strong, the nuclear parameters are all set in fm.
* `input_file` allows a numerically determined nuclear _potential_ to be read of from a file. **note:**
  * The file must contain the _potential_, not the density
  * list, space separated, of r and V(r) values
  * Both r and V(r) must be in atomic units
  * V(r) should go to \f$-Z/r\f$ at lare r
  * They will be interpolated onto the grid used in ampsci, so the grid spacings are not important (so long as they are dense enough to capture any required detail)
  * V(r) will be extrapolated to larger r assuming \f$-Z/r\f$ - but will be not be extrapolated to smaller r than given in the file (set to zero)

Default rms values are taken from:

* Angeli, Marinova, [At. Data Nucl. Data Tables **99**, 69 (2013)](http://dx.doi.org/10.1016/j.adt.2011.12.006).
* It's typically best to leave these as the default parameters, but you're free to update them, e.g., the following will assume a spherically-symmetric (uniformly charged ball) nucleus with rms charge radius of 3.5 fm:

```java
Nucleus{
  rrms = 3.5;
  type = spherical;
}
```

* **Note**: For certain isotopes, the rrms radius may not be in the Angeli tables. In those cases, a default approximate formula is used, and a warning will be printed. The approximate formula is fine in most cases, but it should be checked, particularly in finite-nuclear-size effects are important for the given calculation

-----

## Breit <a name="breit"></a>

The Breit Hamiltonian accounts for magnetic interactions between electrons (also known as the Gaunt interaction), and retardation effects.
It leads to a correction to the electron-electron Coulomb term in the many-body Hamiltonian:

\f[
  \sum_{ij}\frac{1}{r_{ij}}
  \to
  \sum_{ij}\left( \frac{1}{r_{ij}} + \hat h^B_{ij}\right),
\f]

where, in the limit of zero frequency, the two-particle Breit Hamiltonian is

\f[
  h^B_{ij} = - \frac{\vec{\alpha_i}\cdot\vec{\alpha_j} + (\vec{\alpha_i}\cdot\hat{n_{ij}})(\vec{\alpha_j}\cdot\hat{n_{ij}})}{2\, r_{ij}}.
\f]

This can be included at the Hartree-Fock level with the `HartreeFock{Breit = true;}` setting.

```java
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
  Breit = 1.0;
}
```

This setting is a scaling factor; i.e., setting `Breit = 0.5;` will add effective factor of 0.5 in front of \f$h^B\f$. Typically, only 0 or 1 is set; other values are useful for checking for non-physical non-linear-in-Breit effects. Setting to `true` is equivilant to setting to 1.0.

-----

## Radiative QED <a name="qed"></a>

Radiative QED corrections can be included into the wavefunctions using the Flambaum-Ginges radiative potential method.
An effective potential, \f$V_{\rm rad}\f$,
is added to the Hamiltonian before the equations are solved.
The potential can be written as the sum of the Uehling (vacuum polarisation), Wichmann-Kroll (higher-order vacuum polarisation) and self-energy potentials; the self-energy potential itself is written as the sum of the high- and low-frequency electric contributions, and the magnetic contribution:

\f[
  V_{\rm rad}(\vec{r}) = V_{\rm Ueh}(r) + V_{\rm WK}(r) + V_{\rm SE}^{h}(r) +  V_{\rm SE}^{l}(r) + i (\vec{\gamma}\cdot\hat{n}) V^{\rm mag}(r).
\f]

* Flambaum, Ginges, [Phys. Rev. A **72**, 052115 (2005)](http://link.aps.org/doi/10.1103/PhysRevA.72.052115).

This can be included by setting `QED=true;` in the Hartree Fock block:

```java
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
  QED = true;
}
```

* This will use all default values (see below)
* You can also include QED into the valence states only (not the core), by setting `QED=valence;`
  * This is only really for tests, to check to contribution of QED in the core; if you do this, the valence will not be orthogonal to the core.

You can all set detailed QED options within the `RadPot{}` Block. If `QED=true;` in HartreeFock, and the `RadPot{}` block is not explicitely included, then the QED radiative potential will be included with the default parameters. This is equivilant to

```java
RadPot{
  Ueh = 1.0;
  SE_h = 1.0;
  SE_l = 1.0;
  SE_m = 1.0;
  WK = 0.0;
  scale_l;
  scale_rN;
}
```

* Each of the (`Ueh, SE_h, SE_l, SE_m, WK`) options are scaling factors for their corresponding terms in the potential; typically 0 or 1, but can be tuned (for testing).
* `scale_l` takes a list of scaling factors for each l; these will rescale \f$V_{\rm rad}\f$ for each partial wave.
e.g., `scale_l = 0,1,0;` will include \f$V_{\rm rad}\f$ for p states, but not s or d states.
* `scale_rN;` Re-scales the effective nuclear radius; used to test finite-nuclear size effects on \f$V_{\rm rad}\f$. =1 means normal, =0 means assume point-like nucleus (when calculating \f$V_{\rm rad}\f$).
* Usually, the defaults should be used; this block can be used for testing if required.

-----

* See [Modules](\ref tutorial_modules) for tutorial on using the wavefunctions (e.g., calculating matrix elements)
* See [MBPT](\ref tutorial_mbpt) for tutorial on doing more advanced calculations, including correlation corrections
* See physics documentation: [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full physics description of the employed methods
