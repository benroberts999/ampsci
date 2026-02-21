\page tutorial_mbpt Advanced Tutorial: MBPT

\brief MBPT and correlation corrections

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
* [Constructing a basis](#basis)
  * [Spectrum](#spectrum)
* [MBPT: second-order correlations](#secondorder)
* [MBPT: all-order correlations](#allorders)
* [MBPT for matrix elements: Core Polarisation/RPA](#rpa)
* [MBPT for matrix elements: Structure Radiation](#sructrad)

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

### Hartree-Fock <a name="hf"></a>

We focus on case of single-valence systems, and start with so-called \f$V^{N-1}\f$ approximation, in which Hartree-Fock potential is due to the \f$N-1\f$ core electrons:

\f[
  \hat v^{\rm HF}\phi_a(\vec{r_1}) = \sum_{i\neq a}^{N_c}\Bigg(
  \int \frac{\phi_i^\dagger(\vec{r_2})\phi_i(\vec{r_2})}{|r_{12}|}d^3\vec{r_2}\,\phi_a(\vec{r_1})
  -\int \frac{\phi_i^\dagger(\vec{r_2})\phi_a(\vec{r_2})}{|r_{12}|}d^3\vec{r_2}\,\phi_i(\vec{r_1})
  \Bigg),
\f]

First, Hartree-Fock equations are solved self-consistently for all core electrons, then the valence states are found in the Frozen Hartree-Fock potential due to the core.

------

As always, begin by checking the available options:

```bash
 ./ampsci -a HartreeFock
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

If you're unsire of core/valence states, use ampsci's ib-built periodic table

```ampsci -p```
or
```ampsci -p Cs```
which will give:

```text
  H                                                                  He 
  1                                                                   2 
 Li  Be   B                                           C   N   O   F  Ne 
  3   4   5                                           6   7   8   9  10 
 Na  Mg  Al                                          Si   P   S  Cl  Ar 
 11  12  13                                          14  15  16  17  18 
  K  Ca  Sc  Ti   V  Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr 
 19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36 
 Rb  Sr   Y  Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te   I  Xe 
 37  38  39  40  41  42  43  44  45  46  47  48  49  50  51  52  53  54 
 Cs  Ba  *   Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn 
 55  56      72  73  74  75  76  77  78  79  80  81  82  83  84  85  86 
 Fr  Ra  *   Rf  Db  Sg  Bh  Hs  Mt  Ds  Rg  Cn  Nh  Fl  Mc  Lv  Ts  Og 
 87  88     104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 

      *  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu 
         57  58  59  60  61  62  63  64  65  66  67  68  69  70  71 
      *  Ac  Th  Pa   U  Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr 
         89  90  91  92  93  94  95  96  97  98  99 100 101 102 103 


Cs,  cesium.
Z = 55;  A = 133 (default)

Electron config: [Xe],6s1   (guess)
 = 1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 | 6s1 

Isotpe data:
Cs-133 (Z=55, A=133)
r_rms = 4.8041, c = 5.67073, mu = 2.5778, I = 3.5, parity = 1
```

---------------

Besides Hartree-Fock, other methods are available:

* Hartree: exchange term excluded
* KohnSham: Uses Kohn-Sham density functional, which includes a localised exchange potential and Latter correction (nb: valence electron should be included into the 'core' here, i.e. `core=[Xe]:6s1')
  * The colon tells ampsci to include 6s1 into the core potential, but above the 'Fermi leve; (i.e., also allow it in the valence).
* ApproxHF: Approximate (localised) Hartree-Fock potential (only for tests)
* Local: Uses a local parametric potential (Green's potential)

e.g.,

```java
HartreeFock{
  core = [Xe]:6s1;
  valence = 6sp5d;
  method = KohnSham;
}
```

---------------

### Nuclear Potential <a name="nuclear"></a>

If no innputs are given, the code will assume a Fermi-like distribution for the nuclear charge, and look up default values for the nuclear parameters (charge radius and skin thickness). The default values are chosen from the specified isotope in the `Atom{}` block. These may also be specified specifically:

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
* type: Nuclear charge distribution. The available nuclear types are: Fermi, spherical, point-like, and Gaussian, with Fermi being the default.
  * c: nuclear half-density radius (only meaningful for Fermi distribution). Usually this is not set; if it is set, it will _override_ rrms: \f$3c^2=5r_{\rm rms}^2-7\pi^2a^2\f$, and \f$t=4a~{\rm ln}3\f$
  * t: nuclear ``skin thickness'' (90-10% fall-off radius), in fm. Only for Fermi distribution.
  * \f$3c^2=5r_{\rm rms}^2-7\pi^2a^2\f$, and \f$t=4a~{\rm ln}3\f$

Default rms values are taken from:

* Angeli, Marinova, [At. Data Nucl. Data Tables **99**, 69 (2013)](http://dx.doi.org/10.1016/j.adt.2011.12.006).

It's typically best to leave these as the default parameters, but you're free to update them, e.g., the following will assume a spherically-symmetric nucleus with rms charge radius of 3.5 fm:

```java
Nucleus{
  rrms = 3.5;
  type = spherical;
}
```

--------------

### Breit <a name="breit"></a>

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

This can be included at the Hartree-Fock level with the `HartreeFock{Breit = 1.0;}` setting.

```java
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
  Breit = 1.0;
}
```

This setting is a scaling factor; i.e., setting `Breit = 0.5;` will add effective factor of 0.5 in front of \f$h^B\f$. Typically, only 0 or 1 is set; other values are useful for checking for non-physical non-linear-in-Breit effects.

**Note:** it's very important when including Breit effects along with correlations that the Breit Hamiltonian is included in the Dirac equation when forming the basis used to construct the correlation potential! This is to ensure basis is orthogonal to the Hartree-Fock core states
At the moment, there is no way to include the Breit effect into the Hartree-Fock Green's function; as a result, we cannot include Breit at the level of all-order correlations (see below).

---------------

### Radiative QED <a name="qed"></a>

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

You can all set detailed QED options within the `RadPot{}` Block.
If the `RadPot{}` block is included, then the QED radiative potential will be included with the default parameters.

```java
RadPot{}
```

equivalent to:

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

-----------
---------------

## Constructing a basis <a name="basis"></a>

We often require a full "compete" set of solutions to the Hartree-Fock equation in order to perform sums-over-states required in perturbation theory. To form these, the set of atomic orbitals are expanded as

\f[
  \phi_{n\kappa}(\vec{r}) = \sum_i^{2N} p_i S_i(\vec{r}),
\f]

where \f$S_i\f$ are a set of \f$2N\f$ basis orbitals that form an approximately complete set over a sub-domain of the radial grid \f$[0,r_{\rm max}]\f$ (\f$N\f$ is defined this way because of the duel set of positive/negative energy Dirac solutions).
The \f$p_i\f$ expansion coefficients are found by diagonalising the set of basis orbitals with respect to the Hamiltonian matrix, equivalent to solving the eigenvalue problem:
\f[
  \langle{S_i}|h_{\rm HF} |S_j \rangle p_i = \varepsilon\langle S_i|S_j \rangle p_i.
\f]
There are \f$2N\f$ solutions of eigenvalues \f$\varepsilon\f$ with corresponding eigenvectors \f$\vec{p}\f$, which correspond to the spectrum of stationary states; \f$N\f$ of these correspond to negative-energy (\f$\varepsilon<-mc^2\f$) states.

For the \f$S_i\f$ basis orbitals, we use the Duel-Kinetic-Balence basis of Beloy and Dereviano, and \f$S_i\f$ is built from \f$N_{\rm spl}\f$ B-splines of order \f$k\f$.

* Beloy, Derevianko, [Comput. Phys. Commun. **179**, 310 (2008)](https://linkinghub.elsevier.com/retrieve/pii/S0010465508001148).

```java
// Available Basis options/blocks
Basis{
  // Number of splines used in expansion [30]
  number;
  // order of splines ~7-9 [7]
  order;
  // minimum cavity radius (first internal knot) [1.0e-4]
  r0;
  // Select cavity radius r0 for each l by position where |psi(r0)/psi_max|
  // falls below r0_eps [0.0]
  r0_eps;
  // maximum cavity radius [40.0]
  rmax;
  // states to keep (e.g., 30spdf20ghi)
  states;
  // Force orthogonal to core [false]
  orthogonalise;
  // Print all spline energies (for testing) [false]
  print;
  // Include -ve energy states [false]]
  positron;
  // Derevianko (DKB) or Johnson [Derevianko]
  type;
}
```

`r0` is the location of the first 'internal' knot (first actual knot is always placed at r=0).
We may instead set `r0_eps`, which automatically choses `r0` such that the core density \f$\rho_l(r) = \sum_n|\phi_{nl}(r)|^2\f$ drops below given relative value \f$\rho_l(r0)/\rho_{\rm max}<{\rm r0\_{\rm eps}}\f$.

e.g, to calculate basis including up to \f$n=30\f$ for states up to \f$l=6\f$ in a cavity of 40 \f$a_0\f$, using 40 B-splines of order 7:

```java
Atom{
  Z = Cs;
}
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
}
Basis{
  number = 40;
  order = 7;
  r0 = 1.0e-4;
  rmax = 40.0;
  states = 30spdfghi;
}
```

If Breit and/or QED is set (see above), these are included into Hamiltonian and thus into the basis also.

It's important to look at the output when generating the basis.
In the above case (with all other default options for Cs), we get output like:

```text
Constructing B-spline basis with N=40, k=7. Storing: 30spdfghi
Using Derevianko (Duel Kinetic Balance) type splines.
Spline cavity l=0 s: (1.0e-04, 40.0)aB.
Spline cavity l=1 p: (1.0e-04, 40.0)aB.
Spline cavity l=2 d: (1.0e-04, 40.0)aB.
Spline cavity l=3 f: (1.0e-04, 40.0)aB.
Spline cavity l=4 g: (1.0e-04, 40.0)aB.
Spline cavity l=5 h: (1.0e-04, 40.0)aB.
Spline cavity l=6 i: (1.0e-04, 40.0)aB.
Basis/core:
 |<3s+|3s+>-1| = 2.9e-06
 dE/E(4s+)     = 3.6e-05
 <3s+|22s+>    = 6.7e-04
Basis/valence:
 |<5d-|5d->-1| = 9.8e-08
 dE/E(6s+)     = 1.7e-05
 <5d-|17d->    = 1.3e-04
Basis: T = 344.00 ms
 ```

The line `Spline cavity l=0 s: (1.0e-04, 40.0)aB.` tells us the actual point of the first (internal) and final B-spline knot.

The most important output is the `Basis/core:` output - this tells us how orthogonal the generated basis is to the Hartree-Fock core. The MBPT formalism relies on the orthogonality here, so it's important to check. If the results are not good enough (the code will warn you with `**`), try increasing `rmax` and/or `number`.

In this case, the worst normality occurred for the 3s state, where the inner-product of the finite-difference (Hartree Fock) \f$3s\f$ state and the corresponding basis \f$3s\f$ state was different from 1 by parts in \f$10^6\f$.
The worst energy comparison was for \f$4s\f$ state (parts in \f$10^5\f$), and the worst orthogonality comparison was for the \f$3s\f$ finite-difference Hartree-Fock state, and the \f$22s\f$ basis state, which were orthogonal to parts in \f$10^4\f$.

`Basis/valence:` Gives the same info, but for the valence states. In some cases, this is much less important, but on other cases, it matters. This is question of the physics approximation.

### Spectrum <a name="spectrum"></a>

The spectrum is the same as the basis (takes the same options), except that it also includes correlations (see below):

\f[
  \langle{S_i}|\hat h_{\rm HF} + \hat \Sigma|{S_j}\rangle p_i = \varepsilon\langle{S_i|S_j}\rangle p_i.
\f]

We typically use Basis to calculate \f$\Sigma\f$ (or other MBPT corrections), and then use spectrum to calculate atomic properties where a direct sum-over-states is required (e.g., polarisabilities).
Since we often only use basis to calculate MBPT corrections, it doesn't need to be very large.
We directly use the spectrum, so we typically require a larger basis set.
It's now typically also crucial for the spectrum to be orthogonal to the valence states (not just the core states).
e.g.,

```java
Spectrum{
  number = 90;
  order = 9;
  r0_eps = 1.0e-3;
  rmax = 90.0;
  states = 80spd;
}
```

-----------

## MBPT: second-order correlations <a name="secondorder"></a>

\f[
  H = \sum_i h_{\rm HF}(\vec{r_i}) + \delta V_{\rm corr},
\f]

where \f$h_{\rm HF}(\vec{r}_i)\f$ is the single-particle HF Hamiltonian, and

\f[
  \delta V_{\rm corr} = \sum_{i < j}\frac{1}{r_{ij}} - \sum_i V_{\rm HF}(\vec{r_i})
\f]

\f[
  \delta \varepsilon_v =
  \sum_{amn}
    \frac{g_{vamn}\widetilde g_{nmav}}{\varepsilon_v+\varepsilon_a - \varepsilon_m-\varepsilon_n}
  +\sum_{abn}
    \frac{g_{vnab}\widetilde g_{banv}}{\varepsilon_v+\varepsilon_n-\varepsilon_a-\varepsilon_b}  ,
\f]

We define _correlation potential_, \f$\Sigma\f$:

\f[
  \delta \varepsilon_v =\langle v | \Sigma |v \rangle
\f]

We then solve the Hartree-Fock equation, including the correlation potential. The solutions are known as _Brueckner orbitals_, and include correlation corrections.

\f[
  [h_{\rm HF} +  \Sigma(\varepsilon)]\phi^{\rm Br} = \varepsilon^{\rm Br}\phi^{\rm Br}
\f]

There are quite a few options available for Correlation corrections (see `ampsci -a Correlations`) -- we will focus on the most important here.

```java
Correlations{
  n_min_core = 3;
  read;
  write;
  each_valence = true;
}
```

The `read` and `write` option are filenames the correlation potential will be written to and/or read from (this saves time).
If none are given, it will read/write from the default filename.
Set to 'false' to not read/write.

`n_min_core` is minimum n for core states to include in polarisation loops.
Very low core states contribute very little (due to large excitation energy), so this saves much time for almost no degradation in results.

If `each_valence` is true, a new correlation potential will be calculated for each valence state, fully taking energy dependence into account; if false, only 1 \f$\Sigma\f$ for \f$\kappa\f$ will be formed (i.e. assumes the same energy for calculating \f$\Sigma\f$ for each angular symmetry).

```java
Correlations{
  // .... above options
  fitTo_cm;
  lambda_kappa;
}
```

`fitTo_cm` takes a list of experimental energies (in \f${\rm cm}^{-1}\f$). A scaling factor \f$\lambda\f$ will be introduced in front of the correlation potential, and will be tuned so that energies exactly match those given. This is a semi-empirical method for accounting for higher-order correlations in wavefunctions.
`lambda_kappa` is for manually setting the scaling factors.
Note that for both of these, the input list must be in the exact same order as the valence states - check the 'valence' output.

### Excplicit energy dependence

```java
Correlations{
  // .... above options
  // Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127;
  // 7s+=-0.552;}. Blank => HF energies. Takes precedence over each_valence.
  // [blank]
  ek{}
}
```

Will form the correlation potential \f$\Sigma(\varepsilon)\f$ for the listed states at the given energy \f$\varepsilon\f$.

```java
Correlations{
  // .... above options
  // Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127;
  // 7s+=-0.552;}. Blank => HF energies. Takes precidence over each_valence.
  // [blank]
  ek{
    6s+=-0.127;
    6p-=-0.127;
  }
}
```

Will use the correlation potential for \f$\kappa=-1\f$ and \f$\varepsilon = -0.127\f$ for the \f$6s\f$ state, and \f$\kappa=+1\f$ and \f$\varepsilon = -0.127\f$ for the \f$6p_{1/2}\f$ state.

-----------

## MBPT: all-order correlations <a name="allorders"></a>

The most important corrections beyond the second-order correlation potential are the screening of the residual Coulomb interaction by the core electrons and hole-particle interaction.
These are taken into account using the _all-orders correlation potential method_, also called Perturbation Theory in the Screened Coulomb Interaction (PTSCI), developed by Dzuba, Flambaum and Sushkov [1-3].

1. V. A. Dzuba, V. V. Flambaum, P. G. Silvestrov, and O. P. Sushkov, [Phys. Lett. A 131, 461 (1988).](https://dx.doi.org/10.1016/0375-9601(88)90302-7)
2. V. A. Dzuba, V. V. Flambaum, A. Y. Kraftmakher, and O. P. Sushkov, [Phys. Lett. A 142, 373 (1989).](https://dx.doi.org/10.1016/0375-9601(89)90385-X)
3. V. A. Dzuba, V. V. Flambaum, and O. P. Sushkov, [Phys. Lett. A 140, 493 (1989).](https://dx.doi.org/10.1016/0375-9601(89)90129-1)

The starting point uses the Feynman technique, in which the direct part of the correlation potential can be expressed

\f[
  \Sigma_{\rm d} = \int\frac{{\rm d}\omega}{2\pi} G_{12}(\varepsilon+\omega) Q_{1i}\Pi_{ij}(\omega) Q_{j2}(\omega),
\f]

where \f$G_{12} = G(r_1,r_2)\f$ in the Hartee-Fock Feynman Green's function, \f$Q\f$ is the (non-relativistic) Coulomb operator, and \f$\Pi\f$ is the polarisation operator
(subscripts are coordinate indices; integration is assumed over internal \f$i\f$ and \f$j\f$). Note that \f$\Pi\f$ represents polarisation of the atomic core.

The screening of the coulomb interaction can be represented by a series of polarisation corrections

\f[
  \widetilde Q \equiv  Q + \ Q(-i\, \Pi  Q) +  Q(-i\, \Pi  Q)^2+\ldots
\f]

which can be summed exactly:

\f[
  \widetilde Q(\omega) = Q\left[1+i\,\Pi(\omega) Q\right]^{-1}.
\f]

The screening is accounted for in the direct diagrams via \f$Q\to \widetilde Q\f$ for _one_ of the Coulomb operators in \f$\Sigma_{\rm d}\f$.
(Screening for exchange diagrams is taken into account in a simpler fasion, see pdf for details).

The hole-particle interaction arises due to the deviation of the Hartree-Fock potential for the excited core electron in the polarisation loop from that for the non-excited one.
The potential that simultaneously describes the occupied core and excited states is
\f[
  \hat V = V^{N-1} - (1-\hat P_{\rm core})V_{\rm self} (1-\hat P_{\rm core}),
\f]
where \f$P_{\rm core}\f$ is the operator of projection onto the core, and \f$V_{\rm self}\f$ is the self-interaction part of the Hartree-Fock potential for the outgoing electron.
Therefore, hole-particle interaction is accounted for by using this potential when forming the polarisation operator.

To include all-orders correlations, the `Feynman`, `screening`, and `holeParticle` options in the `Correlations{}` block should be set to `true`.
Alternativel, just set `AllOrder=true;`

```java
Correlations{
  n_min_core = 3;
  each_valence = true;
  Feynman = true;
  screening = true;
  holeParticle = true;
}
```

equivilant to:

```java
Correlations{
  n_min_core = 3;
  each_valence = true;
  AllOrder = true;
}
```

More options are available, though they rarely need to be changed from the default.
See [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full details.


-----------

## MBPT for matrix elements: Core Polarisation/RPA <a name="rpa"></a>

Core polarisation (RPA) is included in the matrix elements.

The best method to use is TDHF, which is numerically stable, and includes contribution from negative energy states automatically.

-----------

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
