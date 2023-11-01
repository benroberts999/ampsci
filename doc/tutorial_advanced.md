# Advanced Tutorial: MBPT

\brief Adanced ampsci tutorial: correlations

[[Home](/README.md)]

This assumes you already have ampsci compiled and have a basic understanding of how to run and use it.

* See [doc/tutorial.md](/doc/tutorial.md) for the basic tutorial
* See [ampsci.dev/ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full descrition of the physics

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

The full atomic Hamiltonian for an $N$-electron atom:
$\begin{equation}
  H = \sum_i^N h_0({r}_i) + \sum_{i<j}\frac{1}{{r}_i - {r}_j}
\end{equation}$
$\begin{equation}
  h_0(\vec r_i) = c \vec\alpha_i\cdot\vec p_i + c^2 (\beta_i-1) + V_{\rm nuc}.
\end{equation}$
The starting approximation is the Hartree-Fock method:
$\begin{equation}
  H \approx \sum_i^N [h_0(\vec r_i) + v^{\rm HF} ].
\end{equation}$

### Hartree-Fock <a name="hf"></a>

We focus on case of single-valence systems, and start with so-called $V^{N-1}$ approximation, in which Hartree-Fock potential is due to the $N-1$ core electrons:
$\begin{equation}
  \hat v^{\rm HF}\phi_a(\vec{r}_1) = \sum_{i\neq a}^{N_c}\Bigg(
  \int \frac{\phi_i^\dagger(\vec{r}_2)\phi_i(\vec{r}_2)}{|{r}_{12}|}d^3\vec{r}_2\,\phi_a(\vec{r}_1)
  -\int \frac{\phi_i^\dagger(\vec{r}_2)\phi_a(\vec{r}_2)}{|{r}_{12}|}d^3\vec{r}_2\,\phi_i(\vec{r}_1)
  \Bigg),
\end{equation}$

First, Hartree-Fock equations are solved self-consistantly for all core electrons, then the valence states are found in the Frozen Hartree-Fock potential due to the core.
For example, to run Hartree-Fock for Cs with a $N-1$ (implies Xe-like core),

```java
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
}
```

Besides Hartree-Fock, other methods are available:

* Hartree: exchange term excluded
* KohnSham: Uses Kohn-Sham density functional, which includes a localised exchange potential and Latter correction (nb: valence electron should be included into the 'core' here, i.e. `core=[Xe],6s1')
* ApproxHF: Approximate (localised) Hartree-Fock potential (only for tests)
* Local: Uses a local parametric potential (Green's potential)

e.g.,

```java
HartreeFock{
  core = [Xe],6s1;
  valence = 6sp5d;
  method = KohnSham;
}
```

### Nuclear Potential <a name="nuclear"></a>

If no innputs are given, the code will assume a Fermi-like distribution for the nuclear charge, and look up default values for the nuclear parameters (charge radius and skin thickness). The default values are chosen from the specified isotope in the `Atom{}` block. These may also be specified specifically:

```java
Nucleus{
  rrms;
  type;
  c;
  t;
}
```

* rrms: nuclear root-mean-square charge radius, in fm
* type: Nuclear charge distribution. The available nuclear types are: Fermi, spherical, pointlike, and Gaussian, with Fermi being the default.
  * c: nuclear half-density radius (only meaningful for Fermi distribution). Usually this is not set; if it is set, it will _override_ rrms: $3c^2=5r_{\rm rms}^2-7\pi^2a^2$, and $t=4a~{\rm ln}3$
  * t: nuclear ``skin thickness'' (90-10% fall-off radius), in fm. Only for Fermi distro.
  * $3c^2=5r_{\rm rms}^2-7\pi^2a^2$, and $t=4a~{\rm ln}3$

Default rms values are taken from:

* Angeli, Marinova, [At. Data Nucl. Data Tables **99**, 69 (2013)](http://dx.doi.org/10.1016/j.adt.2011.12.006).

It's typically best to leave these as the default parameters, but you're free to update them, e.g., the following will assume a spherically-symmetric nucleus with rms charge radius of 3.5 fm:

Nucleus{
  rrms = 3.5;
  type = spherical;
}

### Breit <a name="breit"></a>

The Breit Hamiltonian accounts for magnetic interactions between electrons (also known as the Gaunt interaction), and retardation effects.
It leads to a correction to the electron-electron Coulomb term in the many-body Hamiltonian:
$\begin{equation}
  \sum_{ij}\frac{1}{r_{ij}}
  \to
  \sum_{ij}\left( \frac{1}{r_{ij}} + \hat h^B_{ij}\right),
\end{equation}$
where, in the limit of zero frequency, the two-particle Breit Hamiltonian is
$\begin{equation}
  \hat h^B_{ij} = - \frac{\vec{\alpha}_i\cdot\vec{\alpha}_j + (\vec{\alpha}_i\cdot\hat{n}_{ij})(\vec{\alpha}_j\cdot\hat{n}_{ij})}{2\, {r}_{ij}}.
\end{equation}$

This can be included at the Hartree-Fock level with the `HartreeFock{Breit = 1.0;}` setting.

```java
HartreeFock{
  core = [Xe];
  valence = 6sp5d;
  Breit = 1.0;
}
```

This setting is a scaling factor; i.e., setting `Breit = 0.5;` will add effective factor of 0.5 in front of $h^B$. Typically, only 0 or 1 is set; other values are useful for checking for non-physical non-linear-in-Breit effects.

**Note:** it's very important when including Breit effects along with correlations that the breit Hamiltonian is included in the Dirac equation when forming the basis used to construct the correlation potential! This is to ensure basis is orthogonal to the Hartree-Fock core states
At the moment, there is no way to include the Breit effect into the Hartree-Fock Green's function; as a result, we cannot include Breit at the level of all-order correlations (see below).

### Radiative QED <a name="qed"></a>

Radiative QED corrections can be included into the wavefunctions using the Flambaum-Ginges radiative potential method.
An effective potential, $V_{\rm rad}$,
is added to the Hamiltonian before the equations are solved.
The potential can be written as the sum of the Uehling (vacuum polarisation), Wichmann-Kroll (higher-order vacuum polarisation) and self-energy potentials; the self-energy potential itself is written as the sum of the high- and low-frequency electric contributions, and the magnetic contribution:
$\begin{equation}
  V_{\rm rad}(\vec{r}) = V_{\rm Ueh}(r) + V_{\rm WK}(r) + V_{\rm SE}^{h}(r) +  V_{\rm SE}^{l}(r) + i (\vec{\gamma}\cdot\hat{n}) V^{\rm mag}(r).
\end{equation}$

* Flambaum, Ginges, [Phys. Rev. A **72**, 052115 (2005)](http://link.aps.org/doi/10.1103/PhysRevA.72.052115).

If the `RadPot{}` block is inlcuded, then the QED radiative potential will be included with the default parameters.

```java
RadPot{}
```

equivilant to:

```java
RadPot{
  Ueh = 1.0;
  SE_h = 1.0;
  SE_l = 1.0;
  SE_m = 1.0;
  WK = 0.0;
  core_qed = true;
  scale_l;
  scale_rN;
}
```

* Each of the (`Ueh, SE_h, SE_l, SE_m, WK`) options are scaling factors for their corresponding terms in the potential; typically 0 or 1, but can be tuned (for testing).

* `core_qed` means radiative potential will be included into the core states - this gives an important contribution and should only be set to false for testing.
* `scale_l` takes a list of scaling factors for each l; these will rescale $V_{\rm rad}$ for each partial wave.
e.g., `scale_l = 0,1,0;` will include $V_{\rm rad}$ for p states, but not s or d states.
* `scale_rN;` Re-scales the effective nuclear radius; used to test finite-nuclear size effects on $V_{\rm rad}$. =1 means normal, =0 means assume pointlike nucleus (when calculating $V_{\rm rad}$).

-----------

## Constructing a basis <a name="basis"></a>

We often require a full "compete" set of solutions to the Hartree-Fock equation in order to permorm sums-over-states required in perturbation theory. To form these, the set of atomic orbitals are expanded as
$\begin{equation}
  \phi_{n\kappa}(\vec{r}) = \sum_i^{2N} p_i S_i(\vec{r}),
\end{equation}$
where $S_i$ are a set of $2N$ basis orbitals that form an approximately complete set over a sub-domain of the radial grid $[0,r_{\rm max}]$ ($N$ is defined this way because of the duel set of positive/negative energy Dirac solutions).
The $p_i$ expansion coefficients are found by diagonalising the set of basis orbitals with respect to the Hamiltonian matrix, equivalent to solving the eigenvalue problem:
$\begin{equation}
  \langle{S_i}|h_{\rm HF} |S_j \rangle p_i = \varepsilon\langle S_i|S_j \rangle p_i.
\end{equation}$
There are $2N$ solutions of eigenvalues $\varepsilon$ with corresponding eigenvectors $\vec{p}$, which correspond to the spectrum of stationary states; $N$ of these correspond to negative-energy ($\varepsilon<-mc^2$) states.

For the $S_i$ basis orbitals, we use the Duel-Kinetic-Balence basis of Beloy and Dereviano, and $S_i$ is built from $N_{\rm spl}$ B-splines of order $k$.

* Beloy, Derevianko, [Comput. Phys. Commun. **179**, 310 (2008)](https://linkinghub.elsevier.com/retrieve/pii/S0010465508001148).

```java
// Available Basis options/blocks
Basis{
  // Number of splines used in expansion [0]
  number;
  // order of splines ~7-9 [7]
  order;
  // minimum cavity radius (first internal knot) [1.0e-4]
  r0;
  // Select cavity radius r0 for each l by position where |psi(r0)/psi_max|
  // falls below r0_eps [1.0e-3]
  r0_eps;
  // maximum cavity radius [Grid{rmax}]
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
Typically, we instead set `r0_eps`, which automatically choses `r0` such that the core density $\rho_l(r) = \sum_n|\phi_{nl}(r)|^2$ drops below given relative value $\rho_l(r0)/\rho_{\rm max}<{\rm r0\_{\rm eps}}$.

e.g, to calculate basis including up to $n=30$ for states up to $l=6$ in a cavity of 40 $a_0$, using 40 B-splines of order 7:

```java
Basis{
  number = 40;
  order = 7;
  r0_eps = 1.0e-3;
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
Spline cavity l=0 s: (1.4e-04, 40.0)aB.
Spline cavity l=1 p: (2.8e-03, 40.0)aB.
Spline cavity l=2 d: (2.4e-02, 40.0)aB.
Spline cavity l=3 f: (2.4e-02, 40.0)aB.
Spline cavity l=4 g: (2.4e-02, 40.0)aB.
Spline cavity l=5 h: (2.4e-02, 40.0)aB.
Spline cavity l=6 i: (2.4e-02, 40.0)aB.
Basis/core:
 |<3s+|3s+>-1| = 2.4e-06
 dE/E(5s+)     = 4.6e-05
 <3s+|22s+>    = 4.5e-04
Basis/valence:
 |<5d-|5d->-1| = 9.5e-08
 dE/E(6s+)     = 2.9e-05
 <5d-|23d->    = 1.1e-04
 ```

The line `Spline cavity l=0 s: (1.4e-04, 40.0)aB.` tells us the actual point of the first (internal) and final B-spline knot.

The most important output is the `Basis/core:` output - this tells us how orthogonal the generated basis is to the Hartree-Fock core. The MBPT formalism relies on the orthogonality here, so it's important to check. If the results are not good enough (the code will warn you with `**`), try increasing `rmax` and/or `number`.

In this case, the worst normality occured for the 3s state, where the inner-product of the finite-difference (Hartree Fock) $3s$ state and the corresponding basis $3s$ state was different from 1 by parts in $10^6$.
The worst energy comparison was for $5s$ state (parts in $10^5$), and the worst orthogonality comparison was for the $3s$ finite-difference Hartree-Fock state, and the $22s$ basis state, which were orthogonal to parts in $10^4$.

### Spectrum <a name="spectrum"></a>

The spectrum is the same as the basis (takes the same options), except that it also includes correlations (see below):

$\begin{equation}
  \langle{S_i}|\hat h_{\rm HF} + \hat \Sigma|{S_j}\rangle p_i = \varepsilon\langle{S_i|S_j}\rangle p_i.
\end{equation}$

We typically use Basis to calculate $\Sigma$ (or other MBPT corrections), and then use spectrum to calculate atomic properties where a direct sum-over-states is required (e.g., polarisabilities).
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

$\begin{equation}
  H = \sum_i h_{\rm HF}(\vec{r}_i) + \delta V_{\rm corr},
\end{equation}$
where $h_{\rm HF}(\vec{r}_i)$ is the single-particle HF Hamiltonian, and
$\begin{equation}
  \delta V_{\rm corr} = \sum_{i<j}\frac{1}{{r}_{ij}} - \sum_iV_{\rm HF}(\vec{r}_i)
\end{equation}$

$\begin{equation}
  \delta \varepsilon_v =
  \sum_{amn}
    \frac{g_{vamn}\widetilde g_{nmav}}{\varepsilon_v+\varepsilon_a - \varepsilon_m-\varepsilon_n}
  +\sum_{abn}
    \frac{g_{vnab}\widetilde g_{banv}}{\varepsilon_v+\varepsilon_n-\varepsilon_a-\varepsilon_b}  ,
\end{equation}$

We define _correlation potential_, $\Sigma$:
$\begin{equation}
  \delta \varepsilon_v =\langle v | \Sigma |v \rangle
\end{equation}$

We then solve the Hartree-Fock equation, including the correlation potential. The solutions are known as _Brueckner orbitals_, and include correlation corrections.
$\begin{equation}
  [h_{\rm HF} +  \Sigma(\varepsilon)]\phi^{\rm Br} = \varepsilon^{\rm Br}\phi^{\rm Br}
\end{equation}$

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

`n_min_core` is minumum n for core states to include in polarisation loops.
Very low core states contribute very little (due to large excitation energy), so this saves much time for almost no degredation in results.

If `each_valence` is true, a new correlation potential will be calculated for each valence state, fully taking energy dependence into account; if false, only 1 $\Sigma$ for $\kappa$ will be formed (i.e. assumes the same energy for calculating $\Sigma$ for each angular symmetry).

```java
Correlations{
  // .... above options
  fitTo_cm;
  lambda_kappa;
}
```

`fitTo_cm` takes a list of experimental energies (in ${\rm cm}^{-1}$). A scaling factor $\lambda$ will be introduced in front of the correlation potential, and will be tuned so that energies exactly match those given. This is a semi-empirical method for accounting for hiher-order correlations in wavefunctions.
`lambda_kappa` is for manually setting the scaling factors.
Note that for both of these, the input list must be in the exact same order as the valence states - check the 'valence' output.

### Excplicit energy dependence

```java
Correlations{
  // .... above options
  // Block: Explicit list of energies to solve for. e.g., ek{6s+=-0.127;
  // 7s+=-0.552;}. Blank => HF energies. Takes precidence over each_valence.
  // [blank]
  ek{}
}
```

Will form the correlation potential $\Sigma(\varepsilon)$ for the listed states at the given energy $\varepsilon$.

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

Will use the correlation potential for $\kappa=-1$ and $\varepsilon = -0.127$ for the $6s$ state, and $\kappa=+1$ and $\varepsilon = -0.127$ for the $6p_{1/2}$ state.

-----------

## MBPT: all-order correlations <a name="allorders"></a>

The most important corrections beyond the second-order correlation potential are the screening of the residual Coulomb interaction by the core electrons and hole-particle interaction.
These are taken into account using the _all-orders correlation potential method_, also called Perturbation Theory in the Screened Couloumb Interaction (PTSCI), developed by Dzuba, Flambaum and Sushkov [1-3].

1. V. A. Dzuba, V. V. Flambaum, P. G. Silvestrov, and O. P. Sushkov, [Phys. Lett. A 131, 461 (1988).](https://dx.doi.org/10.1016/0375-9601(88)90302-7)
2. V. A. Dzuba, V. V. Flambaum, A. Y. Kraftmakher, and O. P. Sushkov, [Phys. Lett. A 142, 373 (1989).](https://dx.doi.org/10.1016/0375-9601(89)90385-X)
3. V. A. Dzuba, V. V. Flambaum, and O. P. Sushkov, [Phys. Lett. A 140, 493 (1989).](https://dx.doi.org/10.1016/0375-9601(89)90129-1)

The starying point uses the Feynman technique, in which the direct part of the correlation potential can be expressed
$\begin{equation}
  \Sigma_{\rm d} = \int\frac{{\rm d}\omega}{2\pi} G_{12}(\varepsilon+\omega) Q_{1i}\Pi_{ij}(\omega) Q_{j2}(\omega),
\end{equation}$
where $G_{12} = G(r_1,r_2)$ in the Hartee-Fock Feynman Green's function, $Q$ is the (non-relativistic) Coulomb operator, and $\Pi$ is the polarisation operator
(subsripts are coordinate indices; integration is assumed over internal $i$ and $j$). Note that $\Pi$ represents polarisation of the atomic core.

The screening of the coulomb interaction can be represented by a series of polarisation corrections
$\begin{equation}
  \widetilde Q \equiv  Q + \ Q(-i\, \Pi  Q) +  Q(-i\, \Pi  Q)^2+\ldots
\end{equation}$
which can be summed exactly:
$\begin{equation}
  \widetilde Q(\omega) = Q\left[1+i\,\Pi(\omega) Q\right]^{-1}.
\end{equation}$
The screening is accounted for in the direct diagrams via $Q\to \widetilde Q$ for _one_ of the Coulomb operators in $\Sigma_{\rm d}$.
(Screening for exchange diagrams is taken into account in a simpler fasion, see pdf for details).

The hole-particle interaction arises due to the deviation of the Hartree-Fock potential for the excited core electron in the polarisation loop from that for the non-excited one.
The potential that simultaneously describes the occupied core and excited states is
$\begin{equation}
  \hat V = V^{N-1} - (1-\hat P_{\rm core})V_{\rm self} (1-\hat P_{\rm core}),
\end{equation}$
where $P_{\rm core}$ is the operator of projection onto the core, and $V_{\rm self}$ is the self-interaction part of the Hartree-Fock potential for the outgoing electron.
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

Note: currently, there is a numerical problem in calculating polarisation loop for very deep core states. This leads to no issues, unless you attempt to polarise a deep core shell, where large numerical errors are encountered.
Best to set `n_min_core = 2` for Rb, `n_min_core = 3` for Cs, `n_min_core = 4` for Fr etc.
The effect of this is negligable, as can be checked by seeing impact on second-order results.
This should be addressed soon.

Note also: it is currently not possible to correctly account for Breit Hamiltonian within the all-orders method (problem is in construction of Green's function).
Therefore, Breit should only be included at second-order level.

-----------

## MBPT for matrix elements: Core Polarisation/RPA <a name="rpa"></a>

Core polarisation (RPA) is included in the matrix elements.

The best method to use is TDHF, which is numerically stable, and includes contribution from negative energy states automatically.
However, it doesn't work (at the moment) for even parity operators, so for these we must use the diagram method.
The diagram method uses `Basis{}` for internal summations, so one should ensure the basis is complete enough for this to work.

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
