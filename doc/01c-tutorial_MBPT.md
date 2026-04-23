\page tutorial_mbpt MBPT and Correlations

\brief MBPT and correlation corrections

This assumes you already have ampsci compiled and have a basic understanding of how to run and use it, and a basic understanding of the Hartree-Fock starting approximation.

* See [Compilation](\ref compilation) for compilation instructions
* See [Atomic Potentials](\ref tutorial_potentials) for getting started with basic calculations (Hartree-Fock etc.)
* See [ampsci.dev/ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full description of the physics

## Contents

* [Constructing a basis](#basis)
  * [Spectrum](#spectrum)
* [MBPT: second-order correlations](#secondorder)
* [MBPT: all-order correlations](#allorders)

-----

## Constructing a basis <a name="basis"></a>

We often require a full "complete" set of solutions to the Hartree-Fock equation in order to perform sums-over-states required in perturbation theory. To form these, the set of atomic orbitals (Hartree-Fock eigenstates) are expanded as

\f[
  \phi_{n\kappa}(\vec{r}) = \sum_i^{2N} p_i S_i(\vec{r}),
\f]

where \f$S_i\f$ are a set of \f$2N\f$ basis orbitals that form an approximately complete set over a sub-domain of the radial grid \f$[0,r_{\rm max}]\f$ (\f$N\f$ is defined this way because of the duel set of positive/negative energy Dirac solutions).
The \f$p_i\f$ expansion coefficients are found by diagonalising the set of basis orbitals with respect to the Hamiltonian matrix, equivalent to solving the eigenvalue problem:
\f[
  \langle{S_i}|h_{\rm HF} |S_j \rangle p_i = \varepsilon\langle S_i|S_j \rangle p_i.
\f]
There are \f$2N\f$ solutions of eigenvalues \f$\varepsilon\f$ with corresponding eigenvectors \f$\vec{p}\f$, which correspond to the spectrum of stationary states; \f$N\f$ of these correspond to negative-energy (\f$\varepsilon<-mc^2\f$) states.

* The negative energy states are dropped by default but can be kept with `positron=true;`. You can keep a smaller subset of the negative energy states with, for example,  `positron=7spdf;`.
* The negative energy states are always put at the end of the basis, and have \f$\varepsilon < -m_ec^2\f$

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
It's now typically also crucial for the spectrum to be orthogonal to the valence states (since correlations are not included into the HF core, it's expected that the spectrum is not perfectly orthogonal to the core if we are including correlations).
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

-----

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

There are quite a few options available for Correlation corrections (see `ampsci -i Correlations`) -- we will focus on the most important here.

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

-----

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

To include all-orders correlations, the `Feynman`, `screening`, and `hole_particle` options in the `Correlations{}` block should be set to `true`.
Alternativel, just set `all_order=true;`

```java
Correlations{
  n_min_core = 3;
  each_valence = true;
  Feynman = true;
  screening = true;
  hole_particle = true;
}
```

equivilant to:

```java
Correlations{
  n_min_core = 3;
  each_valence = true;
  all_order = true;
}
```

More options are available, though they rarely need to be changed from the default.
See [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full details.

-----

* See [Modules](\ref tutorial_modules) for tutorial on using the wavefunctions (e.g., calculating matrix elements)
* See physics documentation: [ampsci.pdf](https://ampsci.dev/ampsci.pdf) for full physics description of the employed methods
