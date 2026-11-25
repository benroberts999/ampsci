\page module_secondorder Second-Order Amplitudes

\brief Second-order amplitudes for a valence state/transition: polarisabilities, PNC, etc.

The `secondOrder` module (run as `Module::secondOrder{}`) calculates the
second-order (in the external field) amplitude \f$ A^K \f$
for a single-valence transition \f$ A \to B \f$, for two one-body operators
\f$ t \f$ (at frequency \f$ \omega \f$) and \f$ s \f$ (at
\f$ \omega_s = E_B - E_A - \omega \f$):

\f[
  A^K = \sum_n \left[
    c_1 \frac{\langle B||t||n\rangle\langle n||s||A\rangle}{E_A+\omega_s-E_n}
  + c_2 \frac{\langle B||s||n\rangle\langle n||t||A\rangle}{E_A+\omega-E_n}
  \right].
\f]

This single module covers (with the appropriate operators and rank K):

* Scalar polarisability, \f$ \alpha_0 \f$ (t = s = E1, K = 0; static,
  dynamic, or transition)
* Tensor polarisability, \f$ \alpha_2 \f$ (K = 2)
* Vector transition polarisability, \f$ \beta \f$ (K = 1)
* PNC amplitudes (t = E1, s = pnc, K = 1)

and the relevant physical quantity is printed automatically.

The valence sum is evaluated either by sum-over-states (`method = SOS;`,
over the spectrum or basis), or with mixed states (`method = MS;`, solving
the inhomogeneous equation), which is complete: no truncation of the sum.
Run twice to compare the two. The mixed-states method evaluates the sums two
independent ways (with the mixed states of s, and of t), which must agree:
a strong internal check, printed with their relative difference.

The core polarisation contribution is included (for diagonal, K = 0
amplitudes), and the valence sum is decomposed into core-valence, main
(low-n, see `n_main`), and tail parts, for comparison with other
calculations. Structure radiation + normalisation may be added to the
matrix elements of the sum via the `StructureRadiation{}` block
(sum-over-states method only).

This is the single-valence analogue of the CI_secondOrder module.

## Definitions

Same definitions in the CI_secondOrder module, with \f$ j \to J \f$.
\f$ t \f$ has rank \f$ k_t \f$, \f$ s \f$ has rank \f$ k_s \f$, and
\f$ [K] \equiv 2K+1 \f$.

Uncoupled (definite projections \f$ q_1 \f$ of \f$ t \f$, \f$ q_2 \f$ of
\f$ s \f$; the sum over \f$ n \f$ includes magnetic quantum numbers):

\f[
  A^{k_tk_s}_{q_1q_2} = \sum_n \left[
    \frac{\langle B|t_{q_1}|n\rangle\langle n|s_{q_2}|A\rangle}
         {E_A+\omega_s-E_n}
  + \frac{\langle B|s_{q_2}|n\rangle\langle n|t_{q_1}|A\rangle}
         {E_A+\omega-E_n}
  \right].
\f]

Coupled to rank \f$ K \f$, with \f$ Q = q_1+q_2 = m_B-m_A \f$:

\f[
  A^K_Q = \sum_{q_1q_2}\langle k_tq_1\,k_sq_2|KQ\rangle\,A^{k_tk_s}_{q_1q_2}
        = (-1)^{k_t-k_s+Q}\sqrt{[K]}\sum_{q_1q_2}
          \begin{pmatrix} k_t & k_s & K \\ q_1 & q_2 & -Q \end{pmatrix}
          A^{k_tk_s}_{q_1q_2}.
\f]

Reduced, via the Wigner-Eckart theorem:

\f[
  A^K_Q = (-1)^{j_B-m_B}
    \begin{pmatrix} j_B & K & j_A \\ -m_B & Q & m_A \end{pmatrix}
    A^K ,
\f]

\f$ A^K \f$ being the quantity the module reports: the reduced matrix element
of \f$ [t\times s]^K \f$ (both time orderings), with

\f[
  c_1 = (-1)^{K}\sqrt{[K]}\,(-1)^{j_B+j_A}
    \begin{Bmatrix} K & k_s & k_t \\ j_n & j_B & j_A \end{Bmatrix},
  \qquad
  c_2 = (-1)^{k_t+k_s}\sqrt{[K]}\,(-1)^{j_B+j_A}
    \begin{Bmatrix} K & k_t & k_s \\ j_n & j_B & j_A \end{Bmatrix}.
\f]

Uncoupling again gives the z-component (\f$ m_A = m_B = m \f$,
\f$ q_1 = q_2 = 0 \f$), also printed:

\f[
  A_{zz} = \sum_K \langle k_t0\,k_s0|K0\rangle\,(-1)^{j_B-m}
      \begin{pmatrix} j_B & K & j_A \\ -m & 0 & m \end{pmatrix} A^K .
\f]

**Sign convention.** The coupling above is the standard Clebsch-Gordan one,
\f$ t \f$ first. The alternative definition

\f[
  \tilde A^K_Q = (-1)^{Q}\sqrt{[K]}\sum_{q_1q_2}
    \begin{pmatrix} k_t & k_s & K \\ -q_1 & -q_2 & Q \end{pmatrix}
    A^{k_tk_s}_{q_1q_2}
  = (-1)^K A^K_Q
\f]

differs by \f$ (-1)^K \f$. It cancels in \f$ A_{zz} \f$ (so in
\f$ E_{\rm PNC} \f$) and for even \f$ K \f$; it affects only the sign of
\f$ \beta \f$.

## Specific cases

With \f$ t = s = d \f$ (E1) and \f$ [j] \equiv 2j+1 \f$:

\f[
  \alpha_0 = \frac{A^0}{\sqrt{3[j_A]}}
    \quad (K=0,\ B=A),
  \qquad
  \alpha_2 = -\sqrt{\frac{2j(2j-1)}{3(j+1)(2j+1)(2j+3)}}\;A^2
    \quad (K=2,\ j_B=j_A=j\ge1),
\f]
\f[
  \beta = \frac{A^1}{\sqrt{2}\,\langle B||\boldsymbol\sigma||A\rangle}
    \quad (K=1),
\f]

with \f$ \langle B||\boldsymbol\sigma||A\rangle = 2S_{\kappa\kappa'} \f$ for a
single valence electron (radial overlap dropped by convention); for CI states
see \ref CI::sigma_rme "CI::sigma_rme()".

With \f$ t = d \f$ and \f$ s = h_W \f$ (PNC, \f$ k_s = 0 \f$, so \f$ K = 1 \f$
and \f$ \langle 1\,0\,0\,0|1\,0\rangle = 1 \f$), at \f$ m_A = m_B = m \f$:

\f[
  E_{\rm PNC} = A^1_0 = (-1)^{j_B-m}
    \begin{pmatrix} j_B & 1 & j_A \\ -m & 0 & m \end{pmatrix} A^1 .
\f]

Examples -- static polarisability of Cs 6s, and the 6s-7s PNC amplitude:

```java
Module::secondOrder{
  A = 6s+;
}

Module::secondOrder{
  A = 6s+;
  B = 7s+;
  t = E1;
  s = pnc;
}
```

The calculations are performed by the Amplitudes library: see
\ref Amplitudes::sos_valence "Amplitudes::sos_valence()", \ref Amplitudes::ms_valence "Amplitudes::ms_valence()",
\ref Amplitudes::sos_core "Amplitudes::sos_core()", and \ref Amplitudes::ms_core "Amplitudes::ms_core()" for the formulas and conventions (the module,
Module::secondOrder, only parses input and prints).

Available options (from `./ampsci -m secondOrder`):

```java
// Available Module::secondOrder options/blocks
Module::secondOrder{
  // Second-order amplitude A^K for a single-valence
  // transition A -> B, for two operators t and s, at
  // frequencies omega and omega_s:
  // A^K = sum_n [c1 <B||t||n><n||s||A>/(E_A+omega_s-E_n) +
  // c2 <B||s||n><n||t||A>/(E_A+omega-E_n)].
  // Energy conservation fixes omega_s = E_B - E_A - omega,
  // so only omega is an input: at its default, t carries
  // the whole transition frequency and s is static. For a
  // dynamic polarisability, set B = A and omega to the
  // frequency; then omega_s = -omega.
  // The valence sum is evaluated either by sum-over-states
  // (method=SOS), or with mixed states (method=MS), which
  // is complete (no truncation of the sum). Run twice to
  // compare the two.
  A;
    // Initial valence state, e.g., 6s+ [required]
  B;
    // Final valence state [default: same as A]
  t;
    // The operator that carries the frequency omega [E1]
  t_options{}
    // Options for the t operator
  s;
    // The other operator; it carries omega_s = E_b - E_a -
    // omega [E1]
  s_options{}
    // Options for the s operator
  omega;
    // Frequency of t. For a transition, leave as default,
    // so that t carries it all and s is static. For a
    // dynamic polarisability (B = A), set this to the
    // frequency: s then carries -omega. [default: E_b -
    // E_a]
  omega_t;
    // Explicit frequency for t, overriding omega. Energy
    // conservation is then up to you; use with care [rare]
  omega_s;
    // Explicit frequency for s, overriding omega_s = E_b -
    // E_a - omega. Use with care [rare]
  K;
    // Rank K of the amplitude. Requires |kt-ks| <= K <=
    // kt+ks, and the triangle rule for (jb,K,ja) [default:
    // smallest allowed]
  method;
    // SOS (sum-over-states) or MS (mixed states) [MS]
  rpa;
    // Method used for RPA: true(=TDHF), false, TDHF, basis,
    // diagram [true]. MS requires TDHF: the diagram method
    // triggers a warning, and TDHF is used
  use_spectrum;
    // SOS: use the spectrum (which includes Sigma, if it
    // was made with it) for the sum, if available;
    // otherwise the basis is used [true]
  Sigma;
    // MS: include the correlation potential in the mixed
    // states, if it is available (use with Brueckner
    // valence states) [true]
  n_main;
    // The 'main' part of the valence sum: intermediate
    // states up to this n, printed separately (as the pnc
    // module) [max_n_core + 4]
  StructureRadiation{}
    // Options for structure radiation and normalisation. If
    // this block is included, SR+N is added to the
    // single-particle matrix elements of the sum.
    // Sum-over-states method only
}

// Available StructureRadiation options/blocks
StructureRadiation{
  // If this block is included, SR + Normalisation
  // corrections will be included (sum-over-states method
  // only)
  Qk_file;
    // true/false/filename - SR: filename for QkTable file.
    // If blank will not use QkTable; if exists, will read
    // it in; if doesn't exist, will create it and write to
    // disk. If 'true' will use default filename
  n_min_core;
    // Minimum n for core states in the diagrams [1]
  basis;
    // Basis used for the internal lines of the diagrams,
    // e.g., 30spdf [full basis]
  n_max_legs;
    // SR+N is applied to matrix elements whose states both
    // have n <= this (the valence legs, and the low-n
    // intermediate states). SR+N is only meaningful between
    // physical states: the high-n basis states are cavity
    // states [default: max_n_core + 3]
  norm;
    // Include the normalisation of states? If false, only
    // the structure radiation is included [true]
}
```
