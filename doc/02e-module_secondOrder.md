\page module_secondorder Second-Order Amplitudes

\brief Second-order amplitudes for a valence state/transition: polarisabilities, PNC, etc.

The `secondOrder` module (run as `Module::secondOrder{}`) calculates the
second-order (in the external field) amplitude \f$ A^K \f$
for a single-valence transition \f$ v \to w \f$, for two one-body operators
\f$ t \f$ (at frequency \f$ \omega \f$) and \f$ s \f$ (at
\f$ \omega_s = E_w - E_v - \omega \f$):

\f[
  A^K = \sum_n \left[
    c_1 \frac{\langle w||t||n\rangle\langle n||s||v\rangle}{E_v+\omega_s-E_n}
  + c_2 \frac{\langle w||s||n\rangle\langle n||t||v\rangle}{E_v+\omega-E_n}
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
matrix elements of the sum via the `StructureRadiation{}` block. This is
always by sum-over-states (for either method), and printed separately: SR
and Norm, each for t and s, and the small nonlinear part (t and s both
corrected). The basis is used for all of the SR+N part (legs, intermediate
states, and the internal lines of the diagrams).

This is the single-valence analogue of the CI_secondOrder module.

## Definitions

Same definitions in the CI_secondOrder module, with \f$ j \to J \f$.
\f$ t \f$ has rank \f$ k_t \f$, \f$ s \f$ has rank \f$ k_s \f$, and
\f$ [K] \equiv 2K+1 \f$.

Uncoupled (definite projections \f$ q_1 \f$ of \f$ t \f$, \f$ q_2 \f$ of
\f$ s \f$; the sum over \f$ n \f$ includes magnetic quantum numbers):

\f[
  A^{k_tk_s}_{q_1q_2} = \sum_n \left[
    \frac{\langle w|t_{q_1}|n\rangle\langle n|s_{q_2}|v\rangle}
         {E_v+\omega_s-E_n}
  + \frac{\langle w|s_{q_2}|n\rangle\langle n|t_{q_1}|v\rangle}
         {E_v+\omega-E_n}
  \right].
\f]

Coupled to rank \f$ K \f$, with \f$ Q = q_1+q_2 = m_w-m_v \f$:

\f[
  A^K_Q = \sum_{q_1q_2}\langle k_tq_1\,k_sq_2|KQ\rangle\,A^{k_tk_s}_{q_1q_2}
        = (-1)^{k_t-k_s+Q}\sqrt{[K]}\sum_{q_1q_2}
          \begin{pmatrix} k_t & k_s & K \\ q_1 & q_2 & -Q \end{pmatrix}
          A^{k_tk_s}_{q_1q_2}.
\f]

Reduced, via the Wigner-Eckart theorem:

\f[
  A^K_Q = (-1)^{j_w-m_w}
    \begin{pmatrix} j_w & K & j_v \\ -m_w & Q & m_v \end{pmatrix}
    A^K ,
\f]

\f$ A^K \f$ being the quantity the module reports: the reduced matrix element
of \f$ [t\times s]^K \f$ (both time orderings), with

\f[
  c_1 = (-1)^{K}\sqrt{[K]}\,(-1)^{j_w+j_v}
    \begin{Bmatrix} K & k_s & k_t \\ j_n & j_w & j_v \end{Bmatrix},
  \qquad
  c_2 = (-1)^{k_t+k_s}\sqrt{[K]}\,(-1)^{j_w+j_v}
    \begin{Bmatrix} K & k_t & k_s \\ j_n & j_w & j_v \end{Bmatrix}.
\f]

Uncoupling again gives the z-component (\f$ m_v = m_w = m \f$,
\f$ q_1 = q_2 = 0 \f$), also printed:

\f[
  A_{zz} = \sum_K \langle k_t0\,k_s0|K0\rangle\,(-1)^{j_w-m}
      \begin{pmatrix} j_w & K & j_v \\ -m & 0 & m \end{pmatrix} A^K .
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
  \alpha_0 = \frac{A^0}{\sqrt{3[j_v]}}
    \quad (K=0,\ w=v),
  \qquad
  \alpha_2 = -\sqrt{\frac{2j(2j-1)}{3(j+1)(2j+1)(2j+3)}}\;A^2
    \quad (K=2,\ j_w=j_v=j\ge1),
\f]
\f[
  \beta = \frac{A^1}{\sqrt{2}\,\langle w||\boldsymbol\sigma||v\rangle}
    \quad (K=1),
\f]

with \f$ \langle w||\boldsymbol\sigma||v\rangle = 2S_{\kappa\kappa'} \f$ for a
single valence electron (radial overlap dropped by convention); for CI states
see \ref CI::sigma_rme "CI::sigma_rme()".

With \f$ t = d \f$ and \f$ s = h_W \f$ (PNC, \f$ k_s = 0 \f$, so \f$ K = 1 \f$
and \f$ \langle 1\,0\,0\,0|1\,0\rangle = 1 \f$), at \f$ m_v = m_w = m \f$:

\f[
  E_{\rm PNC} = A^1_0 = (-1)^{j_w-m}
    \begin{pmatrix} j_w & 1 & j_v \\ -m & 0 & m \end{pmatrix} A^1 .
\f]

Examples -- static polarisability of Cs 6s, and the 6s-7s PNC amplitude:

```java
Module::secondOrder{
  v = 6s+;
}

Module::secondOrder{
  v = 6s+;
  w = 7s+;
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
  // transition v -> w, for two operators t and s, at
  // frequencies omega(=omega_t) and omega_s:
  // A^K = sum_n [c1 <w||t||n><n||s||v>/(E_v+omega_s-E_n) +
  // c2 <w||s||n><n||t||v>/(E_v+omega-E_n)].
  // Energy conservation fixes omega_s = E_w - E_v - omega,
  // so only omega is an input: at its default, t carries
  // the whole transition frequency and s is static.
  // For a dynamic polarisability, set w = v and omega to
  // the frequency; then omega_s = -omega (only does 1
  // freq.)
  // The valence sum is evaluated either by sum-over-states
  // (method=SOS), or with mixed states (method=MS).
  v;
    // Initial valence state, e.g., 6s+. May be a list (v =
    // 6s+,7s+;): each pair (v_i -> w_i) is calculated in
    // turn, with a summary of the totals at the end
    // [required]
  w;
    // Final valence state(s) [default: same as v]. If the
    // lists differ in length, the shorter is padded with
    // its final state: v = a,b,c; w = d,e; gives (a,d),
    // (b,e), (c,e)
  t;
    // The operator that carries the frequency omega [E1]
  t_options{}
    // Options for the t operator
  s;
    // The other operator; it carries omega_s = E_w - E_v -
    // omega [E1]
  s_options{}
    // Options for the s operator
  omega;
    // Frequency of t. For a transition, leave as default,
    // so that t carries it all and s is static. For a
    // dynamic polarisability (w = v), set this to the
    // frequency: s then carries -omega. [default: E_w -
    // E_v]
  K;
    // Rank K of the amplitude. Requires |kt-ks| <= K <=
    // kt+ks, and the triangle rule for (jw,K,jv) [default:
    // smallest allowed]
  method;
    // SOS (sum-over-states) or MS (mixed states) [MS]
  rpa;
    // Method used for RPA: true(=TDHF), false, TDHF, basis,
    // diagram, cntm [true]. MS requires TDHF.
  DCP;
    // Include the double core polarisation (DCP)
    // correction, by TDHF. With rpa = false: gives
    // lowest-order, with no RPA for t and s. [true]
  replace_w_valence;
    // SOS: replace the basis/spectrum states with the
    // corresponding valence states in the sum, where
    // available (e.g., when the valence states were fitted
    // to experiment) [false]
  n_main;
    // The 'main' part of the valence sum: intermediate
    // states up to this n, printed separately [max_n_core +
    // 4]
  StructureRadiation{}
    // Options for structure radiation and normalisation. If
    // this block is included, SR+N is added to the
    // single-particle matrix elements of the sum, by
    // sum-over-states (for either method). The basis is
    // used for all of it (legs, intermediate states, and
    // the internal lines of the diagrams). Printed
    // separately: SR and Norm, each for t and s, and the
    // nonlinear part (t and s both corrected)

  omega_t;
    // Explicit frequency for t, overriding omega. Energy
    // conservation is then up to you; use with care [rare]
  omega_s;
    // Explicit frequency for s, overriding omega_s = E_w -
    // E_v - omega. Use with care [rare]
}

// Available StructureRadiation options/blocks
StructureRadiation{
  // If this block is included, SR + Normalisation
  // corrections will be included, by sum-over-states (for
  // either method), using the basis for all states
  n_minmax;
    // list; min,max n for core/excited (internal): [1,inf]
  n_max_legs;
    // SR+N is applied to matrix elements whose states both
    // have n <= this (the valence legs, and the low-n
    // intermediate states). SR+N is only meaningful between
    // physical states: the high-n basis states are cavity
    // states. Never applied to core intermediate states
    // [default: max_n_core + 3]
  norm;
    // Include the normalisation of states? If false, only
    // the structure radiation is included [true]
  rpa;
    // Include the core polarisation (RPA) inside SR+N: on
    // the internal lines of the structure radiation
    // diagrams, and in the normalisation. Only if RPA is
    // included in the main block [true]
  omega;
    // Structure radiation denominators at the frequency of
    // each operator (true), or at zero (false): a test of
    // the frequency dependence, usually very weak [true]
  Qk_file;
    // true/false/filename: QkTable of Coulomb integrals;
    // read if it exists, otherwise calculated and written.
    // true: default filename. false: no QkTable (~10x
    // slower, less memory) [true]
}
```
