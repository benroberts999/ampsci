\page using_operators Operators

\brief Available operators and their options: examples for E1, hyperfine, PNC

Modules that calculate matrix elements or amplitudes take an `operator`
option (or `t`/`s` in \ref module_secondorder), naming any of the available
operators, with operator-specific options given in the accompanying
`options{}` block (`t_options{}`/`s_options{}` in secondOrder).

* Get a list of available operators: `./ampsci -o`
* List the options for a given operator: `./ampsci -o <OperatorName>`
  (output is in input-file format: copy+paste it into your input file)
* Operators are constructed by \ref DiracOperator::generate() - see the
  DiracOperator namespace documentation for the definitions.

# Common Operators

\brief Examples for common operators: E1, hyperfine (hfs), and parity violation (pnc).

## E1: electric dipole

Length form by default (operator `E1`); the velocity (momentum) form is a
separate operator, `E1v`, which is frequency dependent (the modules evaluate
it at each transition frequency automatically). No options:

```java
Module::matrixElements{
  operator = E1;
  rpa = true;
}
```

## hfs: magnetic dipole (and higher) hyperfine

Hyperfine structure constants (A, B, ...). By default, nuclear parameters
(moment, spin) are looked up for the default isotope; the nuclear
magnetisation distribution is set with `F`:

```java
Module::matrixElements{
  operator = hfs;
  options{
    // pointlike, ball, shell, SingleParticle, ...
    F = ball;
    // Take mu, I, k from the isotope table by default; k=1: magnetic dipole
    // mu = 2.582;
    // I = 3.5;
    // k = 1;
  }
}
```

By default the module prints hyperfine constants (in MHz) rather than
reduced matrix elements: see the `what` option of
\ref module_matrixelements.

## pnc: nuclear-spin-independent parity violation

The weak-interaction (nuclear spin-independent) PNC operator, with a Fermi
nuclear density. Units are \f$ i (-Q_W/N) 10^{-11} \f$ a.u. Typically used
as the static operator in \ref module_secondorder :

```java
Module::secondOrder{
  v = 6s+;
  w = 7s+;
  t = E1;
  s = pnc;
  // s_options{ c = 5.67073; t = 2.3; }
}
```

Available options (from `./ampsci -o pnc`):

```java
pnc{
  c;
    // Half-density radius for Fermi rho(r). [defaut: from
    // wavefunction]
  t;
    // skin thickness [2.3]
  N;
    // Neutron number, for units [default: from
    // wavefunction]
  print;
    // Write details to screen [true]
}
```

## See also

Writing your own custom operators:

* \subpage modules_custom_operator - \copybrief modules_custom_operator
