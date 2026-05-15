\page modules_custom_operator Writing a custom operator

\brief How to write your own DiracOperator for use in ampsci modules

Operators in ampsci are single-particle tensor operators of rank \f$ k \f$, acting on
\ref DiracSpinor objects. They are represented by the \ref DiracOperator::TensorOperator
virtual base class. All built-in operators (E1, hfs, pnc, ...) derive from this class,
and you can add your own the same way.

* See \ref modules_custom for how to write a module that uses your operator.
* See \ref modules for the list of existing modules and operators.

---

## Why write an operator?

If you need a matrix element that is not already implemented, the operator framework gives you:

* Automatic reduced matrix elements via the Wigner-Eckart theorem.
* Selection rules enforced by the base class (rank, parity).
* Compatibility with all existing modules (e.g., `Module::MatrixElements`, RPA, structure radiation).
* Frequency dependence handled cleanly if needed.

---

## How operators work

See \ref DiracOperator::TensorOperator for the full interface.

The base class computes reduced matrix elements (RMEs) as

\f[
  \langle a \| \hat{h} \| b \rangle = A_{ab} \cdot R_{ab}
\f]

where \f$ A_{ab} \f$ is the angular factor (see `angularF()`) and the default radial integral is

\f[
  R_{ab} = c \int_0^\infty v(r)\left(
    C_{ff}\,f_a f_b + C_{fg}\,f_a g_b +
    C_{gf}\,g_a f_b + C_{gg}\,g_a g_b
  \right)\,{\rm d}r
\f]

where \f$ c \f$ is an overall constant, \f$ v(r) \f$ is an optional radial function,
and the \f$ C_{xy} \f$ are angular coefficients (defaulting to \f$ C_{ff}=C_{gg}=1 \f$,
\f$ C_{fg}=C_{gf}=0 \f$).

* The default radial integral is enough for the vast majority of operators, but for more complicated operators, it too can be overridden.

---

## Creating your own operator

Derive from \ref DiracOperator::TensorOperator and pass the operator properties to the base constructor, for example:

```cpp
class MyOperator : public DiracOperator::TensorOperator {
public:
  MyOperator(const Grid &grid, double constant, std::vector<double> v_of_r /* + any other required options*/)
      : TensorOperator(rank, DiracOperator::Parity::even, constant, v_of_r,
                       DiracOperator::Realness::real, /*freq_dep=*/false) {}

  // Mandatory: angular factor
  double angularF(int kappa_a, int kappa_b) const override {
    return /* ... */;
  }

  // Optiona
  double angularCff(int kappa_a, int kappa_b) const override {
    return /* default returns 1 */;
  }
};
```

Key parameters passed to the base constructor:

* `rank` -- tensor rank \f$ k \f$ (integer).
* `parity` -- `Parity::even` or `Parity::odd`.
* `constant` -- overall multiplicative factor \f$ c \f$.
* `vec` -- radial function \f$ v(r) \f$ as `std::vector<double>` (may be empty, equivalent to 1).
* `RorI` -- `Realness::real` or `Realness::imaginary`.
* `freq_dep` -- set `true` if the operator depends on frequency or momentum transfer.

### Standard case: override angular coefficients

For operators whose radial integral fits the default form above:

* Implement `angularF()` (required)
* optionally override `angularCff()`, `angularCgg()`, `angularCfg()`, `angularCgf()`
  * These return \f$ C_{ff}, C_{gg}, C_{fg}, C_{gf} \f$ respectively; defaults are 1,1,0,0.
  * These are often constant, but may depend on angular quantum numbers too
  * Usually, but not always, operators are either diagonal or off-diagonal in spinor space
    * e.g., diagonal: \f$ C_{ff}, C_{gg} \neq 0 \f$ and \f$ C_{fg} = C_{gf} = 0 \f$
  
The radial integral is then handled automatically.

### Non-standard case: override the radial functions

If the radial integral cannot be expressed in the default form
(e.g., it involves derivatives, non-local terms, or mixing not captured by the \f$ C_{xy} \f$ coefficients),
override **both** `radial_rhs()` and `radialIntegral()` consistently:

* `radial_rhs(ka, Fb)` returns \f$ \delta F_b \f$ such that \f$ F_a \cdot \delta F_b = R_{ab} \f$.
* `radialIntegral(Fa, Fb)` returns \f$ R_{ab} \f$ directly.

* Both must be overridden together, or results will be inconsistent.
* `angularCff()`, `angularCgg()`, `angularCfg()`, `angularCgf()` are not used (unless you explicitely use them in your overridden radial functions)

### Frequency-dependent operators

For frequency-dependent operators, pass `freq_dep=true` to the base constructor and override `updateFrequency(double omega)`.
This is called by modules such as `MatrixElements` whenever the frequency changes
(e.g., when `omega = each` is set in the input). The base implementation aborts at runtime
if not overridden.

---

## Registering your operator

To expose the operator at runtime:

1. Add a `generate_MyOp()` factory function and register it in the `operator_list` in
   `src/DiracOperator/GenerateOperator.hpp`.
2. Include your header in `src/DiracOperator/Operators/include.hpp`.

---

## Existing implementations as examples

Browse `src/DiracOperator/Operators/` for a wide range of examples:

* Electric dipole, \f$k\f$-pole -- see `Ek.hpp`.
* Hyperfine structure with non-trivial angular structure -- see `hfs.hpp`.
* Frequency-dependent operators with spherical Bessel functions -- see the multipole operators.
