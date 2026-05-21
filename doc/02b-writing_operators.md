\page modules_custom_operator Writing a custom operator

\brief How to write your own DiracOperator for use in ampsci modules

Operators in ampsci are single-particle tensor operators of rank \f$ k \f$, acting on
\ref DiracSpinor objects. They are represented by the \ref DiracOperator::TensorOperator
virtual base class. All built-in operators (E1, hfs, pnc, ...) derive from this class,
and you can add your own the same way.

* See \ref modules_custom for how to write a module that uses your operator.
* See \ref modules for the list of existing modules and operators.
* Use `./ampsci -o` command-line option so see a list of available operators
  * `./ampsci -o <OperatorName>` for list of input options for specific operator
  * This will also work for your custom operator once it's compiled

@note If you plan to contribute your operator directly to the ampsci codebase, see the [Contributing guidelines](contributing.html) for how to structure and submit your contribution.

---

## Why write an operator?

If you need a matrix element that is not already implemented, the operator framework gives you:

* Automatic reduced matrix elements via the Wigner-Eckart theorem.
* Selection rules enforced by the base class (rank, parity).
* Compatibility with all existing modules (e.g., `Module::MatrixElements`, RPA, structure radiation).
* Frequency dependence handled cleanly if needed.
* RPA (see \ref ExternalField) and any other many-body effects (Structure Radiation etc.) will all "just work" if you implement the operator in this way

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

Done in three steps (details below):

1. **Write the operator class** -- derive from `TensorOperator`, implement `angularF()` and any needed overrides.
2. **Add a static `generate()` factory** -- required for runtime lookup by name; without it the operator cannot be used from input files or by any module.
3. **Register it** -- add a `Register<T>` entry in a `.cpp` file so ampsci discovers it at startup.

By default, put everything in a single new `.cpp` file (class, `generate()`, and registration), dropped anywhere in the source tree -- no header needed.
Split into `.hpp`/`.cpp` only if the class needs to be visible to other files (i.e., if we plan to add this to the core ampsci operator list).
The `.cpp` file may also live outside the source tree (see "Built-in vs. external operators" below).

### Write the operator

First, derive from \ref DiracOperator::TensorOperator and pass the operator properties to the base constructor, for example:

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

  // Optional (angularCff, angularCfg, angularCgf, angularCgg)
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

You can add any other public/private functions or variables if you require.

@note
It's very important to get the `rank`, `parity`, and `Realness` correct, since they dictate the selection rules and symmetry properties.

#### Standard case: override angular coefficients

For operators whose radial integral fits the default form above:

* Implement `angularF()` (required)
* optionally override `angularCff()`, `angularCgg()`, `angularCfg()`, `angularCgf()`
  * These return \f$ C_{ff}, C_{gg}, C_{fg}, C_{gf} \f$ respectively; defaults are 1,1,0,0.
  * These are often constant, but may depend on angular quantum numbers too
  * Usually, but not always, operators are either diagonal or off-diagonal in spinor space
    * e.g., diagonal: \f$ C_{ff}, C_{gg} \neq 0 \f$ and \f$ C_{fg} = C_{gf} = 0 \f$
  
The radial integral is then handled automatically.

#### Non-standard case: override the radial functions

If the radial integral cannot be expressed in the default form
(e.g., it involves derivatives, non-local terms, or mixing not captured by the \f$ C_{xy} \f$ coefficients),
override **both** `radial_rhs()` and `radialIntegral()` consistently:

* `radial_rhs(ka, Fb)` returns \f$ \delta F_b \f$ such that \f$ F_a \cdot \delta F_b = R_{ab} \f$.
* `radialIntegral(Fa, Fb)` returns \f$ R_{ab} \f$ directly.

* Both must be overridden together, or results will be inconsistent.
* `angularCff()`, `angularCgg()`, `angularCfg()`, `angularCgf()` are not used (unless you explicitely use them in your overridden radial functions)

#### Frequency-dependent operators

For frequency-dependent operators, pass `freq_dep=true` to the base constructor and override `updateFrequency(double omega)`.
This is called by modules such as `MatrixElements` whenever the frequency changes
(e.g., when `omega = each` is set in the input). The base implementation aborts at runtime if this is every required but not overridden to prevent silent failures.

---

### Adding the generate "factory" function: `generate()`

Add this static function to your class so ampsci can find the operator by name (from input files, `ampsci -o`, and modules like `MatrixElements`).
The signature and return value must be exactly like:

```cpp
static std::unique_ptr<DiracOperator::TensorOperator>
generate(const IO::InputBlock &input, const Wavefunction &wf){
  // details
  return std::make_unique<MyOperator>(/*any options your operator takes*/);
}
```

The implementation details are up to you; for most simple operators, they will be very simple (or even blank).
You don't have to use `input` or `wf`, but they are available for reading user options and accessing the wavefunction/nucleaus/grid etc.
For example, an operator that depends on the radial grid, the nuclear charge, as well as two parameters, \f$ c \f$ and \f$ n \f$, may be written like:

```cpp
static std::unique_ptr<TensorOperator>
generate(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"c", "Overall scale factor [default: 1.0]"},
               {"n",    "Power of r [default: 2]"}});
  if (input.has_option("help"))
    return nullptr;

  const auto c = input.get("c", 1.0);
  const auto n = input.get("n", 2);
  return std::make_unique<MyOperator>(wf.grid(), wf.Znuc(), c, n);
}
```

* `input.check(...)` is optional but strongly recommended -- it is what populates `ampsci -o` and catches typos in user input. Return `nullptr` for the `help` case.
* Same for `if (input.has_option("help")) return nullptr;`
* These are what makes the `ampsci -o` command-line helper works, and how users will know what options your operator takes

---

### Registering your operator

Add a `Register<T>` entry, where `T` is your operator class name, into the same .cpp file as your class (inside an anonymous namespace) so ampsci discovers it at startup:

```cpp
namespace {
const DiracOperator::Register<DiracOperator::MyOperator> r_MyOperator{
  "MyOperator", "Short description (shown by ampsci -o)"};
} // namespace
```

(If your operator is written into a .hpp file instead, add this to a seperate .cpp file -- see "Built-in vs. external operators" below).
The first argument is the name used in input files (e.g., `MatrixElements::MyOp{}`); the second is shown by `ampsci -o`.

---

## Built-in vs. external operators

**External (recommended):** put everything in a single `.cpp` file that lives _outside_ the main ampsci source code, and add it to `EXTERNAL_OPERATORS` in your Makefile -- no changes to the ampsci source required:

```makefile
EXTERNAL_OPERATORS = path/to/MyOperator.cpp
```

* This allows you to write and share operators without touching the ampsci source, so ampsci can be updated without conflicts with your code. Can be version-controlled and shared (e.g. via GitHub) independently.
* Downside: the operator can only be accessed through `generate()` -- you cannot use it directly as a concrete C++ type. This is sufficient for the vast majority of use cases.

**Built-in** (contributing to the ampsci repo): define the class in `src/DiracOperator/Operators/MyOperator.hpp`, add a `Register<MyOperator>` line to `src/DiracOperator/RegisterOperators.cpp`, and add the header to `src/DiracOperator/Operators/include.hpp`.

* The full class type is visible to the rest of the codebase, so other code can construct or cast to it directly.

---

## Existing implementations as examples

Browse `src/DiracOperator/Operators/` for a wide range of examples:

* Electric dipole, \f$k\f$-pole -- see `Ek.hpp`.
* Hyperfine structure with non-trivial angular structure -- see `hfs.hpp`.
* Frequency-dependent operators with spherical Bessel functions -- see the multipole operators.

---

## Key API reference

The following links are to the detailed documentation for the specific namespaces/classes that will be most useful for writing your own operator:

* \ref DiracOperator::TensorOperator -- single-particle tensor operators

* \ref DiracSpinor -- single relativistic orbital \f$ F_{n\kappa} = (f, g) \f$;
  radial components, quantum numbers, arithmetic, inner products.

* \ref Wavefunction -- full atomic state: core/valence/basis orbital lists,
  radial grid, nuclear potential, and the HF object. Access valence orbitals
  via `wf.valence()`, core via `wf.core()`, basis via `wf.basis()`.

* \ref Angular -- Angular coeficients, 3j/6j symbols etc.

* \ref Nuclear -- Nuclear data and potentials

* \ref SphericalBessel -- Spherical Bessel functions

* \ref PhysConst -- Physical constants, unit conversions

* \ref Grid -- numerical grids, including Jacobian

* [Namespaces](namespaces.html) -- full ampsci API docs
