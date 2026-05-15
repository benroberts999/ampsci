\page modules_custom Writing custom modules

\brief Instructions for writing your own custom ampsci modules

A module is the standard way to do calculations with ampsci wavefunctions --
if none of the built-in modules do what you need, writing your own is straightforward.
Modules are just C++ functions: you write the physics, ampsci handles the wavefunction.

* The modules system allows the easy calculation of any atomic properties after the wavefunction has been calculated.
* Any number of _modules_ can be run by adding a `Module::moduleName{}` block to the input file.
* Get a list of available modules: `./ampsci -m`
* See \ref modules for details of currently available modules.
* See \ref tutorial_modules for a hands-on introduction to using modules.
* The code is designed so that you can easily create your own modules.
* You can also write your own operators -- see \ref modules_custom_operator.

## Creating your own module

An example module is provided to help you write your own module:

* `src/Modules/exampleModule.hpp`
* `src/Modules/exampleModule.cpp`

Duplicate both files and give them a new name -- that is much easier than starting from scratch.

* Modules are functions with the following signature:

  ```cpp
  namespace Module {
  void exampleModule(const IO::InputBlock &input, const Wavefunction &wf);
  }
  ```

* `input` is an \ref IO::InputBlock that holds any input options supplied by the user.

* `wf` is the \ref Wavefunction object calculated by ampsci -- it gives access to orbitals, the grid, and the HF potential.

* Modules are typically placed in the `Module` namespace, but this is not required.

* Modules typically live in `src/Modules/`, but can live anywhere.

### Including your module into ampsci

* Update `src/Modules/module_list.hpp`:
  * Add the corresponding `#include` at the top.
  * Add a `std::pair` to the `module_list` vector:

    ```cpp
    {"moduleName", &Module::moduleName}
    ```

  * The string is the name used in the input file; the second element is a pointer to the function.

* Recompile ampsci.

* Run the module by adding a `Module::moduleName{}` block to the input file.

### Highly recommended: input checking

* Add an `input.check()` call for all options used in your module:

```cpp
  input.check({{"option1", "Short description of option1 [default1]"},
               {"option2", "Short description of option2 [default2]"}});
```

* Catches spelling mistakes in user input -- a mistyped option is silently ignored otherwise.
* Descriptions are printed when the user requests `help` for the module.

* Return immediately after `check()` if help was requested, to avoid running the module unnecessarily:

```cpp
  if (input.has_option("help")) {
    return;
  }
```

* Minimal example:

```cpp
void exampleModule(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"option1", "Short description of option1 [default1]"},
               {"option2", "Short description of option2 [default2]"}});

  if (input.has_option("help")) {
    return;
  }

  auto option1 = input.get("option1", default1);
  auto option2 = input.get("option2", default2);
  // option1 is set from user input, or default1 if not given
}
```

---

## Key API reference

The following classes form the core API available inside a module.
See the linked API docs for full details.

* \ref DiracSpinor -- single relativistic orbital \f$ F_{n\kappa} = (f, g) \f$; radial components, quantum numbers, arithmetic, inner products.

* \ref Wavefunction -- full atomic state: core/valence/basis orbital lists, radial grid, nuclear potential, and the HF object.
  Access valence orbitals via `wf.valence()`, core via `wf.core()`, basis via `wf.basis()`.

* \ref HF::HartreeFock -- self-consistent HF field; accessible via `wf.vHF()`.

* \ref DiracOperator::TensorOperator -- virtual base class for single-particle tensor operators (E1, hfs, pnc, ...); reduced matrix elements, selection rules, radial integrals.
  See \ref modules_custom_operator for writing your own operator.
