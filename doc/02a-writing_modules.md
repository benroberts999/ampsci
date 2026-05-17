\page modules_custom Writing custom modules

\brief Instructions for writing your own custom ampsci modules

A module is the standard way to do calculations with ampsci wavefunctions --
if none of the built-in modules do what you need, writing your own is
straightforward. Modules are just C++ functions: you write the physics,
ampsci handles the wavefunction.

* Any number of _modules_ can be run by adding `Module::moduleName{}` blocks
  to the input file.
* Get a list of available modules: `./ampsci -m`
* See \ref modules for details of currently available modules.
* See \ref tutorial_modules for a hands-on introduction to using modules.
* You can also write your own operators -- see \ref modules_custom_operator.

## Writing your module

1. **Create the file.** 
 - Copy `src/Modules/exampleModule.cpp` to a new file (e.g. `src/Modules/myModule.cpp`) and rename `exampleModule` to your module name throughout.
 - You can name the module whatever you like (`myModule` in this example), but it must be in the `Module` namespace, and the function signature must be exactly:
  ```cpp
    void myModule(const IO::InputBlock &, const Wavefunction &)
  ```
 - Place the file anywhere under `src/` (though use `src/Modules/` by default); the build system will pick it up automatically.
2. **Write the physics.** 
 - Replace the body of your function with whatever calculation you need. 
 - Inside it, `input` (\ref IO::InputBlock) gives you access to the run-time user options, and `wf` (\ref Wavefunction) gives you the solved atomic wavefunction.
 - See detailed API documentation and the provided examples for details on how to use the \ref Wavefunction
3. **Register the module.** 
 - Add a single `Register` line inside an anonymous namespace (see below for example)
 - This must be in the same .cpp file as your module (_after_ the module has been declared).
 - This is what makes the module visible to ampsci, and allows users to run it.
4. **Recompile and run.** 
 - Add `Module::myModule{}` to an input file, or query its options with `./ampsci -m myModule`. 
 - You must recompile after adding a new module for it to be visible to ampsci.

No other files need to be edited. The module self-registers at program startup.

Each module file declares its function, registers it with the @ref Module
namespace, then defines it -- all in a single `.cpp` file. 
The structure of a module will be like the following example:

```cpp
namespace Module {

// Implement your module:
void myModule(const IO::InputBlock &input, const Wavefunction &wf) {

  input.check({{"value1", "Short description of value1 [3.14]"},
             {"value2", "Short description of value2 [default2]"}});
  
  // parse your input:
  const double value1_default = 3.14;
  const double value1 = input.get("value1", value1_default);

  // Use the wavefunction to do something:
  for (const auto &v : wf.valence()){
    std::cout << v.shortSymbol() << " " << v.en() <<"\n";
  }
}

// Register your module (in an anonymous namespace) so it's visible to ampsci:
namespace {
const Register r_myModule{"myModule",
                          "Short description of myModule",
                          &myModule};
} // namespace

} // namespace Module
```

* **name** (`"myModule"`): the string used to invoke the module from an input file (`Module::myModule{...}`), and as the lookup key for `./ampsci -m myModule`. Conventionally the same as the function name.
* **description**: a short, one-line summary shown by `./ampsci -m`.
* **function** (`&myModule`): pointer to the module function. The function's signature must match exactly `void name(const IO::InputBlock&, const Wavefunction&)` -- otherwise the address won't convert to the @ref Module::ModuleFn function pointer and you'll get a compile error.
* The `r_myModule` name for the @ref Register is just convention -- it can be anything.
  * In fact, the only restriction is that multiple modules within the same .cpp file (compilation unit) have unique names; the names do not even need to be unique across different files. (These names are not used anywhere in the code)
* The registration must appear _after_ the module function has been declared (a forward declaration is enough; the definition can come later).
* Keep the registration in the _same_ `.cpp` file as the module: this isn't a hard language requirement, but it's the right convention -- splitting them buys nothing and makes the module harder to find.

## Highly recommended: input checking

Add an `input.check()` call listing all options used by your module:

```cpp
input.check({{"option1", "Short description of option1 [default1]"},
             {"option2", "Short description of option2 [default2]"}});
```

* Catches spelling mistakes in user input -- a mistyped option would be silently
  ignored otherwise.
* The descriptions are printed when the user runs `./ampsci -m myModule`.
* Without this it is much harder for other users to know how to use your module, and it provides (basic) automatic input checking.

Return immediately after `check()` if help was requested, to avoid running
the module unnecessarily (this is what makes `ampsci -m moduleName` work):

```cpp
if (input.has_option("help"))
  return;
```


## How registration works

Most users don't need to care -- copy the template and it just works.

The `Register` variable's constructor runs once, _before_ `main()`, as part of 
C++ static initialisation. 
The constructor body simply appends a new `ModuleEntry{name, description, &fn}` 
to the singleton `Module::Registry`.
By the time `main()` starts, every compiled-in module has self-registered.
The ampsci driver then iterates over the input file, looks up each 
`Module::name{...}` block in the Registry by name, and dispatches to 
the matching function pointer.

The anonymous namespace surrounding the Register variable gives it 
_internal linkage_; one private copy per `.cpp` file. 
That's why different module files can each use the same identifier 
(`r_myModule`, or whatever) without producing "multiple definition" link errors. 
It's also why the registration must live in a `.cpp`: 
putting it in a header would create one entry in any translation unit that includes the header.

For the full API (Register, Registry, ModuleEntry, runModule, list_modules),
see the @ref Module namespace documentation in
[Modules.hpp](\ref Modules.hpp).

## Key API reference

The following classes form the core API available inside a module.
See the linked API docs for full details.

* \ref DiracSpinor -- single relativistic orbital \f$ F_{n\kappa} = (f, g) \f$;
  radial components, quantum numbers, arithmetic, inner products.

* \ref Wavefunction -- full atomic state: core/valence/basis orbital lists,
  radial grid, nuclear potential, and the HF object. Access valence orbitals
  via `wf.valence()`, core via `wf.core()`, basis via `wf.basis()`.

* \ref HF::HartreeFock -- self-consistent HF field; accessible via `wf.vHF()`.

* \ref DiracOperator::TensorOperator -- virtual base class for single-particle
  tensor operators (E1, hfs, pnc, ...); reduced matrix elements, selection
  rules, radial integrals. See \ref modules_custom_operator for writing your
  own operator.

* \ref IO::InputBlock -- the input block passed to your module; use
  `input.get<T>("key", default_value)` to read user-supplied options, and
  `input.check({...})` to declare and validate them.

* \ref Module -- the namespace containing the registration mechanism
  (@ref Module::Register, @ref Module::Registry); see [Modules.hpp](\ref Modules.hpp)
  for the registration pattern.
