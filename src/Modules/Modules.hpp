#pragma once
#include <string>
#include <vector>

namespace IO {
class InputBlockLegacy;
class InputBlock;
} // namespace IO
class Wavefunction;

/*!
  @brief Modules are user-defined calculations run after the wavefunction has been solved.

  @details
  A _module_ is a self-contained C++ function that takes an already-computed
  atomic wavefunction (Hartree-Fock, MBPT, etc.) and performs some downstream
  calculation: matrix elements, polarisabilities, hyperfine constants, PNC
  amplitudes, lifetimes, and so on. Any number of modules can be run from a
  single input file by adding `Module::ModuleName{ ... }` blocks; the ampsci
  driver looks up each name in the registry below and dispatches to the
  corresponding function.

  See \ref modules_custom for a more detailed tutorial.

  ## Adding a new module

  Adding a module requires editing _exactly one_ new `.cpp` file. No central
  list needs to be updated; the build system picks up new `.cpp` files
  automatically, and the module self-registers at program startup via a static
  initialiser.

  The minimal template:

  \code{.cpp}
  // src/Modules/myModule.cpp
  #include "IO/InputBlock.hpp"
  #include "Modules/Modules.hpp"
  #include "Wavefunction/Wavefunction.hpp"

  namespace Module {

  // Declare, register, then define below.
  void myModule(const IO::InputBlock &input, const Wavefunction &wf);
  namespace {
  const Register r_myModule{"myModule", "One-line description shown in -m",
                            &myModule};
  } // namespace

  void myModule(const IO::InputBlock &input, const Wavefunction &wf) {

    input.check({{"option1", "Description of option1 [default1]"},
                 {"option2", "Description of option2 [default2]"}});
    if (input.has_option("help"))
      return;

    const auto option1 = input.get("option1", default1);
    // ... physics ...
  }

  } // namespace Module
  \endcode

  The first string passed to Register is the name used in input files:
  `Module::myModule{ ... }`. The Register variable itself is never referenced
  by name -- it exists purely so its constructor fires before `main()` and adds
  the entry to the singleton Registry. The anonymous namespace gives it
  internal linkage, so different module files can each use the same identifier
  (`r_myModule`, or whatever you like) without colliding at link time.

  A complete worked example is provided in `src/Modules/exampleModule.cpp`,
  intended as a template for new users to copy.

  ## Running a module

  In the input file:

  \code{.unparsed}
  Module::myModule{
    option1 = 3.14;
    option2 = true;
  }
  \endcode

  From the command line, to list available modules or query a module's
  options:

  \code{.unparsed}
  ./ampsci -m                 # list all modules
  ./ampsci -m myModule        # show options for `myModule`
  \endcode
*/
namespace Module {

//==============================================================================

//! Function-pointer signature shared by every module.
using ModuleFn = void (*)(const IO::InputBlock &, const Wavefunction &);

/*!
  @brief One entry in the module registry.
  @details
  Holds everything required to identify and invoke a module:

  - @c name: the string used in input files (`Module::<name>{ ... }`) and
    when querying with `./ampsci -m <name>`.
  - @c description: one-line summary printed by `./ampsci -m`.
  - @c function: pointer to the module's free function (signature
    @ref ModuleFn).
*/
struct ModuleEntry {
  std::string name;
  std::string description;
  ModuleFn function;
};

//==============================================================================

/*!
  @brief Singleton registry of all compiled-in modules.

  @details
  Populated at program startup via static initialisers -- one per module
  `.cpp` file -- and never modified once `main()` begins.

  The "construct on first use" idiom in get() avoids the static-initialisation
  order fiasco: the Register constructor calls Registry::get(), which forces
  the singleton to be created before any module tries to register itself.
  Order of registration across translation units is unspecified, but for this
  registry that does not matter.

  Module authors do not normally need to interact with the Registry directly;
  they self-register via @ref Register. The ampsci driver uses the Registry
  via @ref runModule, @ref runModules and @ref list_modules.
*/
class Registry {
public:
  /*!
    @brief Access the singleton instance.
    @return Reference to the (only) Registry.
  */
  static Registry &get() {
    static Registry instance;
    return instance;
  }

  /*!
    @brief Append a new entry to the registry.
    @details
    Normally called only by the @ref Register constructor.

    @param name        Module name as used in input files.
    @param description One-line description (shown by `./ampsci -m`).
    @param fn          Pointer to the module function.
  */
  void add(std::string name, std::string description, ModuleFn fn) {
    m_entries.push_back({std::move(name), std::move(description), fn});
  }

  /*!
    @brief All registered modules, in registration order.
    @return Const reference to the underlying vector of @ref ModuleEntry.
  */
  const std::vector<ModuleEntry> &entries() const { return m_entries; }

  /*!
    @brief Look up a module by name.
    @param name  Module name (case sensitive; same string passed to @ref Register).
    @return Pointer to the matching @ref ModuleEntry, or @c nullptr if no
            module with that name is registered.
  */
  const ModuleEntry *find(const std::string &name) const {
    for (const auto &e : m_entries) {
      if (e.name == name)
        return &e;
    }
    return nullptr;
  }

private:
  Registry() : m_entries{} {}
  std::vector<ModuleEntry> m_entries;
};

//==============================================================================

/*!
  @brief Helper struct: constructing a Register adds a module to the Registry.
  @details
  Used at file scope (inside an anonymous namespace) in every module `.cpp`
  file to self-register that module with the @ref Registry. The constructor
  is the only thing of interest -- it runs once, before `main()`, as part of
  static initialisation.

  See the namespace-level documentation above for the usage pattern.

  @note The variable itself is never referenced by name. Pick any identifier
        you like (`registrar`, `r_myModule`, ...); the anonymous namespace
        gives it internal linkage so different files can reuse the same
        identifier without clashing.
*/
struct Register {

  /*!
    @brief Register a module by constructing one of these at file scope.
    @param name        Module name as used in input files: `Module::<name>{}`.
    @param description One-line description, shown by `./ampsci -m`.
    @param fn          Pointer to the module function.
  */
  Register(const char *name, const char *description, ModuleFn fn) {
    Registry::get().add(name, description, fn);
  }
};

//==============================================================================

/*!
  @brief Iterate over the input blocks and run any that are modules.
  @details
  Scans @p input for blocks whose name matches `Module::*` and dispatches
  each to @ref runModule in turn. Called once by the ampsci driver after the
  wavefunction has been solved.

  @param input  Top-level input block.
  @param wf     Solved wavefunction passed to each module.
*/
void runModules(const IO::InputBlockLegacy &input, const Wavefunction &wf);

/*!
  @brief Run a single module.
  @details
  Looks up the module by name in the @ref Registry and calls its function.
  If the name is not found, prints an error message together with a
  nearest-match suggestion and a list of available modules.

  @param input  Input block for this module (typically `Module::Name{}`).
  @param wf     Solved wavefunction.
*/
void runModule(const IO::InputBlock &input, const Wavefunction &wf);

/*!
  @brief Run a single module (legacy InputBlockLegacy overload).
  @details
  Converts @p input to IO::InputBlock and delegates to the InputBlock
  overload. Used by the legacy .in file dispatch path.
*/
void runModule(const IO::InputBlockLegacy &input, const Wavefunction &wf);

/*!
  @brief Iterate over a JSON-format input block and run any modules found.
  @details
  Reads the "Module" array from @p input (new JSON format) and dispatches
  each entry directly to the registered module function as IO::InputBlock.
  No legacy string conversion is performed.

  @param input  Top-level JSON input block (must contain a "Module" array).
  @param wf     Solved wavefunction passed to each module.
*/
void runModules2(const IO::InputBlock &input, const Wavefunction &wf);

/*!
  @brief Print the list of compiled-in modules (name + description).
  @details
  Invoked by `./ampsci -m`. Iterates over Registry::entries() and prints each
  module's name and description (word-wrapped).
*/
void list_modules();

} // namespace Module
