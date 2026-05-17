#pragma once
#include <string>
#include <vector>

namespace IO {
class InputBlock;
}
class Wavefunction;

//! @brief Modules are run using calculated atomic wavefunctions
/*! @details
After the program has generated wavefunctions (Hartree-Fock etc), any number of
modules can then be run (see input documentation for how to run them).
Examples include: calculating matrix elements, lifetimes, PNC etc.

To add a new module:
  - Drop a single .cpp file into `src/Modules/` (or any subdirectory of `src/`
    that is part of the build).
  - In that file, define the module function in `namespace Module`, then
    register it with a `Module::Registrar` in an anonymous namespace.
  - No other files need to be edited; the build system picks up new `.cpp`
    files automatically, and the module self-registers at program startup.
  - See `Modules/exampleModule.cpp` for a template.
*/
namespace Module {

//! Function signature for all modules.
using ModuleFn = void (*)(const IO::InputBlock &, const Wavefunction &);

struct ModuleEntry {
  std::string name;
  std::string description;
  ModuleFn function;
};

/*!
  @brief Singleton registry of all compiled-in modules.
  @details
  Populated at program startup via static initialisers (one per module .cpp
  file). The "construct on first use" idiom in get() avoids the static
  initialisation order problem -- the Registrar constructor below calls
  Registry::get(), forcing the singleton to be created before any module
  tries to register itself.
*/
class Registry {
public:
  static Registry &get() {
    static Registry instance;
    return instance;
  }

  void add(std::string name, std::string description, ModuleFn fn) {
    m_entries.push_back({std::move(name), std::move(description), fn});
  }

  const std::vector<ModuleEntry> &entries() const { return m_entries; }

  //! Returns pointer to entry, or nullptr if not found
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

//! Static-lifetime helper: constructing one of these registers a module.
struct Registrar {
  Registrar(const char *name, const char *description, ModuleFn fn) {
    Registry::get().add(name, description, fn);
  }
};

//! Loops through all given modules, runs them one at a time
void runModules(const IO::InputBlock &input, const Wavefunction &wf);

//! Figures out which module to run, looking it up in the Registry
void runModule(const IO::InputBlock &input, const Wavefunction &wf);

//! Lists all available modules
void list_modules();

} // namespace Module
