#include "Modules/Modules.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <iostream>

namespace Module {

//==============================================================================
void runModules(const IO::InputBlock &input, const Wavefunction &wf) {

  // Check if block is a module (modules state with 'Module::')
  const auto block_is_module = [](auto &block) {
    return qip::ci_wc_compare(block.name(), "Module*");
  };
  // if it is, run it with its input
  for (const auto &block : input.blocks()) {
    if (block_is_module(block)) {
      runModule(block, wf);
      std::cout << std::flush;
    }
  }
}

//==============================================================================
void runModule(const IO::InputBlock &module_input, const Wavefunction &wf) {

  const auto &in_name = module_input.name();
  std::cout << '\n';
  IO::print_line('-');
  std::cout << in_name << "\n" << std::flush;

  const auto &entries = Registry::get().entries();

  // Loop through all registered modules, run correct one
  for (const auto &entry : entries) {
    if (qip::ci_wc_compare("Module::" + entry.name, in_name) ||
        qip::ci_wc_compare(entry.name, in_name))
      return entry.function(module_input, wf);
  }

  // Only reach here if module not found:
  std::cout << "\nCould not find requested module.\n";

  std::vector<std::string> names;
  for (const auto &e : entries) {
    names.push_back(e.name);
  }
  IO::unkown_option(in_name, names);

  std::cout << "\nTo see a list of available modules:  `./ampsci -m`\n";
}

//==============================================================================
void list_modules() {
  for (const auto &entry : Registry::get().entries()) {
    fmt::print(" * {}\n", entry.name);
    if (!entry.description.empty())
      fmt2::styled_print(fg(fmt::color::light_blue), "{}\n",
                         qip::wrap(entry.description, 80, "     "));
  }
}

} // namespace Module
