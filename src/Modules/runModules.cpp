#include "Modules/runModules.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/modules_list.hpp"
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
    }
  }
}

//==============================================================================
void runModule(const IO::InputBlock &module_input, const Wavefunction &wf) {

  const auto &in_name = module_input.name();
  std::cout << '\n';
  IO::print_line('-');
  std::cout << "Module: " << in_name << "\n";

  // Loop through all available modules, run correct one
  for (const auto &[mod_name, mod_func, description] : module_list) {
    if (qip::ci_wc_compare("Module::" + mod_name, in_name) ||
        qip::ci_wc_compare(mod_name, in_name))
      return mod_func(module_input, wf);
  }

  // Only reach here if module not found:
  std::cout << "\nThere is no available module named: " << in_name << "\n";

  // spell-check + nearest suggestion:
  const auto compare_sc = [&in_name](const auto &s1, const auto &s2) {
    return qip::ci_Levenstein(s1.name, in_name) <
           qip::ci_Levenstein(s2.name, in_name);
  };
  const auto guess =
      std::min_element(module_list.cbegin(), module_list.cend(), compare_sc);
  if (guess != module_list.cend()) {
    std::cout << "\nDid you mean: " << guess->name << " ?\n";
  }

  std::cout << "\nAvailable modules:\n";
  list_modules();
  std::cout << "\n";
}

//==============================================================================
void list_modules() {
  for (auto &[name, func, description] : module_list) {
    // fmt::print(" * {:23s} : {}\n", name, description);
    fmt::print(" * {}\n", name);
    if (!description.empty())
      fmt2::styled_print(fg(fmt::color::light_blue), "{}\n",
                         qip::wrap(description, 80, "     "));
  }
}

} // namespace Module

//transitionPolarisability