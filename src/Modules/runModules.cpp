#include "Modules/runModules.hpp"
#include "IO/InputBlock.hpp"
#include "Modules/modules_list.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Module {

//==============================================================================
void runModules(const IO::InputBlock &input, const Wavefunction &wf) {
  for (const auto &module : input.blocks()) {
    if (qip::ci_wildcard_compare(module.name(), "Module*"))
      runModule(module, wf);
  }
}

//==============================================================================
void runModule(const IO::InputBlock &module_input, const Wavefunction &wf) {

  const auto &in_name = module_input.name();
  std::cout << '\n';
  IO::print_line('-');
  std::cout << "Module: " << in_name << "\n";

  // Loop through all available modules, run correct one
  for (const auto &[mod_name, mod_func] : module_list) {
    if (qip::ci_wildcard_compare("Module::" + mod_name, in_name) ||
        qip::ci_wildcard_compare(mod_name, in_name))
      return mod_func(module_input, wf);
  }

  std::cout << "\nThere is no available module named: " << in_name << "\n";

  // spell-check + nearest suggestion:
  const auto compare_sc = [&in_name](const auto &s1, const auto &s2) {
    return qip::ci_Levenstein(s1.first, in_name) <
           qip::ci_Levenstein(s2.first, in_name);
  };
  const auto guess =
      std::min_element(module_list.cbegin(), module_list.cend(), compare_sc);
  if (guess != module_list.cend()) {
    std::cout << "\nDid you mean: " << guess->first << " ?\n";
  }

  std::cout << "\nAvailable modules:\n";
  for (const auto &module : module_list) {
    std::cout << "  " << module.first << "\n";
  }
  std::cout << "\n";
}

} // namespace Module
