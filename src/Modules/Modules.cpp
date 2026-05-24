#include "Modules/Modules.hpp"
#include "IO/InputBlock.hpp"
#include "IO/InputBlockLegacy.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "json/json.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <iostream>

namespace Module {

namespace {
// Convert a legacy IO::InputBlockLegacy to IO::InputBlock by walking the options
// and sub-blocks recursively. Values are parsed as JSON where possible
// (numbers, bools); anything that doesn't parse stays as a string.
IO::InputBlock ib_to_ib2(const IO::InputBlockLegacy &ib) {
  auto j = nlohmann::json::object();
  for (const auto &opt : ib.options()) {
    try {
      j[opt.key] = nlohmann::json::parse(opt.value_str);
    } catch (...) {
      j[opt.key] = opt.value_str;
    }
  }
  for (const auto &block : ib.blocks()) {
    j[std::string{block.name()}] = ib_to_ib2(block).node();
  }
  return IO::InputBlock{std::string{ib.name()}, std::move(j)};
}
} // namespace

//==============================================================================
void runModules(const IO::InputBlockLegacy &input, const Wavefunction &wf) {

  // Check if block is a module (modules state with 'Module::')
  const auto block_is_module = [](auto &block) {
    return qip::ci_wc_compare(block.name(), "Module*");
  };
  // if it is, run it with its input (converting to InputBlock)
  for (const auto &block : input.blocks()) {
    if (block_is_module(block)) {
      runModule(ib_to_ib2(block), wf);
      std::cout << std::flush;
    }
  }
}

//==============================================================================
void runModule(const IO::InputBlock &module_input, const Wavefunction &wf) {

  const auto in_name = std::string{module_input.name()};
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
void runModule(const IO::InputBlockLegacy &module_input,
               const Wavefunction &wf) {
  runModule(ib_to_ib2(module_input), wf);
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

//==============================================================================
void runModules2(const IO::InputBlock &input, const Wavefunction &wf) {
  const auto &node = input.node();
  const auto it = node.find("Module");
  if (it == node.end() || !it->is_array())
    return;

  for (const auto &entry : *it) {
    if (!entry.is_object())
      continue;
    const auto type = entry.value("type", std::string{});
    if (type.empty())
      continue;
    // Strip "type" key (module name metadata) and dispatch directly as
    // InputBlock -- no legacy string conversion.
    auto opts = entry;
    opts.erase("type");
    IO::InputBlock mod_ib2{"Module::" + type, std::move(opts)};
    runModule(mod_ib2, wf);
    std::cout << std::flush;
  }
}

} // namespace Module
