#pragma once
#include "IO/InputBlock.hpp"
#include "IO/InputBlockLegacy.hpp"
#include "json/json.hpp"
#include <string>

namespace IO {

namespace detail {

// Forward declaration for recursion
inline nlohmann::json ib_legacy_to_json(const IO::InputBlockLegacy &ib);

inline nlohmann::json ib_legacy_to_json(const IO::InputBlockLegacy &ib) {
  auto j = nlohmann::json::object();
  auto modules = nlohmann::json::array();

  for (const auto &opt : ib.options()) {
    try {
      j[opt.key] = nlohmann::json::parse(opt.value_str);
    } catch (...) {
      j[opt.key] = opt.value_str;
    }
  }

  for (const auto &block : ib.blocks()) {
    const auto name = std::string{block.name()};
    // Collect Module::* blocks into a "Module" JSON array
    if (name.size() > 8 && name.substr(0, 8) == "Module::") {
      auto entry = ib_legacy_to_json(block);
      entry["type"] = name.substr(8);
      modules.push_back(std::move(entry));
    } else {
      j[name] = ib_legacy_to_json(block);
    }
  }

  if (!modules.empty())
    j["Module"] = std::move(modules);

  return j;
}

} // namespace detail

/*!
  @brief Convert a legacy IO::InputBlockLegacy to IO::InputBlock.
  @details
  Recursively walks the legacy tree of options and sub-blocks and builds an
  equivalent JSON node. Values are parsed as JSON where possible (numbers,
  booleans); anything that does not parse is stored as a string.

  Blocks whose names begin with "Module::" are collected into a "Module" JSON
  array, each with a "type" key set to the module name. This matches the
  format expected by Module::runModules().
*/
inline InputBlock from_legacy(const IO::InputBlockLegacy &ib) {
  return InputBlock{std::string{ib.name()}, detail::ib_legacy_to_json(ib)};
}

} // namespace IO
