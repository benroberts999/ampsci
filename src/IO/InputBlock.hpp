#pragma once
#include "fmt/color.hpp"
#include "json/json.hpp"
#include "qip/String.hpp"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace IO {

//==============================================================================
namespace detail {

//! True if T is std::vector<U> for any U.
template <typename T>
struct is_vector : std::false_type {
  using element_type = T;
};
template <typename T>
struct is_vector<std::vector<T>> : std::true_type {
  using element_type = T;
};

//! True if T is std::array<U,N> for any U and N.
template <typename T>
struct is_std_array : std::false_type {
  using element_type = T;
  static constexpr std::size_t extent = 0;
};
template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type {
  using element_type = T;
  static constexpr std::size_t extent = N;
};

} // namespace detail

//==============================================================================
/*!
  @brief JSON-backed structured input block.
  @details
  Light wrapper around nlohmann::json that provides a similar interface to
  IO::InputBlockLegacy, but backed by a standard JSON object. Designed for reading
  ampsci input files in JSON format.

  Comments (//) and block comments (slash-star) are supported in input files
  via nlohmann's built-in ignore_comments parse option. The file extension
  .json or .jsonc is recommended.

  Top-level blocks are JSON objects. Modules are an array under the
  "Module" key; each entry must have a "type" field giving the module
  name. Nested sub-blocks are nested JSON objects.

  @code{.json}
  {
    "Atom": { "Z": "Cs" },
    "Grid": { "r0": 1e-6, "rmax": 120.0, "num_points": 1000 },
    "Module": [
      { "type": "matrixElements", "operator": "E1", "rpa": true },
      {
        "type": "matrixElements",
        "operator": "hfs",
        "options": {
          "F": "SingleParticle",
          "SingleParticle": { "mu": 3.0, "gs": 1.654 }
        }
      }
    ]
  }
  @endcode

  @note JSON keys are case-sensitive. Adopt lowercase-only keys to match
        existing ampsci conventions.
  @note get<bool>() accepts both native JSON booleans and the strings
        "true"/"yes"/"y" (case-insensitive) for compatibility.
  @note get<std::vector<T>>() accepts both JSON arrays and
        comma-separated strings.
*/
class InputBlock {
public:
  //! Default-construct an empty block with no name.
  InputBlock() = default;

  //! Construct from a name and an existing JSON object node.
  explicit InputBlock(std::string name,
                      nlohmann::json node = nlohmann::json::object())
    : m_name(std::move(name)), m_node(std::move(node)) {}

  /*!
    @brief Parse a JSON (or JSONC) input file from @p path.
    @details
    Both line comments (//) and block comments (slash-star) are stripped
    before parsing.
    Throws std::runtime_error on file-open or parse failure.
  */
  static InputBlock from_file(const std::string &path) {
    std::ifstream f(path);
    if (!f)
      throw std::runtime_error("InputBlock: cannot open file: " + path);
    try {
      // ignore_comments = true: strips // and /* */ before parsing
      auto j = nlohmann::json::parse(f, nullptr,
                                     /*allow_exceptions=*/true,
                                     /*ignore_comments=*/true);
      return InputBlock{path, std::move(j)};
    } catch (const nlohmann::json::parse_error &e) {
      throw std::runtime_error(std::string("InputBlock: JSON parse error in ") +
                               path + ": " + e.what());
    }
  }

  //! Returns the name of this block.
  std::string_view name() const { return m_name; }

  //! Returns a const reference to the underlying JSON node (escape hatch).
  const nlohmann::json &node() const { return m_node; }
  //! Returns a mutable reference to the underlying JSON node (escape hatch).
  nlohmann::json &node() { return m_node; }

  //==============================================================================
  // Value access

  /*!
    @brief Returns the value of @p key, or @p default_value if not found or
           wrong type.
    @details
    Supports T = bool, integral types, floating-point types, std::string,
    std::vector<T>, std::array<T,N>. For bool, also accepts "true"/"yes"/"y"
    strings. For vector/array, accepts either a JSON array or a
    comma-separated string.
    @note Cannot be used with T = const char*; use std::string.
  */
  template <typename T>
  T get(std::string_view key, T default_value) const {
    static_assert(!std::is_same_v<T, const char *>,
                  "Cannot use get with const char* -- use std::string");
    return get<T>(key).value_or(default_value);
  }

  /*!
    @brief Returns an optional value for @p key; empty if not found, null,
           or type mismatch.
  */
  template <typename T = std::string>
  std::optional<T> get(std::string_view key) const;

  //==============================================================================
  // Block access

  /*!
    @brief Returns an optional copy of the child block named @p name.
    @details
    Returns nullopt if the key is absent or its value is not a JSON object.
  */
  std::optional<InputBlock> getBlock(std::string_view name) const {
    const auto it = m_node.find(std::string(name));
    if (it == m_node.end() || !it->is_object())
      return std::nullopt;
    return InputBlock{std::string(name), *it};
  }

  //! Returns the child block named @p name, or an empty block if not found.
  InputBlock get_block(std::string_view name) const {
    auto tmp = getBlock(name);
    return tmp ? *tmp : InputBlock{std::string(name)};
  }

  //! Returns true if a child block (object) named @p name exists.
  bool has_block(std::string_view name) const {
    return getBlock(name).has_value();
  }

  //! Returns true if @p key is present in this block (regardless of value type).
  bool has_option(std::string_view key) const {
    return m_node.contains(std::string(key));
  }

  //! Returns true if @p key is present and its value is not null or empty string.
  bool option_is_set(std::string_view key) const {
    return get<std::string>(key).has_value();
  }

  //==============================================================================
  // Dynamic construction (replaces IO::InputBlockLegacy::add(Option) sticky tape)

  /*!
    @brief Sets (or overwrites) @p key to @p value.
    @details
    Equivalent to `m_node[key] = value`. Supports any type that nlohmann::json
    can serialise (bool, int, double, std::string, etc.).
  */
  template <typename T>
  void set(std::string_view key, T value) {
    m_node[std::string(key)] = std::move(value);
  }

  //! Sets (or overwrites) @p name as a child block.
  void set_block(std::string_view name, InputBlock block) {
    m_node[std::string(name)] = std::move(block.m_node);
  }

  //==============================================================================
  // Validation

  /*!
    @brief Validates the options and sub-blocks in this block against @p list.
    @details
    Checks each key in this block against the allowed list. Keys ending with
    `{}` in the list are treated as allowed sub-blocks; all others are
    treated as allowed options. Unknown keys trigger a warning with a
    nearest-spelling suggestion. If @p print is true, or if a `"help"` key
    is present, the full list of allowed options is printed.

    @param list  Pairs of {key, description}. Sub-blocks are marked with
                 a `{}` suffix on the key, e.g. `{"SubBlock{}", "desc"}`.
    @param print If true, always print the full option list.
    @return true if all keys were found in the list; false otherwise.
  */
  bool check(const std::vector<std::pair<std::string, std::string>> &list,
             bool print = false) const;

  //! Prints the block contents as indented JSON to @p os.
  void print(std::ostream &os = std::cout, int /*indent_depth*/ = 0) const {
    os << m_node.dump(2) << "\n";
  }

  /*!
    @brief Serialise to the legacy ampsci block format (key = value; Block{}).
    @details
    Converts this block's JSON content to the old-style ampsci input string.
    Used as a bridge so that JSON-parsed input can be passed to code that
    still accepts IO::InputBlockLegacy.

    Conversion rules:
    - Scalar values (string, number, bool) become  key = value;
    - Nested objects become  key { <recursive> }
    - A top-level "Module" array of objects expands to
      Module::type { <options> } blocks (one per array entry, skipping "type").
    - Other arrays are serialised as comma-separated lists (scalars only;
      arrays of objects are silently skipped).
  */
  std::string to_ampsci_string() const;

private:
  std::string m_name;
  nlohmann::json m_node = nlohmann::json::object();

  // Used by to_ampsci_string: serialise one key+value to old ampsci format.
  static std::string scalar_or_block(const std::string &key,
                                     const nlohmann::json &val);
};

//==============================================================================
// InputBlock::get<T> implementation
//==============================================================================

template <typename T>
std::optional<T> InputBlock::get(std::string_view key) const {
  const auto it = m_node.find(std::string(key));
  if (it == m_node.end() || it->is_null())
    return std::nullopt;

  const auto &val = *it;

  // --- std::vector<E> ---
  if constexpr (detail::is_vector<T>::value) {
    using E = typename detail::is_vector<T>::element_type;
    if (!val.is_array())
      return std::nullopt;
    std::vector<E> out;
    out.reserve(val.size());
    for (const auto &elem : val) {
      if (elem.is_null())
        continue;
      if constexpr (std::is_arithmetic_v<E>) {
        if (!elem.is_number())
          return std::nullopt;
        out.push_back(static_cast<E>(elem.template get<double>()));
      } else {
        out.push_back(elem.template get<E>());
      }
    }
    return out;
  }

  // --- std::array<E,N> ---
  else if constexpr (detail::is_std_array<T>::value) {
    using E = typename detail::is_std_array<T>::element_type;
    constexpr std::size_t N = detail::is_std_array<T>::extent;
    if (!val.is_array() || val.size() != N)
      return std::nullopt;
    T out;
    for (std::size_t i = 0; i < N; ++i) {
      if constexpr (std::is_arithmetic_v<E>) {
        if (!val[i].is_number())
          return std::nullopt;
        out[i] = static_cast<E>(val[i].template get<double>());
      } else {
        out[i] = val[i].template get<E>();
      }
    }
    return out;
  }

  // --- bool ---
  // Accept native JSON boolean or "true"/"yes"/"y" strings.
  else if constexpr (std::is_same_v<T, bool>) {
    if (val.is_boolean())
      return val.template get<bool>();
    if (val.is_string()) {
      const auto str = val.template get<std::string>();
      if (qip::ci_wc_compare("true", str) || qip::ci_wc_compare("yes", str) ||
          qip::ci_wc_compare("y", str))
        return true;
      if (qip::ci_wc_compare("false", str) || qip::ci_wc_compare("no", str) ||
          qip::ci_wc_compare("n", str))
        return false;
    }
    return std::nullopt;
  }

  // --- std::string ---
  else if constexpr (std::is_same_v<T, std::string>) {
    if (!val.is_string())
      return std::nullopt;
    const auto str = val.template get<std::string>();
    // treat "" or "default" as unset (matches InputBlockLegacy behaviour)
    if (str.empty() || qip::ci_wc_compare("default", str))
      return std::nullopt;
    return str;
  }

  // --- integral types ---
  else if constexpr (std::is_integral_v<T>) {
    if (val.is_number_integer())
      return static_cast<T>(val.template get<int64_t>());
    if (val.is_number_float())
      return static_cast<T>(val.template get<double>());
    return std::nullopt;
  }

  // --- floating-point types ---
  else if constexpr (std::is_floating_point_v<T>) {
    if (val.is_number())
      return static_cast<T>(val.template get<double>());
    return std::nullopt;
  }

  // --- fallback: let nlohmann try ---
  else {
    try {
      return val.template get<T>();
    } catch (...) {
      return std::nullopt;
    }
  }
}

//==============================================================================
// InputBlock::check implementation
//==============================================================================

inline bool
InputBlock::check(const std::vector<std::pair<std::string, std::string>> &list,
                  bool print) const {
  bool all_ok = true;

  // Separate allowed option names from allowed block names
  std::vector<std::string> allowed_opts;
  std::vector<std::string> allowed_blocks;
  for (const auto &[key, _] : list) {
    if (!key.empty() && key.back() == '}') {
      // strip "{}" suffix
      allowed_blocks.push_back(key.substr(0, key.size() - 2));
    } else if (!key.empty()) {
      allowed_opts.push_back(key);
    }
  }

  for (const auto &[key, val] : m_node.items()) {
    // "help" is always silently accepted; triggers full-print
    if (qip::ci_wc_compare("help", key)) {
      print = true;
      continue;
    }

    const bool is_obj = val.is_object();

    if (is_obj) {
      // Check against allowed block names
      const bool found =
        std::any_of(allowed_blocks.cbegin(), allowed_blocks.cend(),
                    [&](const auto &b) { return qip::ci_wc_compare(b, key); });
      if (!found) {
        all_ok = false;
        fmt2::styled_print(fg(fmt::color::orange), "\nWARNING\n");
        std::cout << "Unclear input block within " << m_name << ": " << key
                  << "{}\nBlock and containing options may be ignored!\n"
                  << "Check spelling (or update list of options)\n";
        auto cmp = [&key](const auto &s1, const auto &s2) {
          return qip::ci_Levenstein(s1.first, key) <
                 qip::ci_Levenstein(s2.first, key);
        };
        const auto guess = std::min_element(list.cbegin(), list.cend(), cmp);
        if (guess != list.cend())
          std::cout << "\nDid you mean: " << guess->first << " ?\n";
      }
    } else {
      // Check against allowed option names
      const bool found =
        std::any_of(allowed_opts.cbegin(), allowed_opts.cend(),
                    [&](const auto &o) { return qip::ci_wc_compare(o, key); });
      if (!found) {
        all_ok = false;
        fmt2::styled_print(fg(fmt::color::orange), "\nWARNING\n");
        std::cout << "Unclear input option in " << m_name << ": " << key
                  << " = " << val.dump() << "\n"
                  << "Option may be ignored!\n"
                  << "Check spelling (or update list of options)\n";
        auto cmp = [&key](const auto &s1, const auto &s2) {
          return qip::ci_Levenstein(s1.first, key) <
                 qip::ci_Levenstein(s2.first, key);
        };
        const auto guess = std::min_element(list.cbegin(), list.cend(), cmp);
        if (guess != list.cend())
          std::cout << "\nDid you mean: " << guess->first << " ?\n";
      }
    }
  }

  if (!all_ok || print) {
    fmt2::styled_print(fg(fmt::color::light_blue),
                       "\n// Available {} options/blocks\n", m_name);
    fmt2::styled_print(fmt::emphasis::bold, m_name);
    std::cout << " {\n";
    for (const auto &[opt, desc] : list) {
      if (opt.empty())
        continue;
      const bool is_block = (opt.back() == '}');
      std::cout << "  " << opt << (is_block ? "\n" : ";\n");
      if (!desc.empty()) {
        fmt2::styled_print(fg(fmt::color::light_blue), "{}\n",
                           qip::wrap(desc, 60, "    // "));
      } else {
        std::cout << "\n";
      }
    }
    std::cout << "}\n\n";
  }

  return all_ok;
}

//==============================================================================
// InputBlock::to_ampsci_string implementation
//==============================================================================

inline std::string InputBlock::to_ampsci_string() const {
  std::string out;

  for (const auto &[key, val] : m_node.items()) {
    if (val.is_null())
      continue;

    // Special case: top-level "Module" array -> expand to Module::type{} blocks
    if (key == "Module" && val.is_array()) {
      for (const auto &entry : val) {
        if (!entry.is_object())
          continue;
        const auto type = entry.value("type", std::string{});
        if (type.empty())
          continue;
        out += "Module::" + type + " { ";
        for (const auto &[opt_key, opt_val] : entry.items()) {
          if (opt_key == "type")
            continue;
          out += scalar_or_block(opt_key, opt_val);
        }
        out += "} ";
      }
      continue;
    }

    out += scalar_or_block(key, val);
  }
  return out;
}

// Helper: serialise one key+value pair (used by to_ampsci_string).
// Not a member; defined locally so it can recurse via to_ampsci_string().
inline std::string InputBlock::scalar_or_block(const std::string &key,
                                               const nlohmann::json &val) {
  if (val.is_object()) {
    InputBlock sub{key, val};
    return key + " { " + sub.to_ampsci_string() + "} ";
  }
  if (val.is_string())
    return key + " = " + val.get<std::string>() + "; ";
  if (val.is_boolean())
    return key + " = " + (val.get<bool>() ? "true" : "false") + "; ";
  if (val.is_number())
    return key + " = " + val.dump() + "; ";
  if (val.is_array()) {
    // Scalar-element arrays -> comma-separated list
    std::string list;
    for (const auto &elem : val) {
      if (!elem.is_primitive() || elem.is_null())
        return {}; // array of objects -- skip
      if (!list.empty())
        list += ',';
      list += elem.is_string() ? elem.get<std::string>() : elem.dump();
    }
    return list.empty() ? std::string{} : key + " = " + list + "; ";
  }
  return {};
}

} // namespace IO
