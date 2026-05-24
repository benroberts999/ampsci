#pragma once
#include "fmt/color.hpp"
#include "qip/String.hpp" //for case insensitive
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <istream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace IO {

//==============================================================================
//! Removes all white space (space, tab, newline), except for those in quotes
inline std::string removeSpaces(std::string str);

//! Removes all quote marks
inline std::string removeQuoteMarks(std::string str);

//! Removes all c++ style block comments from a string
inline void removeBlockComments(std::string &input);

//! Removes all c++ style comments from a string (block and line)
inline std::string removeComments(const std::string &input);

//! Expands "#include" files
inline std::string expandIncludes(std::string input);

//! Parses a string to type T by stringstream
template <typename T>
inline T parse_str_to_T(const std::string &value_as_str);

//! Parses entire file into string. Note: v. inefficient
inline std::string file_to_string(const std::istream &file);

//! Compile-time trait: IsVector<T>::v is true if T is std::vector. IsVector<T>::t is the element type.
template <typename T>
struct IsVector {
  constexpr static bool v = false;
  using t = T;
};
template <typename T>
struct IsVector<std::vector<T>> {
  constexpr static bool v = true;
  using t = T;
};

//! Compile-time trait: IsArray<T>::v is true if T is std::array. IsArray<T>::t is the element type; IsArray<T>::size is the extent.
template <typename T>
struct IsArray {
  constexpr static bool v = false;
  using t = T;
  static constexpr std::size_t size = 0;
};
template <typename T, std::size_t N>
struct IsArray<std::array<T, N>> {
  constexpr static bool v = true;
  using t = T;
  static constexpr std::size_t size = N;
};

//==============================================================================
//! Returns current local date and time as a string, e.g. "2026-05-14 13:01:02".
inline std::string time_date() {
  const auto now =
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[30];
  std::strftime(buffer, 30, "%F %T", localtime(&now));
  return buffer;
}
//! Returns current local date as a string, e.g. "2026-05-14".
inline std::string date() {
  const auto now =
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[30];
  std::strftime(buffer, 30, "%F", localtime(&now));
  return buffer;
}
//! Returns current local time as a string, e.g. "13:01:02".
inline std::string time() {
  const auto now =
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[30];
  std::strftime(buffer, 30, "%T", localtime(&now));
  return buffer;
}

//==============================================================================
//! Simple struct; holds key-value pair, both strings. == compares key
struct Option {
  std::string key;
  std::string value_str;

  friend bool operator==(Option option, std::string_view tkey) {
    // return option.key == tkey;
    return qip::ci_wc_compare(tkey, option.key);
  }
  friend bool operator==(std::string_view tkey, Option option) {
    return option == tkey;
  }
  friend bool operator!=(Option option, std::string_view tkey) {
    return !(option == tkey);
  }
  friend bool operator!=(std::string_view tkey, Option option) {
    return !(option == tkey);
  }
};

//==============================================================================
/*!
  @brief Holds a named list of key=value options and nested InputBlocks.
  @details
  Parses and stores structured input in the format:

  ```
  BlockName1 {
    option1 = value1;
    option2 = value2;
    InnerBlock {
      option1 = v3;
    }
  }
  ```

  Can be constructed from a literal option list, a string, or a file stream.
  Block and option name comparison is case-insensitive.

  @note Implementations are inline in the header for single-file header-only use.
*/
// nb: implementations are separated below (inline) for readability only
class InputBlockLegacy {
private:
  std::string m_name{};
  std::vector<Option> m_options{};
  std::vector<InputBlockLegacy> m_blocks{};

public:
  //! Default constructor; name will be blank.
  InputBlockLegacy() {};

  //! Constructs from a literal list of Option structs.
  InputBlockLegacy(std::string_view name,
                   std::initializer_list<Option> options = {})
    : m_name(name), m_options(options) {}

  //! Constructs by parsing @p string_input in Block{option=value;} format.
  InputBlockLegacy(std::string_view name, const std::string &string_input)
    : m_name(name) {
    add(string_input);
  }

  //! Constructs by reading and parsing a plain-text file stream.
  InputBlockLegacy(std::string_view name, const std::istream &file)
    : m_name(name) {
    add(file_to_string(file));
  }

  /*!
    @brief Appends or merges a child InputBlockLegacy.
    @details
    If @p merge is false (default), always appends @p block.
    If @p merge is true and a block with the same name already exists,
    the options from @p block are merged into the existing block instead.
  */
  inline void add(InputBlockLegacy block, bool merge = false);

  //! Merges @p block into an existing block of the same name, or appends it.
  inline void merge(InputBlockLegacy block) { add(block, true); }

  //! Appends a single Option to the option list.
  inline void add(Option option);

  //! Appends a list of Options to the option list.
  inline void add(const std::vector<Option> &options);

  /*!
    @brief Parses @p string and adds its options and blocks.
    @details
    Comments, whitespace, and quote marks are stripped before parsing.
    If @p merge is true, duplicate block names are consolidated rather than
    appended.
  */
  inline void add(const std::string &string, bool merge = false);

  //! Parses @p string and merges any duplicate block names.
  inline void merge(const std::string &string) { add(string, true); }

  //! Returns the name of this block.
  std::string_view name() const { return m_name; }

  //! Returns const reference to the list of options.
  const std::vector<Option> &options() const { return m_options; }

  //! Returns const reference to the list of child blocks.
  const std::vector<InputBlockLegacy> &blocks() const { return m_blocks; }

  //! Equality/inequality compare by block name (case-insensitive).
  friend inline bool operator==(InputBlockLegacy block, std::string_view name);
  friend inline bool operator==(std::string_view name, InputBlockLegacy block);
  friend inline bool operator!=(InputBlockLegacy block, std::string_view name);
  friend inline bool operator!=(std::string_view name, InputBlockLegacy block);

  /*!
    @brief Returns the value of @p key, or @p default_value if not found.
    @details
    If the same key appears more than once, the later occurrence takes
    precedence. For bool, accepts "true"/"yes"/"y" (case-insensitive).
    For std::vector<T> or std::array<T,N>, parses a comma-separated list.
    @note Cannot be used with T = const char*; use std::string instead.
  */
  template <typename T>
  T get(std::string_view key, T default_value) const;

  /*!
    @brief Returns an optional value for @p key; empty if not found or unset.
    @details
    If the same key appears more than once, the later occurrence takes
    precedence. A value of "default" or empty string is treated as unset
    (returns nullopt). Supports T = std::vector<T> or std::array<T,N> for
    comma-separated list values.
  */
  template <typename T = std::string>
  std::optional<T> get(std::string_view key) const;

  //! Returns true if @p key is present in this block's option list, even if unset.
  bool has_option(std::string_view key) const {
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    return !(option == m_options.crend());
  }

  //! Returns true if @p key is present and has a non-default, non-empty value.
  bool option_is_set(std::string_view key) const {
    return !(get(key) == std::nullopt);
  }

  //! Returns value of @p key in a sequence of nested blocks, or @p default_value.
  template <typename T>
  T get(std::initializer_list<std::string> blocks, std::string_view key,
        T default_value) const;

  //! Returns optional value of @p key in a sequence of nested blocks; empty if not found.
  template <typename T>
  std::optional<T> get(std::initializer_list<std::string> blocks,
                       std::string_view key) const;

  //! Returns an optional copy of the child block named @p name; empty if not found.
  inline std::optional<InputBlockLegacy> getBlock(std::string_view name) const;

  //! Returns an optional copy of a block found by traversing @p blocks then looking for @p name.
  inline std::optional<InputBlockLegacy>
  getBlock(std::initializer_list<std::string> blocks,
           std::string_view name) const;

  //! Returns a copy of the child block named @p name, or an empty block if not found.
  inline InputBlockLegacy get_block(std::string_view name) const {
    auto temp = getBlock(name);
    if (temp)
      return *temp;
    return InputBlockLegacy{};
  }

  //! Returns true if a child block named @p name exists in this block.
  bool has_block(std::string_view name) const {
    auto temp = getBlock(name);
    return temp != std::nullopt;
  }

  //! Returns true if a block named @p name exists within the given nested @p blocks.
  bool has_block(std::initializer_list<std::string> blocks,
                 std::string_view name) const {
    auto temp = getBlock(blocks, name);
    return temp != std::nullopt;
  }

  //! Returns the raw Option struct for @p key; rarely needed directly.
  inline std::optional<Option> getOption(std::string_view key) const;

  //! Prints the block contents to @p os in Block{option=value;} format.
  inline void print(std::ostream &os = std::cout, int indent_depth = 0) const;

  /*!
    @brief Validates options and sub-blocks against an allowed list.
    @details
    Checks each option and sub-block in the nested path @p blocks against
    @p list. If any are not found, a warning is printed along with the nearest
    spelling suggestion, and false is returned. If @p print is true, the full
    list of allowed options and descriptions is always printed.

    The list entries are pairs of {option_name, description}. Blocks are
    identified by a trailing `{}` in the name, e.g., "SubBlock{}".
    Descriptions support self-documenting input files.
  */
  inline bool
  check(std::initializer_list<std::string> blocks,
        const std::vector<std::pair<std::string, std::string>> &list,
        bool print = false) const;

  //! Validates options in the current block against @p list. See check() overload for details.
  inline bool
  check(const std::vector<std::pair<std::string, std::string>> &list,
        bool print = false) const {
    return checkBlock(list, print);
  }

private:
  inline bool
  checkBlock(const std::vector<std::pair<std::string, std::string>> &list,
             bool print = false) const;

  //! Return a pointer to a block. Allows editing of blocks
  inline InputBlockLegacy *getBlock_ptr(std::string_view name);
  inline const InputBlockLegacy *getBlock_cptr(std::string_view name) const;

  // Allows returning std::vector: comma-separated list input
  template <typename T>
  std::optional<std::vector<T>> get_vector(std::string_view key) const;

  // Allows returning std::array: comma-separated list input
  template <typename T, std::size_t N>
  std::optional<std::array<T, N>> get_array(std::string_view key) const;

  inline void add_option(std::string_view in_string);
  inline void add_blocks_from_string(std::string_view string, bool merge);
  inline void consolidate();
};

//==============================================================================
//==============================================================================
void InputBlockLegacy::add(InputBlockLegacy block, bool merge) {
  auto existing_block = getBlock_ptr(block.m_name);
  if (merge && existing_block) {
    existing_block->m_options.insert(existing_block->m_options.end(),
                                     block.m_options.cbegin(),
                                     block.m_options.cend());
  } else {
    m_blocks.push_back(block);
  }
}

//==============================================================================
void InputBlockLegacy::add(Option option) { m_options.push_back(option); }
void InputBlockLegacy::add(const std::vector<Option> &options) {
  for (const auto &option : options)
    m_options.push_back(option);
}
//==============================================================================
void InputBlockLegacy::add(const std::string &string, bool merge) {

  add_blocks_from_string(
    removeQuoteMarks(removeSpaces(expandIncludes(removeComments(string)))),
    merge);
}

//==============================================================================
bool operator==(InputBlockLegacy block, std::string_view name) {
  return qip::ci_wc_compare(name, block.m_name);
}
bool operator==(std::string_view name, InputBlockLegacy block) {
  return block == name;
}
bool operator!=(InputBlockLegacy block, std::string_view name) {
  return !(block == name);
}
bool operator!=(std::string_view name, InputBlockLegacy block) {
  return !(block == name);
}

//==============================================================================
template <typename T>
std::optional<T> InputBlockLegacy::get(std::string_view key) const {
  if constexpr (IsVector<T>::v) {
    return get_vector<typename IsVector<T>::t>(key);
  } else if constexpr (IsArray<T>::v) {
    return get_array<typename IsArray<T>::t, IsArray<T>::size>(key);
  } else if constexpr (std::is_same_v<T, bool>) {
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    if (option == m_options.crend())
      return std::nullopt;
    if (qip::ci_wc_compare("default", option->value_str) ||
        option->value_str == "")
      return std::nullopt;
    const auto &str = option->value_str;
    if (qip::ci_wc_compare("true", str) || qip::ci_wc_compare("yes", str) ||
        qip::ci_wc_compare("y", str))
      return true;
    return false;
  } else {
    // Use reverse iterators so that we find _last_ option that matches key
    // i.e., assume later options override earlier ones.
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    if (option == m_options.crend())
      return std::nullopt;
    if (qip::ci_wc_compare("default", option->value_str) ||
        option->value_str == "")
      return std::nullopt;
    return parse_str_to_T<T>(option->value_str);
  }
}

// special function; allows return of std::vector (for comma-separated list
// input). Optional of vector is kind of redundant, but is this way so it
// aligns with the other functions (checks if optional is empty when deciding
// if should return the default value)
template <typename T>
std::optional<std::vector<T>>
InputBlockLegacy::get_vector(std::string_view key) const {
  // Use reverse iterators so that we find _last_ option that matches key
  // i.e., assume later options override earlier ones.
  std::vector<T> out;
  const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
  if (option == m_options.crend())
    return std::nullopt;
  if (option->value_str == "")
    return std::nullopt;
  std::stringstream ss(option->value_str);
  while (ss.good()) {
    // note: *very* innefficient
    std::string substr;
    std::getline(ss, substr, ',');
    out.push_back(parse_str_to_T<T>(substr));
  }
  return out;
}

template <typename T, std::size_t N>
std::optional<std::array<T, N>>
InputBlockLegacy::get_array(std::string_view key) const {
  // Use reverse iterators so that we find _last_ option that matches key
  // i.e., assume later options override earlier ones.
  std::array<T, N> out;
  const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
  if (option == m_options.crend())
    return std::nullopt;
  if (option->value_str == "")
    return std::nullopt;
  std::stringstream ss(option->value_str);
  std::size_t index = 0;
  while (ss.good()) {
    // note: *very* innefficient
    std::string substr;
    std::getline(ss, substr, ',');
    // out.push_back(parse_str_to_T<T>(substr));
    out.at(index) = parse_str_to_T<T>(substr);
    ++index;
  }
  return out;
}

template <typename T>
T InputBlockLegacy::get(std::string_view key, T default_value) const {
  static_assert(!std::is_same_v<T, const char *>,
                "Cannot use get with const char* - use std::string");
  return get<T>(key).value_or(default_value);
}

template <typename T>
T InputBlockLegacy::get(std::initializer_list<std::string> blocks,
                        std::string_view key, T default_value) const {
  return get<T>(blocks, key).value_or(default_value);
}

template <typename T>
std::optional<T>
InputBlockLegacy::get(std::initializer_list<std::string> blocks,
                      std::string_view key) const {
  // Find key in nested blocks
  const InputBlockLegacy *pB = this;
  for (const auto &block : blocks) {
    pB = pB->getBlock_cptr(block);
    if (pB == nullptr)
      return std::nullopt;
  }
  return pB->get<T>(key);
}

//==============================================================================
std::optional<InputBlockLegacy>
InputBlockLegacy::getBlock(std::string_view name) const {
  // note: by copy!
  const auto block = std::find(m_blocks.crbegin(), m_blocks.crend(), name);
  if (block == m_blocks.crend())
    return {};
  return *block;
}

std::optional<InputBlockLegacy>
InputBlockLegacy::getBlock(std::initializer_list<std::string> blocks,
                           std::string_view name) const {
  // note: by copy!
  const InputBlockLegacy *pB = this;
  for (const auto &block : blocks) {
    pB = pB->getBlock_cptr(block);
    if (pB == nullptr)
      return std::nullopt;
  }
  return pB->getBlock(name);
}

//==============================================================================
std::optional<Option> InputBlockLegacy::getOption(std::string_view key) const {
  // Use reverse iterators so that we find _last_ option that matches key
  // i.e., assume later options override earlier ones.
  const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
  if (option != m_options.crend())
    return *option;
  return {};
}

//==============================================================================
void InputBlockLegacy::print(std::ostream &os, int depth) const {

  std::string indent = "";
  for (int i = 1; i < depth; ++i)
    indent += "  ";

  // Don't print outer-most name
  if (depth != 0)
    os << indent << m_name << " { ";

  const auto multi_entry = (!m_blocks.empty() || (m_options.size() > 1));

  if (depth != 0 && multi_entry)
    os << "\n";

  for (const auto &[key, value] : m_options) {
    os << (depth != 0 && multi_entry ? indent + "  " : "");
    if (value == "")
      os << key << ';';
    else
      os << key << " = " << value << ';';
    os << (multi_entry ? '\n' : ' ');
  }

  for (const auto &block : m_blocks)
    block.print(os, depth + 1);

  if (depth != 0 && multi_entry)
    os << indent;

  if (depth != 0)
    os << "}\n";
}

//==============================================================================
bool InputBlockLegacy::checkBlock(
  const std::vector<std::pair<std::string, std::string>> &list,
  bool print) const {
  // Check each option NOT each sub block!
  // For each input option stored, see if it is allowed
  // "allowed" means appears in list
  bool all_ok = true;
  for (const auto &option : m_options) {
    const auto is_optionQ = [&](const auto &l) {
      // return option.key == l.first;
      return qip::ci_wc_compare(l.first, option.key);
    };
    const auto bad_option =
      !std::any_of(list.cbegin(), list.cend(), is_optionQ);
    const auto help = qip::ci_wc_compare("help", option.key) ? true : false;
    if (help)
      print = true;
    if (bad_option && !help) {
      all_ok = false;
      fmt2::styled_print(fg(fmt::color::orange), "\nWARNING\n");
      std::cout << "Unclear input option in " << m_name << ": " << option.key
                << " = " << option.value_str << ";\n"
                << "Option may be ignored!\n"
                << "Check spelling (or update list of options)\n";
      // spell-check + nearest suggestion:
      auto compare_sc = [&option](const auto &s1, const auto &s2) {
        return qip::ci_Levenstein(s1.first, option.key) <
               qip::ci_Levenstein(s2.first, option.key);
      };
      auto guess = std::min_element(list.cbegin(), list.cend(), compare_sc);
      if (guess != list.cend()) {
        std::cout << "\nDid you mean: " << guess->first << " ?\n";
      }
    }
  }

  using namespace std::string_literals;
  for (const auto &block : m_blocks) {
    const auto is_blockQ = [&](const auto &b) {
      return qip::ci_wc_compare(std::string{block.name()} + "{}"s, b.first);
    };
    const auto bad_block = !std::any_of(list.cbegin(), list.cend(), is_blockQ);
    if (bad_block) {
      all_ok = false;
      fmt2::styled_print(fg(fmt::color::orange), "\nWARNING\n");
      std::cout << "Unclear input block within " << m_name << ": "
                << block.name() << "{}\n"
                << "Block and containing options may be ignored!\n"
                << "Check spelling (or update list of options)\n";
      // spell-check + nearest suggestion:
      auto compare_sc = [&block](const auto &s1, const auto &s2) {
        return qip::ci_Levenstein(s1.first, block.name()) <
               qip::ci_Levenstein(s2.first, block.name());
      };
      auto guess = std::min_element(list.cbegin(), list.cend(), compare_sc);
      if (guess != list.cend()) {
        std::cout << "\nDid you mean: " << guess->first << " ?\n";
      }
    }
  }

  if (!all_ok || print) {
    fmt2::styled_print(fg(fmt::color::light_blue),
                       "\n// Available {} options/blocks\n", m_name);
    fmt2::styled_print(fmt::emphasis::bold, m_name);
    std::cout << "{\n";
    std::for_each(list.cbegin(), list.cend(), [](const auto &s) {
      const auto option_is_block = s.first.back() == '}';
      const auto leading_spaces = s.first.empty() ? ""s : "  ";
      if (!s.first.empty()) {
        std::cout << "  " << s.first << (option_is_block ? "\n" : ";\n");
      }
      if (!s.second.empty()) {
        fmt2::styled_print(fg(fmt::color::light_blue), "{}\n",
                           qip::wrap(s.second, 60, leading_spaces + "  // "));
      } else {
        std::cout << "\n";
      }
    });
    std::cout << "}\n\n";
  }
  return all_ok;
}

//! Check one of the sub-blocks
bool InputBlockLegacy::check(
  std::initializer_list<std::string> blocks,
  const std::vector<std::pair<std::string, std::string>> &list,
  bool print) const {
  // Find key in nested blocks
  const InputBlockLegacy *pB = this;
  for (const auto &block : blocks) {
    pB = pB->getBlock_cptr(block);
    if (pB == nullptr) {
      // Did not fund nested block... may be fine
      // Return true, since a missing block is not an issue (or sepparate
      // issue) We are checking to see if blocks exist that shouldn't
      return true;
    }
  }
  return pB->check(list, print);
}

//==============================================================================
void InputBlockLegacy::add_blocks_from_string(std::string_view string,
                                              bool merge) {

  // Expects that string has comments and spaces removed already

  auto start = 0ul;
  while (start < string.length()) {

    // Find the first of either next ';' or open '{'
    // This is the end of the next input option, or start of block
    auto end = std::min(string.find(';', start), string.find('{', start));
    if (end > string.length() || start >= end)
      break;

    if (string.at(end) == ';') {
      // end of option:

      this->add_option(string.substr(start, end - start));

    } else {
      // start of block

      // 'name' directly preceeds "{"
      const auto block_name = string.substr(start, end - start);
      start = end + 1;

      // Now, find *matching* close '}' - ensure balanced
      int depth_count = 1;
      auto next_start = start; // + 1;
      while (depth_count != 0) {
        if (next_start > string.length())
          break;
        const auto next_end =
          std::min(string.find('{', next_start), string.find('}', next_start));
        if (next_end > string.length())
          break;

        // count depth of bracket nesting:
        if (string.at(next_end) == '{')
          ++depth_count;
        else
          --depth_count;

        if (depth_count == 0) {
          end = next_end;
          break;
        }
        if (depth_count > 100) {
          std::cerr
            << "FAIL 271 in InputBlockLegacy::add_blocks_from_string: Depth "
               "error. Check balanced {} in input\n";
          end = next_end;
          break;
        }
        next_start = next_end + 1;
      }

      // Add a new block, populate it with string. Recursive, since blocks may
      // contain blocks
      auto &block = m_blocks.emplace_back(block_name);

      if (end > start)
        block.add_blocks_from_string(string.substr(start, end - start), merge);
    }

    start = end + 1;
  }

  // Merge duplicated blocks.
  if (merge)
    consolidate();
  // No - want ability to have multiple blocks of same name
}

//==============================================================================
void InputBlockLegacy::add_option(std::string_view in_string) {
  const auto pos = in_string.find('=');
  const auto option = in_string.substr(0, pos);
  const auto value = pos < in_string.length() ? in_string.substr(pos + 1) : "";
  m_options.push_back({std::string(option), std::string(value)});
}

//==============================================================================
InputBlockLegacy *InputBlockLegacy::getBlock_ptr(std::string_view name) {
  auto block = std::find(m_blocks.rbegin(), m_blocks.rend(), name);
  if (block == m_blocks.rend())
    return nullptr;
  return &(*block);
}

const InputBlockLegacy *
InputBlockLegacy::getBlock_cptr(std::string_view name) const {
  auto block = std::find(m_blocks.crbegin(), m_blocks.crend(), name);
  if (block == m_blocks.rend())
    return nullptr;
  return &(*block);
}

//==============================================================================
void InputBlockLegacy::consolidate() {
  for (auto bl = m_blocks.end() - 1; bl != m_blocks.begin() - 1; --bl) {
    bl->consolidate();
    auto bl2 = std::find(m_blocks.begin(), bl, bl->name());
    if (bl2 != bl) {
      bl2->m_options.insert(bl2->m_options.end(), bl->m_options.cbegin(),
                            bl->m_options.cend());
      m_blocks.erase(bl);
    }
  }
}

//==============================================================================
//==============================================================================
//==============================================================================
inline std::string expandIncludes(std::string str) {
  const std::string include_text = "#include";

  for (auto ipos = str.find(include_text); ipos != std::string::npos;
       ipos = str.find(include_text)) {
    const auto start = std::min(str.find('"', ipos), str.find('<', ipos));
    const auto end =
      std::min(str.find('"', start + 1), str.find('>', start + 1));
    const auto fname = str.substr(start + 1, end - start - 1);
    str.erase(ipos, end - ipos + 1);
    std::ifstream ifile(fname);
    if (ifile.good()) {
      str.insert(ipos, removeComments(file_to_string(ifile)));
    }
  }

  return str;
}

//==============================================================================
inline std::string removeSpaces(std::string str) {

  bool inside = false;
  auto lambda = [&inside](unsigned char x) {
    if (x == '\"' || x == '\'')
      inside = !inside;
    return ((x == ' ' || x == '\t' || x == '\n') && !inside);
  };

  str.erase(std::remove_if(str.begin(), str.end(), lambda), str.end());

  return str;
}

inline std::string removeQuoteMarks(std::string str) {

  // remove ' and "
  str.erase(std::remove_if(str.begin(), str.end(),
                           [](unsigned char x) { return x == '\''; }),
            str.end());
  str.erase(std::remove_if(str.begin(), str.end(),
                           [](unsigned char x) { return x == '\"'; }),
            str.end());

  return str;
}

//==============================================================================
inline void removeBlockComments(std::string &input) {
  for (auto posi = input.find("/*"); posi != std::string::npos;
       posi = input.find("/*")) {
    auto posf = input.find("*/");
    if (posf != std::string::npos) {
      input = input.substr(0, posi) + input.substr(posf + 2);
    } else {
      input = input.substr(0, posi);
    }
  }
}

//==============================================================================
inline std::string removeComments(const std::string &input) {
  std::string str = "";
  {
    std::string line;
    std::stringstream stream1(input);
    while (std::getline(stream1, line, '\n')) {
      // const auto comm1 = line.find('!');  // nb: char, NOT string literal!
      // auto comm2 = line.find('#');
      const auto comm3 = line.find("//"); // str literal here
      // const auto comm = std::min({comm1, comm3});
      str += line.substr(0, comm3);
      str += '\n';
    }
  }
  removeBlockComments(str);

  return str;
}

//==============================================================================
template <typename T>
T inline parse_str_to_T(const std::string &value_as_str) {
  if constexpr (std::is_same_v<T, std::string>) {
    // already a string, just return value
    return value_as_str;
  } else {
    // T is not a string: convert using stringstream
    T value_T;
    std::stringstream ss(value_as_str);
    ss >> value_T;
    return value_T;
  }
}

//==============================================================================
inline std::string file_to_string(const std::istream &file) {
  std::string out;
  if (!file)
    return "";
  // Horribly inneficient...
  std::ostringstream ss;
  ss << file.rdbuf();
  return ss.str();
}

} // namespace IO
