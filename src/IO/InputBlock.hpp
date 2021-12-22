#pragma once
#include "qip/String.hpp" //for case insensitive
#include <algorithm>
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
//******************************************************************************
//! Removes all white space (space, tab, newline), except for those in quotes
inline std::string removeSpaces(std::string str);

//! Removes all quote marks
inline std::string removeQuoteMarks(std::string str);

//! Removes all c++ style block comments from a string
inline void removeBlockComments(std::string &input);

//! Removes all c++ style comments from a string (block and line)
inline std::string removeComments(const std::string &input);

//! Parses a string to type T by stringstream
template <typename T> inline T parse_str_to_T(const std::string &value_as_str);

//! Parses entire file into string. Note: v. inefficient
inline std::string file_to_string(const std::istream &file);

//! Class to determine if a class template in vector
template <typename T> struct IsVector {
  constexpr static bool v = false;
  using t = T;
};
template <typename T> struct IsVector<std::vector<T>> {
  constexpr static bool v = true;
  // nb: returns conatined type of vector
  using t = T;
};
// e.g.:
// std::cout << std::boolalpha;
// std::cout << IO::IsVector<int>::v << "\n";
// std::cout << IO::IsVector<std::vector<int>>::v << "\n";
// std::cout << IO::IsVector<std::vector<double>>::v << "\n";

//! Prints a line of 'c' characters (dflt '*'), num chars long (dflt 80) to
//! cout
inline void print_line(const char c = '*', const int num = 80) {
  for (int i = 0; i < num; i++)
    std::cout << c;
  std::cout << "\n";
}

//******************************************************************************
inline std::string time_date() {
  const auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[30];
  std::strftime(buffer, 30, "%F %T", localtime(&now));
  return buffer;
}
inline std::string date() {
  const auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[30];
  std::strftime(buffer, 30, "%F", localtime(&now));
  return buffer;
}
inline std::string time() {
  const auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[30];
  std::strftime(buffer, 30, "%T", localtime(&now));
  return buffer;
}

//******************************************************************************
//! Simple struct; holds key-value pair, both strings. == compares key
struct Option {
  std::string key;
  std::string value_str;

  friend bool operator==(Option option, std::string_view tkey) {
    // return option.key == tkey;
    return qip::ci_wildcard_compare(tkey, option.key);
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

//******************************************************************************
//! Holds list of Options, and a list of other InputBlocks. Can be initialised
//! with a list of options, with a string, or from a file (ifstream).
//! Format for input is, e.g.,:
/*!
 BlockName1{
   option1=value1;
   option2=value2;
   InnerBlock{
     option1=v3;
   }
 }

 Note: comparison for block/option names is case insensitive!
*/
// nb: I sepparate the function implementations below (in the header file) and
// mark them as inline. This is for readability only, and ensures this file
// works as a single-file header-only
class InputBlock {
private:
  std::string m_name{};
  std::vector<Option> m_options{};
  std::vector<InputBlock> m_blocks{};

public:
  //! Default constructor: name will be blank
  InputBlock(){};

  //! Construct from literal list of 'Options' (see Option struct)
  InputBlock(std::string_view name, std::initializer_list<Option> options = {})
      : m_name(name), m_options(options) {}

  //! Construct from a string with the correct Block{option=value;} format
  InputBlock(std::string_view name, const std::string &string_input)
      : m_name(name) {
    add(string_input);
  }

  //! Construct from plain text file, in Block{option=value;} format
  InputBlock(std::string_view name, const std::istream &file) : m_name(name) {
    add(file_to_string(file));
  }

  //! Add a new InputBlock (merge: will be merged with existing if names
  //! match)
  inline void add(InputBlock block, bool merge = false);
  inline void merge(InputBlock block) { add(block, true); }
  //! Adds a new option to end of list
  inline void add(Option option);
  inline void add(const std::vector<Option> &options);
  //! Adds options/inputBlocks by parsing a string
  inline void add(const std::string &string, bool merge = false);
  inline void merge(const std::string &string) { add(string, true); }

  std::string_view name() const { return m_name; }
  //! Return const reference to list of options
  const std::vector<Option> &options() const { return m_options; }
  //! Return const reference to list of blocks
  const std::vector<InputBlock> &blocks() const { return m_blocks; }

  //! Comparison of blocks compares the 'name'
  friend inline bool operator==(InputBlock block, std::string_view name);
  friend inline bool operator==(std::string_view name, InputBlock block);
  friend inline bool operator!=(InputBlock block, std::string_view name);
  friend inline bool operator!=(std::string_view name, InputBlock block);

  //! If 'key' exists in the options, returns value. Else, returns
  //! default_value. Note: If two keys with same name, will use the later
  template <typename T> T get(std::string_view key, T default_value) const;

  //! Returns optional value. Contains value if key exists; empty otherwise.
  //! Note: If two keys with same name, will use the later
  template <typename T = std::string>
  std::optional<T> get(std::string_view key) const;

  //! Get value from set of nested blocks. .get({block1,block2},option)
  template <typename T>
  T get(std::initializer_list<std::string> blocks, std::string_view key,
        T default_value) const;
  //! As above, but without default value
  template <typename T>
  std::optional<T> get(std::initializer_list<std::string> blocks,
                       std::string_view key) const;

  //! Returns optional InputBlock. Contains InputBlock if block of given name
  //! exists; empty otherwise.
  inline std::optional<InputBlock> getBlock(std::string_view name) const;
  inline std::optional<InputBlock>
  getBlock(std::initializer_list<std::string> blocks,
           std::string_view name) const;

  //! Get an 'Option' (kay, value) - rarely needed
  inline std::optional<Option> getOption(std::string_view key) const;

  //! Prints options to screen in user-friendly form. Same form as input
  //! string. By default prints to cout, but can be given any ostream
  inline void print(std::ostream &os = std::cout, int indent_depth = 0) const;

  inline bool checkBlock_old(const std::vector<std::string> &list,
                             bool print = false) const;

  //! Check all the options and blocks in this; if any of them are not present
  //! in 'list', then there is likely a spelling error in the input => returns
  //! false, warns user, and prints all options to screen. list is a pair:
  //! {option, description}. Description allws you to explain what each option
  //! is - great for 'self-documenting' code
  //! If print=true - will print all options+descriptions even if all good.
  inline bool
  check(std::initializer_list<std::string> blocks,
        const std::vector<std::pair<std::string, std::string>> &list,
        bool print = false) const;

  //! Override for when condidering current block
  inline bool
  check(const std::vector<std::pair<std::string, std::string>> &list,
        bool print = false) const {
    return checkBlock(list, print);
  }

private:
  inline bool
  checkBlock(const std::vector<std::pair<std::string, std::string>> &list,
             bool print = false) const;

  // Allows returning std::vector: comma-separated list input
  template <typename T>
  std::optional<std::vector<T>> get_vector(std::string_view key) const;

  inline InputBlock *getBlock_ptr(std::string_view name);
  inline const InputBlock *getBlock_cptr(std::string_view name) const;

  inline void add_option(std::string_view in_string);
  inline void add_blocks_from_string(std::string_view string, bool merge);
  inline void consolidate();
};

//******************************************************************************
//******************************************************************************
void InputBlock::add(InputBlock block, bool merge) {
  auto existing_block = getBlock_ptr(block.m_name);
  if (merge && existing_block) {
    existing_block->m_options.insert(existing_block->m_options.end(),
                                     block.m_options.cbegin(),
                                     block.m_options.cend());
  } else {
    m_blocks.push_back(block);
  }
}

//******************************************************************************
void InputBlock::add(Option option) { m_options.push_back(option); }
void InputBlock::add(const std::vector<Option> &options) {
  for (const auto &option : options)
    m_options.push_back(option);
}
//******************************************************************************
void InputBlock::add(const std::string &string, bool merge) {
  add_blocks_from_string(removeQuoteMarks(removeSpaces(removeComments(string))),
                         merge);
}

//******************************************************************************
bool operator==(InputBlock block, std::string_view name) {
  return qip::ci_wildcard_compare(name, block.m_name);
}
bool operator==(std::string_view name, InputBlock block) {
  return block == name;
}
bool operator!=(InputBlock block, std::string_view name) {
  return !(block == name);
}
bool operator!=(std::string_view name, InputBlock block) {
  return !(block == name);
}

//******************************************************************************
template <typename T>
std::optional<T> InputBlock::get(std::string_view key) const {
  if constexpr (IsVector<T>::v) {
    return get_vector<typename IsVector<T>::t>(key);
  } else if constexpr (std::is_same_v<T, bool>) {
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    if (option == m_options.crend())
      return std::nullopt;
    if (qip::ci_wildcard_compare("default", option->value_str) ||
        option->value_str == "")
      return std::nullopt;
    const auto &str = option->value_str;
    if (qip::ci_wildcard_compare("true", str) ||
        qip::ci_wildcard_compare("yes", str) ||
        qip::ci_wildcard_compare("y", str) ||
        qip::ci_wildcard_compare("1", str))
      return true;
    return false;
  } else {
    // Use reverse iterators so that we find _last_ option that matches key
    // i.e., assume later options override earlier ones.
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    if (option == m_options.crend())
      return std::nullopt;
    if (qip::ci_wildcard_compare("default", option->value_str) ||
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
InputBlock::get_vector(std::string_view key) const {
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

template <typename T>
T InputBlock::get(std::string_view key, T default_value) const {
  static_assert(!std::is_same_v<T, const char *>,
                "Cannot use get with const char* - use std::string");
  return get<T>(key).value_or(default_value);
}

template <typename T>
T InputBlock::get(std::initializer_list<std::string> blocks,
                  std::string_view key, T default_value) const {
  return get<T>(blocks, key).value_or(default_value);
}

template <typename T = std::string>
std::optional<T> InputBlock::get(std::initializer_list<std::string> blocks,
                                 std::string_view key) const {
  // Find key in nested blocks
  const InputBlock *pB = this;
  for (const auto &block : blocks) {
    pB = pB->getBlock_cptr(block);
    if (pB == nullptr)
      return std::nullopt;
  }
  return pB->get<T>(key);
}

//******************************************************************************
std::optional<InputBlock> InputBlock::getBlock(std::string_view name) const {
  // note: by copy!
  const auto block = std::find(m_blocks.crbegin(), m_blocks.crend(), name);
  if (block == m_blocks.crend())
    return {};
  return *block;
}

std::optional<InputBlock>
InputBlock::getBlock(std::initializer_list<std::string> blocks,
                     std::string_view name) const {
  // note: by copy!
  const InputBlock *pB = this;
  for (const auto &block : blocks) {
    pB = pB->getBlock_cptr(block);
    if (pB == nullptr)
      return std::nullopt;
  }
  return pB->getBlock(name);
}

//******************************************************************************
std::optional<Option> InputBlock::getOption(std::string_view key) const {
  // Use reverse iterators so that we find _last_ option that matches key
  // i.e., assume later options override earlier ones.
  const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
  if (option != m_options.crend())
    return *option;
  return {};
}

//******************************************************************************
void InputBlock::print(std::ostream &os, int depth) const {

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

//******************************************************************************
bool InputBlock::checkBlock(
    const std::vector<std::pair<std::string, std::string>> &list,
    bool print) const {
  // Check each option NOT each sub block!
  // For each input option stored, see if it is allowed
  // "allowed" means appears in list
  bool all_ok = true;
  for (const auto &option : m_options) {
    const auto is_optionQ = [&](const auto &l) {
      // return option.key == l.first;
      return qip::ci_wildcard_compare(l.first, option.key);
    };
    const auto bad_option =
        !std::any_of(list.cbegin(), list.cend(), is_optionQ);
    const auto help =
        qip::ci_wildcard_compare("help", option.key) ? true : false;
    if (help)
      print = true;
    if (bad_option && !help) {
      all_ok = false;
      std::cout << "\n⚠️  WARNING: Unclear input option in " << m_name
                << ": " << option.key << " = " << option.value_str << ";\n"
                << "Option may be ignored!\n"
                << "Check spelling (or update list of options)\n";
    }
  }

  for (const auto &block : m_blocks) {
    // compares, accounting for wild-cards (*)
    // const auto is_blockQ = [&](const auto &b) {
    //   // check if wildcard:
    //   const auto wc = std::find(b.first.cbegin(), b.first.cend(), '*');
    //   if (wc != b.first.cend()) {
    //     const auto len = std::size_t(std::distance(b.first.cbegin(), wc));
    //     // compare sub-strings (up to wildcard):
    //     return block.name().substr(0, len) == b.first.substr(0, len);
    //   }
    //   return block == b.first;
    // };
    const auto is_blockQ = [&](const auto &b) {
      return qip::ci_wildcard_compare(block.name(), b.first);
    };
    const auto bad_block = !std::any_of(list.cbegin(), list.cend(), is_blockQ);
    if (bad_block) {
      all_ok = false;
      std::cout << "\n⚠️  WARNING: Unclear input block within " << m_name
                << ": " << block.name() << "{}\n"
                << "Block and containing options may be ignored!\n"
                << "Check spelling (or update list of options)\n";
    }
  }

  if (!all_ok || print) {
    std::cout << "\nAvailable " << m_name << " options/blocks are:\n"
              << m_name << "{\n";
    std::for_each(list.cbegin(), list.cend(), [](const auto &s) {
      std::cout << "  " << s.first << ";  // " << s.second << "\n";
    });
    std::cout << "}\n\n";
  }
  return all_ok;
}

bool InputBlock::checkBlock_old(const std::vector<std::string> &list,
                                bool print) const {
  // Check each option NOT each sub block!
  // For each input option stored, see if it is allowed
  // "allowed" means appears in list
  bool all_ok = true;
  for (const auto &option : m_options) {
    // const auto is_optionQ = [&](const auto &l) { return option.key == l; };
    const auto is_optionQ = [&](const auto &l) {
      return qip::ci_wildcard_compare(option.key, l);
    };

    const auto bad_option =
        !std::any_of(list.cbegin(), list.cend(), is_optionQ);
    auto help = qip::ci_wildcard_compare(option.key, "help") ? true : false;
    if (help)
      print = true;
    if (bad_option && !help) {
      all_ok = false;
      std::cout << "\n⚠️  WARNING: Unclear input option in " << m_name
                << ": " << option.key << " = " << option.value_str << ";\n"
                << "Option may be ignored!\n"
                << "Check spelling (or update list of options)\n";
    }
  }

  for (const auto &block : m_blocks) {
    // const auto is_blockQ = [&](const auto &b) { return block == b; };
    const auto is_blockQ = [&](const auto &b) {
      return qip::ci_wildcard_compare(block.name(), b);
    };
    const auto bad_block = !std::any_of(list.cbegin(), list.cend(), is_blockQ);
    if (bad_block) {
      all_ok = false;
      std::cout << "\n⚠️  WARNING: Unclear input block within " << m_name
                << ": " << block.name() << "{}\n"
                << "Block and containing options may be ignored!\n"
                << "Check spelling (or update list of options)\n";
    }
  }

  if (!all_ok || print) {
    std::cout << "\nAvailable " << m_name << " options/blocks are:\n"
              << m_name << "{\n";
    std::for_each(list.cbegin(), list.cend(),
                  [](const auto &s) { std::cout << "  " << s << ";\n"; });
    std::cout << "}\n\n";
  }
  return all_ok;
}

//! Check one of the sub-blocks
bool InputBlock::check(
    std::initializer_list<std::string> blocks,
    const std::vector<std::pair<std::string, std::string>> &list,
    bool print) const {
  // Find key in nested blocks
  const InputBlock *pB = this;
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

//******************************************************************************
void InputBlock::add_blocks_from_string(std::string_view string, bool merge) {

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
        const auto next_end = std::min(string.find('{', next_start),
                                       string.find('}', next_start));
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
          std::cerr << "FAIL 271 in InputBlock::add_blocks_from_string: Depth "
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

//******************************************************************************
void InputBlock::add_option(std::string_view in_string) {
  const auto pos = in_string.find('=');
  const auto option = in_string.substr(0, pos);
  const auto value = pos < in_string.length() ? in_string.substr(pos + 1) : "";
  m_options.push_back({std::string(option), std::string(value)});
}

//******************************************************************************
InputBlock *InputBlock::getBlock_ptr(std::string_view name) {
  auto block = std::find(m_blocks.rbegin(), m_blocks.rend(), name);
  if (block == m_blocks.rend())
    return nullptr;
  return &(*block);
}

const InputBlock *InputBlock::getBlock_cptr(std::string_view name) const {
  auto block = std::find(m_blocks.crbegin(), m_blocks.crend(), name);
  if (block == m_blocks.rend())
    return nullptr;
  return &(*block);
}

//******************************************************************************
void InputBlock::consolidate() {
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

//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
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

//******************************************************************************
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

//******************************************************************************
inline std::string removeComments(const std::string &input) {
  std::string str = "";
  {
    std::string line;
    std::stringstream stream1(input);
    while (std::getline(stream1, line, '\n')) {
      auto comm1 = line.find('!'); // nb: char, NOT string literal!
      auto comm2 = line.find('#');
      auto comm3 = line.find("//"); // str literal here
      auto comm = std::min(comm1, std::min(comm2, comm3));
      str += line.substr(0, comm);
      str += '\n';
    }
  }
  removeBlockComments(str);

  return str;
}

//******************************************************************************
template <typename T> T inline parse_str_to_T(const std::string &value_as_str) {
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

//******************************************************************************
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
