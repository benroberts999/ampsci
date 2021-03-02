#pragma once
#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <istream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace IO {

//******************************************************************************
//! Removes all white space (space, tab, newline), AND quote marks!
inline std::string removeSpaces(std::string lines);

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

//******************************************************************************
//! Prints a line of 'c' characters (dflt '*'), num chars long (dflt 80) to cout
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

//******************************************************************************
//! Simple struct; holds key-value pair, both strings. == compares key
struct Option {
  std::string key;
  std::string value_str;

  friend bool operator==(Option option, std::string_view tkey) {
    return option.key == tkey;
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
*/
class InputBlock {
private:
  std::string m_name{};
  std::vector<Option> m_options{};
  std::vector<InputBlock> m_blocks{};

public:
  InputBlock(){};

  InputBlock(std::string_view name, std::initializer_list<Option> options = {})
      : m_name(name), m_options(options) {}

  InputBlock(std::string_view name, const std::string &string_input)
      : m_name(name) {
    add(string_input);
  }

  InputBlock(std::string_view name, const std::istream &file) : m_name(name) {
    add(file_to_string(file));
  }

  //! Add a new InputBlock (will be merged with existing if names match)
  void add(InputBlock block);
  //! Adds a new option to end of list
  void add(Option option);
  void add(const std::vector<Option> &options);
  //! Adds options/inputBlocks by parsing a string
  void add(const std::string &string);

  std::string_view name() const { return m_name; }
  const std::vector<Option> &options() const { return m_options; }
  const std::vector<InputBlock> &blocks() const { return m_blocks; }

  //! Comparison of blocks compares the 'name'
  friend bool operator==(InputBlock block, std::string_view name);
  friend bool operator==(std::string_view name, InputBlock block);
  friend bool operator!=(InputBlock block, std::string_view name);
  friend bool operator!=(std::string_view name, InputBlock block);

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
  std::optional<InputBlock> getBlock(std::string_view name) const;

  //! Get an 'Option' (kay, value) - rarely needed
  std::optional<Option> getOption(std::string_view key) const;

  //! Prints options to screen in user-friendly form. Same form as input string.
  //! By default prints to cout, but can be given any ostream
  void print(std::ostream &os = std::cout, int indent_depth = 0) const;

  //! Check all the options and blocks in this; if any of them are not present
  //! in 'list', then there is likely a spelling error in the input => returns
  //! false, warns user, and prints all options to screen. list is a pair:
  //! {option, description}. Description allws you to explain what each option
  //! is - great for 'self-documenting' code
  //! If print=true - will print all options+descriptions even if all good.

  // bool checkBlock(const std::vector<std::pair<std::string, std::string>>
  // &list,
  //                 bool print = false) const;
  //
  // bool check(std::initializer_list<std::string> blocks,
  //            const std::vector<std::pair<std::string, std::string>> &list,
  //            bool print = false) const;

  bool checkBlock(const std::vector<std::string> &list,
                  bool print = false) const;

  bool check(std::initializer_list<std::string> blocks,
             const std::vector<std::string> &list, bool print = false) const;

private:
  // Allows returning std::vector: comma-separated list input
  template <typename T>
  std::optional<std::vector<T>> get_vector(std::string_view key) const;

  InputBlock *getBlock_ptr(std::string_view name);
  const InputBlock *getBlock_cptr(std::string_view name) const;

  void add_option(std::string_view in_string);
  void add_blocks_from_string(std::string_view string);
  void consolidate();
};

//******************************************************************************

//******************************************************************************
template <typename T> // typename T = std::string
std::optional<T> InputBlock::get(std::string_view key) const {
  if constexpr (IsVector<T>::v) {
    return get_vector<typename IsVector<T>::t>(key);
  } else if constexpr (std::is_same_v<T, bool>) {
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    if (option == m_options.crend())
      return std::nullopt;
    const auto &str = option->value_str;
    return (str == "True" || str == "true" || str == "Yes" || str == "yes" ||
            str == "1" || str == "Y" || str == "y");
  } else {
    // Use reverse iterators so that we find _last_ option that matches key
    // i.e., assume later options override earlier ones.
    const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
    if (option == m_options.crend())
      return std::nullopt;
    return parse_str_to_T<T>(option->value_str);
  }
}

// special function; allows return of std::vector (for comma-separated list
// input). Optional of vector is kind of redundant, but is this way so it aligns
// with the other functions (checks if optional is empty when deciding if should
// return the default value)
template <typename T>
std::optional<std::vector<T>>
InputBlock::get_vector(std::string_view key) const {
  // Use reverse iterators so that we find _last_ option that matches key
  // i.e., assume later options override earlier ones.
  std::vector<T> out;
  const auto option = std::find(m_options.crbegin(), m_options.crend(), key);
  if (option == m_options.crend())
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
//******************************************************************************
//******************************************************************************

//******************************************************************************
inline std::string removeSpaces(std::string lines) {

  // remove spaces
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == ' '; }),
              lines.end());
  // remove tabs
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == '\t'; }),
              lines.end());
  // remove newlines
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == '\n'; }),
              lines.end());

  // remove ' and "
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == '\''; }),
              lines.end());
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == '\"'; }),
              lines.end());

  return lines;
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
  std::string lines = "";
  {
    std::string line;
    std::stringstream stream1(input);
    while (std::getline(stream1, line, '\n')) {
      auto comm1 = line.find('!'); // nb: char, NOT string literal!
      auto comm2 = line.find('#');
      auto comm3 = line.find("//"); // str literal here
      auto comm = std::min(comm1, std::min(comm2, comm3));
      lines += line.substr(0, comm);
      lines += '\n';
    }
  }
  removeBlockComments(lines);

  return lines;
}

//******************************************************************************
template <typename T> T parse_str_to_T(const std::string &value_as_str) {
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
