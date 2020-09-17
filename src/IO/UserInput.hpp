#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

// Set 'help' option to print out all available options

namespace IO {

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
// Move to qip?:

template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
bool acceptableQ(const T &value, const std::pair<T, T> &list,
                 const bool print = true) {
  // If list.size()==0, and is a numeric type, treat as range
  // Otherwise, treat as exclusive list
  const auto &[min, max] = list;
  const auto ok = value >= min && value <= max;
  if (!ok && print) {
    std::cout << "\n⚠️  WARNING unacceptable value: " << value
              << " not in required range: [" << min << "," << max << "]\n";
  }
  return ok;
}

template <typename T, typename List>
bool acceptableQ(const T &value, const List &list, const bool print = true) {
  // If list.size()==0, and is a numeric type, treat as range
  // Otherwise, treat as exclusive list
  const auto ok = std::find(cbegin(list), cend(list), value) != cend(list);
  if (!ok && print) {
    std::cout << "\n⚠️  WARNING unacceptable value: " << value
              << " not in required set:\n";
    std::for_each(cbegin(list), cend(list),
                  [](auto &x) { std::cout << x << ", "; });
    std::cout << "\n";
  }
  return ok;
}

//******************************************************************************
//! Stores input options for a specific 'block'; usually use 'UserInput'
//! (instead of 'UserInputBlock') - has essentially the same functionality (but
//! don't need to specify block name on '.get')
class UserInputBlock {

public:
  UserInputBlock(const std::string &in_block_name,
                 const std::vector<std::string> &in_input_options)
      : m_block_name(in_block_name), m_input_options(in_input_options) {}
  UserInputBlock(const std::string &in_block_name,
                 const UserInputBlock &in_block)
      : m_block_name(in_block_name),
        m_input_options(in_block.m_input_options) {}

  template <typename T>
  T get(const std::string &option, const T &default_value) const;
  template <typename T> T get(const std::string &option) const;

  template <typename T>
  std::vector<T> get_list(const std::string &option,
                          const std::vector<T> &default_value) const;
  template <typename T>
  std::vector<T> get_list(const std::string &option) const;

  const std::string &name() const { return m_block_name; }

  void print() const;
  bool checkBlock(const std::vector<std::string> &options) const;
  void add(const std::string &new_in) { m_input_options.push_back(new_in); }

  UserInputBlock copy_with(const std::string &new_in) const {
    auto out = *this;
    out.add(new_in);
    return out;
  }

  UserInputBlock subBlock(const std::string &sub_block_name,
                          const std::vector<std::string> &options) const;

private:
  std::string m_block_name;
  std::vector<std::string> m_input_options;

  std::stringstream find_option(const std::string &in_option) const;
};

//******************************************************************************
//! Stores + retrives user input
class UserInput {

public:
  //! Constructor: takes in a plane text file.
  UserInput(const std::string &infile);

  //! Looks up option with default value.
  //! @details If requested option not in input, returns default value. Return
  //! type is specified by the default value (or, explicitely)
  template <typename T>
  T get(const std::string &in_block, const std::string &option,
        const T &default_value) const;

  //! Looks up option with no default. If option not in input, will warn and
  //! abort
  template <typename T>
  T get(const std::string &in_block, const std::string &option) const;

  //! As above, but for a list (vector) - all items in list must be of same
  //! type.
  template <typename T>
  std::vector<T> get_list(const std::string &in_block,
                          const std::string &option,
                          const std::vector<T> &default_value) const;
  //! As above, but for a list (vector) - all items in list must be of same
  //! type.
  template <typename T>
  std::vector<T> get_list(const std::string &in_block,
                          const std::string &option) const;

  //! Returns one of the sub blocks (by const ref)
  const UserInputBlock &get(const std::string &in_block) const;

  //! Returns list of all blocks that are 'Modules' or 'MatrixElements' (copy)
  std::vector<UserInputBlock> module_list() const;

  //! Prints all input options in nice format.
  //! @details Note: out format same as input format, so can be copied into new
  //! input file
  void print() const;

  //! Checker. Returns false if any of the input options in the 'in_block' block
  //! are not listed in the given 'options' list
  //! @details This checks to ensure we aren't trying to set an option that
  //! doesn't exist, or if we've spelled an option wrong. Otherwise, these would
  //! simply be ignored and default values used in the calculations, which is
  //! likely an error (and an easy one to miss)
  bool check(const std::string &in_block,
             const std::vector<std::string> &options) const;

private:
  std::string m_filename;
  std::vector<UserInputBlock> m_blocks = {};
  UserInputBlock m_empty_block{"", {}};
};

//******************************************************************************
//******************************************************************************
namespace UserInputHelper {
template <typename T>
inline T get_impl(std::stringstream &ss, const std::string &in) {
  T val;
  ss >> val;
  if (ss.fail()) {
    std::cerr << "\nWARNING 78 in UserInput: " << in << "=" << ss.str()
              << " invalid?\n";
    std::abort();
  }
  return val;
}

template <typename T>
inline std::vector<T> get_list_impl(std::stringstream &ss,
                                    const std::string &in) {
  std::vector<T> result;
  bool fail = false;

  while (ss.good()) {
    std::string each_str = "";
    getline(ss, each_str, ',');
    std::stringstream ss_each(each_str);
    T val;
    ss_each >> val;
    result.push_back(val);
    if (ss.fail())
      fail = true;
  }

  if (ss.fail() || fail) {
    std::cerr << "\nWARNING 78 in UserInput: " << in << "=" << ss.str()
              << " invalid?\n";
    std::abort();
  }
  return result;
}

template <> inline bool get_impl(std::stringstream &ss, const std::string &) {
  if (ss.str() == "false" || ss.str() == "False" || ss.str() == "0" ||
      ss.str() == "No" || ss.str() == "no")
    return false;
  return true;
}
} // namespace UserInputHelper

//******************************************************************************
template <typename T>
T UserInputBlock::get(const std::string &option, const T &default_value) const {
  auto option_ss = find_option(option);
  if (option_ss.str() == "InputNotFound" || option_ss.str() == "default" ||
      option_ss.str() == "dflt")
    return default_value;
  return UserInputHelper::get_impl<T>(option_ss, m_block_name + '/' + option);
}
template <typename T>
T UserInputBlock::get(const std::string &option) const
// No default value; user input is complulsory
{
  auto option_ss = find_option(option);
  if (option_ss.str() == "InputNotFound") {
    std::cerr << "\nError: Missing required input: " << m_block_name << "/"
              << option << " (compulsory)\n";
    // std::abort();
    std::cout << "Enter input value now, or ctrl+c to quit:\n" << option << "=";
    option_ss.str("");
    std::string tmp;
    std::cin >> tmp;
    std::cout << option << "=" << tmp << "\n";
    option_ss << tmp;
  }
  return UserInputHelper::get_impl<T>(option_ss, m_block_name + '/' + option);
}

//------------------------------------------------------------------------------
template <typename T>
std::vector<T>
UserInputBlock::get_list(const std::string &option,
                         const std::vector<T> &default_value) const {
  auto option_ss = find_option(option);
  if (option_ss.str() == "InputNotFound" || option_ss.str() == "default" ||
      option_ss.str() == "dflt")
    return default_value;
  return UserInputHelper::get_list_impl<T>(option_ss,
                                           m_block_name + '/' + option);
}
template <typename T>
std::vector<T> UserInputBlock::get_list(const std::string &option) const
// No default value; user input is complulsory
{
  auto option_ss = find_option(option);
  if (option_ss.str() == "InputNotFound") {
    std::cerr << "\nError: Missing required input: " << m_block_name << "/"
              << option << " (compulsory)\n";
    std::abort();
  }
  return UserInputHelper::get_list_impl<T>(option_ss,
                                           m_block_name + '/' + option);
}

//******************************************************************************
template <typename T>
T UserInput::get(const std::string &in_block, const std::string &option,
                 const T &default_value) const {
  for (const auto &block : m_blocks) {
    if (in_block == block.name())
      return block.get<T>(option, default_value);
  }
  return default_value;
}

template <typename T>
T UserInput::get(const std::string &in_block, const std::string &option) const {
  for (const auto &block : m_blocks) {
    if (in_block == block.name())
      return block.get<T>(option);
  }
  std::cerr << "\nFAIL: Missing required input: " << in_block << "/" << option
            << " (compulsory)\n";
  std::abort();
}
//------------------------------------------------------------------------------
template <typename T>
std::vector<T> UserInput::get_list(const std::string &in_block,
                                   const std::string &option,
                                   const std::vector<T> &default_value) const {
  for (const auto &block : m_blocks) {
    if (in_block == block.name())
      return block.get_list<T>(option, default_value);
  }
  return default_value;
}

template <typename T>
std::vector<T> UserInput::get_list(const std::string &in_block,
                                   const std::string &option) const {
  for (const auto &block : m_blocks) {
    if (in_block == block.name())
      return block.get_list<T>(option);
  }
  std::cerr << "\nFAIL: Missing required input: " << in_block << "/" << option
            << " (compulsory)\n";
  std::abort();
}

} // namespace IO
