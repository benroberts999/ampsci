#pragma once
// #include "FileIO_fileReadWrite.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

//******************************************************************************
class UserInputBlock {

public:
  UserInputBlock(const std::string &in_block_name,
                 const std::vector<std::string> &in_input_options)
      : m_block_name(in_block_name), m_input_options(in_input_options) {}

  template <typename T>
  T get(const std::string &option, const T &default_value) const;
  template <typename T> T get(const std::string &option) const;

  const std::string &name() const { return m_block_name; }

  void print() const;

private:
  std::string m_block_name;
  std::vector<std::string> m_input_options;

  std::stringstream find_option(const std::string &in_option) const;
};

//******************************************************************************
class UserInput {

public:
  UserInput(const std::string &infile);

  template <typename T>
  T get(const std::string &in_block, const std::string &option,
        const T &default_value) const;

  template <typename T>
  T get(const std::string &in_block, const std::string &option) const;

  std::vector<UserInputBlock> module_list() const;

  void print() const;

private:
  std::string m_filename;
  std::vector<UserInputBlock> m_blocks;
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

template <> inline bool get_impl(std::stringstream &ss, const std::string &in) {
  if (ss.str() == "true" || ss.str() == "True" || ss.str() == "1" ||
      ss.str() == "Yes" || ss.str() == "yes")
    return true;
  if (ss.str() == "false" || ss.str() == "False" || ss.str() == "0" ||
      ss.str() == "No" || ss.str() == "no")
    return false;
  std::cerr << "\nWARNING 44 in UserInput: " << in << "=" << ss.str()
            << " invalid?\n";
  std::abort();
}
} // namespace UserInputHelper

//******************************************************************************
template <typename T>
T UserInputBlock::get(const std::string &option, const T &default_value) const {
  auto option_ss = find_option(option);
  if (option_ss.str() == "InputNotFound")
    return default_value;
  return UserInputHelper::get_impl<T>(option_ss, m_block_name + '/' + option);
}
template <typename T>
T UserInputBlock::get(const std::string &option) const
// No default value; user input is complulsory
{
  auto option_ss = find_option(option);
  if (option_ss.str() == "InputNotFound") {
    std::cerr << "\nFAIL: Missing required input: " << m_block_name << "/"
              << option << " (compulsory)\n";
    std::abort();
  }
  return UserInputHelper::get_impl<T>(option_ss, m_block_name + '/' + option);
}

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
