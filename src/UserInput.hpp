#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

template <typename T> T get_impl(std::stringstream &ss, const std::string &in);

class UserInput {

  using StringPair = std::pair<std::string, std::string>;
  using StringVectorPair = std::pair<std::string, std::vector<std::string>>;

public:
  UserInput(const std::string &infile);

private:
  const std::string m_filename;
  const std::string m_raw_input;
  std::vector<StringVectorPair> m_user_input;

private:
  std::stringstream find(const std::string &block,
                         const std::string &option) const;

public:
  void print() const;

  std::vector<std::string> module_list() const {
    std::vector<std::string> output;
    for (const auto &entry : m_user_input) {
      auto block = entry.first;
      auto pos = block.find("Module::");
      if (pos == 0) {
        output.push_back(block);
      }
    }
    return output;
  }

  template <typename T>
  inline T get(const std::string &block, const std::string &option,
               const T &default_value) const {
    auto a = find(block, option);
    if (a.str() == "InputNotFound")
      return default_value;
    return get_impl<T>(a, block + "/" + option);
  }

  template <typename T>
  inline T get(const std::string &block, const std::string &option) const
  // No default value; user input is complulsory
  {
    auto a = find(block, option);
    if (a.str() == "InputNotFound") {
      std::cerr << "\nFAIL: Missing required input: " << block << "/" << option
                << " (compulsory)\n";
      std::abort();
    }
    return get_impl<T>(a, block + "/" + option);
  }
};

//******************************************************************************
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
