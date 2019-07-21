#include "UserInput.hpp"
#include "FileIO_fileReadWrite.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

UserInput::UserInput(const std::string &infile)                      //
    : m_filename(infile), m_raw_input(FileIO::readInputFile(infile)) //
{
  auto inp = FileIO::splitInput_byBraces(m_raw_input);
  for (const auto &item : inp) {
    auto a = item.first;
    auto b = FileIO::splitInput_bySemiColon(item.second);
    m_user_input.push_back(std::make_pair(a, b));
  }
}

std::stringstream UserInput::find(const std::string &block,
                                  const std::string &option) const {
  auto output = std::stringstream("InputNotFound");
  bool already_found = false; // allows for warning if input given twice
  for (const auto &el : m_user_input) {
    if (block == el.first) {
      for (const auto &op : el.second) {
        auto pos = op.find(option + "=");
        auto len = option.length() + 1;
        if (pos != std::string::npos) {
          output = std::stringstream(op.substr(pos + len));
          if (already_found)
            std::cerr << "Warning: duplicate input for " << block << "/"
                      << option << "\n";
          already_found = true;
        }
      }
    }
  }
  return output;
}

void UserInput::print() const {
  for (const auto &item : m_user_input) {
    std::cout << item.first << ":\n";
    for (const auto &entry : item.second) {
      std::cout << "  " << entry << "\n";
    }
    std::cout << "\n";
  }
}
