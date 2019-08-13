#include "UserInput.hpp"
#include "FileIO_fileReadWrite.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

//******************************************************************************
void UserInputBlock::print() const {
  for (const auto &option : m_input_options) {
    std::cout << option << "\n";
  }
}

//******************************************************************************
std::stringstream
UserInputBlock::find_option(const std::string &in_option) const {
  auto output = std::stringstream("InputNotFound");
  bool already_found = false; // allows for warning if input given twice
  for (const auto &option : m_input_options) {
    auto pos = option.find(in_option + "=");
    if (pos == 0) {
      auto len = in_option.length() + 1;
      output = std::stringstream(option.substr(pos + len));
      if (already_found)
        std::cerr << "Warning: duplicate input for " << m_block_name << "/"
                  << in_option << " (using the latter)" << option << "\n";
      already_found = true;
    }
  }
  return output;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
UserInput::UserInput(const std::string &infile) : m_filename(infile) {
  auto inp = FileIO::splitInput_byBraces(FileIO::readInputFile(infile));
  for (const auto &item : inp) {
    auto block_name = item.first;
    auto option_vector = FileIO::splitInput_bySemiColon(item.second);
    m_blocks.emplace_back(block_name, option_vector);
  }
}

//******************************************************************************
std::vector<UserInputBlock> UserInput::module_list() const {
  std::vector<UserInputBlock> out_module_list;
  for (const auto &block : m_blocks) {
    auto name = block.name();
    auto pos1 = name.find("Module::");
    auto pos2 = name.find("MatrixElements::");
    if (pos1 == 0 || pos2 == 0) {
      out_module_list.push_back(block);
    }
  }
  return out_module_list;
}

//******************************************************************************
void UserInput::print() const {
  for (const auto &block : m_blocks) {
    block.print();
  }
}
