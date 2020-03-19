#include "IO/UserInput.hpp"
#include "IO/FileIO_fileReadWrite.hpp"
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
bool UserInputBlock::checkBlock(const std::vector<std::string> &list) const {

  bool all_ok = true;
  for (const auto &entry : m_input_options) {
    // c==20, use stringview here
    auto pos = entry.find('=');
    auto option = pos < entry.length() ? entry.substr(0, pos) : entry;
    // For each option in
    auto isoption = [&](std::string a) { return option == a; };
    auto bad_option = !std::any_of(list.begin(), list.end(), isoption);
    auto help = (option == "Help" || option == "help") ? true : false;
    if (bad_option && !help) {
      all_ok = false;
      std::cerr << "\nWARNING: Unclear input option in " << m_block_name
                << ": `" << option << "'\n --> " << entry << "\n"
                << "Option may be ignored!\n"
                << "Check spelling (or update list of options)\n";
      help = true;
    }
    if (help) {
      std::cout << "Available " << m_block_name << " options are:\n\n "
                << m_block_name << "{ ";
      std::for_each(list.begin(), list.end(),
                    [](const auto &s) { std::cout << s << "; "; });
      std::cout << " }\n\n";
    }
  }
  return all_ok;
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
const UserInputBlock &UserInput::get(const std::string &in_block) const {
  for (const auto &block : m_blocks) {
    if (in_block == block.name())
      return block;
  }
  std::cerr << "\nFAIL55: Missing required input: " << in_block
            << " (compulsory)\n";
  std::abort();
}

//******************************************************************************
std::vector<UserInputBlock> UserInput::module_list() const {
  std::vector<UserInputBlock> out_module_list;
  for (const auto &block : m_blocks) {
    const auto &name = block.name();
    auto pos1 = name.find("Module::");
    auto pos2 = name.find("MatrixElements::");
    if (pos1 == 0 || pos2 == 0)
      out_module_list.push_back(block);
  }
  return out_module_list;
}

//******************************************************************************
void UserInput::print() const {
  for (const auto &block : m_blocks) {
    block.print();
  }
}
//******************************************************************************
bool UserInput::check(const std::string &in_block,
                      const std::vector<std::string> &options) const {
  for (const auto &block : m_blocks) {
    if (block.name() == in_block)
      return block.checkBlock(options);
  }
  return true;
}
