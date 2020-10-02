#include "IO/UserInput.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace IO {

//******************************************************************************
void UserInputBlock::print() const {
  const bool small_block = m_input_options.size() < 4;
  const bool empty = m_input_options.empty();
  if (m_block_name == "")
    return;
  std::cout << m_block_name << " {";
  if (!empty && small_block)
    std::cout << " ";
  if (!empty && !small_block)
    std::cout << "\n  ";

  int count = 0;
  for (const auto &option : m_input_options) {
    if (count % 3 == 0 && count > 0)
      std::cout << "\n  ";
    else
      std::cout << "";
    std::cout << option << "; ";
    ++count;
  }
  if (!small_block)
    std::cout << "\n";
  std::cout << "}\n";
}
//******************************************************************************
UserInputBlock
UserInputBlock::subBlock(const std::string &sub_block_name,
                         const std::vector<std::string> &options) const {
  UserInputBlock sub_block(sub_block_name, {});
  for (const auto &option : options) {
    for (const auto &input_line : m_input_options) {
      if (input_line.substr(0, option.length() + 1) == option + "=")
        sub_block.add(input_line);
    }
  }
  return sub_block;
}
//******************************************************************************
bool UserInputBlock::checkBlock(const std::vector<std::string> &list) const {

  bool all_ok = true;
  for (const auto &entry : m_input_options) {
    // c==20, use stringview here
    const auto pos = entry.find('=');
    const auto option = pos < entry.length() ? entry.substr(0, pos) : entry;
    // For each option in
    const auto isoption = [&](const auto &a) { return option == a; };
    const auto bad_option = !std::any_of(list.begin(), list.end(), isoption);
    auto help = (option == "Help" || option == "help") ? true : false;
    if (bad_option && !help) {
      all_ok = false;
      std::cout << "\n⚠️  WARNING: Unclear input option in " << m_block_name
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
    const auto pos = option.find(in_option + "=");
    if (pos == 0) {
      const auto len = in_option.length() + 1;
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
  const auto inp = IO::FRW::splitInput_byBraces(IO::FRW::readInputFile(infile));
  for (const auto &item : inp) {
    const auto block_name = item.first;
    const auto option_vector = IO::FRW::splitInput_bySemiColon(item.second);

    // If block already exists, adds options to that block.
    // Otherwise, adds new block:
    // (Unless its a module/matrix element block, in which case add as a new
    // block)... horrible hack...
    const auto exists =
        std::find_if(begin(m_blocks), end(m_blocks), [block_name](auto &block) {
          return block.name() == block_name          //
                 && (block_name.find("Module") == 0) //
                 && (block_name.find("MatrixElements") == 0);
        });

    if (exists != end(m_blocks)) {
      exists->add(option_vector);
    } else {
      m_blocks.emplace_back(block_name, std::move(option_vector));
    }
  }
}

//******************************************************************************
const UserInputBlock &UserInput::get(const std::string &in_block) const {
  for (const auto &block : m_blocks) {
    if (in_block == block.name())
      return block;
  }
  // XXX DUMB use optional..
  return m_empty_block;
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
  std::cout << time_date() << '\n';
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

} // namespace IO
