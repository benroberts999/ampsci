#include "InputBlock.hpp"
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
void InputBlock::add(InputBlock block) {
  auto existing_block = getBlock_ptr(block.m_name);
  if (existing_block) {
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
void InputBlock::add(const std::string &string) {
  add_blocks_from_string(removeSpaces(removeComments(string)));
}

//******************************************************************************
bool operator==(InputBlock block, std::string_view name) {
  return block.m_name == name;
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
std::optional<InputBlock> InputBlock::getBlock(std::string_view name) const {
  // note: by copy!
  const auto block = std::find(m_blocks.crbegin(), m_blocks.crend(), name);
  if (block == m_blocks.crend())
    return {};
  return *block;
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
    // const std::vector<std::pair<std::string, std::string>> &list, bool print)
    const std::vector<std::string> &list, bool print) const {
  // Check each option NOT each sub block!
  // For each input option stored, see if it is allowed
  // "allowde" means appears in list
  bool all_ok = true;
  for (const auto &[key, value] : m_options) {
    // const auto is_optionQ = [&](const auto &l) { return key == l.first; };
    const auto is_optionQ = [&](const auto &l) { return key == l; };
    const auto bad_option =
        !std::any_of(list.cbegin(), list.cend(), is_optionQ);
    auto help = (key == "Help" || key == "help") ? true : false;
    if (bad_option && !help) {
      all_ok = false;
      std::cout << "\n⚠️  WARNING: Unclear input option in " << m_name
                << ": " << key << " = " << value << ";\n"
                << "Option may be ignored!\n"
                << "Check spelling (or update list of options)\n";
    }
  }

  for (const auto &block : m_blocks) {
    // (void)value; // not needed
    // const auto is_blockQ = [&](const auto &b) { return block == b.first; };
    const auto is_blockQ = [&](const auto &b) { return block == b; };
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
      // std::cout << "  " << s.first << ";  // " << s.second << "\n";
      std::cout << "  " << s << ";  // " << s << "\n";
    });
    std::cout << "}\n\n";
  }
  return all_ok;
}

//! Check one of the sub-blocks
bool InputBlock::check(
    std::initializer_list<std::string> blocks,
    // const std::vector<std::pair<std::string, std::string>> &list,
    const std::vector<std::string> &list, bool print) const {
  // Find key in nested blocks
  const InputBlock *pB = this;
  for (const auto &block : blocks) {
    pB = pB->getBlock_cptr(block);
    if (pB == nullptr) {
      // Did not fund nested block... may be fine
      return true;
    }
  }
  return pB->checkBlock(list, print);
}

//******************************************************************************
void InputBlock::add_blocks_from_string(std::string_view string) {

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
      auto next_start = start + 1;
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
      block.add_blocks_from_string(string.substr(start, end - start));
    }

    start = end + 1;
  }
  // Merge duplicated blocks.
  consolidate();
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

} // namespace IO
