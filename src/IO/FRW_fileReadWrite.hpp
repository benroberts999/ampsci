#pragma once
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
/*!
Functions to help read-write data to files, including binary files.
Very old and poorly written..but works.
*/
namespace IO::FRW {

//==============================================================================
//! Reads a text file line-by-line into a pair of vectors.
//! Ignores comments (!, #, //)
inline std::pair<std::vector<double>, std::vector<double>>
readFile_xy_PoV(const std::string &fname) {
  std::pair<std::vector<double>, std::vector<double>> out_list;
  std::ifstream file(fname);
  std::string line = "";
  while (getline(file, line) && (file.is_open())) {
    if (line == "")
      continue;
    if (line.at(0) == '!' || line.at(0) == '#')
      continue;
    if (line.size() >= 2 && line.at(0) == '/' && line.at(1) == '/')
      continue;
    std::stringstream ss(line);
    double x, y;
    ss >> x >> y;
    out_list.first.push_back(x);
    out_list.second.push_back(y);
  }
  return out_list;
}

//==============================================================================
//! Removes block comments from a string
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
//==============================================================================
//! Removes comments, spaces/tabs, quotes and newlines from a string
//! (comments: !, #, //, and block /**/ comments)
inline std::string removeCommentsAndSpaces(const std::string &input) {
  // Note: also squashes lines, except for semi-colons
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
    }
  }

  removeBlockComments(lines);

  // remove spaces
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == ' '; }),
              lines.end());
  // remove tabs
  lines.erase(std::remove_if(lines.begin(), lines.end(),
                             [](unsigned char x) { return x == '\t'; }),
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

//==============================================================================
//! Splits an input string into a vector of pairs of strings; split by '{}'
//! For example, ABC{DE} -> [ABC,DE]
inline std::vector<std::pair<std::string, std::string>>
splitInput_byBraces(const std::string &input) {

  std::vector<std::pair<std::string, std::string>> output;

  auto lines = removeCommentsAndSpaces(input);

  auto find_close = [&lines](auto open) {
    auto next_open = lines.find('{', open + 1);
    auto next_close = lines.find('}', open + 1);
    auto next = std::min(next_open, next_close);
    if (next == next_close || next == std::string::npos)
      return next;
    int depth = 1;
    while (depth > 0) {
      next_open = lines.find('{', next + 1);
      next_close = lines.find('}', next + 1);
      next = std::min(next_open, next_close);
      if (next == next_close)
        --depth;
      if (next == next_open)
        ++depth;
    }
    next = lines.find('}', next + 1);
    return next;
  };

  std::size_t previous_end = 0;
  while (true) {
    auto beg = lines.find('{', previous_end);
    auto end = find_close(beg);
    if (beg == std::string::npos)
      break;
    if (end == std::string::npos) {
      std::cerr << "\nFAIL 114 in FRW: Bad file format (missing '}'?)\n";
      break;
    }
    auto identifier = lines.substr(previous_end, beg - previous_end);
    auto options = lines.substr(beg + 1, end - beg - 1);
    output.push_back(std::make_pair(identifier, options));
    previous_end = end + 1;
  }

  return output;
}

//==============================================================================
//! Splits a string into a vector of strings, by semicolon... with exceptions..
inline std::vector<std::string> splitInput_bySemiColon(const std::string &input)
// ...
{
  std::vector<std::string> entry_list;

  auto lines = removeCommentsAndSpaces(input);

  // BUT, treat anything inside [] brackets as single

  // after removing comments, break up by ';'

  // dumb hack to ingore first braket '['..
  std::size_t start = (input.size() > 0 && input[0] == '{') ? 1 : 0;

  while (true) {
    // Check if we find a ';' or an open braket '[' first:
    const auto bkt = input.find('{', start);
    const auto semic = input.find(';', start);
    const auto bracketQ = bkt < semic;

    // If no braket, 'end' of string is semicolon ';'.
    // If we do have brakets, end is "];"
    const auto end = bracketQ ? input.find("};", start) + 1 : semic;
    if (end == std::string::npos)
      break;
    entry_list.push_back(input.substr(start, end - start));
    start = end + 1;
  }

  return entry_list;
}

//==============================================================================
enum RoW { read, write };

// XXX Wrap these into a class. AND make explicit which types are/aren't
// allowed!

inline void open_binary(std::fstream &stream, const std::string &fname,
                        RoW row) {
  switch (row) {
  case write:
    stream.open(fname, std::ios_base::out | std::ios_base::binary);
    break;
  case read:
    stream.open(fname, std::ios_base::in | std::ios_base::binary);
    break;
  default:
    std::cout << "\nFAIL 16 in FRW\n";
  }
}

inline bool file_exists(const std::string &fileName) {
  if (fileName == "")
    return false;
  std::ifstream infile(fileName);
  return infile.good();
}

//! Function (variadic): reads/writes data from/to binary file. Works for
//! trivial (PoD) types, and std::string only (but not checked)
//! Overload for vectors
template <typename T, typename... Types>
void rw_binary(std::fstream &stream, RoW row, std::vector<T> &value,
               Types &...values) {
  binary_rw_vec(stream, value, row);
  if constexpr (sizeof...(values) != 0)
    rw_binary(stream, row, values...);
}

//! Function (variadic): reads/writes data from/to binary file. Works for
//! trivial (PoD) types, and std::string only (but not checked)
template <typename T, typename... Types>
void rw_binary(std::fstream &stream, RoW row, T &value, Types &...values) {
  if constexpr (std::is_same_v<T, std::string>)
    binary_str_rw(stream, value, row);
  else
    binary_rw(stream, value, row);
  if constexpr (sizeof...(values) != 0)
    rw_binary(stream, row, values...);
}

// make "private" behind 'helper' namespace: used in old code still
template <typename T>
void binary_rw(std::fstream &stream, T &value, RoW row) {
  switch (row) {
  case write:
    stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
    break;
  case read:
    stream.read(reinterpret_cast<char *>(&value), sizeof(T));
    break;
  default:
    std::cout << "\nFAIL 32 in FRW\n";
  }
}

// For vector. {works with Vector of vectors also}
template <typename T>
void binary_rw_vec(std::fstream &stream, std::vector<T> &value, RoW row) {
  std::size_t size = value.size();
  binary_rw(stream, size, row);
  if (row == read)
    value.resize(size);
  for (auto &x : value) {
    rw_binary(stream, row, x);
  }
}

// make "private" behind 'helper' namespace: used in old code still
inline void binary_str_rw(std::fstream &stream, std::string &value, RoW row) {
  if (row == write) {
    std::size_t temp_len = value.length();
    stream.write(reinterpret_cast<const char *>(&temp_len),
                 sizeof(std::size_t));
    stream.write(value.c_str(), long(value.length()));
  } else if (row == read) {
    std::size_t temp_len;
    stream.read(reinterpret_cast<char *>(&temp_len), sizeof(std::size_t));
    char *tvalue = new char[temp_len + 1];
    stream.read(tvalue, long(temp_len));
    tvalue[temp_len] = '\0'; // null 'end of string' character
    value = tvalue;
    delete[] tvalue;
  } else {
    std::cout << "\nFAIL 55 in FRW\n";
  }
}

} // namespace IO::FRW
