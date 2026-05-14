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
  @brief File read/write utilities: text parsing and binary I/O.
  @details
  Low-level helpers for reading and writing text and binary files. Includes
  functions for parsing space- or line-separated input files (with comment
  stripping), reading (x, y) data, and variadic binary read/write for PoD types,
  vectors, and strings.

  @warning Many functions in this namespace are old and some are obsolete.
  Use sparingly and with caution -- prefer higher-level IO interfaces where
  available.
*/
namespace IO::FRW {

//------------------------------------------------------------------------------
//! Reads (x, y) pairs from a file; returns a pair of vectors {xs, ys}.
//! Lines beginning with '#', '!', or '//' are skipped.
inline std::pair<std::vector<double>, std::vector<double>>
readFile_xy_PoV(const std::string &fname) {
  std::pair<std::vector<double>, std::vector<double>> out_list;
  std::ifstream file(fname);
  std::string line = "";
  while (getline(file, line) && (file.is_open())) { // redundant?
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
//! Removes C-style block comments (/* ... */) from a string in-place.
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
//! Strips comments (#, !, //, /* */), spaces, tabs, and quote characters.
//! Lines are squashed together; semicolons are preserved as delimiters.
inline std::string removeCommentsAndSpaces(const std::string &input) {
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
enum RoW { read, write };

// XXX Wrap these into a class. AND make explicit which types are/aren't
// allowed!

//! Opens a binary fstream for reading or writing according to @p row.
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

//! Returns true if the file at @p fileName exists and can be opened.
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
