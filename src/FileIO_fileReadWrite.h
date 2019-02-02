#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
/*
Functions to help read-write data to binary files.
NOTE: Only works for plain data. Have to de/re-construct vectors etc. later!
ALSO: there are no safety checks here! If format is wrong, leads to undefined
behaviour (most likely, a crash)
*/

namespace FileIO {

//******************************************************************************
inline std::vector<std::string> readInputFile(const std::string &fname) {
  std::vector<std::string> entry_list;
  std::ifstream file(fname);
  std::string line;
  while (getline(file, line)) {
    std::stringstream ss(line);
    std::string entry;
    while (ss >> entry) {
      if (entry.at(0) == '!' || entry.at(0) == '#')
        break;
      else
        entry_list.push_back(entry);
    }
  }
  return entry_list;
}

//******************************************************************************
enum ROW { read, write };

inline void open_binary(std::fstream &stream, const std::string &fname,
                        ROW row) {
  switch (row) {
  case write:
    stream.open(fname, std::ios_base::out | std::ios_base::binary);
    break;
  case read:
    stream.open(fname, std::ios_base::in | std::ios_base::binary);
    break;
  default:
    std::cout << "\nFAIL 16 in FileIO\n";
  }
}

template <typename T> void binary_rw(std::fstream &stream, T &value, ROW row) {
  switch (row) {
  case write:
    stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
    break;
  case read:
    stream.read(reinterpret_cast<char *>(&value), sizeof(T));
    break;
  default:
    std::cout << "\nFAIL 32 in FileIO\n";
  }
}

inline void binary_str_rw(std::fstream &stream, std::string &value, ROW row) {
  if (row == write) {
    size_t temp_len = value.length();
    stream.write(reinterpret_cast<const char *>(&temp_len), sizeof(size_t));
    stream.write(value.c_str(), value.length());
  } else if (row == read) {
    size_t temp_len;
    stream.read(reinterpret_cast<char *>(&temp_len), sizeof(size_t));
    char *tvalue = new char[temp_len + 1];
    stream.read(tvalue, temp_len);
    tvalue[temp_len] = '\0'; // null 'end of string' character
    value = tvalue;
    delete[] tvalue;
  } else {
    std::cout << "\nFAIL 55 in FileIO\n";
  }
}

} // namespace FileIO
