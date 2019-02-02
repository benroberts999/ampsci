#include "BRW_binaryReadWrite.h"
#include <fstream>
#include <iostream>
#include <string>
/*
Functions to help read-write data to binary files.
NOTE: Only works for plain data. Have to de/re-construct vectors etc. later!
ALSO: there are no safety checks here! If format is wrong, leads to undefined
behaviour (most likely, a crash)
*/

namespace BRW {

// Open the file for binary read/write
void open_binary(std::fstream &stream, const std::string &fname, ROW row) {
  switch (row) {
  case write:
    stream.open(fname, std::ios_base::out | std::ios_base::binary);
    break;
  case read:
    stream.open(fname, std::ios_base::in | std::ios_base::binary);
    break;
  default:
    std::cout << "\nFAIL 16 in BRW\n";
  }
}

// templates for int, double, float (can add others later)
template void binary_rw(std::fstream &stream, int &value, ROW row);
template void binary_rw(std::fstream &stream, double &value, ROW row);
template void binary_rw(std::fstream &stream, float &value, ROW row);

// Do the actual reading/writing
template <typename T> void binary_rw(std::fstream &stream, T &value, ROW row) {
  switch (row) {
  case write:
    stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
    break;
  case read:
    stream.read(reinterpret_cast<char *>(&value), sizeof(T));
    break;
  default:
    std::cout << "\nFAIL 32 in BRW\n";
  }
}

// For strings:
// void binary_str_rw(std::fstream& stream, std::string& value, ROW row);

// template<typename T>

void binary_str_rw(std::fstream &stream, std::string &value, ROW row) {
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
    std::cout << "\nFAIL 55 in BRW\n";
  }
}

} // namespace BRW
