#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
/*
Functions to help read-write data to binary files.
NOTE: Only works for plain data. Have to de/re-construct vectors etc. later!
ALSO: there are no safety checks here! If format is wrong, leads to undefined
behaviour (most likely, a crash)
*/

namespace FileIO {

//******************************************************************************
/*
Uses compile-time recursion to get access to elements of tuple.
Specifically, string-streams data from a string vector into tuple.
Works with a tuple of references, i.e., std::forward_as_tuple
Idea from:
https://stackoverflow.com/questions/1198260/iterate-over-tuple/23142715
*/
template <std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
stringstreamVectorIntoTuple(std::vector<std::string>, std::tuple<Tp...> &) {}

template <std::size_t I = 0, typename... Tp>
    inline typename std::enable_if <
    I<sizeof...(Tp), void>::type
    stringstreamVectorIntoTuple(std::vector<std::string> lst,
                                std::tuple<Tp...> &t) {
  if (I > lst.size())
    std::cerr << "\nFAIL 34 in FileIO: list shorter than tuple\n";
  std::stringstream(lst[I]) >> std::get<I>(t);
  stringstreamVectorIntoTuple<I + 1, Tp...>(lst, t);
}

//******************************************************************************
inline std::vector<std::string> readInputFile(const std::string &fname) {
  std::vector<std::string> entry_list;
  std::ifstream file(fname);
  std::string line;
  while (getline(file, line) && (file.is_open())) { // redundant?
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
template <typename... Tp>
void setInputParameters(std::string infile, std::tuple<Tp...> &tp) {
  auto input = readInputFile(infile);
  if (sizeof...(Tp) > input.size()) {
    // Note: for now, I allow a longer-than-needed input list.
    // This allows us to not have to comment out all the crap below
    // the input file..
    std::cerr << "\nFail 71 in FileIO: Wrong number of input parameters? "
              << "Reading from file: " << infile << ". Expected "
              << sizeof...(Tp) << " arguments"
              << ", but got " << input.size() << ".\n";
    return;
  }
  stringstreamVectorIntoTuple(input, tp);
}

//******************************************************************************
enum RoW { read, write };

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
    std::cout << "\nFAIL 16 in FileIO\n";
  }
}

template <typename T> void binary_rw(std::fstream &stream, T &value, RoW row) {
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

inline void binary_str_rw(std::fstream &stream, std::string &value, RoW row) {
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
