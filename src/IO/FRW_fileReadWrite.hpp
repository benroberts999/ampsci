#pragma once
#include <algorithm>
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

namespace IO::FRW {

//******************************************************************************
// Uses compile-time recursion to get access to elements of tuple.
// Specifically, string-streams data from a string vector into tuple.
// Works with a tuple of references, i.e., std::forward_as_tuple
// Idea from:
// https://stackoverflow.com/questions/1198260/iterate-over-tuple/23142715
template <std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
stringstreamVectorIntoTuple(const std::vector<std::string> &,
                            std::tuple<Tp...> &) {}

template <std::size_t I = 0, typename... Tp>
    inline typename std::enable_if <
    I<sizeof...(Tp), void>::type
    stringstreamVectorIntoTuple(const std::vector<std::string> &lst,
                                std::tuple<Tp...> &t) {
  if (I > lst.size())
    std::cerr << "\nFAIL 34 in FRW: list shorter than tuple\n";
  std::stringstream(lst[I]) >> std::get<I>(t);
  stringstreamVectorIntoTuple<I + 1, Tp...>(lst, t);
}

//******************************************************************************
inline std::vector<std::string> readInputFile_byEntry(const std::string &fname)
// Reads each item (space separated) from a file into a string vector
// Any text after a '#' or '!' are treated as comments and ignored
// (can come anywhere in the line)
{
  std::vector<std::string> entry_list;
  std::ifstream file(fname);
  std::string line;
  while (getline(file, line) && (file.is_open())) { // redundant?
    std::stringstream ss(line);
    std::string entry;
    while (ss >> entry) {
      if (entry.at(0) == '!' || entry.at(0) == '#')
        break;
      if (entry.size() >= 2 && entry.at(0) == '/' && entry.at(1) == '/')
        break;
      entry_list.push_back(entry);
    }
  }
  return entry_list;
}

//******************************************************************************
inline std::vector<std::string> readInputFile_byLine(const std::string &fname)
// Reads each line from a file into a string vector
// Lines beginning with '!' or '#' are comments
{
  std::vector<std::string> entry_list;
  std::ifstream file(fname);
  std::string line;
  while (getline(file, line) && (file.is_open())) { // redundant?
    if (line == "")
      continue;
    if (line.at(0) == '!' || line.at(0) == '#')
      continue;
    if (line.size() >= 2 && line.at(0) == '/' && line.at(1) == '/')
      continue;
    entry_list.push_back(line);
  }
  return entry_list;
}

//******************************************************************************
inline std::vector<std::pair<double, double>>
readFile_xy_VoP(const std::string &fname)
// Reads each line from a file into a vector of {x,y} points
// Lines beginning with '!' or '#' are comments
// Could generalise this to any number of points?
{
  std::vector<std::pair<double, double>> out_list;
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
    out_list.emplace_back(x, y);
  }
  return out_list;
}
//------------------------------------------------------------------------------
inline std::pair<std::vector<double>, std::vector<double>>
readFile_xy_PoV(const std::string &fname)
// Reads each line from a file;
// Returns a pair of a vectors {{x},{y}}
// Lines beginning with '!' or '#' are comments
{
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

//------------------------------------------------------------------------------
inline void writeFile_xy(const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::string &fname) {
  //
  if (x.size() != y.size()) {
    std::cout << "Warning 139 in FRW: trying to write {x,y} vectors of "
                 "different lengths!\n";
  }
  std::ofstream file(fname);
  auto max = std::min(x.size(), y.size());
  for (auto i = 0ul; i < max; ++i) {
    file << x[i] << " " << y[i] << "\n";
  }
  //
}

//******************************************************************************
inline std::string readInputFile(const std::string &fname) {
  std::ifstream f(fname); // taking file as inputstream
  std::string str;
  if (f) {
    std::ostringstream ss;
    ss << f.rdbuf(); // reading data
    str = ss.str();
  }
  return str;
}

//******************************************************************************
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
//******************************************************************************
inline std::string removeCommentsAndSpaces(const std::string &input)
// Note: also squashes lines, except for semi-colons
{
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

//******************************************************************************
inline std::vector<std::pair<std::string, std::string>>
splitInput_byBraces(const std::string &input) {

  std::vector<std::pair<std::string, std::string>> output;

  auto lines = removeCommentsAndSpaces(input);

  std::size_t previous_end = 0;
  while (true) {
    auto beg = lines.find('{', previous_end);
    auto end = lines.find('}', beg);
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

//******************************************************************************
inline std::vector<std::string> splitInput_bySemiColon(const std::string &input)
// ...
{
  std::vector<std::string> entry_list;

  auto lines = removeCommentsAndSpaces(input);

  // after removing comments, break up by ';'
  std::stringstream stream2(lines);
  std::string entry;
  while (std::getline(stream2, entry, ';')) {
    entry_list.push_back(entry);
  }

  return entry_list;
}

//******************************************************************************
template <typename... Tp>
void setInputParameters(const std::string &infile, std::tuple<Tp...> &tp) {
  auto input = readInputFile_byEntry(infile);
  if (sizeof...(Tp) > input.size()) {
    // Note: for now, I allow a longer-than-needed input list.
    // This allows us to not have to comment out all the crap below
    // the input file..
    std::cerr << "\nFail 71 in FRW: Wrong number of input parameters? "
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
    std::cout << "\nFAIL 16 in FRW\n";
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
    std::cout << "\nFAIL 32 in FRW\n";
  }
}

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
