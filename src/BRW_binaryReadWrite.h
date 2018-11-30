#pragma once
#include <fstream>
#include <string>
/*
Functions to help read-write data to binary files.
NOTE: Only works for plain data. Have to de/re-construct vectors etc. later!
ALSO: there are no safety checks here! If format is wrong, leads to undefined
behaviour (most likely, a crash)
*/

namespace BRW{

enum ROW { read, write };

void open_binary(std::fstream& stream, const std::string &fname, ROW row);

template<typename T>
void binary_rw(std::fstream& stream, T& value, ROW row);

void binary_str_rw(std::fstream& stream, std::string& value, ROW row);

}
