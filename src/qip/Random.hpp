#pragma once
#include <iostream>
#include <random>
#include <string>

namespace qip {

inline std::string random_string(std::size_t length) {
  std::string out;
  out.reserve(length);
  const std::string chars{
      "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"};
  std::uniform_int_distribution<std::size_t> rint(0, chars.size() - 1);
  std::random_device rd;
  std::mt19937 gen(rd());
  for (std::size_t i = 0; i < length; ++i) {
    out.push_back(chars[rint(gen)]);
  }
  return out;
}

} // namespace qip