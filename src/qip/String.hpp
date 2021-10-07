#pragma once
#include <algorithm>
#include <string>
#include <string_view>

namespace qip {

//! Compares two strings, s1 and s2. s2 may contain ONE wildcard ('*') which
//! will match anything
bool wildcard_compare(std::string_view s1, std::string_view s2) {
  // look for wildcard:
  const auto wc = std::find(s2.cbegin(), s2.cend(), '*');
  if (wc == s2.cend())
    return s1 == s2;

  const auto pos_wc = std::size_t(std::distance(s2.cbegin(), wc));

  const auto s1_front = s1.substr(0, pos_wc);
  const auto s2_front = s2.substr(0, pos_wc);

  // number of characters following the '*'
  const auto len_back = std::size_t(std::distance(wc + 1, s2.cend()));

  const auto pos_1_back = s1.length() > len_back ? s1.length() - len_back : 0;
  const auto s1_back = s1.substr(pos_1_back, std::string::npos);
  const auto s2_back = s2.substr(pos_wc + 1, std::string::npos);

  return s1_front == s2_front && s1_back == s2_back;
}

} // namespace qip
