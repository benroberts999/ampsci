#pragma once
#include <algorithm>
#include <cctype>
#include <cstdarg>
#include <string>
#include <string_view>

namespace qip {

//******************************************************************************
//! Returns a formatted std::string, with formatting printf-like commands. Note:
//! maximum string lenth is 256 characters. If longer string required, use
//! provided overload
inline std::string fstring(const std::string format, ...) {
  constexpr std::size_t size = 256;
  std::string fmt_str;
  fmt_str.resize(size + 1); // allow for null

  // C-style variadic param-list, to call c function vsnprintf (varidatic
  // snprintf)
  va_list args;
  // note: format in va_start mist not be a reference type.. so copy the string?
  va_start(args, format);
  vsnprintf(&fmt_str[0], fmt_str.size(), format.c_str(), args);
  va_end(args);

  // resize string, remove part after the buffer (not needed)
  fmt_str.erase(std::find(fmt_str.begin(), fmt_str.end(), '\0'), fmt_str.end());

  return fmt_str;
}

//! Overload: size is maximum string length (buffer size).
inline std::string fstring(const std::size_t size, const std::string format,
                           ...) {
  // nb: cannot just call other overload, since using c-style variadic function
  // (I think?) - so a copy-paste re-implementation
  std::string fmt_str;
  fmt_str.resize(size + 1); // allow for null

  // C-style variadic param-list, to call c function vsnprintf (varidatic
  // snprintf)
  va_list args;
  // note: format in va_start mist not be a reference type.. so copy the string?
  va_start(args, format);
  vsnprintf(&fmt_str[0], fmt_str.size(), format.c_str(), args);
  va_end(args);

  // resize string, remove part after the buffer (not needed)
  fmt_str.erase(std::find(fmt_str.begin(), fmt_str.end(), '\0'), fmt_str.end());

  return fmt_str;
}

//******************************************************************************
//! Compares two strings, s1 and s2. s2 may contain ONE wildcard ('*') which
//! will match anything
inline bool wildcard_compare(std::string_view s1, std::string_view s2) {
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

//******************************************************************************
//! return static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
inline char tolower(char ch) {
  // https://en.cppreference.com/w/cpp/string/byte/tolower
  return static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
}

//******************************************************************************
//! Case insensitive string compare. Essentially: LowerCase(s1)==LowerCase(s2)
inline bool ci_compare(std::string_view s1, std::string_view s2) {
  return std::equal(
      s1.cbegin(), s1.cend(), s2.cbegin(), s2.cend(),
      [](char c1, char c2) { return qip::tolower(c1) == qip::tolower(c2); });
}

//! Compares two strings, s1 and s2. s2 may contain ONE wildcard ('*') which
//! will match anything. Case Insensitive version
inline bool ci_wildcard_compare(std::string_view s1, std::string_view s2) {
  // look for wildcard:
  const auto wc = std::find(s2.cbegin(), s2.cend(), '*');
  if (wc == s2.cend())
    return ci_compare(s1, s2);

  const auto pos_wc = std::size_t(std::distance(s2.cbegin(), wc));

  const auto s1_front = s1.substr(0, pos_wc);
  const auto s2_front = s2.substr(0, pos_wc);

  // number of characters following the '*'
  const auto len_back = std::size_t(std::distance(wc + 1, s2.cend()));

  const auto pos_1_back = s1.length() > len_back ? s1.length() - len_back : 0;
  const auto s1_back = s1.substr(pos_1_back, std::string::npos);
  const auto s2_back = s2.substr(pos_wc + 1, std::string::npos);

  return ci_compare(s1_front, s2_front) && ci_compare(s1_back, s2_back);
}

} // namespace qip
