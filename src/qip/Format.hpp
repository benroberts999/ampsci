#pragma once
#include <algorithm>
#include <cstdarg>
#include <string>

namespace qip {

//! Returns a formatted std::string, with formatting printf-like commands. Note:
//! maximum string lenth is 80 characters. If longer string required, use
//! provided overload
std::string fstring(const std::string &format, ...) {
  static const std::size_t size = 80;
  std::string fmt_str;
  fmt_str.resize(size + 1); // allow for null

  // C-style variadic param-list, to call c function vsnprintf (varidatic
  // snprintf)
  va_list args;
  va_start(args, format);
  vsnprintf(&fmt_str[0], fmt_str.size(), format.c_str(), args);
  va_end(args);

  // resize string, remove part after the buffer (not needed)
  fmt_str.erase(std::find(fmt_str.begin(), fmt_str.end(), '\0'), fmt_str.end());

  return fmt_str;
}

//! Overload: size is maximum string length (buffer size).
std::string fstring(const std::size_t size, const std::string &format, ...) {
  // nb: cannot just call other overload, since using c-style variadic function
  // (I think?) - so a copy-paste re-implementation
  std::string fmt_str;
  fmt_str.resize(size + 1); // allow for null

  // C-style variadic param-list, to call c function vsnprintf (varidatic
  // snprintf)
  va_list args;
  va_start(args, format);
  vsnprintf(&fmt_str[0], fmt_str.size(), format.c_str(), args);
  va_end(args);

  // resize string, remove part after the buffer (not needed)
  fmt_str.erase(std::find(fmt_str.begin(), fmt_str.end(), '\0'), fmt_str.end());

  return fmt_str;
}

} // namespace qip
