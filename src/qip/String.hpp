#pragma once
#include "fmt/color.hpp"
#include <algorithm>
#include <cctype>
#include <cctype> //char from string
#include <cstdarg>
#include <functional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace qip {

//==============================================================================

/*!
  @brief Returns a formatted string using printf-style format specifiers.

  @details Maximum string length is 256 characters.
  Use the size overload for longer strings.
*/
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

//! Overload of fstring with explicit buffer size (maximum string length).
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

//==============================================================================

/*!
  @brief Compares s1 against pattern s2, where s2 may contain one wildcard '*'
  that matches any substring.
*/
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

//==============================================================================

//! Conversion of a single character to lowercase.
inline char tolower(char ch) {
  // https://en.cppreference.com/w/cpp/string/byte/tolower
  return static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
}

//! Returns a lowercase copy of the string.
inline std::string tolower(std::string t_string) {
  for (auto &c : t_string) {
    c = qip::tolower(c);
  }
  return t_string;
}

//==============================================================================

//! Returns true if the_string contains sub_string.
inline bool contains(std::string_view the_string, std::string_view sub_string) {
  return the_string.find(sub_string) != std::string::npos;
}

//! Returns true if the_string contains sub_string (case insensitive).
inline bool ci_contains(const std::string &the_string,
                        const std::string &sub_string) {
  return tolower(the_string).find(tolower(sub_string)) != std::string::npos;
}

//! Returns true if the_string contains any of the sub_strings.
inline bool contains(const std::string &the_string,
                     const std::vector<std::string> &sub_strings) {
  for (const auto &substr : sub_strings) {
    if (contains(the_string, substr))
      return true;
  }
  return false;
}

//! Returns true if the_string contains any of the sub_strings (case insensitive).
inline bool ci_contains(const std::string &the_string,
                        const std::vector<std::string> &sub_strings) {
  for (const auto &substr : sub_strings) {
    if (ci_contains(the_string, substr))
      return true;
  }
  return false;
}

//==============================================================================

//! Case-insensitive string comparison; equivalent to tolower(s1) == tolower(s2).
inline bool ci_compare(std::string_view s1, std::string_view s2) {
  return std::equal(
    s1.cbegin(), s1.cend(), s2.cbegin(), s2.cend(),
    [](char c1, char c2) { return qip::tolower(c1) == qip::tolower(c2); });
}

/*!
  @brief Case-insensitive version of @ref wildcard_compare.

  @details Compares s1 against pattern s2, where s2 may contain one wildcard
  '*' that matches any substring.
*/
inline bool ci_wc_compare(std::string_view s1, std::string_view s2) {
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

//==============================================================================

//! Returns the Levenshtein edit distance between strings a and b.
inline auto Levenstein(std::string_view a, std::string_view b) {
  // https://en.wikipedia.org/wiki/Levenshtein_distance
  // https://stackoverflow.com/a/70237726/8446770
  std::vector<size_t> d_t((a.size() + 1) * (b.size() + 1), size_t(-1));
  auto d = [&](size_t ia, size_t ib) -> size_t & {
    return d_t[ia * (b.size() + 1) + ib];
  };
  std::function<size_t(size_t, size_t)> LevensteinInt =
    [&](size_t ia, size_t ib) -> size_t {
    if (d(ia, ib) != size_t(-1))
      return d(ia, ib);
    size_t dist = 0;
    if (ib >= b.size())
      dist = a.size() - ia;
    else if (ia >= a.size())
      dist = b.size() - ib;
    else if (a[ia] == b[ib])
      dist = LevensteinInt(ia + 1, ib + 1);
    else
      dist = 1 + std::min(std::min(LevensteinInt(ia, ib + 1),
                                   LevensteinInt(ia + 1, ib)),
                          LevensteinInt(ia + 1, ib + 1));
    d(ia, ib) = dist;
    return dist;
  };
  return LevensteinInt(0, 0);
}

//! Case-insensitive version of @ref Levenstein.
inline auto ci_Levenstein(std::string_view a, std::string_view b) {
  std::vector<size_t> d_t((a.size() + 1) * (b.size() + 1), size_t(-1));
  auto d = [&](size_t ia, size_t ib) -> size_t & {
    return d_t[ia * (b.size() + 1) + ib];
  };
  std::function<size_t(size_t, size_t)> LevensteinInt =
    [&](size_t ia, size_t ib) -> size_t {
    if (d(ia, ib) != size_t(-1))
      return d(ia, ib);
    size_t dist = 0;
    if (ib >= b.size())
      dist = a.size() - ia;
    else if (ia >= a.size())
      dist = b.size() - ib;
    else if (qip::tolower(a[ia]) == qip::tolower(b[ib]))
      dist = LevensteinInt(ia + 1, ib + 1);
    else
      dist = 1 + std::min(std::min(LevensteinInt(ia, ib + 1),
                                   LevensteinInt(ia + 1, ib)),
                          LevensteinInt(ia + 1, ib + 1));
    d(ia, ib) = dist;
    return dist;
  };
  return LevensteinInt(0, 0);
}

//! Returns an iterator to the closest match to test_string in list,
//! using @ref Levenstein distance.
inline auto closest_match(std::string_view test_string,
                          const std::vector<std::string> &list) {
  auto compare = [&test_string](const auto &s1, const auto &s2) {
    return qip::Levenstein(s1, test_string) < qip::Levenstein(s2, test_string);
  };
  return std::min_element(list.cbegin(), list.cend(), compare);
}

//! Returns the closest match (case insensitive) to test_string in list,
//! using @ref ci_Levenstein distance.
inline std::string ci_closest_match(const std::string_view test_string,
                                    const std::vector<std::string> &list) {
  auto compare = [&test_string](const auto &s1, const auto &s2) {
    return qip::ci_Levenstein(s1, test_string) <
           qip::ci_Levenstein(s2, test_string);
  };
  using namespace std::string_literals;
  return list.empty() ? ""s :
                        *std::min_element(list.cbegin(), list.cend(), compare);
}

//==============================================================================

/*!
  @brief Returns true if the string represents an integer.

  @details
  Accepts an optional leading '+' or '-'. e.g., "16" and "-12" return true;
  "12x" and "12.5" return false.
*/
inline bool string_is_integer(std::string_view s) {
  return !s.empty() &&
         // checks if all non-leading characters are integer digits
         std::find_if(s.cbegin() + 1, s.cend(),
                      [](auto c) { return !std::isdigit(c); }) == s.end() &&
         // checks if leading character is one of: digit, '+', or '-'
         (std::isdigit(s[0]) || ((s[0] == '-' || s[0] == '+') && s.size() > 1));
}

//==============================================================================

//! Splits a string by delimiter into a vector of substrings.
inline std::vector<std::string> split(const std::string &s, char delim = ' ') {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tmp;
  while (getline(ss, tmp, delim)) {
    out.push_back(tmp);
  }
  return out;
}

//! Concatenates a vector of strings into one, with an optional delimiter.
inline std::string concat(const std::vector<std::string> &v,
                          const std::string &delim = "") {
  std::string out;
  for (std::size_t i = 0; i < v.size(); ++i) {
    out += v[i];
    if (i != v.size() - 1)
      out += delim;
  }
  return out;
}

//==============================================================================

/*!
  @brief Word-wraps input at column at, optionally prefixing each line.

  @details Does not split words unless unavoidable.
*/
inline std::string wrap(const std::string &input, std::size_t at = 80,
                        const std::string &prefix = "") {
  std::string output;
  const auto length = at - prefix.size();
  std::size_t ipos = 0;
  std::size_t fpos = length;
  while (ipos < input.length()) {
    if (!output.empty())
      output += '\n';

    auto temp_pos_nl = input.find('\n', ipos);
    if (temp_pos_nl > ipos && temp_pos_nl < fpos &&
        temp_pos_nl != std::string::npos) {
      output += prefix + input.substr(ipos, temp_pos_nl - ipos);
      ipos = temp_pos_nl + 1;
      fpos = ipos + length;
      continue;
    }

    if (fpos >= input.length()) {
      output += prefix + input.substr(ipos, fpos - ipos);
      break;
    }

    auto temp_pos = input.rfind(' ', fpos);
    if (temp_pos <= ipos || temp_pos == std::string::npos) {
      output += prefix + input.substr(ipos, fpos - ipos);
      ipos = fpos;
      fpos = ipos + length;
    } else {
      output += prefix + input.substr(ipos, temp_pos - ipos);
      ipos = temp_pos + 1;
      fpos = ipos + length;
    }
  }
  return output;
}

//==============================================================================

//! Converts an integer to a Roman numeral string. Assumes |a| <= 3999.
inline std::string int_to_roman(int a) {
  if (a < 0)
    return "-" + int_to_roman(-a);
  if (a > 3999)
    return std::to_string(a);
  static const std::string M[] = {"", "M", "MM", "MMM"};
  static const std::string C[] = {"",  "C",  "CC",  "CCC",  "CD",
                                  "D", "DC", "DCC", "DCCC", "CM"};
  static const std::string X[] = {"",  "X",  "XX",  "XXX",  "XL",
                                  "L", "LX", "LXX", "LXXX", "XC"};
  static const std::string I[] = {"",  "I",  "II",  "III",  "IV",
                                  "V", "VI", "VII", "VIII", "IX"};
  return M[a / 1000] + C[(a % 1000) / 100] + X[(a % 100) / 10] + I[(a % 10)];
}

} // namespace qip
