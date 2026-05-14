#pragma once
#ifdef __GNUC__
// Avoid tons of warnings with root code
#pragma GCC system_header
#endif
#ifdef __clang__
// Avoid tons of warnings with root code
#pragma clang system_header
#endif
#include "fmt/format.hpp"
#include "fmt/include/color.h"
#include <unistd.h> // isatty

namespace fmt2 {

#if defined(__GNUC__) || defined(__clang__)
static const bool enable_fmt_text_style = isatty(STDOUT_FILENO);
#else
static const bool enable_fmt_text_style = false;
#endif

//! Wrapper for text_style formatted fmt::print.
/*!
  @brief Styled terminal print with automatic TTY detection.
  @details
  Calls fmt::print with the given text style, but only applies styling if
  stdout is a terminal. Suppresses ANSI escape codes when output is redirected
  (e.g., piped or written to file with `>>` or `| tee`).
*/
template <typename... Args>
void styled_print(const fmt::text_style &ts, const Args &...args) {
  isatty(STDOUT_FILENO);
  fmt::print(enable_fmt_text_style ? ts : fmt::text_style(), args...);
}

//! Prints orange "WARNING" label to stdout (styled if terminal).
inline void warning() { styled_print(fg(fmt::color::orange), "WARNING"); }
//! Prints red "ERROR" label to stdout (styled if terminal).
inline void error() { styled_print(fg(fmt::color::red), "ERROR"); }
//! Prints red "FAIL" label to stdout (styled if terminal).
inline void fail() { styled_print(fg(fmt::color::red), "FAIL"); }

} // namespace fmt2