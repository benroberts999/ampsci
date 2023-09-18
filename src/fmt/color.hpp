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

//! wrapper for text_style formatted fmt::print. Will only apply styling if
//! output is a terminal (stops ANSI characters written to file with >> or |tee)
template <typename... Args>
void styled_print(const fmt::text_style &ts, const Args &...args) {
  isatty(STDOUT_FILENO);
  fmt::print(enable_fmt_text_style ? ts : fmt::text_style(), args...);
}

inline void warning() { styled_print(fg(fmt::color::orange), "\nWARNING\n"); }
inline void error() { styled_print(fg(fmt::color::red), "\nERROR\n"); }
inline void fail() { styled_print(fg(fmt::color::red), "\nFAIL\n"); }

} // namespace fmt2