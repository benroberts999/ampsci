#pragma once
/*! @file
  @brief Include instead of `<omp.h>` to allow compilation with or without OpenMP.

  @details
  If OpenMP is not available, stub macros are provided so `#pragma omp`
  directives are silently ignored, and the common `omp_*thread` functions return
  0 or 1.
*/
#include <string>

#if defined(_OPENMP)
#include <omp.h>
namespace qip {
//! True if compiled with OpenMP support, false otherwise.
constexpr bool use_omp = true;
} // namespace qip
#else
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace qip {
//! True if compiled with OpenMP support, false otherwise.
constexpr bool use_omp = false;
} // namespace qip
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#endif

namespace qip {
//! Returns a short string describing the threading status, e.g.
//! "Using OpenMP with 8 threads." or "Single-threaded."
inline std::string omp_details() {
  return use_omp ? "Using OpenMP with " +
                     std::to_string(omp_get_max_threads()) + " threads." :
                   "Single-threaded.";
}
} // namespace qip
