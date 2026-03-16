#pragma once
#include <string>

// #include this instead of <omp.h>
// This allows ompenMP code to compile, even if

#if defined(_OPENMP)
#include <omp.h>
namespace qip {
constexpr bool use_omp = true;
} // namespace qip
#else
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
namespace qip {
constexpr bool use_omp = false;
} // namespace qip
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

namespace qip {
inline std::string omp_details() {
  return use_omp ? "Using OpenMP with " +
                       std::to_string(omp_get_max_threads()) + " threads." :
                   "Single-threaded.";
}
} // namespace qip