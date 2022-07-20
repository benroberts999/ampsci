#pragma once

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
