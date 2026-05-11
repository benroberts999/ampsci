#pragma once
#include "qip/omp.hpp"
#include <atomic>
#include <cstdio>
#include <iostream>
// Note: Uses POSIX isatty(). Not portable to Windows without WSL/MinGW.
#include <unistd.h>

namespace qip {

/*! @brief Basic progress bar. Prints new line if (and only if) i==(max-1)
    @details 
    - does not work well in Parellel regions. Use @ref ProgressBar in that case
    - Prints directly to cout - creates a mess if piped to a text file.
      Use @ref ProgressBar if that will be an issue
*/
inline void progbar(int i, int max, int length = 50) {
  const int len = (length - 1);
  const int current = int(len * double(i) / double(max - 1));
  std::cout << "[";
  for (auto j = 0; j < current; ++j) {
    std::cout << "=";
  }
  for (auto j = current; j < len; ++j) {
    std::cout << " ";
  }
  std::cout << "]   \r" << std::flush;
  if (i == max - 1)
    std::cout << "\n";
}

//==============================================================================
/*!
  @brief Thread-safe progress bar for OpenMP parallel loops.

  @details
  Displays a progress bar with percentage. The progress counter uses
  std::atomic for thread-safe updates.
  Each call to update() increments the counter and prints the bar.
  The output is serialised via critical section to prevent garbled output
  from simultaneous writes.

  @warning This adds overhead. 
  Prefer not to use if each OMP task is extremely small, 
  since overhead may become noticable. 
  If each task is large, overhead negligable.

  - If print set to false on construction, does nothing 
  (does not print, does not track progress). 
  Just a simple way of run-time turning off.

  @note
  When stdout is not a normal TTY terminal (i.e., when output is piped), 
  intermediate updates are not printed. The bar will start at 0%, and
  then only jump to 100% when job is completed.

  The @p length template parameter specifies the total character width of the
  output, including '[', the bar content, ']', and the percentage display.
  The default is 55.

  @note If the loop exits before the final iteration (i.e., early `break;`), 
  then final the newline `\n` will not be printed. Output may be messy.
  Cannot break like this in OpenMP loop anyway.

  For non-parallel loops, the simpler @ref qip::progbar() function should work fine.

  Typical usage in a parallel loop: 
  @code
  qip::ProgressBar bar(n_iterations);
  #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < n_iterations; ++i) {
      // ... work ...
      bar.update();
    }
  @endcode

  - Put .update(); _after_ work, to avoid early "100% done" reporting

  @tparam Length Total character width of the output (default 55).

  @note Thread-safe; safe to call from multiple threads simultaneously.
  But may add overhead.
*/
template <std::size_t Length = 55>
class ProgressBar {
private:
  int m_max;
  bool m_print;
  // Atomicly tracks progress
  std::atomic<int> m_progress{0};
  // bar chars printed last time; skip if unchanged
  std::atomic<int> m_last_fill{-1};

  static constexpr std::size_t bar_width = (Length >= 7) ? (Length - 7) : 1;

  void print_prog_bar(int progress) {
    const bool is_initial = (progress == 0);
    const bool is_final = (progress == m_max);

    // when piped, skip intermediate updates to avoid cluttering captured output
    const bool is_tty = isatty(fileno(stdout));
    if (!is_tty && !is_final && !is_initial)
      return;

    // build bar into stack buffer (char array)
    const auto current = int(bar_width * double(progress) / double(m_max));
    char buf[Length + 2];
    char *p = buf;
    *p++ = '[';
    for (int j = 0; j < current; ++j)
      *p++ = '=';
    for (int j = current; j < (int)bar_width; ++j)
      *p++ = ' ';
    *p++ = ']';
    *p++ = ' ';
    p += snprintf(p, 10, "%d%%", int(100.0 * double(progress) / double(m_max)));
    *p++ = is_final ? '\n' : '\r';
    *p = '\0';

    // serialise writes to std::out
#pragma omp critical
    {
      fputs(buf, stdout);
      fflush(stdout);
    }
    m_last_fill.store(current, std::memory_order_relaxed);
  }

public:
  /*!
    @brief Construct progress bar for @p max iterations.
    @param max Total number of iterations (denominator for percentage).
    @param print Runtime switch; if false, class does nothing.
  */
  ProgressBar(int max, bool print) : m_max(max), m_print(print) {
    if (m_print)
      print_prog_bar(0);
  }

  //! Atomically increment progress counter and print updated bar.
  void update() {
    if (!m_print)
      return;
    m_progress++;
    const int progress = m_progress.load();
    if (progress > m_max)
      return;
    const bool is_final = (progress == m_max);

    // skip print if bar fill is unchanged (avoids critical section overhead)
    const int current = int(bar_width * double(progress) / double(m_max));
    if (!is_final && current <= m_last_fill.load(std::memory_order_relaxed))
      return;

    print_prog_bar(progress);
  }
};

} // namespace qip
