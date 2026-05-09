#pragma once
#include "qip/omp.hpp"
#include <atomic>
#include <cstdio>
#include <iostream>
// Note: Uses POSIX functions (isatty, fopen) and /dev/tty.
// Not portable to Windows without WSL/MinGW.
#include <unistd.h>

namespace qip {

/*! @brief Basic progress bar. Prints new line if (and only if) i==(max-1)
    @details 
    - does not work well in Parellel regions. Use @ref ProgressBar in that case
    - Prints directly to cout - creates a mess if piped to a text file.
      Use @ref ProgressBar is that will be an issue
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
  @brief Thread-safe template progress bar for OpenMP parallel loops.
  @details
  Displays a text-based progress bar with percentage. The progress counter uses
  std::atomic for thread-safe updates, eliminating need for explicit
  synchronisation pragmas. Each call to update() increments the counter and prints
  the bar. The output is serialised via critical section to prevent garbled output
  from simultaneous writes.

  @note This adds overhead. Prefer not to use if each OMP task is extremely small,
  since overhead may become noticable. If each task is large, overhead negligable.

  Only prints when stdout is a TTY (interactive terminal). When piped to a file
  or non-interactive stream, update() attemps to write to terminal via `/dev/tty`.
  If this is not available, it returns without printing, avoiding output
  corruption from carriage returns.

  The @p length template parameter specifies the total character width of the
  output, including '[', the bar content, ']', and the percentage display.
  The default is 55.

  @note If the loop exits before the final iteration (i.e., early `break;`), 
  then final the newline `\n` will not be printed. Output may be messy.
  Cannot break like this in OpenMP loop anyway.

  For non-parallel loops, the simpler @ref qipprogbar::() function should work fine.

  Typical usage in a parallel loop: 
  @code
  qip::ProgressBar bar(n_iterations);
  #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < n_iterations; ++i) {
      bar.update();  
      // ... work ...
    }
  @endcode

  @tparam Length Total character width of the output (default 55).

  @note Thread-safe; safe to call from multiple threads simultaneously.
*/
template <std::size_t Length = 55>
class ProgressBar {
private:
  std::atomic<int> m_progress{0};
  int m_max;

  static constexpr std::size_t bar_width = (Length >= 7) ? (Length - 7) : 1;

public:
  //!
  /*!
    @brief Construct progress bar for @p max iterations.
    @param max Total number of iterations (denominator for percentage).
  */
  ProgressBar(int max) : m_max(max) {}

  //! Atomically increment progress counter and print updated bar.
  void update() {

    // Write to stdout if it's a TTY, otherwise try /dev/tty (terminal).
    // This allows progress bar to display when stdout is piped to a file.
    // If can't use /dev/tty, do nothing
    FILE *out = nullptr;
    if (isatty(fileno(stdout))) {
      out = stdout;
    } else {
      out = fopen("/dev/tty", "w");
      // If no TTY available; skip output.
      if (out == nullptr)
        return;
    }

    m_progress++;
    const int progress = m_progress.load();
    const int current = int(bar_width * double(progress) / double(m_max));

    char buf[Length + 2];
    char *p = buf;

    *p++ = '[';
    for (int j = 0; j < current; ++j) {
      *p++ = '=';
    }
    for (int j = current; j < (int)bar_width; ++j) {
      *p++ = ' ';
    }
    *p++ = ']';
    *p++ = ' ';

    const int pct = int(100.0 * double(progress) / double(m_max));
    const int pct_chars = snprintf(p, 10, "%d%%", pct);
    p += pct_chars;

    if (m_progress == m_max) {
      *p++ = '\n';
    } else {
      *p++ = '\r';
    }
    *p = '\0';

#pragma omp critical
    {
      fputs(buf, out);
      fflush(out);
      if (out != stdout)
        fclose(out);
    }
  }
};

} // namespace qip
