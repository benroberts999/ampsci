#pragma once
#include "qip/omp.hpp"
#include <atomic>
#include <cstdio>
#include <iostream>
#include <string>
#include <string_view>
// Note: Uses POSIX isatty(). Not portable to Windows without WSL/MinGW.
#include <unistd.h>

namespace qip {

/*!
  @brief Live status line for iterative loops; overwrites itself on a TTY.
  @details
  Construct before the loop, call operator() each iteration with the current
  status message, then call done() (optionally with a trailing annotation) when
  the loop finishes. The destructor calls done() automatically if it has not
  been called explicitly, so early returns and exceptions are handled cleanly.

  Behaviour depends on whether stdout is a TTY:
  - TTY: each operator() call prints `\r{header}{msg}` and flushes,
    overwriting the current line. done() prints `\r{header}{last_msg}{post}\n`.
  - Non-TTY: the header is printed on construction; operator() calls are
    silent (but the message is buffered). done() prints `{last_msg}{post}\n`.
    This gives one clean output line with no intermediate churn.

  If @p active is false all methods are no-ops (runtime print toggle).

  Typical usage:
  @code
  qip::LiveMessage status("TDHF E1 (w=1.23): ", print);
  for (int it = 0; it < max_its; ++it) {
    // ... work ...
    status(fmt::format("{:2d} {:.1e} [{}]", it, eps, worst));
    if (converged) break;
  }
  status.done("  ***");  // optional trailing annotation; destructor calls done() otherwise
  @endcode

  @param header  Fixed prefix, always printed (on construction for non-TTY,
                 or on every line for TTY).
  @param active  If false, all methods are no-ops. Default true.
*/
class LiveMessage {
public:
  explicit LiveMessage(std::string_view header, bool active = true)
    : m_header(header), m_is_tty(isatty(fileno(stdout))), m_active(active) {
    if (m_active && !m_is_tty) {
      std::fwrite(m_header.data(), 1, m_header.size(), stdout);
      std::fflush(stdout);
    }
  }

  ~LiveMessage() { done(); }

  // Non-copyable: owns the "print header once" invariant
  LiveMessage(const LiveMessage &) = delete;
  LiveMessage &operator=(const LiveMessage &) = delete;

  //! Update the status message. On TTY overwrites the current line; on
  //! non-TTY buffers silently until done() is called.
  void update(std::string_view msg) {
    if (!m_active || m_done)
      return;
    m_last_msg = msg;
    if (m_is_tty) {
      std::fputs("\r", stdout);
      std::fwrite(m_header.data(), 1, m_header.size(), stdout);
      std::fwrite(msg.data(), 1, msg.size(), stdout);
      std::fflush(stdout);
    }
  }

  //! Update the status message. see @ref update()
  void operator()(std::string_view msg) { return update(msg); }

  //! Finalise: print last message + optional @p post, then newline.
  //! Safe to call multiple times (only the first call has effect).
  void done(std::string_view post = {}) {
    if (!m_active || m_done)
      return;
    m_done = true;
    if (m_is_tty) {
      std::fputs("\r", stdout);
      std::fwrite(m_header.data(), 1, m_header.size(), stdout);
    }
    std::fwrite(m_last_msg.data(), 1, m_last_msg.size(), stdout);
    if (!post.empty())
      std::fwrite(post.data(), 1, post.size(), stdout);
    std::fputc('\n', stdout);
    std::fflush(stdout);
  }

private:
  std::string m_header;
  std::string m_last_msg;
  bool m_is_tty;
  bool m_active;
  bool m_done{false};
};

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
  When stdout is not a TTY (i.e., piped), prints comma-separated percentages
  instead of a bar: "0%, 10%, 20%, ... 100%\n". Prints approximately
  @p Length / 5 values, spaced evenly.

  The @p Length template parameter controls output width.
  For TTY mode: total bar width in characters.
  For non-TTY mode: determines how many percentage values are printed (~Length/5).

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

  @tparam Length Total character width of the output (default 60).

  @note Thread-safe; safe to call from multiple threads simultaneously.
  But may add overhead.

  @note Counts converted to `int` internally, so the maximum supported 
  iteration count is INT_MAX.
*/
template <std::size_t Length = 60>
class ProgressBar {
private:
  int m_max;
  bool m_print;
  // Atomicly tracks progress
  std::atomic<int> m_progress{0};
  // TTY: bar fill chars printed last time. Non-TTY: last percentage slot printed.
  std::atomic<int> m_last_fill{-1};
  // Set to true once the final bar (with '\n') has been printed; guarded by
  // omp critical so no further '\r' output can overwrite the newline.
  bool m_done{false};

  static constexpr std::size_t bar_width = (Length >= 7) ? (Length - 7) : 1;

public:
  /*!
    @brief Construct progress bar for @p max iterations.
    @param max Total number of iterations (denominator for percentage).
              Accepts any integral type; converted to int internally.
    @param print Runtime switch; if false, class does nothing.
  */
  template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
  ProgressBar(T max, bool print = true)
    : m_max(static_cast<int>(max)), m_print(print) {
    if (m_print)
      print_prog_bar(0);
  }

  //! Atomically increment progress counter and print updated bar.
  void update() {
    if (!m_print)
      return;
    const int progress = ++m_progress;
    if (progress > m_max)
      return;
    const bool is_final = (progress == m_max);

    // skip print if bar fill is unchanged (avoids critical section overhead)
    const int current = int(bar_width * double(progress) / double(m_max));
    if (!is_final && current <= m_last_fill.load(std::memory_order_relaxed))
      return;

    print_prog_bar(progress);
  }

private:
  //----------------------------------------------------------------------------

  // Prints progress bar (when stdout is a normal tty)
  void print_prog_bar(int progress) {
    const bool is_tty = isatty(fileno(stdout));
    if (!is_tty) {
      print_prog_bar_notty(progress);
      return;
    }

    const bool is_final = (progress == m_max);

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
    p += snprintf(p, sizeof(buf) - std::size_t(p - buf) - 1, "%d%%",
                  int(100.0 * double(progress) / double(m_max)));
    *p++ = is_final ? '\n' : '\r';
    *p = '\0';

    // serialise writes; skip if final bar already printed (prevents a lagging
    // thread's '\r' from overwriting the final '\n' on the terminal)
    bool did_print = false;
#pragma omp critical
    {
      if (!m_done) {
        fputs(buf, stdout);
        fflush(stdout);
        did_print = true;
        if (is_final)
          m_done = true;
      }
    }
    if (did_print)
      m_last_fill.store(current, std::memory_order_relaxed);
  }

  //----------------------------------------------------------------------------

  // Prints progress "bar" (comma separated %) (when stdout is NOT a normal tty)
  void print_prog_bar_notty(int progress) {
    const bool is_initial = (progress == 0);
    const bool is_final = (progress == m_max);

    const int pct =
      is_final ? 100 : int(100.0 * double(progress) / double(m_max));

    // number of entries ~ Length/5 chars each; interval between prints in pct
    static constexpr int n_entries = static_cast<int>(Length) / 5;
    static constexpr int interval =
      (n_entries > 1) ? (100 / (n_entries - 1)) : 100;
    const int slot = pct / interval;

    if (!is_initial && !is_final &&
        slot <= m_last_fill.load(std::memory_order_relaxed))
      return;

    char buf[8];
    if (is_final)
      snprintf(buf, sizeof(buf), "100%%\n");
    else
      snprintf(buf, sizeof(buf), "%d%%, ", pct);

    bool did_print = false;
#pragma omp critical
    {
      if (!m_done) {
        fputs(buf, stdout);
        fflush(stdout);
        did_print = true;
        if (is_final)
          m_done = true;
      }
    }
    if (!is_final && did_print)
      m_last_fill.store(slot, std::memory_order_relaxed);
  }
};

} // namespace qip
