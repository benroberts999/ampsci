#pragma once
#include <chrono>
#include <iomanip> // setprecision
#include <iostream>
#include <sstream> // stringstream
#include <string>

//! In-out (timers, profilers, and read/write data)
namespace IO {

/*!
  @brief Class that uses std::chrono to easily time code
  @details
  Usage:
  Will automatically time a routine, based on scope.
  Automatically starts timing on construction, and returns total time elapsed when
  timer goes out of scope. Also "fancy" commands to start/stop the timer if need
  be (e.g., to not time a certain part of the code)

    - sw.start() -- starts timing
    - sw.stop()  -- 'pauses' timing
    - sw.reset() -- clears timer. Needs to be re-stared. Forgets all times.

    - start/stop lets you time indevidual sections of code, while not timing
    others.
    - Calling start() again will start a new "lap", and save current to total.

    - reading_ms() -- returns total elapsed time as double, in ms
    - lap_reading_ms() -- same, but only returns time since last start()

    - reading_str() and lap_reading_str() -- As above, but outputs as formatted
  string, in units of ms,s,mins, or hours (e.g., "1.56 s" or "2.10 hours")

*/
class ChronoTimer {
public:
  //! Constructs and immediately starts the timer. If @p print is true (default),
  //! prints elapsed time to stdout on destruction.
  ChronoTimer(std::string_view in_name = "", bool print = true);

  /*!
    @brief Starts (or resumes) timing.
    @details
    If the timer is already running, the current lap is saved to the total
    before restarting -- equivalent to calling stop() then start().
  */
  void start();

  //! Pauses timing and accumulates elapsed time into the running total.
  void stop();

  //! Clears all accumulated time and stops the timer. Must call start() to resume.
  void reset();

  //! Equivalent to reset() followed by start().
  void restart();

  //! Prints elapsed time if a name was given at construction.
  ~ChronoTimer();

  //! Total elapsed time in milliseconds (accumulated laps + current lap if running).
  double reading_ms() const;

  //! Elapsed time in milliseconds for the current lap only (since last start()).
  double lap_reading_ms() const;

  //! Total elapsed time as a human-readable string (e.g., "1.56 s", "2.10 hours").
  std::string reading_str() const;

  //! Current lap elapsed time as a human-readable string.
  std::string lap_reading_str() const;

private:
  std::string name;
  bool m_print;
  bool running;
  std::chrono::high_resolution_clock::time_point tstart{};
  double total_time_ms;

  std::string convertHR(double t) const;
};

//==============================================================================
inline ChronoTimer::ChronoTimer(std::string_view in_name, bool print)
  : name(in_name), m_print(print), running(false), total_time_ms(0) {
  start();
}
//==============================================================================
inline ChronoTimer::~ChronoTimer() {
  if (!m_print)
    return;
  if (name.empty()) {
    std::cout << reading_str() << "\n";
  } else {
    std::cout << name << ": T = " << reading_str() << "\n";
  }
}

//==============================================================================
inline void ChronoTimer::start() {
  // note: will over-ride any existing reading!
  if (running)
    stop();
  running = true;
  tstart = std::chrono::high_resolution_clock::now();
}

//==============================================================================
inline void ChronoTimer::stop() {
  if (!running)
    return;

  double current_time = lap_reading_ms();
  running = false; //"turn off" stopwatch

  // update total time
  total_time_ms += current_time;
}

//==============================================================================
inline void ChronoTimer::restart() {
  reset();
  start();
}

//==============================================================================
inline void ChronoTimer::reset() {
  running = false;
  total_time_ms = 0;
}

//==============================================================================
inline double ChronoTimer::lap_reading_ms() const
// Returns value for current riming run (lap)
// Returns double (milliseconds)
{

  if (!running)
    return 0;

  std::chrono::high_resolution_clock::time_point tcurrent =
    std::chrono::high_resolution_clock::now();

  auto duration =
    std::chrono::duration_cast<std::chrono::microseconds>(tcurrent - tstart)
      .count();

  return ((double)duration) * 1.0e-3;
}

//==============================================================================
inline double ChronoTimer::reading_ms() const
// Returns total value for timnig run
// Returns double (milliseconds)
{
  return lap_reading_ms() + total_time_ms;
}
//==============================================================================

inline std::string ChronoTimer::reading_str() const {
  return convertHR(reading_ms());
}
//==============================================================================
inline std::string ChronoTimer::lap_reading_str() const {
  return convertHR(lap_reading_ms());
}

//==============================================================================
inline std::string ChronoTimer::convertHR(double t) const
// Convers double (in ms) into formmated 2 d.p. string in units of either
// ms, s, mins, or hours, depending on size.
{
  double ot;

  std::string un;
  if (t < 1000) {
    ot = t;
    un = "ms";
  } else if (t < 60000) {
    ot = t / 1000.;
    un = "s";
  } else if (t < 3600000) {
    ot = t / 60000.;
    un = "mins";
  } else {
    ot = t / 3600000.;
    un = "hours";
  }

  std::stringstream ss;
  ss << std::fixed << std::setprecision(2) << ot;
  return ss.str() + " " + un;
}

} // namespace IO
