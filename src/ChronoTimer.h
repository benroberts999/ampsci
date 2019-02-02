#pragma once
#include <chrono>
#include <string>
/*
Class that uses std::chrono to easily time code
Usage:
  Create object. By default will not start the timer, unless 'true'
  is given as optional input, e.g.:
    ChronoTimer sw; //will not start the timing.
    ChronoTimer sw(true); will

  sw.start() -- starts timing
  sw.stop()  -- 'pauses' timing
  sw.reset() -- clears timer. Needs to be re-stared. Forgets all times.

  start/stop lets you time indevidual sections of code, while not timing
  others.
  Calling start() again will start a new "lap", and save current to total.

  reading_ms() -- returns total elapsed time as double, in ms
  lap_reading_ms() -- same, but only returns time since last start()

  reading_str() and lap_reading_str() -- As above, but outputs as
  formatted string, in units of ms,s,mins, or hours
  (e.g., "1.56 s" or "2.10 hours")

*/

class ChronoTimer {
public:
  ChronoTimer(bool auto_start = false);
  void start();
  void stop();
  void reset();

  double reading_ms();
  double lap_reading_ms();
  std::string reading_str();
  std::string lap_reading_str();

private:
  bool running;
  std::chrono::high_resolution_clock::time_point tstart;
  std::string convertHR(double t);
  double total_time_ms;
};
