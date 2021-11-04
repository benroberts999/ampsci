#pragma once
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>

// XXX Note: Very slow - should not be used for functions called large number of
// times. Partly issue is due to omp critical - can be fixed

// #define PROFILE(x)
// #ifdef IOPROFILER
// #define PROFILE(x) [[maybe_unused]] auto sp = IO::Profile::safeProfiler(x);
// #endif

//! Simple thread-safe scope-based profiler function.
//! @details Must compile with the pre-processor macro IOPROFILER defined for
//! the pofiler to be used. If not, it does nothing.
//!
//! Use: add following line to scope you wish to profile:
//!  - [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
//!
//! #include "IO/SafeProfiler.hpp"
namespace IO::Profile {

#ifdef IOPROFILER
constexpr bool do_profile = true;
// #define PROFILE(x) [[maybe_unused]] auto sp = IO::Profile::safeProfiler(x);
#else
constexpr bool do_profile = false;
#endif

//******************************************************************************
// Thread-safe profiler (timing) tool
class ProfileLog {

  std::unordered_map<std::string,
                     std::pair<unsigned long, std::chrono::microseconds>>
      m_log = {};

  friend class Profiler;

  ProfileLog(){};

  ~ProfileLog() {
    if (m_log.empty())
      return;
#pragma omp critical(finalise)
    write_log();
  }

  void add(const std::string &name, std::chrono::microseconds duration) {
    // if key already exists, this does nothing and insertedQ==false:
    auto [it, insertedQ] = m_log.insert({name, {1, duration}});
    if (!insertedQ) {
      // an entry already exists, so update it:
      ++(it->second.first);
      it->second.second += duration;
    }
  }

  void write_log() {
    const std::string title = "ProfileLog";
    const std::string ext = ".txt";
    auto file_existsQ = [](const std::string &fileName) {
      std::ifstream infile(fileName);
      return infile.good();
    };
    // std::filesystem::exists(fname + ext) // needs GCC v>8
    // Ensure filename is unique:
    auto fname = title + "_0";
    for (int i = 1; file_existsQ(fname + ext); ++i) {
      fname = title + "_" + std::to_string(i);
    }
    std::ofstream of(fname + ext);
    std::cout << "\nProfile:\n";
    of << "# Function call_count total_time(ms)\n";
    for (const auto &item : m_log) {
      const auto &[name, entry] = item;
      const auto &[count, time] = entry;
      const auto dtime = static_cast<double>(time.count());
      std::cout << name << " x " << count << " = " << dtime / 1000.0 << " ms\n";
      of << name << ", " << count << ", " << dtime / 1000.0 << "\n";
    }
  }
};

//******************************************************************************
class Profiler {
  friend class ProfileLog;
  const std::string name;
  const std::chrono::high_resolution_clock::time_point tstart;
  inline static ProfileLog log;

public:
  Profiler(const char *in_name, const char *extra = "")
      : name([&]() {
          return std::string(extra) == ""
                     ? std::string(in_name)
                     : std::string(in_name) + "_" + std::string(extra);
        }()),
        tstart(std::chrono::high_resolution_clock::now()) {}

  ~Profiler() {
#pragma omp critical(add_entry)
    {
      const auto now = std::chrono::high_resolution_clock::now();
      const auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(now - tstart);
      log.add(name, duration);
    }
  }
};

//******************************************************************************
//! @details This function does the timing/profiling. It returns an object (a
//! profile logger), which must survive to the end of the scope you are trying
//! to profile. Call function once at beginning of scope, store its returned
//! log.
inline auto safeProfiler(const char *in_name, const char *extra = "") {
  if constexpr (do_profile)
    return Profiler(in_name, extra);
  else
    return 0;
}

} // namespace IO::Profile
