#pragma once
#include <chrono>
// #include <filesystem> // needs gcc >8
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

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
#else
constexpr bool do_profile = false;
#endif

struct StrDoubleUnsigned {
  StrDoubleUnsigned(std::string is, double id, unsigned iu)
      : s(std::move(is)), d(id), u(iu) {}
  std::string s;
  double d;
  unsigned u;
};

//******************************************************************************
// Thread-safe profiler (timing) tool
class ProfileLog {
  friend class Profiler;
  ProfileLog(){};
  ~ProfileLog() {
    if (m_log.empty())
      return;
#pragma omp critical(finalise)
    {
      for (const auto &le : m_log) {
        bool new_entry = true;
        auto each_time = static_cast<double>(le.second.count());
        for (auto &[name, time, count] : m_log_sorted) {
          if (name == le.first) {
            new_entry = false;
            time += each_time;
            count++;
            break;
          }
        }
        if (new_entry)
          m_log_sorted.emplace_back(le.first, each_time, 1);
      }
      write_log();
    }
  }
  std::vector<std::pair<std::string, std::chrono::microseconds>> m_log = {};
  std::vector<StrDoubleUnsigned> m_log_sorted = {};

  void write_log() {
    const std::string title = "ProfileLog";
    auto fname = title + "_0";
    std::string ext = ".txt";
    auto file_existsQ = [](const std::string &fileName) {
      std::ifstream infile(fileName);
      return infile.good();
    };
    // std::filesystem::exists(fname + ext) // needs GCC v>8
    for (int i = 1; file_existsQ(fname + ext); ++i) {
      fname = title + "_" + std::to_string(i);
    }
    std::ofstream of(fname + ext);
    std::cout << "\nProfile:\n";
    const auto *p_worst = &(m_log_sorted.front());
    for (const auto &item : m_log_sorted) {
      const auto &[name, time, count] = item;
      std::cout << name << " x " << count << " = " << time / 1000.0 << " ms\n";
      of << name << " x " << count << " = " << time / 1000.0 << " ms\n";
      if (time >= p_worst->d) {
        p_worst = &item;
      }
    }
    const auto &[name, time, count] = *p_worst;
    std::cout << "\n ** Worst: \n";
    std::cout << name << " x " << count << " = " << time / 1000.0 << " ms\n";
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
      log.m_log.emplace_back(name, duration);
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
