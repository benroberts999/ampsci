#pragma once
#include <chrono>
// #include <filesystem> // needs gcc >8
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

namespace SafeProfiler {

constexpr bool do_profile = true;
// auto sp1 = SafeProfiler::profile(__func__);
// #include "IO/SafeProfiler.hpp"

struct StrDoubleUnsigned {
  StrDoubleUnsigned(std::string is, double id, unsigned iu)
      : s(is), d(id), u(iu) {}
  std::string s;
  double d;
  unsigned u;
};

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
  std::vector<std::pair<std::string, std::chrono::microseconds>> m_log;
  // std::vector<std::tuple<std::string, double, unsigned>> m_log_sorted;
  std::vector<StrDoubleUnsigned> m_log_sorted;

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
    for (const auto &[name, time, count] : m_log_sorted) {
      std::cout << name << " x " << count << " = " << time / 1000.0 << " ms\n";
      of << name << " x " << count << " = " << time / 1000.0 << " ms\n";
    }
  }
};

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

struct BlankClass {};

inline auto profile(const char *in_name, const char *extra = "") {
  if constexpr (do_profile)
    return Profiler(in_name, extra);
  else
    return BlankClass();
}

} // namespace SafeProfiler
