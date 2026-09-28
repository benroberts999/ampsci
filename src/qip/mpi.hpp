#pragma once
/*! @file
  @brief MPI support: a plain-data buffer, and a manager/worker task layer.
  Include instead of `<mpi.h>`; compiles with or without MPI.

  @details
  With AMPSCI_USE_MPI defined (the MPI switch in the Makefile; the compiler
  must then be an MPI wrapper such as mpicxx) this wraps MPI. Without it
  every function is the single-process identity: rank 0 of 1, broadcasts and
  reductions leave their data unchanged, and no worker ever exists.

  \par Manager/worker model
  Only rank 0 runs the program. Every other rank waits in worker_loop() until
  rank 0 starts a task (start_task()); the task function then runs on every
  rank, including rank 0, and returns the workers to the loop. So the rest of
  the program is unaware of MPI: only the few functions written as tasks,
  main(), and this file use it. File I/O and printing happen on rank 0 only.

  \par Tasks
  A task is a plain function, registered once by name with a static
  TaskRegistration object. Rank 0 calls start_task(name) and then the same
  function (so that it runs everywhere); inside, rank 0 packs the inputs into
  a Buffer and every rank calls broadcast() to share them. Each rank does its
  share of the work (mine()), and the results are combined with the
  reductions. Tasks must not start other tasks.

  \par Buffer
  Objects that a task needs are sent as their construction parameters (a
  plain-data struct the class provides, e.g. HF::HartreeFock::Params from
  params(), with a constructor taking it back); the task puts the fields and
  gets them in the same order on the workers.

  Run as: mpirun -np N ./ampsci input.in (OpenMP threads within each rank).
*/
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>
#if defined(AMPSCI_USE_MPI)
// Only the C interface is used; the (deprecated) C++ bindings do not compile
// warning-free
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#include <mpi.h>
#pragma GCC diagnostic pop
#include <chrono>
#include <cstdlib>
#include <map>
#include <thread>
#endif

namespace qip::mpi {

#if defined(AMPSCI_USE_MPI)
//! True if compiled with MPI support, false otherwise.
constexpr bool use_mpi = true;
#else
//! True if compiled with MPI support, false otherwise.
constexpr bool use_mpi = false;
#endif

//==============================================================================
/*!
  @brief Plain-data byte stream, used to send objects between ranks.
  @details
  put() appends a value, get() reads the next one; values must be read in the
  order they were written. Supports trivially copyable types (numbers, enums,
  bools, plain structs), std::string, std::vector and std::optional of such.
*/
class Buffer {
public:
  //! Appends a trivially copyable value
  template <typename T>
  void put(const T &value) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "Buffer::put: type must be trivially copyable (or a "
                  "string/vector of such)");
    const auto *bytes = reinterpret_cast<const char *>(&value);
    m_data.insert(m_data.end(), bytes, bytes + sizeof(T));
  }

  //! Appends a string (length, then characters)
  void put(const std::string &value) {
    put(value.size());
    m_data.insert(m_data.end(), value.begin(), value.end());
  }

  //! Appends a vector (length, then elements)
  template <typename T>
  void put(const std::vector<T> &values) {
    put(values.size());
    if constexpr (std::is_trivially_copyable_v<T>) {
      const auto *bytes = reinterpret_cast<const char *>(values.data());
      m_data.insert(m_data.end(), bytes, bytes + sizeof(T) * values.size());
    } else {
      for (const auto &value : values) {
        put(value);
      }
    }
  }

  //! Appends an optional (whether it holds a value, then the value)
  template <typename T>
  void put(const std::optional<T> &value) {
    put(value.has_value());
    if (value) {
      put(*value);
    }
  }

  //! Reads the next value into value
  template <typename T>
  void get(T &value) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "Buffer::get: type must be trivially copyable (or a "
                  "string/vector of such)");
    assert(m_pos + sizeof(T) <= m_data.size() && "Buffer: read past end");
    std::memcpy(&value, m_data.data() + m_pos, sizeof(T));
    m_pos += sizeof(T);
  }

  //! Reads the next string into value
  void get(std::string &value) {
    const auto n = get<std::size_t>();
    assert(m_pos + n <= m_data.size() && "Buffer: read past end");
    value.assign(m_data.data() + m_pos, n);
    m_pos += n;
  }

  //! Reads the next vector into values
  template <typename T>
  void get(std::vector<T> &values) {
    const auto n = get<std::size_t>();
    values.resize(n);
    if constexpr (std::is_trivially_copyable_v<T>) {
      assert(m_pos + sizeof(T) * n <= m_data.size() && "Buffer: read past end");
      std::memcpy(values.data(), m_data.data() + m_pos, sizeof(T) * n);
      m_pos += sizeof(T) * n;
    } else {
      for (auto &value : values) {
        get(value);
      }
    }
  }

  //! Reads the next optional into value
  template <typename T>
  void get(std::optional<T> &value) {
    if (get<bool>()) {
      T held{};
      get(held);
      value = std::move(held);
    } else {
      value.reset();
    }
  }

  //! Reads and returns the next value
  template <typename T>
  T get() {
    T value{};
    get(value);
    return value;
  }

  //! The raw bytes (for sending)
  std::vector<char> &bytes() { return m_data; }
  //! Moves the read position back to the start
  void rewind() { m_pos = 0; }

private:
  std::vector<char> m_data{};
  std::size_t m_pos{0};
};

//==============================================================================
//! Rank of this process (0 without MPI)
inline int rank() {
#if defined(AMPSCI_USE_MPI)
  int r = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
#else
  return 0;
#endif
}

//! Number of ranks (1 without MPI)
inline int size() {
#if defined(AMPSCI_USE_MPI)
  int n = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  return n;
#else
  return 1;
#endif
}

//! True on rank 0: the rank that runs the program
inline bool root() { return rank() == 0; }

//! True if item i of a list shared round robin between the ranks belongs to
//! this rank (i % size == rank)
inline bool mine(std::size_t i) {
  return i % static_cast<std::size_t>(size()) ==
         static_cast<std::size_t>(rank());
}

//! Short description for the screen, e.g. "4 ranks." Empty without MPI.
inline std::string details() {
  if (!use_mpi)
    return "";
  const auto n = size();
  return std::to_string(n) + (n == 1 ? " rank." : " ranks.");
}

//==============================================================================
//! Sends rank 0's buffer to every rank (the others' contents are replaced),
//! and rewinds it for reading. Collective: every rank must call it.
inline void broadcast([[maybe_unused]] Buffer *buffer) {
  assert(buffer != nullptr);
#if defined(AMPSCI_USE_MPI)
  auto &bytes = buffer->bytes();
  unsigned long long n = bytes.size();
  MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  bytes.resize(static_cast<std::size_t>(n));
  // MPI counts are int: send in chunks
  constexpr std::size_t chunk = std::size_t(1) << 30;
  for (std::size_t i0 = 0; i0 < bytes.size(); i0 += chunk) {
    const auto count = static_cast<int>(std::min(chunk, bytes.size() - i0));
    MPI_Bcast(bytes.data() + i0, count, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif
  buffer->rewind();
}

#if defined(AMPSCI_USE_MPI)
namespace detail {
//! In-place reduction of n doubles over all ranks; MPI counts are int, so in
//! chunks
inline void all_reduce(double *data, std::size_t n, MPI_Op op) {
  constexpr std::size_t chunk = std::size_t(1) << 30;
  for (std::size_t i0 = 0; i0 < n; i0 += chunk) {
    const auto count = static_cast<int>(std::min(chunk, n - i0));
    MPI_Allreduce(MPI_IN_PLACE, data + i0, count, MPI_DOUBLE, op,
                  MPI_COMM_WORLD);
  }
}
} // namespace detail
#endif

//! Sums the n values at data over all ranks, in place: every rank gets the
//! total. Collective. No-op without MPI.
inline void all_reduce_sum([[maybe_unused]] double *data,
                           [[maybe_unused]] std::size_t n) {
#if defined(AMPSCI_USE_MPI)
  detail::all_reduce(data, n, MPI_SUM);
#endif
}

//! Element-wise maximum of the n values at data over all ranks, in place.
//! Collective. No-op without MPI.
inline void all_reduce_max([[maybe_unused]] double *data,
                           [[maybe_unused]] std::size_t n) {
#if defined(AMPSCI_USE_MPI)
  detail::all_reduce(data, n, MPI_MAX);
#endif
}

//! Sums a count over all ranks, in place. Collective. No-op without MPI.
inline void all_reduce_sum([[maybe_unused]] std::size_t *count) {
#if defined(AMPSCI_USE_MPI)
  unsigned long long total = *count;
  MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
  *count = static_cast<std::size_t>(total);
#endif
}

//==============================================================================
//! A task: runs on every rank (see file description)
using TaskFunction = void (*)();

#if defined(AMPSCI_USE_MPI)
namespace detail {
//! Registered tasks, by name
inline std::map<std::string, TaskFunction> &task_registry() {
  static std::map<std::string, TaskFunction> registry;
  return registry;
}

//! Commands from rank 0 to the workers: a task name, or the stop command
constexpr std::size_t command_length = 128;
inline const std::string stop_command = "!stop";

//! Sends a command from rank 0 to every worker. Non-blocking collective, so
//! that idle workers can sleep between polls (a blocking wait spins, and
//! would steal cores from rank 0's threads)
inline std::string exchange_command(const std::string &command) {
  char text[command_length] = {};
  if (root()) {
    assert(command.size() < command_length && "task name too long");
    std::copy(command.begin(), command.end(), text);
  }
  MPI_Request request;
  MPI_Ibcast(text, int(command_length), MPI_CHAR, 0, MPI_COMM_WORLD, &request);
  if (root()) {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
  } else {
    int done = 0;
    MPI_Test(&request, &done, MPI_STATUS_IGNORE);
    while (!done) {
      std::this_thread::sleep_for(std::chrono::milliseconds(2));
      MPI_Test(&request, &done, MPI_STATUS_IGNORE);
    }
  }
  return std::string(text);
}

//! Set once the workers have been stopped and MPI finalised
inline bool &finalised() {
  static bool value = false;
  return value;
}

//! On rank 0: stops the workers, then finalises MPI (once)
inline void shutdown() {
  if (finalised())
    return;
  finalised() = true;
  if (root()) {
    exchange_command(stop_command);
  }
  MPI_Finalize();
}
} // namespace detail
#endif

//! Registers a task by name, at static initialisation: declare one of these
//! (at namespace scope) next to the task function
class TaskRegistration {
public:
  TaskRegistration([[maybe_unused]] const std::string &name,
                   [[maybe_unused]] TaskFunction function) {
#if defined(AMPSCI_USE_MPI)
    assert(detail::task_registry().count(name) == 0 && "duplicate task name");
    detail::task_registry()[name] = function;
#endif
  }
};

//! On rank 0: tells the workers to run the named task. Rank 0 must then run
//! the same task function itself. No-op with a single rank.
inline void start_task([[maybe_unused]] const std::string &name) {
#if defined(AMPSCI_USE_MPI)
  assert(root() && "only rank 0 starts tasks");
  assert(detail::task_registry().count(name) == 1 && "unregistered task");
  if (size() > 1) {
    detail::exchange_command(name);
  }
#endif
}

//! On the workers (rank > 0): runs the tasks rank 0 starts, until it stops
//! them. No-op without MPI.
inline void worker_loop() {
#if defined(AMPSCI_USE_MPI)
  assert(!root());
  for (;;) {
    const auto command = detail::exchange_command("");
    if (command == detail::stop_command)
      break;
    const auto task = detail::task_registry().find(command);
    assert(task != detail::task_registry().end() && "unregistered task");
    task->second();
  }
#endif
}

//==============================================================================
/*!
  @brief Initialises MPI on construction; on destruction (or at exit, if the
  program ends through std::exit) rank 0 stops the workers, and MPI is
  finalised.
  @details Create one at the start of main(), before any other MPI call. Then
  the workers call worker_loop(), and only rank 0 runs the program.
*/
class Environment {
public:
  Environment([[maybe_unused]] int *argc, [[maybe_unused]] char ***argv) {
#if defined(AMPSCI_USE_MPI)
    // MPI is called from the main thread only (not inside OpenMP regions)
    int provided = 0;
    MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
    std::atexit(detail::shutdown);
#endif
  }
  ~Environment() {
#if defined(AMPSCI_USE_MPI)
    detail::shutdown();
#endif
  }
  Environment(const Environment &) = delete;
  Environment &operator=(const Environment &) = delete;
};

} // namespace qip::mpi
