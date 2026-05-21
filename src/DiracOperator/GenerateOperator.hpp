#pragma once
#include "DiracOperator/TensorOperator.hpp"
#include "IO/InputBlock.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace DiracOperator {

/*!
  @brief Operators are self-registering: each operator class exposes a static
  generate() factory and registers it via a file-scope Register<T> object.

  ## Adding a new operator

  **Internal operators** (need concrete type access from outside): add a
  `static generate()` to the class in its `.hpp`, and add a `Register<T>`
  entry in `RegisterOperators.cpp`.

  **External operators** (never constructed directly from outside): put the
  class definition, `static generate()`, and `Register<T>` all in a single
  `.cpp` file. Drop it in the build and it self-registers at startup.

  The minimal static factory signature:

  \code{.cpp}
  static std::unique_ptr<TensorOperator>
  generate(const IO::InputBlock &input, const Wavefunction &wf);
  \endcode

  From the command line:
  ```shell
    ampsci -o                 # list all operators
    ampsci -o OperatorName    # show options for OperatorName
  ```
*/

//! Function-pointer signature shared by every operator factory.
using FactoryFn = std::unique_ptr<TensorOperator> (*)(const IO::InputBlock &,
                                                      const Wavefunction &);

/*!
  @brief One entry in the operator registry.
*/
struct OperatorEntry {
  std::string name;
  std::string description;
  FactoryFn factory;
};

/*!
  @brief Singleton registry of all compiled-in operators.
  @details
  Populated at program startup via static initialisers (one per operator) and
  never modified once main() begins. Uses construct-on-first-use to avoid the
  static-initialisation order fiasco.
*/
class Registry {
public:
  //! Access the singleton instance.
  static Registry &get() {
    static Registry instance;
    return instance;
  }

  //! Append a new entry. Normally called only by Register<T>.
  void add(std::string name, std::string description, FactoryFn fn) {
    m_entries.push_back({std::move(name), std::move(description), fn});
  }

  //! All registered operators, in registration order.
  const std::vector<OperatorEntry> &entries() const { return m_entries; }

private:
  Registry() : m_entries{} {}
  std::vector<OperatorEntry> m_entries;
};

/*!
  @brief Constructing a Register<T> adds T::generate to the Registry.
  @details
  Place one of these at file scope (inside an anonymous namespace) in a .cpp
  file to register an operator. T must have a static generate() with the
  FactoryFn signature; the template instantiation enforces this at compile time.

  \code{.cpp}
  namespace {
  const DiracOperator::Register<MyOp> r{"MyOp", "One-line description"};
  }
  \endcode
*/
template <typename T>
struct Register {
  Register(const char *name, const char *description) {
    Registry::get().add(name, description, &T::generate);
  }
};

//! Returns a unique_ptr (polymorphic) to the requested operator.
/*!
  @details
  Looks up @p operator_name in the Registry (case-insensitive) and calls its
  factory. Returns a NullOperator and prints an error if not found.
*/
std::unique_ptr<DiracOperator::TensorOperator>
generate(std::string_view operator_name, const IO::InputBlock &input,
         const Wavefunction &wf);

//! Print the list of compiled-in operators (name + description).
void list_operators();

} // namespace DiracOperator
