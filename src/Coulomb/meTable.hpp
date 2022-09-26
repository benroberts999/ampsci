#pragma once
#include "QkTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cstdint>
#include <unordered_map>

namespace Coulomb {

//! index type for set of 2 orbitals {nk,nk} -> nk4Index
using nk2Index = uint32_t;

//! Look-up table for matrix elements. Note: does not assume any symmetry: (a,b)
//! is stored independantly of (b,a). In general, maps a pair of DiracSpinors to
//! a single value (of any type, T).
template <typename T = double> class meTable {

private:
  std::unordered_map<nk2Index, T> m_data{};

public:
  auto operator-> () { return &m_data; }

  //! Adds new element to table. If already exists, does nothing (does not
  //! update)
  void add(const DiracSpinor &a, const DiracSpinor &b, T value) {
    m_data.insert({FormIndex(a.nk_index(), b.nk_index()), std::move(value)});
  }

  //! Updates given element in table. If element not yet present, adds it
  void update(const DiracSpinor &a, const DiracSpinor &b, T value) {
    m_data.insert_or_assign(FormIndex(a.nk_index(), b.nk_index()),
                            std::move(value));
  }

  //! Checks if given element is in the table
  [[nodiscard]] bool contains(const DiracSpinor &a,
                              const DiracSpinor &b) const {
    return m_data.find(FormIndex(a.nk_index(), b.nk_index())) != m_data.cend();
  }

      //! Gets pointer to const requested element. If element not present,
      //! returns nullptr
      [[nodiscard]] const T *get(const DiracSpinor &a,
                                 const DiracSpinor &b) const {
    const auto map_it = m_data.find(FormIndex(a.nk_index(), b.nk_index()));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }
  //! Gets pointer to mutable requested element. If element not present, returns
  //! nullptr
  [[nodiscard]] T *get(const DiracSpinor &a, const DiracSpinor &b) {
    const auto map_it = m_data.find(FormIndex(a.nk_index(), b.nk_index()));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }

  private :
      // Converts given set of nkIndex's (in any order) to nk4Index
      [[nodiscard]] nk2Index FormIndex(nkIndex a, nkIndex b) const {
        static_assert(sizeof(nk2Index) == 2 * sizeof(nkIndex));
        static_assert(sizeof(nkIndex) * 8 == 16);
        return (nk2Index)b + ((nk2Index)a << 16);
      }
};

} // namespace Coulomb