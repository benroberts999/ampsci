#pragma once
#include "Physics/AtomData.hpp"
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
  auto operator->() { return &m_data; }

  //! Adds new element to table. If already exists, does nothing (does not
  //! update)
  void add(const DiracSpinor &a, const DiracSpinor &b, T value) {
    m_data.insert({FormIndex(a.nk_index(), b.nk_index()), std::move(value)});
  }
  //! Adds new element to table. If already exists, does nothing (does not
  //! update)
  void add(DiracSpinor::Index a, DiracSpinor::Index b, T value) {
    m_data.insert({FormIndex(a, b), std::move(value)});
  }

  //! Adds elements from one Table into another (by copy)
  void add(const meTable<T> &other) {
    m_data.insert(other.cbegin(), other.cend());
  }

  //! Updates given element in table. If element not yet present, adds it
  void update(const DiracSpinor &a, const DiracSpinor &b, T value) {
    m_data.insert_or_assign(FormIndex(a.nk_index(), b.nk_index()),
                            std::move(value));
  }
  //! Updates given element in table. If element not yet present, adds it
  void update(DiracSpinor::Index a, DiracSpinor::Index b, T value) {
    m_data.insert_or_assign(FormIndex(a, b), std::move(value));
  }

  //! Checks if given element is in the table
  [[nodiscard]] bool contains(const DiracSpinor &a,
                              const DiracSpinor &b) const {
    return m_data.find(FormIndex(a.nk_index(), b.nk_index())) != m_data.cend();
  }
  //! Checks if given element is in the table
  [[nodiscard]] bool contains(DiracSpinor::Index a,
                              DiracSpinor::Index b) const {
    return m_data.find(FormIndex(a, b)) != m_data.cend();
  }

  //! Gets pointer to const requested element. If element not present,
  //! returns nullptr
  [[nodiscard]] const T *get(const DiracSpinor &a, const DiracSpinor &b) const {
    const auto map_it = m_data.find(FormIndex(a.nk_index(), b.nk_index()));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }
  //! Gets pointer to const requested element. If element not present,
  //! returns nullptr
  [[nodiscard]] const T *get(DiracSpinor::Index a, DiracSpinor::Index b) const {
    const auto map_it = m_data.find(FormIndex(a, b));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }

  //! Gets pointer to mutable requested element. If element not present, returns
  //! nullptr
  [[nodiscard]] T *get(const DiracSpinor &a, const DiracSpinor &b) {
    const auto map_it = m_data.find(FormIndex(a.nk_index(), b.nk_index()));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }
  //! Gets pointer to mutable requested element. If element not present, returns
  //! nullptr
  [[nodiscard]] T *get(DiracSpinor::Index a, DiracSpinor::Index b) {
    const auto map_it = m_data.find(FormIndex(a, b));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }

  //! Gets value of requested element. If element not present,
  //! returns zero (or default constructed T)
  [[nodiscard]] T getv(const DiracSpinor &a, const DiracSpinor &b) const {
    const auto ptr = get(a, b);
    return ptr ? *ptr : T{};
  }
  //! Gets value of requested element. If element not present,
  //! returns zero (or default constructed T)
  [[nodiscard]] T getv(DiracSpinor::Index a, DiracSpinor::Index b) const {
    const auto ptr = get(a, b);
    return ptr ? *ptr : T{};
  }

  //! Gets pointer to const requested element. If element not present,
  //! returns nullptr. Overload for strings (parses symbol)
  [[nodiscard]] const T *get(const std::string &a, const std::string &b) const {
    const auto [na, ka] = AtomData::parse_symbol(a);
    const auto [nb, kb] = AtomData::parse_symbol(b);
    const auto a_index =
        static_cast<Coulomb::nkIndex>(Angular::nk_to_index(na, ka));
    const auto b_index =
        static_cast<Coulomb::nkIndex>(Angular::nk_to_index(nb, kb));
    const auto map_it = m_data.find(FormIndex(a_index, b_index));
    return (map_it == m_data.cend()) ? nullptr : &(map_it->second);
  }

  //! Provide iterators
  auto begin() { return m_data.begin(); }
  auto end() { return m_data.end(); }
  auto cbegin() const { return m_data.cbegin(); }
  auto cend() const { return m_data.cend(); }

  static std::pair<std::string, std::string> index_to_symbols(nk2Index index) {
    const auto [a, b] = unFormIndex(index);
    const auto [na, ka] = Angular::index_to_nk(int(a));
    const auto [nb, kb] = Angular::index_to_nk(int(b));
    return {AtomData::shortSymbol(na, ka), AtomData::shortSymbol(nb, kb)};
  }

  // private:
  // Converts given set of nkIndex's (in any order) to nk4Index
  [[nodiscard]] static nk2Index FormIndex(nkIndex a, nkIndex b) {
    static_assert(sizeof(nk2Index) == 2 * sizeof(nkIndex));
    static_assert(sizeof(nkIndex) * 8 == 16);
    return (nk2Index)b + ((nk2Index)a << 16);
  }

  [[nodiscard]] static std::pair<nkIndex, nkIndex> unFormIndex(nk2Index index) {
    std::pair<nkIndex, nkIndex> out;
    auto &[a, b] = out;
    b = static_cast<nkIndex>(index);
    // nb: this relies on specific encoding, and may fail?
    a = static_cast<nkIndex>((index - static_cast<nk2Index>(b)) >> 16);
    assert(FormIndex(a, b) == index);
    return out;
  }
};

} // namespace Coulomb