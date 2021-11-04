#pragma once
#include "Angular/SixJTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <array>
#include <unordered_map>

namespace Coulomb {

// ! Symmetry (state index order) for tables.
enum class Symmetry { Qk, Wk, Lk, none };

//******************************************************************************
/*!
@brief
Base (pure virtual) class to store Coulomb integrals, and similar. 3 derived
classes, which account for symmetry.
 @details
Mostly, a wrapper for std::unordered_map.
Stores data in a std::vector of maps, each element (map) for each k
(multipolarity).
The symmetry is taken into account by defining a "Normal Ordering" of set of
orbitals {a,b,c,d}. How this is defined depends on the specific symmetry in
question.
 Qk symmetry:
 {abcd} = cbad = adcb = cdab = badc = bcda = dabc = dcba.
 Normal ordering: i = min(a,b,c,d) is first index.
 Two options (second and fourth may be swapped): choose 2nd to be smallest
 Wk symmetry (same as g):
 {abcd} = badc = cdab = dcba
 Lk symmetry:
 {abcd} = badc
*/
class CoulombTable {

public:
  //! Data type used to store integra;s
  using Real = double;
  //! index tpe for set of 4 orbitals {nk,nk,nk,nk} -> BigIndex
  using BigIndex = uint64_t;
  //! index type for each {nk} (orbital)
  using Index = uint16_t;

protected:
  // each vector element corresponds to a 'k'
  std::vector<std::unordered_map<BigIndex, Real>> m_data{};

public:
  // 'Rule of zero' (except virtual destructor)
  virtual ~CoulombTable() = default;

  //! Gives arrow access to all underlying vector<unordered_map> functions
  auto operator-> () { return &m_data; }

  //! For testing: prints details of coulomb integrals stored
  void summary() const;
  //! Returns number of stored integrals
  int count() const;

  //! adds a new value. Note: does nothing if {k,a,b,c,d} already exists
  // XXX Update to return {it,bool}?
  void add(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d, Real value);
  //! Overload for when 'BigIndex' already known
  void add(int k, BigIndex index, Real value);

  //! Updates value in table. If not present, adds new value
  void update(int k, const DiracSpinor &a, const DiracSpinor &b,
              const DiracSpinor &c, const DiracSpinor &d, Real value);
  void update(int k, BigIndex index, Real value);

  //! Checks if given {k,a,b,c,d} is in the table
  bool contains(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d) const;
  bool contains(int k, BigIndex index) const;

  //! Returns true if table is empty
  bool emptyQ() const {
    return (m_data.size() == 0);
    // nb: possible to be empty if size > 0... (but doesn't matter..)
  }

  //! Retrieve a stored Q. If not present, returns 0. (Returns exactly as
  //! stored in table.)
  Real Q(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;
  //! Overload for when 'BigIndex' already known
  Real Q(int k, BigIndex index) const;

  //! Returns 'R', defined via: R := Q / (angular_coef)
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;

  //! 'Exchange-only', defined via W = Q + P. Optionally, takes pointer to 6J
  //! table (faster eval of 6J symbols)
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;

  //! W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc. Optionally, takes
  //! pointer to 6J table (faster eval of 6J symbols)
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;

  //! Creates single 'BigIndex' corresponding to 'NormalOrder' symmetry of
  //! {a,b,c,d}
  BigIndex NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const;

  //! Checks if set {a,b,c,d} are already in NormalOrder
  bool is_NormalOrdered(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const {
    return NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d);
  }

  //! Writes coulomb integrals to disk
  void write(const std::string &fname) const;
  //! Reads coulomb integrals to disk. Returns false if none read in
  bool read(const std::string &fname);

  //! Operator overload: scale by constant (nb: *= not *, since can't rely on
  //! polymorphism for return type..)
  void operator*=(double x) {
    for (auto &q_k : m_data) {
      for (auto &[key, qk_abcd] : q_k) {
        qk_abcd *= x;
      }
    }
  }
  //! Operator overload: add another table (nb: += not +, since can't rely on
  //! polymorphism for return type..)
  void operator+=(const CoulombTable &Qk) {
    int k = 0;
    for (auto &q_k : m_data) {
      for (auto &[key, qk_abcd] : q_k) {
        // nb: cannot assume each map's entries in exact same order, so need
        // 'find'
        qk_abcd += Qk.Q(k, key);
      }
      ++k;
    }
  }

  auto begin() { return m_data.begin(); }
  auto end() { return m_data.end(); }
  auto cbegin() { return m_data.cbegin(); }
  auto cend() { return m_data.cend(); }
  auto begin(std::size_t k) { return m_data.at(k).begin(); }
  auto end(std::size_t k) { return m_data.at(k).end(); }
  auto cbegin(std::size_t k) { return m_data.at(k).cbegin(); }
  auto cend(std::size_t k) { return m_data.at(k).cend(); }

protected:
  // Creates single 'BigIndex', WITHOUT accounting for 'NormalOrder'. Can be
  // used to check if {a,b,c,d} are already in 'NormalOrder'
  BigIndex CurrentOrder(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const;

  // Converts given set of Index's (in any order) to BigIndex
  BigIndex FormIndex(Index a, Index b, Index c, Index d) const;

  // Breaks BigIndex back into {ia,ib,ic,id}. Not used.
  std::array<Index, 4> UnFormIndex(const BigIndex &index) const;

  // Virtual: returns the "NormalOrder" BigIndex for given set. Overriden by
  // derived classes.
  virtual BigIndex NormalOrder_impl(Index a, Index b, Index c,
                                    Index d) const = 0;
};

//******************************************************************************
//! Derived CoulombTable for the Qk symmetry: {abcd} = cbad = adcb = cdab =
//! badc = bcda = dabc = dcba.
class QkTable : public CoulombTable {

public:
  static constexpr Symmetry symmetry = Symmetry::Qk;

  // //! Takes a constructed YkTable, and fills Coulomb table with all
  // possible
  // //! non-zero Qk elements, accounting for symmetry (only really makes
  // sense for
  // //! QkTable)
  // void fill(const YkTable &yk);
  //! Takes a constructed YkTable, and fills Coulomb table with all possible
  //! non-zero Qk elements, accounting for symmetry (only really makes sense
  //! for QkTable)
  void fill(const std::vector<DiracSpinor> &basis, const YkTable &yk);
  void fill_old(const std::vector<DiracSpinor> &basis, const YkTable &yk);

private:
  virtual BigIndex NormalOrder_impl(Index a, Index b, Index c,
                                    Index d) const override final;
};

//******************************************************************************
//! Derived CoulombTable for the Wk symmetry: {abcd} = badc = cdab = dcba
class WkTable : public CoulombTable {

public:
  static constexpr Symmetry symmetry = Symmetry::Wk;

private:
  virtual BigIndex NormalOrder_impl(Index a, Index b, Index c,
                                    Index d) const override final;
};

//******************************************************************************
//! Derived CoulombTable for the Lk symmetry: {abcd} = badc
class LkTable : public CoulombTable {

public:
  static constexpr Symmetry symmetry = Symmetry::Lk;

private:
  virtual BigIndex NormalOrder_impl(Index a, Index b, Index c,
                                    Index d) const override final;
};

//******************************************************************************
//! Derived CoulombTable for NO symmetry: {abcd} = abcd only
class NkTable : public CoulombTable {

public:
  static constexpr Symmetry symmetry = Symmetry::none;

private:
  virtual BigIndex NormalOrder_impl(Index a, Index b, Index c,
                                    Index d) const override final;
};

} // namespace Coulomb
