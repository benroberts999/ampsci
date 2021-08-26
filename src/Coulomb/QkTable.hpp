#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <array>
#include <unordered_map>

namespace Coulomb {

// ! Symmetry (state index order) for tables.
enum class Symmetry { Qk, Wk, none };

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
 Wk symmetry:
 {abcd} = badc = cdab = dcba
*/
class CoulombTable {

public:
  //! Data type used to store integra;s
  using Real = double;
  //! index tpe for set of 4 orbitals {nk,nk,nk,nk} -> BigIndex
  using BigIndex = uint64_t;
  //! index type for each {nk} (orbital)
  using Index = uint16_t;
  // Type for set of four Index's: {nk,nk,nk,nk}
  using IndexSet = std::array<Index, 4>;

private:
  // each vector element corresponds to a 'k'
  std::vector<std::unordered_map<BigIndex, Real>> m_data{};

public:
  // 'Rule of zero' (except virtual destructor)
  virtual ~CoulombTable() = default;

  //! Gives arrow access to all underlying vector<unordered_map> functions
  auto operator-> () { return &m_data; }

  //! For testing: prints details of coulomb integrals stored
  void count() const;

  //! adds a new value. Note: does nothing if {k,a,b,c,d} already exists
  // XXX Update to return {it,bool}?
  void add(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d, Real value);

  //! Overload for when 'BigIndex' already known
  void add(int k, BigIndex index, Real value);

  //! Updates value in table. If not present, adds new value
  void update(int k, const DiracSpinor &a, const DiracSpinor &b,
              const DiracSpinor &c, const DiracSpinor &d, Real value);

  //! Checks if given {k,a,b,c,d} is in the table
  bool contains(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d) const;

  //! Retrieve a stored Q. If not present, returns 0. (Returns exactly as stored
  //! in table.)
  Real Q(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;

  //! Returns 'R', defined via: R := Q / (angular_coef)
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;

  //! 'Exchange-only', defined via W = Q + P
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;

  //! W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;

  //! Creates single 'BigIndex' corresponding to 'NormalOrder' symmetry of
  //! {a,b,c,d}
  BigIndex NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const;

  //! Checks if set {a,b,c,d} are already in NormalOrder
  bool is_NormalOrdered(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const {
    return NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d);
  }

  //! Takes a constructed YkTable, and fills Coulomb table with all possible
  //! non-zero Qk elements, accounting for symmetry (only really makes sense for
  //! QkTable)
  void fill(const YkTable &yk);

  //! Writes coulomb integrals to disk
  void write(const std::string &fname) const;
  //! Reads coulomb integrals to disk. Returns false if none read in
  bool read(const std::string &fname);

protected:
  // Creates single 'BigIndex', WITHOUT accounting for 'NormalOrder'. Can be
  // used to check if {a,b,c,d} are already in 'NormalOrder'
  BigIndex CurrentOrder(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const;

  // Converts given IndexSet (in any order) to BigIndex
  BigIndex FormIndex(const IndexSet &set) const;

  // Breaks BigIndex back into {ia,ib,ic,id}. Not used.
  IndexSet UnFormIndex(const BigIndex &index) const;

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
//! Derived CoulombTable for NO symmetry: {abcd} = abcd only
class NkTable : public CoulombTable {

public:
  static constexpr Symmetry symmetry = Symmetry::none;

private:
  virtual BigIndex NormalOrder_impl(Index a, Index b, Index c,
                                    Index d) const override final;
};

} // namespace Coulomb
