#pragma once
#include "Angular/SixJTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <array>
#include <cstdint>
#include <unordered_map>

namespace Coulomb {

//! Symmetry (state index order) for tables.
enum class Symmetry { Qk, Wk, Lk, none };
//! Data type used to store integrals
using Real = double;
//! index type for set of 4 orbitals {nk,nk,nk,nk} -> nk4Index
using nk4Index = uint64_t;
//! index type for each {nk} (orbital)
using nkIndex = uint16_t;

//==============================================================================
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
template <Symmetry S> class CoulombTable {

private:
  // each vector element corresponds to a 'k'
  std::vector<std::unordered_map<nk4Index, Real>> m_data{};

public:
  // 'Rule of zero' (except virtual destructor)
  // virtual ~CoulombTable() = default;

  // XXX Should only be for QkTable
  //! Takes a constructed YkTable, and fills Coulomb table with all possible
  //! non-zero Qk elements, accounting for symmetry (only really makes sense
  //! for QkTable), up to maximum k, k_cut (set k_cut to <=0 to use all k)
  void fill(const std::vector<DiracSpinor> &basis, const YkTable &yk,
            int k_cut = -1);

  //! Gives arrow access to all underlying vector<unordered_map> functions
  auto operator->() { return &m_data; }

  //! For testing: prints details of coulomb integrals stored
  void summary() const;
  //! Returns number of stored integrals
  int count() const;

  //! adds a new value. Note: does nothing if {k,a,b,c,d} already exists
  // XXX Update to return {it,bool}?
  void add(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d, Real value);
  //! Overload for when 'nk4Index' already known
  void add(int k, nk4Index index, Real value);

  //! Updates value in table. If not present, adds new value
  void update(int k, const DiracSpinor &a, const DiracSpinor &b,
              const DiracSpinor &c, const DiracSpinor &d, Real value);
  void update(int k, nk4Index index, Real value);

  //! Checks if given {k,a,b,c,d} is in the table
  bool contains(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d) const;
  bool contains(int k, nk4Index index) const;

  //! Returns true if table is empty
  bool emptyQ() const {
    return (m_data.size() == 0);
    // nb: possible to be empty if size > 0... (but doesn't matter..)
  }

  //! Retrieve a stored Q. If not present, returns 0. (Returns exactly as
  //! stored in table.)
  Real Q(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;
  //! Overload for when 'nk4Index' already known
  Real Q(int k, nk4Index index) const;

  //! Returns 'R', defined via: R := Q / (angular_coef)
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;

  //! 'Exchange-only', defined via W = Q + P. Optionally, takes pointer to 6J
  //! table (faster eval of 6J symbols)
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;

  //! 'Exchange-only', defined via W = Q + P. Optionally, takes pointer to 6J
  //! table (faster eval of 6J symbols) - with screening
  Real P2(int k, const DiracSpinor &a, const DiracSpinor &b,
          const DiracSpinor &c, const DiracSpinor &d,
          const Angular::SixJTable &sj, const std::vector<double> &fk) const;

  //! W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc. Optionally, takes
  //! pointer to 6J table (faster eval of 6J symbols)
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;

  //! Creates single 'nk4Index' corresponding to 'NormalOrder' symmetry of
  //! {a,b,c,d}
  nk4Index NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
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

  // //! Operator overload: scale by constant (nb: *= not *, since can't rely on
  // //! polymorphism for return type..)
  // void operator*=(double x) {
  //   for (auto &q_k : m_data) {
  //     for (auto &[key, qk_abcd] : q_k) {
  //       qk_abcd *= x;
  //     }
  //   }
  // }
  // //! Operator overload: add another table (nb: += not +, since can't rely on
  // //! polymorphism for return type..)
  // void operator+=(const CoulombTable &Qk) {
  //   int k = 0;
  //   for (auto &q_k : m_data) {
  //     for (auto &[key, qk_abcd] : q_k) {
  //       // nb: cannot assume each map's entries in exact same order, so need
  //       // 'find'
  //       qk_abcd += Qk.Q(k, key);
  //     }
  //     ++k;
  //   }
  // }

  // auto begin() { return m_data.begin(); }
  // auto end() { return m_data.end(); }
  // auto cbegin() { return m_data.cbegin(); }
  // auto cend() { return m_data.cend(); }
  // auto begin(std::size_t k) { return m_data.at(k).begin(); }
  // auto end(std::size_t k) { return m_data.at(k).end(); }
  // auto cbegin(std::size_t k) { return m_data.at(k).cbegin(); }
  // auto cend(std::size_t k) { return m_data.at(k).cend(); }

private:
  // Creates single 'nk4Index', WITHOUT accounting for 'NormalOrder'. Can be
  // used to check if {a,b,c,d} are already in 'NormalOrder'
  nk4Index CurrentOrder(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const;

  // Converts given set of nkIndex's (in any order) to nk4Index
  static nk4Index FormIndex(nkIndex a, nkIndex b, nkIndex c, nkIndex d);

  // Breaks nk4Index back into {ia,ib,ic,id}. Not used.
  std::array<nkIndex, 4> UnFormIndex(const nk4Index &index) const;

  // // Virtual: returns the "NormalOrder" nk4Index for given set. Overriden by
  // // derived classes.
  // virtual nk4Index NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
  //                                   nkIndex d) const = 0;
  static inline nk4Index NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                          nkIndex d);
};

using QkTable = CoulombTable<Symmetry::Qk>;
using WkTable = CoulombTable<Symmetry::Wk>;
using LkTable = CoulombTable<Symmetry::Lk>;
using NkTable = CoulombTable<Symmetry::none>;

} // namespace Coulomb

#include "QkTable.ipp"
