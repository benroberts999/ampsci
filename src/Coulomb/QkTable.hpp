#pragma once
#include "Angular/SixJTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <array>
#include <cstdint>
#include <functional>
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

static_assert(sizeof(nkIndex) == sizeof(DiracSpinor::Index));

//! Function type for calculating Coulomb(-like) integrals.
//! Takes k and 4 DiracSpinors, returns a double
using CoulombFunction =
    std::function<double(int, const DiracSpinor &a, const DiracSpinor &b,
                         const DiracSpinor &c, const DiracSpinor &d)>;

//! Function type for determining Coulomb(-like) integral selection rules.
//! Takes k and 4 DiracSpinors, returns a true if integral is allowed.
using SelectionRules =
    std::function<bool(int, const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d)>;

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
template <Symmetry S>
class CoulombTable {

private:
  // each vector element corresponds to a 'k'
  std::vector<std::unordered_map<nk4Index, Real>> m_data{};

public:
  // 'Rule of zero' (except virtual destructor)
  // virtual ~CoulombTable() = default;

  //! Takes a constructed YkTable, and fills Coulomb table with all possible
  //! non-zero Qk Coulomb integrals, accounting for symmetry (only allowed
  //! for QkTable), up to maximum k, k_cut (set k_cut to <=0 to use all k)
  void fill(const std::vector<DiracSpinor> &basis, const YkTable &yk,
            int k_cut = -1, bool print = true);

  //! Takes a constructed YkTable, and fills Coulomb table with all possible
  //! non-zero Qk Coulomb integrals that are allowed by the SelectionFunction,
  //! accounting for symmetry (only allowed for QkTable), up to maximum k,
  //! k_cut (set k_cut to <=0 to use all k)
  void fill_if(const std::vector<DiracSpinor> &basis, const YkTable &yk,
               const SelectionRules &SelectionFunction, int k_cut = -1,
               bool print = true);

  void fill(const std::vector<DiracSpinor> &basis, const CoulombFunction &Fk,
            const SelectionRules &Fk_SR, int k_cut = -1, bool print = true);

  //! Gives arrow access to all underlying vector<unordered_map> functions
  auto operator->() { return &m_data; }

  //! For testing: prints details of coulomb integrals stored
  void summary() const;
  //! Returns number of stored integrals
  std::size_t count() const;

  //! adds a new value. Note: does nothing if {k,a,b,c,d} already exists
  void add(int k, const DiracSpinor &a, const DiracSpinor &b,
           const DiracSpinor &c, const DiracSpinor &d, Real value);
  //! adds a new value. Note: does nothing if {k,a,b,c,d} already exists
  void add(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d, Real value);
  //! adds a new value. Note: does nothing if {k,a,b,c,d} already exists
  void add(int k, nk4Index index, Real value);

  //! Updates value in table. If not present, adds new value
  void update(int k, const DiracSpinor &a, const DiracSpinor &b,
              const DiracSpinor &c, const DiracSpinor &d, Real value);
  //! Updates value in table. If not present, adds new value
  void update(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d, Real value);
  //! Updates value in table. If not present, adds new value
  void update(int k, nk4Index index, Real value);

  //! Checks if given {k,a,b,c,d} is in the table
  bool contains(int k, const DiracSpinor &a, const DiracSpinor &b,
                const DiracSpinor &c, const DiracSpinor &d) const;
  //! Checks if given {k,a,b,c,d} is in the table
  bool contains(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d) const;
  //! Checks if given {k,a,b,c,d} is in the table
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
  //! Retrieve a stored Q. If not present, returns 0. (Returns exactly as
  //! stored in table.)

  Real Q(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d) const;
  //! Retrieve a stored Q. If not present, returns 0. (Returns exactly as
  //! stored in table.)
  Real Q(int k, nk4Index index) const;

  //! Returns 'R', defined via: R := Q / (angular_coef)
  Real R(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d) const;
  //! Returns 'R', defined via: R := Q / (angular_coef)
  Real R(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d) const;

  //! 'Exchange-only', defined via W = Q + P. Optionally, takes pointer to 6J
  //! table (faster eval of 6J symbols)
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;
  //! 'Exchange-only', defined via W = Q + P. Optionally, takes pointer to 6J
  //! table (faster eval of 6J symbols)
  Real P(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
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
  //! W^k_abcd = Q^k_abcd + \sum_l [k] 6j * Q^l_abdc. Optionally, takes
  //! pointer to 6J table (faster eval of 6J symbols)
  Real W(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
         const Angular::SixJTable *const sj = nullptr) const;

  //! Returns 'g_abcd'
  Real g(const DiracSpinor &a, const DiracSpinor &b, const DiracSpinor &c,
         const DiracSpinor &d, int tma, int tmb, int tmc, int tmd) const;

  //! Creates single 'nk4Index' corresponding to 'NormalOrder' symmetry of
  //! {a,b,c,d}
  nk4Index NormalOrder(const DiracSpinor &a, const DiracSpinor &b,
                       const DiracSpinor &c, const DiracSpinor &d) const;
  //! Creates single 'nk4Index' corresponding to 'NormalOrder' symmetry of
  //! {a,b,c,d}
  nk4Index NormalOrder(nkIndex a, nkIndex b, nkIndex c, nkIndex d) const;

  //! Checks if set {a,b,c,d} are already in NormalOrder
  bool is_NormalOrdered(nkIndex a, nkIndex b, nkIndex c, nkIndex d) const {
    return NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d);
  }
  //! Checks if set {a,b,c,d} are already in NormalOrder
  bool is_NormalOrdered(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const {
    return NormalOrder(a, b, c, d) == CurrentOrder(a, b, c, d);
  }

  //! Writes coulomb integrals to disk
  void write(const std::string &fname) const;
  //! Reads coulomb integrals to disk. Returns false if none read in
  bool read(const std::string &fname);

  //! Directly gets one of the stored elements, given normal-ordered nk4Index
  double *get(int k, nk4Index index);

  //! Creates single 'nk4Index', WITHOUT accounting for 'NormalOrder'. Can be
  //! used to check if {a,b,c,d} are already in 'NormalOrder'
  nk4Index CurrentOrder(const DiracSpinor &a, const DiracSpinor &b,
                        const DiracSpinor &c, const DiracSpinor &d) const;

  //! Creates single 'nk4Index', WITHOUT accounting for 'NormalOrder'. Can be
  //! used to check if {a,b,c,d} are already in 'NormalOrder'
  nk4Index CurrentOrder(nkIndex a, nkIndex b, nkIndex c, nkIndex d) const;

  //! Converts given set of nkIndex's (in any order) to nk4Index
  static nk4Index FormIndex(nkIndex a, nkIndex b, nkIndex c, nkIndex d);

  //! Breaks nk4Index back into {ia,ib,ic,id}. Not used.
  std::array<nkIndex, 4> UnFormIndex(const nk4Index &index) const;

private:
  // Returns the "NormalOrder" nk4Index for given set. Implemented by
  // specific symmetries
  static inline nk4Index NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                          nkIndex d);
};

//! Stores Coulomb(-like) integrals, assuming Q-like symmetry
using QkTable = CoulombTable<Symmetry::Qk>;
//! Stores Coulomb(-like) integrals, assuming W-like symmetry
using WkTable = CoulombTable<Symmetry::Wk>;
//! Stores Coulomb(-like) integrals, assuming L-like symmetry
using LkTable = CoulombTable<Symmetry::Lk>;
//! Stores Coulomb(-like) integrals, assuming no symmetry
using NkTable = CoulombTable<Symmetry::none>;

} // namespace Coulomb

#include "QkTable.ipp"
