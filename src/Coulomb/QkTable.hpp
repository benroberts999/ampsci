#pragma once
#include "Angular/SixJTable.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "YkTable.hpp"
#include <array>
#include <cstdint>
#include <functional>
#include <string_view>
// #define COULOMB_USE_STD_MAP
#ifdef COULOMB_USE_STD_MAP
#include <unordered_map>
#else
#include "ankerl/unordered_dense.h"
#endif

namespace Coulomb {

//! Symmetry (state index order) for tables.
enum class Symmetry { Qk, Wk, Lk, none };

//! Convert Symmetry to string
inline std::string_view to_string(Symmetry s) {
  switch (s) {
  case Symmetry::Qk:
    return "Qk";
  case Symmetry::Wk:
    return "Wk";
  case Symmetry::Lk:
    return "Lk";
  case Symmetry::none:
    return "none";
  }
  return "";
}

//! Data type used to store integrals
using Real = double;

//! index type for set of 4 orbitals {nk,nk,nk,nk} -> nk4Index [max n: 256]
//! @details @warning Requires 0<n<256, see \ref Angular::nk_to_index
using nk4Index = uint64_t;

//! index type for each {nk} (orbital) [max n: 256]
//! @details @warning Requires 0<n<256, see \ref Angular::nk_to_index
using nkIndex = uint16_t;

static_assert(sizeof(nkIndex) == sizeof(DiracSpinor::Index));

//! Hashmap type used to store integrals (one per k).
//! @details Defaults to ankerl::unordered_dense::map (a fast, densely-stored,
//! open-addressing hashmap). Define COULOMB_USE_STD_MAP at compile time to fall
//! back to std::unordered_map. The two share the API used here, so swapping is
//! transparent to the rest of the code.
#ifdef COULOMB_USE_STD_MAP
using QkMap = std::unordered_map<nk4Index, Real>;
#else
using QkMap = ankerl::unordered_dense::map<nk4Index, Real>;
#endif

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
  Base class template to store Coulomb integrals, and similar. 
  3 specific cases (by template instantiation), account for specific symmetrs.

  @details
  Mostly, a wrapper for std::unordered_map.
  Stores data in a std::vector of maps, each element (map) for each k
  (multipolarity).
  The symmetry is taken into account by defining a "Normal Ordering" of set of
  orbitals {a,b,c,d}, which is defined as the "smallest" arrangment of equivilant indices. The details depend on the specific symmetry in question:

  ### Qk symmetry:
  
  - {abcd} = cbad = adcb = cdab = badc = bcda = dabc = dcba.
  - Normal ordering: "smallest" of above 8 options
  - [ Also assumes Coulomb selection rules, for @ref fill(), and @ref P(), @ref W() ]
  
  ### Wk symmetry (same as g):

  - {abcd} = badc = cdab = dcba
  
  ### Lk symmetry:
  
  - {abcd} = badc

  ### Nk symmetry:

  - No symmetry assumed; each integral treated as unique.

  @note Hashtable lookup is significant bottleneck. Should look into better hashmaps.

  @warning
  Requires 0<n<256, see \ref Angular::nk_to_index

  @note The internal summations in @ref P(), @ref P2(), @ref W(), and @ref g()
  depend on the selection rules, and @p S:
  for @ref Symmetry::Qk, we assume Coulomb angular selection rules, 
  including parity rule, (via @ref k_minmax_Q); 
  for all other symmetries no explicit selection rule is assumed.
  The @ref fill() and @ref fill_if() functions also make this assumption;
  these functions are only available for the @ref Symmetry::Qk symmetry, and
  in that case, assume the regular Coulomb selection rules.
*/
template <Symmetry S>
class CoulombTable {

private:
  // each vector element corresponds to a 'k'
  std::vector<QkMap> m_data{};

public:
  // 'Rule of zero' (except virtual destructor)
  // virtual ~CoulombTable() = default;

  /*!
    @brief Fill table with all non-zero Coulomb Q integrals from a YkTable.
    @details
    Fills the table with Q^k_abcd for all orbitals in @p basis, exploiting
    the 8-fold Qk symmetry to avoid redundant calculations. Uses a precomputed
    @p yk for the radial integrals.

    Filling proceeds in four stages:
    1. Count non-zero integrals per k (for map reservation)
    2. Reserve map capacity
    3. Initialise all entries to zero (parallelised over k)
    4. Fill values in parallel (safe since no new map insertions occur)

    @param basis   Basis set of orbitals.
    @param yk      Precomputed Yk table for radial integrals.
    @param k_cut   Maximum multipolarity k. Set to <=0 for no cut-off.
    @param print   If true, prints timing and summary.

    @note Only valid for QkTable ( @ref Symmetry::Qk ) (enforced via static_assert).
    @warning Does not update existing entries; only adds new ones.
  */
  void fill(const std::vector<DiracSpinor> &basis, const YkTable &yk,
            int k_cut = -1, bool print = true);

  /*!
    @brief Fill table with Coulomb Q integrals satisfying a selection rule.
    @details
    Same as default @ref fill(basis, yk, k_cut, print), but only computes and stores
    integrals for which @p SelectionFunction returns true. Useful for
    restricting to a particular subset of integrals (e.g., CI- or MBPT-relevant).

    Uses the same four-stage filling strategy as fill(basis, yk, ...).

    @param basis               Basis set of orbitals.
    @param yk                  Precomputed Yk table for radial integrals.
    @param SelectionFunction   Returns true for integrals to include.
    @param k_cut               Maximum multipolarity k. Set to <=0 for no cut-off.
    @param print               If true, prints timing and summary.

    @note Only valid for QkTable ( @ref Symmetry::Qk ) (enforced via static_assert).
    @warning Does not update existing entries; only adds new ones.
  */
  void fill_if(const std::vector<DiracSpinor> &basis, const YkTable &yk,
               const SelectionRules &SelectionFunction, int k_cut = -1,
               bool print = true);

  /*!
    @brief Fill table using a general CoulombFunction and selection rules.
    @details
    Fills the table with Fk(k,a,b,c,d) for all orbitals in @p basis,
    restricted to entries for which @p Fk_SR returns true, accounting for
    the symmetry of the table. Unlike the YkTable overloads, this version
    is general-purpose and does not assume Qk angular selection rules.

    Uses the same four-stage filling strategy as fill(basis, yk, ...).

    @param basis    Basis set of orbitals.
    @param Fk       Function to compute the integral value.
    @param Fk_SR    Selection rule; only integrals passing this are stored.
    @param k_cut    Maximum multipolarity k. Set to <=0 for no cut-off.
    @param print    If true, prints timing and summary.

    @warning Does not update existing entries; only adds new ones.
  */
  void fill(const std::vector<DiracSpinor> &basis, const CoulombFunction &Fk,
            const SelectionRules &Fk_SR, int k_cut = -1, bool print = true);

  /*!
    @brief Re-calculate all existing table entries using a CoulombFunction.
    @details
    Updates all integrals currently stored in the table by re-evaluating
    Fk(k,a,b,c,d) using the provided @ref CoulombFunction() for each entry. 
    Unlike @ref fill(), this recalculates all
    existing entries and does not add new ones. Useful for iterating.

    @warning Note: does not calculate any new integrals; 
    only re-calculates existing ones. 
    Table must have corect size/shape before calling.

    Damping: With \f$\eta = \text{damp}\f$, integrals are updated as:
    \f[
      Q^k_{abcd} \to \eta Q^k_{abcd} + (1-\eta) Fk(k,a,b,c,d)
    \f]
    where @p Fk is a @ref CoulombFunction. \f$\eta=0\f$ means no damping.
    Must be between 0 and 1 (0 is allowed, 1 is not).

    @param basis   Basis set of orbitals.
    @param Fk      Function to compute the integral value.
    @param damp    Damping factor, \f$\eta \in [0,1)\f$. 0 means no damping.
    @param print   If true, prints a progress bar.
    @param filter  Optional. If set, entries for which filter(a,b,c,d) returns
                   false are left unchanged (skipped before the k-loop, so no
                   lookup/recalculation is done for them). Empty => update all.
  */
  void update(const std::vector<DiracSpinor> &basis, const CoulombFunction &Fk,
              double damp, bool print = true,
              const std::function<bool(const DiracSpinor &, const DiracSpinor &,
                                       const DiracSpinor &,
                                       const DiracSpinor &)> &filter = {});

  //! Gives arrow access to all underlying vector<unordered_map> functions
  auto operator->() { return &m_data; }

  //! For testing: prints details of coulomb integrals stored
  void summary() const;

  //! Returns number of stored integrals
  std::size_t count() const;

  //!
  /*! Maximum k stored in table. If table is empty, returns '-1'
    @details 
    @note Based on _size_ of map only; doesn't check for values.
    May sometimes have a map with no non-zero integrals for a given k. 
    Guarenteed to be no integrals with k larger than this; 
    no guarentee that there are non-zero integrals at or below this k.
  */
  int max_k() const { return int(m_data.size()) - 1; }

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

  /*! Counts number of non-zero Coulomb integrals, accounting for symmetry
   @details

   @note Specifically assumes Qk selection rules, and therefore only works for Qk.

   In theory, easily updated for general symmetry/selection rules, 
   but this is usually the bottle-neck!
  */
  std::array<std::size_t, 4>
  count_non_zero_integrals(const std::vector<DiracSpinor> &basis,
                           std::size_t max_k = 99, double eF = 0.0) const;

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

  /*!
    @brief Exchange integral P^k_abcd = (2k+1) sum_l {6j} Q^l_abdc.
    @details
    Computes the exchange part of the antisymmetrised W = Q + P integral.
    Optionally uses a precomputed 6J table for faster evaluation.
    @note For @ref Symmetry::Qk, l is iterated using Coulomb angular selection
          rule bounds (k_minmax_Q); for all other symmetries, l runs over
          [0, max_k()]. See @ref CoulombTable.
  */
  Real P(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;

  /*!
    @brief Exchange integral P^k_abcd = (2k+1) sum_l {6j} Q^l_abdc.
    @details
    As @ref P(int, const DiracSpinor&, ...) but takes @ref nkIndex arguments.
    @note For @ref Symmetry::Qk, l is iterated using Coulomb angular selection
          rule bounds (k_minmax_Q); for all other symmetries, l runs over
          [0, max_k()]. See @ref CoulombTable.
  */
  Real P(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
         const Angular::SixJTable *const sj = nullptr) const;

  /*!
    @brief Exchange integral P^k_abcd with effective Coulomb screening.
    @details
    As @ref P(), but each l term is weighted by a screening factor fk[l].
    @note For @ref Symmetry::Qk, l is iterated using Coulomb angular selection
          rule bounds (k_minmax_Q); for all other symmetries, l runs over
          [0, max_k()]. See @ref CoulombTable.
  */
  Real P2(int k, const DiracSpinor &a, const DiracSpinor &b,
          const DiracSpinor &c, const DiracSpinor &d,
          const Angular::SixJTable &sj, const std::vector<double> &fk) const;

  /*!
    @brief Antisymmetrised integral W^k_abcd = Q^k_abcd + P^k_abcd.
    @details
    Returns Q + P. Optionally uses a precomputed 6J table.
    @note l iteration in P() depends on symmetry S; see @ref P() and @ref CoulombTable.
  */
  Real W(int k, const DiracSpinor &a, const DiracSpinor &b,
         const DiracSpinor &c, const DiracSpinor &d,
         const Angular::SixJTable *const sj = nullptr) const;

  /*!
    @brief Antisymmetrised integral W^k_abcd = Q^k_abcd + P^k_abcd.
    @details
    As @ref W(int, const DiracSpinor&, ...) but takes @ref nkIndex arguments.
    @note l iteration in P() depends on symmetry S; see @ref P() and @ref CoulombTable.
  */
  Real W(int k, nkIndex a, nkIndex b, nkIndex c, nkIndex d,
         const Angular::SixJTable *const sj = nullptr) const;

  /*!
    @brief Full matrix element g_abcd with explicit magnetic quantum numbers.
    @details
    Sums over k, weighted by 3j symbols, using the stored Q^k values.
    @note For @ref Symmetry::Qk, k is iterated using Coulomb angular selection
          rule bounds (step 2); for all other symmetries, k runs over
          [0, max_k()] with step 1. See @ref CoulombTable.
  */
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
  void write(const std::string &fname, bool verbose = true) const;
  //! Reads coulomb integrals to disk. Returns false if none read in
  bool read(const std::string &fname, bool verbose = true);

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

  //! Breaks nk4Index back into {ia,ib,ic,id}. Not often used.
  std::array<nkIndex, 4> UnFormIndex(const nk4Index &index) const;

private:
  // Returns the "NormalOrder" nk4Index for given set. Implemented by
  // specific symmetries
  static inline nk4Index NormalOrder_impl(nkIndex a, nkIndex b, nkIndex c,
                                          nkIndex d);

  // True if m_data[k].size() already equals the expected count for every
  // k in [0, count.size()). Used by fill()/fill_if() to skip re-filling a
  // table that already contains all the required entries.
  bool already_filled(const std::vector<std::size_t> &count) const {
    for (auto k = 0ul; k < count.size(); ++k) {
      if (count[k] == 0)
        continue;
      if (k >= m_data.size() || m_data[k].size() < count[k])
        return false;
    }
    return true;
  }
};

/*!
  @brief Coulomb integral table with 8-fold Qk symmetry.
  @details
  Stores Q^k_abcd assuming the symmetry:

   - {abcd} = cbad = adcb = cdab = badc = bcda = dabc = dcba.

  Normal ordering is the lexicographically smallest arrangement among
  these 8 equivalents.

  @note See @ref CoulombTable note regarding selection rules assumptions
*/
using QkTable = CoulombTable<Symmetry::Qk>;

/*!
  @brief Coulomb integral table with 4-fold Wk symmetry.
  @details
  Stores integrals assuming the symmetry:
  {abcd} = badc = cdab = dcba.
*/
using WkTable = CoulombTable<Symmetry::Wk>;

/*!
  @brief Coulomb integral table with 2-fold Lk symmetry.
  @details
  Stores integrals assuming the symmetry:
  {abcd} = badc.
*/
using LkTable = CoulombTable<Symmetry::Lk>;

/*!
  @brief Coulomb integral table with no assumed symmetry.
  @details
  Each {k,a,b,c,d} is treated as a unique entry; no normal ordering is applied.
*/
using NkTable = CoulombTable<Symmetry::none>;

//==============================================================================
//==============================================================================
//==============================================================================

/*!
  @brief Estimate memory required for a QkTable for a given basis.
  @details
  Given the @p basis_string, counts the number of non-zero Coulomb integrals 
  using Qk symmetry and angular selection rules.
  Reports the estimated memory usage in GB for three subsets:
  - All possible Q^k integrals (e.g., structure radiation)
  - Excited-only integrals (all four orbitals above the Fermi level) [CI]
  - Core-Excited integrals (1 or 2 orbitals below Fermi level) [MBPT]

  - And @ref CoulombTable::count_non_zero_integrals() to count integrals
  - Uses @ref estimate_memory_GB() to convert integral counts to memory estimates.

  @note count_non_zero_integrals() takes non-trivial time for large basis

  @param basis_string   Basis specification string (e.g., "20spdf").
  @param core_string    Core configuration string (e.g., "[Xe]").
                        Used to determine the Fermi level.
                        Pass an empty string for no core.
  @param k_cut          Maximum multipolarity k to consider.
*/
void estimate_memory_usage(const std::string &basis_string,
                           const std::string &core_string, int k_cut = 999);

//! Estimate memory required for a Qk table
double estimate_memory_GB(std::size_t number_of_integrals,
                          double efficiency = 0.65);

} // namespace Coulomb

#include "QkTable.ipp"
