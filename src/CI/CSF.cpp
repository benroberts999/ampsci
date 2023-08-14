#include "CSF.hpp"
#include "Angular/Angular.hpp"
#include "Physics/AtomData.hpp"
#include <algorithm>
#include <array>
#include <cassert>

namespace CI {

//==============================================================================
CSF2::CSF2(const DiracSpinor &a, const DiracSpinor &b)
    : states(a <= b ? std::array{&a, &b} : std::array{&b, &a}) {}

const DiracSpinor *CSF2::state(std::size_t i) const { return states.at(i); }

bool operator==(const CSF2 &A, const CSF2 &B) {
  // only works because states is sorted
  return A.states == B.states;
}

bool operator!=(const CSF2 &A, const CSF2 &B) {
  // only works because states is sorted
  return !(A.states == B.states);
}

// Returns number of different orbitals between two CSFs
int CSF2::num_different(const CSF2 &A, const CSF2 &B) {
  if (A == B)
    return 0;
  if (A.state(0) == B.state(0) || A.state(1) == B.state(1) || //
      A.state(0) == B.state(1) || A.state(1) == B.state(0))
    return 1;
  return 2;
}

// returns _different_ orbitals, for case where CSFs differ by 1.
// i.e., returns {n,a} where |A> = |B_a^n> (i.e., A has n, but not a)
std::array<const DiracSpinor *, 2> CSF2::diff_1_na(const CSF2 &A,
                                                   const CSF2 &B) {
  assert(num_different(A, B) == 1); // only valid in this case
  if (A.state(0) == B.state(0))
    return {A.state(1), B.state(1)};
  if (A.state(1) == B.state(1))
    return {A.state(0), B.state(0)};
  if (A.state(0) == B.state(1))
    return {A.state(1), B.state(0)};
  if (A.state(1) == B.state(0))
    return {A.state(0), B.state(1)};
  assert(false); // should be unreachable, for testing
}

//==============================================================================
// Forms list of all possible (2-particle) CSFs with given J and parity
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis) {

  assert((parity == 1 || parity == -1) && "Parity must be +/-1");
  assert(twoJ >= 0 && "2*J must not be negative");

  std::vector<CSF2> CSFs;

  // This is always smaller than actual # of CSFs, but good head-start
  CSFs.reserve(cisp_basis.size());

  for (const auto &v : cisp_basis) {
    for (const auto &w : cisp_basis) {

      // Symmetric: only include unique CSFs once
      if (w < v)
        continue;

      // Parity symmetry:
      if (v.parity() * w.parity() != parity)
        continue;

      // J triangle rule (use M=Jz=J):
      if (v.twoj() + w.twoj() < twoJ || std::abs(v.twoj() - w.twoj()) > twoJ)
        continue;

      // identical particles can only give even J (cannot have mv=mw)
      if (v == w && twoJ % 4 != 0)
        continue;

      CSFs.emplace_back(v, w);
    }
  }

  return CSFs;
}

//==============================================================================
// Takes a subset of input basis according to subset_string.
// Only states *not* included in frozen_core_string are included.
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      const std::string &frozen_core_string) {

  // Form 'subset' from {a} in 'basis', if:
  //    a _is_ in subset_string AND
  //    a _is not_ in basis string

  std::vector<DiracSpinor> subset;
  const auto nmaxk_list = AtomData::n_kappa_list(subset_string);
  const auto core_list = AtomData::core_parser(frozen_core_string);

  for (const auto &a : basis) {

    // Check if a is present in 'subset_string'
    const auto nk =
        std::find_if(nmaxk_list.cbegin(), nmaxk_list.cend(),
                     [&a](const auto &tnk) { return a.kappa() == tnk.second; });
    if (nk == nmaxk_list.cend())
      continue;
    // nk is now max n, for given kappa {max_n, kappa}
    if (a.n() > nk->first)
      continue;

    // assume only filled shells in frozen core
    const auto core = std::find_if(
        core_list.cbegin(), core_list.cend(), [&a](const auto &tcore) {
          return a.n() == tcore.n && a.l() == tcore.l;
        });

    if (core != core_list.cend())
      continue;
    subset.push_back(a);
  }
  return subset;
}

} // namespace CI