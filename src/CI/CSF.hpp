#pragma once
#include "Wavefunction/DiracSpinor.hpp"
#include <array>
#include <vector>

namespace CI {

//==============================================================================
//! Very basic two-electron CSF.
class CSF2 {
public:
  // nb: array of states is always sorted
  std::array<const DiracSpinor *, 2> states;

  CSF2(const DiracSpinor &a, const DiracSpinor &b);

  const DiracSpinor *state(std::size_t i) const;

  friend bool operator==(const CSF2 &A, const CSF2 &B);
  friend bool operator!=(const CSF2 &A, const CSF2 &B);

  //! Returns number of different orbitals between two CSFs
  static int num_different(const CSF2 &A, const CSF2 &B);

  //! returns _different_ orbitals, for case where CSFs differ by 1.
  //! i.e., returns {n,a} where |V> = |X_a^n> (i.e., V has n, but not a)
  static std::array<const DiracSpinor *, 2> diff_1_na(const CSF2 &V,
                                                      const CSF2 &X);

  int parity() const { return states[0]->parity() * states[1]->parity(); }
};

//==============================================================================
//! Forms list of all possible (2-particle) CSFs with given J and parity
std::vector<CSF2> form_CSFs(int twoJ, int parity,
                            const std::vector<DiracSpinor> &cisp_basis);

//==============================================================================
//! Takes a subset of input basis according to subset_string.
//! Only states *not* included in frozen_core_string are included.
std::vector<DiracSpinor> basis_subset(const std::vector<DiracSpinor> &basis,
                                      const std::string &subset_string,
                                      const std::string &frozen_core_string);

} // namespace CI