#include "calcMatrixElements.hpp"
#include "CorePolarisation.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <iostream>
#include <string>
#include <vector>

namespace ExternalField {

std::vector<MEdata> calcMatrixElements(const std::vector<DiracSpinor> &b_orbs,
                                       const std::vector<DiracSpinor> &a_orbs,
                                       DiracOperator::TensorOperator *const h,
                                       CorePolarisation *const dV, double omega,
                                       bool each_freq, bool diagonal,
                                       bool off_diagonal, bool calculate_both) {

  std::vector<MEdata> res;

  const auto pi = h->parity();

  if (&b_orbs != &a_orbs)
    calculate_both = true;

  if (h->freqDependantQ() && !each_freq) {
    h->updateFrequency(omega);
  }

  if (dV && !each_freq) {
    dV->solve_core(omega);
  }

  // First, diagonal:
  if (diagonal) {

    if (each_freq && h->freqDependantQ()) {
      h->updateFrequency(0.0);
    }
    if (each_freq && dV) {
      if (dV->get_eps() > 1.0e-5)
        dV->clear();
      dV->solve_core(0.0);
    }

    for (std::size_t ia = 0; ia < a_orbs.size(); ia++) {
      const auto &Fa = b_orbs.at(ia);

      if (h->isZero(Fa.kappa(), Fa.kappa()))
        continue;

      const auto hab = h->reducedME(Fa, Fa);
      const auto dv = dV ? dV->dV(Fa, Fa) : 0.0;

      const auto w = 0.0;
      res.emplace_back(MEdata{Fa.shortSymbol(), Fa.shortSymbol(), w, hab, dv});
    }
  }

  // Then, off-diagonal:
  if (off_diagonal) {
    for (std::size_t ib = 0; ib < b_orbs.size(); ib++) {
      const auto &Fb = b_orbs.at(ib);
      for (std::size_t ia = 0; ia < a_orbs.size(); ia++) {
        const auto &Fa = b_orbs.at(ia);

        if (Fa == Fb)
          continue;
        if (h->isZero(Fa.kappa(), Fb.kappa()))
          continue;

        // Ensure even-parity state on right for odd-parity operators
        if (pi == -1) {
          if (!calculate_both && Fb.parity() == -1)
            continue;
        } else {
          if (!calculate_both && ib > ia)
            continue;
        }

        const auto ww = std::abs(Fa.en() - Fb.en());
        if (each_freq && h->freqDependantQ()) {
          h->updateFrequency(ww);
        }
        if (each_freq && dV) {
          if (dV->get_eps() > 1.0e-5)
            dV->clear();
          dV->solve_core(ww);
        }

        const auto hab = h->reducedME(Fa, Fb);
        const auto dv = dV ? dV->dV(Fa, Fb) : 0.0;

        const auto w = Fa.en() - Fb.en();
        res.emplace_back(
            MEdata{Fa.shortSymbol(), Fb.shortSymbol(), w, hab, dv});
      }
    }
  }

  return res;
}

} // namespace ExternalField
