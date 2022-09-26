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
                                       bool each_freq, bool diagonal_only,
                                       bool print_both, bool radial_int) {

  std::vector<MEdata> res;

  auto AhfsQ = h->name() == "hfs";

  if (&b_orbs != &a_orbs)
    print_both = true;

  if (h->freqDependantQ() && !each_freq) {
    h->updateFrequency(omega);
  }

  if (dV && !each_freq) {
    dV->solve_core(omega);
  }

  for (const auto &Fb : b_orbs) {
    for (const auto &Fa : a_orbs) {

      if (h->isZero(Fa.kappa(), Fb.kappa()))
        continue;
      if (diagonal_only && Fb != Fa)
        continue;
      if (!print_both && Fb > Fa)
        continue;

      const auto ww = std::abs(Fa.en() - Fb.en());
      if (each_freq && h->freqDependantQ()) {
        h->updateFrequency(ww);
      }
      if (each_freq && dV) {
        if (dV->get_eps() > 1.0e-5)
          dV->clear();
        dV->solve_core(ww);
      }

      // Special case: HFS A:
      const auto a = radial_int ? 1.0 / h->angularF(Fa.kappa(), Fb.kappa()) :
                     AhfsQ ? DiracOperator::HyperfineA::convertRMEtoA(Fa, Fb) :
                             1.0;

      const auto hab = h->reducedME(Fa, Fb) * a;
      const auto dv1 = dV ? dV->dV1(Fa, Fb) * a : 0.0;
      const auto dv = dV ? dV->dV(Fa, Fb) * a : 0.0;

      auto x = res.emplace_back(
          MEdata{Fa.shortSymbol(), Fb.shortSymbol(), hab, dv1, dv});
    }
  }

  return res;
}

} // namespace ExternalField
