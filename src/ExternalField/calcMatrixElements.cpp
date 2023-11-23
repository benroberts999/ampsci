#include "calcMatrixElements.hpp"
#include "CorePolarisation.hpp"
#include "DiagramRPA.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "HF/HartreeFock.hpp"
#include "TDHF.hpp"
#include "TDHFbasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/color.hpp"
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

//==============================================================================
// Required to set omega for freq. dependent operators initially
Coulomb::meTable<double> me_table(const std::vector<DiracSpinor> &a_orbs,
                                  const std::vector<DiracSpinor> &b_orbs,
                                  const DiracOperator::TensorOperator *h,
                                  const CorePolarisation *dV,
                                  const MBPT::StructureRad *srn,
                                  std::optional<double> omega) {

  Coulomb::meTable<double> h_ab;

  const auto a_is_b = &a_orbs == &b_orbs;

  for (const auto &a : a_orbs) {
    for (const auto &b : b_orbs) {
      if (b < a && a_is_b)
        continue;
      if (h->isZero(a, b))
        continue;
      h_ab.add(a, b, 0.0);
      if (a != b) {
        h_ab.add(b, a, 0.0);
      }
    }
  }

#pragma omp parallel for
  for (std::size_t i = 0; i < a_orbs.size(); ++i) {
    const auto &a = a_orbs[i];
    for (const auto &b : b_orbs) {

      if (b < a && a_is_b)
        continue;

      if (h->isZero(a, b))
        continue;

      const auto ww = omega ? *omega : std::abs(b.en() - a.en()); // abs?

      const auto tab = h->reducedME(a, b);
      const auto dv = dV ? dV->dV(a, b) : 0.0;
      const auto sr = srn ? srn->srn(h, a, b, ww, dV).second : 0.0;

      const auto me = tab + dv + sr;

      h_ab.update(a, b, me);
      if (a != b) {
        h_ab.update(b, a, h->symm_sign(a, b) * me);
      }
    }
  }
  return h_ab;
}

//==============================================================================
std::unique_ptr<CorePolarisation>
make_rpa(const std::string &method, const DiracOperator::TensorOperator *h,
         const HF::HartreeFock *vhf, bool print,
         const std::vector<DiracSpinor> &basis, const std::string &identity) {

  // Parse method for RPA:
  auto rpa_method = ExternalField::ParseMethod(method);

  if (rpa_method == ExternalField::Method::Error) {
    fmt2::styled_print(fg(fmt::color::red), "\nError 148: ");
    fmt::print(
        "RPA method {} not found - check spelling? Defaulting to NO rpa\n",
        method);
    rpa_method = ExternalField::Method::none;
  }
  const auto rpaQ = rpa_method != ExternalField::Method::none;

  // do RPA:
  std::unique_ptr<ExternalField::CorePolarisation> dV{nullptr};
  if (rpaQ && print)
    std::cout << "Including RPA: ";
  if (rpa_method == ExternalField::Method::TDHF) {
    if (print)
      std::cout << "TDHF method\n";
    dV = std::make_unique<ExternalField::TDHF>(h, vhf);
  } else if (rpa_method == ExternalField::Method::basis) {
    if (print)
      std::cout << "TDHF/basis method (" << DiracSpinor::state_config(basis)
                << ")\n";
    dV = std::make_unique<ExternalField::TDHFbasis>(h, vhf, basis);
  } else if (rpa_method == ExternalField::Method::diagram) {
    if (print)
      std::cout << "diagram method (" << DiracSpinor::state_config(basis)
                << ")\n";
    dV = std::make_unique<ExternalField::DiagramRPA>(h, basis, vhf, identity);
  }
  return dV;
}

} // namespace ExternalField
