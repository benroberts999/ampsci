#include "ExternalField.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/DiracSpinor.hpp"
#include <algorithm>
#include <vector>

//******************************************************************************
ExternalField::ExternalField(const DiracOperator *const h,
                             const std::vector<DiracSpinor> &core,
                             const double omega)
    : m_h(h), p_core(&core), m_omega(omega), m_rank(h->rank()),
      m_pi(h->parity()), m_imag(h->imaginaryQ()) {}

//******************************************************************************
const std::vector<DiracSpinor> &
ExternalField::get_dPsis(const DiracSpinor &phic, dPsiType XorY) {
  auto index = static_cast<std::size_t>(
      std::find(p_core->cbegin(), p_core->cend(), phic) - p_core->cbegin());
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//******************************************************************************
const DiracSpinor &ExternalField::get_dPsi_x(const DiracSpinor &phic,
                                             dPsiType XorY, const int kappa_x) {
  const auto &dPsis = get_dPsis(phic, XorY);
  auto match_kappa_x = [=](const DiracSpinor &phi) { return phi.k == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//******************************************************************************
void ExternalField::solve_TDHFcore() {}

//******************************************************************************
// does it matter if a or b is in the core?
double ExternalField::dV_ab(const DiracSpinor &phia, const DiracSpinor &phib) {
  return phia * phib; // XX just stand-in
}

//******************************************************************************
DiracSpinor ExternalField::dV_ab_rhs(const DiracSpinor &phia,
                                     const DiracSpinor &phib) {
  return (phia * phib) * phib; // XX just stand-in
}
