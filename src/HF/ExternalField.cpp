#include "ExternalField.hpp"
#include "Angular/Wigner_369j.hpp"
#include "Dirac/DiracOperator.hpp"
#include "Dirac/DiracSpinor.hpp"
#include "HF/HartreeFockClass.hpp"
#include <algorithm>
#include <vector>

//******************************************************************************
ExternalField::ExternalField(const DiracOperator *const h,
                             const std::vector<DiracSpinor> &core,
                             const std::vector<double> &vl, const double alpha,
                             const double omega)
    : m_h(h), p_core(&core), m_vl(vl), m_alpha(alpha), m_omega(omega), //
      static_fieldQ(std::abs(omega) > 0 ? false : true),               //
      m_rank(h->rank()), m_pi(h->parity()), m_imag(h->imaginaryQ())    //
// I think omega must be >0 ?? enforce?
{
  m_X.resize(core.size());
  for (auto ic = 0u; ic < core.size(); ic++) {
    const auto &phic = core[ic];
    const auto pi_ch = phic.parity() * m_pi;
    const auto tj_c = phic.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Wigner::parity_l(l_minus) * pi_ch;
      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      // XXX Add option to restrict l = lc
      const auto kappa = Wigner::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, *(phic.p_rgrid));
    }
  }
  if (!static_fieldQ)
    m_Y = m_X;

  for (const auto &phic : core) {
    auto &dPsi = get_dPsis(phic, dPsiType::X);
    std::cout << phic.symbol() << " : ";
    for (const auto &dX : dPsi) {
      std::cout << dX.symbol() << " ";
    }
    std::cout << "\n";
  }
}
//******************************************************************************
std::size_t ExternalField::core_index(const DiracSpinor &phic) {
  // XXX Note: no bounds checking here!
  return static_cast<std::size_t>(
      std::find(p_core->cbegin(), p_core->cend(), phic) - p_core->cbegin());
}

//******************************************************************************
const std::vector<DiracSpinor> &
ExternalField::get_dPsis(const DiracSpinor &phic, dPsiType XorY) {
  auto index = core_index(phic);
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
void ExternalField::solve_TDHFcore() {

#pragma omp parallel for
  for (auto ic = 0u; ic < p_core->size(); ic++) {
    const auto &phic = (*p_core)[ic];
    auto &dPsis_X = m_X[ic];
    for (auto &Xx : dPsis_X) {
      auto x = Xx.k;
      const auto hPsic = m_h->reduced_rhs(x, phic); // always same!?
      Xx = HartreeFock::solveMixedState(phic, x, m_omega, m_vl, m_alpha,
                                        *p_core, hPsic);
    }
    if (!static_fieldQ) {
      auto &dPsis_Y = m_Y[ic];
      for (auto &Yx : dPsis_Y) {
        auto x = Yx.k;
        const auto hPsic = m_h->reduced_rhs(x, phic); // always same!?
        Yx = HartreeFock::solveMixedState(phic, x, -m_omega, m_vl, m_alpha,
                                          *p_core, hPsic);
      }
    }
  }

  // for (const auto &phic : *p_core) {
  //   auto &dPsi = get_dPsis(phic, dPsiType::X);
  //   for (const auto &Xx : dPsi) {
  //     std::cout << phic.symbol() << " " << Xx.symbol() << " " << Xx * Xx << "
  //     "
  //               << Xx * phic << "\n";
  //   }
  // }
  //
}

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
