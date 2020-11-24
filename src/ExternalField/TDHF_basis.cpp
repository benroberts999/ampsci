#include "TDHF_basis.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <vector>

namespace ExternalField {

//******************************************************************************
TDHFbasis::TDHFbasis(const DiracOperator::TensorOperator *const h,
                     const HF::HartreeFock *const hf,
                     const std::vector<DiracSpinor> &basis)
    : TDHF(h, hf), m_basis(basis) {}

//******************************************************************************
DiracSpinor TDHFbasis::form_dPsi(const DiracSpinor &Fv, const double omega,
                                 dPsiType XorY, const int kappa_beta,
                                 const std::vector<DiracSpinor> &spectrum,
                                 StateType st) const {

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj; //?

  const auto imag = m_h->imaginaryQ();
  const auto s = (imag && conj) ? -1 : 1;

  auto s2 = 1;
  if (st == StateType::bra) {
    // "left-hand-side" : "reduced" ket, so has factor (+ confugate)
    const auto sj =
        Angular::evenQ_2(Fv.twoj() - Angular::twoj_k(kappa_beta)) ? 1 : -1;
    // if conj, extra * => +1
    const auto si = imag && !conj ? -1 : 1;
    s2 = sj * si;
  }

  // nb: Spectrum must constain Sigma to include correlations
  DiracSpinor Xx{0, kappa_beta, Fv.rgrid};
  Xx.p0 = Fv.pinf;
  Xx.pinf = Fv.pinf;
  for (const auto &Fn : spectrum) {
    if (Fn.k != Xx.k || m_h->isZero(Fn.k, Fv.k) || Fv == Fn)
      continue;
    // XXX Check:
    const auto hnc = s2 * (s * m_h->reducedME(Fn, Fv) + dV(Fn, Fv, conj));
    Xx += (hnc / (Fv.en - Fn.en + ww)) * Fn;
  }

  return Xx;
}

//******************************************************************************
//! Forms \delta Psi_v for valence state Fv for all kappas (see solve_dPsi)
std::vector<DiracSpinor>
TDHFbasis::form_dPsis(const DiracSpinor &Fv, const double omega, dPsiType XorY,
                      const std::vector<DiracSpinor> &spectrum,
                      StateType st) const {
  std::vector<DiracSpinor> dFvs;
  const auto tjmin = std::max(1, Fv.twoj() - 2 * m_rank);
  const auto tjmax = Fv.twoj() + 2 * m_rank;
  for (int tjbeta = tjmin; tjbeta <= tjmax; tjbeta += 2) {
    const auto kappa = Angular::kappa_twojpi(tjbeta, Fv.parity() * m_pi);
    dFvs.push_back(form_dPsi(Fv, omega, XorY, kappa, spectrum, st));
  }
  return dFvs;
}

//******************************************************************************
void TDHFbasis::solve_core(const double omega, int max_its, const bool print) {
  const double converge_targ = 1.0e-8;
  const auto damper = HF::rampedDamp(0.75, 0.25, 1, 20);

  const bool staticQ = std::abs(omega) < 1.0e-10;

  auto tmp_X = m_X;
  auto tmp_Y = m_Y;

  auto eps = 0.0;
  int it = 0;
  if (print) {
    printf("TDHFb (w=%.3f): ", omega);
    std::cout << std::flush;
  }
  for (; it < max_its; it++) {
    eps = 0.0;
    std::vector<double> eps_vec(m_core.size(), 0.0);
    const auto a_damp = (it == 0) ? 0.0 : damper(it);
#pragma omp parallel for
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      double eps_c = 0.0;
      const auto &Fc = m_core[ic];

      for (auto j = 0ul; j < m_X[ic].size(); j++) {
        auto &Xx = tmp_X[ic][j];
        const auto &oldX = m_X[ic][j];
        Xx *= a_damp;
        Xx += (1.0 - a_damp) *
              form_dPsi(Fc, omega, dPsiType::X, Xx.k, m_basis, StateType::ket);

        const auto delta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
        if (delta > eps_c)
          eps_c = delta;

        // Xx = a_damp * oldX + (1.0 - a_damp) * Xx;
        if (!staticQ) {
          auto &Yx = tmp_Y[ic][j];
          const auto &oldY = m_Y[ic][j];
          Yx *= a_damp;
          Yx += (1.0 - a_damp) * form_dPsi(Fc, omega, dPsiType::Y, Xx.k,
                                           m_basis, StateType::ket);
        } else {
          tmp_Y[ic][j] = Xx;
        }
      }
      eps_vec[ic] = eps_c;
    }

    m_X = tmp_X;
    m_Y = tmp_Y;
    eps = *std::max_element(cbegin(eps_vec), cend(eps_vec));

    if ((it > 1 && eps < converge_targ))
      break;
  }
  if (print) {
    printf("%2i %.1e\n", it, eps);
  }
  std::cout << std::flush;
  m_core_eps = eps;
  m_core_omega = omega;
}

} // namespace ExternalField
