#include "TDHFbasis.hpp"
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
    // why need multiply dV by s too? Thought I shouldn't??
    const auto hnc = s2 * (s * m_h->reducedME(Fn, Fv) + s * dV(Fn, Fv, conj));
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
    printf("TDHFb %s (w=%.3f): ", m_h->name().c_str(), omega);
    std::cout << std::flush;
  }

  // Store 2D indexes, for more efficient //isation
  std::vector<std::pair<std::size_t, std::size_t>> indexs;
  for (auto ic = 0ul; ic < m_core.size(); ic++) {
    for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
      indexs.emplace_back(ic, beta);
    }
  }

  for (; it < max_its; it++) {
    // eps = 0.0;
    // std::vector<double> eps_vec(m_core.size(), 0.0);
    const auto a_damp = (it == 0) ? 0.0 : damper(it);

#pragma omp parallel for
    for (auto i = 0ul; i < indexs.size(); i++) {
      const auto [ic, beta] = indexs[i];
      const auto &Fc = m_core[ic];

      auto &Xx = tmp_X[ic][beta];
      Xx *= a_damp;
      Xx += (1.0 - a_damp) *
            form_dPsi(Fc, omega, dPsiType::X, Xx.k, m_basis, StateType::ket);

      if (!staticQ) {
        auto &Yx = tmp_Y[ic][beta];
        Yx *= a_damp;
        Yx += (1.0 - a_damp) *
              form_dPsi(Fc, omega, dPsiType::Y, Xx.k, m_basis, StateType::ket);
      } else {
        const auto s = m_h->imaginaryQ() ? -1 : 1;
        tmp_Y[ic][beta] = s * Xx;
      }
    }

    // Find convergance from worst X
    eps = 0.0;
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
        const auto &Xx = tmp_X[ic][beta];
        const auto &oldX = m_X[ic][beta];
        const auto eps_beta = (Xx - oldX) * (Xx - oldX) / (Xx * Xx);
        if (eps_beta > eps) {
          eps = eps_beta;
        }
      }
    }

    m_X = tmp_X;
    m_Y = tmp_Y;

    if ((it > 1 && eps < converge_targ))
      break;
  }
  if (print) {
    printf("%2i %.1e\n", it, eps);
  }
  std::cout << std::flush;
  m_core_eps = eps;
  m_core_omega = omega;
} // namespace ExternalField

} // namespace ExternalField
