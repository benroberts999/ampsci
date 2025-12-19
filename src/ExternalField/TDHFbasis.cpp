#include "TDHFbasis.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <vector>

namespace ExternalField {

//==============================================================================
TDHFbasis::TDHFbasis(const DiracOperator::TensorOperator *const h,
                     const HF::HartreeFock *const hf,
                     const std::vector<DiracSpinor> &basis)
    : TDHF(h, hf) {

  // m_basis = basis;
  // return;

  // find maximum n for each kappa in core:
  const auto max_ki_core = DiracSpinor::max_kindex(hf->core());
  std::vector<int> max_n_ka_core;
  for (int ki = 0; ki <= max_ki_core; ++ki) {
    max_n_ka_core.push_back(
        DiracSpinor::max_n(hf->core(), Angular::kappaFromIndex(ki)));
  }

  // only include excited states in spectrum?
  for (const auto &n : basis) {
    const auto ki = std::size_t(n.k_index());
    const auto n_core = ki < max_n_ka_core.size() ? max_n_ka_core.at(ki) : 0;
    if (n.n() > n_core) {
      m_basis.push_back(n);
    }
  }
}

//==============================================================================
DiracSpinor TDHFbasis::form_dPsi(const DiracSpinor &Fv, const double omega,
                                 dPsiType XorY, const int kappa_beta,
                                 const std::vector<DiracSpinor> &spectrum,
                                 StateType st, bool incl_dV) const {

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj; //?

  const auto imag = m_h->imaginaryQ();
  const auto s = (imag && conj) ? -1 : 1;

  auto s2 = 1;
  if (st == StateType::bra) {
    // "left-hand-side" : "reduced" ket, so has factor (+ conjugate)
    const auto sj =
        Angular::evenQ_2(Fv.twoj() - Angular::twoj_k(kappa_beta)) ? 1 : -1;
    // if conj, extra * => +1
    const auto si = imag && !conj ? -1 : 1;
    s2 = sj * si;
  }

  // nb: Spectrum must constain Sigma to include correlations
  DiracSpinor Xx{0, kappa_beta, Fv.grid_sptr()};
  Xx.min_pt() = Fv.max_pt();
  Xx.max_pt() = Fv.max_pt();
  for (const auto &Fn : spectrum) {
    if (Fn.kappa() != Xx.kappa() || m_h->isZero(Fn.kappa(), Fv.kappa()) ||
        Fv == Fn)
      continue;
    // XXX Check:
    // why need multiply dV by s too? Thought I shouldn't??
    // const auto hnc = s2 * (s * m_h->reducedME(Fn, Fv) + s * dV(Fn, Fv,
    // conj));
    const auto hnc = incl_dV ?
                         s2 * s * (m_h->reducedME(Fn, Fv) + dV(Fn, Fv, conj)) :
                         s2 * s * m_h->reducedME(Fn, Fv);
    Xx += (hnc / (Fv.en() - Fn.en() + ww)) * Fn;
  }

  return Xx;
}

//==============================================================================
//! Forms \delta Psi_v for valence state Fv for all kappas (see solve_dPsi)
std::vector<DiracSpinor>
TDHFbasis::form_dPsis(const DiracSpinor &Fv, const double omega, dPsiType XorY,
                      const std::vector<DiracSpinor> &spectrum, StateType st,
                      bool incl_dV) const {
  std::vector<DiracSpinor> dFvs;
  const auto tjmin = std::max(1, Fv.twoj() - 2 * m_rank);
  const auto tjmax = Fv.twoj() + 2 * m_rank;

  for (int tjbeta = tjmin; tjbeta <= tjmax; tjbeta += 2) {

    const auto pi_vh = Fv.parity() * m_pi;
    const auto tj_v = Fv.twoj();
    const auto l_minus = (tjbeta - 1) / 2;
    const auto pi_chla = Angular::parity_l(l_minus) * pi_vh;
    if (Angular::triangle(tj_v, tjbeta, 2 * m_rank) == 0)
      continue;

    const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
    const auto kappa = Angular::kappa_twojl(tjbeta, l);

    dFvs.push_back(form_dPsi(Fv, omega, XorY, kappa, spectrum, st, incl_dV));
  }
  return dFvs;
}

//==============================================================================
void TDHFbasis::solve_core(const double omega, int max_its, const bool print) {

  const auto converge_targ = m_eps;
  const auto eta_damp = 1.0e-5 * m_eta;
  const auto ww = std::abs(omega);

  const bool staticQ = std::abs(omega) < 1.0e-10;

  auto tmp_X = m_X;
  auto tmp_Y = m_Y;

  auto eps = 0.0;
  std::string s_worst;

  if (print) {
    printf("TDHFb %s (w=%.3f): ", m_h->name().c_str(), ww);
    std::cout << std::flush;
  }

  // Store 2D indexes, for more efficient //isation
  std::vector<std::pair<std::size_t, std::size_t>> indexs;
  for (auto ic = 0ul; ic < m_core.size(); ic++) {
    for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
      indexs.emplace_back(ic, beta);
    }
  }

  int it = 0;
  for (; it < max_its; it++) {
    const auto a_damp = (it == 0) ? 0.0 : eta_damp;

#pragma omp parallel for
    for (auto i = 0ul; i < indexs.size(); i++) {
      const auto [ic, beta] = indexs[i];
      const auto &Fc = m_core[ic];

      auto &Xx = tmp_X[ic][beta];
      Xx *= a_damp;
      Xx += (1.0 - a_damp) *
            form_dPsi(Fc, ww, dPsiType::X, Xx.kappa(), m_basis, StateType::ket);

      if (!staticQ) {
        auto &Yx = tmp_Y[ic][beta];
        Yx *= a_damp;
        Yx += (1.0 - a_damp) * form_dPsi(Fc, ww, dPsiType::Y, Xx.kappa(),
                                         m_basis, StateType::ket);
      } else {
        const auto s = m_h->imaginaryQ() ? -1 : 1;
        tmp_Y[ic][beta] = s * Xx;
      }
    }

    // Find convergance from worst X
    double eps_top = 0.0;
    double eps_bottom = 0.0;
    double worst = 0.0;
    for (auto ic = 0ul; ic < m_core.size(); ic++) {
      for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
        const auto &Xx = tmp_X[ic][beta];
        const auto &oldX = m_X[ic][beta];
        const auto t_eps_top = (Xx - oldX).norm2();
        const auto t_eps_bottom = 0.5 * (Xx + oldX).norm2();
        eps_top += t_eps_top;
        eps_bottom += t_eps_bottom;
        if (t_eps_top / t_eps_bottom > worst) {
          worst = t_eps_top / t_eps_bottom;
          s_worst = m_core[ic].shortSymbol() + "," + Xx.shortSymbol();
        }
      }
    }
    eps = eps_top / eps_bottom;

    m_X = tmp_X;
    m_Y = tmp_Y;

    if (it > 1 && eps < converge_targ)
      break;
  }
  if (print) {
    printf("%2i %.1e [%s]\n", it, eps, s_worst.c_str());
  }
  std::cout << std::flush;
  m_core_eps = eps;
  m_core_its = it;
  m_core_omega = ww; //?
}

} // namespace ExternalField
