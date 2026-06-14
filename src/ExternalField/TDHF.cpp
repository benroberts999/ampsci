#include "TDHF.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/format.hpp"
#include "qip/Widgets.hpp"
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>

namespace ExternalField {

//==============================================================================
TDHF::TDHF(const DiracOperator::TensorOperator *const h_plus,
           const HF::HartreeFock *const hf,
           const DiracOperator::TensorOperator *const h_minus)
  : CorePolarisation((assert(h_plus != nullptr), h_plus)),
    p_hf((assert(hf != nullptr), hf)),
    m_core(hf->core()),
    m_alpha(hf->alpha()),
    p_VBr(hf->vBreit()),
    m_h_minus(h_minus ? h_minus : h_plus) {
  initialise_dPsi();
}

//==============================================================================
void TDHF::initialise_dPsi() {
  // Initialise dPsi vectors, accounting for selection rules
  constexpr bool print = false;
  m_X.resize(m_core.size());
  for (auto ic = 0u; ic < m_core.size(); ic++) {
    const auto &Fc = m_core[ic];
    const auto pi_ch = Fc.parity() * m_pi;
    const auto tj_c = Fc.twoj();
    const auto tjmin_tmp = tj_c - 2 * m_rank;
    const auto tjmin = tjmin_tmp < 1 ? 1 : tjmin_tmp;
    const auto tjmax = tj_c + 2 * m_rank;
    if constexpr (print) {
      std::cout << "|" << Fc.symbol() << ">  -->  ";
    }
    for (int tj = tjmin; tj <= tjmax; tj += 2) {
      const auto l_minus = (tj - 1) / 2;
      const auto pi_chla = Angular::parity_l(l_minus) * pi_ch;

      if (Angular::triangle(tj_c, tj, 2 * m_rank) == 0)
        continue;

      const auto l = (pi_chla == 1) ? l_minus : l_minus + 1;
      const auto kappa = Angular::kappa_twojl(tj, l);
      m_X[ic].emplace_back(0, kappa, Fc.grid_sptr());
      m_X[ic].back().max_pt() = Fc.max_pt();
      if constexpr (print) {
        std::cout << "|" << m_X[ic].back().symbol() << "> + ";
      }
    }
    if constexpr (print) {
      std::cout << "\n";
    }
  }
  m_Y = m_X;
}

//==============================================================================
void TDHF::clear() {
  using namespace qip::overloads;
  m_X *= 0.0;
  m_Y *= 0.0;
}

//==============================================================================
const std::vector<DiracSpinor> &TDHF::get_dPsis(const DiracSpinor &Fc,
                                                dPsiType XorY) const {

  const auto index = static_cast<std::size_t>(
    std::find(m_core.cbegin(), m_core.cend(), Fc) - m_core.cbegin());
  assert(index < m_X.size());
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//==============================================================================
const DiracSpinor &TDHF::get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                    const int kappa_x) const {
  const auto &dPsis = get_dPsis(Fc, XorY);
  auto match_kappa_x = [=](const auto &Fa) { return Fa.kappa() == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//==============================================================================
std::vector<DiracSpinor>
TDHF::solve_dPsis(const DiracSpinor &Fv, const double omega, dPsiType XorY,
                  const MBPT::CorrelationPotential *const Sigma, StateType st,
                  bool incl_dV) const {
  std::vector<DiracSpinor> dFvs;
  const auto tjmin = std::max(1, Fv.twoj() - 2 * m_rank);
  const auto tjmax = Fv.twoj() + 2 * m_rank;
  for (int tjbeta = tjmin; tjbeta <= tjmax; tjbeta += 2) {
    const auto kappa = Angular::kappa_twojpi(tjbeta, Fv.parity() * m_pi);
    dFvs.push_back(solve_dPsi(Fv, omega, XorY, kappa, Sigma, st, incl_dV));
  }
  return dFvs;
}
//==============================================================================
DiracSpinor TDHF::solve_dPsi(const DiracSpinor &Fv, const double omega,
                             dPsiType XorY, const int kappa_x,
                             const MBPT::CorrelationPotential *const Sigma,
                             StateType st, bool incl_dV) const {
  // Solves (H + Sigma - e - w)X = -(h + dV - de)Psi
  // or     (H + Sigma - e + w)Y = -(h^dag + dV^dag - de)Psi

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;

  const auto imag = m_h->imaginaryQ();
  const auto *h_use = conj ? m_h_minus : m_h;

  auto rhs = h_use->reduced_rhs(kappa_x, Fv);
  if (imag && conj)
    rhs *= -1;

  if (incl_dV)
    rhs += dV_rhs(kappa_x, Fv, conj);
  if (kappa_x == Fv.kappa() && !imag) {
    auto de = h_use->reducedME(Fv, Fv);
    if (incl_dV)
      de += dV(Fv, Fv, conj);
    rhs -= de * Fv;
  }

  // Do we want |Y> or <Y| ?
  auto s2 = 1;
  if (st == StateType::bra) {
    // "left-hand-side" : "reduced" ket, so has factor (+ confugate)
    const auto sj =
      Angular::evenQ_2(Fv.twoj() - Angular::twoj_k(kappa_x)) ? 1 : -1;
    // if conj, extra * => +1
    const auto si = imag && !conj ? -1 : 1;
    s2 = sj * si;
  }

  const auto vl = p_hf->vlocal(Angular::l_k(Fv.kappa()));
  const auto &Hmag = p_hf->Hmag(Angular::l_k(Fv.kappa()));
  // The l from X ? or from Fv ?
  return s2 * ExternalField::solveMixedState(Fv, ww, vl, m_alpha, m_core, rhs,
                                             1.0e-12, Sigma, p_VBr, Hmag);
}

//==============================================================================
void TDHF::solve_ms_core(std::vector<DiracSpinor> &dFb, const DiracSpinor &Fb,
                         const std::vector<DiracSpinor> &hFbs,
                         const double omega, dPsiType XorY,
                         double eps_ms) const {
  // Solves (H - e - w)Xb = -(h + dV - de)Psi
  // or     (H - e + w)Y = -(h^dag + dV^dag - de)Psi
  // The diagonal (de) and near-resonant fine-structure partner terms are
  // conditioned inside solveMixedState (see conditioning_states()).

  // Continuum (photoionisation) mode: delegate to the V^{N-1} version
  // (active only via solve_core_cntm()).
  if (m_cntm_mode) {
    solve_ms_core_cntm(dFb, Fb, hFbs, omega, XorY, eps_ms);
    return;
  }

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;

  const auto imag = m_h->imaginaryQ();
  for (auto ibeta = 0ul; ibeta < dFb.size(); ibeta++) {

    auto &dF_beta = dFb[ibeta];
    const int kappa_beta = dF_beta.kappa();

    const auto &hFb = hFbs[ibeta];
    const auto s = (imag && conj) ? -1.0 : 1.0;
    auto rhs = s * hFb + dV_rhs(kappa_beta, Fb, conj);

    const auto vl = p_hf->vlocal(Angular::l_k(Fb.kappa()));
    const auto &Hmag = p_hf->Hmag(Angular::l_k(Fb.kappa()));
    // The l from X ? or from Fv ?
    ExternalField::solveMixedState(dF_beta, Fb, ww, vl, m_alpha, m_core, rhs,
                                   eps_ms, nullptr, p_VBr, Hmag);
  }
}

//==============================================================================
void TDHF::solve_ms_core_b(DiracSpinor &dF_beta, const DiracSpinor &Fb,
                           const DiracSpinor &hFb, const double omega,
                           dPsiType XorY, double eps_ms) const {
  // As solve_ms_core(), but for a single channel (one core orbital Fb, one
  // kappa projection beta, X or Y). Used to parallelise tdhf_core_it() over
  // (orbital x channel x X/Y) for better load balance than per-orbital.
  // Thread-safe: reads only const/shared state (and the previous iteration's
  // m_X/m_Y via dV_rhs); writes only dF_beta.

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;

  const auto imag = m_h->imaginaryQ();
  const int kappa_beta = dF_beta.kappa();
  const auto s = (imag && conj) ? -1.0 : 1.0;
  auto rhs = s * hFb + dV_rhs(kappa_beta, Fb, conj);

  const auto vl = p_hf->vlocal(Angular::l_k(Fb.kappa()));
  const auto &Hmag = p_hf->Hmag(Angular::l_k(Fb.kappa()));
  ExternalField::solveMixedState(dF_beta, Fb, ww, vl, m_alpha, m_core, rhs,
                                 eps_ms, nullptr, p_VBr, Hmag);
}

//==============================================================================
void TDHF::solve_ms_core_cntm(std::vector<DiracSpinor> &dFb,
                              const DiracSpinor &Fb,
                              const std::vector<DiracSpinor> &hFbs,
                              const double omega, dPsiType XorY,
                              double eps_ms) const {
  // Continuum (photoionisation) version of solve_ms_core. Active only via
  // solve_core_cntm(); see there for the V^{N-1} rearrangement.
  // Ionised orbital: en_+ = en_a + omega > 0. The one-electron (spherically
  // averaged) self-interaction V^a_0 = y^0_aa - X_a of the hole (= Fb itself)
  // is subtracted from the static Hamiltonian AND added (lagged) to the
  // source -- an exact rearrangement of the TDHF equations:
  //   (h_HF - V^a_0 - en_pm) phi = -(t + dV - de)phi_a - V^a_0 phi_bar,
  // so the fixed point is the unchanged TDHF one, but:
  //   - the long-ranged hole monopole -y^0_aa is inverted directly (local
  //     potential vl - y^0_aa: residual-ion tail, Z_ion = 1, consistent with
  //     the continuum Coulomb boundary conditions);
  //   - the lagged source is short-ranged: the +y^0_aa*phi_bar term cancels
  //     pointwise the 1/r tail hidden in the b=a part of dV*phi_a;
  //   - the exchange part is carried by the mixed-state solvers (Fhole = &Fb:
  //     vexFa -> vexFa - vexFa_1el on the current iterate).
  // Both partners (X/+ and Y/-) of the ionised orbital are treated this way.
  using namespace qip::overloads;

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;
  const auto imag = m_h->imaginaryQ();

  const bool ionisedQ = !m_suppress_open && (Fb.en() + std::abs(ww) > 0.0);
  // Only the X/+ channel (ww = +omega) reaches en_+ > 0 (continuum solve).
  const bool continuumQ = ionisedQ && (Fb.en() + ww > 0.0);
  // suppress_open: hold the open channels at zero instead.
  const bool suppressQ = m_suppress_open && (Fb.en() + ww > 0.0);

  // One-electron direct potential of the hole orbital (the hole is Fb itself).
  const auto y0aa =
    ionisedQ ? Coulomb::yk_ab(0, Fb, Fb) : std::vector<double>{};

  for (auto ibeta = 0ul; ibeta < dFb.size(); ibeta++) {
    auto &dF_beta = dFb[ibeta];
    const int kappa_beta = dF_beta.kappa();

    if (suppressQ) {
      // Open channel excluded from the response (explicit zero: clears nan)
      dF_beta.f().assign(dF_beta.f().size(), 0.0);
      dF_beta.g().assign(dF_beta.g().size(), 0.0);
      continue;
    }

    const auto &hFb = hFbs[ibeta];
    const auto s = (imag && conj) ? -1.0 : 1.0;
    auto rhs = s * hFb + dV_rhs(kappa_beta, Fb, conj);

    // + V^a_0 phi_bar (lagged): dF_beta still holds the previous iterate
    if (ionisedQ) {
      rhs += (y0aa * dF_beta) + HF::vexFa_1el(dF_beta, Fb);
    }

    if (continuumQ) {
      const auto vl_c = p_hf->vlocal(Angular::l_k(kappa_beta)) - y0aa;
      DiracSpinor Freg{0, kappa_beta, Fb.grid_sptr()};
      DiracSpinor Firr{0, kappa_beta, Fb.grid_sptr()};
      ExternalField::solveContinuumMixedState(
        dF_beta, Freg, Firr, Fb, ww, vl_c, m_alpha, m_core, rhs, eps_ms, &Fb);
      continue;
    }

    if (ionisedQ) {
      // Bound (Y/-) partner of an ionised orbital: same V^{N-1} treatment.
      const auto vl_c = p_hf->vlocal(Angular::l_k(Fb.kappa())) - y0aa;
      const auto &Hmag = p_hf->Hmag(Angular::l_k(Fb.kappa()));
      ExternalField::solveMixedState_cntm(
        dF_beta, Fb, ww, vl_c, m_alpha, m_core, rhs, eps_ms, p_VBr, Hmag, &Fb);
      continue;
    }

    // Closed orbital (not ionised at this omega): bound solve (robust against
    // the spurious-preconditioner divergence at high omega -- see
    // solveMixedState_cntm()).
    const auto vl = p_hf->vlocal(Angular::l_k(Fb.kappa()));
    const auto &Hmag = p_hf->Hmag(Angular::l_k(Fb.kappa()));
    ExternalField::solveMixedState_cntm(dF_beta, Fb, ww, vl, m_alpha, m_core,
                                        rhs, eps_ms, p_VBr, Hmag);
  }
}

//==============================================================================
std::vector<std::vector<DiracSpinor>>
TDHF::form_hFcore(const DiracOperator::TensorOperator *h) const {
  std::vector<std::vector<DiracSpinor>> hFcore(m_core.size());
  for (auto ic = 0ul; ic < m_X.size(); ic++) {
    const auto &Fc = m_core[ic];
    hFcore.reserve(m_X[ic].size()); // each h projection
    for (auto beta = 0ul; beta < m_X[ic].size(); beta++) {
      const auto &Xx = m_X[ic][beta];
      hFcore[ic].push_back(h->reduced_rhs(Xx.kappa(), Fc));
    }
  }
  return hFcore;
}

//==============================================================================
std::pair<double, std::string>
TDHF::eps_dPsi(const std::vector<std::vector<DiracSpinor>> &Xnew,
               bool relative) const {

  // Relative L2 change of the (undamped) X spinors, summed over all core
  // orbitals and channels: ratio = Sum|dX|^2 / Sum|X_new|^2. Returns sqrt(ratio)
  // -- the relative change -- if `relative`, else the squared ratio.
  // (Y need not be included: X and Y are solved with the same dV, so once X is
  // converged dV is converged and Y too.)

  using namespace qip::overloads;
  double DdF2 = 0.0;
  double dF2 = 0.0;
  double worst = 0.0;
  std::string worst_lab;
  for (std::size_t ib = 0; ib < m_core.size(); ++ib) {
    for (std::size_t i = 0; i < Xnew[ib].size(); ++i) {
      const auto &nw = Xnew[ib][i];
      const auto d = (nw - m_X[ib][i]).norm2();
      const auto n = nw.norm2();
      DdF2 += d;
      dF2 += n;
      // per-channel ratio: used only to rank the worst channel for reporting
      const auto e = n == 0.0 ? 0.0 : d / n;
      if (e > worst) {
        worst = e;
        worst_lab = m_core[ib].shortSymbol() + "," + nw.shortSymbol();
      }
    }
  }
  const auto ratio = dF2 == 0.0 ? 0.0 : DdF2 / dF2;
  return {relative ? std::sqrt(ratio) : ratio, worst_lab};
}

//==============================================================================
std::pair<double, std::string> TDHF::tdhf_core_it(double omega,
                                                  double eta_damp) {
  // solve TDHF equation for core - single iteration
  using namespace qip::overloads;
  const bool staticQ = std::abs(omega) < 1.0e-10;
  const bool static_Y = staticQ && !m_h->freqDependantQ();
  const auto s = m_imag ? -1 : 1;

  // make copies (instead of blank spinors) since the solve_mixed_states
  // is faster if it starts from existing solution
  auto Xs = m_X;
  auto Ys = m_Y;

  const auto eps_ms = 1.0e-12;

  // Flatten all (core orbital x channel x X/Y) mixed-states solves into one
  // task list, so the parallel loop has many small, fairly-uniform tasks
  // instead of a few very uneven per-orbital ones (much better load balance).
  // Each task writes a distinct dF and only reads shared state, so they are
  // independent. See solve_ms_core_b().
  struct MsTask {
    DiracSpinor *dF;        // target (in Xs or Ys)
    const DiracSpinor *Fb;  // core orbital
    const DiracSpinor *hFb; // source projection h|Fb>
    dPsiType type;
  };
  std::vector<MsTask> tasks;
  for (auto ib = 0ul; ib < m_core.size(); ib++) {
    for (auto be = 0ul; be < Xs[ib].size(); be++) {
      tasks.push_back(
        {&Xs[ib][be], &m_core[ib], &m_hFcore[ib][be], dPsiType::X});
      if (!static_Y) {
        tasks.push_back(
          {&Ys[ib][be], &m_core[ib], &m_hFcore_minus[ib][be], dPsiType::Y});
      }
    }
  }

  // Solve for the (undamped) dF of this iteration.
#pragma omp parallel for schedule(dynamic)
  for (auto it = 0ul; it < tasks.size(); it++) {
    const auto &t = tasks[it];
    solve_ms_core_b(*t.dF, *t.Fb, *t.hFb, omega, t.type, eps_ms);
  }

  // Measure convergence on the undamped X solution, before damping.
  const auto eps = eps_dPsi(Xs, m_eps_sqrt);

  // Damp the solution (for the next iteration) and store.
#pragma omp parallel for
  for (auto ib = 0ul; ib < m_core.size(); ib++) {
    Xs[ib] = eta_damp * m_X[ib] + (1.0 - eta_damp) * Xs[ib];
    if (static_Y) {
      Ys[ib] = s * Xs[ib];
    } else {
      Ys[ib] = eta_damp * m_Y[ib] + (1.0 - eta_damp) * Ys[ib];
    }
  }
  m_X = std::move(Xs);
  m_Y = std::move(Ys);

  return eps;
}

//==============================================================================
void TDHF::solve_core(double omega, int max_its, bool print) {

  assert(m_h->rank() == m_rank && "Rank must match in solve_core");
  assert(m_h->parity() == m_pi && "Parity must match in solve_core");
  assert(m_h->imaginaryQ() == m_imag && "Imaginarity must match in solve_core");

  const double converge_targ = m_eps;
  const auto eta_damp = m_eta;
  omega = std::abs(omega);

  m_hFcore = form_hFcore(m_h);
  m_hFcore_minus = form_hFcore(m_h_minus);

  std::pair<double, std::string> eps{};
  double best_eps{1.0e30};
  int count_worse = 0;
  int it{0};
  qip::LiveMessage status(
    fmt::format("TDHF {} (w={:.4f}): ", m_h->name(), omega), print);
  for (; it < max_its; it++) {
    const auto eta = it == 0 ? 0.0 : eta_damp;
    eps = tdhf_core_it(omega, eta);

    status(fmt::format("{:2d} {:.1e} [{}]", it, eps.first, eps.second));

    if (eps.first < converge_targ || std::isnan(eps.first))
      break; // converged (or broken)

    // Stalled: no meaningful improvement for several iterations -> give up.
    // (The result is still used; the warning stars flag its quality.)
    if (eps.first < 0.99 * best_eps) {
      best_eps = eps.first;
      count_worse = 0;
    } else if (++count_worse > 5) {
      break;
    }
  }

  // Soft visual warning, relative to the convergence target.
  const auto stars = (max_its <= 1)                      ? "" :
                     (eps.first > 1.0e6 * converge_targ) ? "  ***" :
                     (eps.first > 1.0e4 * converge_targ) ? "  **" :
                     (eps.first > 1.0e2 * converge_targ) ? "  *" :
                                                           "";
  status.done(stars);

  // set last eps (convergance) and frequency (omega)
  m_core_eps = eps.first;
  m_core_its = it;
  m_core_omega = omega;
}

//==============================================================================
void TDHF::solve_core_cntm(double omega, int max_its, bool print,
                           bool suppress_open) {
  // Continuum (photoionisation) RPA. Enables the V^{N-1} continuum branches in
  // solve_ms_core (via m_cntm_mode), then reuses the standard damped solve_core
  // driver -- the same proven iteration as the bound RPA. The continuum
  // mixed-state solve uses forward integration (regular at the origin), so the
  // fixed-point map stays bounded and damped iteration converges. (cf. Dzuba's
  // continuum solver, which likewise matches the inner/outer integrations at a
  // turning point and uses adaptive damping, rather than carrying the
  // origin-divergent irregular solution through a global Green's function.)
  m_cntm_mode = true;
  m_suppress_open = suppress_open;
  solve_core(omega, max_its, print);
  m_cntm_mode = false;
  m_suppress_open = false;
}

//==============================================================================
// does it matter if a or b is in the core?
double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm, bool conj) const {
  const auto s = conj && m_h->imaginaryQ() ? -1 : 1; // careful. OK?
  return s * Fn * dV_rhs(Fn.kappa(), Fm, conj);
}

//==============================================================================
double TDHF::dV(const DiracSpinor &Fn, const DiracSpinor &Fm) const {
  const auto conj = Fm.en() > Fn.en();
  return dV(Fn, Fm, conj);
}

//==============================================================================
double TDHF::dV_cont(const DiracSpinor &Fe, const DiracSpinor &Fa) const {
  // As dV(Fe,Fa), but for a continuum Fe: add the one-electron self-
  // interaction term +V^a_0 phi_pm to the source, matching the V^{N-1}
  // rearrangement used in solve_ms_core_cntm (see solve_core_cntm). Its
  // +y^0_aa*phi term cancels the 1/r tail hidden in the b=a part of dV*phi_a
  // pointwise, so the continuum-continuum overlap is box-independent; Fe must
  // be the V^{N-1} continuum state (hole = Fa, incl. the exchange part).
  using namespace qip::overloads;
  const auto conj = Fa.en() > Fe.en();
  const auto s = conj && m_h->imaginaryQ() ? -1 : 1;
  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto &chi = get_dPsi_x(Fa, ChiType, Fe.kappa());
  const auto y0aa = Coulomb::yk_ab(0, Fa, Fa);
  const auto V0chi = (y0aa * chi) + HF::vexFa_1el(chi, Fa);
  return s * (Fe * (dV_rhs(Fe.kappa(), Fa, conj) + V0chi));
}

//==============================================================================
DiracSpinor TDHF::dV_rhs(int kappa_n, const DiracSpinor &Fa, bool conj) const {

  auto dVFa = DiracSpinor(0, kappa_n, Fa.grid_sptr());
  dVFa.max_pt() = Fa.max_pt();

  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto EtaType = !conj ? dPsiType::Y : dPsiType::X;

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);
  const auto tjn = Angular::twoj_k(kappa_n);

  // nb: faster to not //ize this one
  for (const auto &Fb : m_core) {

    const auto &X_b = get_dPsis(Fb, ChiType);
    const auto &Y_b = get_dPsis(Fb, EtaType);

    for (auto ibeta = 0ul; ibeta < X_b.size(); ++ibeta) {
      const auto &X_beta = X_b[ibeta];
      const auto &Y_beta = Y_b[ibeta];

      const auto sQ = Angular::neg1pow_2(tjn - X_beta.twoj() + 2 * k);
      if (sQ == 1) { // faster than multiplying!
        dVFa += Coulomb::Wkv_bcd(k, kappa_n, Fb, Fa, X_beta);
        dVFa += Coulomb::Wkv_bcd(k, kappa_n, Y_beta, Fa, Fb);
      } else {
        dVFa -= Coulomb::Wkv_bcd(k, kappa_n, Fb, Fa, X_beta);
        dVFa -= Coulomb::Wkv_bcd(k, kappa_n, Y_beta, Fa, Fb);
      }

      // Breit part:
      if (p_VBr) {
        // Makes small difference for E1, huge difference for HFS
        // No frequency-dependence yet. Not sure if we should bother.
        dVFa += tkp1 * p_VBr->dV_Br(kappa_n, k, Fa, Fb, X_beta, Y_beta);
      }
    }
  }

  dVFa *= (1.0 / tkp1);

  return dVFa;
}

} // namespace ExternalField
