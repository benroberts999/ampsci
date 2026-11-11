#include "TDHF_DCP.hpp"
#include "Angular/CkTable.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/Anderson.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/Breit.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/ChronoTimer.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Wavefunction/BSplineBasis.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "fmt/color.hpp"
#include "fmt/format.hpp"
#include "qip/Widgets.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <memory>
#include <vector>

namespace ExternalField {

//==============================================================================

TDHF_DCP::TDHF_DCP(const ExternalField::TDHF *const T,
                   const ExternalField::TDHF *const S, int K,
                   const HF::HartreeFock *const hf)
  // Rank K, parity pi_t*pi_s. Matrix elements of the pair are imaginary if
  // those of exactly one of t, s are
  : CorePolarisation(K, T->parity() * S->parity(), T->imagQ() != S->imagQ()),
    p_hf((assert(hf != nullptr), hf)),
    m_core(hf->core()),
    m_alpha(hf->alpha()),
    p_VBr(hf->vBreit()),
    m_T(T),
    m_S(S) {

  assert(Angular::triangle(T->rank(), S->rank(), K) != 0 &&
         "K must satisfy the triangle rule with the ranks of t and s");

  initialise_dPsi();
}

//==============================================================================

void TDHF_DCP::initialise_dPsi() {
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

void TDHF_DCP::clear() {
  using namespace qip::overloads;
  m_X *= 0.0;
  m_Y *= 0.0;
}

//==============================================================================

const std::vector<DiracSpinor> &TDHF_DCP::get_dPsis(const DiracSpinor &Fc,
                                                    dPsiType XorY) const {

  const auto index = static_cast<std::size_t>(
    std::find(m_core.cbegin(), m_core.cend(), Fc) - m_core.cbegin());
  assert(index < m_X.size());
  return XorY == dPsiType::X ? m_X[index] : m_Y[index];
}

//==============================================================================

const DiracSpinor &TDHF_DCP::get_dPsi_x(const DiracSpinor &Fc, dPsiType XorY,
                                        const int kappa_x) const {
  const auto &dPsis = get_dPsis(Fc, XorY);
  auto match_kappa_x = [=](const auto &Fa) { return Fa.kappa() == kappa_x; };
  return *std::find_if(dPsis.cbegin(), dPsis.cend(), match_kappa_x);
}

//==============================================================================

DiracSpinor TDHF_DCP::source_rhs(const DiracSpinor &Fa, dPsiType XorY,
                                 const int kappa_beta) const {
  // [sigma^a]^K_beta - [lambda^a]^K_beta:
  // sigma: each first-order operator (with its core polarisation) acting on
  // the first-order corrections to Fa of the other operator;
  // lambda: the first-order energy shifts of Fa times those corrections.
  // X: source of the ++ equation, from Xt, Xs (t, s both at positive
  // frequency). Y: that of the -- equation, from Yt, Ys.
  //
  // No explicit factors of iota (+1 real operator, -1 imaginary): the
  // first-order Y corrections are stored with their iota,
  // Y = iota * [phi_-], as is dV_rhs(conj). So the solutions of the --
  // equation are dY = iota_t * iota_s * [phi_--].

  const auto conj = XorY == dPsiType::Y;
  const auto tja = Fa.twoj();
  const auto tjbeta = Angular::twoj_k(kappa_beta);

  // First-order corrections to Fa: (Xt, Xs), or (Yt, Ys)
  const auto &phi_t = m_T->get_dPsis(Fa, XorY);
  const auto &phi_s = m_S->get_dPsis(Fa, XorY);

  DiracSpinor rhs{0, kappa_beta, Fa.grid_sptr()};

  // sigma. Coefficients are zero unless the operator connects alpha to beta
  for (const auto &phi_s_alpha : phi_s) {
    const auto [c_ts, c_st] =
      source_coefs(tja, phi_s_alpha.twoj(), tjbeta, conj);
    if (c_ts != 0.0) {
      rhs += c_ts * h_tilde_rhs(*m_T, kappa_beta, phi_s_alpha, conj);
    }
  }
  for (const auto &phi_t_alpha : phi_t) {
    const auto [c_ts, c_st] =
      source_coefs(tja, phi_t_alpha.twoj(), tjbeta, conj);
    if (c_st != 0.0) {
      rhs += c_st * h_tilde_rhs(*m_S, kappa_beta, phi_t_alpha, conj);
    }
  }

  // lambda. The shift from t multiplies the s correction: t takes a to a,
  // then s takes a to beta; so the intermediate state is a, and the
  // coefficient is that of s acting on a t correction (and t <-> s)
  const auto [c_ts, c_st] = source_coefs(tja, tja, tjbeta, conj);
  const auto de_t = c_st != 0.0 ? h_tilde_de(*m_T, Fa, conj) : 0.0;
  const auto de_s = c_ts != 0.0 ? h_tilde_de(*m_S, Fa, conj) : 0.0;
  if (de_t != 0.0) {
    for (const auto &phi_s_beta : phi_s) {
      if (phi_s_beta.kappa() == kappa_beta) {
        rhs -= (c_st * de_t) * phi_s_beta;
      }
    }
  }
  if (de_s != 0.0) {
    for (const auto &phi_t_beta : phi_t) {
      if (phi_t_beta.kappa() == kappa_beta) {
        rhs -= (c_ts * de_s) * phi_t_beta;
      }
    }
  }

  return rhs;
}

//==============================================================================

std::pair<double, double>
TDHF_DCP::source_coefs(int twoj_a, int twoj_n, int twoj_beta, bool conj) const {
  // c1_wnv = (-1)^{K+jw+jv} * sqrt[K] * {K ks kt; jn jw jv}
  // c2_wnv = (-1)^{kt+ks+jw+jv} * sqrt[K] * {K kt ks; jn jw jv}
  // ++ equation: {c1, c2} for (w, n, v) = (beta, n, a).
  // -- equation: {c2, c1} for (w, n, v) = (a, n, beta)
  const auto K = m_rank;
  const auto kt = m_T->rank();
  const auto ks = m_S->rank();
  const auto twoj_w = conj ? twoj_a : twoj_beta;
  const auto twoj_v = conj ? twoj_beta : twoj_a;

  const auto sqrtK_sign =
    std::sqrt(2.0 * K + 1.0) * Angular::neg1pow_2(twoj_w + twoj_v);

  const auto c1 =
    Angular::neg1pow(K) * sqrtK_sign *
    Angular::sixj_2(2 * K, 2 * ks, 2 * kt, twoj_n, twoj_w, twoj_v);

  const auto c2 =
    Angular::neg1pow(kt + ks) * sqrtK_sign *
    Angular::sixj_2(2 * K, 2 * kt, 2 * ks, twoj_n, twoj_w, twoj_v);

  return conj ? std::pair{c2, c1} : std::pair{c1, c2};
}

//==============================================================================

DiracSpinor TDHF_DCP::h_tilde_rhs(const TDHF &h_tilde, const int kappa,
                                  const DiracSpinor &F, bool conj) const {
  // The emission component of the operator is the stored one, up to the
  // signs of an anti-Hermitian and of an imaginary operator. With these,
  // and dV_rhs(conj), the result follows the convention of the Y
  // corrections
  const auto h = h_tilde.get_operator();
  const auto s_conj =
    conj ? h->h_minus_sign() * (h_tilde.imagQ() ? -1.0 : 1.0) : 1.0;
  auto rhs = s_conj * h->reduced_rhs(kappa, F);
  // Core polarisation of h, unless its first-order TDHF was not run, not
  // iterated (a single iteration: eps = 1, the change from zero), or failed.
  // A loose test: a solve that stops short of its target is still used
  if (h_tilde.last_eps() < 1.0e-2) {
    rhs += h_tilde.dV_rhs(kappa, F, conj);
  }
  return rhs;
}

//==============================================================================

double TDHF_DCP::h_tilde_de(const TDHF &h_tilde, const DiracSpinor &Fa,
                            bool conj) const {
  // Zero by the selection rules, or for an imaginary operator
  if (h_tilde.get_operator()->isZero(Fa, Fa) || h_tilde.imagQ())
    return 0.0;
  return Fa * h_tilde_rhs(h_tilde, Fa.kappa(), Fa, conj);
}

//==============================================================================

std::vector<std::vector<DiracSpinor>>
TDHF_DCP::form_hFcore(dPsiType XorY, bool include_d2V) const {
  const auto conj = XorY == dPsiType::Y;
  std::vector<std::vector<DiracSpinor>> hFcore(m_core.size());
#pragma omp parallel for schedule(dynamic)
  for (auto ic = 0ul; ic < m_core.size(); ic++) {
    const auto &Fc = m_core[ic];
    hFcore[ic].reserve(m_X[ic].size()); // one for each channel
    for (const auto &X_beta : m_X[ic]) {
      const auto kappa_beta = X_beta.kappa();
      hFcore[ic].push_back(source_rhs(Fc, XorY, kappa_beta));
      if (include_d2V) {
        // The parts of d2V that do not depend on dX, dY
        hFcore[ic].back() += d2V_11_rhs(kappa_beta, Fc, conj);
        hFcore[ic].back() += d2V_Norm_rhs(kappa_beta, Fc);
      }
    }
  }
  return hFcore;
}

//==============================================================================

std::pair<double, std::string> TDHF_DCP::tdhf_core_it(double omega,
                                                      bool static_Y) {
  using namespace qip::overloads;

  // Copies (instead of blank spinors): the mixed-state solves are faster
  // if they start from the existing solution
  auto dX_new = m_X;
  auto dY_new = m_Y;

  const auto eps_ms = 1.0e-12;

  // One task for each (core orbital, channel, dX or dY): many small, fairly
  // uniform tasks for the parallel loop. Each writes a distinct spinor, and
  // only reads shared state. See solve_ms_core_b().
  struct MixedStateTask {
    DiracSpinor *dF;           // target (in dX_new or dY_new)
    const DiracSpinor *Fa;     // core orbital
    const DiracSpinor *source; // fixed source of its equation
    dPsiType type;
  };
  std::vector<MixedStateTask> tasks;
  for (auto ia = 0ul; ia < m_core.size(); ia++) {
    for (auto ibeta = 0ul; ibeta < dX_new[ia].size(); ibeta++) {
      tasks.push_back(
        {&dX_new[ia][ibeta], &m_core[ia], &m_hFcore[ia][ibeta], dPsiType::X});
      if (!static_Y) {
        tasks.push_back({&dY_new[ia][ibeta], &m_core[ia],
                         &m_hFcore_minus[ia][ibeta], dPsiType::Y});
      }
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (auto i = 0ul; i < tasks.size(); i++) {
    const auto &task = tasks[i];
    solve_ms_core_b(*task.dF, *task.Fa, *task.source, omega, task.type, eps_ms);
  }

  // dX and dY are solved with the same d2V, so dX alone measures convergence
  const auto eps = eps_dPsi(m_core, m_X, dX_new, false);

  if (static_Y) {
    // Static: the -- equations are those of ++ up to a sign, from the
    // conjugate of the pair coupled to rank K, and of imaginary operators:
    // dY = (-1)^{K+kt+ks} * iota_t * iota_s * dX
    const auto sign_Y =
      Angular::neg1pow(m_rank + m_T->rank() + m_S->rank()) * (m_imag ? -1 : 1);
    for (auto ia = 0ul; ia < m_core.size(); ia++) {
      dY_new[ia] = sign_Y * dX_new[ia];
    }
  }

  m_X = std::move(dX_new);
  m_Y = std::move(dY_new);

  return eps;
}

//==============================================================================

void TDHF_DCP::solve_ms_core_b(DiracSpinor &dF_beta, const DiracSpinor &Fa,
                               const DiracSpinor &source, const double omega,
                               dPsiType XorY, double eps_ms) const {
  // Solves (H - e_a - w) dX = -rhs, or (H - e_a + w) dY = -rhs, with
  // rhs = source + [d2V_2 Fa]_beta. The d2V_2 part is from the corrections
  // of the previous iteration. The solver keeps the solution orthogonal to
  // Fa.
  // Thread-safe: reads only const/shared state; writes only dF_beta.

  const auto conj = XorY == dPsiType::Y;
  const auto ww = conj ? -omega : omega;

  const auto rhs = source + d2V_2_rhs(dF_beta.kappa(), Fa, conj);

  const auto vl = p_hf->vlocal(Angular::l_k(Fa.kappa()));
  const auto &Hmag = p_hf->Hmag(Angular::l_k(Fa.kappa()));
  ExternalField::solveMixedState_Anderson(dF_beta, Fa, ww, vl, m_alpha, m_core,
                                          rhs, eps_ms, p_VBr, Hmag);
}

//==============================================================================

void TDHF_DCP::solve_core(const double omega, int max_its, const bool print) {

  // The second-order corrections are at the sum of the frequencies of t and
  // s, each non-negative: the first-order X are the + components of both.
  // Any other frequency is allowed (it is solved as given), but is not the
  // physical one
  const auto omega_ts = m_T->last_omega() + m_S->last_omega();
  if (std::abs(omega - omega_ts) > 1.0e-10) {
    fmt2::warning();
    fmt::print(": TDHF_DCP: omega = {:.6f}, but omega_t + omega_s = {:.6f} "
               "(first-order solutions, with omega_t, omega_s >= 0)\n",
               omega, omega_ts);
  }

  const double converge_targ = m_eps;

  // Lowest order (a single iteration): d2V is not included in the core
  // equations. So none of it goes in the sources, and the corrections
  // start from zero
  const bool lowest_order = max_its <= 1;
  if (lowest_order) {
    clear();
  }

  // Static (both frequencies zero, and neither operator depends on the
  // frequency): dY follows from dX, so is not solved for; see
  // tdhf_core_it()
  const auto zero = 1.0e-10;
  const bool static_Y = std::abs(omega) < zero &&
                        std::abs(m_T->last_omega()) < zero &&
                        std::abs(m_S->last_omega()) < zero &&
                        !m_T->get_operator()->freqDependantQ() &&
                        !m_S->get_operator()->freqDependantQ();

  // Sources: fixed by the first-order solutions
  m_hFcore = form_hFcore(dPsiType::X, !lowest_order);
  if (!static_Y) {
    m_hFcore_minus = form_hFcore(dPsiType::Y, !lowest_order);
  }

  // Anderson mixing of all the corrections together, dX and dY (see
  // TDHF::solve_core_anderson): the equations are linear in dX, dY
  constexpr std::size_t history_depth = 8;
  AndersonVectorHistory history(history_depth);

  std::pair<double, std::string> eps{};
  int it{0};
  qip::LiveMessage status(
    fmt::format("TDHF_DCP {},{} K={} (w={:.4f}): ", m_T->get_operator()->name(),
                m_S->get_operator()->name(), m_rank, omega),
    print);
  for (; it < max_its; it++) {
    const auto x = dPsi_vector(m_X, m_Y);
    // Undamped map application: m_X, m_Y now hold G(x)
    eps = tdhf_core_it(omega, static_Y);

    status(fmt::format("{:2d} {:.1e} [{}]", it, eps.first, eps.second));

    if (eps.first < converge_targ || std::isnan(eps.first))
      break; // converged (or broken)

    // Anderson extrapolation over the kept history: x = sum_i c_i g_i
    history.push(x, dPsi_vector(m_X, m_Y));
    set_dPsi_vector(history.mix(), m_X, m_Y);
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

double TDHF_DCP::dV(const DiracSpinor &Fn, const DiracSpinor &Fm,
                    bool conj) const {
  // does it matter if a or b is in the core?
  const auto s = conj && m_imag ? -1 : 1; // careful. OK?
  return s * Fn * dV_rhs(Fn.kappa(), Fm, conj);
}

//==============================================================================

DiracSpinor TDHF_DCP::dV_rhs(const int kappa_w, const DiracSpinor &Fv,
                             bool conj) const {
  return d2V_11_rhs(kappa_w, Fv, conj) + d2V_Norm_rhs(kappa_w, Fv) +
         d2V_2_rhs(kappa_w, Fv, conj);
}

//==============================================================================

DiracSpinor TDHF_DCP::d2V_11_rhs(const int kappa_w, const DiracSpinor &Fv,
                                 bool conj) const {
  // Sum over core orbitals a, and pairs of first-order corrections to a,
  // one from each operator: the Y of one in the bra with the X of the other
  // in the ket. No explicit factors of iota: the Y are stored with theirs.
  //
  // conj (the -- component) is the Hermitian conjugate, found from
  // <w||d2V_--||v> = iota_t * iota_s * (-1)^{jw-jv} <v||d2V_++||w>, and
  // W^K_{abcd} = W^K_{cdab}: the roles of X and Y swap, with the phase
  // (-1)^{K+kt+ks} from the conjugate of the pair coupled to rank K.

  auto d2VFv = DiracSpinor(0, kappa_w, Fv.grid_sptr());
  d2VFv.max_pt() = Fv.max_pt();

  const auto K = m_rank;
  const auto kt = m_T->rank();
  const auto ks = m_S->rank();
  const auto tjw = Angular::twoj_k(kappa_w);
  const auto sign_st = Angular::neg1pow(K + ks + kt);
  const auto sign_conj = conj ? sign_st : 1;

  const auto in_ket = conj ? dPsiType::Y : dPsiType::X;
  const auto in_bra = conj ? dPsiType::X : dPsiType::Y;

  for (const auto &Fa : m_core) {
    const auto tja = Fa.twoj();
    const auto sign_wa = sign_conj * Angular::neg1pow_2(tjw + tja);

    // First-order corrections to a. ++ component: X in the ket, Y in the
    // bra; -- component: Y in the ket, X in the bra
    const auto &ket_t = m_T->get_dPsis(Fa, in_ket);
    const auto &bra_t = m_T->get_dPsis(Fa, in_bra);
    const auto &ket_s = m_S->get_dPsis(Fa, in_ket);
    const auto &bra_s = m_S->get_dPsis(Fa, in_bra);

    // s correction in the bra, t correction in the ket
    for (const auto &bra_s_beta : bra_s) {
      for (const auto &ket_t_alpha : ket_t) {
        const auto sixj = Angular::sixj_2(
          2 * kt, 2 * ks, 2 * K, bra_s_beta.twoj(), ket_t_alpha.twoj(), tja);
        if (sixj != 0.0) {
          d2VFv += (sign_wa * sixj) *
                   Coulomb::Wkv_bcd(K, kappa_w, bra_s_beta, Fv, ket_t_alpha);
        }
      }
    }

    // t correction in the bra, s correction in the ket
    for (const auto &bra_t_beta : bra_t) {
      for (const auto &ket_s_alpha : ket_s) {
        const auto sixj = Angular::sixj_2(
          2 * kt, 2 * ks, 2 * K, ket_s_alpha.twoj(), bra_t_beta.twoj(), tja);
        if (sixj != 0.0) {
          d2VFv += (sign_st * sign_wa * sixj) *
                   Coulomb::Wkv_bcd(K, kappa_w, bra_t_beta, Fv, ket_s_alpha);
        }
      }
    }
  }

  d2VFv *= (1.0 / std::sqrt(2.0 * K + 1.0));

  return d2VFv;
}

//==============================================================================

DiracSpinor TDHF_DCP::d2V_Norm_rhs(const int kappa_w,
                                   const DiracSpinor &Fv) const {
  // The second-order corrections are solved orthogonal to their core
  // orbital, a. Normalisation fixes the component along a, from the
  // overlaps of the first-order corrections to a (Y of one operator with X
  // of the other, same channel). It enters as that component times
  // W^K_{wava}. No explicit factors of iota: the Y are stored with theirs.
  //
  // The same for the ++ and -- components: for --, X and Y swap, with the
  // phase (-1)^{K+kt+ks} (see d2V_11_rhs), which leaves the sum unchanged.

  auto d2VFv = DiracSpinor(0, kappa_w, Fv.grid_sptr());
  d2VFv.max_pt() = Fv.max_pt();

  const auto K = m_rank;
  const auto kt = m_T->rank();
  const auto ks = m_S->rank();
  const auto tjw = Angular::twoj_k(kappa_w);
  const auto sign_ts = Angular::neg1pow(K + ks + kt);

  for (const auto &Fa : m_core) {
    const auto tja = Fa.twoj();

    // W^K_{wava} is zero unless (ja, ja, K) is a triangle
    if (Angular::triangle(tja, tja, 2 * K) == 0)
      continue;

    const auto &Xt = m_T->get_dPsis(Fa, dPsiType::X);
    const auto &Yt = m_T->get_dPsis(Fa, dPsiType::Y);
    const auto &Xs = m_S->get_dPsis(Fa, dPsiType::X);
    const auto &Ys = m_S->get_dPsis(Fa, dPsiType::Y);

    // Sum over the channels, nu, common to the corrections from t and s
    double overlaps = 0.0;
    for (auto it = 0ul; it < Xt.size(); ++it) {
      for (auto is = 0ul; is < Xs.size(); ++is) {
        if (Xt[it].kappa() != Xs[is].kappa())
          continue;
        const auto tjnu = Xt[it].twoj();
        const auto sixj =
          Angular::sixj_2(2 * kt, 2 * ks, 2 * K, tja, tja, tjnu);
        if (sixj == 0.0)
          continue;
        overlaps += Angular::neg1pow_2(tjw - tjnu) * sixj *
                    (sign_ts * (Ys[is] * Xt[it]) + Yt[it] * Xs[is]);
      }
    }

    if (overlaps != 0.0) {
      d2VFv += overlaps * Coulomb::Wkv_bcd(K, kappa_w, Fa, Fv, Fa);
    }
  }

  d2VFv *= (1.0 / std::sqrt(2.0 * K + 1.0));

  return d2VFv;
}

//==============================================================================

DiracSpinor TDHF_DCP::d2V_2_rhs(const int kappa_w, const DiracSpinor &Fv,
                                bool conj) const {
  // The first-order formula for dV, at rank K, with the second-order
  // corrections dX, dY in place of X, Y

  auto d2VFv = DiracSpinor(0, kappa_w, Fv.grid_sptr());
  d2VFv.max_pt() = Fv.max_pt();

  const auto K = m_rank;
  const auto tKp1 = double(2 * K + 1);
  const auto tjw = Angular::twoj_k(kappa_w);

  // nb: faster to not //ize this one
  for (auto ia = 0ul; ia < m_core.size(); ++ia) {
    const auto &Fa = m_core[ia];

    // conj (the -- component): the roles of dX and dY swap
    const auto &dX_a = conj ? m_Y[ia] : m_X[ia];
    const auto &dY_a = conj ? m_X[ia] : m_Y[ia];

    for (auto ibeta = 0ul; ibeta < dX_a.size(); ++ibeta) {
      const auto &dX_beta = dX_a[ibeta];
      const auto &dY_beta = dY_a[ibeta];

      const auto sQ = Angular::neg1pow_2(tjw - dX_beta.twoj() + 2 * K);
      if (sQ == 1) { // faster than multiplying!
        d2VFv += Coulomb::Wkv_bcd(K, kappa_w, Fa, Fv, dX_beta);
        d2VFv += Coulomb::Wkv_bcd(K, kappa_w, dY_beta, Fv, Fa);
      } else {
        d2VFv -= Coulomb::Wkv_bcd(K, kappa_w, Fa, Fv, dX_beta);
        d2VFv -= Coulomb::Wkv_bcd(K, kappa_w, dY_beta, Fv, Fa);
      }

      // Breit part:
      if (p_VBr) {
        d2VFv += tKp1 * p_VBr->dV_Br(kappa_w, K, Fv, Fa, dX_beta, dY_beta);
      }
    }
  }

  d2VFv *= (1.0 / tKp1);

  return d2VFv;
}

} // namespace ExternalField
