#include "TDHFcomplex.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/ContinuumState.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/HartreeFock.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "qip/Vector.hpp"
#include "qip/Widgets.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <utility>
#include <vector>

namespace ExternalField {

//==============================================================================
TDHFcntm::TDHFcntm(const DiracOperator::TensorOperator *const h_plus,
                   const HF::HartreeFock *const hf,
                   const DiracOperator::TensorOperator *const h_minus)
  : TDHF(h_plus, hf, h_minus) {
  zero_imaginary_sets();
}

//==============================================================================
void TDHFcntm::zero_imaginary_sets() {
  using namespace qip::overloads;
  m_Xi = m_X;
  m_Yi = m_Y;
  m_Xi *= 0.0;
  m_Yi *= 0.0;
}

//==============================================================================
void TDHFcntm::clear() {
  TDHF::clear();
  zero_imaginary_sets();
  m_channels.clear();
  m_omega = -1.0;
}

//==============================================================================
void TDHFcntm::prepare_channels(double omega) {
  // For each open channel (en_+ = en_b + omega > 0):
  //  1. the homogeneous pair Freg/Firr at en_+, in the fixed local channel
  //     potential v = vlocal(l) - y^0_bb + U_KS (the potential
  //     solveContinuumMixedState uses, so the pair is reused across the
  //     iterations);
  //  2. the exchange-dressed regular solution F^HF: Freg plus the standing
  //     response to its exchange deficit (vexFa - vexFa_1el - U_KS) Freg,
  //     F^HF ~ Freg + K_ex Firr. Diagonal channels (kappa = kappa_b: even
  //     parity operators) orthogonalise the incident wave against phi_b,
  //     with the matching deficit term c (omega + V^b_0) phi_b [from
  //     h_HF phi_b = en_b phi_b], so the outgoing solution stays orthogonal
  //     to phi_b; phi_b decays, so the asymptotics are unchanged.
  // A channel whose continuum solve returns zero (the grid does not resolve
  // en_+) is left closed, and is zeroed in the iteration.
  using namespace qip::overloads;

  m_channels.clear();
  m_channels.resize(m_X.size());
  m_omega = omega;

  const auto Ux = HF::vex_KS(m_core);

  std::vector<std::pair<std::size_t, std::size_t>> open_list;
  for (std::size_t ib = 0; ib < m_X.size(); ++ib) {
    const auto &Fb = m_core[ib];
    const auto en_plus = Fb.en() + omega;
    const bool ionised = en_plus > 0.0;
    const auto y0bb =
      ionised ? Coulomb::yk_ab(0, Fb, Fb) : std::vector<double>{};

    for (const auto &X_beta : m_X[ib]) {
      const auto kappa = X_beta.kappa();
      ContinuumChannel channel{false,
                               DiracSpinor{0, kappa, Fb.grid_sptr()},
                               DiracSpinor{0, kappa, Fb.grid_sptr()},
                               DiracSpinor{0, kappa, Fb.grid_sptr()},
                               0.0,
                               0.0};
      if (ionised) {
        const auto v = p_hf->vlocal(Angular::l_k(kappa)) - y0bb + Ux;
        DiracODE::solveContinuum(channel.Freg, en_plus, v, m_alpha);
        if (channel.Freg.norm2() == 0.0) {
          fmt::print("Warning: TDHFcntm: continuum solve failed for {} "
                     "(kappa = {}) at en = {:.4f}; grid too coarse. Channel "
                     "excluded.\n",
                     Fb.shortSymbol(), kappa, en_plus);
        } else {
          DiracODE::solveContinuumIrregular(channel.Firr, channel.Freg, en_plus,
                                            v, m_alpha);
          channel.open = true;
          open_list.emplace_back(ib, m_channels[ib].size());
        }
      }
      m_channels[ib].push_back(std::move(channel));
    }
  }

  const auto eps_ms = 1.0e-12;
#pragma omp parallel for schedule(dynamic)
  for (std::size_t i = 0; i < open_list.size(); ++i) {
    const auto [ib, be] = open_list[i];
    auto &channel = m_channels[ib][be];
    const auto &Fa = m_core[ib];
    const auto &Freg = channel.Freg;
    const auto kappa = Freg.kappa();

    auto deficit =
      HF::vexFa(Freg, m_core) - HF::vexFa_1el(Freg, Fa) - (Ux * Freg);
    auto incident = Freg;
    if (kappa == Fa.kappa()) {
      const auto c_a = Fa * incident;
      incident -= c_a * Fa;
      deficit += c_a * (omega * Fa + hole_compensation(Fa, Fa));
    }

    const auto y0aa = Coulomb::yk_ab(0, Fa, Fa);
    const auto vl_c = p_hf->vlocal(Angular::l_k(kappa)) - y0aa;
    DiracSpinor correction{0, kappa, Fa.grid_sptr()};
    double K_ex = 0.0;
    // The cached pair is reused (same energy), not rebuilt
    ExternalField::solveContinuumMixedState(
      &correction, &channel.Freg, &channel.Firr, &K_ex, Fa, omega, vl_c,
      m_alpha, m_core, deficit, eps_ms, &Fa);
    channel.Fhf = incident + correction;
    channel.K_ex = K_ex;
    channel.K = 0.0;
  }
}

//==============================================================================
void TDHFcntm::solve_core(double omega, int max_its, bool print) {

  assert(m_h->rank() == m_rank && "Rank must match in solve_core");
  assert(m_h->parity() == m_pi && "Parity must match in solve_core");
  assert(m_h->imaginaryQ() == m_imag && "Imaginarity must match in solve_core");

  const double converge_targ = m_eps;
  omega = std::abs(omega);

  m_hFcore = form_hFcore(m_h);
  m_hFcore_minus = form_hFcore(m_h_minus);

  // Continuum data of the open channels; kept while omega is unchanged (a
  // re-solve at the same omega warm-starts from the previous corrections)
  if (omega != m_omega || m_channels.empty()) {
    prepare_channels(omega);
  }

  // Anderson history of the whole state: map outputs g_k = G(x_k)
  // and residuals r_k = g_k - x_k. Each entry holds two full copies of the
  // state (the dominant memory use); depth 4 converges as well as deeper
  // histories here.
  constexpr std::size_t history_depth = 4;
  std::vector<std::vector<double>> g_hist;
  std::vector<std::vector<double>> r_hist;

  std::pair<double, std::string> eps{};
  double best_eps{1.0e30};
  int count_worse = 0;
  int it{0};
  qip::LiveMessage status(
    fmt::format("TDHFcntm {} (w={:.4f}): ", m_h->name(), omega), print);
  for (; it < max_its; it++) {
    const auto x = state_vector();
    // Undamped map application: the sets now hold G(x); eps measures G(x)
    // against x, the residual in the physical measure
    eps = tdhf_core_it_complex(omega);

    status(fmt::format("{:2d} {:.1e} [{}]", it, eps.first, eps.second));

    if (eps.first < converge_targ || std::isnan(eps.first))
      break; // converged (or broken)

    auto g = state_vector();
    auto r = g;
    for (std::size_t i = 0; i < r.size(); ++i) {
      r[i] -= x[i];
    }
    g_hist.push_back(std::move(g));
    r_hist.push_back(std::move(r));
    if (g_hist.size() > history_depth) {
      g_hist.erase(g_hist.begin());
      r_hist.erase(r_hist.begin());
    }

    // Anderson extrapolation over the kept history: x = sum_i c_i g_i
    const auto m = r_hist.size();
    std::vector<std::vector<double>> gram(m, std::vector<double>(m, 0.0));
    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        double sum = 0.0;
        for (std::size_t k = 0; k < r_hist[i].size(); ++k) {
          sum += r_hist[i][k] * r_hist[j][k];
        }
        gram[i][j] = sum;
        gram[j][i] = sum;
      }
    }
    const auto c = anderson_coefficients(gram);
    const auto n_drop = long(m - c.size());
    g_hist.erase(g_hist.begin(), g_hist.begin() + n_drop);
    r_hist.erase(r_hist.begin(), r_hist.begin() + n_drop);
    // Nothing to mix (a single entry is the plain step, which the sets
    // already hold)
    if (c.size() <= 1)
      continue;
    std::vector<double> x_new(g_hist[0].size(), 0.0);
    for (std::size_t i = 0; i < c.size(); ++i) {
      for (std::size_t k = 0; k < x_new.size(); ++k) {
        x_new[k] += c[i] * g_hist[i][k];
      }
    }
    set_state(x_new);

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

  m_core_eps = eps.first;
  m_core_its = it;
  m_core_omega = omega;
}

//==============================================================================
std::vector<double> TDHFcntm::state_vector() const {
  std::vector<double> state;
  for (const auto *set : {&m_X, &m_Y, &m_Xi, &m_Yi}) {
    for (const auto &Fs : *set) {
      for (const auto &F : Fs) {
        state.insert(state.end(), F.f().cbegin(), F.f().cend());
        state.insert(state.end(), F.g().cbegin(), F.g().cend());
      }
    }
  }
  for (const auto &orbital_channels : m_channels) {
    for (const auto &channel : orbital_channels) {
      state.push_back(channel.K.real());
      state.push_back(channel.K.imag());
    }
  }
  return state;
}

//==============================================================================
void TDHFcntm::set_state(const std::vector<double> &state) {
  std::size_t i = 0;
  for (auto *set : {&m_X, &m_Y, &m_Xi, &m_Yi}) {
    for (auto &Fs : *set) {
      for (auto &F : Fs) {
        std::copy_n(state.cbegin() + long(i), F.f().size(), F.f().begin());
        i += F.f().size();
        std::copy_n(state.cbegin() + long(i), F.g().size(), F.g().begin());
        i += F.g().size();
      }
    }
  }
  for (auto &orbital_channels : m_channels) {
    for (auto &channel : orbital_channels) {
      channel.K = {state[i], state[i + 1]};
      i += 2;
    }
  }
  assert(i == state.size());
}

//==============================================================================
std::pair<double, std::string> TDHFcntm::tdhf_core_it_complex(double omega) {
  // One undamped iteration of the complex TDHF map. The dV sources of the
  // real and imaginary sets are built separately (the builder is
  // real-linear; the sets couple only through the outgoing addition in the
  // open channels). The external field drives the real part only.
  assert(omega >= 0.0 && "solve_core() passes omega = |omega|");
  assert(!m_channels.empty() && "iteration requires prepared channels");

  auto Xs = m_X;
  auto Ys = m_Y;
  auto Xis = m_Xi;
  auto Yis = m_Yi;

  const auto eps_ms = 1.0e-12;

  // Previous iteration's amplitudes (for the convergence measure); the
  // channel solves write the new K into m_channels in place
  std::vector<std::vector<std::complex<double>>> K_old(m_channels.size());
  for (std::size_t ib = 0; ib < m_channels.size(); ++ib) {
    for (const auto &channel : m_channels[ib]) {
      K_old[ib].push_back(channel.K);
    }
  }

  // Flatten all (core orbital x channel x X/Y) solves into one task list
  // (load balance). Only the X tasks of an ionised orbital carry a
  // continuum channel; the Y partners are always bound.
  struct Task {
    DiracSpinor *dF_re;        // targets (in Xs/Ys and Xis/Yis)
    DiracSpinor *dF_im;        //
    ContinuumChannel *channel; // continuum data (nullptr if bound)
    const DiracSpinor *Fb;     // core orbital
    const DiracSpinor *hFb;    // external source h|Fb> (real part only)
    dPsiType type;
  };
  std::vector<Task> tasks;
  for (std::size_t ib = 0; ib < m_core.size(); ++ib) {
    const bool ionised = m_core[ib].en() + omega > 0.0;
    for (std::size_t be = 0; be < Xs[ib].size(); ++be) {
      tasks.push_back({&Xs[ib][be], &Xis[ib][be],
                       ionised ? &m_channels[ib][be] : nullptr, &m_core[ib],
                       &m_hFcore[ib][be], dPsiType::X});
      tasks.push_back({&Ys[ib][be], &Yis[ib][be], nullptr, &m_core[ib],
                       &m_hFcore_minus[ib][be], dPsiType::Y});
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t it = 0; it < tasks.size(); ++it) {
    const auto &t = tasks[it];
    const auto conj = t.type == dPsiType::Y;
    const auto kappa = t.dF_re->kappa();
    // [dV phi_b]_beta from the previous iterate: real and imaginary sets
    const auto dV_re = dV_rhs_sets(kappa, *t.Fb, conj, m_X, m_Y);
    const auto dV_im = dV_rhs_sets(kappa, *t.Fb, conj, m_Xi, m_Yi);
    if (t.channel == nullptr) {
      // Bound channel: two independent real solves, the external field in
      // the real part only
      solve_channel_bound(t.dF_re, *t.Fb, t.hFb, dV_re, omega, t.type, eps_ms);
      solve_channel_bound(t.dF_im, *t.Fb, nullptr, dV_im, omega, t.type,
                          eps_ms);
    } else if (t.channel->open) {
      solve_channel_outgoing(t.dF_re, t.dF_im, t.channel, *t.Fb, *t.hFb, dV_re,
                             dV_im, omega, eps_ms);
    } else {
      // Open channel the grid cannot resolve: excluded (explicit zero)
      t.dF_re->f().assign(t.dF_re->f().size(), 0.0);
      t.dF_re->g().assign(t.dF_re->g().size(), 0.0);
      t.dF_im->f().assign(t.dF_im->f().size(), 0.0);
      t.dF_im->g().assign(t.dF_im->g().size(), 0.0);
      t.channel->K = 0.0;
    }
  }

  const auto eps = eps_complex(Xs, Xis, K_old);

  m_X = std::move(Xs);
  m_Y = std::move(Ys);
  m_Xi = std::move(Xis);
  m_Yi = std::move(Yis);

  return eps;
}

//==============================================================================
void TDHFcntm::solve_channel_outgoing(
  DiracSpinor *X_re, DiracSpinor *X_im, ContinuumChannel *channel,
  const DiracSpinor &Fb, const DiracSpinor &hFb, const DiracSpinor &dV_re,
  const DiracSpinor &dV_im, double omega, double eps_ms) const {
  // Open X channel of ionised orbital Fb, en_+ = en_b + omega > 0, with the
  // outgoing-wave boundary condition (see class description). The
  // one-electron self-interaction V^b_0 = y^0_bb + X_b of the hole is
  // subtracted from the static Hamiltonian (local potential vl - y^0_bb;
  // exchange part via the solver's Fhole option) and compensated in the
  // source from the previous iterate: an exact rearrangement, whose source
  // is short ranged.
  using namespace qip::overloads;
  assert(channel != nullptr && channel->open);
  const int kappa_beta = X_re->kappa();

  // Physical sources: [(t + dV)phi_b]_beta. The external field is real in
  // the code's convention, so it drives the real part only.
  auto rhs_re = dV_re + hFb;
  auto rhs_im = dV_im;
  // Diagonal channel (even-parity operators only): subtract the
  // norm-conservation term de*phi_b, de = <b|(t + dV)phi_b>. Applies to the
  // physical source only: the V^b_0 compensation below is part of the
  // operator rearrangement, and projecting it too would over-subtract.
  if (kappa_beta == Fb.kappa()) {
    rhs_re -= (Fb * rhs_re) * Fb;
    rhs_im -= (Fb * rhs_im) * Fb;
  }
  // The compensation from the full previous complex iterate (X_re, X_im
  // still hold it)
  rhs_re += hole_compensation(Fb, *X_re);
  rhs_im += hole_compensation(Fb, *X_im);

  const auto y0bb = Coulomb::yk_ab(0, Fb, Fb);
  const auto vl_c = p_hf->vlocal(Angular::l_k(kappa_beta)) - y0bb;

  // Warm start of the inner (exchange) iterations from the standing parts
  // of the previous iterate: X = phi - i K F^HF, so
  // phi_re = X_re - Im(K) F^HF and phi_im = X_im + Re(K) F^HF
  const auto &Fhf = channel->Fhf;
  *X_re -= channel->K.imag() * Fhf;
  *X_im += channel->K.real() * Fhf;

  // Two standing-wave solves on the cached pair (K amplitudes out)
  double K_re = 0.0;
  double K_im = 0.0;
  ExternalField::solveContinuumMixedState(X_re, &channel->Freg, &channel->Firr,
                                          &K_re, Fb, omega, vl_c, m_alpha,
                                          m_core, rhs_re, eps_ms, &Fb);
  ExternalField::solveContinuumMixedState(X_im, &channel->Freg, &channel->Firr,
                                          &K_im, Fb, omega, vl_c, m_alpha,
                                          m_core, rhs_im, eps_ms, &Fb);

  // Outgoing amplitude and addition: K_+ = (K_r + i K_i)/(1 + i K_ex),
  // X = phi - i K_+ F^HF
  const std::complex<double> K_hf{K_re, K_im};
  const std::complex<double> i_unit{0.0, 1.0};
  channel->K = K_hf / (1.0 + i_unit * channel->K_ex);
  *X_re += channel->K.imag() * Fhf;
  *X_im -= channel->K.real() * Fhf;
}

//==============================================================================
void TDHFcntm::solve_channel_bound(DiracSpinor *dF_beta, const DiracSpinor &Fb,
                                   const DiracSpinor *hFb,
                                   const DiracSpinor &dV_src, double omega,
                                   dPsiType XorY, double eps_ms) const {
  // Bound channel (closed orbital, or the Y partner of an ionised orbital):
  // plain V^N bound solve, as TDHF, with Anderson mixing of the mixed-state
  // iteration (solveMixedState_cntm). The V^{N-1} rearrangement applies to
  // the open channels only.
  using namespace qip::overloads;
  const auto ww = XorY == dPsiType::X ? omega : -omega;
  const auto conj = XorY == dPsiType::Y;
  const auto s = (m_h->imaginaryQ() && conj) ? -1.0 : 1.0;
  auto rhs = dV_src;
  if (hFb) {
    rhs += s * (*hFb);
  }
  const auto vl = p_hf->vlocal(Angular::l_k(Fb.kappa()));
  const auto &Hmag = p_hf->Hmag(Angular::l_k(Fb.kappa()));
  ExternalField::solveMixedState_cntm(*dF_beta, Fb, ww, vl, m_alpha, m_core,
                                      rhs, eps_ms, p_VBr, Hmag);
}

//==============================================================================
DiracSpinor TDHFcntm::hole_compensation(const DiracSpinor &Fa,
                                        const DiracSpinor &chi) const {
  using namespace qip::overloads;
  const auto y0aa = Coulomb::yk_ab(0, Fa, Fa);
  return (y0aa * chi) + HF::vexFa_1el(chi, Fa);
}

//==============================================================================
std::pair<double, std::string> TDHFcntm::eps_complex(
  const std::vector<std::vector<DiracSpinor>> &Xs_re,
  const std::vector<std::vector<DiracSpinor>> &Xs_im,
  const std::vector<std::vector<std::complex<double>>> &K_old) const {
  double DdF2 = 0.0;
  double dF2 = 0.0;
  double worst = 0.0;
  double worst_K = 0.0;
  std::string worst_lab;

  // Floor for the open-channel measure: the largest |K_+|^2 (a channel
  // passing through zero must not gate the iteration)
  double K_max2 = 0.0;
  for (const auto &orbital_channels : m_channels) {
    for (const auto &channel : orbital_channels) {
      K_max2 = std::max(K_max2, std::norm(channel.K));
    }
  }
  const auto K_floor2 = 1.0e-4 * K_max2;

  for (std::size_t ib = 0; ib < m_core.size(); ++ib) {
    for (std::size_t i = 0; i < Xs_re[ib].size(); ++i) {
      const auto &nw_re = Xs_re[ib][i];
      const auto &nw_im = Xs_im[ib][i];
      // Norm scale: ALL X channels contribute to the denominator
      dF2 += nw_re.norm2() + nw_im.norm2();
      const auto lab = m_core[ib].shortSymbol() + "," + nw_re.shortSymbol();
      if (m_channels[ib][i].open) {
        const auto Kn = m_channels[ib][i].K;
        const auto dK2 = std::norm(Kn - K_old[ib][i]);
        const auto den = std::max(std::norm(Kn), K_floor2);
        const auto e = den == 0.0 ? 0.0 : dK2 / den;
        worst_K = std::max(worst_K, e);
        if (e > worst) {
          worst = e;
          worst_lab = lab;
        }
      } else {
        const auto d =
          (nw_re - m_X[ib][i]).norm2() + (nw_im - m_Xi[ib][i]).norm2();
        const auto n = nw_re.norm2() + nw_im.norm2();
        DdF2 += d;
        // per-channel ratio: used only to rank the worst channel
        const auto e = n == 0.0 ? 0.0 : d / n;
        if (e > worst) {
          worst = e;
          worst_lab = lab;
        }
      }
    }
  }
  const auto ratio = dF2 == 0.0 ? 0.0 : DdF2 / dF2;
  const auto eps = std::max(ratio, worst_K);
  return {m_eps_sqrt ? std::sqrt(eps) : eps, worst_lab};
}

//==============================================================================
std::vector<TDHFcntm::Channel> TDHFcntm::channel_list() const {
  std::vector<Channel> channels;
  for (std::size_t ib = 0; ib < m_channels.size(); ++ib) {
    for (std::size_t be = 0; be < m_channels[ib].size(); ++be) {
      if (m_channels[ib][be].open) {
        channels.push_back(
          {ib, m_X[ib][be].kappa(), m_core[ib].en() + m_omega});
      }
    }
  }
  return channels;
}

//==============================================================================
std::pair<std::size_t, std::size_t>
TDHFcntm::channel_index(const DiracSpinor &Fb, int kappa) const {
  const auto ib = static_cast<std::size_t>(
    std::find(m_core.cbegin(), m_core.cend(), Fb) - m_core.cbegin());
  assert(ib < m_X.size() && "Fb must be a core orbital");
  const auto be = static_cast<std::size_t>(
    std::find_if(m_X[ib].cbegin(), m_X[ib].cend(),
                 [kappa](const auto &X) { return X.kappa() == kappa; }) -
    m_X[ib].cbegin());
  assert(be < m_X[ib].size() && "kappa must be one of Fb's channels");
  return {ib, be};
}

//==============================================================================
std::vector<std::complex<double>> TDHFcntm::A_phys() const {
  std::vector<std::complex<double>> A;
  for (const auto &orbital_channels : m_channels) {
    for (const auto &channel : orbital_channels) {
      if (channel.open) {
        A.push_back(std::conj(channel.K) / M_PI);
      }
    }
  }
  return A;
}

//==============================================================================
std::complex<double> TDHFcntm::A_phys(const DiracSpinor &Fa, int kappa) const {
  if (m_channels.empty()) {
    return 0.0;
  }
  const auto [ib, be] = channel_index(Fa, kappa);
  if (!m_channels[ib][be].open) {
    return 0.0;
  }
  return std::conj(m_channels[ib][be].K) / M_PI;
}

//==============================================================================
std::complex<double> TDHFcntm::dV_complex(const DiracSpinor &Fa,
                                          const DiracSpinor &Fb) const {
  // As TDHF::dV, on the real and imaginary sets. A continuum bra (en_a > 0)
  // gets the V^{N-1} treatment of the ejected electron: the source carries
  // the hole term V^b_0 chi, whose y^0_bb chi part cancels pointwise the
  // 1/r tail of the b-diagonal part of dV phi_b, so the continuum overlap
  // is box independent.
  using namespace qip::overloads;
  assert(Fb.en() <= 0.0 && "the ket must be a core orbital: dV(Fe, Fa)");
  const auto conj = Fb.en() > Fa.en();
  const auto s = conj && m_h->imaginaryQ() ? -1.0 : 1.0;
  auto S_re = dV_rhs(Fa.kappa(), Fb, conj);
  auto S_im = dV_rhs_sets(Fa.kappa(), Fb, conj, m_Xi, m_Yi);
  if (Fa.en() > 0.0) {
    assert(!conj);
    const auto [ib, be] = channel_index(Fb, Fa.kappa());
    S_re += hole_compensation(Fb, m_X[ib][be]);
    S_im += hole_compensation(Fb, m_Xi[ib][be]);
  }
  return s * std::complex<double>{Fa * S_re, Fa * S_im};
}

//==============================================================================
double TDHFcntm::dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  assert(Fa.en() <= 0.0 && Fb.en() <= 0.0 &&
         "TDHFcntm::dV is real: for a continuum state use dV_complex");
  return dV_complex(Fa, Fb).real();
}

} // namespace ExternalField
