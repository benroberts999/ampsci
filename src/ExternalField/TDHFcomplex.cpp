#include "TDHFcomplex.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/ContinuumState.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "LinAlg/Solvers.hpp"
#include "LinAlg/Vector.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "qip/Maths.hpp"
#include "qip/Vector.hpp"
#include "qip/Widgets.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <utility>
#include <vector>

namespace ExternalField {

//==============================================================================
TDHFcomplex::TDHFcomplex(const DiracOperator::TensorOperator *const h_plus,
                         const HF::HartreeFock *const hf,
                         const DiracOperator::TensorOperator *const h_minus)
  : TDHF(h_plus, hf, h_minus),
    // Fixed core-core y^l(Fi, Fj) table for the exchange part of
    // dV_rhs_all(); the core never changes, so build it once here (copies
    // share it)
    m_ycc(std::make_shared<const Coulomb::YkTable>(m_core)) {
  zero_imaginary_sets();
}

//==============================================================================
void TDHFcomplex::zero_imaginary_sets() {
  using namespace qip::overloads;
  m_Xi = m_X;
  m_Yi = m_Y;
  m_Xi *= 0.0;
  m_Yi *= 0.0;
}

//==============================================================================
void TDHFcomplex::clear() {
  TDHF::clear();
  zero_imaginary_sets();
  m_channels.clear();
  m_have_solution = false;
  m_omega = -1.0;
  m_KpiD = 0.0;
  m_KpiD_lab.clear();
  m_excluded.clear();
}

//==============================================================================
void TDHFcomplex::set_operator(const DiracOperator::TensorOperator *h_plus,
                               const DiracOperator::TensorOperator *h_minus) {
  assert(h_plus != nullptr);
  // The omega-only caches depend on the operator only through the channel
  // structure (rank, parity): keep them when that is unchanged, else
  // rebuild the channel sets and drop the caches
  const bool same_channels =
    h_plus->rank() == m_rank && h_plus->parity() == m_pi;
  m_h = h_plus;
  m_h_minus = h_minus ? h_minus : h_plus;
  m_rank = h_plus->rank();
  m_pi = h_plus->parity();
  m_imag = h_plus->imaginaryQ();
  if (same_channels) {
    TDHF::clear();
    zero_imaginary_sets();
    m_have_solution = false;
  } else {
    m_X.clear();
    m_Y.clear();
    initialise_dPsi();
    clear();
  }
}

//==============================================================================
void TDHFcomplex::prepare(double omega) {
  omega = std::abs(omega);
  if (omega != m_omega || m_channels.empty()) {
    prepare_channels(omega);
  }
}

//==============================================================================
std::vector<TDHFcomplex::Channel> TDHFcomplex::channel_list() const {
  std::vector<Channel> channels;
  if (m_channels.empty() || m_omega < 0.0) {
    return channels;
  }
  list_open_channels(&channels);
  return channels;
}

//==============================================================================
void TDHFcomplex::prepare_channels(double omega) {
  // (Re)build the per-channel continuum caches for this omega. For each
  // open channel (en_+ = en_b + omega > 0):
  //  1. the homogeneous pair Freg/Firr at en_+, built ONCE in the fixed
  //     conditioning potential v = vlocal(l) - y^0_bb + U_KS (the same
  //     potential solveContinuumMixedState() uses, so it reuses them
  //     across the iterations; barely-open channels are handled by the
  //     outward-extension construction of the pair, high-energy channels
  //     by the auxiliary-grid solve and truncation at Freg.max_pt());
  //  2. the exchange-dressed incident wave F^HF: Freg (an eigenfunction of
  //     the local conditioning potential) plus the standing response to
  //     its exchange deficit (vexFa - vexFa_1el - U_KS) Freg, with
  //     F^HF ~ Freg + K_ex Firr. Diagonal channels (even-parity operators)
  //     orthogonalise the incident wave against phi_b, with the matching
  //     deficit term c*(omega + V^b_0) phi_b (from h_HF phi_b = en_b phi_b),
  //     so the whole outgoing solution stays orthogonal to phi_b; phi_b
  //     decays, so the asymptotics are unchanged.
  // The only exclusion is the backstop: solveContinuum returning zero
  // (nothing resolvable on this grid) -- a zero pair must never reach the
  // extraction, which divides by the pair Wronskian.
  using namespace qip::overloads;

  m_channels.clear();
  m_channels.resize(m_X.size());
  m_omega = omega;
  m_have_solution = false;
  m_excluded.clear();

  const auto Ux = HF::vex_KS(m_core);

  std::vector<std::pair<std::size_t, std::size_t>> open_list;
  for (auto ib = 0ul; ib < m_X.size(); ib++) {
    const auto &Fb = m_core[ib];
    const auto en_plus = Fb.en() + omega;
    const bool open = en_plus > 0.0;

    const auto y0bb = open ? Coulomb::yk_ab(0, Fb, Fb) : std::vector<double>{};

    std::string failed{};
    for (const auto &X_beta : m_X[ib]) {
      const auto kappa = X_beta.kappa();
      ContinuumChannel channel{false,
                               DiracSpinor{0, kappa, Fb.grid_sptr()},
                               DiracSpinor{0, kappa, Fb.grid_sptr()},
                               DiracSpinor{0, kappa, Fb.grid_sptr()},
                               0.0,
                               0.0};
      if (open) {
        const auto v = p_hf->vlocal(Angular::l_k(kappa)) - y0bb + Ux;
        DiracODE::solveContinuum(channel.Freg, en_plus, v, m_alpha);
        if (channel.Freg.norm2() == 0.0) {
          failed += (failed.empty() ? "" : ",") + X_beta.shortSymbol();
        } else {
          DiracODE::solveContinuumIrregular(channel.Firr, channel.Freg, en_plus,
                                            v, m_alpha);
          channel.open = true;
          open_list.emplace_back(ib, m_channels[ib].size());
        }
      }
      // nb: "open" in the cache means open AND usable: failed channels are
      // dispatched exactly like suppress_open (zeroed)
      m_channels[ib].push_back(std::move(channel));
    }

    if (!failed.empty()) {
      m_excluded +=
        (m_excluded.empty() ? "" : "; ") + Fb.shortSymbol() + "," + failed;
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
TDHFcomplex::dVrhs_XY
TDHFcomplex::dV_rhs_all(const std::vector<std::vector<DiracSpinor>> &X,
                        const std::vector<std::vector<DiracSpinor>> &Y,
                        bool include_Y) const {
  // See hpp. Follows TDHF::dV_rhs term-by-term (both W = Q + P calls, same
  // angular factors), with the radial yk functions hoisted out of the
  // per-task loop; verified against dV_rhs directly in the tests.
  const auto Nc = m_core.size();
  const auto num_points = m_core.front().grid().num_points();

  dVrhs_XY rhs;
  rhs.X.resize(Nc);
  rhs.Y.resize(include_Y ? Nc : 0);
  for (std::size_t ib = 0; ib < Nc; ++ib) {
    for (const auto &X_beta : X[ib]) {
      rhs.X[ib].emplace_back(0, X_beta.kappa(), m_core[ib].grid_sptr());
      if (include_Y) {
        rhs.Y[ib].emplace_back(0, X_beta.kappa(), m_core[ib].grid_sptr());
      }
    }
  }

  // Breit: the dV_Br term is not restructured -- per-task fallback
  if (p_VBr) {
#pragma omp parallel for schedule(dynamic)
    for (std::size_t ib = 0; ib < Nc; ++ib) {
      for (std::size_t be = 0; be < X[ib].size(); ++be) {
        rhs.X[ib][be] = dV_rhs_sets(X[ib][be].kappa(), m_core[ib], false, X, Y);
        if (include_Y) {
          rhs.Y[ib][be] =
            dV_rhs_sets(X[ib][be].kappa(), m_core[ib], true, X, Y);
        }
      }
    }
    return rhs;
  }

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  // Flat channel list, for the shared yk passes
  std::vector<std::pair<std::size_t, std::size_t>> all_channels;
  for (std::size_t jb = 0; jb < Nc; ++jb) {
    for (std::size_t be = 0; be < X[jb].size(); ++be) {
      all_channels.emplace_back(jb, be);
    }
  }

  // Nonzero flags, once per channel function (fresh solves start from 0)
  std::vector<bool> non_zero_X(all_channels.size());
  std::vector<bool> non_zero_Y(all_channels.size());
  for (std::size_t ic = 0; ic < all_channels.size(); ++ic) {
    const auto [jb, be] = all_channels[ic];
    non_zero_X[ic] = X[jb][be].norm2() != 0.0;
    non_zero_Y[ic] = Y[jb][be].norm2() != 0.0;
  }

  // Rank-k channel-core yk pass: y^k(Fb, X_beta) and y^k(Fb, Y_beta) for
  // every channel (left empty while the channel is still zero), combined
  // into the direct screening functions
  //   S_X = sum_beta u_beta [tCk(k, kb, k_beta) y^k(Fb, X_beta)
  //                          + tCk(k, k_beta, kb) y^k(Y_beta, Fb)]
  // (S_Y with X <-> Y), where u(2j) = (-1)^((2j - 1)/2) is the channel
  // factor of the separated dV_rhs phase (-1)^((2j_n - 2j_beta)/2 + k).
  std::vector<std::vector<double>> ykX(all_channels.size()),
    ykY(all_channels.size());
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ic = 0; ic < all_channels.size(); ++ic) {
    const auto [jb, be] = all_channels[ic];
    const auto &Fb = m_core[jb];
    if (non_zero_X[ic]) {
      ykX[ic] = Coulomb::yk_ab(k, Fb, X[jb][be]);
    }
    if (non_zero_Y[ic]) {
      ykY[ic] = Coulomb::yk_ab(k, Fb, Y[jb][be]);
    }
  }
  std::vector<double> S_X(num_points, 0.0), S_Y(num_points, 0.0);
  bool have_S = false;
  for (std::size_t ic = 0; ic < all_channels.size(); ++ic) {
    const auto [jb, be] = all_channels[ic];
    const auto kb = m_core[jb].kappa();
    const auto kbeta = X[jb][be].kappa();
    const auto u_beta = Angular::neg1pow_2(X[jb][be].twoj() - 1);
    const auto c1 = u_beta * Angular::tildeCk_kk(k, kb, kbeta);
    const auto c2 = u_beta * Angular::tildeCk_kk(k, kbeta, kb);
    if (!ykX[ic].empty()) {
      have_S = true;
      for (std::size_t i = 0; i < num_points; ++i) {
        S_X[i] += c1 * ykX[ic][i];
        S_Y[i] += c2 * ykX[ic][i];
      }
    }
    if (!ykY[ic].empty()) {
      have_S = true;
      for (std::size_t i = 0; i < num_points; ++i) {
        S_X[i] += c2 * ykY[ic][i];
        S_Y[i] += c1 * ykY[ic][i];
      }
    }
  }

  // Per core orbital a: the exchange yk block, then all of a's tasks
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ib = 0; ib < Nc; ++ib) {
    const auto &Fa = m_core[ib];
    const auto ka = Fa.kappa();
    const auto tja = Fa.twoj();
    const auto &Fa_f = Fa.f();
    const auto &Fa_g = Fa.g();

    // Exchange block for this a: y^l(chi, Fa) for every channel function
    // chi (s = 0: X_beta, 1: Y_beta) and every l allowed by Ck(l, k_beta,
    // ka), computed up to the partner orbital's extent (it multiplies Fb).
    // Indexed y_exchange[ic][s][(l - lmin)/2].
    std::vector<std::array<std::vector<std::vector<double>>, 2>> y_exchange(
      all_channels.size());
    for (std::size_t ic = 0; ic < all_channels.size(); ++ic) {
      const auto [jb, be] = all_channels[ic];
      const auto max_pt = m_core[jb].max_pt();
      const auto [lmin, lmax] = Angular::kminmax_Ck(X[jb][be].kappa(), ka);
      for (int l = lmin; l <= lmax; l += 2) {
        y_exchange[ic][0].push_back(non_zero_X[ic] ?
                                      Coulomb::yk_ab(l, X[jb][be], Fa, max_pt) :
                                      std::vector<double>{});
        y_exchange[ic][1].push_back(non_zero_Y[ic] ?
                                      Coulomb::yk_ab(l, Y[jb][be], Fa, max_pt) :
                                      std::vector<double>{});
      }
    }

    // One task: rhs = dV_rhs(kappa_n, Fa, conj), assembled from the shared
    // radial functions with dV_rhs's angular factors
    const auto assemble = [&](int kappa_n, bool conj, DiracSpinor *out) {
      const auto tjn = Angular::twoj_k(kappa_n);
      auto &out_f = out->f();
      auto &out_g = out->g();

      // Direct (Q) parts of both W terms: u_n tCk(k, kappa_n, ka) S(r) Fa
      if (have_S) {
        const auto &S = conj ? S_Y : S_X;
        const auto cQ =
          Angular::neg1pow_2(tjn - 1) * Angular::tildeCk_kk(k, kappa_n, ka);
        if (cQ != 0.0) {
          for (auto i = Fa.min_pt(); i < Fa.max_pt(); ++i) {
            out_f[i] += cQ * S[i] * Fa_f[i];
            out_g[i] += cQ * S[i] * Fa_g[i];
          }
        }
      }

      // Exchange (P) parts
      for (std::size_t ic = 0; ic < all_channels.size(); ++ic) {
        const auto [jb, be] = all_channels[ic];
        const auto &Fb = m_core[jb];
        const auto kb = Fb.kappa();
        const auto tjb = Fb.twoj();
        const auto &chi = conj ? Y[jb][be] : X[jb][be];
        const auto kbeta = chi.kappa();
        const auto tjbeta = chi.twoj();
        const auto sQ = Angular::neg1pow_2(tjn - tjbeta + 2 * k);
        const auto lmin_exchange = Angular::kminmax_Ck(kbeta, ka).first;

        // P of W(k, n, Fb, Fa, chi) = sum_l 6j Q^l(n, Fb, chi, Fa):
        // fixed y^l(Fb, Fa), multiplies chi
        if (conj ? non_zero_Y[ic] : non_zero_X[ic]) {
          const auto min_twol =
            std::max(std::abs(tjbeta - tjn), std::abs(tja - tjb));
          const auto max_twol = std::min(tjbeta + tjn, tja + tjb);
          const auto &chi_f = chi.f();
          const auto &chi_g = chi.g();
          for (int tl = min_twol; tl <= max_twol; tl += 2) {
            const auto l = tl / 2;
            if (!Angular::Ck_kk_SR(l, kb, ka) ||
                !Angular::Ck_kk_SR(l, kappa_n, kbeta)) {
              continue;
            }
            const auto sixj = Angular::sixj_2(tja, tjn, 2 * k, tjbeta, tjb, tl);
            if (sixj == 0.0) {
              continue;
            }
            const auto m1tl = Angular::evenQ(l) ? 1 : -1;
            const auto coef = sQ * tkp1 * sixj * m1tl *
                              Angular::tildeCk_kk(l, kappa_n, kbeta) *
                              Angular::tildeCk_kk(l, kb, ka);
            const auto *y = m_ycc->get(l, Fb, Fa);
            assert(y && "Ck-allowed l is always stored in the YkTable");
            for (auto i = chi.min_pt(); i < chi.max_pt(); ++i) {
              out_f[i] += coef * (*y)[i] * chi_f[i];
              out_g[i] += coef * (*y)[i] * chi_g[i];
            }
          }
        }

        // P of W(k, n, eta, Fa, Fb) = sum_l 6j Q^l(n, eta, Fb, Fa):
        // y^l(eta, Fa) from the block, multiplies Fb
        if (conj ? non_zero_X[ic] : non_zero_Y[ic]) {
          const auto s_eta = conj ? 0 : 1;
          const auto min_twol =
            std::max(std::abs(tjb - tjn), std::abs(tja - tjbeta));
          const auto max_twol = std::min(tjb + tjn, tja + tjbeta);
          const auto &Fb_f = Fb.f();
          const auto &Fb_g = Fb.g();
          for (int tl = min_twol; tl <= max_twol; tl += 2) {
            const auto l = tl / 2;
            if (!Angular::Ck_kk_SR(l, kbeta, ka) ||
                !Angular::Ck_kk_SR(l, kappa_n, kb)) {
              continue;
            }
            const auto sixj = Angular::sixj_2(tja, tjn, 2 * k, tjb, tjbeta, tl);
            if (sixj == 0.0) {
              continue;
            }
            const auto m1tl = Angular::evenQ(l) ? 1 : -1;
            const auto coef = sQ * tkp1 * sixj * m1tl *
                              Angular::tildeCk_kk(l, kappa_n, kb) *
                              Angular::tildeCk_kk(l, kbeta, ka);
            const auto &y = y_exchange[ic][std::size_t(s_eta)]
                                      [std::size_t(l - lmin_exchange) / 2];
            for (auto i = Fb.min_pt(); i < Fb.max_pt(); ++i) {
              out_f[i] += coef * y[i] * Fb_f[i];
              out_g[i] += coef * y[i] * Fb_g[i];
            }
          }
        }
      }

      (*out) *= (1.0 / tkp1);
    };

    for (std::size_t be = 0; be < X[ib].size(); ++be) {
      assemble(X[ib][be].kappa(), false, &rhs.X[ib][be]);
      if (include_Y) {
        assemble(X[ib][be].kappa(), true, &rhs.Y[ib][be]);
      }
    }
  }

  return rhs;
}

//==============================================================================
void TDHFcomplex::solve_core(double omega, int max_its, bool print) {
  assert(m_h->rank() == m_rank && "Rank must match in solve_core");
  assert(m_h->parity() == m_pi && "Parity must match in solve_core");
  assert(m_h->imaginaryQ() == m_imag && "Imaginarity must match in solve_core");

  omega = std::abs(omega);

  m_hFcore = form_hFcore(m_h);
  m_hFcore_minus = form_hFcore(m_h_minus);

  prepare(omega);
  assert(m_Xi.size() == m_X.size() && m_Yi.size() == m_Y.size());

  // Starts from the current (complex) corrections: a warm start at the
  // same omega, else the previous omega's solution as the initial guess
  solve_core_anderson(omega, max_its, print);
  m_have_solution = true;

  // K_+ = pi*D consistency of the converged open channels (see KpiD_dev):
  // the one measure that reliably detects an unreliable channel solve
  // (marginal grid resolution, invalid extraction window, bad pair), which
  // can otherwise LOOK converged. Absolute deviation relative to the
  // largest amplitude, so a weak channel with a large relative but
  // negligible absolute error cannot fire the alarm. Stored, never printed.
  m_KpiD = 0.0;
  m_KpiD_lab.clear();
  std::vector<std::pair<const DiracSpinor *, OutgoingChannel>> all_channels;
  double K_scale = 0.0;
  for (const auto &Fb : m_core) {
    for (const auto &channel : outgoing_channels(Fb)) {
      K_scale =
        std::max({K_scale, std::abs(channel.K), M_PI * std::abs(channel.D)});
      all_channels.emplace_back(&Fb, channel);
    }
  }
  for (const auto &[pFb, channel] : all_channels) {
    if (K_scale == 0.0) {
      break;
    }
    const auto dev = std::abs(channel.K - M_PI * channel.D) / K_scale;
    if (dev > m_KpiD) {
      m_KpiD = dev;
      m_KpiD_lab =
        pFb->shortSymbol() + "," + AtomData::kappa_symbol(channel.kappa);
    }
  }
}

//==============================================================================
void TDHFcomplex::solve_core_anderson(double omega, int max_its, bool print) {
  // Anderson/Pulay (DIIS) accelerated driver. One "map application" is a
  // single UNDAMPED tdhf_core_it_complex (all channels, X and Y, real and
  // imaginary): the TDHF fixed-point problem is linear, x = A x + b, so
  // DIIS extrapolation over the map outputs (equivalent to preconditioned
  // GMRES) converges wherever (1 - A) is non-singular -- in particular
  // through the resonance windows where damped iteration diverges. The
  // state = all corrections (real and imaginary parts) and the complex
  // channel amplitudes, which are linear in the state and so extrapolate
  // with the same coefficients.

  const double converge_targ = m_eps;

  // Flatten the full state into a single real vector, and back
  const auto flatten = [this]() {
    std::vector<double> v;
    for (const auto *set : {&m_X, &m_Y, &m_Xi, &m_Yi}) {
      for (const auto &Fs : *set) {
        for (const auto &F : Fs) {
          v.insert(v.end(), F.f().cbegin(), F.f().cend());
          v.insert(v.end(), F.g().cbegin(), F.g().cend());
        }
      }
    }
    for (const auto &orbital_channels : m_channels) {
      for (const auto &channel : orbital_channels) {
        v.push_back(channel.K.real());
        v.push_back(channel.K.imag());
      }
    }
    return v;
  };
  const auto unflatten = [this](const std::vector<double> &v) {
    std::size_t i = 0;
    for (auto *set : {&m_X, &m_Y, &m_Xi, &m_Yi}) {
      for (auto &Fs : *set) {
        for (auto &F : Fs) {
          std::copy_n(v.cbegin() + long(i), F.f().size(), F.f().begin());
          i += F.f().size();
          std::copy_n(v.cbegin() + long(i), F.g().size(), F.g().begin());
          i += F.g().size();
        }
      }
    }
    for (auto &orbital_channels : m_channels) {
      for (auto &channel : orbital_channels) {
        channel.K = {v[i], v[i + 1]};
        i += 2;
      }
    }
    assert(i == v.size());
  };
  const auto dot = [](const std::vector<double> &a,
                      const std::vector<double> &b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
      sum += a[i] * b[i];
    }
    return sum;
  };

  // Anderson history depth: the dominant memory user -- each history entry
  // holds TWO full copies of the state (all real and imaginary spinor
  // components). Depth 4 converges essentially as well as deeper
  // histories here (typical solves take a handful of iterations).
  constexpr std::size_t and_dim = 4;
  constexpr double cond_floor = 1.0e-14;   // drop-oldest threshold on B
  std::vector<std::vector<double>> g_hist; // map outputs g_k = G(x_k)
  std::vector<std::vector<double>> r_hist; // residuals r_k = g_k - x_k

  std::pair<double, std::string> eps{};
  double best_eps{1.0e30};
  int count_worse = 0;
  int it{0};
  qip::LiveMessage status(
    fmt::format("TDHFcomplex {} (w={:.4f}): ", m_h->name(), omega), print);
  for (; it < max_its; it++) {
    auto x = flatten();
    // Undamped map application; eps measures G(x) vs x -- exactly the
    // residual, in the physical (K-weighted) measure
    eps = tdhf_core_it_complex(omega);

    status(fmt::format("{:2d} {:.1e} [{}]", it, eps.first, eps.second));

    if (std::isnan(eps.first))
      break; // broken
    if (eps.first < converge_targ)
      break; // converged (the sets hold the latest map output)

    auto g = flatten();
    auto r = g;
    for (std::size_t i = 0; i < r.size(); ++i) {
      r[i] -= x[i];
    }
    g_hist.push_back(std::move(g));
    r_hist.push_back(std::move(r));
    if (g_hist.size() > and_dim) {
      g_hist.erase(g_hist.begin());
      r_hist.erase(r_hist.begin());
    }

    // First step: nothing to mix yet -- plain fixed-point step (the sets
    // already hold g; from a fresh start this is the first-order result,
    // so max_its == 1 gives the exact first-order correction)
    auto m = g_hist.size();
    if (m == 1)
      continue;

    // Anderson coefficients: minimise || sum_i c_i r_i ||^2 s.t.
    // sum_i c_i = 1, via the bordered system [B 1; 1' 0][c; lam] = [0; 1],
    // B_ij = <r_i|r_j>. If B is ill-conditioned (saturated history), drop
    // the oldest and re-build.
    LinAlg::Vector<double> c;
    while (true) {
      LinAlg::Matrix<double> B(m + 1, m + 1);
      LinAlg::Vector<double> bvec(m + 1);
      double dmax = 0.0, dmin = 1.0e300;
      for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
          B(i, j) = dot(r_hist[i], r_hist[j]);
        }
        B(i, m) = 1.0;
        B(m, i) = 1.0;
        bvec(i) = 0.0;
        dmax = std::max(dmax, B(i, i));
        dmin = std::min(dmin, B(i, i));
      }
      B(m, m) = 0.0;
      bvec(m) = 1.0;
      if (m > 2 && dmin < cond_floor * dmax) {
        g_hist.erase(g_hist.begin());
        r_hist.erase(r_hist.begin());
        --m;
        continue;
      }
      c = LinAlg::solve_Axeqb<double>(B, bvec);
      // Guard a singular/near-singular LU (no exceptions in this codebase):
      // drop the oldest entry and re-build if the coefficients are bad
      bool finite = true;
      for (std::size_t i = 0; i < m; ++i) {
        if (!std::isfinite(c(i))) {
          finite = false;
        }
      }
      if (finite)
        break;
      g_hist.erase(g_hist.begin());
      r_hist.erase(r_hist.begin());
      --m;
      if (m == 0)
        break;
    }
    if (m == 0)
      continue;

    // New iterate: x = sum_i c_i g_i (DIIS extrapolation)
    auto x_new = std::vector<double>(g_hist[0].size(), 0.0);
    for (std::size_t i = 0; i < m; ++i) {
      const auto ci = c(i);
      const auto &gi = g_hist[i];
      for (std::size_t k = 0; k < x_new.size(); ++k) {
        x_new[k] += ci * gi[k];
      }
    }
    unflatten(x_new);

    // Stalled: no meaningful improvement for several iterations -> give up
    // (the result is still used; the warning stars flag its quality)
    if (eps.first < 0.99 * best_eps) {
      best_eps = eps.first;
      count_worse = 0;
    } else if (++count_worse > 5) {
      break;
    }
  }

  // Soft visual warning, relative to the convergence target
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
std::pair<double, std::string> TDHFcomplex::tdhf_core_it_complex(double omega) {
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

  const auto dV_re = dV_rhs_all(m_X, m_Y, true);
  const auto dV_im = dV_rhs_all(m_Xi, m_Yi, true);

  // Previous iteration's amplitudes (for the convergence measure); the
  // channel solves write the new K into m_channels in place
  std::vector<std::vector<std::complex<double>>> K_old(m_channels.size());
  for (std::size_t ib = 0; ib < m_channels.size(); ++ib) {
    for (const auto &channel : m_channels[ib]) {
      K_old[ib].push_back(channel.K);
    }
  }

  // Flatten all (core orbital x channel x X/Y) solves into one task list
  // (load balance). Only the X/+ tasks carry the continuum channel cache
  // (Y/- partners are always bound).
  struct MsTask {
    DiracSpinor *dF_re;        // targets (in Xs/Ys and Xis/Yis)
    DiracSpinor *dF_im;        //
    ContinuumChannel *channel; // continuum cache (nullptr for Y)
    const DiracSpinor *Fb;     // core orbital
    const DiracSpinor *hFb;    // external source h|Fb> (real part only)
    const DiracSpinor *dV_re;  // this task's [dV phi]_beta, real set
    const DiracSpinor *dV_im;  // and imaginary set
    dPsiType type;
  };
  std::vector<MsTask> tasks;
  for (std::size_t ib = 0; ib < m_core.size(); ++ib) {
    for (std::size_t be = 0; be < Xs[ib].size(); ++be) {
      tasks.push_back({&Xs[ib][be], &Xis[ib][be], &m_channels[ib][be],
                       &m_core[ib], &m_hFcore[ib][be], &dV_re.X[ib][be],
                       &dV_im.X[ib][be], dPsiType::X});
      tasks.push_back({&Ys[ib][be], &Yis[ib][be], nullptr, &m_core[ib],
                       &m_hFcore_minus[ib][be], &dV_re.Y[ib][be],
                       &dV_im.Y[ib][be], dPsiType::Y});
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (std::size_t it = 0; it < tasks.size(); ++it) {
    const auto &t = tasks[it];
    // Every X/+ channel of an ionised orbital is open (en_b + omega > 0)
    const bool ionised = t.type == dPsiType::X && t.Fb->en() + omega > 0.0;
    if (ionised && (m_suppress_open || !t.channel->open)) {
      // Open channel excluded from the response (explicit zero): by the
      // suppress_open option, or the continuum solve returned zero
      // (prepare_channels backstop)
      t.dF_re->f().assign(t.dF_re->f().size(), 0.0);
      t.dF_re->g().assign(t.dF_re->g().size(), 0.0);
      t.dF_im->f().assign(t.dF_im->f().size(), 0.0);
      t.dF_im->g().assign(t.dF_im->g().size(), 0.0);
      t.channel->K = 0.0;
    } else if (ionised) {
      solve_channel_outgoing(t.dF_re, t.dF_im, t.channel, *t.Fb, *t.hFb,
                             *t.dV_re, *t.dV_im, omega, eps_ms);
    } else {
      // Bound channel: two independent real solves, the external field in
      // the real part only
      solve_channel_bound(t.dF_re, *t.Fb, t.hFb, *t.dV_re, omega, t.type,
                          eps_ms);
      solve_channel_bound(t.dF_im, *t.Fb, nullptr, *t.dV_im, omega, t.type,
                          eps_ms);
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
void TDHFcomplex::solve_channel_outgoing(
  DiracSpinor *X_re, DiracSpinor *X_im, ContinuumChannel *channel,
  const DiracSpinor &Fb, const DiracSpinor &hFb, const DiracSpinor &dV_re,
  const DiracSpinor &dV_im, double omega, double eps_ms) const {
  // Open X/+ channel of ionised orbital Fb, en_+ = en_b + omega > 0, with
  // the outgoing-wave boundary condition (see class description). The
  // one-electron self-interaction V^b_0 = y^0_bb + X_b of the hole is
  // subtracted from the static Hamiltonian (local potential vl - y^0_bb:
  // residual-ion tail; exchange part via the solver's Fhole option) AND
  // compensated (lagged) in the source, an exact rearrangement whose
  // lagged source is short-ranged.
  using namespace qip::overloads;
  assert(channel != nullptr && channel->open);
  const int kappa_beta = X_re->kappa();

  // Physical sources: [(t + dV)phi_b]_beta. The external field is real in
  // the code's convention, so it drives the real part only.
  auto rhs_re = dV_re + hFb;
  auto rhs_im = dV_im;
  // Diagonal channel (even-parity operators only): subtract the
  // norm-conservation (Lagrange multiplier) term de*phi_b, de = <b|(t +
  // dV)phi_b>. Applies to the physical source ONLY: the V^b_0 compensation
  // below is part of the operator rearrangement, and projecting it too
  // would over-subtract.
  if (kappa_beta == Fb.kappa()) {
    rhs_re -= (Fb * rhs_re) * Fb;
    rhs_im -= (Fb * rhs_im) * Fb;
  }
  // The lagged compensation from the FULL previous complex iterate (X_re,
  // X_im still hold it)
  rhs_re += hole_compensation(Fb, *X_re);
  rhs_im += hole_compensation(Fb, *X_im);

  const auto y0bb = Coulomb::yk_ab(0, Fb, Fb);
  const auto vl_c = p_hf->vlocal(Angular::l_k(kappa_beta)) - y0bb;

  // Warm start of the inner (exchange) iterations from the STANDING parts
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
void TDHFcomplex::solve_channel_bound(DiracSpinor *dF_beta,
                                      const DiracSpinor &Fb,
                                      const DiracSpinor *hFb,
                                      const DiracSpinor &dV_src, double omega,
                                      dPsiType XorY, double eps_ms) const {
  // Bound channel (closed orbital, or the Y/- partner of an ionised
  // orbital): plain V^N bound solve, as bound TDHF (solveMixedState_cntm:
  // Anderson mixing, robust against the spurious-preconditioner divergence
  // at high omega). The V^{N-1} rearrangement is applied per channel and
  // ONLY to the open channels; the bound Y/- partners of an ionised
  // orbital must not run the lagged compensation through their
  // near-resonant conditioning (transient mismatches between
  // near-degenerate X/Y partner channels are amplified by ~1/(e0 - em)).
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
DiracSpinor TDHFcomplex::hole_compensation(const DiracSpinor &Fa,
                                           const DiracSpinor &chi) const {
  // + V^a_0 chi = (y^0_aa + X_a) chi, the lagged compensation for the
  // one-electron self-interaction moved into the static V^{N-1} Hamiltonian
  using namespace qip::overloads;
  const auto y0aa = Coulomb::yk_ab(0, Fa, Fa);
  return (y0aa * chi) + HF::vexFa_1el(chi, Fa);
}

//==============================================================================
std::pair<double, std::string> TDHFcomplex::eps_complex(
  const std::vector<std::vector<DiracSpinor>> &Xs_re,
  const std::vector<std::vector<DiracSpinor>> &Xs_im,
  const std::vector<std::vector<std::complex<double>>> &K_old) const {
  double DdF2 = 0.0;
  double dF2 = 0.0;
  double worst = 0.0;
  double worst_K = 0.0;
  std::string worst_lab;

  // Floor for the open-channel measure: the largest |K_+|^2 (a channel
  // through a zero must not gate the SCF)
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
std::vector<std::pair<std::size_t, std::size_t>>
TDHFcomplex::list_open_channels(std::vector<Channel> *channels) const {
  std::vector<std::pair<std::size_t, std::size_t>> index_list;
  for (auto ib = 0ul; ib < m_channels.size(); ib++) {
    for (auto be = 0ul; be < m_channels[ib].size(); be++) {
      if (m_channels[ib][be].open) {
        index_list.emplace_back(ib, be);
        channels->push_back(
          {ib, m_X[ib][be].kappa(), m_core[ib].en() + m_omega});
      }
    }
  }
  return index_list;
}

//==============================================================================
std::pair<std::size_t, std::size_t>
TDHFcomplex::channel_index(const DiracSpinor &Fb, int kappa) const {
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
std::vector<std::complex<double>> TDHFcomplex::A_phys() const {
  std::vector<std::complex<double>> A;
  if (m_channels.empty() || m_omega < 0.0) {
    return A;
  }
  std::vector<Channel> channels;
  const auto index_list = list_open_channels(&channels);
  for (const auto &[ib, be] : index_list) {
    A.push_back(std::conj(m_channels[ib][be].K) / M_PI);
  }
  return A;
}

//==============================================================================
std::complex<double> TDHFcomplex::A_phys(const DiracSpinor &Fa,
                                         int kappa_e) const {
  if (m_channels.empty() || m_omega < 0.0) {
    return 0.0;
  }
  const auto [ib, be] = channel_index(Fa, kappa_e);
  if (!m_channels[ib][be].open) {
    return 0.0;
  }
  return std::conj(m_channels[ib][be].K) / M_PI;
}

//==============================================================================
double TDHFcomplex::D_phys(const DiracSpinor &Fa, int kappa_e) const {
  return std::abs(A_phys(Fa, kappa_e));
}

//==============================================================================
std::vector<TDHFcomplex::OutgoingChannel>
TDHFcomplex::outgoing_channels(const DiracSpinor &Fa) const {
  // Complex K_+ and D = <Freg|S> from the TOTAL effective source of the
  // converged solve,
  //   S = (t + dV')phi_a + (V^nl - X_a - U_KS) chi ,
  // real and imaginary parts separately (with the diagonal de projection,
  // mirroring the channel solve). The radial solve satisfies
  // K = pi <Freg|S> identically (surface-term identity for the eigenpair
  // of the local operator), so the residual measures only the solve /
  // extraction accuracy, not the exchange.
  using namespace qip::overloads;
  std::vector<OutgoingChannel> out;
  if (m_channels.empty() || m_hFcore.empty() || m_omega < 0.0) {
    return out;
  }
  const auto ib = static_cast<std::size_t>(
    std::find(m_core.cbegin(), m_core.cend(), Fa) - m_core.cbegin());
  assert(ib < m_channels.size());
  const auto Ux = HF::vex_KS(m_core);

  const auto total_source = [&](int kappa, const DiracSpinor &S_phys,
                                const DiracSpinor &chi) {
    auto S_t = S_phys;
    if (kappa == Fa.kappa()) {
      S_t -= (Fa * S_t) * Fa;
    }
    return S_t + hole_compensation(Fa, chi) + HF::vexFa(chi, m_core) -
           HF::vexFa_1el(chi, Fa) - (Ux * chi);
  };

  for (std::size_t be = 0; be < m_channels[ib].size(); ++be) {
    const auto &channel = m_channels[ib][be];
    if (!channel.open) {
      continue;
    }
    const auto kappa = channel.Freg.kappa();
    const auto S_re = total_source(
      kappa, m_hFcore[ib][be] + dV_rhs(kappa, Fa, false), m_X[ib][be]);
    const auto S_im = total_source(
      kappa, dV_rhs_sets(kappa, Fa, false, m_Xi, m_Yi), m_Xi[ib][be]);
    const std::complex<double> D{channel.Freg * S_re, channel.Freg * S_im};
    out.push_back({kappa, Fa.en() + m_omega, channel.K, D});
  }
  return out;
}

//==============================================================================
std::complex<double> TDHFcomplex::dV_complex(const DiracSpinor &Fa,
                                             const DiracSpinor &Fb) const {
  // As TDHF::dV, on the real and imaginary sets. A continuum bra (en_a > 0)
  // gets the V^{N-1} treatment of the ionised electron: the source carries
  // the hole compensation V^b_0 chi, whose +y^0_bb*chi term cancels
  // pointwise the 1/r tail hidden in the b-diagonal part of dV*phi_b, so
  // the continuum-continuum overlap is box-independent; Fa must then be
  // the V^{N-1} continuum state of hole Fb (incl. the exchange part).
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
double TDHFcomplex::dV(const DiracSpinor &Fa, const DiracSpinor &Fb) const {
  assert(Fa.en() <= 0.0 && Fb.en() <= 0.0 &&
         "TDHFcomplex::dV is real: for a continuum state use dV_complex");
  return dV_complex(Fa, Fb).real();
}

//==============================================================================
double TDHFcomplex::imaginary_norm2() const {
  double sum = 0.0;
  for (const auto *set : {&m_Xi, &m_Yi}) {
    for (const auto &Fs : *set) {
      for (const auto &F : Fs) {
        sum += F.norm2();
      }
    }
  }
  return sum;
}

} // namespace ExternalField
