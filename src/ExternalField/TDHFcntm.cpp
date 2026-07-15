#include "TDHFcntm.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "DiracODE/ContinuumState.hpp"
#include "DiracOperator/TensorOperator.hpp"
#include "ExternalField/MixedStates.hpp"
#include "HF/HartreeFock.hpp"
#include "LinAlg/Matrix.hpp"
#include "LinAlg/Solvers.hpp"
#include "LinAlg/Vector.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "fmt/format.hpp"
#include "qip/Maths.hpp"
#include "qip/Widgets.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

namespace ExternalField {

//==============================================================================
TDHFcntm::TDHFcntm(const DiracOperator::TensorOperator *const h_plus,
                   const HF::HartreeFock *const hf,
                   const DiracOperator::TensorOperator *const h_minus)
  : TDHF(h_plus, hf, h_minus) {}

//==============================================================================
void TDHFcntm::clear() {
  TDHF::clear();
  m_ch.clear();
  m_omega = -1.0;
  m_resc.reset();
}

//==============================================================================
void TDHFcntm::prepare_channels(double omega, bool print) {
  // (Re)build the per-channel continuum caches for this omega: openness, and
  // for each open channel the homogeneous pair Freg/Firr at en_+ = en_b +
  // omega, built ONCE in the fixed conditioning potential
  //   v = vlocal(l) - y^0_bb + U_KS
  // -- the same potential solveContinuumMixedState() uses, so it can reuse
  // them across all the outer TDHF iterations (this both saves the two ODE
  // solves per channel per iteration, and removes the per-rebuild F_reg
  // renormalisation noise, which otherwise sets the TDHF convergence floor).
  //
  // An open channel is kept only if it is RESOLVABLE on the radial box: the
  // standing-wave boundary condition and the energy normalisation both need
  // the outer (constant-amplitude) envelope, i.e. a few asymptotic
  // wavelengths inside r_max. Just above a shell's threshold (en_+ -> 0+)
  // that fails: the continuum solve returns noise, which feeds noise into
  // dV and the TDHF cannot converge (the map is no longer reproducibly
  // linear). Such channels are excluded from the response (zeroed, as
  // suppress_open does) -- an approximation ONLY for omega within
  // ~ (3*2*pi/r_max)^2/2 of that shell's threshold; enlarge the box to
  // resolve them.
  using namespace qip::overloads;

  m_ch.clear();
  m_ch.resize(m_X.size());
  m_omega = omega;

  const auto Ux = HF::vex_KS(m_core);

  const auto r_max = m_core.front().grid().r().back();
  // at least ~one asymptotic wavelength must fit in the box (nonrel k is
  // fine for this estimate)
  const auto n_waves_min = 1.0;
  const auto ec_min = 0.5 * qip::pow<2>(n_waves_min * 2.0 * M_PI / r_max);

  for (auto ib = 0ul; ib < m_X.size(); ib++) {
    const auto &Fb = m_core[ib];
    const auto en_plus = Fb.en() + omega;
    const bool open = en_plus > 0.0;
    const bool resolvable = en_plus > ec_min;

    if (open && !resolvable && print) {
      fmt::print(
        "\nWarning: TDHFcntm: {} channels open but unresolvable on this box "
        "(en_+ = {:.1e} < {:.1e} au): excluded from the response. Enlarge "
        "r_max to include them.\n",
        Fb.shortSymbol(), en_plus, ec_min);
    }

    const auto y0bb = open ? Coulomb::yk_ab(0, Fb, Fb) : std::vector<double>{};

    for (const auto &X_beta : m_X[ib]) {
      const auto kappa = X_beta.kappa();
      DiracSpinor Freg{0, kappa, Fb.grid_sptr()};
      DiracSpinor Firr{0, kappa, Fb.grid_sptr()};
      if (open && resolvable) {
        const auto v = p_hf->vlocal(Angular::l_k(kappa)) - y0bb + Ux;
        DiracODE::solveContinuum(Freg, en_plus, v, m_alpha);
        DiracODE::solveContinuumIrregular(Firr, Freg, en_plus, v, m_alpha);
      }
      // nb: "open" in the cache means open AND resolvable: unresolvable
      // channels are dispatched exactly like suppress_open (zeroed)
      m_ch[ib].emplace_back(open && resolvable, std::move(Freg),
                            std::move(Firr));
    }
  }
}

//==============================================================================
void TDHFcntm::build_ycc() {
  // See hpp: fixed core-core y^l(Fi, Fj) table for the exchange part of
  // dV_rhs_all(). Full grid: the y^l tail region multiplies the (continuum)
  // channel spinors.
  if (m_ycc) {
    return;
  }
  const auto Nc = m_core.size();
  auto ycc = std::make_shared<YccTable>(Nc * (Nc + 1) / 2);
  std::vector<std::pair<std::size_t, std::size_t>> prs;
  for (std::size_t i = 0; i < Nc; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      prs.emplace_back(i, j);
    }
  }
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ij = 0; ij < prs.size(); ++ij) {
    const auto &Fi = m_core[prs[ij].first];
    const auto &Fj = m_core[prs[ij].second];
    const auto [lmin, lmax] = Angular::kminmax_Ck(Fi.kappa(), Fj.kappa());
    auto &yij = (*ycc)[ij];
    for (int l = lmin; l <= lmax; l += 2) {
      yij.push_back(Coulomb::yk_ab(l, Fi, Fj));
    }
  }
  m_ycc = ycc;
}

//------------------------------------------------------------------------------
const std::vector<double> &TDHFcntm::ycc_get(std::size_t i, std::size_t j,
                                             int l) const {
  // y^l(core_i, core_j); the caller guarantees l is Ck-allowed for (i, j)
  if (j > i) {
    std::swap(i, j);
  }
  const auto lmin =
    Angular::kminmax_Ck(m_core[i].kappa(), m_core[j].kappa()).first;
  return (*m_ycc)[i * (i + 1) / 2 + j][std::size_t(l - lmin) / 2];
}

//==============================================================================
TDHFcntm::DVRhsAll TDHFcntm::dV_rhs_all(bool include_Y) const {
  // See hpp. Follows TDHF::dV_rhs term-by-term (both W = Q + P calls, same
  // angular factors), with the radial yk functions hoisted out of the
  // per-task loop; verified against dV_rhs directly in the [cntmrpa] tests.
  assert(m_ycc && "solve_core() builds the core-core yk table first");

  const auto Nc = m_core.size();
  const auto num_points = m_core.front().grid().num_points();

  DVRhsAll rhs;
  rhs.X.resize(Nc);
  rhs.Y.resize(include_Y ? Nc : 0);
  for (std::size_t ib = 0; ib < Nc; ++ib) {
    for (const auto &X_beta : m_X[ib]) {
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
      for (std::size_t be = 0; be < m_X[ib].size(); ++be) {
        rhs.X[ib][be] = dV_rhs(m_X[ib][be].kappa(), m_core[ib], false);
        if (include_Y) {
          rhs.Y[ib][be] = dV_rhs(m_X[ib][be].kappa(), m_core[ib], true);
        }
      }
    }
    return rhs;
  }

  const auto k = m_h->rank();
  const auto tkp1 = double(2 * k + 1);

  // Flat channel list, for the shared yk passes
  std::vector<std::pair<std::size_t, std::size_t>> chans;
  for (std::size_t jb = 0; jb < Nc; ++jb) {
    for (std::size_t be = 0; be < m_X[jb].size(); ++be) {
      chans.emplace_back(jb, be);
    }
  }

  // Nonzero flags, once per channel function (fresh solves start from 0)
  std::vector<char> nzX(chans.size()), nzY(chans.size());
  for (std::size_t ic = 0; ic < chans.size(); ++ic) {
    const auto [jb, be] = chans[ic];
    nzX[ic] = m_X[jb][be].norm2() != 0.0;
    nzY[ic] = m_Y[jb][be].norm2() != 0.0;
  }

  // Rank-k channel-core yk pass: y^k(Fb, X_beta) and y^k(Fb, Y_beta) for
  // every channel (left empty while the channel is still zero), combined
  // into the direct screening functions
  //   S_X = sum_beta u_beta [tCk(k, kb, k_beta) y^k(Fb, X_beta)
  //                          + tCk(k, k_beta, kb) y^k(Y_beta, Fb)]
  // (S_Y with X <-> Y), where u(2j) = (-1)^((2j - 1)/2) is the channel
  // factor of the separated dV_rhs phase (-1)^((2j_n - 2j_beta)/2 + k).
  std::vector<std::vector<double>> ykX(chans.size()), ykY(chans.size());
#pragma omp parallel for schedule(dynamic)
  for (std::size_t ic = 0; ic < chans.size(); ++ic) {
    const auto [jb, be] = chans[ic];
    const auto &Fb = m_core[jb];
    if (nzX[ic]) {
      ykX[ic] = Coulomb::yk_ab(k, Fb, m_X[jb][be]);
    }
    if (nzY[ic]) {
      ykY[ic] = Coulomb::yk_ab(k, Fb, m_Y[jb][be]);
    }
  }
  std::vector<double> S_X(num_points, 0.0), S_Y(num_points, 0.0);
  bool have_S = false;
  for (std::size_t ic = 0; ic < chans.size(); ++ic) {
    const auto [jb, be] = chans[ic];
    const auto kb = m_core[jb].kappa();
    const auto kbeta = m_X[jb][be].kappa();
    const auto u_beta = Angular::neg1pow_2(m_X[jb][be].twoj() - 1);
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
    const auto &faf = Fa.f();
    const auto &fag = Fa.g();

    // Exchange block for this a: y^l(chi, Fa) for every channel function
    // chi (s = 0: X_beta, 1: Y_beta) and every l allowed by Ck(l, k_beta,
    // ka), computed up to the partner orbital's extent (it multiplies Fb).
    // Indexed yblk[ic][s][(l - lmin)/2].
    std::vector<std::array<std::vector<std::vector<double>>, 2>> yblk(
      chans.size());
    for (std::size_t ic = 0; ic < chans.size(); ++ic) {
      const auto [jb, be] = chans[ic];
      const auto maxi = m_core[jb].max_pt();
      const auto [lmin, lmax] = Angular::kminmax_Ck(m_X[jb][be].kappa(), ka);
      for (int l = lmin; l <= lmax; l += 2) {
        yblk[ic][0].push_back(nzX[ic] ?
                                Coulomb::yk_ab(l, m_X[jb][be], Fa, maxi) :
                                std::vector<double>{});
        yblk[ic][1].push_back(nzY[ic] ?
                                Coulomb::yk_ab(l, m_Y[jb][be], Fa, maxi) :
                                std::vector<double>{});
      }
    }

    // One task: rhs = dV_rhs(kappa_n, Fa, conj), assembled from the shared
    // radial functions with dV_rhs's angular factors
    const auto assemble = [&](int kappa_n, bool conj, DiracSpinor *out) {
      const auto tjn = Angular::twoj_k(kappa_n);
      auto &of = out->f();
      auto &og = out->g();

      // Direct (Q) parts of both W terms: u_n tCk(k, kappa_n, ka) S(r) Fa
      if (have_S) {
        const auto &S = conj ? S_Y : S_X;
        const auto cQ =
          Angular::neg1pow_2(tjn - 1) * Angular::tildeCk_kk(k, kappa_n, ka);
        if (cQ != 0.0) {
          for (auto i = Fa.min_pt(); i < Fa.max_pt(); ++i) {
            of[i] += cQ * S[i] * faf[i];
            og[i] += cQ * S[i] * fag[i];
          }
        }
      }

      // Exchange (P) parts
      for (std::size_t ic = 0; ic < chans.size(); ++ic) {
        const auto [jb, be] = chans[ic];
        const auto &Fb = m_core[jb];
        const auto kb = Fb.kappa();
        const auto tjb = Fb.twoj();
        const auto &chi = conj ? m_Y[jb][be] : m_X[jb][be];
        const auto kbeta = chi.kappa();
        const auto tjbeta = chi.twoj();
        const auto sQ = Angular::neg1pow_2(tjn - tjbeta + 2 * k);
        const auto lmin_blk = Angular::kminmax_Ck(kbeta, ka).first;

        // P of W(k, n, Fb, Fa, chi) = sum_l 6j Q^l(n, Fb, chi, Fa):
        // fixed y^l(Fb, Fa), multiplies chi
        if (conj ? nzY[ic] : nzX[ic]) {
          const auto min_twol =
            std::max(std::abs(tjbeta - tjn), std::abs(tja - tjb));
          const auto max_twol = std::min(tjbeta + tjn, tja + tjb);
          const auto &cf = chi.f();
          const auto &cg = chi.g();
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
            const auto &y = ycc_get(jb, ib, l);
            for (auto i = chi.min_pt(); i < chi.max_pt(); ++i) {
              of[i] += coef * y[i] * cf[i];
              og[i] += coef * y[i] * cg[i];
            }
          }
        }

        // P of W(k, n, eta, Fa, Fb) = sum_l 6j Q^l(n, eta, Fb, Fa):
        // y^l(eta, Fa) from the block, multiplies Fb
        if (conj ? nzX[ic] : nzY[ic]) {
          const auto s_eta = conj ? 0 : 1;
          const auto min_twol =
            std::max(std::abs(tjb - tjn), std::abs(tja - tjbeta));
          const auto max_twol = std::min(tjb + tjn, tja + tjbeta);
          const auto &bf = Fb.f();
          const auto &bg = Fb.g();
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
            const auto &y =
              yblk[ic][std::size_t(s_eta)][std::size_t(l - lmin_blk) / 2];
            for (auto i = Fb.min_pt(); i < Fb.max_pt(); ++i) {
              of[i] += coef * y[i] * bf[i];
              og[i] += coef * y[i] * bg[i];
            }
          }
        }
      }

      (*out) *= (1.0 / tkp1);
    };

    for (std::size_t be = 0; be < m_X[ib].size(); ++be) {
      assemble(m_X[ib][be].kappa(), false, &rhs.X[ib][be]);
      if (include_Y) {
        assemble(m_X[ib][be].kappa(), true, &rhs.Y[ib][be]);
      }
    }
  }

  return rhs;
}

//==============================================================================
void TDHFcntm::solve_core(double omega, int max_its, bool print) {
  // Shared set-up, then dispatch to the Anderson (default) or plain damped
  // self-consistency driver; both use the continuum-aware channel dispatch
  // (tdhf_core_it_cntm) for each iteration.

  assert(m_h->rank() == m_rank && "Rank must match in solve_core");
  assert(m_h->parity() == m_pi && "Parity must match in solve_core");
  assert(m_h->imaginaryQ() == m_imag && "Imaginarity must match in solve_core");

  omega = std::abs(omega);

  m_hFcore = form_hFcore(m_h);
  m_hFcore_minus = form_hFcore(m_h_minus);

  // Fixed core-core yk table for dV_rhs_all (first call only)
  build_ycc();

  // New solve: the cached rescattering (D_phys) is stale
  m_resc.reset();
  m_max_its = max_its;

  // Warm start: continuing a previous solve at the SAME omega (e.g. after a
  // first-order run). The first damped iteration must then be damped: an
  // undamped step from an already-large near-resonant X can diverge. From a
  // fresh start (X = 0), the first iteration is undamped as usual (it just
  // builds the first-order correction).
  const bool warm_start = omega == m_omega && !m_ch.empty();

  // (Re)build the continuum channel caches when omega changes
  if (omega != m_omega || m_ch.empty()) {
    prepare_channels(omega, print);
  }

  if (m_anderson) {
    solve_core_anderson(omega, max_its, print);
  } else {
    solve_core_damped(omega, max_its, print, warm_start);
  }
}

//==============================================================================
void TDHFcntm::solve_core_damped(double omega, int max_its, bool print,
                                 bool warm_start) {
  // Damped fixed-point driver, as TDHF::solve_core (with Johnson staging).
  // nb: diverges wherever the Picard multiplier exceeds 1 (autoionising
  // resonances, occupied-occupied near-degeneracies) -- see set_anderson.

  const double converge_targ = m_eps;
  const auto eta_damp = m_eta;

  // Staged iteration (Johnson): it 0 solves both X and Y (dV = 0, so
  // max_its == 1 still gives the exact first-order correction); in stage A
  // the Y/- corrections are then held frozen at first order while X is
  // iterated; stage B (full X+Y) begins once X has roughly converged, a
  // third of the iteration budget is used, or stage A stalls.
  // Only applies in the photoionisation regime (open channels present): for
  // an all-bound core the X-only map can diverge where the full X+Y damped
  // iteration is fine, and there is nothing to gain from staging.
  const bool any_open =
    std::any_of(m_ch.cbegin(), m_ch.cend(), [](const auto &chs) {
      return std::any_of(chs.cbegin(), chs.cend(),
                         [](const auto &ch) { return ch.open; });
    });
  bool stageA = m_staged_Y && max_its > 2 && any_open && !m_suppress_open;

  std::pair<double, std::string> eps{};
  double best_eps{1.0e30};
  int count_worse = 0;
  int it{0};
  qip::LiveMessage status(
    fmt::format("TDHFcntm {} (w={:.4f}): ", m_h->name(), omega), print);
  for (; it < max_its; it++) {
    const auto eta = (it == 0 && !warm_start) ? 0.0 : eta_damp;
    const bool include_Y = it == 0 || !stageA;
    eps = tdhf_core_it_cntm(omega, eta, include_Y);

    status(fmt::format("{:2d} {:.1e} [{}]", it, eps.first, eps.second));

    if (std::isnan(eps.first))
      break; // broken

    if (stageA) {
      // Track stage-A progress; on stall, release Y rather than give up.
      bool stalled = false;
      if (eps.first < 0.99 * best_eps) {
        best_eps = eps.first;
        count_worse = 0;
      } else if (++count_worse > 3) {
        stalled = true;
      }
      // Release Y (stage B); reset the stall detection for the new stage.
      // (No convergence exit while Y is frozen.)
      if (eps.first < std::sqrt(converge_targ) || it + 1 >= max_its / 3 ||
          stalled) {
        stageA = false;
        best_eps = 1.0e30;
        count_worse = 0;
      }
      continue;
    }

    if (eps.first < converge_targ)
      break; // converged

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
void TDHFcntm::solve_core_anderson(double omega, int max_its, bool print) {
  // Anderson/Pulay (DIIS) accelerated driver. One "map application" is a
  // single UNDAMPED tdhf_core_it_cntm (all channels, X and Y): the TDHF
  // fixed-point problem is linear, x = A x + b, so DIIS extrapolation over
  // the map outputs (equivalent to preconditioned GMRES) converges wherever
  // (1 - A) is non-singular -- in particular through the resonance windows
  // where the damped driver diverges. Mirrors the Anderson implementation
  // in solveMixedState_cntm(), with the state = all (X, Y) corrections
  // (and the channel K amplitudes, which are linear in the state, so they
  // extrapolate with the same coefficients).

  const double converge_targ = m_eps;

  // Flatten the full state (all X, all Y spinor components, all K) into a
  // single vector, and back.
  const auto flatten = [this]() {
    std::vector<double> v;
    for (const auto *set : {&m_X, &m_Y}) {
      for (const auto &Fs : *set) {
        for (const auto &F : Fs) {
          v.insert(v.end(), F.f().cbegin(), F.f().cend());
          v.insert(v.end(), F.g().cbegin(), F.g().cend());
        }
      }
    }
    for (const auto &chs : m_ch) {
      for (const auto &ch : chs) {
        v.push_back(ch.K);
      }
    }
    return v;
  };
  const auto unflatten = [this](const std::vector<double> &v) {
    std::size_t i = 0;
    for (auto *set : {&m_X, &m_Y}) {
      for (auto &Fs : *set) {
        for (auto &F : Fs) {
          std::copy_n(v.cbegin() + long(i), F.f().size(), F.f().begin());
          i += F.f().size();
          std::copy_n(v.cbegin() + long(i), F.g().size(), F.g().begin());
          i += F.g().size();
        }
      }
    }
    for (auto &chs : m_ch) {
      for (auto &ch : chs) {
        ch.K = v[i];
        ++i;
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
  // holds TWO full copies of the state (all X, Y spinor components), i.e.
  // 2 * and_dim * n_channels * 2 * num_points doubles in total. Depth 4
  // converges essentially as well as deeper histories here (typical solves
  // take only a handful of iterations).
  constexpr std::size_t and_dim = 4;
  constexpr double cond_floor = 1.0e-14;   // drop-oldest threshold on B
  std::vector<std::vector<double>> g_hist; // map outputs g_k = G(x_k)
  std::vector<std::vector<double>> r_hist; // residuals r_k = g_k - x_k

  std::pair<double, std::string> eps{};
  double best_eps{1.0e30};
  int count_worse = 0;
  int it{0};
  // Label seeded (homogeneous) solves by their seed channel
  const auto seed_lab =
    m_seed ? fmt::format(" seed {},{}", m_core[m_seed->ib].shortSymbol(),
                         m_X[m_seed->ib][m_seed->be].shortSymbol()) :
             std::string{};
  qip::LiveMessage status(
    fmt::format("TDHFcntm {} (w={:.4f}) {}: ", m_h->name(), omega, seed_lab),
    print);
  for (; it < max_its; it++) {
    auto x = flatten();
    // Undamped map application; eps (from eps_cntm) measures G(x) vs x --
    // exactly the residual, in the physical (K-weighted) measure.
    eps = tdhf_core_it_cntm(omega, 0.0, true);

    status(fmt::format("{:2d} {:.1e} [{}]", it, eps.first, eps.second));

    if (std::isnan(eps.first))
      break; // broken
    if (eps.first < converge_targ)
      break; // converged (m_X/m_Y hold the latest map output)

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

    // First step: nothing to mix yet -- plain fixed-point step (m_X/m_Y
    // already hold g; from a fresh start this is the first-order result,
    // so max_its == 1 gives the exact first-order correction).
    auto m = g_hist.size();
    if (m == 1)
      continue;

    // Anderson coefficients: minimise || sum_i c_i r_i ||^2 s.t.
    // sum_i c_i = 1, via the bordered system [B 1; 1' 0][c; lam] = [0; 1],
    // B_ij = <r_i|r_j>. If B is ill-conditioned (saturated history), drop
    // the oldest and re-build. (As solveMixedState_cntm().)
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
      // drop the oldest entry and re-build if the coefficients are bad.
      bool finite = true;
      for (std::size_t i = 0; i < m; ++i) {
        if (!std::isfinite(c(i)))
          finite = false;
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

    // New iterate: x = sum_i c_i g_i (DIIS extrapolation).
    auto x_new = std::vector<double>(g_hist[0].size(), 0.0);
    for (std::size_t i = 0; i < m; ++i) {
      const auto ci = c(i);
      const auto &gi = g_hist[i];
      for (std::size_t k = 0; k < x_new.size(); ++k) {
        x_new[k] += ci * gi[k];
      }
    }
    unflatten(x_new);

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
std::pair<double, std::string>
TDHFcntm::tdhf_core_it_cntm(double omega, double eta_damp, bool include_Y) {
  // Single iteration; mirrors TDHF::tdhf_core_it, but each channel solve
  // goes through the continuum dispatch (solve_channel_cntm). The bound
  // static-Y shortcut is dropped: photoionisation implies omega > 0.
  using namespace qip::overloads;
  assert(omega >= 0.0 && "solve_core() passes omega = |omega|");

  auto Xs = m_X;
  auto Ys = m_Y;

  const auto eps_ms = 1.0e-12;

  // Batched dV sources for every task, one shared-yk pass (reads m_X/m_Y,
  // which are fixed for this iteration: the solves write into Xs/Ys).
  const auto dv = dV_rhs_all(include_Y);

  // Previous iteration's K amplitudes (for the convergence measure); the
  // channel solves write the new K into m_ch in place.
  std::vector<std::vector<double>> K_old(m_ch.size());
  for (auto ib = 0ul; ib < m_ch.size(); ib++) {
    for (const auto &ch : m_ch[ib]) {
      K_old[ib].push_back(ch.K);
    }
  }

  // Flatten all (core orbital x channel x X/Y) solves into one task list
  // (load balance; see TDHF::tdhf_core_it). Only the X/+ tasks carry the
  // continuum channel cache (Y/- partners are always bound). Y tasks are
  // skipped when Y is frozen (staged iteration). In seeded (homogeneous)
  // mode there is no external field: hFb is nullptr for every task.
  struct MsTask {
    DiracSpinor *dF;           // target (in Xs or Ys)
    CntmChannel *ch;           // continuum cache (nullptr for Y)
    const DiracSpinor *Fb;     // core orbital
    const DiracSpinor *hFb;    // source projection h|Fb> (nullptr if seeded)
    const DiracSpinor *dV_src; // this task's [dV phi]_beta (dV_rhs_all)
    dPsiType type;
  };
  std::vector<MsTask> tasks;
  for (auto ib = 0ul; ib < m_core.size(); ib++) {
    for (auto be = 0ul; be < Xs[ib].size(); be++) {
      tasks.push_back({&Xs[ib][be], &m_ch[ib][be], &m_core[ib],
                       m_seed ? nullptr : &m_hFcore[ib][be], &dv.X[ib][be],
                       dPsiType::X});
      if (include_Y) {
        tasks.push_back({&Ys[ib][be], nullptr, &m_core[ib],
                         m_seed ? nullptr : &m_hFcore_minus[ib][be],
                         &dv.Y[ib][be], dPsiType::Y});
      }
    }
  }

  // Solve for the (undamped) dF of this iteration.
#pragma omp parallel for schedule(dynamic)
  for (auto it = 0ul; it < tasks.size(); it++) {
    const auto &t = tasks[it];
    solve_channel_cntm(t.dF, t.ch, *t.Fb, t.hFb, *t.dV_src, omega, t.type,
                       eps_ms);
  }

  // Measure convergence on the undamped X solution (and undamped K), before
  // damping.
  const auto eps = eps_cntm(Xs, K_old);

  // Damp the solution (for the next iteration) and store. K is linear in the
  // spinor, so it is damped with the same weights (keeps the stored K
  // consistent with the stored X).
  for (auto ib = 0ul; ib < m_ch.size(); ib++) {
    for (auto be = 0ul; be < m_ch[ib].size(); be++) {
      m_ch[ib][be].K =
        eta_damp * K_old[ib][be] + (1.0 - eta_damp) * m_ch[ib][be].K;
    }
  }
#pragma omp parallel for
  for (auto ib = 0ul; ib < m_core.size(); ib++) {
    Xs[ib] = eta_damp * m_X[ib] + (1.0 - eta_damp) * Xs[ib];
    if (include_Y) {
      Ys[ib] = eta_damp * m_Y[ib] + (1.0 - eta_damp) * Ys[ib];
    }
  }
  m_X = std::move(Xs);
  if (include_Y) {
    m_Y = std::move(Ys);
  }

  return eps;
}

//==============================================================================
std::pair<double, std::string>
TDHFcntm::eps_cntm(const std::vector<std::vector<DiracSpinor>> &Xnew,
                   const std::vector<std::vector<double>> &K_old) const {
  // Bound channels: relative L2 change of the (undamped) X spinors, summed:
  // ratio = Sum|dX|^2 / Sum|X_new|^2 (as TDHF::eps_dPsi), where the
  // denominator runs over ALL X channels (bound AND open): when only
  // physically-negligible bound channels remain (e.g. just the tiny 1s
  // response once the outer shells are open), their noise-level relative
  // change must not gate the TDHF -- their dV impact is norm-weighted, so
  // norm-weighting the measure is the consistent choice. Open channels:
  // per-channel squared relative change of the standing-wave amplitude,
  // (dK/K)^2 -- box-independent and directly physical (K = pi*D).
  // Returns the worse of the two families (sqrt if m_eps_sqrt), with the
  // worst channel's label.
  double DdF2 = 0.0;
  double dF2 = 0.0;
  double worst = 0.0;
  double worst_K = 0.0;
  std::string worst_lab;

  // Scale for the open-channel measure: the largest |K| across all open
  // channels. Floors the per-channel denominator, so a channel whose K
  // passes through zero (e.g. near a Cooper-type minimum) cannot blow up
  // the relative measure while the physics is converged.
  double K_max2 = 0.0;
  for (const auto &chs : m_ch) {
    for (const auto &ch : chs) {
      if (ch.open && ch.K * ch.K > K_max2) {
        K_max2 = ch.K * ch.K;
      }
    }
  }
  const auto K_floor2 = 1.0e-4 * K_max2;

  for (std::size_t ib = 0; ib < m_core.size(); ++ib) {
    for (std::size_t i = 0; i < Xnew[ib].size(); ++i) {
      const auto &nw = Xnew[ib][i];
      // Norm scale: ALL X channels contribute to the denominator
      dF2 += nw.norm2();
      if (m_ch[ib][i].open) {
        const auto Kn = m_ch[ib][i].K;
        const auto dK = Kn - K_old[ib][i];
        const auto den = std::max(Kn * Kn, K_floor2);
        const auto e = den == 0.0 ? 0.0 : (dK * dK) / den;
        if (e > worst_K) {
          worst_K = e;
        }
        if (e > worst) {
          worst = e;
          worst_lab = m_core[ib].shortSymbol() + "," + nw.shortSymbol();
        }
      } else {
        const auto d = (nw - m_X[ib][i]).norm2();
        const auto n = nw.norm2();
        DdF2 += d;
        // per-channel ratio: used only to rank the worst channel
        const auto e = n == 0.0 ? 0.0 : d / n;
        if (e > worst) {
          worst = e;
          worst_lab = m_core[ib].shortSymbol() + "," + nw.shortSymbol();
        }
      }
    }
  }
  const auto ratio = dF2 == 0.0 ? 0.0 : DdF2 / dF2;
  const auto eps = std::max(ratio, worst_K);
  return {m_eps_sqrt ? std::sqrt(eps) : eps, worst_lab};
}

//==============================================================================
std::vector<TDHFcntm::OpenChannel>
TDHFcntm::open_channels(const DiracSpinor &Fa) const {
  // Per-open-channel K and the internal dressed amplitude D (see hpp).
  // D is built from the solver's own (local-KS) F_reg and the TOTAL
  // effective source of the converged solve,
  //   S = (t + dV')phi_a + (V^nl - X_a - U_KS) chi ,
  // i.e. including the iterated nonlocal-exchange remainder: the radial
  // solve satisfies K = pi * <F_reg|S> identically (surface-term identity
  // for the (F_reg, F_irr) eigenpair of the local operator), so the K = pi*D
  // residual measures only the solve/extraction accuracy, not the exchange.
  using namespace qip::overloads;
  std::vector<OpenChannel> out;
  if (m_ch.empty() || m_hFcore.empty() || m_omega < 0.0) {
    return out;
  }
  const auto ib = static_cast<std::size_t>(
    std::find(m_core.cbegin(), m_core.cend(), Fa) - m_core.cbegin());
  assert(ib < m_ch.size());
  const auto Ux = HF::vex_KS(m_core);
  for (std::size_t be = 0; be < m_ch[ib].size(); ++be) {
    const auto &ch = m_ch[ib][be];
    if (!ch.open || ch.Freg.norm2() == 0.0) {
      continue;
    }
    const auto kappa = ch.Freg.kappa();
    const auto &chi = m_X[ib][be];
    // Physical source (with the diagonal de projection, mirroring the
    // channel solve), plus compensation and exchange remainder:
    auto S_t = m_hFcore[ib][be] + dV_rhs(kappa, Fa, false);
    if (kappa == Fa.kappa()) {
      S_t -= (Fa * S_t) * Fa;
    }
    const auto S = S_t + hole_compensation(Fa, chi) + HF::vexFa(chi, m_core) -
                   HF::vexFa_1el(chi, Fa) - (Ux * chi);
    const auto D = ch.Freg * S;
    out.push_back({kappa, Fa.en() + m_omega, ch.K, D});
  }
  return out;
}

//==============================================================================
void TDHFcntm::solve_channel_cntm(DiracSpinor *dF_beta, CntmChannel *ch,
                                  const DiracSpinor &Fb, const DiracSpinor *hFb,
                                  const DiracSpinor &dV_src, const double omega,
                                  dPsiType XorY, double eps_ms) const {
  // Continuum (V^{N-1}) dispatch for a single channel; see class docs.
  // Ionised orbital: en_+ = en_b + omega > 0. The one-electron (spherically
  // averaged) self-interaction V^a_0 = y^0_aa + X_a of the hole (= Fb) is
  // subtracted from the static Hamiltonian AND compensated (lagged) in the
  // source (hole_compensation, via dVprime_rhs) -- an exact rearrangement
  // of the TDHF equations:
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

  const auto ww = XorY == dPsiType::X ? omega : -omega;
  auto conj = XorY == dPsiType::Y;
  if (omega < 0.0)
    conj = !conj;
  const auto imag = m_h->imaginaryQ();

  // Only channels with e0 = en_b + ww > 0 are open; for an ionised orbital
  // that is (all of) its X/+ channels. The V^{N-1} rearrangement is applied
  // per channel, and ONLY to the open (continuum) channels, where it is
  // needed for the residual-ion (Z_ion = 1) boundary condition. Every bound
  // channel -- closed orbitals AND the bound Y/- partners of ionised
  // orbitals -- is solved in plain V^N form, exactly as bound TDHF: the
  // rearrangement is channel-wise exact, so the fixed point is unchanged,
  // but the lagged compensation must not run through the bound channels'
  // near-resonant conditioning (transient mismatches between near-degenerate
  // X/Y partner channels are amplified by ~1/(e0-em) and can diverge the
  // outer TDHF where plain bound TDHF converges).
  const bool openQ = Fb.en() + ww > 0.0;

  if (openQ && (m_suppress_open || ch == nullptr || !ch->open)) {
    // Open channel excluded from the response (explicit zero: clears nan):
    // either by the suppress_open option, or unresolvable on this radial
    // box (en_+ too close to threshold; see prepare_channels).
    dF_beta->f().assign(dF_beta->f().size(), 0.0);
    dF_beta->g().assign(dF_beta->g().size(), 0.0);
    if (ch) {
      ch->K = 0.0;
    }
    return;
  }

  const auto s = (imag && conj) ? -1.0 : 1.0;

  if (openQ) {
    using namespace qip::overloads;
    assert(ch != nullptr && "open X channel requires its CntmChannel cache");
    const int kappa_beta = dF_beta->kappa();
    // Physical source: [(t + dV)phi_a]_beta
    auto rhs = dV_src;
    if (hFb) {
      rhs += s * (*hFb);
    }
    // Diagonal channel (even-parity operators only): subtract the
    // norm-conservation (Lagrange multiplier) term de*phi_a, with
    // de = <a|(t + dV)phi_a>. Applies to the physical source ONLY -- the
    // V^a_0 compensation below is part of the operator rearrangement, and
    // projecting it too would over-subtract by <a|V^a_0 chi>*phi_a. The
    // exact counterpart of the bound solver's conditioning projection;
    // guarantees e.g. zero response for t -> identity (j0(qr), q -> 0).
    if (kappa_beta == Fb.kappa()) {
      rhs -= (Fb * rhs) * Fb;
    }
    // dF_beta still holds the previous iterate: the lagged chi in dV'
    // (in seeded mode, the TOTAL seed + correction for the seeded channel)
    rhs += hole_compensation(Fb, *dF_beta);
    const auto y0aa = Coulomb::yk_ab(0, Fb, Fb);
    const auto vl_c = p_hf->vlocal(Angular::l_k(kappa_beta)) - y0aa;
    // Seeded (homogeneous) channel: solve for the CORRECTION only -- the
    // seed is a solution of the conditioning-potential equation, so its
    // exchange deficit (vnl - U_KS)*seed joins the source; the total
    // (seed + correction) is stored back so dV sees the seed. The K
    // written to the cache is the correction's cos-amplitude = Kbar_ij.
    const bool seeded = m_seed && ch == &m_ch[m_seed->ib][m_seed->be];
    if (seeded) {
      rhs += m_seed->deficit;
      auto corr = *dF_beta - m_seed->seed;
      ExternalField::solveContinuumMixedState(&corr, &ch->Freg, &ch->Firr,
                                              &ch->K, Fb, ww, vl_c, m_alpha,
                                              m_core, rhs, eps_ms, &Fb);
      *dF_beta = corr + m_seed->seed;
      return;
    }
    // Freg/Firr reused from the cache (built in prepare_channels); the
    // standing-wave K amplitude of this channel is written back to the cache.
    ExternalField::solveContinuumMixedState(dF_beta, &ch->Freg, &ch->Firr,
                                            &ch->K, Fb, ww, vl_c, m_alpha,
                                            m_core, rhs, eps_ms, &Fb);
    return;
  }

  // Bound channel (closed orbital, or the Y/- partner of an ionised
  // orbital): plain V^N bound solve, as bound TDHF. (solveMixedState_cntm:
  // Anderson mixing, robust against the spurious-preconditioner divergence
  // at high omega.)
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
DiracSpinor TDHFcntm::dVprime_rhs(int kappa_beta, const DiracSpinor &Fa,
                                  bool conj, const DiracSpinor &chi) const {
  return dV_rhs(kappa_beta, Fa, conj) + hole_compensation(Fa, chi);
}

//==============================================================================
DiracSpinor TDHFcntm::hole_compensation(const DiracSpinor &Fa,
                                        const DiracSpinor &chi) const {
  // + V^a_0 chi = (y^0_aa + X_a) chi, the lagged compensation for the
  // one-electron self-interaction moved into the static V^{N-1} Hamiltonian.
  using namespace qip::overloads;
  const auto y0aa = Coulomb::yk_ab(0, Fa, Fa);
  return (y0aa * chi) + HF::vexFa_1el(chi, Fa);
}

//==============================================================================
void TDHFcntm::solve_homogeneous(std::size_t ib, std::size_t be, int max_its,
                                 bool print) {
  // Seeded homogeneous TDHF (Johnson and Cheng 1979, appendix): no external
  // field; channel (ib, be) is seeded with its regular continuum orbital
  // (the cached Freg -- KS phase reference). m_X stores the TOTAL for the
  // seeded channel, so dV is built from (seed + correction); the channel's
  // own solve is for the correction, with the seed's exchange deficit
  // (vnl - U_KS)*seed carried in the source (the seed solves the local
  // conditioning-potential equation exactly, not the full HF one).
  // At convergence, the per-channel K amplitudes in m_ch are column
  // (ib, be) of the rescattering matrix Kbar. Mutates *this; call on a copy.
  using namespace qip::overloads;
  assert(!m_ch.empty() && m_omega >= 0.0 &&
         "solve_homogeneous requires prepared channels (run solve_core)");
  assert(m_ch[ib][be].open && m_ch[ib][be].Freg.norm2() != 0.0);

  const auto &Freg = m_ch[ib][be].Freg;
  const auto Ux = HF::vex_KS(m_core);
  auto deficit =
    HF::vexFa(Freg, m_core) - HF::vexFa_1el(Freg, m_core[ib]) - (Ux * Freg);
  auto seed = Freg;
  // Diagonal seeded channel (even-parity operators): orthogonalise the seed
  // against phi_a, so the converged total w = seed + correction satisfies
  // the norm-conservation constraint <a|w> = 0 exactly (every solved piece
  // is kept orthogonal to phi_a; the raw F_reg is not). The seed then no
  // longer solves the local equation: the extra term joins the deficit,
  //   (h^{N-1} - en_+)(F_reg - c*phi_a) = deficit + c*(omega + V^a_0)phi_a,
  // using h_HF phi_a = en_a phi_a. phi_a decays, so the asymptotics (and
  // the Kbar decomposition) are unchanged.
  if (seed.kappa() == m_core[ib].kappa()) {
    const auto c_a = m_core[ib] * seed;
    seed -= c_a * m_core[ib];
    deficit +=
      c_a * (m_omega * m_core[ib] + hole_compensation(m_core[ib], m_core[ib]));
  }

  // Fresh start: zero all corrections and K amplitudes (keep the channel
  // caches), then insert the seed as the seeded channel's total.
  TDHF::clear();
  for (auto &chs : m_ch) {
    for (auto &ch : chs) {
      ch.K = 0.0;
    }
  }
  m_X[ib][be] = seed;
  m_seed = SeedInfo{ib, be, seed, deficit};

  solve_core_anderson(m_omega, max_its, print);

  m_seed.reset();
}

//==============================================================================
TDHFcntm::Rescattering TDHFcntm::rescattering(int max_its, bool print,
                                              bool parallel) const {
  // On-shell rescattering matrix and unitarised (physical) amplitudes; see
  // hpp. One seeded homogeneous solve per open channel gives one column of
  // Kbar; the driven amplitudes D = K/pi are read from the converged
  // channel caches of *this (which is left untouched: the homogeneous
  // solves run on copies).
  Rescattering out;
  if (m_ch.empty() || m_omega < 0.0) {
    return out;
  }

  // Open (and usable) channels, in fixed (ib, be) order
  std::vector<std::pair<std::size_t, std::size_t>> idx;
  for (auto ib = 0ul; ib < m_ch.size(); ib++) {
    for (auto be = 0ul; be < m_ch[ib].size(); be++) {
      const auto &ch = m_ch[ib][be];
      if (ch.open && ch.Freg.norm2() != 0.0) {
        idx.emplace_back(ib, be);
        out.channels.push_back(
          {ib, m_X[ib][be].kappa(), m_core[ib].en() + m_omega});
        out.D.push_back(ch.K / M_PI);
      }
    }
  }
  const auto Np = idx.size();
  if (Np == 0) {
    return out;
  }

  // Column j of Kbar: seeded solve for open channel j (on a copy). The
  // columns are independent, so they can run in parallel (printing forces
  // serial, so the labelled convergence lines stay readable). Inside an
  // active parallel region (e.g. omega-parallel caller) the nested pragma
  // is inert and the seeds run serially on that thread, as before.
  const bool par_seeds = parallel && !print;
  out.Kbar = LinAlg::Matrix<double>(Np, Np);
#pragma omp parallel for schedule(dynamic) if (par_seeds)
  for (auto j = 0ul; j < Np; j++) {
    auto homog = *this;
    homog.solve_homogeneous(idx[j].first, idx[j].second, max_its, print);
    for (auto i = 0ul; i < Np; i++) {
      out.Kbar(i, j) = homog.m_ch[idx[i].first][idx[i].second].K;
    }
  }

  // Symmetry diagnostic: Kbar is symmetric in exact arithmetic
  double kmax = 0.0, asym = 0.0;
  for (auto i = 0ul; i < Np; i++) {
    for (auto ij = 0ul; ij < Np; ij++) {
      kmax = std::max(kmax, std::abs(out.Kbar(i, ij)));
      asym = std::max(asym, std::abs(out.Kbar(i, ij) - out.Kbar(ij, i)));
    }
  }
  out.asymmetry = kmax == 0.0 ? 0.0 : asym / kmax;

  // A = (1 - i*Kbar)^{-1} pi*D, in real arithmetic:
  // (1 + Kbar*Kbar) Re(A) = pi*D,  Im(A) = Kbar * Re(A)
  LinAlg::Matrix<double> B = out.Kbar * out.Kbar;
  for (auto i = 0ul; i < Np; i++) {
    B(i, i) += 1.0;
  }
  LinAlg::Vector<double> piD(Np);
  for (auto i = 0ul; i < Np; i++) {
    piD(i) = M_PI * out.D[i];
  }
  const auto Ar = LinAlg::solve_Axeqb(B, piD);
  out.D_phys.resize(Np);
  for (auto i = 0ul; i < Np; i++) {
    double Ai = 0.0;
    for (auto j = 0ul; j < Np; j++) {
      Ai += out.Kbar(i, j) * Ar(j);
    }
    out.D_phys[i] = std::sqrt(Ar(i) * Ar(i) + Ai * Ai) / M_PI;
  }

  return out;
}

//==============================================================================
double TDHFcntm::D_phys(const DiracSpinor &Fa, int kappa_e) {
  // Unitarised (physical) amplitude |A|/pi for the (Fa, kappa_e) open
  // channel; see hpp. All channels are computed together (the seeded solves
  // couple them), so the rescattering runs once and is cached.
  if (!m_resc) {
    m_resc = rescattering(m_max_its, false, true);
  }
  const auto ib = static_cast<std::size_t>(
    std::find(m_core.cbegin(), m_core.cend(), Fa) - m_core.cbegin());
  for (std::size_t i = 0; i < m_resc->channels.size(); ++i) {
    if (m_resc->channels[i].i_core == ib &&
        m_resc->channels[i].kappa == kappa_e) {
      return m_resc->D_phys[i];
    }
  }
  return 0.0;
}

//==============================================================================
double TDHFcntm::dV_cntm(const DiracSpinor &Fe, const DiracSpinor &Fa) const {
  // As TDHF::dV(Fe,Fa), but for a continuum Fe: the source carries the
  // hole compensation, matching the V^{N-1} rearrangement (see class docs).
  // Its +y^0_aa*chi term cancels the 1/r tail hidden in the b=a part of
  // dV*phi_a pointwise, so the continuum-continuum overlap is
  // box-independent; Fe must be the V^{N-1} continuum state (hole = Fa,
  // incl. the exchange part).
  const auto conj = Fa.en() > Fe.en();
  const auto s = conj && m_h->imaginaryQ() ? -1 : 1;
  const auto ChiType = !conj ? dPsiType::X : dPsiType::Y;
  const auto &chi = get_dPsi_x(Fa, ChiType, Fe.kappa());
  return s * (Fe * dVprime_rhs(Fe.kappa(), Fa, conj, chi));
}

} // namespace ExternalField
