#include "HF/HartreeFock.hpp"
#include "Angular/Angular_369j.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracODE/Adams_Greens.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/Breit.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/radiativePotential.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

namespace HF {

//******************************************************************************
Method parseMethod(const std::string &in_method) {
  if (in_method == "HartreeFock")
    return Method::HartreeFock;
  if (in_method == "ApproxHF")
    return Method::ApproxHF;
  if (in_method == "Hartree")
    return Method::Hartree;
  std::cout << "Warning: HF Method: " << in_method << " ?? Defaulting to HF\n";
  return Method::HartreeFock;
}

//******************************************************************************
//******************************************************************************
HartreeFock::HartreeFock(const Grid &in_grid,
                         const std::vector<double> &in_vnuc,
                         std::vector<DiracSpinor> *in_core,
                         const RadiativePotential::Vrad *const in_vrad,
                         double in_alpha, Method method, bool incl_Breit,
                         double in_eps)
    : p_rgrid(&in_grid),
      p_vnuc(&in_vnuc), // or, just have a copy?
      p_vrad(in_vrad),
      m_alpha(in_alpha),
      m_method(method),
      m_include_Breit(incl_Breit),
      m_eps_HF(std::abs(in_eps) < 1.0 ? in_eps : std::pow(10, -in_eps)),
      p_core(in_core),
      m_vdir(in_grid.num_points),
      m_Yab(p_rgrid, p_core),
      m_excludeExchange(
          !(m_method == Method::HartreeFock || m_method == Method::ApproxHF)) {}

//------------------------------------------------------------------------------
HartreeFock::HartreeFock(Wavefunction *wf, Method method, bool incl_Breit,
                         double eps)
    : HartreeFock(wf->rgrid, wf->vnuc, &wf->core, &wf->vrad, wf->alpha, method,
                  incl_Breit, eps) {}

//******************************************************************************
const std::vector<double> &HartreeFock::solveCore() {

  // Core orbs must already be solutions... this OK?

  if (p_core->empty())
    return m_vdir;

  switch (m_method) {
  case Method::HartreeFock:
    hf_core_approx(m_eps_HF);
    hf_core_refine();
    break;
  case Method::ApproxHF:
    hf_core_approx(m_eps_HF);
    break;
  case Method::Hartree:
    hf_core_approx(m_eps_HF);
    break;
  default:
    m_Yab.update_y_ints(); // needed?
  }

  // Frozen core Breit:
  // Once core HF done, core "Fronen", Breit operator created
  // (Temporary VBr used in HF routine)
  if (m_include_Breit)
    m_VBr = std::make_unique<HF::Breit>(*p_core);

  return m_vdir;
}

//******************************************************************************
std::vector<double> HartreeFock::get_vlocal(int l) const {
  const auto &vrad_el = get_Hrad_el(l);
  return NumCalc::add_vectors(*p_vnuc, m_vdir, vrad_el);
}

//******************************************************************************
DiracSpinor HartreeFock::VBr(const DiracSpinor &Fv) const {
  if (m_VBr)
    return (*m_VBr)(Fv);
  else
    return 0.0 * Fv;
}

//******************************************************************************
void HartreeFock::hf_core_approx(const double eps_target_HF) {
  auto sp = IO::Profile::safeProfiler(__func__);
  if (p_core->empty()) {
    return;
  }

  auto damper = rampedDamp(0.7, 0.2, 3, 10);
  // don't include all pts in PT for new e guess:
  static const std::size_t de_stride = 5;

  // initialise 'old' potentials
  auto vdir_old = m_vdir;
  std::vector<std::vector<double>> appr_vex_core(p_core->size());
  auto vex_old = appr_vex_core;

  // Start the HF itterative procedure:
  int hits = 1;
  double t_eps = 1.0;
  auto t_eps_prev = 1.0;
  for (; hits < m_max_hf_its; hits++) {
    auto eta = damper(hits);

    // Store old vdir/vex
    vdir_old = m_vdir;
    vex_old = appr_vex_core;

    // Form new v_dir and v_ex:
    m_Yab.update_y_ints();
    form_vdir(m_vdir, false);
    form_approx_vex_core(appr_vex_core);
    if (hits == 1)
      vex_old = appr_vex_core; // We didn't have old vex before

    for (std::size_t j = 0; j < p_rgrid->num_points; j++) {
      m_vdir[j] = (1.0 - eta) * m_vdir[j] + eta * vdir_old[j];
      for (std::size_t i = 0; i < p_core->size(); i++) {
        appr_vex_core[i][j] =
            (1.0 - eta) * appr_vex_core[i][j] + eta * vex_old[i][j];
      }
    }

    const auto v_local = NumCalc::add_vectors(*p_vnuc, m_vdir);

    // Solve Dirac Eq. for each state in core, using Vdir+Vex:
    t_eps = 0;
    for (std::size_t i = 0; i < p_core->size(); i++) {
      auto &Fa = (*p_core)[i];
      double en_old = Fa.en;
      // calculate de from PT
      double dEa = 0;
      for (std::size_t j = 0; j < Fa.pinf; j += de_stride) {
        double dv =
            (m_vdir[j] - vdir_old[j]) + (appr_vex_core[i][j] - vex_old[i][j]);
        dEa += dv * Fa.f[j] * Fa.f[j] * p_rgrid->drdu[j];
      }
      dEa *= p_rgrid->du * de_stride;
      double en_guess = (en_old < -dEa) ? en_old + dEa : en_old;
      const auto &vrad_el = get_Hrad_el(Fa.l());
      const auto &vrad_mag = get_Hrad_mag(Fa.l());
      const auto v = NumCalc::add_vectors(v_local, vrad_el, appr_vex_core[i]);
      DiracODE::boundState(Fa, en_guess, v, vrad_mag, m_alpha, 6);
      double state_eps = fabs((Fa.en - en_old) / en_old);
      // convergance based on worst orbital:
      t_eps = (state_eps > t_eps) ? state_eps : t_eps;
      if constexpr (print_each_eps) {
        std::cout << __LINE__ << "| ";
        printf(" --- %2i,%2i: en=%11.5f  HFeps = %.0e;  Adams = %.0e[%2i]  "
               "(%4i)\n",
               Fa.n, Fa.k, Fa.en, state_eps, Fa.eps, Fa.its, (int)Fa.pinf);
      }
    } // core states
    if constexpr (print_each_eps) {
      std::cerr << __LINE__ << "| "
                << "HF core it: " << hits << ": eps=" << t_eps << "\n";
    }

    // Force all core orbitals to be orthogonal to each other
    if (m_explicitOrthog_cc)
      Wavefunction::orthonormaliseOrbitals((*p_core), 1);
    auto getting_worse = (hits > 20 && t_eps > t_eps_prev && t_eps < 1.e-5);
    auto converged = (t_eps < eps_target_HF);
    if (converged || getting_worse)
      break;
    t_eps_prev = t_eps;
  } // hits
  if (verbose && m_method != Method::HartreeFock)
    printf("HF core      it:%3i eps=%6.1e\n", hits, t_eps);
  if constexpr (print_final_eps) {
    printf("HF core (approx)  it:%3i eps=%6.1e\n", hits, t_eps);
  }

  // Now, re-solve core orbitals with higher precission
  for (std::size_t i = 0; i < p_core->size(); i++) {
    auto &Fa = (*p_core)[i];
    const auto &vrad_el = get_Hrad_el(Fa.l());
    const auto &vrad_mag = get_Hrad_mag(Fa.l());
    const auto v =
        NumCalc::add_vectors(*p_vnuc, m_vdir, vrad_el, appr_vex_core[i]);
    DiracODE::boundState(Fa, Fa.en, v, vrad_mag, m_alpha, 15);
  }
  if (m_explicitOrthog_cc)
    Wavefunction::orthonormaliseOrbitals((*p_core), 2);
}

//******************************************************************************
void HartreeFock::solveValence(std::vector<DiracSpinor> *valence,
                               const bool print) {
  auto sp = IO::Profile::safeProfiler(__func__);

  if (valence == nullptr || valence->empty())
    return;

  const auto Nval = valence->size();
  auto do_refine = (m_method == Method::HartreeFock && !p_core->empty());

  std::vector<EpsIts> eis(Nval);

#pragma omp parallel for
  for (std::size_t i = 0; i < Nval; i++) {
    auto &Fa = (*valence)[i];
    eis[i] = hf_valence_approx(Fa, m_eps_HF);
    if (do_refine)
      eis[i] = hf_valence_refine(Fa);
  }

  double eps_worst = 0.0, eps_best = 10.0;
  std::size_t i_worst = 0, i_best = 0;
  for (std::size_t i = 0; i < Nval; i++) {
    const auto &ei = eis[i];
    if (ei.eps >= eps_worst) {
      eps_worst = ei.eps;
      i_worst = i;
    }
    if (ei.eps < eps_best) {
      eps_best = ei.eps;
      i_best = i;
    }
  }
  if (verbose && print)
    printf("HF valence:  %3i eps=%6.1e for %s  [%6.1e for %s w/%3i]\n", //
           eis[i_worst].its, eis[i_worst].eps,
           valence->at(i_worst).symbol().c_str(), eis[i_best].eps,
           valence->at(i_best).symbol().c_str(), eis[i_best].its);
}
//******************************************************************************
void HartreeFock::solveBrueckner(std::vector<DiracSpinor> *valence,
                                 const MBPT::CorrelationPotential &Sigma2,
                                 const bool print) {
  auto sp = IO::Profile::safeProfiler(__func__);

  if (valence->empty())
    return;

  if (print) {
    std::cout << "\nSolving for Brueckner orbitals (correlation potential)\n";
    Sigma2.print_scaling();
  }
  for (auto &Fv : *valence) {
    const auto en_old = Fv.en;
    if (print)
      std::cout << Fv.symbol() << ":" << std::flush;
    const auto eis = hf_Brueckner(Fv, Sigma2);
    const auto delta = Fv.en - en_old;
    if (print) {
      printf(" delta=%8.5f; eps=%6.1e [its=%3i]", delta, eis.eps, eis.its);
      if (eis.eps > m_eps_HF && eis.eps > 1.0e-6) {
        std::cout << " *****";
      }
      std::cout << "\n";
    }
  }
}

//******************************************************************************
EpsIts HartreeFock::hf_valence_approx(DiracSpinor &Fa, double eps_target_HF)
// Solves HF for given orbital Fa, in frozen core.
// Does not store vex (must be done outside)
// Can be used to generate a set of virtual/basis orbitals
{
  auto sp = IO::Profile::safeProfiler(__func__);
  Fa.occ_frac = 1.0 / Fa.twojp1();

  auto damper = rampedDamp(0.7, 0.2, 2, 6);
  // don't include all pts in PT for new e guess
  static const std::size_t de_stride = 5;

  std::vector<double> vexa(p_rgrid->num_points, 0);
  auto vexa_old = vexa;

  const auto &vrad_el = get_Hrad_el(Fa.l());
  const auto &Hmag = get_Hrad_mag(Fa.l());
  const auto v_local = NumCalc::add_vectors(*p_vnuc, m_vdir, vrad_el);

  double eps = -1, eps_prev = -1;
  int hits = 1;
  for (; hits < m_max_hf_its; hits++) {
    auto eta = damper(hits);

    double en_old = Fa.en;
    vexa_old = vexa;
    if (!m_excludeExchange)
      vexa = vex_approx(Fa, *p_core);

    for (std::size_t i = 0; i < p_rgrid->num_points; i++) {
      vexa[i] = (1.0 - eta) * vexa[i] + eta * vexa_old[i];
    }
    // Use P.T. to calculate energy change:
    double en_new_guess = 0;
    for (std::size_t i = 0; i < Fa.pinf; i += de_stride) {
      en_new_guess +=
          (vexa[i] - vexa_old[i]) * Fa.f[i] * Fa.f[i] * p_rgrid->drdu[i];
    }
    en_new_guess = en_old + en_new_guess * p_rgrid->du * de_stride;
    // Solve Dirac using new potential:
    DiracODE::boundState(Fa, en_new_guess, NumCalc::add_vectors(v_local, vexa),
                         Hmag, m_alpha, 15);
    eps = fabs((Fa.en - en_old) / en_old);
    // Force valence state to be orthogonal to core:
    if constexpr (m_explicitOrthog_cv)
      Wavefunction::orthonormaliseWrt(Fa, (*p_core));

    auto getting_worse = (hits > 20 && eps >= eps_prev && eps < 1.e-5);
    auto converged = (eps <= eps_target_HF);
    if (converged || getting_worse)
      break;
    eps_prev = eps;
  }
  if constexpr (print_final_eps) {
    printf("HF val: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", Fa.n, Fa.k, hits,
           eps, Fa.en);
  }

  if constexpr (m_explicitOrthog_cv)
    Wavefunction::orthonormaliseWrt(Fa, (*p_core));
  return {eps, hits};
}

//******************************************************************************
double HartreeFock::calculateCoreEnergy() const
// Calculates the total HF core energy:
//   E = \sum_a [ja]e_a - 0.5 \sum_(ab) (R^0_abab - \sum_k L^k_ab R^k_abba)
// where:
//   R^k_abcd = Integral [f_a*f_c + g_a*g_c] * v^k_bd
//   R^0_abab is not absymmetric
//   R^k_abba _is_ ab symmetric
{
  double Etot = 0.0;
  for (const auto &Fa : (*p_core)) {
    const auto tja = Fa.twoj();

    double e1 = 0.0, e2 = 0.0, e3 = 0.0;
    const double xtjap1 = (tja + 1) * Fa.occ_frac;
    e1 += xtjap1 * Fa.en;
    for (const auto &Fb : *p_core) {
      const auto tjb = Fb.twoj();
      const double xtjbp1 = (tjb + 1) * Fb.occ_frac;
      const auto &v0bb = m_Yab.get_yk_ab(0, Fb, Fb);
      const auto R0fg2 = Fa * (v0bb * Fa);
      e2 += xtjap1 * xtjbp1 * R0fg2;
      // take advantage of symmetry for third term:
      if (Fb > Fa)
        continue;
      const auto y = (Fa == Fb) ? 1.0 : 2.0; // symmetry
      const int kmin = std::abs(tja - tjb) / 2;
      const int kmax = (tja + tjb) / 2;
      const auto &vabk = m_Yab.get_y_ab(Fa, Fb);
      for (auto k = kmin; k <= kmax; k++) {
        const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.k, Fb.k);
        if (Labk == 0)
          continue;
        const auto ik = std::size_t(k - kmin);
        const auto R0fg3 = Fa * (vabk[ik] * Fb);
        e3 += y * xtjap1 * xtjbp1 * Labk * (R0fg3);
      }
    }
    Etot += e1 - 0.5 * (e2 - e3);
  }
  return Etot;
}

//******************************************************************************
void HartreeFock::form_vdir(std::vector<double> &vdir, bool re_scale) const
// Forms the direct part of the potential.
// Must call either form_vbb0 or form_vabk_core first!
// Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
// If re_scale==true, will scale by (N-1)/N. This then given the averaged
// Hartree potential (local, same each state, no exchange).
// re_scale=false by default
{
  auto sp = IO::Profile::safeProfiler(__func__);
  for (auto &v_dir : vdir) {
    v_dir = 0.0;
  }
  const auto Ncore = std::accumulate(
      p_core->cbegin(), p_core->cend(), 0,
      [](int n, const auto &Fa) { return n + Fa.num_electrons(); });
  const double sf = re_scale ? (1.0 - 1.0 / Ncore) : 1.0;
  for (const auto &Fb : (*p_core)) {
    const double f_sf = sf * (Fb.twoj() + 1) * Fb.occ_frac;
    const auto &v0bb = m_Yab.get_yk_ab(0, Fb, Fb);
    for (std::size_t i = 0; i < p_rgrid->num_points; i++) {
      vdir[i] += v0bb[i] * f_sf;
    }
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_core(
    std::vector<std::vector<double>> &vex) const
// Forms the 2D "approximate" exchange potential for each core state, a.
// NOTE: Must call form_vabk_core first!
// Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
{
  auto sp = IO::Profile::safeProfiler(__func__);
  vex.resize(p_core->size());
#pragma omp parallel for
  for (std::size_t a = 0; a < p_core->size(); a++) {
    form_approx_vex_core_a((*p_core)[a], vex[a]);
  }
}

std::vector<double>
HartreeFock::form_approx_vex_core_a(const DiracSpinor &Fa) const {
  std::vector<double> vex_a;
  form_approx_vex_core_a(Fa, vex_a);
  return vex_a;
}

//******************************************************************************
void HartreeFock::form_approx_vex_core_a(const DiracSpinor &Fa,
                                         std::vector<double> &vex_a) const
// XXX Fa MUST be in core
// Forms the 2D "approximate" exchange potential for given core state, a.
// Does the a=b case seperately, since it's a little simpler
// Approximate:
// In order to approximate solution to HF equations, I form "local" ex.
// potential
//   [v_ex*Fa](r) = \sum_b v_ex^(a,b)(r) * Fb(r)
// v_ex is non-local; cannot write: [v_ex*Fa](r) =/= v_ex(r)*Fa(r)
// Instead: define local approx: vex_a
//   vex_a = [v_ex*Fa](r) *(Fa/Fa^2)
//         = \sum_b v_ex^(a,b)(r)*Fb(r) * (Fa/Fa^2)
//         = \sum_b v_ex^(a,b)(r)*(Fb(r)*Fa) / Fa^2
// This vex_a is then a local potential (different for each state!) that can
// be used as an addition to local direct potential to solve Dirac Eq. as
// normal. In theory, this is exact. Clearly, however, there is an issue when
// Fa is small. Luckily, however, we don't care as much when Fa is
// small! Also, since v_ex is already small (compared to vdir), we can make
// good approximation. Therefore, I only calculate vex_a when a=b, or when
// |Fa| > 1.e3 Further, largest part of v_ex is when a=b. In this case, the
// factor=1 is exact!
{
  auto sp = IO::Profile::safeProfiler(__func__);
  vex_a.clear();
  vex_a.resize(p_rgrid->num_points);

  // don't subtract self-potential term for Hartree {match Core-Hartree}
  if (m_excludeExchange)
    return;

  const auto twoj_a = Fa.twoj();

  static auto max_abs = [](double a, double b) {
    return (std::abs(a) < std::abs(b));
  };
  const auto max =
      std::abs(*std::max_element(Fa.f.begin(), Fa.f.end(), max_abs));
  const auto cut_off = 0.003 * max;

  // bool a_in_coreQ = false;
  if (!m_excludeExchange) {
    for (const auto &Fb : (*p_core)) { // b!=a
      if (Fb == Fa) {
        // a_in_coreQ = true;
        continue;
      }
      const auto tjb = Fb.twoj();
      const double x_tjbp1 = (tjb + 1) * Fb.occ_frac;
      const auto irmax = std::min(Fa.pinf, Fb.pinf);
      const int kmin = std::abs(twoj_a - tjb) / 2;
      const int kmax = (twoj_a + tjb) / 2;
      const auto &vabk = m_Yab.get_y_ab(Fb, Fa);

      // hold "fraction" Fa*Fb/(Fa^2):
      std::vector<double> v_Fab(p_rgrid->num_points);
      for (std::size_t i = 0; i < irmax; i++) {
        // This is the approximte part! Divides by Fa
        if (std::abs(Fa.f[i]) < cut_off)
          continue;
        const auto fac_top = Fa.f[i] * Fb.f[i] + Fa.g[i] * Fb.g[i];
        const auto fac_bot = Fa.f[i] * Fa.f[i] + Fa.g[i] * Fa.g[i];
        v_Fab[i] = -x_tjbp1 * fac_top / fac_bot;
      } // r
      for (int k = kmin; k <= kmax; k++) {
        const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.k, Fb.k);
        if (Labk == 0)
          continue;
        const auto ik = std::size_t(k - kmin);
        for (std::size_t i = 0; i < irmax; i++) {
          vex_a[i] += Labk * vabk[ik][i] * v_Fab[i];
        } // r
      }   // k
    }     // b
  }

  // now, do a=b, ONLY if a is in the core!
  // if (a_in_coreQ)
  {
    const double x_tjap1 = (twoj_a + 1); // no occ_frac here ?
    const int kmax = twoj_a;
    const auto &vaak = m_Yab.get_y_ab(Fa, Fa);
    const auto irmax = Fa.pinf; // p_rgrid->num_points;
    for (int k = 0; k <= kmax; k++) {
      const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.k, Fa.k);
      if (Labk == 0)
        continue;
      for (std::size_t i = 0; i < irmax; i++) {
        // nb: If I don't 'cut' here, or fails w/ f states... ?? XX
        // Of course, cutting is fine. But WHY FAIL??
        vex_a[i] += -Labk * vaak[std::size_t(k)][i] * x_tjap1;
      }
    } // k
  }   // if a in core
}

//------------------------------------------------------------------------------
std::vector<double> vex_approx(const DiracSpinor &Fa,
                               const std::vector<DiracSpinor> &core, int k_cut)
// Free function
{
  auto sp = IO::Profile::safeProfiler(__func__);

  std::vector<double> vex(Fa.p_rgrid->num_points);
  std::vector<double> vabk;

  const auto tja = Fa.twoj();
  const auto la = Fa.l();

  static auto max_abs = [](double a, double b) {
    return (std::abs(a) < std::abs(b));
  };
  const auto max =
      std::abs(*std::max_element(Fa.f.begin(), Fa.f.end(), max_abs));
  // const auto cut_off = 0.01 * max; //why was this different?
  const auto cut_off = 0.003 * max;

  for (const auto &Fb : core) {
    const auto tjb = Fb.twoj();
    const auto lb = Fb.l();
    const double x_tjbp1 = (tjb + 1) * Fb.occ_frac; // when in core??
    const auto irmax = std::min(Fa.pinf, Fb.pinf);
    const int kmin = std::abs(tja - tjb) / 2;
    int kmax = (tja + tjb) / 2;
    if (kmax > k_cut)
      kmax = k_cut;

    // hold "fraction" Fa*Fb/(Fa^2):
    std::vector<double> v_Fab(Fa.p_rgrid->num_points);
    for (std::size_t i = 0; i < irmax; i++) {
      // This is the approximate part! Divides by Fa
      if (std::abs(Fa.f[i]) < cut_off)
        continue;
      const auto fac_top = Fa.f[i] * Fb.f[i] + Fa.g[i] * Fb.g[i];
      const auto fac_bot = Fa.f[i] * Fa.f[i] + Fa.g[i] * Fa.g[i];
      v_Fab[i] = -x_tjbp1 * fac_top / fac_bot;
    } // r

    for (int k = kmin; k <= kmax; k++) {
      const auto parity = Angular::parity(la, lb, k);
      if (parity == 0)
        continue;
      const auto tjs = Angular::threej_2(tjb, tja, 2 * k, -1, 1, 0);
      if (tjs == 0)
        continue;
      const auto tjs2 = tjs * tjs;
      Coulomb::yk_ab(Fb, Fa, k, vabk);

      for (std::size_t i = 0; i < irmax; i++) {
        if (v_Fab[i] == 0)
          continue;
        vex[i] += tjs2 * vabk[i] * v_Fab[i];
      } // r
    }   // k
  }     // b

  return vex;
}

//******************************************************************************
DiracSpinor HartreeFock::calc_vexFa_core(const DiracSpinor &Fa) const
// calculates V_ex Fa (returns new Dirac Spinor)
// Fa can be any orbital (so long as coulomb integrals exist!)
// ..... i.e., must be core orbital
{
  DiracSpinor VxFa(Fa.n, Fa.k, *(Fa.p_rgrid));
  vex_psia_core(Fa, VxFa);
  return VxFa;
}
//------------------------------------------------------------------------------
void HartreeFock::vex_psia_core(const DiracSpinor &Fa, DiracSpinor &VxFa) const
// calculates V_ex Fa
// Fa can be any orbital (so long as coulomb integrals exist!)
// ..... i.e., must be core orbital
{
  auto sp = IO::Profile::safeProfiler(__func__);
  VxFa.pinf = Fa.f.size(); // silly hack. Make sure VxFa = 0 after pinf
  VxFa *= 0.0;
  VxFa.pinf = Fa.pinf; // updated below

  if (m_excludeExchange)
    return;

  auto twoj_a = Fa.twoj();
  for (const auto &Fb : (*p_core)) {
    VxFa.pinf = std::max(VxFa.pinf, Fb.pinf);
    auto tjb = Fb.twoj();
    double x_tjbp1 = (Fa == Fb) ? (tjb + 1) : (tjb + 1) * Fb.occ_frac;
    int kmin = std::abs(twoj_a - tjb) / 2;
    int kmax = (twoj_a + tjb) / 2;
    const auto &vabk = m_Yab.get_y_ab(Fb, Fa);
    for (int k = kmin; k <= kmax; k++) {
      const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.k, Fb.k);
      if (Labk == 0)
        continue;
      for (auto i = 0u; i < Fb.pinf; i++) {
        auto v = -x_tjbp1 * Labk * vabk[std::size_t(k - kmin)][i];
        VxFa.f[i] += v * Fb.f[i];
        VxFa.g[i] += v * Fb.g[i];
      } // r
    }   // k
  }     // b
}

// -----------------------------------------------------------------------------
DiracSpinor vexFa(const DiracSpinor &Fa, const std::vector<DiracSpinor> &core,
                  int k_cut)
// Free Function
// calculates V_ex Fa (returns new Dirac Spinor)
// Fa can be any orbital (Calculates coulomb integrals here!)
{
  auto sp = IO::Profile::safeProfiler(__func__);
  DiracSpinor VxFa(Fa.n, Fa.k, *(Fa.p_rgrid));
  VxFa.pinf = Fa.pinf; // nb: updated below!
  // note: VxFa.pinf can be larger than psi.pinf!

  std::vector<double> vabk(Fa.p_rgrid->num_points);

  const auto tja = Fa.twoj();
  const auto la = Fa.l();
  for (const auto &Fb : core) {
    VxFa.pinf = std::max(Fb.pinf, VxFa.pinf);
    const auto tjb = Fb.twoj();
    const auto lb = Fb.l();
    const double x_tjbp1 = (Fa == Fb) ? (tjb + 1) : (tjb + 1) * Fb.occ_frac;
    const int kmin = std::abs(tja - tjb) / 2;
    const int kmax = std::min((tja + tjb) / 2, k_cut);
    for (int k = kmin; k <= kmax; k++) {
      if (Angular::parity(la, lb, k) == 0)
        continue;
      const auto tjs = Angular::threej_2(tjb, tja, 2 * k, -1, 1, 0);
      if (tjs == 0)
        continue;
      Coulomb::yk_ab(Fb, Fa, k, vabk, Fb.pinf);
      const auto factor = -x_tjbp1 * tjs * tjs;
      for (auto i = 0u; i < Fb.pinf; i++) {
        VxFa.f[i] += factor * vabk[i] * Fb.f[i];
        VxFa.g[i] += factor * vabk[i] * Fb.g[i];
      } // r
    }   // k
  }     // b

  return VxFa;
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
void HartreeFock::hf_orbital(DiracSpinor &Fa, double en,
                             const std::vector<double> &vl,
                             const std::vector<double> &H_mag,
                             const DiracSpinor &vx_phi,
                             const std::vector<DiracSpinor> &static_core,
                             const std::vector<double> &v0,
                             const HF::Breit *const VBr) const
// Solve Dirac Equation (Eigenvalue): (move to DiracODE??)
//  (H0 + Vl + Vx)Fa = 0
//  (H0 + Vl)Fa = -VxFa
// Vl is local (e.g., Vnuc + fVdir), Vx is non-local (e.g., (1-f)Vdir + Vex)
// where v0 = (1-f)Vdir  [f=1 for valence states!, so v0 may be empty]
// Small energy adjustmenets (and wfs), solve:
// (Hl - e) dF = de * F -VxFa
// e -> e+de, F->F+dF
// Core is input so can call in a thread-safe way! (with a 'old_core' copy)
// Only used in dE from dF
{
  auto sp = IO::Profile::safeProfiler(__func__);
  // pull these outside? But make sure thread safe!
  DiracSpinor Gzero(Fa.n, Fa.k, *(Fa.p_rgrid));
  DiracSpinor Ginf(Fa.n, Fa.k, *(Fa.p_rgrid));
  DiracSpinor VxFh(Fa.n, Fa.k, *(Fa.p_rgrid));
  DiracSpinor dFa(Fa.n, Fa.k, *(Fa.p_rgrid));
  const auto eps_target = 1.0e-17; // m_eps_HF;
  const auto k_max = 1;            // max k for Vex into dEa

  // Note: vx_phi includes VBr*Fa [if Breit].
  // VBr used in the 'small energy correction' portion

  const auto alpha = m_alpha;
  DiracODE::solve_inhomog(Fa, Gzero, Ginf, en, vl, H_mag, alpha, -1.0 * vx_phi);

  // make small adjustments to energy to normalise Fa:
  DiracODE::Adams::GreenSolution(dFa, Ginf, Gzero, alpha, Fa);
  // should dFa = dEa * dFa, but makes it worse?
  // nb: after first it, becomes correct.
  auto dEa = 0.5 * (Fa * Fa - 1.0) / (Fa * dFa);
  auto eps = std::abs(dEa / en);
  int tries = 0;
  for (; tries <= m_max_hf_its; ++tries) {
    if (eps < eps_target)
      break;

    VxFh = v0 * dFa + vexFa(dFa, static_core, k_max);
    if (VBr)
      VxFh += (*VBr)(dFa);
    DiracODE::Adams::GreenSolution(dFa, Ginf, Gzero, alpha, dEa * Fa - VxFh);

    const auto delta_Norm = Fa * Fa - 1.0;
    const auto de0 = dEa;
    dEa *= 0.5 * (delta_Norm) / (Fa * dFa);
    eps = std::abs(dEa / en);
    en += dEa;
    dFa.pinf = Fa.pinf;
    Fa -= (dEa / de0) * dFa;
  }
  Fa.en = en;
  Fa.eps = eps;
  Fa.its = tries;
  if (tries == 0 || tries == m_max_hf_its)
    Fa.normalise(); //? Not needed
}
//******************************************************************************
void HartreeFock::brueckner_orbital(DiracSpinor &Fa, double en,
                                    const std::vector<double> &vl,
                                    const std::vector<double> &H_mag,
                                    const DiracSpinor &VxF,
                                    const MBPT::CorrelationPotential &Sigma,
                                    const std::vector<DiracSpinor> &static_core,
                                    const HF::Breit *const VBr) const
// ..................
{
  auto sp = IO::Profile::safeProfiler(__func__);
  // pull these outside? But make sure thread safe!
  DiracSpinor Gzero(Fa.n, Fa.k, *(Fa.p_rgrid));
  DiracSpinor Ginf(Fa.n, Fa.k, *(Fa.p_rgrid));
  DiracSpinor dFa(Fa.n, Fa.k, *(Fa.p_rgrid));
  const auto eps_target = 1.0e-16; // m_eps_HF;

  const auto SigmaF = Sigma(Fa);

  // nb: VxF includes Breit. VBr used for energy correction

  const auto alpha = m_alpha;
  DiracODE::solve_inhomog(Fa, Gzero, Ginf, en, vl, H_mag, alpha,
                          -1.0 * VxF - SigmaF);

  // make small adjustments to energy to normalise Fa:
  DiracODE::Adams::GreenSolution(dFa, Ginf, Gzero, alpha, Fa);
  // should dFa = dEa * dFa, but makes it worse?
  // nb: after first it, becomes correct.
  auto dEa = 0.5 * (Fa * Fa - 1.0) / (Fa * dFa);
  auto eps = std::abs(dEa / en);
  int tries = 0;
  for (; tries <= m_max_hf_its; ++tries) {
    if (eps < eps_target)
      break;
    {
      auto VxF_tilde = vexFa(dFa, static_core);
      if (VBr)
        VxF_tilde += (*VBr)(dFa);
      const auto SigmaF_tilde = Sigma(dFa);
      DiracODE::Adams::GreenSolution(dFa, Ginf, Gzero, alpha,
                                     dEa * Fa - VxF_tilde - SigmaF_tilde);
    }
    const auto de0 = dEa;
    dEa *= 0.5 * (Fa * Fa - 1.0) / (Fa * dFa);
    eps = std::abs(dEa / en);
    en += dEa;
    dFa.pinf = Fa.pinf;
    Fa -= (dEa / de0) * dFa;
  }
  Fa.en = en;
  Fa.eps = eps;
  Fa.its = tries;
  // if (tries == 0 || tries == m_max_hf_its)
  Fa.normalise(); //? Not needed
}

//******************************************************************************
const static std::vector<double> tmp_empty_vector{};
const std::vector<double> &HartreeFock::get_Hrad_el(int l) const {
  return p_vrad ? p_vrad->get_Hel(l) : tmp_empty_vector;
}
const std::vector<double> &HartreeFock::get_Hrad_mag(int l) const {
  return p_vrad ? p_vrad->get_Hmag(l) : tmp_empty_vector;
}
//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
EpsIts HartreeFock::hf_valence_refine(DiracSpinor &Fa) {
  auto sp = IO::Profile::safeProfiler(__func__);
  if (p_core->empty())
    return {0, 0};

  const auto eps_target = m_eps_HF;

  const auto damper = rampedDamp(0.8, 0.2, 5, 25);
  double extra_damp = 0.0;

  const auto &vrad_el = get_Hrad_el(Fa.l());
  const auto &Hmag = get_Hrad_mag(Fa.l());
  const auto vl = NumCalc::add_vectors(*p_vnuc, m_vdir, vrad_el);

  const auto Fzero = Fa;
  const auto vexFzero = vex_approx(Fa, *p_core) * Fa;

  auto prev_en = Fa.en;
  double best_eps = 1.0;
  auto VxFa = DiracSpinor(Fa.n, Fa.k, *p_rgrid);
  int it = 0;
  double eps = 1.0;
  int worse_count = 0;
  for (; it <= m_max_hf_its; ++it) {
    const auto a_damp = damper(it) + extra_damp;

    VxFa = calc_vexFa(Fa);
    if (m_VBr) { // Breit
      VxFa += (*m_VBr)(Fa);
    }
    const auto oldphi = Fa;
    const auto en = Fzero.en + (Fzero * VxFa - Fa * vexFzero) / (Fa * Fzero);
    hf_orbital(Fa, en, vl, Hmag, VxFa, (*p_core));
    eps = std::fabs((prev_en - Fa.en) / Fa.en);
    prev_en = Fa.en;

    if (it > 20 && eps > 1.5 * best_eps) {
      ++worse_count;
      extra_damp = extra_damp > 0 ? 0 : 0.1;
    } else {
      worse_count = 0;
    }
    const bool converged = (eps <= eps_target && it > 0);
    if (converged || worse_count > 2)
      break;

    if (eps < best_eps)
      best_eps = eps;

    if constexpr (print_each_eps) {
      std::cout << __LINE__ << "| " << it << " " << eps << " " << Fa.en << " "
                << en - Fzero.en << " " << Fa * Fa << "\n";
    }

    Fa = (1.0 - a_damp) * Fa + a_damp * oldphi;
    if constexpr (m_explicitOrthog_cv) {
      Wavefunction::orthonormaliseWrt(Fa, (*p_core));
    } else {
      Fa.normalise();
    }
  } // End HF its

  if constexpr (m_explicitOrthog_cv)
    Wavefunction::orthonormaliseWrt(Fa, (*p_core));

  if constexpr (print_final_eps) {
    printf("refine: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", Fa.n, Fa.k, it, eps,
           Fa.en);
  }
  return {eps, it};
}

//******************************************************************************
EpsIts HartreeFock::hf_Brueckner(DiracSpinor &Fa,
                                 const MBPT::CorrelationPotential &Sigma) {
  auto sp = IO::Profile::safeProfiler(__func__);
  if (p_core->empty())
    return {0, 0};

  const auto eps_target = m_eps_HF;
  const auto max_br_its = m_max_hf_its;

  auto damper = rampedDamp(0.33, 0.05, 1, 15);

  const auto &vrad_el = get_Hrad_el(Fa.l());
  const auto &Hmag = get_Hrad_mag(Fa.l());
  const auto vl = NumCalc::add_vectors(*p_vnuc, m_vdir, vrad_el);

  const auto Fzero = Fa;
  auto vexFzero = calc_vexFa(Fa);
  if (m_VBr) { // Breit
    vexFzero += (*m_VBr)(Fa);
  }

  auto prev_en = Fa.en;
  double best_eps = 1.0;
  int it = 0;
  double eps = 1.0;
  int worse_count = 0;
  for (; it <= max_br_its; ++it) {
    const auto a_damp = damper(it);

    auto VxFa = calc_vexFa(Fa);
    if (m_VBr) { // Breit
      VxFa += (*m_VBr)(Fa);
    }
    const auto SigmaFa = Sigma(Fa);
    const auto oldphi = Fa;

    const auto en_new =
        Fzero.en + (Fzero * (VxFa + SigmaFa) - Fa * vexFzero) / (Fa * Fzero);
    const auto en = (Fa.en + 2.0 * en_new) / 3.0;

    brueckner_orbital(Fa, en, vl, Hmag, VxFa, Sigma, (*p_core), m_VBr.get());

    eps = std::fabs((prev_en - Fa.en) / Fa.en);
    prev_en = Fa.en;

    if (it > 20 && eps > 2.0 * best_eps) {
      ++worse_count;
    } else {
      worse_count = 0;
    }
    const bool converged = (eps <= eps_target && it > 0);
    if (converged || worse_count > 3)
      break;

    if (eps < best_eps)
      best_eps = eps;

    if constexpr (print_each_eps) {
      std::cout << __LINE__ << "| " << it << " " << eps << " " << Fa.en << " "
                << en - Fzero.en << " " << Fa * Fa << "\n";
    }

    Fa = (1.0 - a_damp) * Fa + a_damp * oldphi;
    Fa.normalise();
  } // End HF its

  if constexpr (print_final_eps) {
    printf("Br2: %2i %2i | %3i eps=%6.1e  en=%11.8f\n", Fa.n, Fa.k, it, eps,
           Fa.en);
  }
  return {eps, it};
}

//******************************************************************************
inline void HartreeFock::hf_core_refine() {
  auto sp = IO::Profile::safeProfiler(__func__);
  if (p_core->empty()) {
    return;
  }

  const double eps_target = m_eps_HF;
  m_Yab.update_y_ints(); // only needed if not already done!
  auto damper = rampedDamp(0.8, 0.2, 5, 30);
  double extra_damp = 0;

  std::vector<double> vl(p_rgrid->num_points); // Vnuc + fVd
  std::vector<double> v0(p_rgrid->num_points); // (1-f)Vd
  const auto Ncore = std::accumulate(
      p_core->cbegin(), p_core->cend(), 0,
      [](int sum, const auto &Fa) { return sum + Fa.num_electrons(); });
  const auto f_core_tmp = double(Ncore - 1) / double(Ncore);
  const auto f_core = 0.5 * (1.0 + f_core_tmp);
  const auto &vd = m_vdir;

  // Store arrays of intitial Fa and VexPsi, and VdirPsi (for En guess)
  // And allocate arrays for VexPsi, so can //-ise it loop (over orbs)!
  const auto core_zero = (*p_core);
  auto core_prev = (*p_core);
  std::vector<DiracSpinor> vexCore_zero;
  const auto vd0 = m_vdir;
  std::vector<DiracSpinor> vexF_list;
  const auto num_core_states = p_core->size();
  std::vector<double> eps_lst(num_core_states, 0.0);
  for (std::size_t i = 0; i < num_core_states; ++i) {
    auto &Fa = (*p_core)[i];
    vexCore_zero.push_back(form_approx_vex_core_a(Fa) * Fa);
    vexF_list.push_back(DiracSpinor(Fa.n, Fa.k, *(Fa.p_rgrid)));
  }

  double eps = 0.0;
  double best_eps = 1.0;
  double best_worst_eps = 1.0;
  std::size_t worst_index = 0;
  std::size_t best_index = 0;
  int worse_count = 0;
  int it = 0;
  for (; it <= m_max_hf_its; it++) {
    auto a_damp = damper(it) + extra_damp;

    // re-calculate each Vl = vnuc + fvdir, v0 = (1-f)vdir:
    for (auto i = 0ul; i < p_rgrid->num_points; i++) {
      vl[i] = (*p_vnuc)[i] + f_core * vd[i];
      v0[i] = (1.0 - f_core) * vd[i];
    }

    // re-calculate each VexPsi:
    for (std::size_t i = 0; i < num_core_states; ++i) {
      const auto &Fa = (*p_core)[i];
      vex_psia_core(Fa, vexF_list[i]);
    }

    core_prev = (*p_core);

    // Temporary Breit operator (with 'static' core [frozen single iteration])
    const std::unique_ptr VBr =
        m_include_Breit ? std::make_unique<HF::Breit>(core_prev) : nullptr;

#pragma omp parallel for
    for (std::size_t i = 0; i < num_core_states; ++i) {
      auto &Fa = (*p_core)[i];
      const auto &Fzero = core_zero[i];
      const auto &vexFzero = vexCore_zero[i];

      const auto oldphi = core_prev[i];
      const auto &VxFa = vexF_list[i];

      auto en = Fzero.en + (Fzero * VxFa - Fa * vexFzero + Fzero * (vd * Fa) -
                            Fa * (vd0 * Fzero)) /
                               (Fa * Fzero);
      auto v_nonlocal = v0 * Fa + VxFa;
      if (m_include_Breit) {
        const auto VBrFa = (*VBr)(Fa); // depends on previous core!
        en += (Fzero * VBrFa) / (Fa * Fzero);
        v_nonlocal += VBrFa;
      }

      const auto &Hrad_el = get_Hrad_el(Fa.l());
      const auto &Hmag = get_Hrad_mag(Fa.l());
      const auto &VlVr = NumCalc::add_vectors(vl, Hrad_el);
      hf_orbital(Fa, en, VlVr, Hmag, v_nonlocal, core_prev, v0, VBr.get());
      Fa = (1.0 - a_damp) * Fa + a_damp * oldphi;
      Fa.normalise();
      auto d_eps = std::fabs((oldphi.en - Fa.en) / Fa.en);
      eps_lst[i] = d_eps;
    }

    eps = eps_lst[0];
    best_eps = eps_lst[0];
    for (std::size_t i = 1; i < num_core_states; ++i) {
      auto t_eps = eps_lst[i];
      if (t_eps >= eps) {
        eps = t_eps;
        worst_index = i;
      }
      if (t_eps < best_eps) {
        best_eps = t_eps;
        best_index = i;
      }
    }
    if constexpr (print_each_eps) {
      std::cout << __LINE__ << "| " << it << " " << eps << " "
                << (*p_core)[worst_index].symbol() << " -- "
                << " " << best_eps << " " << (*p_core)[best_index].symbol()
                << "\n";
    }
    if constexpr (print_final_eps) {
      if (std::isnan(eps))
        std::cout << __LINE__ << "| eps is NaN: " << it << " " << eps << " "
                  << (*p_core)[worst_index].symbol() << " -- "
                  << " " << best_eps << " " << (*p_core)[best_index].symbol()
                  << "\n";
    }

    if (it > 20 && eps > 1.5 * best_worst_eps) {
      ++worse_count;
      extra_damp = extra_damp > 0 ? 0 : 0.4;
    } else {
      worse_count = 0;
    }
    const bool converged = (eps <= eps_target && it > 0);
    if (converged || worse_count > 3)
      break;

    if (eps < best_worst_eps)
      best_worst_eps = eps;
    if (m_explicitOrthog_cc)
      Wavefunction::orthonormaliseOrbitals((*p_core));
    m_Yab.update_y_ints();
    form_vdir(m_vdir);
  }
  if (m_explicitOrthog_cc)
    Wavefunction::orthonormaliseOrbitals((*p_core), 2);

  if (verbose)
    printf("HF core:  it:%3i eps=%6.1e for %s  [%6.1e for %s]\n", //
           it, eps, (*p_core)[worst_index].symbol().c_str(), best_eps,
           (*p_core)[best_index].symbol().c_str());
}

} // namespace HF
