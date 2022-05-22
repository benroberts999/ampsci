#include "HF/HartreeFock.hpp"
#include "Angular/Wigner369j.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracODE/Adams_Greens.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/Breit.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "Maths/Grid.hpp"
#include "Physics/Parametric_potentials.hpp"
#include "Physics/RadPot.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include "qip/Vector.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

namespace HF {

//==============================================================================
Method parseMethod(const std::string &in_method) {
  if (in_method == "HartreeFock")
    return Method::HartreeFock;
  if (in_method == "ApproxHF")
    return Method::ApproxHF;
  if (in_method == "Hartree")
    return Method::Hartree;
  if (in_method == "KohnSham")
    return Method::KohnSham;
  if (in_method == "Local")
    return Method::Local;
  std::cout << "Warning: HF Method: " << in_method << " ?? Defaulting to HF\n";
  return Method::HartreeFock;
}

std::string parseMethod(const Method &in_method) {
  if (in_method == Method::HartreeFock)
    return "HartreeFock";
  if (in_method == Method::ApproxHF)
    return "ApproxHF";
  if (in_method == Method::Hartree)
    return "Hartree";
  if (in_method == Method::KohnSham)
    return "KohnSham";
  if (in_method == Method::Local)
    return "Local";
  return "HartreeFock";
}

//==============================================================================
//==============================================================================
HartreeFock::HartreeFock(std::shared_ptr<const Grid> in_grid,
                         const std::vector<double> &in_vnuc,
                         std::vector<DiracSpinor> *in_core,
                         std::optional<QED::RadPot> in_vrad, double in_alpha,
                         Method method, double x_Breit, double in_eps)
    : rgrid(in_grid),
      p_core(in_core),
      m_vnuc(in_vnuc),
      m_vrad(in_vrad),
      m_VBr(x_Breit == 0.0 ? std::optional<HF::Breit>{} :
                             std::optional<HF::Breit>(x_Breit)),
      m_alpha(in_alpha),
      m_method(method),
      m_eps_HF(std::abs(in_eps) < 1.0 ? in_eps : std::pow(10, -in_eps)),
      m_vdir(rgrid->num_points()),
      m_Yab(*p_core) {}

//------------------------------------------------------------------------------
HartreeFock::HartreeFock(Wavefunction *wf, Method method, double x_Breit,
                         double eps)
    : HartreeFock(wf->rgrid, wf->vnuc, &wf->core,
                  wf->qed ? *(wf->qed) : std::optional<QED::RadPot>{},
                  wf->alpha, method, x_Breit, eps) {}

//==============================================================================
const std::vector<double> &HartreeFock::solve_core() {

  // Core orbs must already be solutions... this OK?

  if (p_core->empty())
    return m_vdir;

  switch (m_method) {
  case Method::HartreeFock:
    hf_approx_core(1.0e-8);
    hf_core_refine();
    break;
  case Method::ApproxHF:
    hf_approx_core(m_eps_HF);
    break;
  case Method::Hartree:
    hf_approx_core(0.1 * m_eps_HF);
    break;
  case Method::KohnSham:
    KohnSham_core(0.1 * m_eps_HF);
    break;
  case Method::Local:
    break;
  }

  return m_vdir;
}

//==============================================================================
std::vector<double> HartreeFock::get_vlocal(int l) const {
  const auto &vrad_el = get_Hrad_el(l);
  return qip::add(m_vnuc, m_vdir, vrad_el);
}

//==============================================================================
DiracSpinor HartreeFock::VBr(const DiracSpinor &Fv) const {
  if (m_VBr)
    return m_VBr->VbrFa(Fv, *p_core);
  else
    return 0.0 * Fv;
}

//==============================================================================
// Solves "approximate" (localised) HF equations for the core
void HartreeFock::hf_approx_core(const double eps_target_HF) {
  if (p_core->empty()) {
    return;
  }

  using namespace qip::overloads;

  const auto eta_damp = 0.55;

  // Start the HF itterative procedure:
  int hits = 1;
  double eps = 1.0;
  std::vector<std::vector<double>> vex_core(p_core->size());
  for (; hits < m_max_hf_its; hits++) {

    // Store old vdir/vex
    const auto vdir_old = m_vdir;
    const auto vex_old = vex_core;

    // Form new v_dir and v_ex:
    m_Yab.calculate(*p_core);
    form_vdir(&m_vdir);
    vex_core = form_approx_vex_core();

    // damp vdir and vex
    m_vdir = (1.0 - eta_damp) * m_vdir + eta_damp * vdir_old;
    for (std::size_t i = 0; i < p_core->size(); i++) {
      vex_core[i] = (1.0 - eta_damp) * vex_core[i] + eta_damp * vex_old[i];
    }

    // Solve Dirac Eq. for each state in core, using Vdir+Vex:
    double t_eps = 0.0;
    for (std::size_t i = 0; i < p_core->size(); i++) {
      auto &Fa = (*p_core)[i];
      const double en_old = Fa.en();
      // Guess energy correction [speeds up boundState()]
      const double dEa =
          Fa * ((m_vdir - vdir_old + vex_core[i] - vex_old[i]) * Fa);
      const double en_guess = (en_old < -dEa) ? en_old + dEa : en_old;
      // Solve Dirac equation:
      const auto Hmag = get_Hrad_mag(Fa.l());
      const auto v = get_vlocal(Fa.l()) + vex_core[i];
      DiracODE::boundState(Fa, en_guess, v, Hmag, m_alpha, 7);
      // Find largest epsilon:
      const double state_eps = std::abs((Fa.en() - en_old) / en_old);
      t_eps = std::max(state_eps, t_eps);
    }
    eps = t_eps;

    if (eps < eps_target_HF)
      break;
  }

  if (verbose && m_method != Method::HartreeFock)
    printf("HF(aprx) core      it:%3i eps=%6.1e\n", hits, eps);

  // Now, re-solve core orbitals with higher precission
  for (std::size_t i = 0; i < p_core->size(); i++) {
    auto &Fa = (*p_core)[i];
    const auto &vrad_el = get_Hrad_el(Fa.l());
    const auto &vrad_mag = get_Hrad_mag(Fa.l());
    const auto v = qip::add(m_vnuc, m_vdir, vrad_el, vex_core[i]);
    DiracODE::boundState(Fa, Fa.en(), v, vrad_mag, m_alpha, 14);
  }
}

//==============================================================================
void HartreeFock::KohnSham_core(const double eps_target_HF) {

  if (p_core->empty()) {
    return;
  }

  const auto eta_damp = 0.5;

  // Start the itterative procedure:
  int hits = 1;
  double t_eps = 1.0;
  for (; hits < m_max_hf_its; hits++) {

    // Store old vdir/vex
    const auto vdir_old = m_vdir;

    using namespace qip::overloads;

    // Form new v_dir, and damp:
    m_Yab.calculate(*p_core);
    form_vdir(&m_vdir);
    KohnSham_addition(&m_vdir);
    m_vdir = (1.0 - eta_damp) * m_vdir + eta_damp * vdir_old;

    // Solve Dirac Eq. for each state in core, using Vdir (incl KS):
    t_eps = 0;
    for (std::size_t i = 0; i < p_core->size(); i++) {
      auto &Fa = (*p_core)[i];

      const double en_old = Fa.en();
      const auto dEa = Fa * ((m_vdir - vdir_old) * Fa);
      const double en_guess = (en_old < -dEa) ? en_old + dEa : en_old;

      const auto vrad_mag = get_Hrad_mag(Fa.l());
      const auto v = get_vlocal(Fa.l());
      DiracODE::boundState(Fa, en_guess, v, vrad_mag, m_alpha, 7);
      const double state_eps = std::abs((Fa.en() - en_old) / en_old);
      // convergance based on worst orbital:
      t_eps = (state_eps > t_eps) ? state_eps : t_eps;

    } // core states

    if (t_eps < eps_target_HF)
      break;
  } // hits

  if (verbose) {
    printf("KS core      it:%3i eps=%6.1e\n", hits, t_eps);
  }

  // Now, re-solve core orbitals with higher precission
  for (auto &Fa : *p_core) {
    const auto &vrad_mag = get_Hrad_mag(Fa.l());
    const auto v = get_vlocal(Fa.l());
    DiracODE::boundState(Fa, Fa.en(), v, vrad_mag, m_alpha, 15);
  }
}
//==============================================================================
void HartreeFock::KohnSham_addition(std::vector<double> *vdir) const {

  const auto f =
      -(2.0 / 3.0) * std::pow(81.0 / (32.0 * M_PI * M_PI), 1.0 / 3.0);

  std::vector<double> rho(rgrid->num_points());
  for (const auto &Fc : *p_core) {
    rho = qip::add(rho, Fc.rho());
  }

  for (std::size_t i = 0; i < rgrid->num_points(); ++i) {
    const auto r = rgrid->r(i);
    vdir->at(i) += (f / r) * std::pow(r * rho[i], 1.0 / 3.0);
  }

  // KS is a V^N (not V^N-1) approximation
  const auto z_ion = zion() + 1;
  // Latter correction:
  // Enforce V(r) = Vnuc(r)+Vel(r) =~ -1/r at large r
  for (std::size_t i = rgrid->num_points() - 1; i != 0; --i) {
    // nb: miss i=0, but fine. Only applies large r[i]
    const auto vn = (m_vnuc)[i];
    const auto r = rgrid->r(i);
    if (r * std::abs(vn + vdir->at(i)) > z_ion)
      break;
    vdir->at(i) = -z_ion / r - vn;
  }
}

//==============================================================================
void HartreeFock::solveValence(std::vector<DiracSpinor> *valence,
                               const bool print) {

  if (valence == nullptr || valence->empty())
    return;

  const auto Nval = valence->size();
  auto do_refine = (m_method == Method::HartreeFock && !p_core->empty());

  std::vector<EpsIts> eis(Nval);

#pragma omp parallel for
  for (std::size_t i = 0; i < Nval; i++) {
    auto &Fa = (*valence)[i];
    if (do_refine)
      eis[i] = hf_valence(Fa);
    else
      eis[i] = hf_valence_approx(Fa, m_eps_HF);
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
//==============================================================================
void HartreeFock::solveBrueckner(std::vector<DiracSpinor> *valence,
                                 const MBPT::CorrelationPotential &Sigma2,
                                 const bool print) {

  if (valence->empty())
    return;

  if (print) {
    std::cout << "Solving for Brueckner orbitals (correlation potential)\n";
    Sigma2.print_scaling();
  }
  for (auto &Fv : *valence) {
    const auto en_old = Fv.en();
    if (print)
      std::cout << Fv.symbol() << ":" << std::flush;
    const auto eis = hf_valence(Fv, &Sigma2);
    const auto delta = Fv.en() - en_old;
    if (print) {
      printf(" delta=%8.5f; eps=%6.1e [its=%3i]", delta, eis.eps, eis.its);
      if (eis.eps > m_eps_HF && eis.eps > 1.0e-6) {
        std::cout << " ====*";
      }
      std::cout << "\n";
    }
  }
}

//==============================================================================
EpsIts HartreeFock::hf_valence_approx(DiracSpinor &Fa, double eps_target_HF)
// Solves HF for given orbital Fa, in frozen core.
// Does not store vex (must be done outside)
// Can be used to generate a set of virtual/basis orbitals
{
  Fa.occ_frac() = 1.0 / Fa.twojp1();

  auto damper = rampedDamp(0.7, 0.3, 2, 6);
  // don't include all pts in PT for new e guess
  const std::size_t de_stride = 5;

  std::vector<double> vexa(rgrid->num_points(), 0);
  auto vexa_old = vexa;

  const auto &vrad_el = get_Hrad_el(Fa.l());
  const auto &Hmag = get_Hrad_mag(Fa.l());
  const auto v_local = qip::add(m_vnuc, m_vdir, vrad_el);

  double eps = -1, eps_prev = -1;
  int hits = 1;
  for (; hits < m_max_hf_its; hits++) {
    auto eta = damper(hits);

    double en_old = Fa.en();
    vexa_old = vexa;
    if (!excludeExchangeQ())
      vexa = vex_approx(Fa, *p_core);

    for (std::size_t i = 0; i < rgrid->num_points(); i++) {
      vexa[i] = (1.0 - eta) * vexa[i] + eta * vexa_old[i];
    }
    // Use P.T. to calculate energy change:
    double en_new_guess = 0;
    for (std::size_t i = 0; i < Fa.max_pt(); i += de_stride) {
      en_new_guess +=
          (vexa[i] - vexa_old[i]) * Fa.f(i) * Fa.f(i) * rgrid->drdu(i);
    }
    en_new_guess = en_old + en_new_guess * rgrid->du() * de_stride;
    // Solve Dirac using new potential:
    DiracODE::boundState(Fa, en_new_guess, qip::add(v_local, vexa), Hmag,
                         m_alpha, 15);
    eps = std::abs((Fa.en() - en_old) / en_old);

    auto getting_worse = (hits > 20 && eps >= eps_prev && eps < 1.e-5);
    auto converged = (eps <= eps_target_HF);
    if (converged || getting_worse)
      break;
    eps_prev = eps;
  }

  return {eps, hits};
}

//==============================================================================
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
    const double xtjap1 = (tja + 1) * Fa.occ_frac();
    e1 += xtjap1 * Fa.en();
    for (const auto &Fb : *p_core) {
      const auto tjb = Fb.twoj();
      const double xtjbp1 = (tjb + 1) * Fb.occ_frac();
      const auto &v0bb = *m_Yab.get(0, Fb, Fb); // XXX may be null
      const auto R0fg2 = Fa * (v0bb * Fa);
      e2 += xtjap1 * xtjbp1 * R0fg2;
      // take advantage of symmetry for third term:
      if (Fb > Fa)
        continue;
      const auto y = (Fa == Fb) ? 1.0 : 2.0; // symmetry
      const int kmin = std::abs(tja - tjb) / 2;
      const int kmax = (tja + tjb) / 2;
      // const auto &vabk = m_Yab.get_y_ab(Fa, Fb);
      for (auto k = kmin; k <= kmax; k++) {
        const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.kappa(), Fb.kappa());
        if (Labk == 0)
          continue;
        const auto vabk = m_Yab.get(k, Fa, Fb);
        if (vabk == nullptr)
          continue;
        const auto R0fg3 = Fa * (*vabk * Fb);
        e3 += y * xtjap1 * xtjbp1 * Labk * (R0fg3);
      }
    }
    Etot += e1 - 0.5 * (e2 - e3);
  }
  return Etot;
}

//==============================================================================
int HartreeFock::num_core_electrons() const {
  return std::accumulate(
      p_core->cbegin(), p_core->cend(), 0,
      [](int n, const auto &Fa) { return n + Fa.num_electrons(); });
}

//==============================================================================
void HartreeFock::form_vdir(std::vector<double> *vdir, ReScale re_scale) const
// Forms the direct part of the potential.
// Must call either form_vbb0 or form_vabk_core first!
// Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
// If re_scale==true, will scale by (N-1)/N. This then given the averaged
// Hartree potential (local, same each state, no exchange).
// re_scale=false by default
{
  using namespace qip::overloads;
  (*vdir) *= 0.0;

  // Rescale so that V_nuc + Vdir = -1/r at large r
  // (Only done in initial approx., before exchange)
  const double scale =
      re_scale == ReScale::yes ? (1.0 - 1.0 / num_core_electrons()) : 1.0;

  for (const auto &Fb : *p_core) {
    const double f_sf = scale * Fb.twojp1() * Fb.occ_frac();
    const auto v0bb = m_Yab.get(0, Fb, Fb);
    assert(v0bb != nullptr && "Missing Fb from m_Yab in HF::form_vdir?");
    (*vdir) += f_sf * (*v0bb);
  }
}

//==============================================================================
const std::vector<double> &
HartreeFock::get_vdir_single(const DiracSpinor &Fa) const {
  return *m_Yab.get(0, Fa, Fa);
}
//----------
std::vector<double> HartreeFock::calc_vdir_single(const DiracSpinor &Fa) {
  return Coulomb::yk_ab(Fa, Fa, 0);
}

//==============================================================================
void HartreeFock::form_approx_vex_core(
    std::vector<std::vector<double>> &vex) const
// Forms the 2D "approximate" exchange potential for each core state, a.
// NOTE: Must call form_vabk_core first!
// Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
{
  vex.resize(p_core->size());
#pragma omp parallel for
  for (std::size_t a = 0; a < p_core->size(); a++) {
    form_approx_vex_core_a((*p_core)[a], vex[a]);
  }
}

std::vector<std::vector<double>> HartreeFock::form_approx_vex_core() const {
  std::vector<std::vector<double>> vex_core(p_core->size());
  form_approx_vex_core(vex_core);
  return vex_core;
}

std::vector<double>
HartreeFock::form_approx_vex_core_a(const DiracSpinor &Fa) const {
  std::vector<double> vex_a;
  form_approx_vex_core_a(Fa, vex_a);
  return vex_a;
}

//==============================================================================
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
  vex_a.clear();
  vex_a.resize(rgrid->num_points());

  // don't subtract self-potential term for Hartree {match Core-Hartree}
  if (excludeExchangeQ())
    return;

  const auto twoj_a = Fa.twoj();

  const auto max_abs = [](double a, double b) {
    return (std::abs(a) < std::abs(b));
  };
  const auto max =
      std::abs(*std::max_element(Fa.f().begin(), Fa.f().end(), max_abs));
  const auto cut_off = 0.003 * max;

  if (!excludeExchangeQ()) {
    for (const auto &Fb : (*p_core)) {
      // b!=a
      if (Fb == Fa)
        continue;
      const auto tjb = Fb.twoj();
      const double x_tjbp1 = (tjb + 1) * Fb.occ_frac();
      const auto irmax = std::min(Fa.max_pt(), Fb.max_pt());
      const int kmin = std::abs(twoj_a - tjb) / 2;
      const int kmax = (twoj_a + tjb) / 2;
      // const auto &vabk = m_Yab.get_y_ab(Fb, Fa);

      // hold "fraction" Fa*Fb/(Fa^2):
      std::vector<double> v_Fab(rgrid->num_points());
      for (std::size_t i = 0; i < irmax; i++) {
        // This is the approximte part! Divides by Fa
        if (std::abs(Fa.f(i)) < cut_off)
          continue;
        const auto fac_top = Fa.f(i) * Fb.f(i) + Fa.g(i) * Fb.g(i);
        const auto fac_bot = Fa.f(i) * Fa.f(i) + Fa.g(i) * Fa.g(i);
        v_Fab[i] = -x_tjbp1 * fac_top / fac_bot;
      } // r
      for (int k = kmin; k <= kmax; k++) {
        const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.kappa(), Fb.kappa());
        if (Labk == 0)
          continue;
        // const auto ik = std::size_t(k - kmin);
        const auto vabk = m_Yab.get(k, Fb, Fa);
        if (vabk == nullptr)
          continue;
        for (std::size_t i = 0; i < irmax; i++) {
          vex_a[i] += Labk * (*vabk)[i] * v_Fab[i];
        } // r
      }   // k
    }     // b
  }

  // now, do a=b, ONLY if a is in the core!
  // if (a_in_coreQ)
  {
    const double x_tjap1 = (twoj_a + 1); // no occ_frac here ?
    const int kmax = twoj_a;
    // const auto &vaak = m_Yab.get_y_ab(Fa, Fa);
    const auto irmax = Fa.max_pt(); // rgrid->num_points();
    for (int k = 0; k <= kmax; k++) {
      const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.kappa(), Fa.kappa());
      if (Labk == 0)
        continue;
      const auto vaak = m_Yab.get(k, Fa, Fa);
      if (vaak == nullptr)
        continue;
      for (std::size_t i = 0; i < irmax; i++) {
        // nb: If I don't 'cut' here, or fails w/ f states... ?? XX
        // Of course, cutting is fine. But WHY FAIL??
        vex_a[i] += -Labk * (*vaak)[i] * x_tjap1;
      }
    } // k
  }   // if a in core
}

//------------------------------------------------------------------------------
std::vector<double> vex_approx(const DiracSpinor &Fa,
                               const std::vector<DiracSpinor> &core, int k_cut,
                               const double lambda_cut)
// Free function
{

  std::vector<double> vex(Fa.grid().num_points());
  std::vector<double> vabk;

  const auto tja = Fa.twoj();
  const auto la = Fa.l();

  const auto max_abs = [](double a, double b) {
    return (std::abs(a) < std::abs(b));
  };
  const auto max =
      std::abs(*std::max_element(Fa.f().begin(), Fa.f().end(), max_abs));
  const auto cut_off = lambda_cut * max;

  for (const auto &Fb : core) {
    const auto tjb = Fb.twoj();
    const auto lb = Fb.l();
    const double x_tjbp1 = (tjb + 1) * Fb.occ_frac(); // when in core??
    const auto irmax = std::min(Fa.max_pt(), Fb.max_pt());
    const int kmin = std::abs(tja - tjb) / 2;
    int kmax = (tja + tjb) / 2;
    if (kmax > k_cut)
      kmax = k_cut;

    // hold "fraction" Fa*Fb/(Fa^2):
    std::vector<double> v_Fab(Fa.grid().num_points());
    for (std::size_t i = 0; i < irmax; i++) {
      // This is the approximate part! Divides by Fa
      if (std::abs(Fa.f(i)) < cut_off)
        continue;
      const auto fac_top = Fa.f(i) * Fb.f(i) + Fa.g(i) * Fb.g(i);
      const auto fac_bot = Fa.f(i) * Fa.f(i) + Fa.g(i) * Fa.g(i);
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

//==============================================================================
DiracSpinor HartreeFock::calc_vexFa_core(const DiracSpinor &Fa) const
// calculates V_ex Fa (returns new Dirac Spinor)
// Fa can be any orbital (so long as coulomb integrals exist!)
// ..... i.e., must be core orbital
{
  DiracSpinor VxFa(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  vex_psia_core(Fa, VxFa);
  return VxFa;
}
//------------------------------------------------------------------------------
void HartreeFock::vex_psia_core(const DiracSpinor &Fa, DiracSpinor &VxFa) const
// calculates V_ex Fa
// Fa can be any orbital (so long as coulomb integrals exist!)
// ..... i.e., must be core orbital
{
  VxFa.max_pt() = Fa.f().size(); // silly hack. Make sure VxFa = 0 after pinf
  VxFa *= 0.0;
  VxFa.max_pt() = Fa.max_pt(); // updated below

  if (excludeExchangeQ())
    return;

  const auto twoj_a = Fa.twoj();
  for (const auto &Fb : (*p_core)) {
    VxFa.max_pt() = std::max(VxFa.max_pt(), Fb.max_pt());
    const auto tjb = Fb.twoj();
    const double x_tjbp1 = (Fa == Fb) ? (tjb + 1) : (tjb + 1) * Fb.occ_frac();
    const int kmin = std::abs(twoj_a - tjb) / 2;
    const int kmax = (twoj_a + tjb) / 2;
    // const auto &vabk = m_Yab.get_y_ab(Fb, Fa);
    for (int k = kmin; k <= kmax; k++) {
      const auto Labk = m_Yab.Ck().get_Lambdakab(k, Fa.kappa(), Fb.kappa());
      if (Labk == 0)
        continue;
      const auto vabk = m_Yab.get(k, Fb, Fa);
      if (vabk == nullptr)
        continue;
      for (auto i = 0u; i < Fb.max_pt(); i++) {
        const auto v = -x_tjbp1 * Labk * (*vabk)[i];
        VxFa.f(i) += v * Fb.f(i);
        VxFa.g(i) += v * Fb.g(i);
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
  DiracSpinor VxFa(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  VxFa.max_pt() = Fa.max_pt(); // nb: updated below!
  // note: VxFa.max_pt() can be larger than psi.max_pt()!

  std::vector<double> vabk(Fa.grid().num_points());

  const auto tja = Fa.twoj();
  const auto la = Fa.l();
  for (const auto &Fb : core) {
    VxFa.max_pt() = std::max(Fb.max_pt(), VxFa.max_pt());
    const auto tjb = Fb.twoj();
    const auto lb = Fb.l();
    const double x_tjbp1 = (Fa == Fb) ? (tjb + 1) : (tjb + 1) * Fb.occ_frac();
    const int kmin = std::abs(tja - tjb) / 2;
    const int kmax = std::min((tja + tjb) / 2, k_cut);
    for (int k = kmin; k <= kmax; k++) {
      if (Angular::parity(la, lb, k) == 0)
        continue;
      const auto tjs = Angular::threej_2(tjb, tja, 2 * k, -1, 1, 0);
      if (tjs == 0)
        continue;
      Coulomb::yk_ab(Fb, Fa, k, vabk, Fb.max_pt());
      const auto factor = -x_tjbp1 * tjs * tjs;
      for (auto i = 0u; i < Fb.max_pt(); i++) {
        VxFa.f(i) += factor * vabk[i] * Fb.f(i);
        VxFa.g(i) += factor * vabk[i] * Fb.g(i);
      } // r
    }   // k
  }     // b

  return VxFa;
}

//==============================================================================
std::vector<double> HartreeFock::get_Hrad_el(int l) const {
  return m_vrad ? m_vrad->Vel(l) : std::vector<double>{};
}
std::vector<double> HartreeFock::get_Hrad_mag(int l) const {
  return m_vrad ? m_vrad->Hmag(l) : std::vector<double>{};
}
//==============================================================================
//==============================================================================
//==============================================================================

//==============================================================================
EpsIts HartreeFock::hf_valence(DiracSpinor &Fa,
                               const MBPT::CorrelationPotential *const Sigma) {
  if (p_core->empty())
    return {0, 0};

  const auto eps_target = 0.001 * m_eps_HF;
  const auto damper = 0.4;

  const auto Hmag = get_Hrad_mag(Fa.l());
  const auto vl = get_vlocal(Fa.l());

  auto prev_en = Fa.en();
  int it = 0;
  double eps = 1.0;
  for (; it <= m_max_hf_its; ++it) {
    const auto a_damp = damper;

    auto VxFa = calc_vexFa(Fa);
    if (m_VBr) {
      VxFa += m_VBr->VbrFa(Fa, *p_core);
    }
    if (Sigma) {
      VxFa += (*Sigma)(Fa);
    }
    const auto Fa_prev = Fa;
    DiracODE::boundState(Fa, Fa.en(), vl, Hmag, m_alpha, 16, &VxFa, &Fa_prev,
                         zion());
    eps = std::abs((prev_en - Fa.en()) / Fa.en());
    prev_en = Fa.en();

    const bool converged = (eps <= eps_target && it > 0);
    if (converged || it == m_max_hf_its)
      break;

    Fa = (1.0 - a_damp) * Fa + a_damp * Fa_prev;
    Fa.normalise();

  } // End HF its

  return {eps, it};
}

//==============================================================================
inline void HartreeFock::hf_core_refine() {
  if (p_core->empty()) {
    return;
  }

  const double eps_target = m_eps_HF;
  // m_Yab.update_y_ints(); // only needed if not already done!
  m_Yab.calculate(*p_core);
  auto damper = 0.35; // rampedDamp(0.35, 0.35, 5, 30);
  // double extra_damp = 0;

  std::vector<double> vl(rgrid->num_points()); // Vnuc + fVd
  std::vector<double> v0(rgrid->num_points()); // (1-f)Vd
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
    vexF_list.push_back(DiracSpinor(Fa.n(), Fa.kappa(), Fa.grid_sptr()));
  }

  double eps = 0.0;
  double best_eps = 1.0;
  double best_worst_eps = 1.0;
  std::size_t worst_index = 0;
  std::size_t best_index = 0;
  // int worse_count = 0;
  int it = 0;
  for (; it <= m_max_hf_its; it++) {
    const auto a_damp = damper; //(it);// + extra_damp;

    // re-calculate each Vl = vnuc + fvdir, v0 = (1-f)vdir:
    for (auto i = 0ul; i < rgrid->num_points(); i++) {
      vl[i] = (m_vnuc)[i] + f_core * vd[i];
      v0[i] = (1.0 - f_core) * vd[i];
    }

    // re-calculate each VexPsi:
    for (std::size_t i = 0; i < num_core_states; ++i) {
      const auto &Fa = (*p_core)[i];
      vex_psia_core(Fa, vexF_list[i]);
    }

    core_prev = (*p_core);

#pragma omp parallel for
    for (std::size_t i = 0; i < num_core_states; ++i) {
      auto &Fa = (*p_core)[i];
      const auto &Fzero = core_zero[i];
      const auto &vexFzero = vexCore_zero[i];

      const auto Fa_prev = core_prev[i];
      const auto &VxFa = vexF_list[i];

      auto en = Fzero.en() + (Fzero * VxFa - Fa * vexFzero + Fzero * (vd * Fa) -
                              Fa * (vd0 * Fzero)) /
                                 (Fa * Fzero);
      auto v_nonlocal = v0 * Fa + VxFa;
      if (m_VBr) {
        const auto VbrFa =
            m_VBr->VbrFa(Fa, core_prev); // depends on previous core!
        en += (Fzero * VbrFa) / (Fa * Fzero);
        v_nonlocal += VbrFa;
      }

      const auto &Hrad_el = get_Hrad_el(Fa.l());
      const auto &Hmag = get_Hrad_mag(Fa.l());
      const auto &VlVr = qip::add(vl, Hrad_el);
      // hf_orbital(Fa, en, VlVr, Hmag, v_nonlocal, core_prev, v0, VBr.get());
      hf_orbital_green(Fa, en, VlVr, Hmag, v_nonlocal, core_prev, v0,
                       get_Breit());
      Fa = (1.0 - a_damp) * Fa + a_damp * Fa_prev;
      Fa.normalise();
      auto d_eps = std::abs((Fa_prev.en() - Fa.en()) / Fa.en());
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

    const bool converged = (eps <= eps_target && it > 0);
    if (converged || it == m_max_hf_its)
      break;

    if (eps < best_worst_eps)
      best_worst_eps = eps;

    m_Yab.calculate(*p_core);
    form_vdir(&m_vdir);
  }

  if (verbose)
    printf("HF core:  it:%3i eps=%6.1e for %s  [%6.1e for %s]\n", //
           it, eps, (*p_core)[worst_index].symbol().c_str(), best_eps,
           (*p_core)[best_index].symbol().c_str());
}

//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================

//==============================================================================
void HartreeFock::hf_orbital_green(
    DiracSpinor &Fa, double en0, const std::vector<double> &vl,
    const std::vector<double> &H_mag, const DiracSpinor &VnlFa0,
    const std::vector<DiracSpinor> &static_core, const std::vector<double> &dv0,
    const HF::Breit *const tVBr,
    const MBPT::CorrelationPotential *const Sigma) const
// Solve Dirac Equation (Eigenvalue):
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
  // should work for core, valence, w/wo (Br, Sigma, H_mag)

  // pull these outside? But make sure thread safe!
  DiracSpinor Gzero(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  DiracSpinor Ginf(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  DiracSpinor dFa(Fa.n(), Fa.kappa(), Fa.grid_sptr());
  const auto eps_target = 1.0e-18; // m_eps_HF;
  const auto k_max = 1;            // max k for Vex into dEa

  // const auto SigmaF = (*Sigma)(Fa);

  // nb: VnlFa0 includes Breit + Sigma etc. VBr used for energy correction

  const auto alpha = m_alpha;
  DiracODE::solve_inhomog(Fa, Gzero, Ginf, en0, vl, H_mag, alpha,
                          -1.0 * VnlFa0);

  // make small adjustments to energy to normalise Fa:
  DiracODE::Adams::GreenSolution(dFa, Ginf, Gzero, alpha, Fa);
  // should dFa = dEa * dFa, but makes it worse?
  // nb: after first it, becomes correct.
  double en = en0;
  double dEa = 0.5 * (Fa * Fa - 1.0) / (Fa * dFa);
  double eps = std::abs(dEa / en0);
  int tries = 0;
  for (; tries <= m_max_hf_its; ++tries) {
    if (eps < eps_target && tries > 1)
      break;
    {
      auto VnlF_tilde = vexFa(dFa, static_core, k_max);
      if (tVBr)
        VnlF_tilde += tVBr->VbrFa(dFa, static_core);
      if (Sigma)
        VnlF_tilde += (*Sigma)(dFa);
      if (!dv0.empty())
        VnlF_tilde += dv0 * dFa;
      // const auto SigmaF_tilde = Sigma(dFa);
      DiracODE::Adams::GreenSolution(dFa, Ginf, Gzero, alpha,
                                     dEa * Fa - VnlF_tilde);
    }
    const auto delta_Norm = Fa * Fa - 1.0;
    const auto de0 = dEa;
    dEa *= 0.5 * (delta_Norm) / (Fa * dFa);
    eps = std::abs(dEa / en);
    en += dEa;
    dFa.max_pt() = Fa.max_pt();
    Fa -= (dEa / de0) * dFa;
  }
  Fa.en() = en;
  Fa.eps() = eps;
  Fa.its() = tries;
  Fa.normalise();
}

//==============================================================================
double HartreeFock::zion() const {
  const auto z = -1.0 * m_vnuc.back() * rgrid->r().back();
  return z - double(num_core_electrons());
}

} // namespace HF
