#include "HartreeFockClass.hpp"
#include "ATI_atomInfo.hpp"
#include "DiracSpinor.hpp"
#include "ElectronOrbitals.hpp"
#include "Grid.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include "PRM_localPotentials.hpp"
#include "Wigner_369j.hpp"
#include <cmath>
#include <vector>
/*
Calculates self-consistent Hartree-Fock potential, including exchange.
Solves all core and valence states.
*/

#define DO_DEBUG false
#if DO_DEBUG
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif // DEBUG

//******************************************************************************
HartreeFock::HartreeFock(ElectronOrbitals &wf, const std::string &in_core,
                         double eps_HF, bool in_ExcludeExchange)
    : p_wf(&wf), p_rgrid(&wf.rgrid), m_excludeExchange(in_ExcludeExchange) {

  m_eps_HF = eps_HF;
  if (eps_HF > 1)
    m_eps_HF = pow(10, -1 * eps_HF); // can give as log..

  starting_approx_core(in_core);
  m_num_core_states = p_wf->coreIndexList.size();

  // store l, 2j, and "kappa_index" in arrays for faster/easier access
  twoj_list.reserve(m_num_core_states);
  kappa_index_list.reserve(m_num_core_states);
  for (auto &psi : p_wf->orbitals) {
    twoj_list.push_back(psi.twoj());
    kappa_index_list.push_back(ATI::indexFromKappa(psi.k));
  }

  // Run HF for all core states
  hartree_fock_core();
}

//******************************************************************************
void HartreeFock::hartree_fock_core() {

  static const double eta1 = 0.35;
  static const double eta2 = 0.7; // this value after 4 its
  // don't include all pts in PT for new e guess:
  static const std::size_t de_stride = 5;

  form_core_Lambda_abk();

  vex.resize(m_num_core_states, std::vector<double>(p_rgrid->ngp)); // Also
  initialise_m_arr_v_abk_r_core();

  std::vector<double> vdir_old(p_rgrid->ngp);
  std::vector<std::vector<double>> vex_old(m_num_core_states,
                                           std::vector<double>(p_rgrid->ngp));

  // Start the HF itterative procedure:
  int hits = 1;
  double t_eps;
  double eta = eta1;
  for (; hits < MAX_HART_ITS; hits++) {
    DEBUG(std::cerr << "HF core it: " << hits << "\n";)
    if (hits == 4)
      eta = eta2;
    if (hits == 16)
      eta = 0.5 * (eta1 + eta2);
    if (hits == 32)
      eta = eta1;

    // Form new v_dir and v_ex:
    vdir_old = p_wf->vdir;
    vex_old = vex;

    form_vabk_core();
    form_vdir(p_wf->vdir, false);
    form_approx_vex_core(vex);

    for (std::size_t j = 0; j < p_rgrid->ngp; j++) {
      p_wf->vdir[j] = eta * p_wf->vdir[j] + (1. - eta) * vdir_old[j];
      for (std::size_t i = 0; i < m_num_core_states; i++) {
        vex[i][j] = eta * vex[i][j] + (1. - eta) * vex_old[i][j];
      }
    }

    // Solve Dirac Eq. for each state in core, using Vdir+Vex:
    t_eps = 0;
    for (std::size_t i = 0; i < m_num_core_states; i++) {
      auto &phi = p_wf->orbitals[i];
      double en_old = phi.en;
      // calculate de from PT
      double del_e = 0;
      for (std::size_t j = 0; j < phi.pinf; j += de_stride) {
        double dv = (p_wf->vdir[j] - vdir_old[j]) + (vex[i][j] - vex_old[i][j]);
        del_e += dv * phi.f[j] * phi.f[j] * p_rgrid->drdu[j];
      }
      del_e *= p_rgrid->du * de_stride;
      double en_guess = (en_old < -del_e) ? en_old + del_e : en_old;
      p_wf->reSolveDirac(phi, en_guess, vex[i], 3);
      double state_eps = fabs((phi.en - en_old) / en_old);
      // convergance based on worst orbital:
      DEBUG(printf(" --- %2i,%2i: en=%11.5f  HFeps = %.0e;  Adams = %.0e[%2i]  "
                   "(%4i)\n",
                   phi.n, phi.k, phi.en, state_eps, phi.eps, phi.its,
                   (int)phi.pinf);)
      if (state_eps > t_eps)
        t_eps = state_eps;
    } // core states
    DEBUG(std::cerr << "HF core it: " << hits << ": eps=" << t_eps << "\n\n";
          std::cin.get();)

    // Force all core orbitals to be orthogonal to each other
    p_wf->orthonormaliseOrbitals(1);
    if (t_eps < m_eps_HF)
      break;
  } // hits
  printf("\rHF core        it:%3i eps=%6.1e              \n", hits, t_eps);

  // Now, re-solve core orbitals with higher precission
  for (std::size_t i = 0; i < m_num_core_states; i++) {
    p_wf->reSolveDirac(p_wf->orbitals[i], p_wf->orbitals[i].en, vex[i], 14);
  }
  p_wf->orthonormaliseOrbitals(2);
}

//******************************************************************************
void HartreeFock::solveValence(int n, int kappa) {

  // XXX This is slow! Much better to have another function that takes in
  // a _list_ of valence states, and reserves all space first!
  // Could also be a little clever, and do 'initial solve' outside,
  // then could run the HF part in parallel.
  // Note: probs not necisary: the 'form vex' part is // anyway, and fairly
  // efficient
  twoj_list.push_back(ATI::twoj_k(kappa));
  kappa_index_list.push_back(ATI::indexFromKappa(kappa));
  extend_Lambda_abk(kappa);
  extend_m_arr_v_abk_r_valence(kappa);

  // just use direct to solve initial
  p_wf->solveLocalDirac(n, kappa, 0, 3);
  auto a = p_wf->valenceIndexList.back();
  int twoJplus1 = ATI::twoj_k(kappa) + 1;
  p_wf->orbitals.back().occ_frac = 1. / twoJplus1;

  auto &phi = p_wf->orbitals.back();

  static const double eta1 = 0.35;
  static const double eta2 = 0.7; // this value after 4 its
  // don't include all pts in PT for new e guess
  static const std::size_t de_stride = 5;

  std::vector<double> vexa(p_rgrid->ngp);
  std::vector<double> vexa_old;

  int hits = 1;
  double eps = -1;
  double eta = eta1;
  for (; hits < MAX_HART_ITS; hits++) {
    if (hits == 4)
      eta = eta2;

    double en_old = phi.en;
    vexa_old = vexa;

    form_vabk_valence(a);
    form_approx_vex_a(a, vexa);
    for (std::size_t i = 0; i < p_rgrid->ngp; i++) {
      vexa[i] = eta * vexa[i] + (1. - eta) * vexa_old[i];
    }
    // Use P.T. to calculate energy change:
    double en_new_guess = 0;
    for (std::size_t i = 0; i < phi.pinf; i += de_stride) {
      en_new_guess +=
          (vexa[i] - vexa_old[i]) * phi.f[i] * phi.f[i] * p_rgrid->drdu[i];
    }
    en_new_guess = en_old + en_new_guess * p_rgrid->du * de_stride;
    // Solve Dirac using new potential:
    p_wf->reSolveDirac(phi, en_new_guess, vexa, 3);
    eps = fabs((phi.en - en_old) / en_old);
    // Force valence states to be orthogonal to core:
    p_wf->orthonormaliseValence(phi, 1, true);
    if (eps < m_eps_HF)
      break;
  }
  printf("\rHF val:%2lu %2i %2i | %3i eps=%6.1e  en=%11.8f\n", a, n, kappa,
         hits, eps, phi.en);

  // Re-solve w/ higher precission
  p_wf->reSolveDirac(phi, phi.en, vexa, 15);
  p_wf->orthonormaliseValence(phi, 2);
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
  double Etot = 0;
  for (std::size_t a = 0; a < m_num_core_states; a++) {
    const auto &phi_a = p_wf->orbitals[a];

    double E1 = 0, E2 = 0, E3 = 0;
    double xtjap1 = (twoj_list[a] + 1) * phi_a.occ_frac;
    E1 += xtjap1 * phi_a.en;
    for (std::size_t b = 0; b < m_num_core_states; b++) {
      const auto &phi_b = p_wf->orbitals[b];
      double xtjbp1 = (twoj_list[b] + 1) * phi_b.occ_frac;
      auto irmax = std::min(phi_a.pinf, phi_b.pinf);
      auto &v0bb = get_v_aa0(b);
      double R0f2 = NumCalc::integrate(phi_a.f, phi_a.f, v0bb, p_rgrid->drdu, 1,
                                       0, irmax);
      double R0g2 = NumCalc::integrate(phi_a.g, phi_a.g, v0bb, p_rgrid->drdu, 1,
                                       0, irmax);
      E2 += xtjap1 * xtjbp1 * (R0f2 + R0g2);
      // take advantage of symmetry for third term:
      if (b > a)
        continue;
      double y = (a == b) ? 1 : 2;
      int kmin = abs(twoj_list[a] - twoj_list[b]) / 2;
      int kmax = (twoj_list[a] + twoj_list[b]) / 2;
      auto &vabk = get_v_abk(a, b);
      for (int k = kmin; k <= kmax; k++) {
        double L_abk = get_Lambda_abk(a, b, k);
        if (L_abk == 0)
          continue;
        int ik = k - kmin;
        double R0f3 =
            NumCalc::integrate(phi_a.f, phi_b.f, vabk[ik], p_rgrid->drdu);
        double R0g3 =
            NumCalc::integrate(phi_a.g, phi_b.g, vabk[ik], p_rgrid->drdu);
        E3 += y * xtjap1 * xtjbp1 * L_abk * (R0f3 + R0g3);
      }
    }
    {
      Etot += E1 - 0.5 * (E2 - E3) * p_rgrid->du; // update running total
    }
  }
  return Etot;
}

//******************************************************************************
void HartreeFock::starting_approx_core(const std::string &in_core)
// Starting approx for HF. Uses Green parametric
// Later, can put other options if you want.
{
  // make sure correct size (should already)
  p_wf->vdir.clear();
  p_wf->vdir.reserve(p_rgrid->ngp);
  int Z = p_wf->Znuc();
  // Get default values for Green potential [should be fitted to HF core orbs]
  double Gh, Gd; // Green potential parameters
  PRM::defaultGreenCore(Z, Gh, Gd);
  // Fill the the potential, using Greens PRM
  for (auto r : p_rgrid->r) {
    p_wf->vdir.emplace_back(PRM::green(Z, r, Gh, Gd));
  }
  p_wf->solveInitialCore(in_core, 3);
}

//******************************************************************************
void HartreeFock::form_core_Lambda_abk()
// Calculate + store the angular coefifienct:
//   Lambda^k_ab := 3js((ja,jb,k),(-1/2,1/2,0))^2*parity(la+lb+k)
// Lambda^k_ab = Lambda^k_ka,kb :: i.e., only depends on kappa!
// This routine re-sizes the m_arr_Lambda_nmk array
// New routine for valence? Or make so can re-call this one??
// XXX This should call the 'extend lambda' function !! XXX
{
  m_arr_Lambda_nmk.clear(); // should already be empty!

  // Find largest existing kappa index
  int max_kappa_index = 0;
  for (auto ki : kappa_index_list) {
    if (ki > max_kappa_index)
      max_kappa_index = ki;
  }
  m_max_kappa_index_so_far = max_kappa_index;

  m_arr_Lambda_nmk.reserve(max_kappa_index + 1);
  for (int n = 0; n <= max_kappa_index; n++) {
    int tja = ATI::twojFromIndex(n);
    int la = ATI::lFromIndex(n);
    std::vector<std::vector<double>> Lmk;
    Lmk.reserve(n + 1);
    for (int m = 0; m <= n; m++) {
      int tjb = ATI::twojFromIndex(m);
      int lb = ATI::lFromIndex(m);
      int kmin = (tja - tjb) / 2; // don't need abs, as m\leq n => ja\geq jb
      int kmax = (tja + tjb) / 2;
      std::vector<double> Lk(kmax - kmin + 1, 0);
      for (int k = kmin; k <= kmax; k++) {
        int ik = k - kmin;
        if (Wigner::parity(la, lb, k) == 0)
          continue;
        double tjs = Wigner::threej_2(tja, tjb, 2 * k, -1, 1, 0);
        Lk[ik] = tjs * tjs;
      } // k
      Lmk.push_back(Lk);
    } // m
    m_arr_Lambda_nmk.push_back(Lmk);
  } // n
}

//******************************************************************************
void HartreeFock::extend_Lambda_abk(int kappa_a)
// Note: there is code overlap with: form_core_Lambda_abk
// Could create a new function? Or just leave it?
// Or, merge this with above (if statements etc??)
// Note: don't just add this one, because might have skipped indexes!
// Add all we might need, keep order matchine index!
{
  int n_a = ATI::indexFromKappa(kappa_a);
  if (n_a <= m_max_kappa_index_so_far)
    return; // already done

  for (int n = m_max_kappa_index_so_far + 1; n <= n_a; n++) {
    int tja = ATI::twojFromIndex(n);
    int la = ATI::lFromIndex(n);
    std::vector<std::vector<double>> Lmk;
    Lmk.reserve(n + 1);
    for (int m = 0; m <= n; m++) {
      int tjb = ATI::twojFromIndex(m);
      int lb = ATI::lFromIndex(m);
      int kmin = (tja - tjb) / 2; // don't need abs, as m\leq n => ja\geq jb
      int kmax = (tja + tjb) / 2;
      std::vector<double> Lk(kmax - kmin + 1, 0);
      for (int k = kmin; k <= kmax; k++) {
        int ik = k - kmin;
        if (Wigner::parity(la, lb, k) == 0)
          continue;
        double tjs = Wigner::threej_2(tja, tjb, 2 * k, -1, 1, 0);
        Lk[ik] = tjs * tjs;
      }
      Lmk.push_back(Lk);
    }
    m_arr_Lambda_nmk.push_back(Lmk);
  }
  m_max_kappa_index_so_far = n_a;
}

//******************************************************************************
double HartreeFock::get_Lambda_abk(std::size_t a, std::size_t b, int k) const
// Simple routine to (semi-)safely return Lambda_abk
// Note: input a and b are regular ElectronOrbitals state indexes
// No checks. given k must be ok.
{
  // Get the kappa index for each state:
  int n = kappa_index_list[a];
  int m = kappa_index_list[b];
  int kmin = abs(twoj_list[a] - twoj_list[b]) / 2;

  return (n > m) ? m_arr_Lambda_nmk[n][m][k - kmin]
                 : m_arr_Lambda_nmk[m][n][k - kmin];
}

//******************************************************************************
void HartreeFock::initialise_m_arr_v_abk_r_core()
// Initialise (re-size) array to store CORE HF screening functions, v^k_ab(r)
// Note: only for core. These are stored in m_arr_v_abk_r array (class member)
{
  m_arr_v_abk_r.clear();
  m_arr_v_abk_r.resize(m_num_core_states); // reserve + push_back ?!
  for (std::size_t a = 0; a < m_num_core_states; a++) {
    m_arr_v_abk_r[a].resize(a + 1);
    std::size_t tja = (std::size_t)twoj_list[a];
    for (std::size_t b = 0; b <= a; b++) {
      std::size_t tjb = (std::size_t)twoj_list[b];
      std::size_t num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
      m_arr_v_abk_r[a][b].resize(num_k);
      for (std::size_t ik = 0; ik < num_k; ik++) {
        m_arr_v_abk_r[a][b][ik].resize(p_rgrid->ngp);
      } // k
    }   // b
  }     // a
}
//******************************************************************************
void HartreeFock::extend_m_arr_v_abk_r_valence(int kappa_a)
// This enlargens the m_arr_v_abk_r to make room for the valence states
{
  std::vector<std::vector<std::vector<double>>> v_abk_tmp(m_num_core_states);
  int tja = 2 * abs(kappa_a) - 1; // |2k|=2j+1
  for (std::size_t b = 0; b < m_num_core_states; b++) {
    int tjb = twoj_list[b];
    int num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
    v_abk_tmp[b].resize(num_k, std::vector<double>(p_rgrid->ngp));
  } // b
  m_arr_v_abk_r.push_back(v_abk_tmp);
}

//******************************************************************************
void HartreeFock::calculate_v_abk(const DiracSpinor &phi_a,
                                  const DiracSpinor &phi_b, int k,
                                  std::vector<double> &vabk) const
// Calculalates v^k_ab screening function.
// Note: should only call for a>=b, and for k's with non-zero angular coefs
// (nothing bad will happen othersie, but no point!)
// Since v_ab = v_ba
//
// Stores in vabk (reference to whatever) - must already be sized corectly!
//
// r_min := min(r,r')
// rho(r') := fa(r')*fb(r') + ga(r')gb(r')
// v^k_ab(r) = Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
//           = Int_0^r [r'^k/r^(k+1)]*rho(r') dr'
//             + Int_r^inf [r^k/r'^(k+1)]*rho(r') dr'
//          := A(r)/r^(k+1) + B(r)*r^k
// A(r0)  = 0
// B(r0)  = Int_0^inf [r^k/r'^(k+1)]*rho(r') dr'
// A(r_n) = A(r_{n-1}) + (rho(r_{n-1})*r_{n-1}^k)*dr
// B(r_n) = A(r_{n-1}) + (rho(r_{n-1})/r_{n-1}^(k+1))*dr
// v^k_ab(rn) = A(rn)/rn^(k+1) + B(rn)*rn^k
{
  auto irmax = std::min(phi_a.pinf, phi_b.pinf);

  double Ax = 0, Bx = 0;
  for (std::size_t i = 0; i < irmax; i++) {
    Bx += p_rgrid->drdu[i] *
          (phi_a.f[i] * phi_b.f[i] + phi_a.g[i] * phi_b.g[i]) /
          pow(p_rgrid->r[i], k + 1);
  }

  // For "direct" part, can't cut!
  if (phi_a == phi_b)
    irmax = (std::size_t)p_rgrid->ngp;

  vabk[0] = Bx * p_rgrid->du;
  for (std::size_t i = 1; i < irmax; i++) {
    double Fdr = p_rgrid->drdu[i - 1] * (phi_a.f[i - 1] * phi_b.f[i - 1] +
                                         phi_a.g[i - 1] * phi_b.g[i - 1]);
    Ax = Ax + Fdr * pow(p_rgrid->r[i - 1], k);
    Bx = Bx - Fdr / pow(p_rgrid->r[i - 1], k + 1);
    vabk[i] = p_rgrid->du *
              (Ax / pow(p_rgrid->r[i], k + 1) + Bx * pow(p_rgrid->r[i], k));
  }
  for (std::size_t i = irmax; i < (std::size_t)p_rgrid->ngp; i++) {
    vabk[i] = 0; // maybe not needed?
  }
}

//******************************************************************************
void HartreeFock::form_vbb0()
// When doing Hartree (no exchange) only need v^0_bb
// Don't call this as well as form_vabk_core, not needed (won't break though)
{
  for (std::size_t b = 0; b < m_num_core_states; b++) {
    calculate_v_abk(p_wf->orbitals[b], p_wf->orbitals[b], 0,
                    m_arr_v_abk_r[b][b][0]);
  }
}

//******************************************************************************
void HartreeFock::form_vabk_core()
// Calculates [calls calculate_v_abk] and stores the v^k_ab Hartree-Fock
// sreening functions for (a,b) in the core.
// Takes advantage of a/b symmetry.
// Skips if Lambda=0 (integral=0 from angles) Note: only for core-core states!
// (for now?)
{
#pragma omp parallel for
  for (std::size_t a = 0; a < m_num_core_states; a++) {
    for (std::size_t b = 0; b <= a; b++) {
      int kmin = abs(twoj_list[a] - twoj_list[b]) / 2;
      int kmax = (twoj_list[a] + twoj_list[b]) / 2;
      for (int k = kmin; k <= kmax; k++) {
        if (get_Lambda_abk(a, b, k) == 0)
          continue;
        calculate_v_abk(p_wf->orbitals[a], p_wf->orbitals[b], k,
                        m_arr_v_abk_r[a][b][k - kmin]);
      } // k
    }   // b
  }     // a
}

//******************************************************************************
void HartreeFock::form_vabk_valence(std::size_t w)
// Calculates [calls calculate_v_abk] and stores the Hartree-Fock screening
// functions v^k_wb for a single (given) valence state (w=valence, b=core).
// Stores in m_arr_v_abk_r
{
#pragma omp parallel for
  for (std::size_t b = 0; b < m_num_core_states; b++) {
    int kmin = abs(twoj_list[w] - twoj_list[b]) / 2;
    int kmax = (twoj_list[w] + twoj_list[b]) / 2;
    for (int k = kmin; k <= kmax; k++) {
      if (get_Lambda_abk(w, b, k) == 0)
        continue;
      calculate_v_abk(p_wf->orbitals[w], p_wf->orbitals[b], k,
                      m_arr_v_abk_r[w][b][k - kmin]);
    } // k
  }   // b
}

//******************************************************************************
const std::vector<std::vector<double>> &
HartreeFock::get_v_abk(std::size_t a, std::size_t b) const
// Returns a reference to a 2D-array (a subset of the m_arr_v_abk_r array)
// Returned array is of form: array[ik][r]; ik runs from 0 -> |kmax-kmin+1|
//   array.size() = |kmax-kmin+1| = number of k's
//   array[0].size() = ngp
// Allows to call for any a,b, even though only calculated for a>=b (symmetry)
{
  return (a > b) ? m_arr_v_abk_r[a][b] : m_arr_v_abk_r[b][a];
}
//******************************************************************************
const std::vector<double> &HartreeFock::get_v_aa0(std::size_t a) const
// Same as above, but for v^0_aa, only need to return 1D array: array[r]
// array.size()=ngp
{
  return m_arr_v_abk_r[a][a][0];
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
  for (auto &v_dir : vdir) {
    v_dir = 0;
  }
  double sf = re_scale ? (1. - 1. / p_wf->Ncore()) : 1;
  for (std::size_t b = 0; b < m_num_core_states; b++) {
    double f = (twoj_list[b] + 1) * p_wf->orbitals[b].occ_frac;
    const std::vector<double> &v0bb = get_v_aa0(b);
    for (std::size_t i = 0; i < p_rgrid->ngp; i++) {
      vdir[i] += f * v0bb[i] * sf;
    }
  } // b
}

//******************************************************************************
void HartreeFock::form_approx_vex_core(
    std::vector<std::vector<double>> &vex) const
/*
Forms the 2D "approximate" exchange potential for each core state, a.
NOTE: Must call form_vabk_core first!
Doesn't calculate, assumes m_arr_v_abk_r array exists + is up-to-date
*/
{
#pragma omp parallel for
  for (std::size_t a = 0; a < m_num_core_states; a++) {
    form_approx_vex_a(a, vex[a]);
  }
}

//******************************************************************************
void HartreeFock::form_approx_vex_a(std::size_t a,
                                    std::vector<double> &vex_a) const
// Forms the 2D "approximate" exchange potential for given core state, a.
// Does the a=b case seperately, since it's a little simpler
// Approximate:
// In order to approximate solution to HF equations, I form "local" ex.
// potential
//   [v_ex*psi_a](r) = \sum_b v_ex^(a,b)(r) * psi_b(r)
// v_ex is non-local; cannot write: [v_ex*psi_a](r) =/= v_ex(r)*psi_a(r)
// Instead: define local approx: vex_a
//   vex_a = [v_ex*psi_a](r) *(psi_a/psi_a^2)
//         = \sum_b v_ex^(a,b)(r)*psi_b(r) * (psi_a/psi_a^2)
//         = \sum_b v_ex^(a,b)(r)*(psi_b(r)*psi_a) / psi_a^2
// This vex_a is then a local potential (different for each state!) that can
// be used as an addition to local direct potential to solve Dirac Eq. as
// normal. In theory, this is exact. Clearly, however, there is an issue when
// psi_a is small. Luckily, however, we don't care as much when psi_a is small!
// Also, since v_ex is already small (compared to vdir), we can make good
// approximation. Therefore, I only calculate vex_a when a=b, or when |psi_a|
// > 1.e3 Further, largest part of v_ex is when a=b. In this case, the factor=1
// is exact!
{
  for (auto &va : vex_a) {
    va = 0;
  }

  auto &phi_a = p_wf->orbitals[a];

  if (!m_excludeExchange) {
    for (std::size_t b = 0; b < m_num_core_states; b++) { // b!=a
      if (b == a)
        continue;
      auto &phi_b = p_wf->orbitals[b];
      double x_tjbp1 = (twoj_list[b] + 1) * phi_b.occ_frac;
      auto irmax = std::min(phi_a.pinf, phi_b.pinf);
      int kmin = abs(twoj_list[a] - twoj_list[b]) / 2;
      int kmax = (twoj_list[a] + twoj_list[b]) / 2;
      const std::vector<std::vector<double>> &vabk = get_v_abk(a, b);
      // hold "fraction" psi_a*psi_b/(psi_a^2):
      std::vector<double> v_Fab(p_rgrid->ngp);
      for (std::size_t i = 0; i < irmax; i++) {
        // This is the approximte part! Divides by psi_a
        if (fabs(phi_a.f[i]) < 1.e-3)
          continue;
        double fac_top = phi_a.f[i] * phi_b.f[i] + phi_a.g[i] * phi_b.g[i];
        double fac_bot = phi_a.f[i] * phi_a.f[i] + phi_a.g[i] * phi_a.g[i];
        v_Fab[i] = -1. * x_tjbp1 * fac_top / fac_bot;
      } // r
      for (int k = kmin; k <= kmax; k++) {
        double L_abk = get_Lambda_abk(a, b, k);
        if (L_abk == 0)
          continue;
        for (std::size_t i = 0; i < irmax; i++) {
          if (v_Fab[i] == 0)
            continue;
          vex_a[i] += L_abk * vabk[k - kmin][i] * v_Fab[i];
        } // r
      }   // k
    }     // b
  }

  // now, do a=b, ONLY if a is in the core!
  if (a < m_num_core_states) {
    double x_tjap1 = (twoj_list[a] + 1); // no occ_frac here
    int kmax = twoj_list[a];
    const auto &vaak = get_v_abk(a, a);
    auto irmax = phi_a.pinf;
    for (int k = 0; k <= kmax; k++) {
      double L_abk = get_Lambda_abk(a, a, k);
      if (L_abk == 0)
        continue;
      for (std::size_t i = 0; i < irmax; i++) {
        // nb: need to 'cut' here, or fails w/ f states...
        vex_a[i] += -1 * L_abk * vaak[k][i] * x_tjap1;
      }
    } // k
  }   // if a in core
}
