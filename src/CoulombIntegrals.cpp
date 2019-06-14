#include "CoulombIntegrals.hpp"
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <vector>

//******************************************************************************
void Coulomb::calculate_v_abk(const DiracSpinor &phi_a,
                              const DiracSpinor &phi_b, const int k,
                              std::vector<double> &vabk)
// This is static
// Calculalates v^k_ab screening function.
// Note: should only call for a>=b, and for k's with non-zero angular coefs
// (nothing bad will happen otherwise, but no point!)
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
  auto du = phi_a.p_rgrid->du;

  double Ax = 0, Bx = 0;
  for (std::size_t i = 0; i < irmax; i++) {
    Bx += phi_a.p_rgrid->drdu[i] *
          (phi_a.f[i] * phi_b.f[i] + phi_a.g[i] * phi_b.g[i]) /
          pow(phi_a.p_rgrid->r[i], k + 1);
  }

  // For "direct" part, can't cut!
  if (phi_a == phi_b)
    irmax = phi_a.p_rgrid->ngp;

  vabk[0] = Bx * du;
  for (std::size_t i = 1; i < irmax; i++) {
    auto Fdr = phi_a.p_rgrid->drdu[i - 1] * (phi_a.f[i - 1] * phi_b.f[i - 1] +
                                             phi_a.g[i - 1] * phi_b.g[i - 1]);
    Ax = Ax + Fdr * pow(phi_a.p_rgrid->r[i - 1], k);
    Bx = Bx - Fdr / pow(phi_a.p_rgrid->r[i - 1], k + 1);
    vabk[i] = du * (Ax / pow(phi_a.p_rgrid->r[i], k + 1) +
                    Bx * pow(phi_a.p_rgrid->r[i], k));
  }
  for (std::size_t i = irmax; i < phi_a.p_rgrid->ngp; i++) {
    vabk[i] = 0; // maybe not needed?
  }
}

//******************************************************************************
std::vector<std::vector<double>> &Coulomb::get_vab_kr(const DiracSpinor &phi_a,
                                                      const DiracSpinor &phi_b)
// returns v_ab  (v_ab[k][r] - 2D vector)
// Uses symmetry, v_ab = v_ba
// nb: this not const, cos used in construction!
{

  if (std::find(nka_list.begin(), nka_list.end(), State(phi_a.n, phi_a.k)) ==
      nka_list.end()) {
    // phi_a not in nka list (means I should swap phi_a and phi_b)
  }

  auto ia =
      std::find(nka_list.begin(), nka_list.end(), State(phi_a.n, phi_a.k)) -
      nka_list.begin();
  auto ib =
      std::find(nkb_list.begin(), nkb_list.end(), State(phi_b.n, phi_b.k)) -
      nkb_list.begin();
  if (ia == (int)nka_list.size()) {
    ia = std::find(nkb_list.begin(), nkb_list.end(), State(phi_a.n, phi_a.k)) -
         nkb_list.begin();
    ib = std::find(nka_list.begin(), nka_list.end(), State(phi_b.n, phi_b.k)) -
         nka_list.begin();
  }

  return (phi_a > phi_b) ? m_v_abkr[ia][ib] : m_v_abkr[ib][ia];
}

//******************************************************************************
const std::vector<double> &Coulomb::get_vabk_r(const DiracSpinor &phi_a,
                                               const DiracSpinor &phi_b, int k)
// returns v^k_ab  (v_abk[r] - 1D vector)
// Uses symmetry, v_ab = v_ba
{
  auto ia =
      std::find(nka_list.begin(), nka_list.end(), State(phi_a.n, phi_a.k)) -
      nka_list.begin();
  auto ib =
      std::find(nkb_list.begin(), nkb_list.end(), State(phi_b.n, phi_b.k)) -
      nkb_list.begin();
  if (ia == (int)nka_list.size()) {
    ia = std::find(nkb_list.begin(), nkb_list.end(), State(phi_a.n, phi_a.k)) -
         nkb_list.begin();
    ib = std::find(nka_list.begin(), nka_list.end(), State(phi_b.n, phi_b.k)) -
         nka_list.begin();
  }

  int kmin = abs(phi_a.twoj() - phi_b.twoj()) / 2; // kmin
  return (phi_a > phi_b) ? m_v_abkr[ia][ib][k - kmin]
                         : m_v_abkr[ib][ia][k - kmin];
}

//******************************************************************************
void Coulomb::form_v_abk(const std::vector<DiracSpinor> &a_orbitals,
                         const std::vector<DiracSpinor> &b_orbitals) {
  for (const auto &phi_a : a_orbitals) {
    form_v_abk(phi_a, b_orbitals);
  }
}
//******************************************************************************
void Coulomb::form_v_abk(const DiracSpinor &phi_a,
                         const std::vector<DiracSpinor> &b_orbitals) {
  auto tja = phi_a.twoj();
  for (const auto &phi_b : b_orbitals) {
    if (phi_b > phi_a)
      continue;
    auto tjb = phi_b.twoj();
    auto &vabk = get_vab_kr(phi_a, phi_b);
    int k = abs(tja - tjb) / 2; // kmin
    for (auto &vab_k : vabk) {
      calculate_v_abk(phi_a, phi_b, k, vab_k);
      k++;
    }
  }
}

//******************************************************************************
void Coulomb::initialise_v_abkr(const std::vector<DiracSpinor> &a_orbitals,
                                const std::vector<DiracSpinor> &b_orbitals) {
  for (const auto &phi_b : b_orbitals) {
    nkb_list.emplace_back(phi_b.n, phi_b.k);
  }
  for (const auto &phi_a : a_orbitals) {
    extend_v_abkr(phi_a, b_orbitals);
  }
}
//******************************************************************************
void Coulomb::initialise_v_abkr(const std::vector<DiracSpinor> &orbitals) {
  for (const auto &phi : orbitals) {
    nkb_list.emplace_back(phi.n, phi.k);
    extend_v_abkr(phi, orbitals);
  }
}

//******************************************************************************
void Coulomb::extend_v_abkr(const DiracSpinor &phi_a,
                            const std::vector<DiracSpinor> &b_orbitals)
// Prepare new orbital
{

  // // XXX This adds all the phi_a's to the list.. but not the phi_b's!
  // // If phi_a is in {b_orbs}, then all fine!
  // // If phi_a is not in {b_orbs}, means we need to add b_obs first!
  // // This is annoyingly slow....
  // if (std::find(b_orbitals.begin(), b_orbitals.end(), phi_a) ==
  //     b_orbitals.end()) {
  //   // means phi_a is NOT in phi_b
  //
  // }

  // First, make sure not already in list (SAFETY)
  State nk{phi_a.n, phi_a.k};
  if (std::find(nka_list.begin(), nka_list.end(), nk) != nka_list.end())
    return;

  // add orbital info to list:
  nka_list.push_back(nk);
  auto ki = phi_a.k_index();
  ki_list.push_back(ki);
  auto tja = AtomInfo::twojFromIndex(ki);
  twoj_list.push_back(tja);

  // extend m_v_abkr array
  auto ngp = phi_a.p_rgrid->ngp;
  std::vector<std::vector<std::vector<double>>> va_bkr;
  for (const auto &phi_b : b_orbitals) {
    if (phi_b > phi_a)
      continue; // not break, may not be in order
    auto tjb = phi_b.twoj();
    std::size_t num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
    va_bkr.emplace_back(
        std::vector<std::vector<double>>(num_k, std::vector<double>(ngp)));
  }
  m_v_abkr.push_back(va_bkr);

  calculate_angular(ki);
}

//******************************************************************************
void Coulomb::calculate_angular(int ki)
// Extend the angular coeficient (AND calculate it!)
{
  if (ki <= m_largest_ki)
    return;
  auto prev_largest_ki = m_largest_ki;
  m_largest_ki = ki;
  for (auto kia = prev_largest_ki + 1; kia <= m_largest_ki; kia++) {
    auto tja = AtomInfo::twojFromIndex(kia);
    auto la = AtomInfo::lFromIndex(kia);
    std::vector<std::vector<double>> C_ka_kbk;
    std::vector<std::vector<double>> L_ka_kbk;
    for (auto kib = 0; kib <= kia; kib++) {
      auto tjb = AtomInfo::twojFromIndex(kib);
      auto lb = AtomInfo::lFromIndex(kib);
      auto kmin = (tja - tjb) / 2; // don't need abs, as m\leq n => ja\geq jb
      auto kmax = (tja + tjb) / 2;
      std::vector<double> C_k(kmax - kmin + 1, 0);
      std::vector<double> L_k(kmax - kmin + 1, 0);
      for (auto k = kmin; k <= kmax; k++) {
        if (Wigner::parity(la, lb, k) == 0)
          continue;
        int ik = k - kmin;
        auto tjs = Wigner::threej_2(tja, tjb, 2 * k, -1, 1, 0);
        // auto sign = ((tja + 1) / 2 % 2 == 0) ? 1 : -1;
        C_k[ik] = sqrt((tja + 1) * (tjb + 1)) * tjs; // XXX no sign!
        L_k[ik] = tjs * tjs;
      } // k
      C_ka_kbk.push_back(C_k);
      L_ka_kbk.push_back(L_k);
    }
    m_angular_C_kakbk.push_back(C_ka_kbk);
    m_angular_L_kakbk.push_back(L_ka_kbk);
  }
}

//******************************************************************************
double Coulomb::get_angular_C_kiakibk(int kia, int kib, int k) {
  int kmin =
      abs(AtomInfo::twojFromIndex(kia) - AtomInfo::twojFromIndex(kib)) / 2;
  return kia > kib ? m_angular_C_kakbk[kia][kib][k - kmin]
                   : m_angular_C_kakbk[kib][kia][k - kmin];
}
//******************************************************************************
const std::vector<double> &Coulomb::get_angular_C_kiakib_k(int kia, int kib) {
  // note:output is of-set by k_min!
  return kia > kib ? m_angular_C_kakbk[kia][kib] : m_angular_C_kakbk[kib][kia];
}

//******************************************************************************
double Coulomb::get_angular_L_kiakibk(int kia, int kib, int k) {
  int kmin =
      abs(AtomInfo::twojFromIndex(kia) - AtomInfo::twojFromIndex(kib)) / 2;
  return kia > kib ? m_angular_L_kakbk[kia][kib][k - kmin]
                   : m_angular_L_kakbk[kib][kia][k - kmin];
}
//******************************************************************************
const std::vector<double> &Coulomb::get_angular_L_kiakib_k(int kia, int kib) {
  // note:output is of-set by k_min!
  return kia > kib ? m_angular_L_kakbk[kia][kib] : m_angular_L_kakbk[kib][kia];
}
