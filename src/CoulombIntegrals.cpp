#include "CoulombIntegrals.hpp"
#include "DiracSpinor.hpp"
#include "NumCalc_quadIntegrate.hpp"
#include <vector>

//******************************************************************************
Coulomb::Coulomb(const std::vector<DiracSpinor> &in_core,
                 const std::vector<DiracSpinor> &in_valence = {})
    : c_orbs_ptr(&in_core), v_orbs_ptr(&in_valence),
      rgrid_ptr((in_core.size() > 0) ? in_core.front().p_rgrid
                                     : in_valence.front().p_rgrid) {
  // Call initialise here?
  std::cout << "Called A\n";
  std::cout << "core: " << c_orbs_ptr->size() << "\n";
  std::cout << "Val : " << v_orbs_ptr->size() << "\n";
  std::cout << "core: *" << c_orbs_ptr << "\n";
  std::cout << "Val : *" << v_orbs_ptr << "\n";
  initialise_core_core();
}

//******************************************************************************
void Coulomb::initialise_core_core() {
  // XXX should only be able to call this once?
  auto ngp = rgrid_ptr->ngp;
  m_y_abkr.reserve(c_orbs_ptr->size());
  for (std::size_t ia = 0; ia < c_orbs_ptr->size(); ia++) {
    auto tja = (*c_orbs_ptr)[ia].twoj();
    std::vector<std::vector<std::vector<double>>> ya_bkr;
    for (std::size_t ib = 0; ib <= ia; ib++) {
      auto tjb = (*c_orbs_ptr)[ib].twoj();
      std::size_t num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
      ya_bkr.emplace_back(
          std::vector<std::vector<double>>(num_k, std::vector<double>(ngp)));
    }
    m_y_abkr.push_back(ya_bkr);
    calculate_angular((*c_orbs_ptr)[ia].k_index());
  }
}
//------------------------------------------------------------------------------
void Coulomb::initialise_core_valence() {
  // can call this more than once..

  auto ngp = rgrid_ptr->ngp;
  for (std::size_t iv = num_initialised_vc; iv < v_orbs_ptr->size(); iv++) {
    auto tjv = (*v_orbs_ptr)[iv].twoj();
    std::vector<std::vector<std::vector<double>>> va_bkr;
    for (const auto &phi_c : *c_orbs_ptr) {
      auto tjc = phi_c.twoj();
      std::size_t num_k = (tjv > tjc) ? (tjc + 1) : (tjv + 1);
      va_bkr.emplace_back(
          std::vector<std::vector<double>>(num_k, std::vector<double>(ngp)));
    }
    m_y_vckr.push_back(va_bkr);
    calculate_angular((*v_orbs_ptr)[iv].k_index());
  }
  num_initialised_vc = v_orbs_ptr->size();
}
//------------------------------------------------------------------------------
void Coulomb::initialise_valence_valence() {
  auto ngp = rgrid_ptr->ngp;

  for (std::size_t iv = num_initialised_vv; iv < v_orbs_ptr->size(); iv++) {
    auto tjv = (*v_orbs_ptr)[iv].twoj();
    std::vector<std::vector<std::vector<double>>> vv_wkr;
    for (std::size_t iw = 0; iw <= iv; iw++) {
      auto tjw = (*v_orbs_ptr)[iw].twoj();
      std::size_t num_k = (tjv > tjw) ? (tjw + 1) : (tjv + 1);
      vv_wkr.emplace_back(
          std::vector<std::vector<double>>(num_k, std::vector<double>(ngp)));
    }
    m_y_vwkr.push_back(vv_wkr);
    calculate_angular((*v_orbs_ptr)[iv].k_index());
  }
  num_initialised_vv = v_orbs_ptr->size();
}

//******************************************************************************
std::size_t Coulomb::find_core_index(const DiracSpinor &phi) const {
  // This is very slow..is there another way?
  auto ia = std::find(c_orbs_ptr->begin(), c_orbs_ptr->end(), phi);
  return (std::size_t)std::distance(c_orbs_ptr->begin(), ia);
}
//------------------------------------------------------------------------------
std::size_t Coulomb::find_valence_index(const DiracSpinor &phi) const {
  auto ia = std::find(v_orbs_ptr->begin(), v_orbs_ptr->end(), phi);
  return (std::size_t)std::distance(v_orbs_ptr->begin(), ia);
}

//******************************************************************************
void Coulomb::calculate_core_core() {

  for (std::size_t ia = 0; ia < c_orbs_ptr->size(); ia++) {
    const auto &phi_a = (*c_orbs_ptr)[ia];
    auto tja = phi_a.twoj();
    auto kia = phi_a.k_index();
    for (std::size_t ib = 0; ib <= ia; ib++) {
      const auto &phi_b = (*c_orbs_ptr)[ib];
      auto tjb = phi_b.twoj();
      auto kib = phi_b.k_index();
      auto kmin = abs(tja - tjb) / 2; // kmin
      auto num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
      const auto &Lk = get_angular_L_kiakib_k(kia, kib);
      for (int ik = 0; ik < num_k; ik++) {
        if (Lk[ik] == 0)
          continue;
        calculate_y_ijk(phi_a, phi_b, kmin + ik, m_y_abkr[ia][ib][ik]);
      }
    }
  }
}
//******************************************************************************
void Coulomb::calculate_valence_valence() {
  // XXX this is same as core-core; just call it, passing in different orbs!
  // (internally; externally, call seperate functions, without passing!)
  initialise_valence_valence(); // call this each time?
  for (std::size_t iv = 0; iv < v_orbs_ptr->size(); iv++) {
    const auto &phi_v = (*v_orbs_ptr)[iv];
    auto tjv = phi_v.twoj();
    auto kiv = phi_v.k_index();
    for (std::size_t iw = 0; iw <= iv; iw++) {
      const auto &phi_w = (*v_orbs_ptr)[iw];
      auto tjw = phi_w.twoj();
      auto kiw = phi_w.k_index();
      auto kmin = abs(tjv - tjw) / 2; // kmin
      auto num_k = (tjv > tjw) ? (tjw + 1) : (tjv + 1);
      const auto &Lk = get_angular_L_kiakib_k(kiv, kiw);
      for (int ik = 0; ik < num_k; ik++) {
        if (Lk[ik] == 0)
          continue;
        calculate_y_ijk(phi_v, phi_w, kmin + ik, m_y_vwkr[iv][iw][ik]);
      }
    }
  }
}
//******************************************************************************
void Coulomb::calculate_core_valence() {
  // XXX this is same as core-core; just call it, passing in different orbs!
  // (internally; externally, call seperate functions, without passing!)
  initialise_core_valence(); // call this each time?

  for (std::size_t iv = 0; iv < v_orbs_ptr->size(); iv++) {
    const auto &phi_v = (*v_orbs_ptr)[iv];
    auto tjv = phi_v.twoj();
    auto kiv = phi_v.k_index();
    for (std::size_t ic = 0; ic < c_orbs_ptr->size(); ic++) {
      const auto &phi_c = (*c_orbs_ptr)[ic];
      auto tjc = phi_c.twoj();
      auto kic = phi_c.k_index();
      auto kmin = abs(tjc - tjv) / 2; // kmin
      auto num_k = (tjc > tjv) ? (tjv + 1) : (tjc + 1);
      const auto &Lk = get_angular_L_kiakib_k(kic, kiv);
      for (int ik = 0; ik < num_k; ik++) {
        if (Lk[ik] == 0)
          continue;
        calculate_y_ijk(phi_c, phi_v, kmin + ik, m_y_vckr[iv][ic][ik]);
      }
    }
  }
}

//******************************************************************************
const std::vector<std::vector<double>> &
Coulomb::get_y_abk(std::size_t a, std::size_t b) const {
  return (a > b) ? m_y_abkr[a][b] : m_y_abkr[a][b];
}
//------------------------------------------------------------------------------
const std::vector<std::vector<double>> &
Coulomb::get_y_vwk(std::size_t v, std::size_t w) const {
  return (v > w) ? m_y_abkr[v][w] : m_y_abkr[w][v];
}
//------------------------------------------------------------------------------
const std::vector<std::vector<double>> &
Coulomb::get_y_vck(std::size_t v, std::size_t c) const {
  return m_y_vckr[v][c];
}
//******************************************************************************
const std::vector<std::vector<double>> &
Coulomb::get_y_ijk(const DiracSpinor &phi_i, const DiracSpinor &phi_j) const {

  auto i = find_core_index(phi_i);
  bool ival = false;
  if (i == c_orbs_ptr->size()) {
    i = find_valence_index(phi_i);
    ival = true;
  }
  auto j = find_core_index(phi_j);
  bool jval = false;
  if (j == c_orbs_ptr->size()) {
    j = find_valence_index(phi_j);
    jval = true;
  }
  if (!ival && !jval)
    return (i > j) ? m_y_abkr[i][j] : m_y_abkr[i][j];
  if (ival && !jval)
    return m_y_vckr[i][j];
  if (!ival && jval)
    return m_y_vckr[j][i];
  return (i > j) ? m_y_vckr[i][j] : m_y_vckr[i][j];
}
//------------------------------------------------------------------------------
const std::vector<double> &Coulomb::get_y_ijk(const DiracSpinor &phi_i,
                                              const DiracSpinor &phi_j,
                                              int k) const {
  auto tji = phi_i.twoj();
  auto tjj = phi_j.twoj();
  auto kmin = abs(tji - tjj) / 2; // kmin
  auto kmax = (tji + tjj) / 2;    // kmax
  if (k > kmax || k < kmin) {
    std::cerr << "FAIL 214 in CI; bad k\n";
    std::abort();
  }
  const auto &tmp = get_y_ijk(phi_i, phi_j);
  return tmp[k - kmin];
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
double Coulomb::get_angular_C_kiakibk(const DiracSpinor &phi_a,
                                      const DiracSpinor &phi_b, int k) const {
  auto kia = phi_a.k_index();
  auto kib = phi_b.k_index();
  int kmin =
      abs(AtomInfo::twojFromIndex(kia) - AtomInfo::twojFromIndex(kib)) / 2;
  int kmax =
      abs(AtomInfo::twojFromIndex(kia) + AtomInfo::twojFromIndex(kib)) / 2;
  if (k < kmin || k > kmax)
    return 0;
  return kia > kib ? m_angular_C_kakbk[kia][kib][k - kmin]
                   : m_angular_C_kakbk[kib][kia][k - kmin];
}
//******************************************************************************
const std::vector<double> &Coulomb::get_angular_C_kiakib_k(int kia,
                                                           int kib) const {
  // note:output is of-set by k_min!
  return kia > kib ? m_angular_C_kakbk[kia][kib] : m_angular_C_kakbk[kib][kia];
}

//******************************************************************************
double Coulomb::get_angular_L_kiakibk(const DiracSpinor &phi_a,
                                      const DiracSpinor &phi_b, int k) const {
  auto kia = phi_a.k_index();
  auto kib = phi_b.k_index();
  int kmin =
      abs(AtomInfo::twojFromIndex(kia) - AtomInfo::twojFromIndex(kib)) / 2;
  int kmax =
      abs(AtomInfo::twojFromIndex(kia) + AtomInfo::twojFromIndex(kib)) / 2;
  if (k < kmin || k > kmax)
    return 0;
  return kia > kib ? m_angular_L_kakbk[kia][kib][k - kmin]
                   : m_angular_L_kakbk[kib][kia][k - kmin];
}
//******************************************************************************
const std::vector<double> &Coulomb::get_angular_L_kiakib_k(int kia,
                                                           int kib) const {
  // note:output is off-set by k_min!
  return kia > kib ? m_angular_L_kakbk[kia][kib] : m_angular_L_kakbk[kib][kia];
}

//******************************************************************************
void Coulomb::calculate_y_ijk(const DiracSpinor &phi_a,
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
