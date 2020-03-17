#include "HF/CoulombIntegrals.hpp"
#include "Angular/Angular_369j.hpp"
#include "Angular/Angular_tables.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

//******************************************************************************
double Coulomb::Zk_abcd_any(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const int k) {
  // Z^k_abcd = s ( Q^k_abcd + sum_l [k] 6j * Q^l_abdc)

  auto s = Angular::evenQ_2(Fa.twoj() + Fb.twoj() + 2) ? 1 : -1;
  auto Qk_abcd = Qk_abcd_any(Fa, Fb, Fc, Fd, k);

  auto tkp1 = 2 * k + 1;

  auto min_twol = std::max(std::abs(Fd.twoj() - Fa.twoj()),
                           std::abs(Fc.twoj() - Fb.twoj()));

  auto max_twol = std::min(Fd.twoj() + Fa.twoj(), Fc.twoj() + Fb.twoj());

  double sum = 0.0;
  for (int tl = min_twol; tl <= max_twol; tl += 2) {
    auto sixj = Angular::sixj_2(Fc.twoj(), Fa.twoj(), 2 * k, //
                                Fd.twoj(), Fb.twoj(), tl);
    if (sixj == 0)
      continue;
    auto Ql_abdc = Qk_abcd_any(Fa, Fb, Fd, Fc, tl / 2);
    sum += sixj * Ql_abdc;
  }

  return s * (Qk_abcd + tkp1 * sum);
}

//******************************************************************************
double Coulomb::Qk_abcd_any(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const int k) {

  auto tCac = Angular::tildeCk_kk(k, Fa.k, Fc.k);
  if (tCac == 0.0)
    return 0.0;
  auto tCbd = Angular::tildeCk_kk(k, Fb.k, Fd.k);
  if (tCbd == 0.0)
    return 0.0;
  auto Rkabcd = Rk_abcd_any(Fa, Fb, Fc, Fd, k);
  auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rkabcd;
}

//******************************************************************************
double Coulomb::Rk_abcd_any(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd,
                            const int k) //
{
  //
  auto imax = Angular::max4(Fa.pinf, Fb.pinf, Fc.pinf, Fd.pinf);
  const auto &drdu = Fa.p_rgrid->drdu; // save typing
  const auto du = Fa.p_rgrid->du;

  std::vector<double> yk_bd; // XXX Can already exist!
  calculate_y_ijk(Fb, Fd, k, yk_bd);
  auto Rff = NumCalc::integrate(1.0, 0, imax, Fa.f, Fc.f, yk_bd, drdu);
  auto Rgg = NumCalc::integrate(1.0, 0, imax, Fa.g, Fc.g, yk_bd, drdu);
  return (Rff + Rgg) * du;
}

//******************************************************************************
DiracSpinor Coulomb::Rk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                 const DiracSpinor &Fc, const DiracSpinor &Fd,
                                 const int k) //
{
  // make rhs an in/out parameter??
  // does 2 allocations [3v's: 2 for DS, onr for y]
  std::vector<double> yk_bd; // XXX Can already exist!
  calculate_y_ijk(Fb, Fd, k, yk_bd, Fc.pinf);
  auto out = DiracSpinor(0, Fa.k, *(Fa.p_rgrid));
  out.pinf = Fc.pinf;
  // out.f = Fc.f;
  // out.g = Fc.g;
  // out *= yk_bd;
  out.f = NumCalc::mult_vectors(Fc.f, yk_bd);
  out.g = NumCalc::mult_vectors(Fc.g, yk_bd);
  return out;
}

//******************************************************************************
DiracSpinor Coulomb::Qk_abcd_rhs(const DiracSpinor &Fa, const DiracSpinor &Fb,
                                 const DiracSpinor &Fc, const DiracSpinor &Fd,
                                 const int k) {
  // XXX have bool! can skip!
  auto tCac = Angular::tildeCk_kk(k, Fa.k, Fc.k);
  if (tCac == 0.0)
    return 0.0 * Fa;
  auto tCbd = Angular::tildeCk_kk(k, Fb.k, Fd.k);
  if (tCbd == 0.0)
    return 0.0 * Fa;
  auto Rk_bcd = Rk_abcd_rhs(Fa, Fb, Fc, Fd, k);
  auto m1tk = Angular::evenQ(k) ? 1 : -1;
  return m1tk * tCac * tCbd * Rk_bcd;
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//******************************************************************************
Coulomb::Coulomb(const std::vector<DiracSpinor> &in_core,
                 const std::vector<DiracSpinor> &in_valence = {})
    : c_orbs_ptr(&in_core), v_orbs_ptr(&in_valence),
      rgrid_ptr(in_core.empty() ? in_valence.front().p_rgrid
                                : in_core.front().p_rgrid) {
  initialise_core_core();
  // Valence orbitals may change (new orbitals added) - but pointer to val.
  // orbitals remains const
}
//******************************************************************************
Coulomb::Coulomb(const Grid &in_grid, const std::vector<DiracSpinor> &in_core,
                 const std::vector<DiracSpinor> &in_valence = {})
    : c_orbs_ptr(&in_core), v_orbs_ptr(&in_valence), rgrid_ptr(&in_grid)
// This should probably be the default one...
{
  initialise_core_core();
}

//******************************************************************************
void Coulomb::initialise_core_core()
// Initialises memory (sizes arays) used to store core-core y^k_ab C ints
{
  auto num_points = rgrid_ptr->num_points;
  m_y_abkr.clear();
  m_y_abkr.reserve(c_orbs_ptr->size());
  for (std::size_t ia = 0; ia < c_orbs_ptr->size(); ia++) {
    auto tja = (*c_orbs_ptr)[ia].twoj();
    std::vector<std::vector<std::vector<double>>> ya_bkr;
    for (std::size_t ib = 0; ib <= ia; ib++) {
      auto tjb = (*c_orbs_ptr)[ib].twoj();
      std::size_t num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
      ya_bkr.emplace_back(std::vector<std::vector<double>>(
          num_k, std::vector<double>(num_points)));
    }
    m_y_abkr.push_back(ya_bkr);
    calculate_angular((*c_orbs_ptr)[ia].k_index());
  }
}
//------------------------------------------------------------------------------
void Coulomb::initialise_core_valence()
// Initialises memory (sizes arays) used to store core-valence y^k_vc C ints
// May be called many times.
// Must be called after new valence orbitals are constructed before can
// calculate the new integrals (for now, this happens automatically)
{
  auto num_points = rgrid_ptr->num_points;
  for (std::size_t iv = num_initialised_vc; iv < v_orbs_ptr->size(); iv++) {
    auto tjv = (*v_orbs_ptr)[iv].twoj();
    std::vector<std::vector<std::vector<double>>> va_bkr;
    for (const auto &Fc : *c_orbs_ptr) {
      auto tjc = Fc.twoj();
      std::size_t num_k = (tjv > tjc) ? (tjc + 1) : (tjv + 1);
      va_bkr.emplace_back(std::vector<std::vector<double>>(
          num_k, std::vector<double>(num_points)));
    }
    m_y_vckr.push_back(va_bkr);
    calculate_angular((*v_orbs_ptr)[iv].k_index());
  }
  num_initialised_vc = v_orbs_ptr->size();
}
//------------------------------------------------------------------------------
void Coulomb::initialise_valence_valence()
// Initialises memory (sizes arays) used to store valence-valence y^k_vw C ints
// May be called many times.
// Must be called after new valence orbitals are constructed before can
// calculate the new integrals (for now, this happens automatically)
{
  auto num_points = rgrid_ptr->num_points;
  for (std::size_t iv = num_initialised_vv; iv < v_orbs_ptr->size(); iv++) {
    auto tjv = (*v_orbs_ptr)[iv].twoj();
    std::vector<std::vector<std::vector<double>>> vv_wkr;
    for (std::size_t iw = 0; iw <= iv; iw++) {
      auto tjw = (*v_orbs_ptr)[iw].twoj();
      std::size_t num_k = (tjv > tjw) ? (tjw + 1) : (tjv + 1);
      vv_wkr.emplace_back(std::vector<std::vector<double>>(
          num_k, std::vector<double>(num_points)));
    }
    m_y_vwkr.push_back(vv_wkr);
    calculate_angular((*v_orbs_ptr)[iv].k_index());
  }
  num_initialised_vv = v_orbs_ptr->size();
}

//******************************************************************************
std::size_t Coulomb::find_core_index(const DiracSpinor &Fa) const
// This finds the array index of a particular orbital
// Note: if orbital is not in the list (e.g., if called with a valence orbital
// instead of a core orbital), will return [size_of_orbitals], which is an
// invalid index (can cause undefined behaviour!) This is very slow..is there
// another way?
{
  auto ia = std::find(c_orbs_ptr->begin(), c_orbs_ptr->end(), Fa);
  return (std::size_t)std::distance(c_orbs_ptr->begin(), ia);
}
//------------------------------------------------------------------------------
std::size_t Coulomb::find_valence_index(const DiracSpinor &Fa) const {
  auto ia = std::find(v_orbs_ptr->begin(), v_orbs_ptr->end(), Fa);
  return (std::size_t)std::distance(v_orbs_ptr->begin(), ia);
}
//------------------------------------------------------------------------------
std::size_t Coulomb::find_either_index(const DiracSpinor &Fa,
                                       bool &valenceQ) const {
  auto sp1 = SafeProfiler::profile(__func__);
  auto ia = std::find(c_orbs_ptr->begin(), c_orbs_ptr->end(), Fa);
  auto i = (std::size_t)std::distance(c_orbs_ptr->begin(), ia);
  valenceQ = false;
  if (i == c_orbs_ptr->size()) {
    valenceQ = true;
    ia = std::find(v_orbs_ptr->begin(), v_orbs_ptr->end(), Fa);
    i = (std::size_t)std::distance(v_orbs_ptr->begin(), ia);
  }
  return i;
}

//******************************************************************************
void Coulomb::form_core_core()
// Calls calculate_y_ijk, fills the core-core C int arrays
// Note: symmety: y_ij = y_ji, therefore only calculates y_ij with i >= j
{
  auto sp1 = SafeProfiler::profile(__func__, "a");
  auto Ncore = c_orbs_ptr->size();
#pragma omp parallel for
  for (std::size_t ia = 0; ia < Ncore; ia++) {
    const auto &Fa = (*c_orbs_ptr)[ia];
    auto tja = Fa.twoj();
    auto kia = Fa.k_index();
    for (std::size_t ib = 0; ib <= ia; ib++) {
      const auto &Fb = (*c_orbs_ptr)[ib];
      auto tjb = Fb.twoj();
      auto kib = Fb.k_index();
      auto kmin = std::abs(tja - tjb) / 2;
      auto num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
      const auto &Lk = get_angular_L_kiakib_k(kia, kib);
      for (int ik = 0; ik < num_k; ik++) {
        if (Lk[ik] == 0)
          continue;
        calculate_y_ijk(Fa, Fb, kmin + ik, m_y_abkr[ia][ib][ik]);
      }
    }
  }
}
//******************************************************************************
void Coulomb::form_core_core(const DiracSpinor &Fa)
// Calls calculate_y_ijk, only for terms involving Fa!
// Note: symmety: y_ij = y_ji, therefore only calculates y_ij with i >= j
// NOTE: ONLY call this if we only need to update one!
// Inneficient to call it for each orbital, since only need for a<=b !
{
  auto sp1 = SafeProfiler::profile(__func__, "b");
  auto Ncore = c_orbs_ptr->size();
  auto ia = find_core_index(Fa); // slow?
  auto tja = Fa.twoj();
  auto kia = Fa.k_index();

#pragma omp parallel for
  for (std::size_t ib = 0; ib < Ncore; ib++) {
    const auto &Fb = (*c_orbs_ptr)[ib];
    auto tjb = Fb.twoj();
    auto kib = Fb.k_index();
    auto kmin = std::abs(tja - tjb) / 2;
    auto num_k = (tja > tjb) ? (tjb + 1) : (tja + 1);
    const auto &Lk = get_angular_L_kiakib_k(kia, kib);
    for (int ik = 0; ik < num_k; ik++) {
      if (Lk[ik] == 0)
        continue;
      auto &yab = ia > ib ? m_y_abkr[ia][ib][ik] : m_y_abkr[ib][ia][ik];
      calculate_y_ijk(Fa, Fb, kmin + ik, yab);
    }
  }
}

//******************************************************************************
void Coulomb::form_valence_valence()
// Calls calculate_y_ijk, fills the valence-valence C int arrays
// Note: symmety: y_ij = y_ji, therefore only calculates y_ij with i >= j
{
  auto sp1 = SafeProfiler::profile(__func__);
  initialise_valence_valence(); // call this each time?
  auto Nval = v_orbs_ptr->size();
#pragma omp parallel for
  for (std::size_t iv = 0; iv < Nval; iv++) {
    const auto &Fv = (*v_orbs_ptr)[iv];
    auto tjv = Fv.twoj();
    auto kiv = Fv.k_index();
    for (std::size_t iw = 0; iw <= iv; iw++) {
      const auto &Fw = (*v_orbs_ptr)[iw];
      auto tjw = Fw.twoj();
      auto kiw = Fw.k_index();
      auto kmin = std::abs(tjv - tjw) / 2;
      auto num_k = (tjv > tjw) ? (tjw + 1) : (tjv + 1);
      const auto &Lk = get_angular_L_kiakib_k(kiv, kiw);
      for (int ik = 0; ik < num_k; ik++) {
        if (Lk[ik] == 0)
          continue;
        calculate_y_ijk(Fv, Fw, kmin + ik, m_y_vwkr[iv][iw][ik]);
      }
    }
  }
}
//******************************************************************************
void Coulomb::form_core_valence(const DiracSpinor &Fn)
// Calls calculate_y_ijk, only for terms involving Fn!
// If Fn is valence, calculates Fn-core integrals
// If Fn is core, calculates valence-Fn integrals
{
  auto sp1 = SafeProfiler::profile(__func__, "a");
  initialise_core_valence();

  bool n_valenceQ{};
  auto in = find_either_index(Fn, n_valenceQ);
  auto tjn = Fn.twoj();
  auto kin = Fn.k_index();

  //"other" orbital set (opposite to Fn):
  const auto &orbs = n_valenceQ ? *c_orbs_ptr : *v_orbs_ptr;

#pragma omp parallel for
  for (std::size_t im = 0; im < orbs.size(); im++) {
    const auto &Fm = orbs[im];
    auto tjm = Fm.twoj();
    auto kim = Fm.k_index();
    auto kmin = std::abs(tjm - tjn) / 2;
    auto num_k = (tjm > tjn) ? (tjn + 1) : (tjm + 1);
    const auto &Lk = get_angular_L_kiakib_k(kim, kin);
    for (int ik = 0; ik < num_k; ik++) {
      if (Lk[ik] == 0)
        continue;
      auto &yvc = n_valenceQ ? m_y_vckr[in][im][ik] : m_y_vckr[im][in][ik];
      calculate_y_ijk(Fm, Fn, kmin + ik, yvc);
    }
  }
}

//******************************************************************************
void Coulomb::form_core_valence()
// Calls calculate_y_ijk, fills the core-valence C int arrays
// Note: no symmetry here! y_ij != y_ji [j and i same index, NOT same orbital!]
{
  auto sp1 = SafeProfiler::profile(__func__, "b");
  initialise_core_valence(); // call this each time?
  auto Nval = v_orbs_ptr->size();
#pragma omp parallel for // two-level?
  for (std::size_t iv = 0; iv < Nval; iv++) {
    const auto &Fv = (*v_orbs_ptr)[iv];
    auto tjv = Fv.twoj();
    auto kiv = Fv.k_index();
    for (std::size_t ic = 0; ic < c_orbs_ptr->size(); ic++) {
      const auto &Fc = (*c_orbs_ptr)[ic];
      auto tjc = Fc.twoj();
      auto kic = Fc.k_index();
      auto kmin = std::abs(tjc - tjv) / 2;
      auto num_k = (tjc > tjv) ? (tjv + 1) : (tjc + 1);
      const auto &Lk = get_angular_L_kiakib_k(kic, kiv);
      for (int ik = 0; ik < num_k; ik++) {
        if (Lk[ik] == 0)
          continue;
        calculate_y_ijk(Fc, Fv, kmin + ik, m_y_vckr[iv][ic][ik]);
      }
    }
  }
}

//******************************************************************************
const std::vector<std::vector<double>> &
Coulomb::get_y_abk(std::size_t a, std::size_t b) const {
  return (a > b) ? m_y_abkr[a][b] : m_y_abkr[b][a];
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
Coulomb::get_y_ijk(const DiracSpinor &Fi, const DiracSpinor &Fj) const {
  auto sp1 = SafeProfiler::profile(__func__, "a");

  bool ival{}, jval{};
  auto i = find_either_index(Fi, ival);
  auto j = find_either_index(Fj, jval);

  if (!ival && !jval)
    return (i > j) ? m_y_abkr[i][j] : m_y_abkr[j][i];
  if (ival && !jval)
    return m_y_vckr[i][j];
  if (!ival && jval)
    return m_y_vckr[j][i];
  return (i > j) ? m_y_vckr[i][j] : m_y_vckr[j][i];
}
//------------------------------------------------------------------------------
const std::vector<double> &
Coulomb::get_y_ijk(const DiracSpinor &Fi, const DiracSpinor &Fj, int k) const {
  auto sp1 = SafeProfiler::profile(__func__, "b");
  auto tji = Fi.twoj();
  auto tjj = Fj.twoj();
  auto kmin = std::abs(tji - tjj) / 2; // kmin
  auto kmax = (tji + tjj) / 2;         // kmax
  if (k > kmax || k < kmin) {
    std::cerr << "FAIL 214 in CI; bad k\n";
    std::abort();
  }
  const auto &tmp = get_y_ijk(Fi, Fj);
  return tmp[k - kmin];
}

//******************************************************************************
void Coulomb::calculate_angular(int ki)
// Calculated the angular coeficients C and L.
// Automatically allocates memory (sizes arrays); can be called any number of
// times. [Called automatically]
// Stores them in arrays, indexed directly by kappa_index (ki) for 'k',
// only store between kmin and kmax, so off-set by kmin For kappa_index (kia,
// kib), and k: C[kia][kib][k-kmin] Due to symmetry, C_ab = C_ba only store for
// kia >= kib Always use getter functions to access array. Definitions: L =
// Lambda^k_ij := 3js((ji,jj,k),(-1/2,1/2,0))^2 * parity(li+lj+k) C = |
// <k||C^k||k'> |
//   = std::sqrt([ji][jj]) * 3js((ji,jj,k),(-1/2,1/2,0)) * parity(li+lj+k)
// Note: C is abs value) - if sign needed, do seperately [sign NOT symmetric!]
// Also:
// k_min = |j - j'|; k_max = |j + j'|
// num_k = (j' + 1) if j>j', (j + 1) if j'>j;
{
  auto sp1 = SafeProfiler::profile(__func__);
  if (ki <= m_largest_ki)
    return;
  auto prev_largest_ki = m_largest_ki;
  m_largest_ki = ki;
  for (auto kia = prev_largest_ki + 1; kia <= m_largest_ki; kia++) {
    auto tja = AtomData::twojFromIndex(kia);
    auto la = AtomData::lFromIndex(kia);
    std::vector<std::vector<double>> C_ka_kbk;
    std::vector<std::vector<double>> L_ka_kbk;
    for (auto kib = 0; kib <= kia; kib++) {
      auto tjb = AtomData::twojFromIndex(kib);
      auto lb = AtomData::lFromIndex(kib);
      auto kmin = (tja - tjb) / 2; // don't need abs, as b\leq a => ja\geq jb
      auto kmax = (tja + tjb) / 2;
      std::vector<double> C_k(kmax - kmin + 1, 0);
      std::vector<double> L_k(kmax - kmin + 1, 0);
      for (auto k = kmin; k <= kmax; k++) {
        if (Angular::parity(la, lb, k) == 0)
          continue;
        int ik = k - kmin;
        auto tjs = Angular::threej_2(tja, tjb, 2 * k, -1, 1, 0);
        C_k[ik] = std::sqrt((tja + 1) * (tjb + 1)) * tjs; // nb: no sign!
        L_k[ik] = tjs * tjs;
      } // k
      C_ka_kbk.push_back(C_k);
      L_ka_kbk.push_back(L_k);
    }
    m_C_kakbk.push_back(C_ka_kbk);
    m_L_kakbk.push_back(L_ka_kbk);
  }
}

//******************************************************************************
const std::vector<double> &Coulomb::get_angular_C_kiakib_k(int kia,
                                                           int kib) const {
  // note:output is of-set by k_min!
  return kia > kib ? m_C_kakbk[kia][kib] : m_C_kakbk[kib][kia];
}

//******************************************************************************
const std::vector<double> &Coulomb::get_angular_L_kiakib_k(int kia,
                                                           int kib) const {
  // note:output is off-set by k_min!
  return kia > kib ? m_L_kakbk[kia][kib] : m_L_kakbk[kib][kia];
}

//******************************************************************************
std::vector<double>
Coulomb::calculate_R_abcd_k(const DiracSpinor &Fa, const DiracSpinor &Fb,
                            const DiracSpinor &Fc, const DiracSpinor &Fd) const
// R^k_abcd = Int_0^inf [fa*fc + ga*gc]*y^k_bd(r) dr
// Symmetry: a<->c, and b<->d
// NOTE: NOT offset by k_min, so will calculate for k=0,1,2,...,k_max
{
  auto sp1 = SafeProfiler::profile(__func__);
  auto kmin = std::abs(Fb.twoj() - Fd.twoj()) / 2;
  auto kmax = std::abs(Fb.twoj() + Fd.twoj()) / 2;
  const auto &drdu = Fa.p_rgrid->drdu; // save typing
  const auto du = Fa.p_rgrid->du;

  auto pinf = std::min(Fa.pinf, Fc.pinf); // XXX Check!
  // auto p0 = std::max(Fa.p0, Fc.p0);       // XXX Check!

  // For now, this returns. Later, might be faster to swap to in/out param!
  // (To avoid huge amount of re-alocating memory)
  // Actually, typically only need to call this once (for each a,b,c,d)
  // So, will be equally as fast with nRVO
  std::vector<double> Rabcd(kmax + 1, 0);

  const auto &ybd_kr = get_y_ijk(Fb, Fd);
  for (int k = kmin; k <= kmax; k++) {
    const auto &ybdk_r = ybd_kr[k - kmin];
    auto ffy = NumCalc::integrate(1.0, 0, pinf, Fa.f, Fc.f, ybdk_r, drdu);
    auto ggy = NumCalc::integrate(1.0, 0, pinf, Fa.g, Fc.g, ybdk_r, drdu);
    Rabcd[k] = (ffy + ggy) * du;
  }
  return Rabcd;
}

//******************************************************************************
//******************************************************************************
// Templates for faster method to calculate r^k
template <int k> static inline double powk_new(const double x) {
  return x * powk_new<k - 1>(x);
}
template <> inline double powk_new<0>(const double) { return 1.0; }

//******************************************************************************
template <int k>
static inline void yk_ijk(const int l, const DiracSpinor &Fa,
                          const DiracSpinor &Fb, std::vector<double> &vabk,
                          const std::size_t maxi)
// Calculalates y^k_ab screening function.
// Note: is symmetric: y_ab = y_ba
//
// Stores in vabk (in/out parameter, reference to whatever)
//
// r_min := min(r,r')
// rho(r') := fa(r')*fb(r') + ga(r')gb(r')
// y^k_ab(r) = Int_0^inf [r_min^k/r_max^(k+1)]*rho(f') dr'
//           = Int_0^r [r'^k/r^(k+1)]*rho(r') dr'
//             + Int_r^inf [r^k/r'^(k+1)]*rho(r') dr'
//          := A(r)/r^(k+1) + B(r)*r^k
// A(r0)  = 0
// B(r0)  = Int_0^inf [r^k/r'^(k+1)]*rho(r') dr'
// A(r_n) = A(r_{n-1}) + (rho(r_{n-1})*r_{n-1}^k)*dr
// B(r_n) = A(r_{n-1}) + (rho(r_{n-1})/r_{n-1}^(k+1))*dr
// y^k_ab(rn) = A(rn)/rn^(k+1) + B(rn)*rn^k
//
// Also uses Quadrature integer rules! (Defined in NumCalc)
{
  const auto &gr = Fa.p_rgrid; // just save typing
  const auto du = gr->du;
  const auto num_points = gr->num_points;
  vabk.resize(num_points); // for safety
  const auto irmax = (maxi == 0 || maxi > num_points) ? num_points : maxi;

  // faster method to calculate r^k
  auto powk = [=](double x) {
    if constexpr (k < 0)
      return std::pow(x, l);
    else
      return powk_new<k>(x);
  };

  // Quadrature integration weights:
  auto w = [=](std::size_t i) {
    if (i < NumCalc::Nquad)
      return NumCalc::dq_inv * NumCalc::cq[i];
    if (i < num_points - NumCalc::Nquad)
      return 1.0;
    return NumCalc::dq_inv * NumCalc::cq[num_points - i - 1];
  };

  double Ax = 0.0, Bx = 0.0;

  const auto bmax = std::min(Fa.pinf, Fb.pinf);
  for (std::size_t i = 0; i < bmax; i++) {
    Bx += w(i) * gr->drduor[i] * (Fa.f[i] * Fb.f[i] + Fa.g[i] * Fb.g[i]) /
          powk(gr->r[i]);
  }

  vabk[0] = Bx * du * powk(gr->r[0]);
  for (std::size_t i = 1; i < irmax; i++) {
    const auto rm1_to_k = powk(gr->r[i - 1]);
    const auto inv_rm1_to_kp1 = 1.0 / (rm1_to_k * gr->r[i - 1]);
    const auto r_to_k = powk(gr->r[i]);
    const auto inv_r_to_kp1 = 1.0 / (r_to_k * gr->r[i]);
    const auto Fdr = gr->drdu[i - 1] *
                     (Fa.f[i - 1] * Fb.f[i - 1] + Fa.g[i - 1] * Fb.g[i - 1]);
    const auto wi = w(i - 1);
    Ax += wi * Fdr * rm1_to_k;
    Bx -= wi * Fdr * inv_rm1_to_kp1;
    vabk[i] = du * (Ax * inv_r_to_kp1 + Bx * r_to_k);
  }
  for (std::size_t i = irmax; i < num_points; i++) {
    vabk[i] = 0.0;
  }
}

//******************************************************************************
void Coulomb::calculate_y_ijk(const DiracSpinor &Fa, const DiracSpinor &Fb,
                              const int k, std::vector<double> &vabk,
                              const std::size_t maxi)
// This is static
{
  auto sp1 = SafeProfiler::profile(__func__);

  if (k == 0)
    yk_ijk<0>(k, Fa, Fb, vabk, maxi);
  else if (k == 1)
    yk_ijk<1>(k, Fa, Fb, vabk, maxi);
  else if (k == 2)
    yk_ijk<2>(k, Fa, Fb, vabk, maxi);
  else if (k == 3)
    yk_ijk<3>(k, Fa, Fb, vabk, maxi);
  else if (k == 4)
    yk_ijk<4>(k, Fa, Fb, vabk, maxi);
  else if (k == 5)
    yk_ijk<5>(k, Fa, Fb, vabk, maxi);
  else if (k == 6)
    yk_ijk<6>(k, Fa, Fb, vabk, maxi);
  else if (k == 7)
    yk_ijk<7>(k, Fa, Fb, vabk, maxi);
  else if (k == 8)
    yk_ijk<8>(k, Fa, Fb, vabk, maxi);
  else
    yk_ijk<-1>(k, Fa, Fb, vabk, maxi);
}
