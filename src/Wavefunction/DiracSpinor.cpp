#include "Wavefunction/DiracSpinor.hpp"
#include "Angular/Wigner369j.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include "Physics/DiracHydrogen.hpp"
#include "Physics/PhysConst_constants.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
DiracSpinor::DiracSpinor(int in_n, int in_k,
                         std::shared_ptr<const Grid> in_rgrid)
  : m_rgrid(in_rgrid),
    m_n(in_n),
    m_kappa(in_k),
    m_f(std::vector<double>(m_rgrid ? m_rgrid->num_points() : 0, 0.0)),
    m_g(m_f),
    m_pinf(m_rgrid ? m_rgrid->num_points() : 0),
    m_twoj(Angular::twoj_k(in_k)),
    m_l(Angular::l_k(in_k)),
    m_parity(Angular::parity_k(in_k)),
    m_kappa_index(Angular::indexFromKappa(in_k)),
    m_nkappa_index(static_cast<Index>(Angular::nk_to_index(in_n, in_k))) {}

//==============================================================================
std::string DiracSpinor::symbol(bool gnuplot) const {
  // Readable symbol (s_1/2, p_{3/2} etc.).
  // gnuplot-firndly '{}' braces optional.
  std::string ostring1 = (m_n != 0) ?
                           std::to_string(m_n) + AtomData::l_symbol(m_l) :
                           AtomData::l_symbol(m_l);
  std::string ostring2 = gnuplot ? "_{" + std::to_string(m_twoj) + "/2}" :
                                   "_" + std::to_string(m_twoj) + "/2";
  return ostring1 + ostring2;
}

std::string DiracSpinor::shortSymbol() const {
  return shortSymbol(m_n, m_kappa);
}

std::string DiracSpinor::shortSymbol(int n, int kappa) {
  const std::string pm = (kappa < 0) ? "+" : "-";
  int l = Angular::l_k(kappa);
  return (n != 0) ? std::to_string(n) + AtomData::l_symbol(l) + pm :
                    AtomData::l_symbol(l) + pm;
}

//==============================================================================
void DiracSpinor::set_new_kappa(int new_kappa) {
  m_kappa = new_kappa;
  m_twoj = Angular::twoj_k(new_kappa);
  m_l = Angular::l_k(new_kappa);
  m_parity = Angular::parity_k(new_kappa);
  m_kappa_index = Angular::indexFromKappa(new_kappa);
  m_nkappa_index = static_cast<Index>(Angular::nk_to_index(m_n, new_kappa));
}

//==============================================================================
double DiracSpinor::norm() const { return std::sqrt((*this) * (*this)); }
double DiracSpinor::norm2() const { return (*this) * (*this); }

//==============================================================================
const DiracSpinor &DiracSpinor::scale(const double factor) {
  for (std::size_t i = m_p0; i < m_pinf; ++i)
    m_f[i] *= factor;
  for (std::size_t i = m_p0; i < m_pinf; ++i)
    m_g[i] *= factor;
  return *this;
}
//------------------------------------------------------------------------------
const DiracSpinor &DiracSpinor::scale(const std::vector<double> &v) {
  const auto max = std::min(m_pinf, v.size());
  for (std::size_t i = m_p0; i < max; ++i) {
    m_f[i] *= v[i];
    m_g[i] *= v[i];
  }
  for (std::size_t i = max; i < m_pinf; ++i) {
    m_f[i] = 0;
    m_g[i] = 0;
  }
  return *this;
}

//==============================================================================
void DiracSpinor::normalise(double norm_to) { scale(norm_to / norm()); }

//------------------------------------------------------------------------------
void DiracSpinor::zero_boundaries() {
  for (std::size_t i = 0; i < m_p0; ++i) {
    m_f[i] = 0.0;
    m_g[i] = 0.0;
  }
  for (std::size_t i = m_pinf; i < m_f.size(); ++i) {
    m_f[i] = 0.0;
    m_g[i] = 0.0;
  }
}

//==============================================================================
std::pair<double, double> DiracSpinor::r0pinfratio() const {
  auto max_abs_compare = [](double a, double b) {
    return std::fabs(a) < std::fabs(b);
  };
  const auto max_pos =
    std::max_element(m_f.begin(), m_f.begin() + long(m_pinf), max_abs_compare);
  const auto r0_ratio = m_f[m_p0] / *max_pos;
  const auto pinf_ratio = m_f[m_pinf - 1] / *max_pos;
  return std::make_pair(r0_ratio, pinf_ratio);
  // nb: do i care about ratio to max? or just value?
}

//==============================================================================
std::vector<double> DiracSpinor::rho() const {
  std::vector<double> psi2;
  psi2.reserve(m_rgrid->num_points());
  const auto factor = twojp1() * m_occ_frac;
  for (auto i = 0ul; i < m_rgrid->num_points(); ++i) {
    psi2.emplace_back(factor * (m_f[i] * m_f[i] + m_g[i] * m_g[i]));
  }
  return psi2;
}

//==============================================================================
int DiracSpinor::num_electrons() const {
  return static_cast<int>(std::round((twoj() + 1) * m_occ_frac));
}

//==============================================================================
double DiracSpinor::r0() const { return m_rgrid->r()[m_p0]; }
double DiracSpinor::rinf() const { return m_rgrid->r()[m_pinf - 1]; }

//==============================================================================
//==============================================================================
double operator*(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  // Note: ONLY radial part ("F" radial spinor)
  const auto imin = std::max(lhs.min_pt(), rhs.min_pt());
  const auto imax = std::min(lhs.max_pt(), rhs.max_pt());
  const auto &dr = lhs.m_rgrid->drdu();
  return (NumCalc::integrate(1, imin, imax, lhs.m_f, rhs.m_f, dr) +
          NumCalc::integrate(1, imin, imax, lhs.m_g, rhs.m_g, dr)) *
         lhs.m_rgrid->du();
}

DiracSpinor &DiracSpinor::operator+=(const DiracSpinor &rhs) {
  // assert(this->m_kappa == rhs.m_kappa);

  if (rhs.max_pt() > m_pinf)
    m_pinf = rhs.max_pt();
  if (rhs.min_pt() < m_p0)
    m_p0 = rhs.min_pt();

  for (std::size_t i = rhs.min_pt(); i < rhs.max_pt(); i++)
    m_f[i] += rhs.m_f[i];
  for (std::size_t i = rhs.min_pt(); i < rhs.max_pt(); i++)
    m_g[i] += rhs.m_g[i];
  return *this;
}
DiracSpinor &DiracSpinor::operator-=(const DiracSpinor &rhs) {
  // assert(this->m_kappa == rhs.m_kappa);

  if (rhs.max_pt() > m_pinf)
    m_pinf = rhs.max_pt();
  if (rhs.min_pt() < m_p0)
    m_p0 = rhs.min_pt();

  for (std::size_t i = rhs.min_pt(); i < rhs.max_pt(); i++)
    m_f[i] -= rhs.m_f[i];
  for (std::size_t i = rhs.min_pt(); i < rhs.max_pt(); i++)
    m_g[i] -= rhs.m_g[i];

  return *this;
}

DiracSpinor operator+(DiracSpinor lhs, const DiracSpinor &rhs) {
  return lhs += rhs;
}
DiracSpinor operator-(DiracSpinor lhs, const DiracSpinor &rhs) {
  return lhs -= rhs;
}

DiracSpinor &DiracSpinor::operator*=(const double x) {
  scale(x);
  return *this;
}
DiracSpinor operator*(DiracSpinor lhs, const double x) { return lhs *= x; }
DiracSpinor operator*(const double x, DiracSpinor rhs) { return rhs *= x; }

DiracSpinor &DiracSpinor::operator*=(const std::vector<double> &v) {
  scale(v);
  return *this;
}
DiracSpinor operator*(const std::vector<double> &v, DiracSpinor rhs) {
  return rhs *= v;
}

//==============================================================================
//==============================================================================
// comparitor overloads:

bool operator==(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return lhs.n() == rhs.n() && lhs.kappa() == rhs.kappa();
}

bool operator!=(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return !(lhs == rhs);
}

bool operator<(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  if (lhs.n() == rhs.n())
    return lhs.m_kappa_index < rhs.m_kappa_index;
  return lhs.n() < rhs.n();
}

bool operator>(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return rhs < lhs;
}

bool operator<=(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return !(lhs > rhs);
}

bool operator>=(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return !(lhs < rhs);
}
//==============================================================================
// static
std::string DiracSpinor::state_config(const std::vector<DiracSpinor> &orbs) {
  std::string result = "";
  if (orbs.empty())
    return result;

  // find max l
  const auto maxl =
    std::max_element(orbs.cbegin(), orbs.cend(), DiracSpinor::comp_l)->l();

  // for each l, count num, add to string
  int prev_max_n = 0;
  for (int l = 0; l <= maxl; ++l) {

    auto find_max_n_given_l = [l](int max_n, const auto &Fn) {
      return (Fn.l() == l && Fn.n() > max_n) ? Fn.n() : max_n;
    };
    const auto max_n =
      std::accumulate(orbs.cbegin(), orbs.cend(), 0, find_max_n_given_l);

    // format 'state string' into required notation:
    if (max_n == prev_max_n && max_n != 0)
      result += AtomData::l_symbol(l);
    else if (max_n != 0)
      result += std::to_string(max_n) + AtomData::l_symbol(l);

    prev_max_n = max_n;
  }

  return result;
}

//==============================================================================
// static
int DiracSpinor::max_tj(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty() ?
           0 :
           std::max_element(cbegin(orbs), cend(orbs), comp_j)->twoj();
}
// static
int DiracSpinor::max_l(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty() ? 0 :
                        std::max_element(cbegin(orbs), cend(orbs), comp_l)->l();
}
// static
int DiracSpinor::max_n(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty() ? 0 :
                        std::max_element(cbegin(orbs), cend(orbs), comp_n)->n();
}
// static
int DiracSpinor::max_kindex(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty() ?
           0 :
           std::max_element(cbegin(orbs), cend(orbs), comp_ki)->k_index();
}

// static
const DiracSpinor *DiracSpinor::find(int n, int k,
                                     const std::vector<DiracSpinor> &orbs) {
  const auto find_nk = [n, k](const auto Fa) {
    return Fa.n() == n && Fa.kappa() == k;
  };
  auto Fnk = std::find_if(cbegin(orbs), cend(orbs), find_nk);
  if (Fnk != cend(orbs)) {
    return &*Fnk;
  }
  // otherwise, return nope
  return nullptr;
}

//==============================================================================
std::pair<double, std::string>
DiracSpinor::check_ortho(const std::vector<DiracSpinor> &a,
                         const std::vector<DiracSpinor> &b) {
  double worst_del = 0.0;
  std::string worst_F = "";
  for (const auto &Fa : a) {
    for (const auto &Fb : b) {
      if (Fb.kappa() != Fa.kappa())
        continue;
      const auto del =
        Fa == Fb ? std::abs(std::abs(Fa * Fb) - 1.0) : std::abs(Fa * Fb);
      // nb: sometimes sign of Fb is wrong. Perhaps this is an issue??
      if (del > worst_del) {
        worst_del = del;
        worst_F = "<" + Fa.shortSymbol() + "|" + Fb.shortSymbol() + ">";
      }
    }
  }
  return {worst_del, worst_F};
}

//==============================================================================
void DiracSpinor::orthog(const DiracSpinor &rhs) {
  if (rhs.kappa() != m_kappa)
    return;
  const auto k = *this * rhs;
  const auto p0 = std::max(m_p0, rhs.m_p0);
  const auto pi = std::min(m_pinf, rhs.m_pinf);
  for (std::size_t i = p0; i < pi; ++i) {
    m_f[i] -= k * rhs.m_f[i];
    m_g[i] -= k * rhs.m_g[i];
  }
}

//==============================================================================
void DiracSpinor::orthonormaliseOrbitals(std::vector<DiracSpinor> &in_orbs,
                                         int num_its)
// Note: this function is static
// Forces ALL orbitals to be orthogonal to each other, and normal
// Note: workes best if run twice!
// |a> ->  |a> - \sum_{b!=a} |b><b|a>
// Then:
// |a> -> |a> / <a|a>
// c_ba = c_ab = <a|b>
// num_its is optional parameter. Repeats that many times!
// Note: I force all orthog to each other - i.e. double count.
// {force <2|1>=0 and then <1|2>=0}
// Would be 2x faster not to do this
//  - but that would treat some orbitals special!
// Hence factor of 0.5
// Note: For HF, should never be called after core is frozen!
//
// Note: This allows wfs to extend past pinf!
// ==> This causes the possible orthog issues..
{
  const auto Ns = in_orbs.size();

  // Calculate c_ab = <a|b>  [only for b>a -- symmetric]
  std::vector<std::vector<double>> c_ab(Ns, std::vector<double>(Ns));
  for (std::size_t a = 0; a < Ns; a++) {
    const auto &Fa = in_orbs[a];
    for (auto b = a + 1; b < Ns; b++) {
      const auto &Fb = in_orbs[b];
      if (Fa.kappa() != Fb.kappa())
        continue;
      c_ab[a][b] = 0.5 * (Fa * Fb);
    }
  }
  // note: above loop executes psia*psib half as many times as below would

  // Orthogonalise + re-norm orbitals:
  for (std::size_t a = 0; a < Ns; a++) {
    auto &Fa = in_orbs[a];
    for (std::size_t b = 0; b < Ns; b++) {
      const auto &Fb = in_orbs[b];
      if (Fa.kappa() != Fb.kappa() || Fa.n() == Fb.n())
        continue;
      double cab = (a < b) ? c_ab[a][b] : c_ab[b][a];
      Fa -= cab * Fb;
    }
    Fa.normalise();
  }

  // If necisary: repeat
  if (num_its > 1)
    orthonormaliseOrbitals(in_orbs, num_its - 1);
}

//==============================================================================
DiracSpinor
DiracSpinor::orthogonaliseWrt(const DiracSpinor &Fin,
                              const std::vector<DiracSpinor> &orbs) {
  auto Fv = Fin; //psi_v may be in in_orbs
  for (const auto &Fn : orbs) {
    if (Fv.kappa() != Fn.kappa())
      continue;
    if (Fv == Fn) {
      Fv = Fn; // also copies energy
      break;
    } else {
      Fv -= (Fin * Fn) * Fn; // does not update energy
    }
  }
  return Fv;
}
//------------------------------------------------------------------------------
DiracSpinor
DiracSpinor::orthonormaliseWrt(const DiracSpinor &psi_v,
                               const std::vector<DiracSpinor> &in_orbs) {
  auto Fv = orthogonaliseWrt(psi_v, in_orbs);
  Fv.normalise();
  return Fv;
}

//==============================================================================
DiracSpinor DiracSpinor::exactHlike(int n, int kappa,
                                    std::shared_ptr<const Grid> rgrid,
                                    double zeff, double alpha) {
  if (alpha <= 0.0)
    alpha = PhysConst::alpha;
  DiracSpinor Fa(n, kappa, rgrid);
  using namespace DiracHydrogen;
  Fa.m_en = enk(PrincipalQN(n), DiracQN(kappa), Zeff(zeff), AlphaFS(alpha));
  for (std::size_t i = 0; i < rgrid->num_points(); ++i) {
    Fa.m_f[i] = DiracHydrogen::f(RaB(rgrid->r(i)), PrincipalQN(n),
                                 DiracQN(kappa), Zeff(zeff), AlphaFS(alpha));
    Fa.m_g[i] = DiracHydrogen::g(RaB(rgrid->r(i)), PrincipalQN(n),
                                 DiracQN(kappa), Zeff(zeff), AlphaFS(alpha));
  }
  return Fa;
}

//==============================================================================
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>>
DiracSpinor::split_by_energy(const std::vector<DiracSpinor> &orbitals,
                             double energy, int n_min_core) {
  std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> out;
  auto &[below, above] = out;
  for (const auto &n : orbitals) {
    if (n.en() <= energy) {
      if (n.n() >= n_min_core) {
        below.push_back(n);
      }
    } else {
      above.push_back(n);
    }
  }
  return out;
}

//==============================================================================
std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>>
DiracSpinor::split_by_core(const std::vector<DiracSpinor> &orbitals,
                           const std::vector<DiracSpinor> &core,
                           int n_min_core) {
  std::pair<std::vector<DiracSpinor>, std::vector<DiracSpinor>> out;
  auto &[below, above] = out;
  for (const auto &n : orbitals) {
    bool in_core = std::find(core.cbegin(), core.cend(), n) != core.cend();
    if (in_core) {
      if (n.n() >= n_min_core) {
        below.push_back(n);
      }
    } else {
      above.push_back(n);
    }
  }
  return out;
}

//==============================================================================
std::vector<DiracSpinor>
DiracSpinor::subset(const std::vector<DiracSpinor> &basis,
                    const std::string &subset_string) {

  // Form 'subset' from {a} in 'basis', if:
  //    a _is_ in subset_string AND
  //    a _is not_ in exclude_string

  std::vector<DiracSpinor> subset;
  const auto nmaxk_list = AtomData::n_kappa_list(subset_string);

  for (const auto &a : basis) {

    // Check if a is present in 'subset_string'
    const auto nk =
      std::find_if(nmaxk_list.cbegin(), nmaxk_list.cend(),
                   [&a](const auto &tnk) { return a.kappa() == tnk.second; });
    if (nk == nmaxk_list.cend())
      continue;
    // nk is now max n, for given kappa {max_n, kappa}
    if (a.n() > nk->first)
      continue;

    subset.push_back(a);
  }
  return subset;
}
