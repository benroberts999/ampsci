#include "Wavefunction/DiracSpinor.hpp"
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

//******************************************************************************
DiracSpinor::DiracSpinor(int in_n, int in_k,
                         std::shared_ptr<const Grid> in_rgrid)
    : rgrid(in_rgrid),
      n(in_n),
      k(in_k),
      m_f(std::vector<double>(rgrid ? rgrid->num_points() : 0, 0.0)),
      m_g(m_f),
      m_pinf(rgrid ? rgrid->num_points() : 0),
      m_twoj(AtomData::twoj_k(in_k)),
      m_l(AtomData::l_k(in_k)),
      m_parity(AtomData::parity_k(in_k)),
      m_k_index(AtomData::indexFromKappa(in_k)),
      m_nk_index(static_cast<Index>(AtomData::nk_to_index(in_n, in_k))) {}

//******************************************************************************
std::string DiracSpinor::symbol(bool gnuplot) const {
  // Readable symbol (s_1/2, p_{3/2} etc.).
  // gnuplot-firndly '{}' braces optional.
  std::string ostring1 = (n != 0) ?
                             std::to_string(n) + AtomData::l_symbol(m_l) :
                             AtomData::l_symbol(m_l);
  std::string ostring2 = gnuplot ? "_{" + std::to_string(m_twoj) + "/2}" :
                                   "_" + std::to_string(m_twoj) + "/2";
  return ostring1 + ostring2;
}

std::string DiracSpinor::shortSymbol() const {
  const std::string pm = (k < 0) ? "+" : "-";
  return (n != 0) ? std::to_string(n) + AtomData::l_symbol(m_l) + pm :
                    AtomData::l_symbol(m_l) + pm;
}

//******************************************************************************
double DiracSpinor::norm() const { return std::sqrt((*this) * (*this)); }

//******************************************************************************
const DiracSpinor &DiracSpinor::scale(const double factor) {
  for (std::size_t i = m_p0; i < m_pinf; ++i)
    m_f[i] *= factor;
  for (std::size_t i = m_p0; i < m_pinf; ++i)
    m_g[i] *= factor;
  // zero_boundaries(); // shouln't be needed!
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

//******************************************************************************
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

//******************************************************************************
std::pair<double, double> DiracSpinor::r0pinfratio() const {
  auto max_abs_compare = [](double a, double b) {
    return std::fabs(a) < std::fabs(b);
  };
  const auto max_pos = std::max_element(m_f.begin(), m_f.begin() + long(m_pinf),
                                        max_abs_compare);
  const auto r0_ratio = m_f[m_p0] / *max_pos;
  const auto pinf_ratio = m_f[m_pinf - 1] / *max_pos;
  return std::make_pair(r0_ratio, pinf_ratio);
  // nb: do i care about ratio to max? or just value?
}

//******************************************************************************
std::vector<double> DiracSpinor::rho() const {
  std::vector<double> psi2;
  psi2.reserve(rgrid->num_points());
  const auto factor = twojp1() * m_occ_frac;
  for (auto i = 0ul; i < rgrid->num_points(); ++i) {
    psi2.emplace_back(factor * (m_f[i] * m_f[i] + m_g[i] * m_g[i]));
  }
  return psi2;
}

//******************************************************************************
int DiracSpinor::num_electrons() const {
  return static_cast<int>(std::round((twoj() + 1) * m_occ_frac));
}

//******************************************************************************
double DiracSpinor::r0() const { return rgrid->r()[m_p0]; }
double DiracSpinor::rinf() const { return rgrid->r()[m_pinf - 1]; }

//******************************************************************************
//******************************************************************************
double operator*(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  // Note: ONLY radial part ("F" radial spinor)
  const auto imin = std::max(lhs.min_pt(), rhs.min_pt());
  const auto imax = std::min(lhs.max_pt(), rhs.max_pt());
  const auto &dr = lhs.rgrid->drdu();
  return (NumCalc::integrate(1, imin, imax, lhs.m_f, rhs.m_f, dr) +
          NumCalc::integrate(1, imin, imax, lhs.m_g, rhs.m_g, dr)) *
         lhs.rgrid->du();
}

DiracSpinor &DiracSpinor::operator+=(const DiracSpinor &rhs) {
  assert(this->k == rhs.k);

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
  assert(this->k == rhs.k);

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
  lhs += rhs;
  return lhs;
}
DiracSpinor operator-(DiracSpinor lhs, const DiracSpinor &rhs) {
  lhs -= rhs;
  return lhs;
}

DiracSpinor &DiracSpinor::operator*=(const double x) {
  scale(x);
  return *this;
}
DiracSpinor operator*(DiracSpinor lhs, const double x) {
  lhs *= x;
  return lhs;
}
DiracSpinor operator*(const double x, DiracSpinor rhs) {
  rhs *= x;
  return rhs;
}

DiracSpinor &DiracSpinor::operator*=(const std::vector<double> &v) {
  scale(v);
  return *this;
}
DiracSpinor operator*(const std::vector<double> &v, DiracSpinor rhs) {
  rhs *= v;
  return rhs;
}

//******************************************************************************
DiracSpinor &DiracSpinor::operator=(const DiracSpinor &other) {
  assert(*this == other); // same n and kappa
  if (this != &other) {   // same memory location
    m_en = other.m_en;
    m_f = other.m_f;
    m_g = other.m_g;
    m_p0 = other.min_pt();
    m_pinf = other.max_pt();
    m_occ_frac = other.m_occ_frac;
    m_its = other.m_its;
    m_eps = other.m_eps;
  }
  return *this;
}

//******************************************************************************
//******************************************************************************
// comparitor overloads:

bool operator==(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return lhs.n == rhs.n && lhs.k == rhs.k;
}

bool operator!=(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  return !(lhs == rhs);
}

bool operator<(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  if (lhs.n == rhs.n)
    return lhs.m_k_index < rhs.m_k_index;
  return lhs.n < rhs.n;
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
//******************************************************************************
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
      return (Fn.l() == l && Fn.n > max_n) ? Fn.n : max_n;
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

//******************************************************************************
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
int DiracSpinor::max_kindex(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty() ?
             0 :
             std::max_element(cbegin(orbs), cend(orbs), comp_ki)->k_index();
}

//******************************************************************************
std::pair<double, std::string>
DiracSpinor::check_ortho(const std::vector<DiracSpinor> &a,
                         const std::vector<DiracSpinor> &b) {
  double worst_del = 0.0;
  std::string worst_F = "";
  for (const auto &Fa : a) {
    for (const auto &Fb : b) {
      if (Fb.k != Fa.k)
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

//******************************************************************************
DiracSpinor DiracSpinor::exactHlike(int n, int k,
                                    std::shared_ptr<const Grid> rgrid,
                                    double zeff, double alpha) {
  if (alpha <= 0.0)
    alpha = PhysConst::alpha;
  DiracSpinor Fa(n, k, rgrid);
  using namespace DiracHydrogen;
  Fa.m_en = enk(PrincipalQN(n), DiracQN(k), Zeff(zeff), AlphaFS(alpha));
  for (std::size_t i = 0; i < rgrid->num_points(); ++i) {
    Fa.m_f[i] = DiracHydrogen::f(RaB(rgrid->r()[i]), PrincipalQN(n), DiracQN(k),
                                 Zeff(zeff), AlphaFS(alpha));
    Fa.m_g[i] = DiracHydrogen::g(RaB(rgrid->r()[i]), PrincipalQN(n), DiracQN(k),
                                 Zeff(zeff), AlphaFS(alpha));
  }
  return Fa;
}
