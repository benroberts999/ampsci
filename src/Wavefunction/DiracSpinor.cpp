#include "Wavefunction/DiracSpinor.hpp"
#include "Maths/Grid.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
#include "Physics/AtomData.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

constexpr bool update_pinf = true; // for psi += psi'
// XXX If true, sets all to num_points! [when damping orbitals!?]

//******************************************************************************
DiracSpinor::DiracSpinor(int in_n, int in_k, const Grid &rgrid)
    : p_rgrid(&rgrid), //
      n(in_n),
      k(in_k),
      en(0.0),                                       //
      f(std::vector<double>(rgrid.num_points, 0.0)), //
      g(f),                                          //
      p0(0),                                         //
      pinf(rgrid.num_points),                        //
      its(-1),
      eps(-1),
      occ_frac(0),                        //
      m_twoj(AtomData::twoj_k(in_k)),     //
      m_l(AtomData::l_k(in_k)),           //
      m_parity(AtomData::parity_k(in_k)), //
      m_k_index(AtomData::indexFromKappa(in_k)) {}

//******************************************************************************
std::string DiracSpinor::symbol(bool gnuplot) const {
  // Readable symbol (s_1/2, p_{3/2} etc.).
  // gnuplot-firndly '{}' braces optional.
  std::string ostring1 = (n != 0) ? std::to_string(n) + AtomData::l_symbol(m_l)
                                  : AtomData::l_symbol(m_l);
  std::string ostring2 = gnuplot ? "_{" + std::to_string(m_twoj) + "/2}"
                                 : "_" + std::to_string(m_twoj) + "/2";
  return ostring1 + ostring2;
}

std::string DiracSpinor::shortSymbol() const {
  std::string pm = (k < 0) ? "+" : "-";
  return n > 0 ? std::to_string(n) + AtomData::l_symbol(m_l) + pm
               : AtomData::l_symbol(m_l) + pm;
}

//******************************************************************************
double DiracSpinor::norm() const { return std::sqrt((*this) * (*this)); }

//******************************************************************************
const DiracSpinor &DiracSpinor::scale(const double factor) {
  for (std::size_t i = p0; i < pinf; ++i)
    f[i] *= factor;
  for (std::size_t i = p0; i < pinf; ++i)
    g[i] *= factor;
  // // XXX Need this for some reason!??
  // Means something beyond pinf is hapenning!?!? XXX XXX
  for (std::size_t i = pinf; i < f.size(); ++i) { // shouln't be needed!
    f[i] = 0;
    g[i] = 0;
  }
  for (std::size_t i = 0; i < p0; ++i) { // shouln't be needed!
    f[i] = 0;
    g[i] = 0;
  }
  // if (!df.empty()) {
  //   for (std::size_t i = p0; i < pinf; ++i) {
  //     df[i] *= factor;
  //     dg[i] *= factor;
  //   }
  // }
  return *this;
}
//------------------------------------------------------------------------------
const DiracSpinor &DiracSpinor::scale(const std::vector<double> &v) {
  const auto max = std::min(pinf, v.size());
  for (std::size_t i = p0; i < max; ++i) {
    f[i] *= v[i];
    g[i] *= v[i];
  }
  for (std::size_t i = max; i < pinf; ++i) {
    f[i] *= 0;
    g[i] *= 0;
  }
  // if (!df.empty()) {
  //   for (std::size_t i = p0; i < pinf; ++i) {
  //     df[i] *= v[i];
  //     dg[i] *= v[i];
  //   }
  // }
  return *this;
}

//******************************************************************************
void DiracSpinor::normalise(double norm_to) {
  double rescale_factor = norm_to / norm();
  scale(rescale_factor);
}

//******************************************************************************
std::pair<double, double> DiracSpinor::r0pinfratio() const {
  auto max_abs_compare = [](double a, double b) {
    return std::fabs(a) < std::fabs(b);
  };
  auto max_pos =
      std::max_element(f.begin(), f.begin() + long(pinf), max_abs_compare);
  auto r0_ratio = f[p0] / *max_pos;
  auto pinf_ratio = f[pinf - 1] / *max_pos;
  return std::make_pair(r0_ratio, pinf_ratio);
  // nb: do i care about ratio to max? or just value?
}

//******************************************************************************
std::vector<double> DiracSpinor::rho() const {
  std::vector<double> psi2;
  psi2.reserve(p_rgrid->num_points);
  const auto factor = twojp1() * occ_frac;
  for (auto i = 0ul; i < p_rgrid->num_points; ++i) {
    psi2.emplace_back(factor * (f[i] * f[i] + g[i] * g[i]));
  }
  return psi2;
}

//******************************************************************************
int DiracSpinor::num_electrons() const {
  return static_cast<int>(std::round((twoj() + 1) * occ_frac));
}

//******************************************************************************
//******************************************************************************
double operator*(const DiracSpinor &lhs, const DiracSpinor &rhs) {
  // Note: ONLY radial part ("F" radial spinor)
  const auto imin = std::max(lhs.p0, rhs.p0);
  const auto imax = std::min(lhs.pinf, rhs.pinf);
  const auto &dr = lhs.p_rgrid->drdu;
  const auto ff = NumCalc::integrate(1.0, imin, imax, lhs.f, rhs.f, dr);
  const auto gg = NumCalc::integrate(1.0, imin, imax, lhs.g, rhs.g, dr);
  return (ff + gg) * lhs.p_rgrid->du;
}

DiracSpinor &DiracSpinor::operator+=(const DiracSpinor &rhs) {
  // // XXX Here: pinf update_pinf
  if (this->k != rhs.k) {
    std::cerr
        << "Fail 100 in DiracSpinor: cannot add spinors of different kappa!\n";
    std::abort();
  }
  if (update_pinf) {
    pinf = std::max(pinf, rhs.pinf);
    p0 = std::min(p0, rhs.p0);
  }
  auto imax = std::min(pinf, rhs.pinf); // shouldn't be needed, but safer
  auto imin = std::max(p0, rhs.p0);     // shouldn't be needed, but safer
  for (std::size_t i = imin; i < imax; i++)
    f[i] += rhs.f[i];
  for (std::size_t i = imin; i < imax; i++)
    g[i] += rhs.g[i];
  // if (!df.empty() && !rhs.df.empty()) { // XXX Dangerous?
  //   for (std::size_t i = imin; i < imax; ++i) {
  //     df[i] += rhs.df[i];
  //     dg[i] += rhs.dg[i];
  //   }
  // }
  return *this;
}
DiracSpinor operator+(DiracSpinor lhs, const DiracSpinor &rhs) {
  lhs += rhs;
  return lhs;
}
DiracSpinor &DiracSpinor::operator-=(const DiracSpinor &rhs) {
  // XXX Here: pinf update_pinf
  if (this->k != rhs.k) {
    std::cerr
        << "Fail 121 in DiracSpinor: cannot add spinors of different kappa!\n";
    std::abort();
  }
  if (update_pinf) {
    pinf = std::max(pinf, rhs.pinf);
    p0 = std::min(p0, rhs.p0);
  }
  auto imax = std::min(pinf, rhs.pinf); // shouldn't be needed, but safer
  auto imin = std::max(p0, rhs.p0);     // shouldn't be needed, but safer
  for (std::size_t i = imin; i < imax; i++)
    f[i] -= rhs.f[i];
  for (std::size_t i = imin; i < imax; i++)
    g[i] -= rhs.g[i];
  // if (!df.empty() && !rhs.df.empty()) { // XXX Dangerous?
  //   for (std::size_t i = imin; i < imax; ++i) {
  //     df[i] -= rhs.df[i];
  //     dg[i] -= rhs.dg[i];
  //   }
  // }
  return *this;
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

DiracSpinor &DiracSpinor::operator=(const DiracSpinor &other) {
  if (*this != other) { // same n and kappa
    std::cerr
        << "Fail 152 in DiracSpinor: cannot assign spinors of different nk\n";
    std::abort();
  }
  if (this != &other) { // same memory location
    en = other.en;
    f = other.f;
    g = other.g;
    p0 = other.p0;
    pinf = other.pinf;
    occ_frac = other.occ_frac;
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
double DiracSpinor::r0() const { return p_rgrid->r[p0]; }
double DiracSpinor::rinf() const { return p_rgrid->r[pinf - 1]; }
