#include "Anderson.hpp"
#include "LinAlg/Matrix.hpp"
#include "LinAlg/Solvers.hpp"
#include "LinAlg/Vector.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

namespace ExternalField {

//==============================================================================
std::vector<double>
anderson_coefficients(const std::vector<std::vector<double>> &gram) {
  // Bordered system [B 1; 1^T 0][c; lambda] = [0; 1] over the kept entries.
  // Drops the oldest entries while B is ill-conditioned (saturated history)
  // or the solve fails; there are no exceptions in this codebase, so a
  // failed LU is detected from non-finite coefficients.
  constexpr double cond_floor = 1.0e-14;
  const auto n = gram.size();
  for (auto m = n; m > 0; --m) {
    const auto i0 = n - m; // oldest kept entry
    LinAlg::Matrix<double> B(m + 1, m + 1);
    LinAlg::Vector<double> rhs(m + 1);
    double dmax = 0.0, dmin = 1.0e300;
    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < m; ++j) {
        B(i, j) = gram[i0 + i][i0 + j];
      }
      B(i, m) = 1.0;
      B(m, i) = 1.0;
      rhs(i) = 0.0;
      dmax = std::max(dmax, B(i, i));
      dmin = std::min(dmin, B(i, i));
    }
    B(m, m) = 0.0;
    rhs(m) = 1.0;
    if (m > 2 && dmin < cond_floor * dmax)
      continue;
    const auto c = LinAlg::solve_Axeqb<double>(B, rhs);
    std::vector<double> coefs(m);
    bool finite = true;
    for (std::size_t i = 0; i < m; ++i) {
      coefs[i] = c(i);
      finite = finite && std::isfinite(coefs[i]);
    }
    if (finite)
      return coefs;
  }
  return {};
}

//==============================================================================
std::vector<double>
anderson_coefficients(const std::vector<DiracSpinor> &residuals) {
  const auto n = residuals.size();
  std::vector<std::vector<double>> gram(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      gram[i][j] = residuals[i] * residuals[j];
      gram[j][i] = gram[i][j];
    }
  }
  return anderson_coefficients(gram);
}

//==============================================================================
void AndersonHistory::push(DiracSpinor g, DiracSpinor r, double r2, double k) {
  if (m_g.size() >= m_depth) {
    drop_oldest(m_g.size() + 1 - m_depth);
  }
  // New row and column of the Gram matrix: one inner product per kept entry
  const auto n = m_r.size();
  std::vector<double> row(n + 1);
  for (std::size_t i = 0; i < n; ++i) {
    row[i] = m_r[i] * r;
    m_gram[i].push_back(row[i]);
  }
  row[n] = r2;
  m_gram.push_back(std::move(row));
  m_g.push_back(std::move(g));
  m_r.push_back(std::move(r));
  m_k.push_back(k);
}

//------------------------------------------------------------------------------
void AndersonHistory::drop_oldest(std::size_t n) {
  n = std::min(n, m_g.size());
  if (n == 0)
    return;
  const auto ln = long(n);
  m_g.erase(m_g.begin(), m_g.begin() + ln);
  m_r.erase(m_r.begin(), m_r.begin() + ln);
  m_k.erase(m_k.begin(), m_k.begin() + ln);
  m_gram.erase(m_gram.begin(), m_gram.begin() + ln);
  for (auto &row : m_gram) {
    row.erase(row.begin(), row.begin() + ln);
  }
}

//------------------------------------------------------------------------------
void AndersonHistory::mix(DiracSpinor *x, double *k) {
  assert(x != nullptr && !m_g.empty());
  const auto c = anderson_coefficients(m_gram);
  if (c.empty()) {
    // No usable entry: the plain step, and the history starts afresh
    *x = std::move(m_g.back());
    if (k) {
      *k = m_k.back();
    }
    drop_oldest(m_g.size());
    return;
  }
  drop_oldest(m_g.size() - c.size());
  *x = m_g[0];
  *x *= c[0];
  double k_mixed = c[0] * m_k[0];
  for (std::size_t i = 1; i < c.size(); ++i) {
    x->add_scaled(c[i], m_g[i]);
    k_mixed += c[i] * m_k[i];
  }
  if (k) {
    *k = k_mixed;
  }
}

//==============================================================================
void AndersonVectorHistory::push(const std::vector<double> &x,
                                 std::vector<double> g) {
  assert(x.size() == g.size());
  if (m_g.size() >= m_depth) {
    drop_oldest(m_g.size() + 1 - m_depth);
  }
  auto r = g;
  for (std::size_t k = 0; k < r.size(); ++k) {
    r[k] -= x[k];
  }
  // New row and column of the Gram matrix: one inner product per kept entry
  const auto n = m_r.size();
  std::vector<double> row(n + 1);
  for (std::size_t i = 0; i < n; ++i) {
    row[i] =
      std::inner_product(m_r[i].cbegin(), m_r[i].cend(), r.cbegin(), 0.0);
    m_gram[i].push_back(row[i]);
  }
  row[n] = std::inner_product(r.cbegin(), r.cend(), r.cbegin(), 0.0);
  m_gram.push_back(std::move(row));
  m_g.push_back(std::move(g));
  m_r.push_back(std::move(r));
}

//------------------------------------------------------------------------------
void AndersonVectorHistory::drop_oldest(std::size_t n) {
  n = std::min(n, m_g.size());
  if (n == 0)
    return;
  const auto ln = long(n);
  m_g.erase(m_g.begin(), m_g.begin() + ln);
  m_r.erase(m_r.begin(), m_r.begin() + ln);
  m_gram.erase(m_gram.begin(), m_gram.begin() + ln);
  for (auto &row : m_gram) {
    row.erase(row.begin(), row.begin() + ln);
  }
}

//------------------------------------------------------------------------------
std::vector<double> AndersonVectorHistory::mix() {
  assert(!m_g.empty());
  const auto c = anderson_coefficients(m_gram);
  if (c.empty()) {
    // No usable entry: the plain step, and the history starts afresh
    auto x = std::move(m_g.back());
    drop_oldest(m_g.size());
    return x;
  }
  drop_oldest(m_g.size() - c.size());
  std::vector<double> x(m_g[0].size(), 0.0);
  for (std::size_t i = 0; i < c.size(); ++i) {
    for (std::size_t k = 0; k < x.size(); ++k) {
      x[k] += c[i] * m_g[i][k];
    }
  }
  return x;
}

} // namespace ExternalField
