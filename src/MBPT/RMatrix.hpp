#pragma once
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <iostream>
#include <type_traits>

// XXX NOTE: with multiplication defined as int... identity is NOT 1!

namespace MBPT {

/*! Defines RMatrix, Radial matrix. Designed to store radial coordinate operators 
(like the Coulomb and polarisation operators). The matrix is stored on a sub-grid (between r0
and rmax), with a specified stride.

@details

*/

template <typename T>
class RMatrix {

  std::size_t m_i0, m_stride;
  std::size_t m_size;
  std::size_t m_g_size;
  LinAlg::Matrix<T> m_Rmatrix;
  std::shared_ptr<const Grid> m_rgrid; // "full" grid
  std::vector<double> sub_r{};         // sub grid

public:
  //============================================================================

  RMatrix(std::size_t i0, std::size_t stride, std::size_t size,
          std::shared_ptr<const Grid> rgrid)
      : m_i0(i0),
        m_stride(stride),
        m_size(size),
        m_Rmatrix(m_size),
        m_rgrid(rgrid) {
    //------------------
    // create vector of r on sub-grid, used to interpolate values onto full
    const auto &r = m_rgrid->r();
    sub_r.reserve(m_size);
    assert(m_i0 + m_stride * m_size <= r.size());
    for (std::size_t i = 0; i < m_size; ++i) {
      sub_r.push_back(r[index_to_fullgrid(i)]);
    }
    assert(m_i0 + m_stride * m_size <= r.size());
    assert(sub_r[1] == r[index_to_fullgrid(1)]);
    assert(sub_r[m_size - 1] == r[index_to_fullgrid(m_size - 1)]);
    //------------------
  }

  //============================================================================
  //! direct access to matrix elements
  T &operator()(std::size_t i, std::size_t j) { return m_Rmatrix(i, j); }
  const T operator()(std::size_t i, std::size_t j) const {
    return m_Rmatrix(i, j);
  }

  //! direct access to radial matrix
  const LinAlg::Matrix<T> &Rmatrix() const { return m_Rmatrix; }
  LinAlg::Matrix<T> &Rmatrix() { return m_Rmatrix; }

  std::size_t size() const { return m_size; }
  std::size_t i0() const { return m_i0; }
  std::size_t stride() const { return m_stride; }

  //============================================================================
  //! Sets the radial matrix to have all zero's
  void zero() { m_Rmatrix.zero(); }

  // Makes 1
  [[deprecated]] void make_identity() { m_Rmatrix.make_identity(); }

  // only for transition, kill!
  [[deprecated]] RMatrix<T> &plusIdent(T a = T{1.0}) {
    (*this) += a;
    return *this;
  }

  //============================================================================
  //! Matrix adition +,-
  RMatrix<T> &operator+=(const RMatrix<T> &rhs) {
    m_Rmatrix += rhs.m_Rmatrix;
    return *this;
  }
  //! Matrix adition +,-
  RMatrix<T> &operator-=(const RMatrix<T> &rhs) {
    m_Rmatrix -= rhs.m_Rmatrix;
    return *this;
  }
  //! Scalar multiplication
  RMatrix<T> &operator*=(const T x) {
    m_Rmatrix *= x;
    return *this;
  }

  //! Matrix adition +,-
  [[nodiscard]] friend RMatrix<T> operator+(RMatrix<T> lhs,
                                            const RMatrix<T> &rhs) {
    return (lhs += rhs);
  }
  //! Matrix adition +,-
  [[nodiscard]] friend RMatrix<T> operator-(RMatrix<T> lhs,
                                            const RMatrix<T> &rhs) {
    return (lhs -= rhs);
  }
  //! Scalar multiplication
  [[nodiscard]] friend RMatrix<T> operator*(const T x, RMatrix<T> rhs) {
    return (rhs *= x);
  }

  //! Adition of identity: Matrix<T> += T : T assumed to be *Identity!
  RMatrix<T> &operator+=(T aI) {
    m_Rmatrix += aI;
    return *this;
  }
  //! Adition of identity: Matrix<T> -= T : T assumed to be *Identity!
  RMatrix<T> &operator-=(T aI) {
    m_Rmatrix -= aI;
    return *this;
  }

  //! Adition of identity: Matrix<T> + T : T assumed to be *Identity!
  [[nodiscard]] friend RMatrix<T> operator+(RMatrix<T> M, T aI) {
    return (M += aI);
  }
  //! Adition of identity: Matrix<T> - T : T assumed to be *Identity!
  [[nodiscard]] friend RMatrix<T> operator-(RMatrix<T> M, T aI) {
    return (M -= aI);
  }

  //============================================================================

  //! Matrix multplication: C=A*B := Cij = \sum_k Aik*Bkj.
  //! Note: integration measure not included: call .drj() first to include it!
  [[nodiscard]] friend RMatrix<T> operator*(const RMatrix<T> &a,
                                            const RMatrix<T> &b) {

    RMatrix<T> out(a.m_i0, a.m_stride, a.m_size, a.m_rgrid);

    out.Rmatrix() = a.Rmatrix() * b.Rmatrix();

    return out;
  }

  //============================================================================
  //! Multiply coordinate elements (in place): Gij -> Gij*Bij
  RMatrix<T> &mult_elements_by(const RMatrix<T> &rhs) {
    m_Rmatrix.mult_elements_by(rhs.ff());
    return *this;
  }
  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[nodiscard]] friend RMatrix<T> mult_elements(RMatrix<T> lhs,
                                                const RMatrix<T> &rhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //============================================================================

  //! Returns conjugate of matrix
  [[nodiscard]] RMatrix<T> conj() const {
    auto out = *this;
    out.Rmatrix().conj_in_place();
    return out;
  }
  //! Returns real part of complex matrix (changes type; returns a real
  //! matrix)
  [[nodiscard]] RMatrix<double> real() const {
    RMatrix<double> out(m_i0, m_stride, m_size, m_rgrid);
    out.Rmatrix() = m_Rmatrix.real();
    return out;
  }
  //! Returns imag part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] RMatrix<double> imag() const {
    RMatrix<double> out(m_i0, m_stride, m_size, m_rgrid);
    out.Rmatrix() = m_Rmatrix.imag();
    return out;
  }
  //! Converts a real to complex matrix (changes type; returns a complex
  //! matrix)
  [[nodiscard]] RMatrix<std::complex<double>> complex() const {
    RMatrix<std::complex<double>> out(m_i0, m_stride, m_size, m_rgrid);
    out.Rmatrix() = m_Rmatrix.complex();
    return out;
  }

  //============================================================================
  //! Inversion (in place)
  RMatrix<T> &invert_in_place() {
    m_Rmatrix.invert_in_place();
    return *this;
  }
  //! Returns inverse of matrix; original matrix unchanged
  [[nodiscard]] RMatrix<T> inverse() const {
    auto out = *this; //
    return out.invert_in_place();
  }

  //============================================================================
  //! Multiplies by drj: Q_ij -> Q_ij*dr_j, in place
  RMatrix<T> &drj_in_place() {
    const auto dus = m_rgrid->du() * double(m_stride);
    for (auto i = 0ul; i < m_size; ++i) {
      for (auto j = 0ul; j < m_size; ++j) {
        const auto sj = index_to_fullgrid(j);
        const auto dr = m_rgrid->drdu(sj) * dus;
        m_Rmatrix[i][j] *= dr;
      }
    }
    return *this;
  }
  //! Multiplies by dri: Q_ij -> Q_ij*dr_i, in place
  RMatrix<T> &dri_in_place() {
    const auto dus = m_rgrid->du() * double(m_stride);
    for (auto i = 0ul; i < m_size; ++i) {
      const auto si = index_to_fullgrid(i);
      const auto dr = m_rgrid->drdu(si) * dus;
      for (auto j = 0ul; j < m_size; ++j) {
        m_Rmatrix[i][j] *= dr;
      }
    }
    return *this;
  }
  //! Multiplies by drj: Q_ij -> Q_ij*dr_j. Returns new matrix (orig unchanged)
  RMatrix<T> drj() const {
    auto out = *this;
    return out.drj_in_place();
  }
  //! Multiplies by dri: Q_ij -> Q_ij*dr_i. Returns new matrix (orig unchanged)
  RMatrix<T> dri() const {
    auto out = *this;
    return out.dri_in_place();
  }

  //! returns dr at position along sub grid
  double dr(std::size_t sub_index) const {
    const auto full_index = index_to_fullgrid(sub_index);
    return m_rgrid->drdu(full_index) * m_rgrid->du() * double(m_stride);
  }

  //============================================================================
  //! Converts an index on the sub-grid to the full grid.
  std::size_t index_to_fullgrid(std::size_t i) const {
    return m_i0 + i * m_stride;
  }

  //============================================================================

  /*
  //! Adds k*|ket><bra| to matrix (used for building Green's functions)
  void add(const DiracSpinor &ket, const DiracSpinor &bra, T k = T(1.0)) {
    // Adds (k)*|ket><bra| to G matrix
    // G_ij = f * Q_i * W_j
    // Q = Q(1) = ket, W = W(2) = bra
    // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
    for (auto i = 0ul; i < m_size; ++i) {
      const auto si = index_to_fullgrid(i);
      for (auto j = 0ul; j < m_size; ++j) {
        const auto sj = index_to_fullgrid(j);
        m_ff[i][j] += k * ket.f(si) * bra.f(sj);
      } // j
    } // i

    if (m_incl_g) {
      for (auto i = 0ul; i < m_size; ++i) {
        const auto si = index_to_fullgrid(i);
        for (auto j = 0ul; j < m_size; ++j) {
          const auto sj = index_to_fullgrid(j);
          // XXX Double check fg/gf right way!
          m_fg[i][j] += k * ket.f(si) * bra.g(sj);
          m_gf[i][j] += k * ket.g(si) * bra.f(sj); // symmetric, transpose?
          m_gg[i][j] += k * ket.g(si) * bra.g(sj);
        } // j
      } // i
    }
  }

  //============================================================================
  //! Action of RDMatrix operator on DiracSpinor. Inludes Integration:
  //! G*F = Int[G(r,r')*F(r') dr'] = sum_j G_ij*F_j*drdu_j*du
  DiracSpinor operator*(const DiracSpinor &Fn) const {

    const auto &r = Fn.grid().r();
    // const auto &drdu = Fn.grid().drdu();
    // const double s_du = double(m_stride) * Fn.grid().du();

    // include dr?? No, not by default?
    std::vector<double> f(m_size), g;
    for (auto i = 0ul; i < m_size; ++i) {
      for (auto j = 0ul; j < m_size; ++j) {
        const auto j_f = index_to_fullgrid(j);
        f[i] += m_ff(i, j) * Fn.f(j_f); // * drdu[j_f] * s_du;
      }
    }
    if (m_incl_g) {
      g.resize(m_size);
      for (auto i = 0ul; i < m_size; ++i) {
        for (auto j = 0ul; j < m_size; ++j) {
          const auto j_f = index_to_fullgrid(j);
          f[i] += m_fg(i, j) * Fn.g(j_f); // * drdu[j_f] * s_du;
          g[i] += (m_gf(i, j) * Fn.f(j_f) + m_gg(i, j) * Fn.g(j_f)); // *
          // drdu[j_f] * s_du;
        }
      }
    }

    DiracSpinor out = Fn * 0.0;
    // Interpolate from sub-grid to full grid
    out.f() = Interpolator::interpolate(sub_r, f, r);
    if (m_incl_g) {
      out.g() = Interpolator::interpolate(sub_r, g, r);
    }
    return out;
  }

  //============================================================================
  // For testing only:
  friend std::ostream &operator<<(std::ostream &os, const RDMatrix<T> &a) {
    os << "FF:\n";
    os << a.m_ff;
    if (a.m_incl_g) {
      os << "FG:\n";
      os << a.m_fg;
      os << "GF:\n";
      os << a.m_gf;
      os << "GG:\n";
      os << a.m_gg;
    }
    return os;
  }
*/
};

//! Checks if two matrix's are equal (to within parts in 10^12)
template <typename T>
bool equal(const RMatrix<T> &lhs, const RMatrix<T> &rhs) {
  return equal(lhs.Rmatrix(), rhs.Rmatrix());
}

//! returns maximum element (by abs)
template <typename T>
double max_element(const RMatrix<T> &a) {
  double xmax = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      if (std::abs(a.ff(i, j)) > xmax) {
        xmax = std::abs(a.m_Rmatrix(i, j));
      }
    }
  }
  return xmax;
}

//! returns maximum difference (abs) between two matrixs
template <typename T>
double max_delta(const RMatrix<T> &a, const RMatrix<T> &b) {
  double delta = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      if (std::abs(a.Rmatrix(i, j) - b.Rmatrix(i, j)) > delta) {
        delta = std::abs(a.Rmatrix(i, j) - b.Rmatrix(i, j));
      }
    }
  }
  return delta;
}

//! returns maximum relative diference [aij-bij/(aij+bij)] (abs) between two
//! matrices
template <typename T>
double max_epsilon(const RMatrix<T> &a, const RMatrix<T> &b) {
  double eps = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      const auto eps_temp = std::abs((a.Rmatrix(i, j) - b.Rmatrix(i, j)) /
                                     (a.Rmatrix(i, j) + b.Rmatrix(i, j)));
      if (eps_temp > eps)
        eps = eps_temp;
    }
  }
  return eps;
}

//==============================================================================

} // namespace MBPT
