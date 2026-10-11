#pragma once
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "RadialMatrix.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <type_traits>

namespace MBPT {

/*!
  @brief Radial Dirac spinor matrix.

  @details
  Designed to store Greens-function like operators: \f$ \ket{F_a}\bra{F_b} \f$
  (where Fa, Fb are radial Dirac spinors), as a radial matrix. The matrix is
  stored on a sub-grid (between r0 and rmax), with a specified stride.

  SpinorMatrix is a 2*2 matrix in spinor space {ff, fg, gf, gg} - the g blocks
  are small and are optional. Each block is an N*N radial matrix, where N is a
  subset of the number of points along the full radial grid. May store doubles
  or complex doubles.

  \f[
  \begin{align}
    G &= \begin{pmatrix} G_{ff} & G_{fg} \\ G_{gf} & G_{gg} \end{pmatrix},\\
    G F &= \begin{pmatrix} G_{ff} & G_{fg} \\ G_{gf} & G_{gg} \end{pmatrix}
          \begin{pmatrix} f \\ g \end{pmatrix}
        = \begin{pmatrix}
            G_{ff}(r,r')\,f(r') + G_{fg}(r,r')\,g(r') \\
            G_{gf}(r,r')\,f(r') + G_{gg}(r,r')\,g(r')
          \end{pmatrix}.
  \end{align}
  \f]

  Careful to distinguish SpinorMatrix multiplication/integration:
  \f[
  \begin{align}
    G_1 * G_2 &= \int G_1(r_a,r_b)\,G_2(r_b,r_c)
               = \sum_j G_1(i,j)\,G_2(j,k),\\
    G_1.{\rm drj}() * G_2 &= \int G_1(r_a,r_b)\,G_2(r_b,r_c)\,dr_b
                           = \sum_j G_1(i,j)\,G_2(j,k)\,dr_j
                           = G_1 * G_2.{\rm dri}(),
  \end{align}
  \f]
  where \f$ dr_j = (dr/du)_j \times du \times {\rm stride} \f$ (see dr()).

  While almost always symmetric, this doesn't assume that.
*/
template <typename T>
class SpinorMatrix {

  std::size_t m_i0, m_stride;
  std::size_t m_size;
  std::size_t m_g_size;
  LinAlg::Matrix<T> m_ff, m_fg, m_gf, m_gg;
  bool m_incl_g;
  std::shared_ptr<const Grid> m_rgrid; // "full" grid
  std::vector<double> sub_r{};         // sub grid

public:
  //============================================================================

  /*!
    @brief Constructs a zero matrix on the given radial sub-grid.

    @details
    Sub-grid point i is full-grid point i0 + i*stride (see index_to_fullgrid).
    The last sub-grid point must lie on the full grid (checked by assert).

    @param i0      Full-grid index of the first sub-grid point
    @param stride  Sub-grid stride: every stride-th point of the full grid
    @param size    Number of sub-grid points (each block is size*size)
    @param incl_g  Include the g blocks (fg, gf, gg)
    @param rgrid   Full radial grid
  */
  SpinorMatrix(std::size_t i0, std::size_t stride, std::size_t size,
               bool incl_g, std::shared_ptr<const Grid> rgrid)
    : m_i0(i0),
      m_stride(stride),
      m_size(size),
      m_g_size(incl_g ? size : 0),
      m_ff(m_size),
      m_fg(m_g_size),
      m_gf(m_g_size),
      m_gg(m_g_size),
      m_incl_g(incl_g),
      m_rgrid(rgrid) {
    //------------------
    // create vector of r on sub-grid, used to interpolate values onto full
    const auto &r = m_rgrid->r();
    sub_r.reserve(m_size);
    // last sub-grid point, i0 + stride*(size-1), must lie on the full grid
    assert(m_size > 0 && index_to_fullgrid(m_size - 1) < r.size());
    for (std::size_t i = 0; i < m_size; ++i) {
      sub_r.push_back(r[index_to_fullgrid(i)]);
    }
    assert(m_size < 2 || sub_r[1] == r[index_to_fullgrid(1)]);
    assert(sub_r[m_size - 1] == r[index_to_fullgrid(m_size - 1)]);
    //------------------
  }

  //============================================================================
  //! direct access to matrix elements
  T &ff(std::size_t i, std::size_t j) { return m_ff(i, j); }
  T &fg(std::size_t i, std::size_t j) { return m_fg(i, j); }
  T &gf(std::size_t i, std::size_t j) { return m_gf(i, j); }
  T &gg(std::size_t i, std::size_t j) { return m_gg(i, j); }
  const T ff(std::size_t i, std::size_t j) const { return m_ff(i, j); }
  const T fg(std::size_t i, std::size_t j) const { return m_fg(i, j); }
  const T gf(std::size_t i, std::size_t j) const { return m_gf(i, j); }
  const T gg(std::size_t i, std::size_t j) const { return m_gg(i, j); }

  //! direct access to matrices
  const LinAlg::Matrix<T> &ff() const { return m_ff; }
  const LinAlg::Matrix<T> &fg() const { return m_fg; }
  const LinAlg::Matrix<T> &gf() const { return m_gf; }
  const LinAlg::Matrix<T> &gg() const { return m_gg; }
  LinAlg::Matrix<T> &ff() { return m_ff; }
  LinAlg::Matrix<T> &fg() { return m_fg; }
  LinAlg::Matrix<T> &gf() { return m_gf; }
  LinAlg::Matrix<T> &gg() { return m_gg; }

  /*!
    @brief Spinor block by index: 0 = f, 1 = g.

    @details
    (mu, nu) = (0, 0): ff; (0, 1): fg; (1, 0): gf; (1, 1): gg.

    @note The g blocks are empty (0*0) if g parts are not included.
    @warning @p mu and @p nu must be 0 or 1; aborts otherwise.
  */
  LinAlg::Matrix<T> &sp(std::size_t mu, std::size_t nu) {
    assert(mu < 2 && nu < 2);
    if (mu == 0 && nu == 0)
      return m_ff;
    if (mu == 0 && nu == 1)
      return m_fg;
    if (mu == 1 && nu == 0)
      return m_gf;
    if (mu == 1 && nu == 1)
      return m_gg;
    assert(false && "SpinorMatrix::sp: mu, nu must be 0 or 1");
    std::abort();
  }

  //! Spinor block by index: 0 = f, 1 = g (const)
  const LinAlg::Matrix<T> &sp(std::size_t mu, std::size_t nu) const {
    assert(mu < 2 && nu < 2);
    if (mu == 0 && nu == 0)
      return m_ff;
    if (mu == 0 && nu == 1)
      return m_fg;
    if (mu == 1 && nu == 0)
      return m_gf;
    if (mu == 1 && nu == 1)
      return m_gg;
    assert(false && "SpinorMatrix::sp: mu, nu must be 0 or 1");
    std::abort();
  }

  //! Number of sub-grid points (dimension of each block)
  std::size_t size() const { return m_size; }
  //! Dimension of the g blocks: size(), or 0 if g parts are not included
  std::size_t g_size() const { return m_g_size; }
  //! True if g parts (fg, gf, gg) are included
  bool includes_g() const { return m_g_size == m_size; };
  //! Full-grid index of the first sub-grid point
  std::size_t i0() const { return m_i0; }
  //! Sub-grid stride: every stride-th point of the full grid
  std::size_t stride() const { return m_stride; }

  //============================================================================
  //! Sets all matrix elements to zero
  void zero() {
    m_ff.zero();
    m_fg.zero();
    m_gf.zero();
    m_gg.zero();
  }

  //============================================================================
  //! Kills g parts of spinor matrix, in place!
  SpinorMatrix<T> &drop_g() {
    m_g_size = 0;
    m_incl_g = false;
    m_fg.resize(0, 0);
    m_gf.resize(0, 0);
    m_gg.resize(0, 0);
    return *this;
  }

  /*!
    @brief Creates g parts of spinor matrix - will have value 0.

    @note If g parts already exist, they are reset to 0.
  */
  SpinorMatrix<T> &create_g() {
    m_g_size = m_size;
    m_incl_g = true;
    m_fg.resize(m_size, m_size);
    m_gf.resize(m_size, m_size);
    m_gg.resize(m_size, m_size);
    return *this;
  }

  //============================================================================
  /*!
    @brief Matrix addition +,-.

    @details
    A matrix without g parts is treated as having zero g parts: adding one to
    a matrix with g parts leaves its g parts unchanged; adding a matrix with g
    parts to one without creates them.
  */
  SpinorMatrix<T> &operator+=(const SpinorMatrix<T> &rhs) {
    m_ff += rhs.m_ff;
    if (rhs.m_incl_g) {
      if (!m_incl_g) {
        create_g();
      }
      m_fg += rhs.m_fg;
      m_gf += rhs.m_gf;
      m_gg += rhs.m_gg;
    }
    return *this;
  }
  //! Matrix addition +,- (see operator+= for matrices without g parts)
  SpinorMatrix<T> &operator-=(const SpinorMatrix<T> &rhs) {
    m_ff -= rhs.m_ff;
    if (rhs.m_incl_g) {
      if (!m_incl_g) {
        create_g();
      }
      m_fg -= rhs.m_fg;
      m_gf -= rhs.m_gf;
      m_gg -= rhs.m_gg;
    }
    return *this;
  }
  //! Scalar multiplication
  SpinorMatrix<T> &operator*=(const T x) {
    m_ff *= x;
    m_fg *= x;
    m_gf *= x;
    m_gg *= x;
    return *this;
  }

  //! Matrix addition +,-
  [[nodiscard]] friend SpinorMatrix<T> operator+(SpinorMatrix<T> lhs,
                                                 const SpinorMatrix<T> &rhs) {
    return (lhs += rhs);
  }
  //! Matrix addition +,-
  [[nodiscard]] friend SpinorMatrix<T> operator-(SpinorMatrix<T> lhs,
                                                 const SpinorMatrix<T> &rhs) {
    return (lhs -= rhs);
  }
  //! Scalar multiplication
  [[nodiscard]] friend SpinorMatrix<T> operator*(const T x,
                                                 SpinorMatrix<T> rhs) {
    return (rhs *= x);
  }

  //! Addition of identity: Matrix<T> += T : T assumed to be *Identity!
  SpinorMatrix<T> &operator+=(T aI) {
    m_ff += aI;
    m_gg += aI;
    return *this;
  }
  //! Addition of identity: Matrix<T> -= T : T assumed to be *Identity!
  SpinorMatrix<T> &operator-=(T aI) {
    m_ff -= aI;
    m_gg -= aI;
    return *this;
  }

  //! Addition of identity: Matrix<T> + T : T assumed to be *Identity!
  [[nodiscard]] friend SpinorMatrix<T> operator+(SpinorMatrix<T> M, T aI) {
    return (M += aI);
  }
  //! Addition of identity: Matrix<T> - T : T assumed to be *Identity!
  [[nodiscard]] friend SpinorMatrix<T> operator-(SpinorMatrix<T> M, T aI) {
    return (M -= aI);
  }

  //============================================================================

  /*!
    @brief Matrix multiplication.

    @details
    \f$ C=A\times B \equiv C_{ij} = \sum_k A_{ik}\,B_{kj}, \f$ including the
    sum over spinor blocks (e.g., ff = a.ff*b.ff + a.fg*b.gf).

    A matrix without g parts is treated as having zero g parts (as in
    operator+=); the result includes g parts if either @p a or @p b does.

    @note Integration measure not automatically included: call .drj() first
          to include it!
  */
  [[nodiscard]] friend SpinorMatrix<T> operator*(const SpinorMatrix<T> &a,
                                                 const SpinorMatrix<T> &b) {

    SpinorMatrix<T> out(a.m_i0, a.m_stride, a.m_size, a.m_incl_g || b.m_incl_g,
                        a.m_rgrid);

    // FF = FF*FF + FG*GF
    // FG = FF*FG + FG*GG
    // GF = GF*FF + GG*GF
    // GG = GF*FG + GG*GG
    out.ff() = a.ff() * b.ff();
    if (a.m_incl_g && b.m_incl_g) {
      out.ff() += a.fg() * b.gf();
      out.fg() = a.ff() * b.fg() + a.fg() * b.gg();
      out.gf() = a.gf() * b.ff() + a.gg() * b.gf();
      out.gg() = a.gf() * b.fg() + a.gg() * b.gg();
    } else if (a.m_incl_g) {
      // b has zero g parts: only GF = GF*FF survives
      out.gf() = a.gf() * b.ff();
    } else if (b.m_incl_g) {
      // a has zero g parts: only FG = FF*FG survives
      out.fg() = a.ff() * b.fg();
    }
    return out;
  }

  //============================================================================
  //! Multiply elements (in place): Gij -> Gij*Bij
  [[deprecated]] SpinorMatrix<T> &mult_elements_by(const SpinorMatrix<T> &rhs) {
    m_ff.mult_elements_by(rhs.ff());
    if (this->m_incl_g) {
      if (rhs.m_incl_g) {
        m_fg.mult_elements_by(rhs.fg());
        m_gf.mult_elements_by(rhs.gf());
        m_gg.mult_elements_by(rhs.gg());
      } else {
        // rhs has zero g parts
        m_fg.zero();
        m_gf.zero();
        m_gg.zero();
      }
    }
    return *this;
  }
  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[deprecated]] [[nodiscard]] friend SpinorMatrix<T>
  mult_elements(SpinorMatrix<T> lhs, const SpinorMatrix<T> &rhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //============================================================================
  //! Multiply elements (in place): Gij -> Gij*Bij
  SpinorMatrix<T> &mult_elements_by(const RadialMatrix<T> &rhs) {
    m_ff.mult_elements_by(rhs.Rmatrix());
    if (this->m_incl_g) {
      m_fg.mult_elements_by(rhs.Rmatrix());
      m_gf.mult_elements_by(rhs.Rmatrix());
      m_gg.mult_elements_by(rhs.Rmatrix());
    }
    return *this;
  }

  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[nodiscard]] friend SpinorMatrix<T>
  mult_elements(SpinorMatrix<T> lhs, const RadialMatrix<T> &rhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[nodiscard]] friend SpinorMatrix<T> mult_elements(const RadialMatrix<T> &rhs,
                                                     SpinorMatrix<T> lhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //============================================================================

  //! Returns conjugate of matrix
  [[nodiscard]] SpinorMatrix<T> conj() const {
    auto out = *this;
    out.ff().conj_in_place();
    out.fg().conj_in_place();
    out.gf().conj_in_place();
    out.gg().conj_in_place();
    return out;
  }
  //! Returns real part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] SpinorMatrix<double> real() const {
    SpinorMatrix<double> out(m_i0, m_stride, m_size, m_incl_g, m_rgrid);
    out.ff() = m_ff.real();
    out.fg() = m_fg.real();
    out.gf() = m_gf.real();
    out.gg() = m_gg.real();
    return out;
  }
  //! Returns imag part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] SpinorMatrix<double> imag() const {
    SpinorMatrix<double> out(m_i0, m_stride, m_size, m_incl_g, m_rgrid);
    out.ff() = m_ff.imag();
    out.fg() = m_fg.imag();
    out.gf() = m_gf.imag();
    out.gg() = m_gg.imag();
    return out;
  }
  //! Converts a real to complex matrix (changes type; returns a complex matrix)
  [[nodiscard]] SpinorMatrix<std::complex<double>> complex() const {
    SpinorMatrix<std::complex<double>> out(m_i0, m_stride, m_size, m_incl_g,
                                           m_rgrid);
    out.ff() = m_ff.complex();
    out.fg() = m_fg.complex();
    out.gf() = m_gf.complex();
    out.gg() = m_gg.complex();
    return out;
  }

  //============================================================================
  //! Inversion (in place)
  SpinorMatrix<T> &invert_in_place() {
    m_ff.invert_in_place();
    if (m_incl_g) {
      const auto &ai = m_ff; // already inverted
      const auto &b = m_fg;
      const auto &c = m_gf;
      const auto &d = m_gg;
      const auto cai = c * ai;
      const auto dmcaib = (d - cai * b).invert_in_place();
      const auto aib_dmcaib = ai * b * dmcaib;
      m_ff += aib_dmcaib * cai;
      m_fg = -1.0 * aib_dmcaib;
      m_gf = -1.0 * dmcaib * cai;
      m_gg = dmcaib;
    }
    return *this;
  }
  //! Returns inverse of matrix; original matrix unchanged
  [[nodiscard]] SpinorMatrix<T> inverse() const {
    auto out = *this; //
    return out.invert_in_place();
  }

  //============================================================================
  //! Multiplies by drj: Q_ij -> Q_ij*dr_j, in place
  SpinorMatrix<T> &drj_in_place() {
    const auto dus = m_rgrid->du() * double(m_stride);
    for (auto i = 0ul; i < m_size; ++i) {
      for (auto j = 0ul; j < m_size; ++j) {
        const auto sj = index_to_fullgrid(j);
        const auto dr = m_rgrid->drdu(sj) * dus;
        m_ff[i][j] *= dr;
      }
    }
    if (m_incl_g) {
      for (auto i = 0ul; i < m_size; ++i) {
        for (auto j = 0ul; j < m_size; ++j) {
          const auto sj = index_to_fullgrid(j);
          const auto dr = m_rgrid->drdu(sj) * dus;
          m_fg[i][j] *= dr;
          m_gf[i][j] *= dr;
          m_gg[i][j] *= dr;
        }
      }
    }
    return *this;
  }
  //! Multiplies by dri: Q_ij -> Q_ij*dr_i, in place
  SpinorMatrix<T> &dri_in_place() {
    const auto dus = m_rgrid->du() * double(m_stride);
    for (auto i = 0ul; i < m_size; ++i) {
      const auto si = index_to_fullgrid(i);
      const auto dr = m_rgrid->drdu(si) * dus;
      for (auto j = 0ul; j < m_size; ++j) {
        m_ff[i][j] *= dr;
      }
    }
    if (m_incl_g) {
      for (auto i = 0ul; i < m_size; ++i) {
        const auto si = index_to_fullgrid(i);
        const auto dr = m_rgrid->drdu(si) * dus;
        for (auto j = 0ul; j < m_size; ++j) {
          m_fg[i][j] *= dr;
          m_gf[i][j] *= dr;
          m_gg[i][j] *= dr;
        }
      }
    }
    return *this;
  }
  //! Multiplies by drj: Q_ij -> Q_ij*dr_j. Returns new matrix (orig unchanged)
  SpinorMatrix<T> drj() const {
    auto out = *this;
    return out.drj_in_place();
  }
  //! Multiplies by dri: Q_ij -> Q_ij*dr_i. Returns new matrix (orig unchanged)
  SpinorMatrix<T> dri() const {
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
      }
    }

    if (m_incl_g) {
      for (auto i = 0ul; i < m_size; ++i) {
        const auto si = index_to_fullgrid(i);
        for (auto j = 0ul; j < m_size; ++j) {
          const auto sj = index_to_fullgrid(j);
          // XXX Double check fg/gf right way!
          m_fg[i][j] += k * ket.f(si) * bra.g(sj);
          m_gf[i][j] += k * ket.g(si) * bra.f(sj); // symmetric, transpose?
          m_gg[i][j] += k * ket.g(si) * bra.g(sj);
        }
      }
    }
  }

  //============================================================================
  /*!
    @brief Action of SpinorMatrix operator on DiracSpinor.

    @details
    Assumes matrix already includes integration measure. Formed on the
    sub-grid, then interpolated onto the full grid of @p Fn. Result is zero
    outside the sub-grid.

    @note If g parts are not included, the g component of the result is zero
          (and the g component of @p Fn is not used).
  */
  DiracSpinor operator*(const DiracSpinor &Fn) const {

    const auto &r = Fn.grid().r();

    std::vector<double> f(m_size), g;
    for (auto i = 0ul; i < m_size; ++i) {
      for (auto j = 0ul; j < m_size; ++j) {
        const auto j_f = index_to_fullgrid(j);
        f[i] += m_ff(i, j) * Fn.f(j_f);
      }
    }
    if (m_incl_g) {
      g.resize(m_size);
      for (auto i = 0ul; i < m_size; ++i) {
        for (auto j = 0ul; j < m_size; ++j) {
          const auto j_f = index_to_fullgrid(j);
          f[i] += m_fg(i, j) * Fn.g(j_f);
          g[i] += (m_gf(i, j) * Fn.f(j_f) + m_gg(i, j) * Fn.g(j_f));
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
  friend std::ostream &operator<<(std::ostream &os, const SpinorMatrix<T> &a) {
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
};

//! Checks if two matrices are equal (to within parts in 10^12)
template <typename T>
bool equal(const SpinorMatrix<T> &lhs, const SpinorMatrix<T> &rhs) {
  return equal(lhs.ff(), rhs.ff()) && equal(lhs.fg(), rhs.fg()) &&
         equal(lhs.gf(), rhs.gf()) && equal(lhs.gg(), rhs.gg());
}

//! returns maximum element (by abs)
template <typename T>
double max_element(const SpinorMatrix<T> &a) {
  double xff = 0.0, xfg = 0.0, xgf = 0.0, xgg = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      if (std::abs(a.ff(i, j)) > xff)
        xff = std::abs(a.ff(i, j));
      if (a.g_size() != 0) {
        if (std::abs(a.fg(i, j)) > xfg)
          xfg = std::abs(a.fg(i, j));
        if (std::abs(a.gf(i, j)) > xgf)
          xgf = std::abs(a.gf(i, j));
        if (std::abs(a.gg(i, j)) > xgg)
          xgg = std::abs(a.gg(i, j));
      }
    }
  }
  return std::max({xff, xfg, xgf, xgg});
}

/*!
  @brief Returns maximum difference (abs) between two matrices.

  @details
  A matrix without g parts is treated as having zero g parts.
*/
template <typename T>
double max_delta(const SpinorMatrix<T> &a, const SpinorMatrix<T> &b) {
  if (a.includes_g() != b.includes_g()) {
    return a.includes_g() ? max_delta(a, SpinorMatrix<T>(b).create_g()) :
                            max_delta(SpinorMatrix<T>(a).create_g(), b);
  }
  double xff = 0.0, xfg = 0.0, xgf = 0.0, xgg = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      if (std::abs(a.ff(i, j) - b.ff(i, j)) > xff)
        xff = std::abs(a.ff(i, j) - b.ff(i, j));
      if (a.g_size() != 0) {
        if (std::abs(a.fg(i, j) - b.fg(i, j)) > xfg)
          xfg = std::abs(a.fg(i, j) - b.fg(i, j));
        if (std::abs(a.gf(i, j) - b.gf(i, j)) > xgf)
          xgf = std::abs(a.gf(i, j) - b.gf(i, j));
        if (std::abs(a.gg(i, j) - b.gg(i, j)) > xgg)
          xgg = std::abs(a.gg(i, j) - b.gg(i, j));
      }
    }
  }
  return std::max({xff, xfg, xgf, xgg});
}

/*!
  @brief Returns maximum relative difference (abs) between two matrices.

  @details
  Maximum over all elements of \f$ |(a_{ij}-b_{ij})/(a_{ij}+b_{ij})| \f$.
  A matrix without g parts is treated as having zero g parts.
*/
template <typename T>
double max_epsilon(const SpinorMatrix<T> &a, const SpinorMatrix<T> &b) {
  if (a.includes_g() != b.includes_g()) {
    return a.includes_g() ? max_epsilon(a, SpinorMatrix<T>(b).create_g()) :
                            max_epsilon(SpinorMatrix<T>(a).create_g(), b);
  }
  double xff = 0.0, xfg = 0.0, xgf = 0.0, xgg = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      const auto eps_ff =
        std::abs((a.ff(i, j) - b.ff(i, j)) / (a.ff(i, j) + b.ff(i, j)));
      if (eps_ff > xff)
        xff = eps_ff;
      if (a.g_size() != 0) {
        const auto eps_fg =
          std::abs((a.fg(i, j) - b.fg(i, j)) / (a.fg(i, j) + b.fg(i, j)));
        const auto eps_gf =
          std::abs((a.gf(i, j) - b.gf(i, j)) / (a.gf(i, j) + b.gf(i, j)));
        const auto eps_gg =
          std::abs((a.gg(i, j) - b.gg(i, j)) / (a.gg(i, j) + b.gg(i, j)));
        if (eps_fg > xfg)
          xfg = eps_fg;
        if (eps_gf > xgf)
          xgf = eps_gf;
        if (eps_gg > xgg)
          xgg = eps_gg;
      }
    }
  }
  return std::max({xff, xfg, xgf, xgg});
}

//==============================================================================
using GMatrix = SpinorMatrix<double>;
using ComplexGMatrix = SpinorMatrix<std::complex<double>>;
using ComplexDouble = std::complex<double>;

} // namespace MBPT
