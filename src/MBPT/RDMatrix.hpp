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

/*! Defines RDMatrix, Radial Dirac (Spinor) matrix. Designed to store
  Greens-function like operators: |Fa><Fb| (where Fa, Fb are radial Dirac
  spinors), as a radial matrix. The matrix is stored on a sub-grid (between r0
  and rmax), with a specified stride.

@details

RDMatrix is a 4*4 matrix in spinor space {ff, fg, gf, gg} - the g blocks are
small and are optional. Each block is an N*N radial matrix, where N is a subset
of the number of points along the full radial grid. May store doubles or complex
doubles.

   RDMatrix = {ff fg}
              {gf gg}
   RDMatrix * F = {ff fg} * (f)
                  {gf gg}   (g)
                = (ff(r,r')*f(r') + fg(r,r')*g(r'))
                  (gf(r,r')*f(r') + gg(r,r')*g(r'))

Note: Careful to distinguish RDMatrix multiplication/integration:
  G1 * G2 = Int G1(ra,rb)*G2(rb,rc)
          = Sum_j G1(i,j)*G2(j,k)

G1.drj() * G2 = Int G1(ra,rb)*G2(rb,rc)*dr_b
              = Sum_j G1(i,j)*G2(j,k)*drdu_j*du
              = G1 * G2.dri()

While almost always symmetric, this doesn't assume that.
*/
template <typename T> class RDMatrix {

  std::size_t m_i0, m_stride;
  std::size_t m_size;
  std::size_t m_g_size;
  LinAlg::Matrix<T> m_ff, m_fg, m_gf, m_gg;
  bool m_incl_g;
  std::shared_ptr<const Grid> m_rgrid; // "full" grid
  std::vector<double> sub_r{};         // sub grid

public:
  //****************************************************************************

  RDMatrix(std::size_t i0, std::size_t stride, std::size_t size, bool incl_g,
           std::shared_ptr<const Grid> rgrid)
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
    assert(m_i0 + m_stride * m_size <= r.size());
    for (std::size_t i = 0; i < m_size; ++i) {
      sub_r.push_back(r[index_to_fullgrid(i)]);
    }
    assert(m_i0 + m_stride * m_size <= r.size());
    assert(sub_r[1] == r[index_to_fullgrid(1)]);
    assert(sub_r[m_size - 1] == r[index_to_fullgrid(m_size - 1)]);
    //------------------
  }

  //****************************************************************************
  //! direct access to matrix elements
  T &ff(std::size_t i, std::size_t j) { return m_ff(i, j); }
  T &fg(std::size_t i, std::size_t j) { return m_fg(i, j); }
  T &gf(std::size_t i, std::size_t j) { return m_gf(i, j); }
  T &gg(std::size_t i, std::size_t j) { return m_gg(i, j); }
  const T ff(std::size_t i, std::size_t j) const { return m_ff(i, j); }
  const T fg(std::size_t i, std::size_t j) const { return m_fg(i, j); }
  const T gf(std::size_t i, std::size_t j) const { return m_gf(i, j); }
  const T gg(std::size_t i, std::size_t j) const { return m_gg(i, j); }

  //! direct access to matrix's
  const LinAlg::Matrix<T> &ff() const { return m_ff; }
  const LinAlg::Matrix<T> &fg() const { return m_fg; }
  const LinAlg::Matrix<T> &gf() const { return m_gf; }
  const LinAlg::Matrix<T> &gg() const { return m_gg; }
  LinAlg::Matrix<T> &set_ff() { return m_ff; }
  LinAlg::Matrix<T> &set_fg() { return m_fg; }
  LinAlg::Matrix<T> &set_gf() { return m_gf; }
  LinAlg::Matrix<T> &set_gg() { return m_gg; }

  std::size_t size() const { return m_size; }
  std::size_t g_size() const { return m_g_size; }

  //****************************************************************************
  //! Sets all matrix elements to zero
  void zero() {
    m_ff.zero();
    m_fg.zero();
    m_gf.zero();
    m_gg.zero();
  }

  // Makes 1
  [[deprecated]] void make_identity() {
    m_ff.make_identity();
    m_fg.make_identity();
    m_gf.make_identity();
    m_gg.make_identity();
  }

      // only for transition, kill!
      [[deprecated]] RDMatrix<T> &plusIdent(T a = T{1.0}) {
    (*this) += a;
    return *this;
  }

  //****************************************************************************
  //! Matrix adition +,-
  RDMatrix<T> &operator+=(const RDMatrix<T> &rhs) {
    m_ff += rhs.m_ff;
    m_fg += rhs.m_fg;
    m_gf += rhs.m_gf;
    m_gg += rhs.m_gg;
    return *this;
  }
  //! Matrix adition +,-
  RDMatrix<T> &operator-=(const RDMatrix<T> &rhs) {
    m_ff -= rhs.m_ff;
    m_fg -= rhs.m_fg;
    m_gf -= rhs.m_gf;
    m_gg -= rhs.m_gg;
    return *this;
  }
  //! Scalar multiplication
  RDMatrix<T> &operator*=(const T x) {
    m_ff *= x;
    m_fg *= x;
    m_gf *= x;
    m_gg *= x;
    return *this;
  }

  //! Matrix adition +,-
  [[nodiscard]] friend RDMatrix<T> operator+(RDMatrix<T> lhs,
                                             const RDMatrix<T> &rhs) {
    return (lhs += rhs);
  }
  //! Matrix adition +,-
  [[nodiscard]] friend RDMatrix<T> operator-(RDMatrix<T> lhs,
                                             const RDMatrix<T> &rhs) {
    return (lhs -= rhs);
  }
  //! Scalar multiplication
  [[nodiscard]] friend RDMatrix<T> operator*(const T x, RDMatrix<T> rhs) {
    return (rhs *= x);
  }

  //! Adition of identity: Matrix<T> += T : T assumed to be *Identity!
  RDMatrix<T> &operator+=(T aI) {
    m_ff += aI;
    m_gg += aI;
    return *this;
  }
  //! Adition of identity: Matrix<T> -= T : T assumed to be *Identity!
  RDMatrix<T> &operator-=(T aI) {
    m_ff -= aI;
    m_gg -= aI;
    return *this;
  }

  //! Adition of identity: Matrix<T> + T : T assumed to be *Identity!
  [[nodiscard]] friend RDMatrix<T> operator+(RDMatrix<T> M, T aI) {
    return (M += aI);
  }
  //! Adition of identity: Matrix<T> - T : T assumed to be *Identity!
  [[nodiscard]] friend RDMatrix<T> operator-(RDMatrix<T> M, T aI) {
    return (M -= aI);
  }

  //****************************************************************************

  //! Matrix multplication: C=A*B := Cij = \sum_k Aik*Bkj.
  //! Note: integration measure not included: call .drj() first to include it!
  [[nodiscard]] friend RDMatrix<T> operator*(const RDMatrix<T> &a,
                                             const RDMatrix<T> &b) {

    RDMatrix<T> out(a.m_i0, a.m_stride, a.m_size, a.m_incl_g, a.m_rgrid);

    // FF = FF*FF + FG*GF
    // FG = FF*FG + FG*GG
    // GF = GF*FF + GG*GF
    // GG = GF*FG + GG*GG
    out.set_ff() = a.ff() * b.ff();
    if (a.m_incl_g) {
      out.set_ff() += a.fg() * b.gf();
      out.set_fg() = a.ff() * b.fg() + a.fg() * b.gg();
      out.set_gf() = a.gf() * b.ff() + a.gg() * b.gf();
      out.set_gg() = a.gf() * b.fg() + a.gg() * b.gg();
    }
    return out;
  }

  //****************************************************************************
  //! Multiply elements (in place): Gij -> Gij*Bij
  RDMatrix<T> &mult_elements_by(const RDMatrix<T> &rhs) {
    m_ff.mult_elements_by(rhs.ff());
    m_fg.mult_elements_by(rhs.fg());
    m_gf.mult_elements_by(rhs.gf());
    m_gg.mult_elements_by(rhs.gg());
    return *this;
  }
  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[nodiscard]] friend RDMatrix<T> mult_elements(RDMatrix<T> lhs,
                                                 const RDMatrix<T> &rhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //****************************************************************************

  //! Returns conjugate of matrix
  [[nodiscard]] RDMatrix<T> conj() const {
    auto out = *this;
    out.set_ff().conj_in_place();
    out.set_fg().conj_in_place();
    out.set_gf().conj_in_place();
    out.set_gg().conj_in_place();
    return out;
  }
      //! Returns real part of complex matrix (changes type; returns a real
      //! matrix)
          [[nodiscard]] RDMatrix<double> real() const {
    RDMatrix<double> out(m_i0, m_stride, m_size, m_incl_g, m_rgrid);
    out.set_ff() = m_ff.real();
    out.set_fg() = m_fg.real();
    out.set_gf() = m_gf.real();
    out.set_gg() = m_gg.real();
    return out;
  }
  //! Returns imag part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] RDMatrix<double> imag() const {
    RDMatrix<double> out(m_i0, m_stride, m_size, m_incl_g, m_rgrid);
    out.set_ff() = m_ff.imag();
    out.set_fg() = m_fg.imag();
    out.set_gf() = m_gf.imag();
    out.set_gg() = m_gg.imag();
    return out;
  }
      //! Converts a real to complex matrix (changes type; returns a complex
      //! matrix)
          [[nodiscard]] RDMatrix<std::complex<double>> complex() const {
    RDMatrix<std::complex<double>> out(m_i0, m_stride, m_size, m_incl_g,
                                       m_rgrid);
    out.set_ff() = m_ff.complex();
    out.set_fg() = m_fg.complex();
    out.set_gf() = m_gf.complex();
    out.set_gg() = m_gg.complex();
    return out;
  }

  //****************************************************************************
  //! Inversion (in place)
  RDMatrix<T> &invert_in_place() {
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
  [[nodiscard]] RDMatrix<T> inverse() const {
    auto out = *this; //
    return out.invert_in_place();
  }

  //****************************************************************************
  //! Multiplies by drj: Q_ij -> Q_ij*dr_j, in place
  RDMatrix<T> &drj_in_place() {
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
  RDMatrix<T> &dri_in_place() {
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
  RDMatrix<T> drj() const {
    auto out = *this;
    return out.drj_in_place();
  }
  //! Multiplies by dri: Q_ij -> Q_ij*dr_i. Returns new matrix (orig unchanged)
  RDMatrix<T> dri() const {
    auto out = *this;
    return out.dri_in_place();
  }

  //! returns dr at position along sub grid
  double dr(std::size_t sub_index) const {
    const auto full_index = index_to_fullgrid(sub_index);
    return m_rgrid->drdu(full_index) * m_rgrid->du() * double(m_stride);
  }

  //****************************************************************************
  //! Converts an index on the sub-grid to the full grid.
  std::size_t index_to_fullgrid(std::size_t i) const {
    return m_i0 + i * m_stride;
  }

  //****************************************************************************
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
    }   // i

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
      }   // i
    }
  }

  //****************************************************************************
  //! Action of RDMatrix operator on DiracSpinor. Inludes Integration:
  //! G*F = Int[G(r,r')*F(r') dr'] = sum_j G_ij*F_j*drdu_j*du
  DiracSpinor operator*(const DiracSpinor &Fn) const {

    const auto &r = Fn.rgrid->r();
    // const auto &drdu = Fn.rgrid->drdu();
    // const double s_du = double(m_stride) * Fn.rgrid->du();

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
          g[i] += (m_gf(i, j) * Fn.f(j_f) +
                   m_gg(i, j) * Fn.g(j_f)); // *
                                            // drdu[j_f] * s_du;
        }
      }
    }

    DiracSpinor out = Fn;
    // Interpolate from sub-grid to full grid
    out.set_f() = Interpolator::interpolate(sub_r, f, r);
    if (m_incl_g) {
      out.set_g() = Interpolator::interpolate(sub_r, g, r);
    }
    return out;
  }

  //****************************************************************************
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
};

//! Checks if two matrix's are equal (to within parts in 10^12)
template <typename T>
bool equal(const RDMatrix<T> &lhs, const RDMatrix<T> &rhs) {
  return equal(lhs.ff(), rhs.ff()) && equal(lhs.fg(), rhs.fg()) &&
         equal(lhs.gf(), rhs.gf()) && equal(lhs.gg(), rhs.gg());
}

//! returns maximum element (by abs)
template <typename T> double max_element(const RDMatrix<T> &a) {
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

//! returns maximum difference (abs) between two matrixs
template <typename T>
double max_delta(const RDMatrix<T> &a, const RDMatrix<T> &b) {
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

//! returns maximum relative diference [aij-bij/(aij+bij)] (abs) between two
//! matrices
template <typename T>
double max_epsilon(const RDMatrix<T> &a, const RDMatrix<T> &b) {
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
        if (eps_ff > xfg)
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

//******************************************************************************
using GMatrix = RDMatrix<double>;
using ComplexGMatrix = RDMatrix<std::complex<double>>;
using ComplexDouble = std::complex<double>;

} // namespace MBPT
