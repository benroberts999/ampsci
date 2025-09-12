#pragma once
#include "LinAlg/Matrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <cassert>
#include <iostream>
#include <type_traits>

namespace MBPT {

/*! 
Radial matrix: stored on a sub-grid (defined by i0, stride, size).

XXX Should actually just derive/specialise SpinorMatrix!!

@details
 - i0 : the first grid-point included in subgrid
 - stride: stride between grid-points used in subgrid
 - size: total number of points used in subgrid
 - rgrid: pointer to full grid
 - T: type. Should usually be double or complex<double>
*/
template <typename T>
class RadialMatrix {

  std::size_t m_i0, m_stride;
  std::size_t m_size;
  LinAlg::Matrix<T> m_Rmatrix;
  std::shared_ptr<const Grid> m_rgrid; // "full" grid
  std::vector<double> sub_r{};         // sub grid (required for interpolation)

public:
  //============================================================================

  RadialMatrix(std::size_t i0, std::size_t stride, std::size_t size,
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

  T &at(std::size_t i, std::size_t j) { return m_Rmatrix(i, j); }
  T at(std::size_t i, std::size_t j) const { return m_Rmatrix(i, j); }

  T &operator()(std::size_t i, std::size_t j) { return at(i, j); }
  T operator()(std::size_t i, std::size_t j) const { return at(i, j); }

  //! direct access to radial matrix
  const LinAlg::Matrix<T> &Rmatrix() const { return m_Rmatrix; }
  LinAlg::Matrix<T> &Rmatrix() { return m_Rmatrix; }

  std::size_t size() const { return m_size; }
  std::size_t i0() const { return m_i0; }
  std::size_t stride() const { return m_stride; }

  //! First point along full grid for which subgrid is defined
  double r0() const { return sub_r.front(); }
  //! Last point along full grid for which subgrid is defined
  double rmax() const { return sub_r.back(); }

  //! returns r at position sub_i along sub grid
  double r(std::size_t sub_i) const { return sub_r.at(sub_i); }

  //! returns dr at position along sub grid
  double dr(std::size_t sub_i) const {
    const auto full_index = index_to_fullgrid(sub_i);
    return m_rgrid->drdu(sub_i) * m_rgrid->du() * double(m_stride);
  }

  //! Converts an index on the sub-grid to the full grid.
  std::size_t index_to_fullgrid(std::size_t i) const {
    return m_i0 + i * m_stride;
  }

  //============================================================================
  //! Sets the radial matrix to have all zero's
  void zero() { m_Rmatrix.zero(); }

  //============================================================================
  //! Matrix adition +,-
  RadialMatrix<T> &operator+=(const RadialMatrix<T> &rhs) {
    m_Rmatrix += rhs.m_Rmatrix;
    return *this;
  }
  //! Matrix adition +,-
  RadialMatrix<T> &operator-=(const RadialMatrix<T> &rhs) {
    m_Rmatrix -= rhs.m_Rmatrix;
    return *this;
  }
  //! Scalar multiplication
  RadialMatrix<T> &operator*=(const T x) {
    m_Rmatrix *= x;
    return *this;
  }

  //! Matrix adition +,-
  [[nodiscard]] friend RadialMatrix<T> operator+(RadialMatrix<T> lhs,
                                                 const RadialMatrix<T> &rhs) {
    return lhs += rhs;
  }
  //! Matrix adition +,-
  [[nodiscard]] friend RadialMatrix<T> operator-(RadialMatrix<T> lhs,
                                                 const RadialMatrix<T> &rhs) {
    return lhs -= rhs;
  }
  //! Scalar multiplication
  [[nodiscard]] friend RadialMatrix<T> operator*(const T x,
                                                 RadialMatrix<T> rhs) {
    return rhs *= x;
  }

  //! Adition of identity: Matrix<T> += T : T assumed to be *Identity!
  RadialMatrix<T> &operator+=(T aI) {
    m_Rmatrix += aI;
    return *this;
  }
  //! Adition of identity: Matrix<T> -= T : T assumed to be *Identity!
  RadialMatrix<T> &operator-=(T aI) {
    m_Rmatrix -= aI;
    return *this;
  }

  //! Adition of identity: Matrix<T> + T : T assumed to be *Identity!
  [[nodiscard]] friend RadialMatrix<T> operator+(RadialMatrix<T> M, T aI) {
    return (M += aI);
  }
  //! Adition of identity: Matrix<T> - T : T assumed to be *Identity!
  [[nodiscard]] friend RadialMatrix<T> operator-(RadialMatrix<T> M, T aI) {
    return (M -= aI);
  }

  //============================================================================
  //! Multiply coordinate elements (in place): Gij -> Gij*Bij
  RadialMatrix<T> &mult_elements_by(const RadialMatrix<T> &rhs) {
    m_Rmatrix.mult_elements_by(rhs.m_Rmatrix);
    return *this;
  }

  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[nodiscard]] friend RadialMatrix<T>
  mult_elements(RadialMatrix<T> lhs, const RadialMatrix<T> &rhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //! Matrix multiplication:  Gij = Aik*Bkj
  //! Note: integration measure not included: call .drj() first to include it!
  [[nodiscard]] friend RadialMatrix<T>
  matrix_multiply(const RadialMatrix<T> &a, const RadialMatrix<T> &b) {
    RadialMatrix<T> out(a.m_i0, a.m_stride, a.m_size, a.m_rgrid);
    out.Rmatrix() = a.Rmatrix() * b.Rmatrix();
    return out;
  }

  //! Matrix multiplication:  Gij = Aik*Bkj
  [[nodiscard]] friend RadialMatrix<T> operator*(const RadialMatrix<T> &a,
                                                 const RadialMatrix<T> &b) {
    return matrix_multiply(a, b);
  }

  //============================================================================

  //! Returns conjugate of matrix
  [[nodiscard]] RadialMatrix<T> conj() const {
    auto out = *this;
    out.Rmatrix().conj_in_place();
    return out;
  }

  //! Conjuagtes current matrix, in place
  RadialMatrix<T> &conj_in_place() {
    m_Rmatrix.conj_in_place();
    return *this;
  }

  //! Returns real part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] RadialMatrix<double> real() const {
    RadialMatrix<double> out(m_i0, m_stride, m_size, m_rgrid);
    out.Rmatrix() = m_Rmatrix.real();
    return out;
  }

  //! Returns imag part of complex matrix (changes type; returns a real matrix)
  [[nodiscard]] RadialMatrix<double> imag() const {
    RadialMatrix<double> out(m_i0, m_stride, m_size, m_rgrid);
    out.Rmatrix() = m_Rmatrix.imag();
    return out;
  }

  //! Converts a real to complex matrix (changes type; returns a complex matrix)
  [[nodiscard]] RadialMatrix<std::complex<double>> complex() const {
    RadialMatrix<std::complex<double>> out(m_i0, m_stride, m_size, m_rgrid);
    out.Rmatrix() = m_Rmatrix.complex();
    return out;
  }

  //============================================================================
  //! Inversion (in place)
  RadialMatrix<T> &invert_in_place() {
    m_Rmatrix.invert_in_place();
    return *this;
  }
  //! Returns inverse of matrix; original matrix unchanged
  [[nodiscard]] RadialMatrix<T> inverse() const {
    auto out = *this; //
    return out.invert_in_place();
  }

  //============================================================================
  //! Returns transpose of matrix; original matrix unchanged
  [[nodiscard]] RadialMatrix<T> transpose() const {
    RadialMatrix<T> out(m_i0, m_stride, m_size, m_rgrid);
    out.m_Rmatrix = m_Rmatrix.transpose();
    return out;
  }

  //============================================================================
  //! Multiplies by drj: Q_ij -> Q_ij*dr_j, in place
  RadialMatrix<T> &drj_in_place() {
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
  RadialMatrix<T> &dri_in_place() {
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
  [[nodiscard]] RadialMatrix<T> drj() const {
    auto out = *this;
    return out.drj_in_place();
  }

  //! Multiplies by dri: Q_ij -> Q_ij*dr_i. Returns new matrix (orig unchanged)
  [[nodiscard]] RadialMatrix<T> dri() const {
    auto out = *this;
    return out.dri_in_place();
  }
};

//==============================================================================

//! Checks if two matrix's are equal (to within parts in 10^12)
template <typename T>
bool equal(const RadialMatrix<T> &lhs, const RadialMatrix<T> &rhs) {
  return LinAlg::equal(lhs.Rmatrix(), rhs.Rmatrix());
}

//! returns maximum element (by abs)
template <typename T>
double max_element(const RadialMatrix<T> &a) {
  double xmax = 0.0;
  for (auto i = 0ul; i < a.size(); ++i) {
    for (auto j = 0ul; j < a.size(); ++j) {
      if (std::abs(a(i, j)) > xmax) {
        xmax = std::abs(a(i, j));
      }
    }
  }
  return xmax;
}

//! returns maximum difference (abs) between two matrixs
template <typename T>
double max_delta(const RadialMatrix<T> &a, const RadialMatrix<T> &b) {
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
double max_epsilon(const RadialMatrix<T> &a, const RadialMatrix<T> &b) {
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

using RMatrix = RadialMatrix<double>;

} // namespace MBPT
