#pragma once
#include "LinAlg/LinAlg.hpp"
#include <iostream>
#include <type_traits>
namespace MBPT {

//******************************************************************************
//! Holds Green's fn operator of form: |ket><bra| [4x4 matrix of NxN matrix]
template <typename T> class GreenMatrix {
public:
  bool m_include_G;
  std::size_t size;
  std::size_t G_size;

public:
  T ff, fg, gf, gg;

public:
  GreenMatrix(std::size_t in_size, bool in_include_G)
      : m_include_G(in_include_G),
        size(in_size),
        G_size(m_include_G ? size : 0),
        ff(size),
        fg(G_size),
        gf(G_size),
        gg(G_size) {
    zero();
  }

  //! Sets all matrix elements to zero
  void zero() {
    ff.zero();
    if (m_include_G) {
      fg.zero();
      gf.zero();
      gg.zero();
    }
  }

  //! Makes 1
  void make_identity() {
    ff.make_identity();
    if (m_include_G) {
      fg.make_identity();
      gf.make_identity();
      gg.make_identity();
    }
  }

  GreenMatrix<T> &plusIdent(double a = 1.0) {
    ff.plusIdent(a);
    if (m_include_G) {
      gg.plusIdent(a);
    }
    return *this;
  }

  GreenMatrix<T> &plusIdent(double re, double im) {
    static_assert(std::is_same<T, LinAlg::Matrix<std::complex<double>>>::value,
                  "Can only call plusIdent(x,y) from Complex GMatrix!");
    ff.plusIdent(re, im);
    if (m_include_G) {
      gg.plusIdent(re, im);
    }
    return *this;
  }

  void checkNaN() const {
    ff.checkNaN();
    if (m_include_G) {
      fg.checkNaN();
      gf.checkNaN();
      gg.checkNaN();
    }
  }

  // temp! XXX
  std::complex<double> ffc(std::size_t i, std::size_t j) const {
    static_assert(std::is_same<T, LinAlg::Matrix<std::complex<double>>>::value,
                  "Can only call ffc from Complex GMatrix!");
    return ff.get_copy(i, j);
  }
  gsl_complex &ffc2(std::size_t i, std::size_t j) const {
    static_assert(std::is_same<T, LinAlg::Matrix<std::complex<double>>>::value,
                  "Can only call ffc from Complex GMatrix!");
    return ff[i][j]; // pointer? reference? None! KILL
  }

  //! Can add/subtract matrices (in place)
  GreenMatrix<T> &operator+=(const GreenMatrix<T> &rhs) {
    ff += rhs.ff;
    if (m_include_G) {
      fg += rhs.fg;
      gf += rhs.gf;
      gg += rhs.gg;
    }
    return *this;
  }
  [[nodiscard]] friend GreenMatrix<T> operator+(GreenMatrix<T> lhs,
                                                const GreenMatrix<T> &rhs) {
    return lhs += rhs;
  }
  GreenMatrix<T> &operator-=(const GreenMatrix<T> &rhs) {
    ff -= rhs.ff;
    if (m_include_G) {
      fg -= rhs.fg;
      gf -= rhs.gf;
      gg -= rhs.gg;
    }
    return *this;
  }
  [[nodiscard]] friend GreenMatrix<T> operator-(GreenMatrix<T> lhs,
                                                const GreenMatrix<T> &rhs) {
    return lhs -= rhs;
  }

  GreenMatrix<T> &operator*=(double x) {
    ff *= x;
    if (m_include_G) {
      fg *= x;
      gf *= x;
      gg *= x;
    }
    return *this;
  }
  [[nodiscard]] friend GreenMatrix<T> operator*(double x, GreenMatrix<T> rhs) {
    // rhs *= x;
    return rhs *= x;
  }

  GreenMatrix<T> &operator*=(const std::complex<double> &x) {
    static_assert(std::is_same<T, LinAlg::Matrix<std::complex<double>>>::value,
                  "Can only call *=Complex from Complex GMatrix!");
    ff *= x;
    if (m_include_G) {
      fg *= x;
      gf *= x;
      gg *= x;
    }
    return *this;
  }

  // [[nodiscard]] inline friend GreenMatrix<T>
  // operator*(const std::complex<double> &x, GreenMatrix<T> rhs) {
  //   return rhs *= x;
  // }

  [[nodiscard]] inline friend auto operator*(const std::complex<double> &x,
                                             GreenMatrix<T> rhs) {
    if constexpr (std::is_same<T,
                               LinAlg::Matrix<std::complex<double>>>::value) {
      return rhs *= x;
    } else {
      return rhs.make_complex({1.0, 0.0}) *= x;
      // return c_rhs *= x;
    }
  }

  //! Matrix multplication (in place): Gij -> \sum_k Gik*Bkj
  GreenMatrix<T> &operator*=(const GreenMatrix<T> &b) {
    if (m_include_G) {
      auto tmp = ff;
      ff = ff * b.ff + fg * b.gf;
      fg = tmp * b.fg + fg * b.gg; // ff
      tmp = gf;
      gf = gf * b.ff + gg * b.gf;
      gg = tmp * b.fg + gg * b.gg; // gf
    } else {
      ff = ff * b.ff;
    }
    return *this;
  }
  [[nodiscard]] friend GreenMatrix<T> operator*(GreenMatrix<T> lhs,
                                                const GreenMatrix<T> &rhs) {
    return lhs *= rhs;
  }

  //! Inversion (in place)
  GreenMatrix<T> &invert() {
    ff.invert();
    if (m_include_G) {
      // std::cout << "\n XXX \n Warning: Inversion not tesed!\n";
      const auto &ai = ff; // already inverted
      const auto &b = fg;
      const auto &c = gf;
      const auto &d = gg;
      const auto cai = c * ai;
      const auto dmcaib = (d - cai * b).invert();
      const auto aib_dmcaib = ai * b * dmcaib;
      ff += aib_dmcaib * cai;
      if constexpr (std::is_same<T,
                                 LinAlg::Matrix<std::complex<double>>>::value) {
        const auto neg_one = std::complex<double>{-1.0, 0.0};
        fg = neg_one * aib_dmcaib;
        gf = neg_one * dmcaib * cai;
      } else {
        fg = -1.0 * aib_dmcaib;
        gf = -1.0 * dmcaib * cai;
      }
      gg = dmcaib;
    }
    return *this;
  }
  [[nodiscard]] GreenMatrix<T> inverse() const {
    auto out = *this; //
    return out.invert();
  }

  auto max_el() const {
    double maxre = 0.0;
    double maxim = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
      for (std::size_t j = 0; j < size; ++j) {
        // const auto [re, im] = ff[i][j];
        const auto re = ff[i][j].real();
        const auto im = ff[i][j].imag();
        if (std::abs(re) > std::abs(maxre))
          maxre = re;
        if (std::abs(im) > std::abs(maxim))
          maxim = im;
      }
    }
    if (m_include_G) {
      for (const auto tmp : {&fg, &gf, &gg}) {
        for (std::size_t i = 0; i < G_size; ++i) {
          for (std::size_t j = 0; j < G_size; ++j) {
            // const auto [re, im] = (*tmp)[i][j];
            const auto re = (*tmp)[i][j].real();
            const auto im = (*tmp)[i][j].imag();
            if (std::abs(re) > std::abs(maxre))
              maxre = re;
            if (std::abs(im) > std::abs(maxim))
              maxim = im;
          }
        }
      }
    }
    return std::pair{maxre, maxim};
  };

  //! Multiply elements (in place): Gij -> Gij*Bij
  GreenMatrix<T> &mult_elements_by(const GreenMatrix<T> &rhs) {
    ff.mult_elements_by(rhs.ff);
    if (m_include_G) {
      fg.mult_elements_by(rhs.fg);
      gf.mult_elements_by(rhs.gf);
      gg.mult_elements_by(rhs.gg);
    }
    return *this;
  }
  //! Multiply elements (new matrix): Gij = Aij*Bij
  [[nodiscard]] friend GreenMatrix<T> mult_elements(GreenMatrix<T> lhs,
                                                    const GreenMatrix<T> &rhs) {
    lhs.mult_elements_by(rhs);
    return lhs;
  }

  //! Return the real-part of a complex GreenMatrix (by copy)
  [[nodiscard]] GreenMatrix<LinAlg::Matrix<double>> get_real() const {
    static_assert(std::is_same<T, LinAlg::Matrix<std::complex<double>>>::value,
                  "Can only call get_real from Complex GMatrix!");
    auto gmat = GreenMatrix<LinAlg::Matrix<double>>(size, m_include_G);
    gmat.ff = ff.real();
    if (m_include_G) {
      gmat.fg = fg.real();
      gmat.gf = gf.real();
      gmat.gg = gg.real();
    }
    return gmat;
  }

  //! Return the imaginary-part of a complex GreenMatrix (by copy)
  [[nodiscard]] GreenMatrix<LinAlg::Matrix<double>> get_imaginary() const {
    static_assert(std::is_same<T, LinAlg::Matrix<std::complex<double>>>::value,
                  "Can only call get_imaginary from Complex GMatrix!");
    auto gmat = GreenMatrix<LinAlg::Matrix<double>>(size, m_include_G);
    gmat.ff = ff.imag();
    if (m_include_G) {
      gmat.fg = fg.imag();
      gmat.gf = gf.imag();
      gmat.gg = gg.imag();
    }
    return gmat;
  }

  //! Construct a complex GreenMatrix (C) from a Real one (R), by C = x*R
  [[nodiscard]] GreenMatrix<LinAlg::Matrix<std::complex<double>>>
  make_complex(const std::complex<double> &x = {1.0, 0.0}) const {
    static_assert(std::is_same<T, LinAlg::Matrix<double>>::value,
                  "Can only call make_complex from Real GMatrix!");
    auto gmat =
        GreenMatrix<LinAlg::Matrix<std::complex<double>>>(size, m_include_G);
    gmat.ff = x * ff.complex();
    if (m_include_G) {
      gmat.fg = x * fg.complex();
      gmat.gf = x * gf.complex();
      gmat.gg = x * gg.complex();
    }
    return gmat;
  }
};

//******************************************************************************

using GMatrix = GreenMatrix<LinAlg::Matrix<double>>;
using ComplexGMatrix = GreenMatrix<LinAlg::Matrix<std::complex<double>>>;
using ComplexDouble = std::complex<double>;

} // namespace MBPT
