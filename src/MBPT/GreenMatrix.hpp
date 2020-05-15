#pragma once
#include "Maths/LinAlg_MatrixVector.hpp"
#include <type_traits>
namespace MBPT {

// enum class Gtype { Real, Complex };
// enum class Gincl { Fonly, FandG };

//******************************************************************************
//! Holds Green's fn operator of form: |ket><bra| [4x4 matrix of NxN matrix]
template <typename T> class GreenMatrix {
  bool include_G;
  std::size_t size;
  std::size_t G_size;

public:
  T ff, fg, gf, gg;

public:
  GreenMatrix(std::size_t in_size, bool in_include_G)
      : include_G(in_include_G),
        size(in_size),
        G_size(include_G ? size : 0),
        ff(size),
        fg(G_size),
        gf(G_size),
        gg(G_size) {
    zero();
  }

  //! Sets all matrix elements to zero
  void zero() {
    ff.zero();
    if (include_G) {
      fg.zero();
      gf.zero();
      gg.zero();
    }
  }

  //! Makes 1
  void make_identity() {
    ff.make_identity();
    if (include_G) {
      fg.make_identity();
      gf.make_identity();
      gg.make_identity();
    }
  }

  GreenMatrix<T> &plusIdent(double a = 1.0) {
    ff.plusIdent(a);
    if (include_G) {
      gg.plusIdent(a);
    }
    return *this;
  }

  GreenMatrix<T> &plusIdent(double re, double im) {
    static_assert(std::is_same<T, LinAlg::ComplexSqMatrix>::value,
                  "Can only call plusIdent(x,y) from Complex GMatrix!");
    ff.plusIdent(re, im);
    if (include_G) {
      gg.plusIdent(re, im);
    }
    return *this;
  }

  void checkNaN() const {
    ff.checkNaN();
    if (include_G) {
      fg.checkNaN();
      gf.checkNaN();
      gg.checkNaN();
    }
  }

  // temp! XXX
  LinAlg::Complex<double> ffc(std::size_t i, std::size_t j) const {
    static_assert(std::is_same<T, LinAlg::ComplexSqMatrix>::value,
                  "Can only call ffc from Complex GMatrix!");
    return ff.get_copy(i, j);
  }

  //! Can add/subtract matrices (in place)
  GreenMatrix<T> &operator+=(const GreenMatrix<T> &rhs) {
    ff += rhs.ff;
    if (include_G) {
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
    if (include_G) {
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
    if (include_G) {
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

  GreenMatrix<T> &operator*=(const LinAlg::Complex<double> &x) {
    static_assert(std::is_same<T, LinAlg::ComplexSqMatrix>::value,
                  "Can only call *=Complex from Complex GMatrix!");
    ff *= x;
    if (include_G) {
      fg *= x;
      gf *= x;
      gg *= x;
    }
    return *this;
  }
  [[nodiscard]] friend GreenMatrix<T>
  operator*(const LinAlg::Complex<double> &x, GreenMatrix<T> rhs) {
    // rhs *= x;
    return rhs *= x;
  }

  //! Matrix multplication (in place): Gij -> \sum_k Gik*Bkj
  GreenMatrix<T> &operator*=(const GreenMatrix<T> &b) {
    if (include_G) {
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
    if (include_G) {
      std::cout << "\n XXX \n Warning: Inversion not tesed!\n";
      const auto &ai = ff; // already inverted
      const auto &b = fg;
      const auto &c = gf;
      const auto &d = gg;
      const auto cai = c * ai;
      const auto dmcaib = (d - cai * b).invert();
      const auto aib_dmcaib = ai * b * dmcaib;
      ff += aib_dmcaib * cai;
      if constexpr (std::is_same<T, LinAlg::ComplexSqMatrix>::value) {
        const auto neg_one = LinAlg::Complex<double>{-1.0, 0.0};
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

  //! Multiply elements (in place): Gij -> Gij*Bij
  GreenMatrix<T> &mult_elements_by(const GreenMatrix<T> &rhs) {
    ff.mult_elements_by(rhs.ff);
    if (include_G) {
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
  [[nodiscard]] GreenMatrix<LinAlg::SqMatrix> get_real() const {
    static_assert(std::is_same<T, LinAlg::ComplexSqMatrix>::value,
                  "Can only call get_real from Complex GMatrix!");
    auto gmat = GreenMatrix<LinAlg::SqMatrix>(size, include_G);
    gmat.ff = ff.real();
    if (include_G) {
      gmat.fg = fg.real();
      gmat.gf = gf.real();
      gmat.gg = gg.real();
    }
    return gmat;
  }

  //! Return the imaginary-part of a complex GreenMatrix (by copy)
  [[nodiscard]] GreenMatrix<LinAlg::SqMatrix> get_imaginary() const {
    static_assert(std::is_same<T, LinAlg::ComplexSqMatrix>::value,
                  "Can only call get_imaginary from Complex GMatrix!");
    auto gmat = GreenMatrix<LinAlg::SqMatrix>(size, include_G);
    gmat.ff = ff.imaginary();
    if (include_G) {
      gmat.fg = fg.imaginary();
      gmat.gf = gf.imaginary();
      gmat.gg = gg.imaginary();
    }
    return gmat;
  }

  //! Construct a complex GreenMatrix (C) from a Real one (R), by C = x*R
  [[nodiscard]] GreenMatrix<LinAlg::ComplexSqMatrix>
  make_complex(const LinAlg::Complex<double> &x) const {
    static_assert(std::is_same<T, LinAlg::SqMatrix>::value,
                  "Can only call make_complex from Real GMatrix!");
    auto gmat = GreenMatrix<LinAlg::ComplexSqMatrix>(size, include_G);
    gmat.ff = LinAlg::ComplexSqMatrix::make_complex(x, ff);
    if (include_G) {
      gmat.fg = LinAlg::ComplexSqMatrix::make_complex(x, fg);
      gmat.gf = LinAlg::ComplexSqMatrix::make_complex(x, gf);
      gmat.gg = LinAlg::ComplexSqMatrix::make_complex(x, gg);
    }
    return gmat;
  }
}; // namespace MBPT

//******************************************************************************

} // namespace MBPT
