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
  GreenMatrix<T> &operator-=(const GreenMatrix<T> &rhs) {
    ff -= rhs.ff;
    if (include_G) {
      fg -= rhs.fg;
      gf -= rhs.gf;
      gg -= rhs.gg;
    }
    return *this;
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
      fg = -1.0 * aib_dmcaib;
      gf = -1.0 * dmcaib * cai;
      gg = dmcaib;
    }
    return *this;
  }

  //! Multiply elements (in place): Gij -> Gij*Bij
  void mult_elements_by(const GreenMatrix<T> &rhs) {
    ff.mult_elements_by(rhs.ff);
    if (include_G) {
      fg.mult_elements_by(rhs.fg);
      gf.mult_elements_by(rhs.gf);
      gg.mult_elements_by(rhs.gg);
    }
  }

  //! Return the real-part of a complex GreenMatrix (by copy)
  GreenMatrix<LinAlg::SqMatrix> get_real() const {
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
  GreenMatrix<LinAlg::SqMatrix> get_imaginary() const {
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
  GreenMatrix<LinAlg::ComplexSqMatrix>
  make_complex(const LinAlg::Complex<double> &x) const {
    static_assert(std::is_same<T, LinAlg::SqMatrix>::value,
                  "Can only call make_complex from Real GMatrix!");
    auto gmat = GreenMatrix<LinAlg::ComplexSqMatrix>(size, include_G);
    gmat.ff = make_complex(x, ff);
    if (include_G) {
      gmat.fg = make_complex(x, fg);
      gmat.gf = make_complex(x, gf);
      gmat.gg = make_complex(x, gg);
    }
    return gmat;
  }
};

//******************************************************************************

} // namespace MBPT
