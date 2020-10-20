#include "StructureRad.hpp"
#include "DiracOperator/DiracOperator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <numeric>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
constexpr bool use_omp = true;
#else
constexpr bool use_omp = true;
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

namespace MBPT {

//******************************************************************************
StructureRad::StructureRad(const std::vector<DiracSpinor> &basis,
                           double en_core)
    : mBasis(basis),    // store a local copy (?)
      mCoreEn(en_core), // Fa.en < en_core ==> core state!
      mY(basis.front().rgrid, &mBasis) {
  // Sort into core/excited; store pointers
  std::cout << __LINE__ << std::endl;
  for (const auto &Fn : mBasis) {
    if (Fn.en < mCoreEn)
      mCore.push_back(Fn);
    else
      mExcited.push_back(Fn);
  }
}

//******************************************************************************
double StructureRad::srTB(const DiracOperator::TensorOperator *const h,
                          const DiracSpinor &w, const DiracSpinor &v,
                          double omega) const {

  if (h->isZero(w.k, v.k))
    return 0.0;

  const auto k = h->rank();

  std::vector<double> srt(mCore.size());

#pragma omp parallel for
  for (auto ic = 0ul; ic < mCore.size(); ++ic) {
    const auto &c = mCore[ic];
    for (const auto &r : mExcited) {

      if (h->isZero(c.k, r.k))
        continue;

      const auto inv_erc_pw = 1.0 / (r.en - c.en + omega);
      const auto inv_erc_mw = 1.0 / (r.en - c.en - omega);

      const auto t_cr = h->reducedME(c, r);
      const auto t_rc = h->reducedME(r, c); // blah
      // XXX Can include RPA here also?

      const auto T_wrvc = t1(k, w, r, v, c) + t2(k, w, r, v, c) +
                          t3(k, w, r, v, c) + t4(k, w, r, v, c);

      const auto B_wcvr = b1(k, w, c, v, r) + b2(k, w, c, v, r) +
                          b3(k, w, c, v, r) + b4(k, w, c, v, r);

      srt[ic] += (t_cr * T_wrvc * inv_erc_pw + 0 * t_rc * B_wcvr * inv_erc_mw);
    }
  }

  return std::accumulate(cbegin(srt), cend(srt), 0.0);
}

//******************************************************************************
double StructureRad::t1(int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  // (-1)^{w-v+K}
  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto inv_e_rwab = 1.0 / (r.en + w.en - a.en - b.en);

      const auto [minL, maxL] = mY.k_minmax(a, c);

      for (int u = 0; u <= max_u; ++u) {

        const auto z = mY.Zk(u, w, r, b, a);
        if (Angular::zeroQ(z))
          continue;

        for (int l = minL; l <= maxL; ++l) {
          const auto sj1 = mY.sixj(w, v, K, l, u, b);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, a);
          if (Angular::zeroQ(sj2))
            continue;

          const auto x = mY.Xk(l, b, a, v, c); // *
          if (Angular::zeroQ(x))
            continue;

          t += sj1 * sj2 * z * x * inv_e_rwab;
        }
      }
    }
  }
  return s * t;

  // * note: It would seem that it should be faster to calculate x
  // and store in array first. Turns out, this is not the case.
  // Probably, because we take advantage of selection rules, and only calc x
  // when the other terms are all non-zero!
}

//******************************************************************************
double StructureRad::t2(int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnav = 1.0 / (r.en + n.en - a.en - v.en);

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, r, c, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2fu =
            Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u) / (2 * u + 1);

        const auto z1 = mY.Zk(u, a, w, n, c);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, n, r, a, v);
        if (Angular::zeroQ(z2))
          continue;

        t += s2fu * sj * z1 * z2 * inv_e_rnav;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::t3(int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nwac = 1.0 / (n.en + w.en - a.en - c.en);

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, r, c, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2fu =
            Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u) / (2 * u + 1);

        const auto z1 = mY.Zk(u, n, w, a, c);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, a, r, n, v);
        if (Angular::zeroQ(z2))
          continue;

        t += s2fu * sj * z1 * z2 * inv_e_nwac;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::t4(int K, const DiracSpinor &w, const DiracSpinor &r,
                        const DiracSpinor &v, const DiracSpinor &c) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_l = mY.Ck().max_k();

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto inv_e_nmcv = 1.0 / (n.en + m.en - c.en - v.en);

      const auto [minU, maxU] = mY.k_minmax(r, m);

      for (int u = minU; u <= maxU; ++u) {

        const auto x = mY.Xk(u, w, r, n, m);
        if (Angular::zeroQ(x))
          continue;

        for (int l = 0; l <= max_l; ++l) {

          const auto z = mY.Zk(l, n, m, v, c);
          if (Angular::zeroQ(z))
            continue;

          const auto sj1 = mY.sixj(w, v, K, l, u, n);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          t += sj1 * sj2 * x * z * inv_e_nmcv;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
//******************************************************************************
double StructureRad::b1(int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &b : mCore) {

      const auto inv_e_rvab = 1.0 / (r.en + v.en - a.en - b.en);

      const auto [minL, maxL] = mY.k_minmax(c, a);

      for (int u = 0; u <= max_u; ++u) {

        const auto z = mY.Zk(u, b, a, v, r);
        if (Angular::zeroQ(z))
          continue;

        for (int l = minL; l <= maxL; ++l) {
          const auto sj1 = mY.sixj(w, v, K, u, l, b);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, a);
          if (Angular::zeroQ(sj2))
            continue;

          const auto x = mY.Xk(l, w, c, b, a);
          if (Angular::zeroQ(x))
            continue;

          t += sj1 * sj2 * z * x * inv_e_rvab;
        }
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b2(int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_rnaw = 1.0 / (r.en + n.en - a.en - w.en);

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, c, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2fu =
            Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u) / (2 * u + 1);

        const auto z1 = mY.Zk(u, n, c, a, v);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, a, w, n, r);
        if (Angular::zeroQ(z2))
          continue;

        t += s2fu * sj * z1 * z2 * inv_e_rnaw;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b3(int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(r.twoj() - c.twoj() + 2 * K);
  const auto max_u = mY.Ck().max_k();

  for (const auto &a : mCore) {
    for (const auto &n : mExcited) {

      const auto inv_e_nvac = 1.0 / (n.en + v.en - a.en - c.en);

      for (int u = 0; u <= max_u; ++u) {

        const auto sj = mY.sixj(w, v, K, c, r, u);
        if (Angular::zeroQ(sj))
          continue;

        const auto s2fu =
            Angular::neg1pow_2(a.twoj() - n.twoj() + 2 * u) / (2 * u + 1);

        const auto z1 = mY.Zk(u, a, c, n, v);
        if (Angular::zeroQ(z1))
          continue;

        const auto z2 = mY.Zk(u, n, w, a, r);
        if (Angular::zeroQ(z2))
          continue;

        t += s2fu * sj * z1 * z2 * inv_e_nvac;
      }
    }
  }
  return s * t;
}

//******************************************************************************
double StructureRad::b4(int K, const DiracSpinor &w, const DiracSpinor &c,
                        const DiracSpinor &v, const DiracSpinor &r) const {

  double t = 0.0;

  const auto s = Angular::neg1pow_2(w.twoj() - v.twoj() + 2 * K);
  const auto max_l = mY.Ck().max_k();

  for (const auto &n : mExcited) {
    for (const auto &m : mExcited) {

      const auto inv_e_nmcw = 1.0 / (n.en + m.en - c.en - w.en);

      const auto [minU, maxU] = mY.k_minmax(m, r);

      for (int u = minU; u <= maxU; ++u) {

        const auto x = mY.Xk(u, n, m, v, r);
        if (Angular::zeroQ(x))
          continue;

        for (int l = 0; l <= max_l; ++l) {

          const auto z = mY.Zk(l, w, c, n, m);
          if (Angular::zeroQ(z))
            continue;

          const auto sj1 = mY.sixj(w, v, K, u, l, n);
          if (Angular::zeroQ(sj1))
            continue;
          const auto sj2 = mY.sixj(r, c, K, l, u, m);
          if (Angular::zeroQ(sj2))
            continue;

          t += sj1 * sj2 * x * z * inv_e_nmcw;
        }
      }
    }
  }
  return s * t;
}

} // namespace MBPT
