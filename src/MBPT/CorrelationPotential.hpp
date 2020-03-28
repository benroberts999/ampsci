#pragma once
#include "Angular/Angular_tables.hpp"
#include "Coulomb/YkTable.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <vector>
class Grid;
// namespace Angular {
// class Ck_ab;
// class SixJ;
// } // namespace Angular
// namespace Coulomb {
// class YkTable;
// } // namespace Coulomb

namespace MBPT {

struct GMatrix {
  GMatrix(int size) : ff(size), fg(size), gf(size), gg(size) {
    ff.zero();
    fg.zero();
    gf.zero();
    gg.zero();
  }
  LinAlg::SqMatrix ff;
  LinAlg::SqMatrix fg;
  LinAlg::SqMatrix gf;
  LinAlg::SqMatrix gg;

  // Check SqMatrix impl, might be slow?
  void operator+=(const GMatrix &rhs) {
    ff += rhs.ff;
    fg += rhs.fg;
    gf += rhs.gf;
    gg += rhs.gg;
  }
  void operator-=(const GMatrix &rhs) {
    ff -= rhs.ff;
    fg -= rhs.fg;
    gf -= rhs.gf;
    gg -= rhs.gg;
  }
};

//******************************************************************************
class CorrelationPotential {
public:
  // CorrelationPotential(const std::vector<DiracSpinor> &core,
  //                      const std::vector<DiracSpinor> &excited);

  CorrelationPotential(const std::vector<DiracSpinor> &core,
                       const std::vector<DiracSpinor> &excited,
                       const std::vector<double> &en_list = {});
  // XX Needs to know which are the core states!

  // No, give the class () operator!
  // One for each energy? One each kappa?
  // Probably: each kappa stored in here, each at single energy!

  DiracSpinor operator()(const DiracSpinor &Fv) const { return Sigma2Fv(Fv); }
  double operator()(const DiracSpinor &Fv, const DiracSpinor &Fw) const {
    return Sigma2vw(Fv, Fw);
  }

  // XXX Make these all const!
  DiracSpinor Sigma2Fv(const DiracSpinor &Fv) const;

  void fill_Gkappa(GMatrix *Gmat, const int kappa, const double en);

  double Sigma2vw(const DiracSpinor &Fv, const DiracSpinor &Fw) const;
  // double Sigma2vw(const DiracSpinor &Fv) const; // kill this, will be
  // confusing!

  void addto_Gff(GMatrix *Gmat, const DiracSpinor &ket, const DiracSpinor &bra,
                 const double f = 1.0) const {
    // G_ij = f * Q_i * [W_j drdu_j] * du
    const auto &gr = *(bra.p_rgrid);
    for (int i = 0; i < max_stride; ++i) {
      const auto si = (imin + i) * stride;
      for (int j = 0; j < max_stride; ++j) {
        const auto sj = (imin + j) * stride;
        Gmat->ff[i][j] += f * ket.f[si] * bra.f[sj];
        // Gmat->fg[i][j] += f * ket.f[si] * bra.g[sj];
        // Gmat->gf[i][j] += f * ket.g[si] * bra.f[sj];
        // Gmat->gg[i][j] += f * ket.g[si] * bra.g[sj];
      }
    }
  }

  DiracSpinor Sigma_G_Fv(const GMatrix &Gmat, const DiracSpinor &Fv) const {
    const auto &gr = *(Fv.p_rgrid);
    auto SigmaFv = DiracSpinor(0, Fv.k, gr);
    std::vector<double> f(max_stride);
    std::vector<double> g(max_stride);
    for (int i = 0; i < max_stride; ++i) {
      for (int j = 0; j < max_stride; ++j) {
        const auto sj = (imin + j) * stride;
        const auto dr = gr.drdu[sj] * gr.du * double(stride);
        f[i] += Gmat.ff[i][j] * Fv.f[sj] * dr;
        //
        // g[i] += Gmat.gf[i][j] * Fv.f[sj] * dr;
        //
        // f[i] += (Gmat.ff[i][j] * Fv.f[sj] + Gmat.fg[i][j] * Fv.g[sj]) * dr;
        // g[i] += (Gmat.gf[i][j] * Fv.f[sj] + Gmat.gg[i][j] * Fv.g[sj]) * dr;

        //
        // f[i] += (Gmat.ff[i][j] * Fv.f[sj] + Gmat.fg[i][j] * Fv.g[sj]) * dr;
        // g[i] += (Gmat.gf[i][j] * Fv.f[sj]) * dr;
      }
    }
    // SigmaFv.f = f;
    SigmaFv.f = Interpolator::interpolate(r_stride, f, gr.r);
    SigmaFv.g = Interpolator::interpolate(r_stride, g, gr.r);

    return SigmaFv;
  }

private:
  const int stride = 4;
  int max_stride{};
  std::size_t imin{};
  std::vector<double> r_stride = {};
  std::vector<DiracSpinor> m_core;
  std::vector<DiracSpinor> m_excited;
  Coulomb::YkTable m_yec; // constains Ck
  int m_maxk;
  Angular::SixJ m_6j;

  std::vector<GMatrix> G_kappa = {};
  std::vector<double> en_kappa;
};

} // namespace MBPT
