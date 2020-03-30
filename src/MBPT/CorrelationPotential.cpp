#include "MBPT/CorrelationPotential.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

static inline int find_max_tj(const std::vector<DiracSpinor> &core,
                              const std::vector<DiracSpinor> &excited) {
  // returns maximum value of 2*j in {core,excited}
  auto maxtj1 =
      std::max_element(core.cbegin(), core.cend(), DiracSpinor::comp_j)->twoj();
  auto maxtj2 =
      std::max_element(excited.cbegin(), excited.cend(), DiracSpinor::comp_j)
          ->twoj();
  return std::max(maxtj1, maxtj2);
}

//******************************************************************************
//******************************************************************************
namespace MBPT {

CorrelationPotential::CorrelationPotential(
    const std::vector<DiracSpinor> &core,
    const std::vector<DiracSpinor> &excited, const std::vector<double> &en_list)
    : m_core(core),
      m_excited(excited),
      m_yec(m_excited.front().p_rgrid, &m_excited, &m_core),
      m_maxk(find_max_tj(core, excited)),
      m_6j(m_maxk, m_maxk) {
  auto sp = IO::Profile::safeProfiler(__func__);

  // Form the "Sigma sub-grid"
  const auto &rvec = m_excited.front().p_rgrid->r;
  const double rmin = m_core.front().r0();   // 2.0e-4;
  const double rmax = m_core.front().rinf(); // 30.0;
  imin = 0;
  for (auto i = 0; i < (int)rvec.size(); i += stride) {
    auto r = rvec[std::size_t(i)];
    if (r < rmin) {
      imin++;
      continue;
    }
    if (r > rmax)
      break;
    r_stride.push_back(r);
  }

  stride_points = int(r_stride.size());
  std::cout << "\nCorrelation potential (Sigma)\n";
  printf(
      "Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, stride=%i]\n",
      r_stride.front(), r_stride.back(), stride_points, imin, stride);

  std::cout << "(Including FF";
  if (include_FG)
    std::cout << ", FG";
  if (include_GG)
    std::cout << ", GG";
  std::cout << ")\n";

  if (!en_list.empty())
    form_Sigma(en_list);
}

//******************************************************************************
void CorrelationPotential::addto_G(GMatrix *Gmat, const DiracSpinor &ket,
                                   const DiracSpinor &bra,
                                   const double f) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Adds (f)*|ket><bra| to G matrix
  // G_ij = f * Q_i * W_j      // Q = Q(1) = ket, W = W(2) = bra
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  for (int i = 0; i < stride_points; ++i) {
    const auto si = static_cast<std::size_t>((imin + i) * stride);
    for (int j = 0; j < stride_points; ++j) {
      const auto sj = std::size_t((imin + j) * stride);
      Gmat->ff[i][j] += f * ket.f[si] * bra.f[sj];
      if constexpr (include_FG) {
        Gmat->fg[i][j] += f * ket.f[si] * bra.g[sj];
        Gmat->gf[i][j] += f * ket.g[si] * bra.f[sj];
      }
      if constexpr (include_GG) {
        Gmat->gg[i][j] += f * ket.g[si] * bra.g[sj];
      }
    } // j
  }   // i
}

//******************************************************************************
DiracSpinor CorrelationPotential::Sigma_G_Fv(const GMatrix &Gmat,
                                             const DiracSpinor &Fv,
                                             const double lambda) const {
  // lambda is fitting factor (just scales Sigma|v>)
  // Sigma|v> = int G(r1,r2)*v(r2) dr2
  // (S|v>)_i = sum_j G_ij v_j drdu_j du
  // nb: G is on sub-grid, |v> and S|v> on full-grid. Use interpolation
  auto sp = IO::Profile::safeProfiler(__func__);
  const auto &gr = *(Fv.p_rgrid);
  auto SigmaFv = DiracSpinor(0, Fv.k, gr);
  std::vector<double> f(r_stride.size());
  std::vector<double> g;
  if constexpr (include_FG || include_GG) {
    g.resize(r_stride.size());
  }
  for (int i = 0; i < stride_points; ++i) {
    const auto si = std::size_t(i);
    for (int j = 0; j < stride_points; ++j) {
      const auto sj = std::size_t((imin + j) * stride);
      const auto dr = gr.drdu[sj] * gr.du * double(stride);
      f[si] += Gmat.ff[i][j] * Fv.f[sj] * dr * lambda;

      if constexpr (include_FG) {
        f[si] += Gmat.fg[i][j] * Fv.g[sj] * dr * lambda;
        g[si] += Gmat.gf[i][j] * Fv.f[sj] * dr * lambda;
      }
      if constexpr (include_GG) {
        g[si] += Gmat.gg[i][j] * Fv.g[sj] * dr * lambda;
      }
    }
  }
  // Interpolate from sub-grid to fill-grid
  SigmaFv.f = Interpolator::interpolate(r_stride, f, gr.r);
  if constexpr (include_FG || include_GG) {
    SigmaFv.g = Interpolator::interpolate(r_stride, g, gr.r);
  }

  return SigmaFv;
}

//******************************************************************************
void CorrelationPotential::form_Sigma(const std::vector<double> &en_list) {
  auto sp = IO::Profile::safeProfiler(__func__);

  G_kappa.resize(en_list.size(), stride_points);

  std::cout << "Forming correlation potential for:\n";
  for (auto ki = 0ul; ki < en_list.size(); ki++) {
    auto kappa = Angular::kappaFromIndex(int(ki));
    printf(" k=%2i %6s at en=%8.5f.. ", kappa,
           AtomData::kappa_symbol(kappa).c_str(), en_list[ki]);
    std::cout << std::flush;
    fill_Gkappa(&G_kappa[ki], kappa, en_list[ki]);
    // find lowest excited state:
    auto find_kappa = [=](const auto &a) { return a.k == kappa; };
    const auto vk =
        std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
    std::cout << "de=" << *vk * Sigma2Fv(*vk) << "\n";
  }
}

//******************************************************************************
DiracSpinor CorrelationPotential::Sigma2Fv(const DiracSpinor &v,
                                           const double lambda) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Find correct G matrix (corresponds to kappa_v), return Sigma|v>
  // If G_kappa doesn't exist, returns |0>
  auto kappa_index = std::size_t(Angular::indexFromKappa(v.k));
  if (kappa_index >= G_kappa.size())
    return 0.0 * v;
  return Sigma_G_Fv(G_kappa[kappa_index], v, lambda);
}

//******************************************************************************
void CorrelationPotential::fill_Gkappa(GMatrix *Gmat, const int kappa,
                                       const double en) {
  auto sp = IO::Profile::safeProfiler(__func__);

  // Four second-order diagrams:
  // Diagram (a):
  // |Q^k_amn><Q^k_amn| / de_amn / [k][j]
  // Diagram (b) (exchange):
  // |Q^k_amn><P^k_amn| / de_amn / [k][j]
  // Diagram (c):
  // |Q^k_nba><Q^k_nba| / de_nba / [k][j]
  // Diagram (d) (exchange):
  // |Q^k_nba><P^k_nba| / de_nba / [k][j]
  // where:
  // All indeces are summed over,
  // a & b are core states, n & m are virtual excited states,
  // k is multipolarity [Coloulmb expansion]
  // de_xyz = e_v + e_x - e_y - e_z

  const auto &Ck = m_yec.Ck();

  const auto &gr = *(m_core.front().p_rgrid);

  // Just for safety, should already be zero (unless re-calcing G)
  Gmat->ff.zero();
  Gmat->fg.zero();
  Gmat->gf.zero();
  Gmat->gg.zero();

  // auto fk = [&](int k) { return 1.0; };

  // Note: get_yk_ab() must only be called with k for which y^k_ab exists!
  // Therefore, must use the k_minmax() function provided
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_core.size(); ia++) {
    const auto &a = m_core[ia];
    GMatrix G_a(stride_points);
    auto Qkv = DiracSpinor(0, kappa, gr); // re-use to reduce alloc'ns
    auto Pkv = DiracSpinor(0, kappa, gr); // re-use to reduce alloc'ns
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = m_yec.k_minmax(n, a);
      const auto max_k = std::min(m_maxk, kmax_nb);
      for (int k = kmin_nb; k <= max_k; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * (Angular::twoj_k(kappa) + 1);
        const auto &yknb = m_yec(k, n, a);

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          if (Ck(k, kappa, m.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, a, m, n, k, yknb, Ck);
          Coulomb::Pkv_bcd(&Pkv, a, m, n, k, m_yec(m, a), Ck, m_6j);
          const auto dele = en + a.en - m.en - n.en;
          const auto factor = 1.0 / (f_kkjj * dele);
          addto_G(&G_a, Qkv, Qkv + Pkv, factor);
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_core) {
          if (Ck(k, kappa, b.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, n, b, a, k, yknb, Ck);
          Coulomb::Pkv_bcd(&Pkv, n, b, a, k, m_yec(n, b), Ck, m_6j);
          const auto dele = en + n.en - b.en - a.en;
          const auto factor = 1.0 / (f_kkjj * dele);
          addto_G(&G_a, Qkv, Qkv + Pkv, factor);
        } // b

      } // k
    }   // n
#pragma omp critical(sumG)
    { *Gmat += G_a; }
  } // a
}

//******************************************************************************
double CorrelationPotential::Sigma2vw(const DiracSpinor &v,
                                      const DiracSpinor &w) const {
  if (v.k != w.k)
    return 0.0;

  const auto &Ck = m_yec.Ck();

  std::vector<double> delta_a(m_core.size());
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_core.size(); ia++) {
    const auto &a = m_core[ia];
    auto &del_a = delta_a[ia];
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = m_yec.k_minmax(n, a);
      const auto max_k = std::min(m_maxk, kmax_nb);
      for (int k = kmin_nb; k <= max_k; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * v.twojp1();
        const auto &yknb = m_yec.get_yk_ab(k, n, a);

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          const auto Qk = Coulomb::Qk_abcd(v, a, m, n, k, yknb, Ck);
          if (Qk == 0.0)
            continue;
          const auto &ybm = m_yec.get_y_ab(m, a);
          const auto Pk = Coulomb::Pk_abcd(v, a, m, n, k, ybm, Ck, m_6j);
          const auto dele = v.en + a.en - m.en - n.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qk + Pk)) * Qk;
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_core) {
          const auto Qk = Coulomb::Qk_abcd(v, n, b, a, k, yknb, Ck);
          if (Qk == 0.0)
            continue;
          const auto &yna = m_yec.get_y_ab(n, b);
          const auto Pk = Coulomb::Pk_abcd(v, n, b, a, k, yna, Ck, m_6j);
          const auto dele = v.en + n.en - b.en - a.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qk + Pk)) * Qk;
        } // b

      } // k
    }   // n
  }     // a

  return std::accumulate(delta_a.cbegin(), delta_a.cend(), 0.0);
}

} // namespace MBPT
