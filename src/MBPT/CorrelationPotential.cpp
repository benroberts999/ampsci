#include "MBPT/CorrelationPotential.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Physics/AtomData.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "Wavefunction/Wavefunction.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

//******************************************************************************
// Helper function:
static inline int find_max_tj(const std::vector<DiracSpinor> &core,
                              const std::vector<DiracSpinor> &excited) {
  if (core.empty() || excited.empty())
    return 0;
  // returns maximum value of 2*j in {core,excited}
  auto maxtj1 =
      std::max_element(core.cbegin(), core.cend(), DiracSpinor::comp_j)->twoj();
  auto maxtj2 =
      std::max_element(excited.cbegin(), excited.cend(), DiracSpinor::comp_j)
          ->twoj();
  return std::max(maxtj1, maxtj2);
}
//******************************************************************************

namespace MBPT {

//******************************************************************************
//******************************************************************************
GMatrix::GMatrix(int in_size)
    : size(in_size), ff(size), fg(size), gf(size), gg(size) {
  zero();
}
void GMatrix::zero() {
  ff.zero();
  fg.zero();
  gf.zero();
  gg.zero();
}
GMatrix &GMatrix::operator+=(const GMatrix &rhs) {
  ff += rhs.ff;
  fg += rhs.fg;
  gf += rhs.gf;
  gg += rhs.gg;
  return *this;
}
GMatrix &GMatrix::operator-=(const GMatrix &rhs) {
  ff -= rhs.ff;
  fg -= rhs.fg;
  gf -= rhs.gf;
  gg -= rhs.gg;
  return *this;
}

//******************************************************************************
//******************************************************************************
CorrelationPotential::CorrelationPotential(
    const Grid &gr, const std::vector<DiracSpinor> &core,
    const std::vector<DiracSpinor> &excited, const int in_stride,
    const std::vector<double> &en_list, const std::string &in_fname)
    : p_gr(&gr),
      m_core(core),
      m_excited(excited),
      m_yec(&gr, &m_excited, &m_core),
      m_maxk(find_max_tj(core, excited)),
      m_6j(m_maxk, m_maxk),
      stride(in_stride) {
  auto sp = IO::Profile::safeProfiler(__func__);

  std::cout << "\nCorrelation potential (Sigma)\n";

  std::cout << "(Including FF";
  if (include_FG)
    std::cout << ", FG";
  if (include_GG)
    std::cout << ", GG";
  std::cout << ")\n";

  const auto fname = in_fname == "" ? "" : in_fname + ".Sigma";
  if (fname != "" && IO::FRW::file_exists(fname)) {
    read_write(fname, IO::FRW::read);
  } else if (!en_list.empty()) {
    setup_subGrid();
    form_Sigma(en_list, fname);
  }
}

//******************************************************************************
void CorrelationPotential::setup_subGrid() {
  // Form the "Sigma sub-grid"
  const auto &rvec = p_gr->r;
  // const double rmin = m_core.empty() ? 1.0e-4 : m_core.front().r0();
  // const double rmax = m_core.empty() ? 30.0 : m_core.front().rinf();
  const double rmin = 1.0e-4;
  const double rmax = 30.0;

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
  printf(
      "Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, stride=%i]\n",
      r_stride.front(), r_stride.back(), stride_points, imin, stride);
}

//******************************************************************************
void CorrelationPotential::addto_G(GMatrix *Gmat, const DiracSpinor &ket,
                                   const DiracSpinor &bra,
                                   const double f) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Adds (f)*|ket><bra| to G matrix
  // G_ij = f * Q_i * W_j
  // Q = Q(1) = ket, W = W(2) = bra
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
                                             const DiracSpinor &Fv) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // lambda is fitting factor (just scales Sigma|v>)
  // Sigma|v> = int G(r1,r2)*v(r2) dr2
  // (S|v>)_i = sum_j G_ij v_j drdu_j du
  // nb: G is on sub-grid, |v> and S|v> on full-grid. Use interpolation

  const auto ki = std::size_t(Fv.k_index());
  const auto lambda = ki >= m_lambda_kappa.size() ? 1.0 : m_lambda_kappa[ki];

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
  // Interpolate from sub-grid to full grid
  SigmaFv.f = Interpolator::interpolate(r_stride, f, gr.r);
  if constexpr (include_FG || include_GG) {
    SigmaFv.g = Interpolator::interpolate(r_stride, g, gr.r);
  }

  return SigmaFv;
}

//******************************************************************************
void CorrelationPotential::form_Sigma(const std::vector<double> &en_list,
                                      const std::string &fname) {
  auto sp = IO::Profile::safeProfiler(__func__);

  Sigma_kappa.resize(en_list.size(), stride_points);

  if (m_core.empty() || m_excited.empty()) {
    std::cerr << "\nERROR 162 in form_Sigma: No basis! Sigma will just be 0!\n";
    return;
  }

  std::cout << "Forming correlation potential for:\n";
  for (auto ki = 0ul; ki < en_list.size(); ki++) {
    const auto kappa = Angular::kappaFromIndex(int(ki));

    // if v.kappa > basis, then Ck angular factor won't exist!
    auto tj = Angular::twojFromIndex(int(ki));
    if (tj > m_yec.Ck().max_tj())
      continue;

    printf(" k=%2i %6s at en=%8.5f.. ", kappa,
           AtomData::kappa_symbol(kappa).c_str(), en_list[ki]);
    std::cout << std::flush;
    fill_Sigma_k_Gold(&Sigma_kappa[ki], kappa, en_list[ki]);
    // find lowest excited state, output <v|S|v> energy shift:
    auto find_kappa = [=](const auto &a) { return a.k == kappa; };
    const auto vk =
        std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
    if (vk != cend(m_excited))
      std::cout << "de=" << *vk * Sigma2Fv(*vk);
    std::cout << "\n";
  }

  // write to disk
  if (fname != "")
    read_write(fname, IO::FRW::write);
}

//******************************************************************************
DiracSpinor CorrelationPotential::Sigma2Fv(const DiracSpinor &v) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Find correct G matrix (corresponds to kappa_v), return Sigma|v>
  // If Sigma_kappa doesn't exist, returns |0>
  auto kappa_index = std::size_t(Angular::indexFromKappa(v.k));
  if (kappa_index >= Sigma_kappa.size())
    return 0.0 * v;
  return Sigma_G_Fv(Sigma_kappa[kappa_index], v);
}

//******************************************************************************
void CorrelationPotential::fill_Sigma_k_Gold(GMatrix *Gmat, const int kappa,
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

  // Just for safety, should already be zero (unless re-calcing G)
  Gmat->ff.zero();
  Gmat->fg.zero();
  Gmat->gf.zero();
  Gmat->gg.zero();

  if (m_core.empty())
    return;
  const auto &gr = *(m_core.front().p_rgrid);

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
  // Calculates <Fv|Sigma|Fw> from scratch, at Fv energy [full grid + fg+gg]
  if (v.k != w.k)
    return 0.0;

  const auto &Ck = m_yec.Ck();

  // if v.kappa > basis, then Ck angular factor won't exist!
  if (v.twoj() > Ck.max_tj())
    return 0.0;

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
          const auto Qkv = Coulomb::Qk_abcd(v, a, m, n, k, yknb, Ck);
          if (Qkv == 0.0)
            continue;
          const auto Qkw =
              (&v == &w) ? Qkv : Coulomb::Qk_abcd(w, a, m, n, k, yknb, Ck);
          const auto &ybm = m_yec.get_y_ab(m, a);
          const auto Pkw = Coulomb::Pk_abcd(w, a, m, n, k, ybm, Ck, m_6j);
          const auto dele = v.en + a.en - m.en - n.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qkw + Pkw)) * Qkv;
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_core) {
          const auto Qkv = Coulomb::Qk_abcd(v, n, b, a, k, yknb, Ck);
          if (Qkv == 0.0)
            continue;
          const auto Qkw =
              (&v == &w) ? Qkv : Coulomb::Qk_abcd(w, n, b, a, k, yknb, Ck);
          const auto &yna = m_yec.get_y_ab(n, b);
          const auto Pkw = Coulomb::Pk_abcd(w, n, b, a, k, yna, Ck, m_6j);
          const auto dele = v.en + n.en - b.en - a.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qkw + Pkw)) * Qkv;
        } // b

      } // k
    }   // n
  }     // a

  return std::accumulate(delta_a.cbegin(), delta_a.cend(), 0.0);
}

//******************************************************************************
void CorrelationPotential::read_write(const std::string &fname,
                                      IO::FRW::RoW rw) {
  auto rw_str = rw == IO::FRW::write ? "Writing to " : "Reading from ";
  std::cout << rw_str << "Sigma file: " << fname << " ... " << std::flush;

  std::fstream iofs;
  IO::FRW::open_binary(iofs, fname, rw);

  // // write/read some grid parameters - just to check
  {
    double r0 = rw == IO::FRW::write ? p_gr->r0 : 0;
    double rmax = rw == IO::FRW::write ? p_gr->rmax : 0;
    double b = rw == IO::FRW::write ? p_gr->b : 0;
    std::size_t pts = rw == IO::FRW::write ? p_gr->num_points : 0;
    rw_binary(iofs, rw, r0, rmax, b, pts);
    if (rw == IO::FRW::read) {
      const bool grid_ok = std::abs((r0 - p_gr->r0) / r0) < 1.0e-6 &&
                           std::abs(rmax - p_gr->rmax) < 0.001 &&
                           (b - p_gr->b) < 0.001 && pts == p_gr->num_points;
      if (!grid_ok) {
        std::cerr << "\nFAIL 335 in read_write Sigma: Grid mismatch\n"
                  << "Have in file:\n"
                  << r0 << ", " << rmax << " w/ N=" << pts << ", b=" << b
                  << ", but expected:\n"
                  << p_gr->r0 << ", " << p_gr->rmax
                  << " w/ N=" << p_gr->num_points << ", b=" << p_gr->b << "\n";
        std::abort(); // abort?
      }
    }
  }

  // Sub-grid:
  rw_binary(iofs, rw, stride_points, imin, stride);
  if (rw == IO::FRW::read) {
    r_stride.resize(std::size_t(stride_points));
  }
  for (auto i = 0ul; i < r_stride.size(); ++i) {
    rw_binary(iofs, rw, r_stride[i]);
  }

  // Number of kappas (number of Sigma/G matrices)
  std::size_t num_kappas = rw == IO::FRW::write ? Sigma_kappa.size() : 0;
  rw_binary(iofs, rw, num_kappas);
  if (rw == IO::FRW::read) {
    Sigma_kappa.resize(num_kappas, stride_points);
  }

  // Check if include FG/GG written. Note: doesn't matter if mis-match?!
  auto incl_fg = rw == IO::FRW::write ? include_FG : 0;
  auto incl_gg = rw == IO::FRW::write ? include_GG : 0;
  rw_binary(iofs, rw, incl_fg, incl_gg);

  // Read/Write G matrices
  for (auto &Gk : Sigma_kappa) {
    for (int i = 0; i < stride_points; ++i) {
      for (int j = 0; j < stride_points; ++j) {
        rw_binary(iofs, rw, Gk.ff[i][j]);
        if (incl_fg) {
          rw_binary(iofs, rw, Gk.fg[i][j]);
          rw_binary(iofs, rw, Gk.gf[i][j]);
        }
        if (incl_gg) {
          rw_binary(iofs, rw, Gk.gg[i][j]);
        }
      }
    }
  }
  std::cout << "... done.\n";
  printf(
      "Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, stride=%i]\n",
      r_stride.front(), r_stride.back(), stride_points, imin, stride);
}

//******************************************************************************
void CorrelationPotential::print_scaling() const {
  if (!m_lambda_kappa.empty()) {
    std::cout << "Scaled Sigma, with: lambda_kappa = ";
    for (const auto &l : m_lambda_kappa) {
      std::cout << l << ", ";
    }
    std::cout << "\n";
  }
}

} // namespace MBPT
