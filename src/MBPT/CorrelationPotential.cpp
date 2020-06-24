#include "MBPT/CorrelationPotential.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

namespace MBPT {

//******************************************************************************
//******************************************************************************
CorrelationPotential::CorrelationPotential(
    const HF::HartreeFock *const in_hf, const std::vector<DiracSpinor> &basis,
    const Sigma_params &sigp, const rgrid_params &subgridp)
    : p_gr(in_hf->rgrid),
      m_holes(copy_holes(basis, in_hf->get_core(), sigp.min_n_core)),
      m_excited(copy_excited(basis, in_hf->get_core())),
      m_yeh(p_gr, &m_excited, &m_holes),
      m_maxk(find_max_tj(in_hf->get_core(), basis)),
      m_6j(m_maxk, m_maxk),
      m_stride(subgridp.stride) {
  setup_subGrid(subgridp.r0, subgridp.rmax);
}

//******************************************************************************
void CorrelationPotential::setup_subGrid(double rmin, double rmax) {
  // Form the "Sigma sub-grid"
  const auto &rvec = p_gr->r;

  m_imin = 0;
  for (auto i = 0ul; i < rvec.size(); i += m_stride) {
    auto r = rvec[i];
    if (r < rmin) {
      m_imin++;
      continue;
    }
    if (r > rmax)
      break;
    m_subgrid_r.push_back(r);
  }

  m_subgrid_points = m_subgrid_r.size();
}

//******************************************************************************
std::vector<DiracSpinor>
CorrelationPotential::copy_holes(const std::vector<DiracSpinor> &basis,
                                 const std::vector<DiracSpinor> &core,
                                 int n_min) const {
  std::vector<DiracSpinor> t_holes;
  for (const auto &Fi : basis) {
    const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
    if (inCore && Fi.n >= n_min) {
      t_holes.push_back(Fi);
    }
  }
  return t_holes;
}
std::vector<DiracSpinor>
CorrelationPotential::copy_excited(const std::vector<DiracSpinor> &basis,
                                   const std::vector<DiracSpinor> &core) const {
  std::vector<DiracSpinor> t_excited;
  for (const auto &Fi : basis) {
    const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
    if (!inCore) {
      t_excited.push_back(Fi);
    }
  }
  return t_excited;
}

//******************************************************************************
void CorrelationPotential::addto_G(GMatrix *Gmat, const DiracSpinor &ket,
                                   const DiracSpinor &bra,
                                   const double f) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Adds (f)*|ket><bra| to G matrix
  // G_ij = f * Q_i * W_j
  // Q = Q(1) = ket, W = W(2) = bra
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto sj = ri_subToFull(j);
      Gmat->ff[i][j] += f * ket.f[si] * bra.f[sj];
      if constexpr (m_include_G) {
        Gmat->fg[i][j] += f * ket.f[si] * bra.g[sj];
        Gmat->gf[i][j] += f * ket.g[si] * bra.f[sj];
        Gmat->gg[i][j] += f * ket.g[si] * bra.g[sj];
      }
    } // j
  }   // i
}

//******************************************************************************
DiracSpinor CorrelationPotential::Sigma_G_Fv(const GMatrix &Gmat,
                                             const DiracSpinor &Fv) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // lambda is fitting factor (just scales Sigma|v>)
  // Sigma|v> = int G(r1,r2)*v(r2) dr2
  // (S|v>)_i = sum_j G_ij v_j drdu_j du
  // nb: G is on sub-grid, |v> and S|v> on full-grid. Use interpolation

  const auto ki = std::size_t(Fv.k_index());
  const auto lambda = ki >= m_lambda_kappa.size() ? 1.0 : m_lambda_kappa[ki];

  const auto &gr = *(Fv.rgrid);
  auto SigmaFv = DiracSpinor(0, Fv.k, gr);
  std::vector<double> f(m_subgrid_r.size());
  std::vector<double> g;
  if constexpr (m_include_G) {
    g.resize(m_subgrid_r.size());
  }
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto sj = ri_subToFull(j);
      const auto dr = gr.drdu[sj] * gr.du * double(m_stride);
      f[i] += Gmat.ff[i][j] * Fv.f[sj] * dr * lambda;

      if constexpr (m_include_G) {
        f[i] += Gmat.fg[i][j] * Fv.g[sj] * dr * lambda;
        g[i] += Gmat.gf[i][j] * Fv.f[sj] * dr * lambda;
        g[i] += Gmat.gg[i][j] * Fv.g[sj] * dr * lambda;
      }
    }
  }
  // Interpolate from sub-grid to full grid
  SigmaFv.f = Interpolator::interpolate(m_subgrid_r, f, gr.r);
  if constexpr (m_include_G) {
    SigmaFv.g = Interpolator::interpolate(m_subgrid_r, g, gr.r);
  }

  return SigmaFv;
}

//******************************************************************************
void CorrelationPotential::form_Sigma(const std::vector<double> &en_list,
                                      const std::string &fname) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  if (en_list.empty())
    return;

  m_Sigma_kappa.resize(en_list.size(), {m_subgrid_points, m_include_G});

  print_subGrid();

  std::cout << "Forming correlation potential (";
  std::cout << DiracSpinor::state_config(m_holes) << "/"
            << DiracSpinor::state_config(m_excited) << ") for:\n";
  for (auto ki = 0ul; ki < en_list.size(); ki++) {
    const auto kappa = Angular::kappaFromIndex(int(ki));

    // if v.kappa > basis, then Ck angular factor won't exist!
    auto tj = Angular::twojFromIndex(int(ki));
    if (tj > m_yeh.Ck().max_tj())
      continue;

    printf("k=%2i at en=%8.5f.. ", kappa, en_list[ki]);
    std::cout << std::flush;
    this->fill_Sigma_k(&m_Sigma_kappa[ki], kappa, en_list[ki]);

    // find lowest excited state, output <v|S|v> energy shift:
    auto find_kappa = [=](const auto &a) { return a.k == kappa; };
    const auto vk =
        std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
    if (vk != cend(m_excited)) {
      printf("%.6f", *vk * SigmaFv(*vk));
    }
    std::cout << "\n";
  }

  // write to disk
  if (fname != "")
    read_write(fname, IO::FRW::write);
}

//******************************************************************************
DiracSpinor CorrelationPotential::SigmaFv(const DiracSpinor &v) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Find correct G matrix (corresponds to kappa_v), return Sigma|v>
  // If m_Sigma_kappa doesn't exist, returns |0>
  auto kappa_index = std::size_t(Angular::indexFromKappa(v.k));
  if (kappa_index >= m_Sigma_kappa.size())
    return 0.0 * v;
  return Sigma_G_Fv(m_Sigma_kappa[kappa_index], v);
}

//******************************************************************************
double CorrelationPotential::SOEnergyShift(const DiracSpinor &v,
                                           const DiracSpinor &w) const {
  // Calculates <Fv|Sigma|Fw> from scratch, at Fw energy [full grid + fg+gg]
  if (v.k != w.k)
    return 0.0;

  const auto &Ck = m_yeh.Ck();

  // if v.kappa > basis, then Ck angular factor won't exist!
  // XXX Fix! extend Ck! (and sixj)
  if (v.twoj() > Ck.max_tj())
    return 0.0;

  std::vector<double> delta_a(m_holes.size());
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_holes.size(); ia++) {
    const auto &a = m_holes[ia];
    auto &del_a = delta_a[ia];
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = m_yeh.k_minmax(n, a);
      const auto max_k = std::min(m_maxk, kmax_nb);
      for (int k = kmin_nb; k <= max_k; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * v.twojp1();
        const auto &yknb = m_yeh.get_yk_ab(k, n, a);

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          const auto Qkv = Coulomb::Qk_abcd(v, a, m, n, k, yknb, Ck);
          if (Qkv == 0.0)
            continue;
          const auto Qkw =
              (&v == &w) ? Qkv : Coulomb::Qk_abcd(w, a, m, n, k, yknb, Ck);
          const auto &ybm = m_yeh.get_y_ab(m, a);
          const auto Pkw = Coulomb::Pk_abcd(w, a, m, n, k, ybm, Ck, m_6j);
          const auto dele = v.en + a.en - m.en - n.en;
          del_a += ((1.0 / dele / f_kkjj) * (Qkw + Pkw)) * Qkv;
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_holes) {
          const auto Qkv = Coulomb::Qk_abcd(v, n, b, a, k, yknb, Ck);
          if (Qkv == 0.0)
            continue;
          const auto Qkw =
              (&v == &w) ? Qkv : Coulomb::Qk_abcd(w, n, b, a, k, yknb, Ck);
          const auto &yna = m_yeh.get_y_ab(n, b);
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
bool CorrelationPotential::read_write(const std::string &fname,
                                      IO::FRW::RoW rw) {

  if (rw == IO::FRW::read && !IO::FRW::file_exists(fname))
    return false;

  const auto rw_str = rw == IO::FRW::write ? "Writing to " : "Reading from ";
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
                           std::abs(b - p_gr->b) < 0.001 &&
                           pts == p_gr->num_points;
      if (!grid_ok) {
        std::cout << "\nCannot read from:" << fname << ". Grid mismatch\n"
                  << "Read: " << r0 << ", " << rmax << " w/ N=" << pts
                  << ", b=" << b << ",\n but expected: " << p_gr->r0 << ", "
                  << p_gr->rmax << " w/ N=" << p_gr->num_points
                  << ", b=" << p_gr->b << "\n";
        std::cout << "Will calculate from scratch, + over-write file.\n";
        return false;
      }
    }
  }

  // read/write basis config:
  std::string basis_config =
      (rw == IO::FRW::write) ? DiracSpinor::state_config(m_excited) : "";
  rw_binary(iofs, rw, basis_config);
  // XXX Add info on method! (or, different filename!)

  // Sub-grid:
  rw_binary(iofs, rw, m_subgrid_points, m_imin, m_stride);
  if (rw == IO::FRW::read) {
    m_subgrid_r.resize(m_subgrid_points);
  }
  for (auto i = 0ul; i < m_subgrid_r.size(); ++i) {
    rw_binary(iofs, rw, m_subgrid_r[i]);
  }

  // Number of kappas (number of Sigma/G matrices)
  std::size_t num_kappas = rw == IO::FRW::write ? m_Sigma_kappa.size() : 0;
  rw_binary(iofs, rw, num_kappas);
  if (rw == IO::FRW::read) {
    m_Sigma_kappa.resize(num_kappas, {m_subgrid_points, m_include_G});
  }

  // Check if include FG/GG written. Note: doesn't matter if mis-match?!
  auto incl_g = rw == IO::FRW::write ? m_include_G : 0;
  rw_binary(iofs, rw, incl_g);

  // Read/Write G matrices
  for (auto &Gk : m_Sigma_kappa) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        rw_binary(iofs, rw, Gk.ff[i][j]);
        if (incl_g) {
          rw_binary(iofs, rw, Gk.fg[i][j]);
          rw_binary(iofs, rw, Gk.gf[i][j]);
          rw_binary(iofs, rw, Gk.gg[i][j]);
        }
      }
    }
  }
  std::cout << "done.\n";
  if (rw == IO::FRW::read) {
    std::cout << "Sigma basis: " << basis_config << "\n";
    print_subGrid();
  }
  return true;
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

//******************************************************************************
void CorrelationPotential::print_subGrid() const {
  if (m_subgrid_r.empty()) {
    std::cout << "No sub-grid??\n";
    return;
  }
  printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
         "stride=%i]\n",
         m_subgrid_r.front(), m_subgrid_r.back(), int(m_subgrid_points),
         int(m_imin), int(m_stride));
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
std::size_t CorrelationPotential::ri_subToFull(std::size_t i) const {
  return ((m_imin + i) * m_stride);
}
double CorrelationPotential::dr_subToFull(std::size_t i) const {
  return p_gr->drdu[ri_subToFull(i)] * p_gr->du * double(m_stride);
}

} // namespace MBPT
