#include "MBPT/CorrelationPotential.hpp"
#include "Angular/CkTable.hpp"
#include "Coulomb/CoulombIntegrals.hpp"
#include "Coulomb/QkTable.hpp"
#include "Coulomb/YkTable.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "MBPT/Ladder.hpp" //?
#include "MBPT/RDMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Wavefunction/DiracSpinor.hpp"
#include "qip/omp.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

namespace MBPT {

//==============================================================================
//==============================================================================
CorrelationPotential::CorrelationPotential(
    const HF::HartreeFock *const in_hf, const std::vector<DiracSpinor> &basis,
    const Sigma_params &sigp, const rgrid_params &subgridp)
    : p_gr(in_hf->grid_sptr()),
      m_holes(copy_holes(basis, in_hf->core(), sigp.min_n_core)),
      m_excited(copy_excited(basis, in_hf->core())),
      m_yeh(m_excited, m_holes),
      m_maxk(std::max(DiracSpinor::max_tj(in_hf->core()),
                      DiracSpinor::max_tj(basis))),
      m_6j(2 * m_maxk),
      m_ladder_file(sigp.ladder_file),
      m_stride(subgridp.stride),
      m_include_G(sigp.include_G),
      m_fk(std::move(sigp.fk)),
      m_eta(std::move(sigp.etak)) {
  setup_subGrid(subgridp.r0, subgridp.rmax);
  if (m_ladder_file != "") {
    // read in ladder
    // m_lk = std::make_unique<Coulomb::LkTable>();
    m_lk = Coulomb::LkTable{};
    const auto read_lad = m_lk->read(m_ladder_file);
    if (!read_lad) {
      std::cout << "WARNING: couln't find ladder file: " << m_ladder_file
                << " - NOT READ ****\n";
    } else {
      std::cout << "Including ladder diagrams from: " << m_ladder_file << "\n";
    }
  }
}

//==============================================================================
void CorrelationPotential::setup_subGrid(double rmin, double rmax) {
  // Form the "Sigma sub-grid"
  // const auto &rvec = p_gr->r();
  m_imin = p_gr->getIndex(rmin);
  auto t_i_max = p_gr->getIndex(rmax);
  // ensure multiple of stride:
  m_subgrid_points = (t_i_max - m_imin) / m_stride + 1;
  m_subgrid_r.resize(m_subgrid_points); // kill! no longer used!
}

//==============================================================================
std::size_t CorrelationPotential::getSigmaIndex(int n, int kappa) const {

  //
  // If n=0, find first Sigma of that kappa

  assert(m_Sigma_kappa.size() == m_nk.size());
  std::size_t index = 0;
  for (const auto &[nn, kk, en] : m_nk) {
    (void)en;
    if (kk == kappa && (nn == n || n == 0))
      break;
    ++index;
  }
  const bool not_found = (index == m_nk.size());
  if (not_found) {
    // If not in list, return index for lowest n
    // Return largest instead?
    index = 0;
    for (const auto &[nn, kk, en] : m_nk) {
      (void)en;
      if (kk == kappa)
        break;
      ++index;
    }
  }
  return index;
  // nb: will return  index = m_nk.size() if kappa not found
}

const GMatrix *CorrelationPotential::getSigma(int n, int kappa) const {
  const auto is = getSigmaIndex(n, kappa);
  return (is < m_Sigma_kappa.size()) ? &m_Sigma_kappa[is] : nullptr;
}

//******************************************************************************
DiracSpinor CorrelationPotential::SigmaFv(const DiracSpinor &v,
                                          bool lad) const {
  // Find correct G matrix (corresponds to kappa_v), return Sigma|v>
  // If m_Sigma_kappa doesn't exist, returns |0>

  DiracSpinor SFv = 0.0 * v;
  const auto is = getSigmaIndex(v.n(), v.kappa());

  if (is < m_Sigma_kappa.size())
    SFv = act_G_Fv(m_Sigma_kappa[is], v);

  if (lad && m_lk && !m_ratio_ladder_method) {
    // SFv += Sigmal_Fv(v, m_yeh, *m_lk, m_holes, m_excited, m_fk, m_eta);
    SFv += Sigmal_Fv(v, m_yeh, *m_lk, m_holes, m_excited);
  }

  // Aply lambda, if exists:
  const auto lambda = is >= m_lambda_kappa.size() ? 1.0 : m_lambda_kappa[is];
  if (lambda != 1.0)
    SFv *= lambda;

  // if (is < m_Sigma_kappa.size())
  //   return lambda == 1.0 ? act_G_Fv(m_Sigma_kappa[is], v) :
  //                          lambda * act_G_Fv(m_Sigma_kappa[is], v);

  return SFv;
}

//==============================================================================
void CorrelationPotential::scale_Sigma(int n, int kappa, double lambda) {
  // XXX Careful; likely to be incorrect?? if given too-large n...
  const auto is = getSigmaIndex(n, kappa);
  if (is >= m_lambda_kappa.size()) {
    m_lambda_kappa.resize(is + 1, 1.0);
  }
  m_lambda_kappa[is] = lambda;
}

//==============================================================================
std::vector<DiracSpinor>
CorrelationPotential::copy_holes(const std::vector<DiracSpinor> &basis,
                                 const std::vector<DiracSpinor> &core,
                                 int n_min) const {
  std::vector<DiracSpinor> t_holes;
  for (const auto &Fi : basis) {
    const bool inCore = std::find(cbegin(core), cend(core), Fi) != cend(core);
    if (inCore && Fi.n() >= n_min) {
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

//==============================================================================
void CorrelationPotential::addto_G(GMatrix *Gmat, const DiracSpinor &ket,
                                   const DiracSpinor &bra,
                                   const double f) const {
  // XXX kill
  Gmat->add(ket, bra, f);
}

//==============================================================================
DiracSpinor CorrelationPotential::act_G_Fv(const GMatrix &Gmat,
                                           const DiracSpinor &Fv) const {

  // lambda is fitting factor (just scales Sigma|v>)
  // Sigma|v> = int G(r1,r2)*v(r2) dr2
  // (S|v>)_i = sum_j G_ij v_j drdu_j du
  // nb: G is on sub-grid, |v> and S|v> on full-grid. Use interpolation

  // XXX Update!
  return Gmat.drj() * Fv;

  // const auto &gr = *(Fv.grid_sptr());
  // auto SigmaFv = DiracSpinor(0, Fv.kappa(), Fv.grid_sptr());
  // std::vector<double> f(m_subgrid_r.size());
  // std::vector<double> g;
  //
  // for (auto i = 0ul; i < m_subgrid_points; ++i) {
  //   for (auto j = 0ul; j < m_subgrid_points; ++j) {
  //     const auto sj = ri_subToFull(j);
  //     const auto dr = gr.drdu()[sj] * gr.du() * double(m_stride);
  //     f[i] += Gmat.ff[i][j] * Fv.f(sj) * dr;
  //   }
  // }
  //
  // if (m_include_G) {
  //   g.resize(m_subgrid_r.size());
  //   for (auto i = 0ul; i < m_subgrid_points; ++i) {
  //     for (auto j = 0ul; j < m_subgrid_points; ++j) {
  //       const auto sj = ri_subToFull(j);
  //       const auto dr = gr.drdu()[sj] * gr.du() * double(m_stride);
  //       f[i] += Gmat.fg[i][j] * Fv.g(sj) * dr;
  //       g[i] += Gmat.gf[i][j] * Fv.f(sj) * dr;
  //       g[i] += Gmat.gg[i][j] * Fv.g(sj) * dr;
  //     }
  //   }
  // }
  //
  // // Interpolate from sub-grid to full grid
  // SigmaFv.f() = Interpolator::interpolate(m_subgrid_r, f, gr.r());
  // if (m_include_G) {
  //   SigmaFv.g() = Interpolator::interpolate(m_subgrid_r, g, gr.r());
  // }
  // return SigmaFv;
}

//==============================================================================
double CorrelationPotential::act_G_Fv_2(const DiracSpinor &Fa,
                                        const GMatrix &Gmat,
                                        const DiracSpinor &Fb) const {
  // Dores not include Jacobian (assumed already in Gmat)
  return Fa * (Gmat * Fb);
  //
  // // Dores not include Jacobian (assumed already in Gmat)
  //
  // auto aGb = 0.0;
  // for (auto i = 0ul; i < m_subgrid_points; ++i) {
  //   const auto si = ri_subToFull(i);
  //   for (auto j = 0ul; j < m_subgrid_points; ++j) {
  //     const auto sj = ri_subToFull(j);
  //     aGb += Fa.f(si) * Gmat.ff[i][j] * Fb.f(sj);
  //   }
  // }
  //
  // if (m_include_G) {
  //   for (auto i = 0ul; i < m_subgrid_points; ++i) {
  //     const auto si = ri_subToFull(i);
  //     for (auto j = 0ul; j < m_subgrid_points; ++j) {
  //       const auto sj = ri_subToFull(j);
  //       aGb += Fa.f(si) * Gmat.fg[i][j] * Fb.g(sj);
  //       aGb += Fa.g(si) * Gmat.gf[i][j] * Fb.f(sj);
  //       aGb += Fa.g(si) * Gmat.gg[i][j] * Fb.g(sj);
  //     }
  //   }
  // }
  //
  // return aGb;
}

//==============================================================================
double CorrelationPotential::Sigma_vw(const DiracSpinor &v,
                                      const DiracSpinor &w, int max_l) const {
  // Calculates <Fv|Sigma|Fw> from scratch, at Fw energy [full grid + fg+gg]
  if (v.kappa() != w.kappa())
    return 0.0;

  const auto &Ck = m_yeh.Ck();

  // if v.kappa > basis, then Ck angular factor won't exist!
  // Don't extend, since this is a const function
  if (v.twoj() > Ck.max_tj()) {
    std::cout << "\nError: J too large for valence state " << v.symbol()
              << "\n";
    return 0.0;
  }

  if (max_l < 0)
    max_l = 99;

  std::vector<double> delta_a(m_holes.size());
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_holes.size(); ia++) {
    const auto &a = m_holes[ia];
    if (a.l() > max_l)
      continue;
    auto &del_a = delta_a[ia];
    for (const auto &n : m_excited) {
      if (n.l() > max_l)
        continue;
      const auto [kmin_nb, kmax_nb] = Coulomb::k_minmax_Ck(n, a);
      const auto max_k = std::min(m_maxk, kmax_nb);
      for (int k = kmin_nb; k <= max_k; ++k) {
        if (Ck(k, a.kappa(), n.kappa()) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * v.twojp1();

        // Diagrams (a) [direct] and (b) [exchange]
        for (const auto &m : m_excited) {
          if (m.l() > max_l)
            continue;
          const auto Qkv = m_yeh.Q(k, v, a, m, n);
          if (Qkv == 0.0)
            continue;
          const auto Qkw = (&v == &w) ? Qkv : m_yeh.Q(k, w, a, m, n);

          const auto Pkw = m_yeh.P(k, w, a, m, n);
          const auto dele = v.en() + a.en() - m.en() - n.en();
          del_a += ((1.0 / dele / f_kkjj) * (Qkw + Pkw)) * Qkv;
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_holes) {
          if (b.l() > max_l)
            continue;
          const auto Qkv = m_yeh.Q(k, v, n, b, a);
          if (Qkv == 0.0)
            continue;
          const auto Qkw = (&v == &w) ? Qkv : m_yeh.Q(k, w, n, b, a);
          const auto Pkw = m_yeh.P(k, w, n, b, a);
          const auto dele = v.en() + n.en() - b.en() - a.en();
          del_a += ((1.0 / dele / f_kkjj) * (Qkw + Pkw)) * Qkv;
        } // b

      } // k
    }   // n
  }     // a

  return std::accumulate(delta_a.cbegin(), delta_a.cend(), 0.0);
}

//==============================================================================
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
    double r0 = rw == IO::FRW::write ? p_gr->r0() : 0;
    double rmax = rw == IO::FRW::write ? p_gr->rmax() : 0;
    double b = rw == IO::FRW::write ? p_gr->loglin_b() : 0;
    std::size_t pts = rw == IO::FRW::write ? p_gr->num_points() : 0;
    rw_binary(iofs, rw, r0, rmax, b, pts);
    if (rw == IO::FRW::read) {
      const bool grid_ok = std::abs((r0 - p_gr->r0()) / r0) < 1.0e-6 &&
                           std::abs(rmax - p_gr->rmax()) < 0.001 &&
                           std::abs(b - p_gr->loglin_b()) < 0.001 &&
                           pts == p_gr->num_points();
      if (!grid_ok) {
        std::cout << "\nCannot read from:" << fname << ". Grid mismatch\n"
                  << "Read: " << r0 << ", " << rmax << " w/ N=" << pts
                  << ", b=" << b << ",\n but expected: " << p_gr->r0() << ", "
                  << p_gr->rmax() << " w/ N=" << p_gr->num_points()
                  << ", b=" << p_gr->loglin_b() << "\n";
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
    // m_Sigma_kappa.resize(num_kappas, {m_subgrid_points, m_include_G});
    m_Sigma_kappa.resize(
        num_kappas, {m_imin, m_stride, m_subgrid_points, m_include_G, p_gr});
  }

  // Check if include FG/GG written. Note: doesn't matter if mis-match?!
  auto incl_g = rw == IO::FRW::write ? m_include_G : 0;
  rw_binary(iofs, rw, incl_g);

  if (rw == IO::FRW::read) {
    m_nk.resize(num_kappas);
  }
  for (auto &[n, k, en] : m_nk) {
    rw_binary(iofs, rw, n, k, en);
  }

  // Read/Write G matrices
  for (auto &Gk : m_Sigma_kappa) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        rw_binary(iofs, rw, Gk.ff(i, j));
        if (incl_g) {
          rw_binary(iofs, rw, Gk.fg(i, j));
          rw_binary(iofs, rw, Gk.gf(i, j));
          rw_binary(iofs, rw, Gk.gg(i, j));
        }
      }
    }
  }
  std::cout << "done.\n";
  if (rw == IO::FRW::read) {
    std::cout << "Sigma basis: " << basis_config << "\n";
    print_info();
  }
  return true;
}

//==============================================================================
void CorrelationPotential::print_scaling() const {
  if (!m_lambda_kappa.empty()) {
    std::cout << "Scaled Sigma, with: lambda_kappa = ";
    for (const auto &l : m_lambda_kappa) {
      std::cout << l << ", ";
    }
    std::cout << "\n";
  }
}

//==============================================================================
void CorrelationPotential::print_subGrid() const {
  if (m_subgrid_r.empty()) {
    std::cout << "No sub-grid??\n";
    return;
  }
  printf("Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, "
         "stride=%i]\n",
         p_gr->r(m_imin), p_gr->r(m_imin + m_stride * (m_subgrid_points - 1)),
         int(m_subgrid_points), int(m_imin), int(m_stride));
}

//******************************************************************************
//******************************************************************************
// static
DiracSpinor CorrelationPotential::Sigmal_Fv(
    const DiracSpinor &v, const Coulomb::YkTable &yk,
    const Coulomb::LkTable &lk, const std::vector<DiracSpinor> &core,
    const std::vector<DiracSpinor> &excited,
    const std::vector<double> & /*vfk*/, const std::vector<double> & /*veta*/) {

  DiracSpinor SlFv{v.n(), v.kappa(), v.grid_sptr()};
  std::vector<DiracSpinor> SlFv_s(std::size_t(omp_get_max_threads()),
                                  {v.n(), v.kappa(), v.grid_sptr()});

  // auto get_fk = [&vfk](int k) {
  //   if (k < int(vfk.size())) {
  //     return vfk[std::size_t(k)];
  //   }
  //   return 1.0;
  // };
  // auto get_eta = [&veta](int k) {
  //   if (k < int(veta.size())) {
  //     return veta[std::size_t(k)];
  //   }
  //   return 1.0;
  // };

// reduction(+ : SlFv)
#pragma omp parallel for
  for (auto in = 0ul; in < excited.size(); ++in) {
    const auto &n = excited[in];
    const auto ithr = std::size_t(omp_get_thread_num());
    auto &SlFv_ithr = SlFv_s[ithr];
    for (auto &a : core) {

      for (auto &m : excited) {
        const auto [k0, kI] = Coulomb::k_minmax_Q(v, a, m, n);
        for (int k = k0; k <= kI; k += 2) {
          const auto de = v.en() + a.en() - m.en() - n.en();
          const auto tkp1 = (2 * k + 1);
          // #pragma omp critical
          SlFv_ithr += (lk.W(k, m, n, v, a) / de / tkp1) *
                       yk.Qkv_bcd(k, v.kappa(), a, m, n);
          // const auto fk = get_fk(k);
          // const auto etak = get_eta(k);
          // SlFv_ithr += (lk.Q(k, m, n, v, a) * etak / de / tkp1) *
          //              yk.Qkv_bcd(v.k, a, m, n, k);
          // SlFv_ithr += (lk.P(k, m, n, v, a) * fk / de / tkp1) *
          //              yk.Qkv_bcd(v.k, a, m, n, k);
        }
      }

      for (auto &b : core) {
        const auto [k0, kI] = Coulomb::k_minmax_Q(v, n, a, b);
        for (int k = k0; k <= kI; k += 2) {
          const auto de = v.en() + n.en() - a.en() - b.en();
          const auto tkp1 = (2 * k + 1);
          // #pragma omp critical
          SlFv_ithr += (lk.W(k, v, n, a, b) / de / tkp1) *
                       yk.Qkv_bcd(k, v.kappa(), n, a, b);
        }
      }

      //
    }
  }

  SlFv = std::accumulate(SlFv_s.begin(), SlFv_s.end(), SlFv);
  const auto tjp1 = v.twoj() + 1;
  SlFv *= (1.0 / tjp1);

  // // re-scale to account for difference between Fv(BO) and Fv(HF)
  // const auto &v0 = *std::find(excited.begin(), excited.end(), v);
  // for (auto i = 0ul; i < v.rgrid->num_points(); ++i) {
  //   const auto fac_f = std::abs(v0.f(i)) > 1.0e-3 ? v.f(i) / v0.f(i) : 0.0;
  //   const auto fac_g = std::abs(v0.g(i)) > 1.0e-5 ? v.g(i) / v0.g(i) : 0.0;
  //   SlFv.set_f(i) *= fac_f;
  //   SlFv.set_g(i) *= fac_g;
  // }

  return SlFv;
}

//==============================================================================
RDMatrix<double>
CorrelationPotential::Sigma_l(const DiracSpinor &v, const Coulomb::YkTable &yk,
                              const Coulomb::LkTable &lk,
                              const std::vector<DiracSpinor> &core,
                              const std::vector<DiracSpinor> &excited) const {

  RDMatrix<double> Sigma{m_imin, m_stride, m_subgrid_points, m_include_G, p_gr};

  // XXX nb: not 100% tested...

  // XXX Note: careful w/ YkTable - should have same basis as Qk?

  std::vector<RDMatrix<double>> Sigma_n(excited.size(), Sigma);
#pragma omp parallel for
  for (auto in = 0ul; in < excited.size(); ++in) {
    const auto &n = excited[in];
    for (const auto &a : core) {
      // Diagrams (a) + (b)
      for (const auto &m : excited) {

        const auto inv_de = 1.0 / (v.en() + a.en() - m.en() - n.en());
        const auto [k0, kI] = Coulomb::k_minmax_Q(v, a, m, n);
        for (int k = k0; k <= kI; k += 2) {

          // Effective screening parameter:
          const auto fk = 1.0;   // get_fk(k);
          const auto etak = 1.0; // get_eta(k);

          // form Lkv_amn
          const auto Qkv_amn = yk.Qkv_bcd(k, v.kappa(), a, m, n);
          const auto Q_kvamn = yk.Q(k, v, a, m, n);
          const auto L_kmnva = lk.Q(k, m, n, v, a);
          const auto ratio1 = Q_kvamn == 0.0 ? 0.0 : (L_kmnva / Q_kvamn);
          const auto Lkv_amn = ratio1 * Qkv_amn;

          // form Lambdakv_amn
          DiracSpinor Lambdakv_amn{0, v.kappa(), v.grid_sptr()};
          const auto [l0, lI] = Coulomb::k_minmax_Q(m, n, a, v);
          for (int l = l0; l <= lI; l += 2) {
            const auto Lambda_kmnav = lk.Q(l, m, n, a, v);
            const auto Q_kvanm = yk.Q(l, v, a, n, m);
            const auto ratio = Q_kvanm == 0.0 ? 0.0 : Lambda_kmnav / Q_kvanm;
            const auto sj = m_6j.get(m, v, k, n, a, l);
            Lambdakv_amn += sj * ratio * yk.Qkv_bcd(l, v.kappa(), a, n, m); //?
          }
          Lambdakv_amn *= (2 * k + 1);

          // Only include fk on Qk, not Lk integrals - already included there
          // const auto ratio = Omega_kvamn / W_kvamn;
          const auto f = inv_de / (2 * k + 1);
          // XXX Pretty sure fk should only appear in (b), not (a)
          // XXX But maybe etak should not appear either..
          // Sigma.add(Qkv_amn, Lkv_amn, /*fk **/ etak * f); //(a) //fk both?
          // Sigma.add(Qkv_amn, Lambdakv_amn, fk * f);       //(b)
          Sigma_n[in].add(Qkv_amn, Lkv_amn,
                          /*fk **/ etak * f);             //(a) //fk both?
          Sigma_n[in].add(Qkv_amn, Lambdakv_amn, fk * f); //(b)

        } // k
      }   // m

      // Diagrams (c) + (d)
      for (const auto &b : core) {
        const auto inv_de = 1.0 / (v.en() + n.en() - a.en() - b.en());
        const auto [k0, kI] = Coulomb::k_minmax_Q(v, n, a, b);
        for (int k = k0; k <= kI; k += 2) {

          // Effective screening parameter:
          const auto fk = 1.0;   // get_fk(k);
          const auto etak = 1.0; // get_eta(k);

          // form Lkv_nab
          const auto Qkv_nab = yk.Qkv_bcd(k, v.kappa(), n, a, b);
          const auto Q_kvnab = yk.Q(k, v, n, a, b);
          const auto L_kvnab = lk.Q(k, v, n, a, b);
          const auto ratio1 = Q_kvnab == 0.0 ? 0.0 : (L_kvnab / Q_kvnab);
          const auto Lkv_nab = ratio1 * Qkv_nab;

          // form Lambdakv_amn
          DiracSpinor Lambdakv_nab{0, v.kappa(), v.grid_sptr()};
          const auto [l0, lI] = Coulomb::k_minmax_Q(v, n, b, a);
          for (int l = l0; l <= lI; l += 2) {
            const auto Lambda_kvnba = lk.Q(l, v, n, b, a);
            if (Lambda_kvnba == 0.0)
              continue;
            const auto Q_kvnba = yk.Q(l, v, n, b, a);
            const auto ratio = Q_kvnab == 0.0 ? 0.0 : Lambda_kvnba / Q_kvnba;
            const auto sj = m_6j.get(v, a, k, n, b, l);
            Lambdakv_nab += sj * ratio * yk.Qkv_bcd(l, v.kappa(), n, b, a);
          }
          Lambdakv_nab *= (2 * k + 1);

          const auto f = inv_de / (2 * k + 1);
          // XXX Pretty sure fk should only appear in (d), not (c)
          // XXX But maybe etak should not appear either..
          Sigma_n[in].add(Qkv_nab, Lkv_nab, /*fk **/ etak * f); //(c)
          Sigma_n[in].add(Qkv_nab, Lambdakv_nab, fk * f);       //(d)
          //

        } // k
      }   // b

      //
    } // a
  }   // n

  Sigma = std::accumulate(Sigma_n.cbegin(), Sigma_n.cend(), Sigma);
  return (1.0 / v.twojp1()) * Sigma;
}

} // namespace MBPT
