#include "MBPT/CorrelationPotential.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/Interpolator.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include "Maths/NumCalc_quadIntegrate.hpp"
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
  // returns maximum value of 2*j in {core,excited}
  auto maxtj1 = core.empty() ? 0
                             : std::max_element(core.cbegin(), core.cend(),
                                                DiracSpinor::comp_j)
                                   ->twoj();
  auto maxtj2 = excited.empty()
                    ? 0
                    : std::max_element(excited.cbegin(), excited.cend(),
                                       DiracSpinor::comp_j)
                          ->twoj();
  return std::max(maxtj1, maxtj2);
}
static inline int find_max_l(const std::vector<DiracSpinor> &orbs) {
  return orbs.empty()
             ? 0
             : std::max_element(orbs.cbegin(), orbs.cend(), DiracSpinor::comp_l)
                   ->l();
}
//******************************************************************************

namespace MBPT {

//******************************************************************************
//******************************************************************************
CorrelationPotential::CorrelationPotential(
    const HF::HartreeFock *const in_hf,      //
    const std::vector<DiracSpinor> &core,    //
    const std::vector<DiracSpinor> &excited, //
    const Sigma_params &sigp,                //
    const rgrid_params &subgridp,            //
    const std::vector<double> &en_list,      //
    const std::string &in_fname)
    : p_gr(in_hf->p_rgrid),
      m_core(core),
      m_excited(excited),
      m_yec(p_gr, &m_excited, &m_core),
      stride(subgridp.stride),
      method(sigp.method),
      m_omre(sigp.real_omega),
      p_hf(in_hf) {
  auto sp = IO::Profile::safeProfiler(__func__);

  std::cout << "\nCorrelation potential (Sigma)\n";
  setup_subGrid(subgridp.r0, subgridp.rmax);

  m_min_core_n = sigp.min_n_core;
  // Only used for Feynman:
  m_maxkindex_core = 2 * find_max_l(p_hf->get_core());
  // m_maxkindex = std::max(2 * sigp.max_l_excited, m_maxkindex_core);
  m_maxkindex = 2 * sigp.max_l_excited;

  // find max k and 2j:
  auto max_tj = 0;
  if (method == Method::Feynman) {
    auto kmax_c = find_max_tj(p_hf->get_core(), excited);
    auto kmax_ex = 2 * sigp.max_l_excited + 1;
    max_tj = std::max(kmax_c, kmax_ex);
    m_maxk = std::min(max_tj, sigp.max_k);
  } else {
    max_tj = find_max_tj(core, excited);
    m_maxk = std::min(max_tj, sigp.max_k);
  }
  // fill sixj and Ck:
  m_6j.fill(m_maxk, max_tj);

  const auto fname = in_fname == "" ? "" : in_fname + ".Sigma";
  if (fname != "" && IO::FRW::file_exists(fname)) {
    read_write(fname, IO::FRW::read);
    return;
  }

  std::cout << "Form correlation potential: ";
  if (method == Method::Feynman) {
    std::cout << "Feynman method\n";
    // only needed for Feynman
    prep_Feynman();
    m_yec.extend_Ck(m_maxk, max_tj);
  } else {
    std::cout << "Goldstone method\n";
  }
  if (include_G)
    std::cout << "(Including FG/GF and GG)\n";
  std::cout << "\n";

  form_Sigma(en_list, fname);
}

//******************************************************************************
void CorrelationPotential::setup_subGrid(double rmin, double rmax) {
  // Form the "Sigma sub-grid"
  const auto &rvec = p_gr->r;

  imin = 0;
  for (auto i = 0ul; i < rvec.size(); i += stride) {
    auto r = rvec[i];
    if (r < rmin) {
      imin++;
      continue;
    }
    if (r > rmax)
      break;
    r_stride.push_back(r);
  }

  stride_points = r_stride.size();
  printf(
      "Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, stride=%i]\n",
      r_stride.front(), r_stride.back(), int(stride_points), int(imin),
      int(stride));
}

//******************************************************************************
ComplexGMatrix CorrelationPotential::G_single(const DiracSpinor &ket,
                                              const DiracSpinor &bra,
                                              const ComplexDouble f) const {
  ComplexGMatrix Gmat(stride_points, include_G);
  const auto [x, iy] = f;
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = ri_subToFull(j);
      Gmat.ff[i][j] = gsl_complex_rect(x * ket.f[si] * bra.f[sj],
                                       iy * ket.f[si] * bra.f[sj]);
      if constexpr (include_G) {
        Gmat.fg[i][j] = gsl_complex_rect(x * ket.f[si] * bra.g[sj],
                                         iy * ket.f[si] * bra.g[sj]);
        Gmat.gf[i][j] = gsl_complex_rect(x * ket.g[si] * bra.f[sj],
                                         iy * ket.g[si] * bra.f[sj]);
        Gmat.gg[i][j] = gsl_complex_rect(x * ket.g[si] * bra.g[sj],
                                         iy * ket.g[si] * bra.g[sj]);
      }
    } // j
  }   // i
  return Gmat;
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
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = ri_subToFull(j);
      Gmat->ff[i][j] += f * ket.f[si] * bra.f[sj];
      if constexpr (include_G) {
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
  if constexpr (include_G) {
    g.resize(r_stride.size());
  }
  for (auto i = 0ul; i < stride_points; ++i) {
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = ri_subToFull(j);
      const auto dr = gr.drdu[sj] * gr.du * double(stride);
      f[i] += Gmat.ff[i][j] * Fv.f[sj] * dr * lambda;

      if constexpr (include_G) {
        f[i] += Gmat.fg[i][j] * Fv.g[sj] * dr * lambda;
        g[i] += Gmat.gf[i][j] * Fv.f[sj] * dr * lambda;
        g[i] += Gmat.gg[i][j] * Fv.g[sj] * dr * lambda;
      }
    }
  }
  // Interpolate from sub-grid to full grid
  SigmaFv.f = Interpolator::interpolate(r_stride, f, gr.r);
  if constexpr (include_G) {
    SigmaFv.g = Interpolator::interpolate(r_stride, g, gr.r);
  }

  return SigmaFv;
}

//******************************************************************************
void CorrelationPotential::form_Sigma(const std::vector<double> &en_list,
                                      const std::string &fname) {
  auto sp = IO::Profile::safeProfiler(__func__);

  if (en_list.empty())
    return;

  Sigma_kappa.resize(en_list.size(), {stride_points, include_G});

  std::cout << "Forming correlation potential for:\n";
  for (auto ki = 0ul; ki < en_list.size(); ki++) {
    const auto kappa = Angular::kappaFromIndex(int(ki));

    // if v.kappa > basis, then Ck angular factor won't exist!
    auto tj = Angular::twojFromIndex(int(ki));
    if (tj > m_yec.Ck().max_tj())
      continue;

    printf("k=%2i %6s at en=%8.5f..\n", kappa,
           AtomData::kappa_symbol(kappa).c_str(), en_list[ki]);
    std::cout << std::flush;
    if (method == Method::Goldstone)
      fill_Sigma_k_Gold(&Sigma_kappa[ki], kappa, en_list[ki]);
    else
      fill_Sigma_k_Feyn(&Sigma_kappa[ki], kappa, en_list[ki]);
    // find lowest excited state, output <v|S|v> energy shift:
    auto find_kappa = [=](const auto &a) { return a.k == kappa; };
    const auto vk =
        std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
    if (vk != cend(m_excited))
      std::cout << "  de=" << *vk * Sigma2Fv(*vk) << "\n";
    // std::cout << "\n";
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
void CorrelationPotential::fill_Sigma_k_Feyn(GMatrix *Gmat, const int kappa,
                                             const double en) {
  auto sp = IO::Profile::safeProfiler(__func__);

  Gmat->zero();

  auto find_kappa = [=](const auto &f) { return f.k == kappa; };
  const auto Fk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);

  const auto dir = FeynmanDirect(kappa, en);

  if (Fk != cend(m_excited)) {
    std::cout << "  <d>=" << *Fk * Sigma_G_Fv(dir, *Fk)
              << ", <e1>=" << std::flush;
  }

  const auto exch = FeynmanEx_1(kappa, en);

  if (Fk != cend(m_excited)) {
    std::cout << *Fk * Sigma_G_Fv(exch, *Fk) << "." << std::flush;
  }

  *Gmat = dir + exch;
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
  Gmat->zero();
  auto Sdir = *Gmat; // blank copies!
  auto Sexch = *Gmat;

  if (m_core.empty())
    return;
  const auto &gr = *p_gr;

  // auto fk = [&](int k) { return 1.0; };

  // Note: get_yk_ab() must only be called with k for which y^k_ab exists!
  // Therefore, must use the k_minmax() function provided
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_core.size(); ia++) {
    const auto &a = m_core[ia];
    GMatrix Ga_d(stride_points, include_G);
    GMatrix Ga_x(stride_points, include_G);
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
          addto_G(&Ga_d, Qkv, Qkv, factor);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // m

        // Diagrams (c) [direct] and (d) [exchange]
        for (const auto &b : m_core) {
          if (Ck(k, kappa, b.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, n, b, a, k, yknb, Ck);
          Coulomb::Pkv_bcd(&Pkv, n, b, a, k, m_yec(n, b), Ck, m_6j);
          const auto dele = en + n.en - b.en - a.en;
          const auto factor = 1.0 / (f_kkjj * dele);
          addto_G(&Ga_d, Qkv, Qkv, factor);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // b

      } // k
    }   // n
#pragma omp critical(sum_D)
    { Sdir += Ga_d; }
#pragma omp critical(sum_X)
    { Sexch += Ga_x; }
    // #pragma omp critical(sumG)
    //     { *Gmat += G_a; }
  } // a

  *Gmat = Sdir + Sexch;

  auto find_kappa = [=](const auto &f) { return f.k == kappa; };
  const auto Fk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);
  if (Fk != cend(m_excited)) {
    std::cout << "  <d>=" << *Fk * Sigma_G_Fv(Sdir, *Fk) << ", ";
    std::cout << "<x>=" << *Fk * Sigma_G_Fv(Sexch, *Fk) << ".";
  }
}

//******************************************************************************
double CorrelationPotential::SOEnergyShift(const DiracSpinor &v,
                                           const DiracSpinor &w) const {
  // Calculates <Fv|Sigma|Fw> from scratch, at Fw energy [full grid + fg+gg]
  if (v.k != w.k)
    return 0.0;

  const auto &Ck = m_yec.Ck();

  // if v.kappa > basis, then Ck angular factor won't exist!
  // XXX Fix! extend Ck! (and sixj)
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
    Sigma_kappa.resize(num_kappas, {stride_points, include_G});
  }

  // Check if include FG/GG written. Note: doesn't matter if mis-match?!
  auto incl_g = rw == IO::FRW::write ? include_G : 0;
  rw_binary(iofs, rw, incl_g);

  // Read/Write G matrices
  for (auto &Gk : Sigma_kappa) {
    for (auto i = 0ul; i < stride_points; ++i) {
      for (auto j = 0ul; j < stride_points; ++j) {
        rw_binary(iofs, rw, Gk.ff[i][j]);
        if (incl_g) {
          rw_binary(iofs, rw, Gk.fg[i][j]);
          rw_binary(iofs, rw, Gk.gf[i][j]);
          rw_binary(iofs, rw, Gk.gg[i][j]);
        }
      }
    }
  }
  std::cout << "... done.\n";
  printf(
      "Sigma sub-grid: r=(%.1e, %.1f)aB with %i points. [i0=%i, stride=%i]\n",
      r_stride.front(), r_stride.back(), int(stride_points), int(imin),
      int(stride));
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
//******************************************************************************

//******************************************************************************
ComplexGMatrix CorrelationPotential::Green_core(int kappa, double en_re,
                                                double en_im) const {
  auto sp = IO::Profile::safeProfiler(__func__, "complex");
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.k=k
  ComplexGMatrix Gcore(stride_points, include_G);

  // loop over HF core, not Sigma core (used in subtraction to get G^excited)
  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    // for (const auto &a : core) {
    if (a.k != kappa)
      continue;
    // factor = 1/(a+ib) = (a-ib)/(a^2+b^2), a= enre-e_a, b=enim
    const auto de_re = en_re - a.en;
    const auto a2pb2 = (de_re * de_re + en_im * en_im);
    const ComplexDouble factor{de_re / a2pb2, -en_im / a2pb2};
    Gcore += factor * m_Pa[ia];
  }
  return Gcore;
}

//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::Green_ex(int kappa, double en_re,
                                              double en_im) const {
  if constexpr (basis_for_Green) {
    return Green_hf_basis(kappa, en_re, en_im, true);
  }
  return Green_hf(kappa, en_re, en_im) - Green_core(kappa, en_re, en_im);
}

//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::Green_hf_basis(int kappa, double en_re,
                                                    double en_im,
                                                    bool ex_only) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Test Green fn using basis:
  // XXX THIS MAKES A DIFFERENCE!!! WHY!? Only in Polarisation Op.
  ComplexGMatrix Gc(stride_points, include_G);
  const auto &core = p_hf->get_core();
  const auto &ex = m_excited;
  for (const auto orbs : {&core, &ex}) {
    if (ex_only && orbs == &core)
      continue;
    for (const auto &a : *orbs) {
      if (a.k != kappa)
        continue;
      const auto de_re = en_re - a.en;
      const auto a2pb2 = (de_re * de_re + en_im * en_im);
      const ComplexDouble factor{de_re / a2pb2, -en_im / a2pb2};
      Gc += G_single(a, a, factor);
    }
  }
  return Gc;
}

//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::Green_hf(int kappa, double en_re,
                                              double en_im) const {
  auto sp = IO::Profile::safeProfiler(__func__, "complexE");
  return ComplexG(Green_hf(kappa, en_re), en_im);
}

//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::Green_hf(int kappa, double en) const {
  auto sp = IO::Profile::safeProfiler(__func__);

  if constexpr (basis_for_Green) {
    return Green_hf_basis(kappa, en, 0.0, false);
  }

  DiracSpinor x0(0, kappa, *p_gr);
  DiracSpinor xI(0, kappa, *p_gr);

  const auto &Hmag = p_hf->get_Hrad_mag(x0.l());
  const auto alpha = p_hf->get_alpha();
  const auto vl = p_hf->get_vlocal(x0.l());
  DiracODE::regularAtOrigin(x0, en, vl, Hmag, alpha);
  DiracODE::regularAtInfinity(xI, en, vl, Hmag, alpha);

  const auto pp = std::size_t(0.65 * double(xI.pinf));
  const auto w = -1.0 * (xI.f[pp] * x0.g[pp] - x0.f[pp] * xI.g[pp]) / alpha;
  // not sure why -ve sign.. is sign even defined by DiracODE?

  const auto g0 = MakeGreensG(x0, xI, w);
  const auto &Vx = m_Vxk[std::size_t(Angular::indexFromKappa(kappa))];
  // nb: much faster to invert _before_ make complex!
  return ((-1.0 * g0 * Vx).plusIdent(1.0).invert() * g0)
      .make_complex({1.0, 0.0});
}
//------------------------------------------------------------------------------
GMatrix CorrelationPotential::MakeGreensG(const DiracSpinor &x0,
                                          const DiracSpinor &xI,
                                          const double w) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  // G(r1,r2) = x0(rmin)*xI(imax)/w
  GMatrix g0I(stride_points, include_G);
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = ri_subToFull(j);
      const auto irmin = std::min(sj, si);
      const auto irmax = std::max(sj, si);
      g0I.ff[i][j] = x0.f[irmin] * xI.f[irmax] / w;
      if constexpr (include_G) {
        g0I.fg[i][j] = x0.f[irmin] * xI.g[irmax] / w;
        g0I.gf[i][j] = x0.g[irmin] * xI.f[irmax] / w;
        g0I.gg[i][j] = x0.g[irmin] * xI.g[irmax] / w;
      }
    } // j
  }   // i
  return g0I;
}
//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::ComplexG(const ComplexGMatrix &Gr,
                                              double om_imag) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Given G(wr) and wi, returns G(wr+i*wi)
  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = Gr * [1 + i*wi*Gr]^-1 = Gr * iX
  const ComplexDouble iw{0.0, om_imag};
  return ((iw * Gr).mult_elements_by(*m_drj).plusIdent(1.0).invert()) * Gr;
}
//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::ComplexG(const ComplexGMatrix &Gr,
                                              double de_re,
                                              double de_im) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Given G(wr) and wi, returns G(wr+i*wi)
  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = Gr * [1 + i*wi*Gr]^-1 = Gr * iX
  // Note: only works when de_re=0.0 ??
  const ComplexDouble dele{de_re, de_im};
  return ((dele * Gr).mult_elements_by(*m_drj).plusIdent(1.0).invert()) * Gr;
}

//******************************************************************************
ComplexGMatrix CorrelationPotential::Polarisation(int k_a, int k_alpha,
                                                  double om_re,
                                                  double om_im) const {
  auto sp = IO::Profile::safeProfiler(__func__);

  ComplexGMatrix pi(stride_points, include_G);

  const auto &core = p_hf->get_core();
  const auto Iunit = ComplexDouble{0.0, 1.0};

  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.n < m_min_core_n)
      continue;
    if (a.k != k_a)
      continue;
    const auto &Ga = m_Pa[ia];

    // // Basis for real part, invert for imag: Works
    // // Trying to see if we can shift (e_a +/ wr) to 'rhs'.. no?
    // pi += (ComplexG2(Green_hf_basis(k_alpha, +om_re, 0.0, false),
    //                  a.en - 2 * om_re, -om_im) -
    //        Green_core(k_alpha, a.en - om_re, -om_im) +
    //        ComplexG2(Green_hf_basis(k_alpha, om_re, 0.0, false), a.en, om_im)
    //        - Green_core(k_alpha, a.en + om_re, om_im))
    //           .mult_elements_by(Ga);

    if constexpr (basis_for_Pol) {
      // Full basis: Works
      pi += (Green_hf_basis(k_alpha, a.en - om_re, -om_im, true) +
             Green_hf_basis(k_alpha, a.en + om_re, om_im, true))
                .mult_elements_by(Ga);
    } else {
      // No basis: doesn't work ???
      // What if I solve at omre, but then shift by 'de' for each?
      pi += (Green_ex(k_alpha, a.en - om_re, -om_im) +
             Green_ex(k_alpha, a.en + om_re, om_im))
                .mult_elements_by(Ga);
    }
  }
  return Iunit * pi;
}
//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::Polarisation_a(const ComplexGMatrix &pa,
                                                    double ena, int kA,
                                                    double omre,
                                                    double omim) const {
  auto sp = IO::Profile::safeProfiler(__func__);

  const auto Iunit = ComplexDouble{0.0, 1.0};

  if constexpr (basis_for_Pol) {
    // Full basis: Works
    return Iunit * (Green_hf_basis(kA, ena - omre, -omim, true) +
                    Green_hf_basis(kA, ena + omre, omim, true))
                       .mult_elements_by(pa);
  } else {
    // No basis: no work
    return Iunit *
           (Green_ex(kA, ena - omre, -omim) + Green_ex(kA, ena + omre, omim))
               .mult_elements_by(pa);
  }
}
// //------------------------------------------------------------------------------
// void CorrelationPotential::sumPol(const ComplexGMatrix &pi_aA) const {
//   // this is just to check Pol
//   auto rePi = pi_aA.get_real();
//   auto imPi = pi_aA.get_imaginary();
//   double sum_re = 0.0;
//   double sum_im = 0.0;
//   for (auto i = 0ul; i < stride_points; ++i) {
//     const auto dri = dr_subToFull(i);
//     for (auto j = 0ul; j < stride_points; ++j) {
//       const auto drj = dr_subToFull(j);
//       sum_re += rePi.ff[i][j] * dri * drj;
//       sum_im += imPi.ff[i][j] * dri * drj;
//     }
//   }
//   std::cout << sum_re << "+" << sum_im << "i\n";
// }
//******************************************************************************
std::size_t CorrelationPotential::ri_subToFull(std::size_t i) const {
  return ((imin + i) * stride);
}
double CorrelationPotential::dr_subToFull(std::size_t i) const {
  return p_gr->drdu[ri_subToFull(i)] * p_gr->du * double(stride);
}
//******************************************************************************
void CorrelationPotential::prep_Feynman() {

  std::vector<LinAlg::SqMatrix> qhat(std::size_t(m_maxk + 1), stride_points);

  LinAlg::SqMatrix t_dri{stride_points};

  // calculate (real) Qhat and Jacobian
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = ri_subToFull(i);
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = ri_subToFull(j);
      const auto drj = dr_subToFull(j);
      const auto rl = p_gr->r[std::min(si, sj)];
      const auto rm = p_gr->r[std::max(si, sj)];
      const auto ratio = rl / rm; // = r_< / r_>
      // q0 = (1.0 / rm) = r_<^k / r_>^k+1 ,  for k=0
      qhat[0][i][j] = (1.0 / rm) * dri * drj;
      t_dri[i][j] = dri;
      for (auto k = 1ul; k <= std::size_t(m_maxk); ++k) {
        qhat[k][i][j] = qhat[k - 1][i][j] * ratio;
      }
    } // j
  }   // i

  // fill (complex) Qhat
  m_qhat.resize(std::size_t(m_maxk) + 1, {stride_points, include_G});
  for (auto k = 0ul; k <= std::size_t(m_maxk); ++k) {
    m_qhat[k].ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, qhat[k]);
    if (include_G) {
      m_qhat[k].gg = m_qhat[k].ff;
    }
  }

  // fill (complex) Jacobian
  m_dri = std::make_unique<ComplexGMatrix>(stride_points, include_G);
  m_drj = std::make_unique<ComplexGMatrix>(stride_points, include_G);
  m_dri->ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, t_dri);
  m_drj->ff = m_dri->ff.transpose();

  // Fill Vx:
  const auto max_kappa_index = std::size_t(m_maxkindex);
  m_Vxk.resize(std::size_t(max_kappa_index + 1), {stride_points, include_G});
  for (auto ik = 0ul; ik <= max_kappa_index; ik++) {
    m_Vxk[ik] = Make_Vx(Angular::kappaFromIndex(int(ik)));
  }

  // Fill core |a><a|
  const auto &core = p_hf->get_core();
  m_Pa.resize(core.size(), {stride_points, include_G});
  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_Pa[ia] = G_single(core[ia], core[ia], {1.0, 0.0});
  }

  // Set up imaginary frequency grid:
  std::cout << "Re(w) = " << m_omre << "\n";
  {
    const auto w0 = 0.03;
    const auto wmax = 700.0;
    const double wratio = 1.3;
    const std::size_t wsteps = Grid::calc_num_points_from_du(
        w0, wmax, std::log(wratio), GridType::logarithmic);
    // const auto wgrid = Grid(w0, wmax, wsteps, GridType::logarithmic);
    m_wgridD = std::make_unique<Grid>(w0, wmax, wsteps, GridType::logarithmic);
    std::cout << "Dir | Im(w) " << m_wgridD->gridParameters();
    printf(". r=%.2f\n", m_wgridD->r[1] / m_wgridD->r[0]);
  }
  {
    const auto w0 = 0.08;
    const auto wmax = 35.0;
    const double wratio = 2.5;
    const std::size_t wsteps = Grid::calc_num_points_from_du(
        w0, wmax, std::log(wratio), GridType::logarithmic);
    // const auto wgrid = Grid(w0, wmax, wsteps, GridType::logarithmic);
    m_wgridX = std::make_unique<Grid>(w0, wmax, wsteps, GridType::logarithmic);
    std::cout << "Exch| Im(w) " << m_wgridX->gridParameters();
    printf(". r=%.2f\n", m_wgridX->r[1] / m_wgridX->r[0]);
  }

  // // Fill green:
  // const auto tmp_max_kindex = std::max(m_maxkindex, m_maxkindex_core);
  // m_Gr.reserve(std::size_t(tmp_max_kindex + 1));
  // for (int ik = 0; ik <= tmp_max_kindex; ++ik) {
  //   const auto kB = Angular::kappaFromIndex(ik);
  //   gBetas.push_back(Green_hf(kB, env + omre));
  // }
}

//******************************************************************************
GMatrix CorrelationPotential::Make_Vx(int kappa) const {
  auto sp = IO::Profile::safeProfiler(__func__);

  const auto tj = Angular::twoj_k(kappa);
  GMatrix Vx(stride_points, include_G);
  const auto &core = p_hf->get_core();

  const auto kmax = Angular::twojFromIndex(m_maxkindex_core);
  for (int k = 0; k <= kmax; ++k) {
    GMatrix Vx_k(stride_points, include_G);
    for (const auto &a : core) {
      const auto ck = Angular::Ck_kk(k, kappa, a.k); // XXX lookup
      if (ck == 0.0)
        continue;
      const auto c_ang = -1.0 * ck * ck / double(tj + 1);
      addto_G(&Vx_k, a, a, c_ang);
    }
    Vx_k.mult_elements_by(m_qhat[std::size_t(k)].get_real());
    Vx += Vx_k;
  }

  // // Just for testing Vx matrix.
  // for (const auto &a : core) {
  //   if (a.k != kappa)
  //     continue;
  //   double vxmat = 0.0;
  //   for (auto i = 0ul; i < stride_points; i++) {
  //     const auto si = ri_subToFull(i);
  //     for (auto j = 0ul; j < stride_points; j++) {
  //       const auto sj = ri_subToFull(j);
  //       vxmat += a.f[si] * Vx.ff[i][j] * a.f[sj]; // * double(stride);
  //     }
  //   }
  //   auto vxhf = a * (p_hf->calc_vexFa_core(a));
  //   std::cout << a.symbol() << " " << vxmat << " " << vxhf << "\n";
  // }
  // std::cin.get();

  return Vx;
}

//******************************************************************************
GMatrix CorrelationPotential::FeynmanDirect(int kv, double env) {
  auto sp = IO::Profile::safeProfiler(__func__);

  GMatrix Sigma(stride_points, include_G);

  // // Set up imaginary frequency grid:
  const double omre = m_omre; // does seem to depend on this..
  const auto &wgrid = *m_wgridD;

  // Greens function (at om_re) remains same, so calculate it once only:
  std::vector<ComplexGMatrix> gBetas;
  gBetas.reserve(std::size_t(m_maxkindex + 1));
  for (int ibeta = 0; ibeta <= m_maxkindex; ++ibeta) {
    const auto kB = Angular::kappaFromIndex(ibeta);
    gBetas.push_back(Green_hf(kB, env + omre));
  }

  for (int ia = 0; ia <= m_maxkindex_core; ++ia) {
    // printf("Sigma F: %3i/%3i\r", ia + 1, m_maxkindex_core + 1);
    // std::cout << std::flush;
    const auto ka = Angular::kappaFromIndex(ia);
#pragma omp parallel for
    for (int ialpha = 0; ialpha <= m_maxkindex; ++ialpha) {
      const auto kA = Angular::kappaFromIndex(ialpha);

      // nb: do w loop outside beta, since Pi depends only on a,A, and w
      for (auto iw = 0ul; iw < wgrid.num_points; iw++) { // for omega integral
        const auto omim = wgrid.r[iw];
        const auto dw = wgrid.drdu[iw] * wgrid.du;

        const auto pi_aalpha = Polarisation(ka, kA, omre, omim);

        for (int ibeta = 0; ibeta <= m_maxkindex; ++ibeta) {
          const auto kB = Angular::kappaFromIndex(ibeta);

          const auto &g_beta_re = gBetas[std::size_t(ibeta)];
          const auto g_beta = ComplexG(g_beta_re, omim);

          // sum over k:
          const auto gqpq = sumk_cGQPQ(kv, ka, kA, kB, g_beta, pi_aalpha);
#pragma omp critical(sumSig)
          { Sigma += (dw / M_PI) * gqpq; } // 2.0 for -ve

        } // beta
      }   // w
    }     // alpha
  }       // a

  // devide through by dri, drj [these included in q's, but want differential
  // operator for sigma]
  // or.. include one of these in definition of opertion S|v> ?
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sigma.ff[i][j] /= (dri * drj);
      if (include_G) {
        Sigma.fg[i][j] /= (dri * drj);
        Sigma.gf[i][j] /= (dri * drj);
        Sigma.gg[i][j] /= (dri * drj);
      }
    }
  }

  return Sigma;
}

//******************************************************************************
GMatrix CorrelationPotential::FeynmanEx_1(int kv, double env) {
  auto sp = IO::Profile::safeProfiler(__func__);

  // XXX Sort of works? Out by factor? Unstable?

  GMatrix Sx_e1(stride_points, include_G);

  // // Set up imaginary frequency grid:
  const double omre = m_omre; // does seem to depend on this..
  const auto &wgrid = *m_wgridX;

  // Greens function (at om_re) remains same, so calculate it once only:
  std::vector<ComplexGMatrix> gws;
  for (int ik = 0; ik <= m_maxkindex; ++ik) {
    const auto kappa = Angular::kappaFromIndex(ik);
    gws.push_back(Green_hf(kappa, env + omre));
  }

  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    // printf("Sigma F(e1): %3i/%3i\r", int(ia) + 1, (int)core.size());
    // std::cout << std::flush;
    const auto &Fa = core[ia];
    if (Fa.n < m_min_core_n)
      continue;
    const auto ka = Fa.k;
#pragma omp parallel for
    for (int ialpha = 0; ialpha <= m_maxkindex; ++ialpha) {
      const auto kA = Angular::kappaFromIndex(ialpha);
      const auto &gA_re = gws[std::size_t(ialpha)];
      for (int ibeta = 0; ibeta <= m_maxkindex; ++ibeta) {
        const auto kB = Angular::kappaFromIndex(ibeta);

        for (auto iw = 0ul; iw < wgrid.num_points; iw++) {
          const auto omim = wgrid.r[iw];
          const auto dw1 = wgrid.drdu[iw] * wgrid.du;

          const auto gA = ComplexG(gA_re, omim);
          const auto gxBm = Green_ex(kB, Fa.en - omre, -omim);
          const auto gxBp = Green_ex(kB, Fa.en + omre, omim);
          const auto &pa = m_Pa[ia];

          const auto gqpg = sumkl_GQPGQ(gA, gxBm, gxBp, pa, kv, kA, kB, ka);

#pragma omp critical(sumSigX)
          {
            const auto fudge_factor = 1.0 / (2 * M_PI); //???
            Sx_e1 += (fudge_factor * dw1 / M_PI) * gqpg;
          }

        } // w1
      }   // beta
    }     // alpha
  }       // a

  return Sx_e1;
}

//******************************************************************************
GMatrix
CorrelationPotential::sumk_cGQPQ(int kv, int ka, int kA, int kB,
                                 const ComplexGMatrix &g_beta,
                                 const ComplexGMatrix &pi_aalpha) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // sum over k:
  // sum_k [ck qk * pi(w) * qk], ck angular factor

  auto gqpq = ComplexGMatrix(stride_points, include_G);
  // min/max k to include in sum [just to save calculating terms=0]
  // based on two Ck angular factors
  const auto [kmin1, kmax1] = Angular::kminmax_Ck(kv, kB); // Ck_vB
  const auto [kmin2, kmax2] = Angular::kminmax_Ck(ka, kA); // Ck_aA
  const auto kmin = std::max(kmin1, kmin2);
  const auto kmax = std::min({kmax1, kmax2, m_maxk});
  const auto tjvp1 = Angular::twoj_k(kv) + 1; //[jv]
  for (int k = kmin; k <= kmax; ++k) {
    const auto c1 = Angular::Ck_kk(k, kv, kB) * Angular::Ck_kk(k, ka, kA);
    if (c1 == 0.0)
      continue;
    const double c_ang = c1 * c1 / double(tjvp1 * (2 * k + 1));
    const auto &q = m_qhat[std::size_t(k)];
    gqpq += c_ang * q * pi_aalpha * q;
  } // k
  // return gqpq;
  return gqpq.mult_elements_by(g_beta).get_real();
}

//******************************************************************************
GMatrix CorrelationPotential::sumkl_GQPGQ(const ComplexGMatrix &gA,
                                          const ComplexGMatrix &gxBm,
                                          const ComplexGMatrix &gxBp,
                                          const ComplexGMatrix &pa, int kv,
                                          int kA, int kB, int ka) const {
  auto sp = IO::Profile::safeProfiler(__func__);

  auto sum_GQPG = GMatrix(stride_points, include_G);

  const auto [kmin1, kmax1] = Angular::kminmax_Ck(kv, kA); // Ck_vA
  const auto [kmin2, kmax2] = Angular::kminmax_Ck(ka, kB); // Ck_aB
  const auto kmin = std::max(kmin1, kmin2);
  const auto kmax = std::min({kmax1, kmax2, m_maxk});

  const auto [lmin1, lmax1] = Angular::kminmax_Ck(kv, kB); // Ck_vB
  const auto [lmin2, lmax2] = Angular::kminmax_Ck(ka, kA); // Ck_aA
  const auto lmin = std::max(lmin1, lmin2);
  const auto lmax = std::min({lmax1, lmax2, m_maxk});

  const auto tjv = Angular::twoj_k(kv);
  const auto tjA = Angular::twoj_k(kA);
  const auto tjB = Angular::twoj_k(kB);
  const auto tja = Angular::twoj_k(ka);
  const auto tjvp1 = tjv + 1; //[jv]
  const auto s2 = Angular::evenQ_2(tja - tjB) ? 1 : -1;

  const ComplexDouble iI{0.0, 1.0}; // i

  for (auto k = kmin; k <= kmax; ++k) {
    const auto c0 = Angular::Ck_kk(k, kv, kA) * Angular::Ck_kk(k, ka, kB);
    if (Angular::zeroQ(c0))
      continue;

    const auto &qk = m_qhat[std::size_t(k)];

    for (auto l = lmin; l <= lmax; ++l) {
      const auto c1 = Angular::Ck_kk(l, kv, kB) * Angular::Ck_kk(l, ka, kA);
      const auto c2 = Angular::Ck_kk(l, kv, ka) * Angular::Ck_kk(l, kB, kA);
      if (Angular::zeroQ(c1) && Angular::zeroQ(c2))
        continue;

      const auto sj1 =
          c1 == 0.0 ? 0.0 : Angular::sixj_2(tjv, tjA, 2 * k, tja, tjB, 2 * l);
      const auto sj2 =
          Angular::zeroQ(c2)
              ? 0.0
              : tja == tjB ? sj1
                           : Angular::sixj_2(tjv, tjA, 2 * k, tjB, tja, 2 * l);
      if (Angular::zeroQ(sj1) && Angular::zeroQ(sj2))
        continue;

      const auto &ql = m_qhat[std::size_t(l)];

      const auto s0 = Angular::evenQ(k + l) ? 1 : -1;
      const double cang1 = s0 * c0 * c1 * sj1 / tjvp1;
      const double cang2 = s0 * c0 * c2 * s2 * sj2 / tjvp1;

      for (auto r1 = 0ul; r1 < stride_points; ++r1) {
        const auto dr1 =
            p_gr->drdu[ri_subToFull(r1)] * p_gr->du * double(stride);
        for (auto r2 = 0ul; r2 < stride_points; ++r2) {
          const auto dr2 =
              p_gr->drdu[ri_subToFull(r2)] * p_gr->du * double(stride);
          //
          double ij_sum1 = 0.0; // real
          double ij_sum2 = 0.0; // real
          if (!Angular::zeroQ(cang1)) {
            for (auto i = 0ul; i < stride_points; ++i) {
              ComplexDouble jsum{0.0, 0.0}; // complex!
              for (auto j = 0ul; j < stride_points; ++j) {
                jsum += qk.ffc(r1, j) * pa.ffc(i, j) * gxBm.ffc(j, r2);
              }
              ij_sum1 += (iI * jsum * gA.ffc(r1, i) * ql.ffc(i, r2)).re;
            }
          }
          if (!Angular::zeroQ(cang2)) {
            for (auto i = 0ul; i < stride_points; ++i) {
              ComplexDouble jsum{0.0, 0.0}; // complex!
              for (auto j = 0ul; j < stride_points; ++j) {
                jsum += qk.ffc(r1, j) * gxBp.ffc(i, j) * pa.ffc(j, r2);
              }
              ij_sum2 += (iI * jsum * gA.ffc(r1, i) * ql.ffc(i, r2)).re;
            }
          }
          //
          sum_GQPG.ff[r1][r2] +=
              (cang1 * ij_sum1 + cang2 * ij_sum2) / (dr1 * dr2);
        }
      }
      //
    } // l
  }   // k

  return sum_GQPG;
}

} // namespace MBPT
