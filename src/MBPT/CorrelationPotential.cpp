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
CorrelationPotential::CorrelationPotential(
    const Grid &gr, const std::vector<DiracSpinor> &core,
    const std::vector<DiracSpinor> &excited, const int in_stride,
    const std::vector<double> &en_list, const std::string &in_fname,
    const HF::HartreeFock *const in_hf)
    : p_gr(&gr),
      m_core(core),
      m_excited(excited),
      m_yec(&gr, &m_excited, &m_core),
      m_maxk(find_max_tj(core, excited)),
      m_6j(m_maxk, m_maxk),
      stride(std::size_t(in_stride)),
      p_hf(in_hf) {
  auto sp = IO::Profile::safeProfiler(__func__);

  std::cout << "\nCorrelation potential (Sigma)\n";

  std::cout << "(Including FF";
  if (include_G)
    std::cout << ", FG/GF, and GG";
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
  // Work out from 'r0eps' for min/max (or, overwrite with doubles?)
  const double rmin = 1.0e-4;
  const double rmax = 30.0;
  // const double rmin = 1.0e-6;
  // const double rmax = 100.0;

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
GMatrix CorrelationPotential::G_single(const DiracSpinor &ket,
                                       const DiracSpinor &bra,
                                       const double f) const {
  GMatrix Gmat(stride_points, include_G);
  addto_G(&Gmat, ket, bra, f);
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
    const auto si = onto_fullGrid(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = onto_fullGrid(j);
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
      const auto sj = onto_fullGrid(j);
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

  Sigma_kappa.resize(en_list.size(), {stride_points, include_G});

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
    GMatrix G_a(stride_points, include_G);
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
GMatrix CorrelationPotential::Green_core(int kappa, double en) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // G_core = \sum_a |a><a|/(e-ea), for all a with a.k=k
  GMatrix Gcore(stride_points, include_G);

  // loop over HF core, not Sigma core (used in subtraction to get G^excited)
  const auto &core = p_hf->get_core();
  for (const auto &a : core) {
    if (a.k == kappa)
      addto_G(&Gcore, a, a, 1.0 / (en - a.en));
  }
  return Gcore;
}
//------------------------------------------------------------------------------
ComplexGMatrix CorrelationPotential::Green_core(int kappa, double en_re,
                                                double en_im) const {
  auto sp = IO::Profile::safeProfiler(__func__, "complex");
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.k=k
  ComplexGMatrix Gcore(stride_points, include_G);

  // loop over HF core, not Sigma core (used in subtraction to get G^excited)
  const auto &core = p_hf->get_core();
  for (const auto &a : core) {
    if (a.k != kappa)
      continue;
    // factor = 1/(a+ib) = (a-ib)/(a^2+b^2), a= enre-e_a, b=enim
    const auto de_re = en_re - a.en;
    const auto a2pb2 = (de_re * de_re + en_im * en_im);
    const ComplexDouble factor{de_re / a2pb2, -en_im / a2pb2};
    Gcore += G_single(a, a, 1.0).make_complex(factor); // expensive!
  }
  return Gcore;
}

//******************************************************************************
GMatrix CorrelationPotential::Green_hf(int kappa, double en) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  DiracSpinor x0(0, kappa, *p_gr);
  DiracSpinor xI(0, kappa, *p_gr);

  const auto &Hmag = p_hf->get_Hrad_mag(x0.l());
  const auto alpha = p_hf->m_alpha;

  auto vl = p_hf->get_vlocal(x0.l());

  // // // This seems to not help...
  // const auto Nc = p_hf->num_core_electrons();
  // const auto eta = -0.5 / Nc;
  // auto eta_vd = p_hf->get_vdir();
  // NumCalc::scaleVec(eta_vd, eta);
  // NumCalc::add_to_vector(vl, eta_vd);

  DiracODE::regularAtOrigin(x0, en, vl, Hmag, alpha);
  DiracODE::regularAtInfinity(xI, en, vl, Hmag, alpha);

  const auto pp = std::size_t(0.65 * double(xI.pinf));
  const auto w = -1.0 * (xI.f[pp] * x0.g[pp] - x0.f[pp] * xI.g[pp]) / alpha;
  // not sure why -ve sign.. is sign even defined by DiracODE?

  const auto g0 = MakeGreensG(x0, xI, w);

  GMatrix Ident(stride_points, include_G);
  Ident.make_identity();
  // const auto Vx = Make_Vx(kappa, eta_vd);
  const auto Vx = Make_Vx(kappa, {});

  return ((Ident - g0 * Vx).invert()) * g0;
}

//******************************************************************************
ComplexGMatrix CorrelationPotential::Polarisation(int k_a, int k_alpha,
                                                  double om_re,
                                                  double om_im) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // XXX This is very inefficient! Can be much improved!
  ComplexGMatrix pi(stride_points, include_G);
  for (const auto &a : m_core) {
    if (a.k != k_a)
      continue;

    auto gre_tmp = Green_hf(k_alpha, a.en - om_re); // real part
    auto g_ex = ComplexG(gre_tmp, -om_im);          // complex g (all states)
    // Now include +w term [re-use gre_tmp allocation]
    gre_tmp = Green_hf(k_alpha, a.en + om_re); // real part
    g_ex += ComplexG(gre_tmp, om_im);          // Full Complex [G(e-w)+G(e+w)]
    // subtract core part: gives [G^ex(e-w) + G^ex(e+w)]
    g_ex -= Green_core(k_alpha, a.en - om_re, -om_im);
    g_ex -= Green_core(k_alpha, a.en + om_re, om_im);

    // ketbra_a = |a><a|
    const auto ketbra_a = G_single(a, a, 1.0).make_complex({0.0, 1.0});
    pi += g_ex.mult_elements_by(ketbra_a);
  }
  return pi;
}

//******************************************************************************
ComplexGMatrix CorrelationPotential::ComplexG(const GMatrix &Gr,
                                              double om_imag) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Given G(wr) and wi, returns G(wr+i*wi)

  // // g = g / [1 + iwI*g0] = 1/[g^-1 + iwI]
  // // G(w) = [g^-1 + iwI]^-1 = G(re(w)+im(w))
  // ComplexGMatrix cG = Gre.make_complex({1.0, 0.0});
  // (cG.invert().plusIdent(0.0, om_imag)).invert();
  // return cG;
  // note: simpler method doesn't work, g^-1 doesn't exist?

  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = Gr * [1 + i*wi*Gr]^-1 = Gr * iX

  auto iX = Gr.make_complex({0.0, om_imag})
                .mult_elements_by(*m_drj) // XXX ???
                .plusIdent(1.0)
                .invert();
  // include dri?

  return iX * Gr.make_complex({1.0, 0.0});

  // Doesn't work. G inverse contains NaN? But GSL doesn' complain?
  // return Gr.make_complex({1.0, 0.0}).invert().plusIdent(0.0,
  // om_imag).invert();
}

//******************************************************************************
GMatrix CorrelationPotential::MakeGreensG(const DiracSpinor &x0,
                                          const DiracSpinor &xI,
                                          const double w) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  // G(r1,r2) = x0(rmin)*xI(imax)/w
  GMatrix g0I(stride_points, include_G);
  // XXX Take advantage of symmetry!?
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = onto_fullGrid(i);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = onto_fullGrid(j);
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

//******************************************************************************
void CorrelationPotential::fill_qhat() {

  std::vector<LinAlg::SqMatrix> qhat(std::size_t(m_maxk + 1), stride_points);

  LinAlg::SqMatrix t_dri{stride_points};

  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = onto_fullGrid(i);
    const auto dri = p_gr->drdu[si] * p_gr->du * double(stride);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = onto_fullGrid(j);
      const auto drj = p_gr->drdu[sj] * p_gr->du * double(stride);
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

  m_qhat.resize(std::size_t(m_maxk) + 1, {stride_points, include_G});
  for (auto k = 0ul; k <= std::size_t(m_maxk); ++k) {
    m_qhat[k].ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, qhat[k]);
    if (include_G) {
      m_qhat[k].gg = m_qhat[k].ff;
    }
  }

  m_dri = new ComplexGMatrix{stride_points, include_G}; // XXX TEMP! XXX
  m_drj = new ComplexGMatrix{stride_points, include_G}; // XXX TEMP! XXX
  m_dri->ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, t_dri);
  m_drj->ff = m_dri->ff.transpose();
}

//******************************************************************************
GMatrix CorrelationPotential::Make_Vx(int kappa,
                                      const std::vector<double> vx0) const {
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
  //     const auto si = onto_fullGrid(i);
  //     for (auto j = 0ul; j < stride_points; j++) {
  //       const auto sj = onto_fullGrid(j);
  //       vxmat += a.f[si] * Vx.ff[i][j] * a.f[sj]; // * double(stride);
  //     }
  //   }
  //   auto vxhf = a * (p_hf->calc_vexFa_core(a));
  //   std::cout << a.symbol() << " " << vxmat << " " << vxhf << "\n";
  // }
  // std::cin.get();

  // subtract off "effective/approx" exchange term (eta*v_dir)
  if (!vx0.empty()) {
    for (auto i = 0ul; i < stride_points; ++i) {
      const auto si = onto_fullGrid(i);
      const auto dri = p_gr->drdu[si] * p_gr->du * double(stride);
      // subtract off "effective/approx" exchange term (eta*v_dir)
      // Think there should only be single dri here (already one integral??)
      Vx.ff[i][i] -= vx0[si] * dri;
      if constexpr (include_G) {
        Vx.fg[i][i] -= vx0[si] * dri;
        Vx.gf[i][i] -= vx0[si] * dri;
        Vx.gg[i][i] -= vx0[si] * dri;
      }
    }
  }

  return Vx;
}

//******************************************************************************
void CorrelationPotential::FeynmanDirect(int kv) {
  auto sp = IO::Profile::safeProfiler(__func__);

  fill_qhat();
  Sigma_kappa[0].zero();

  const double env = -0.127368;
  const double omre = -0.3;

  // Set up imaginary frequency grid:
  std::vector<double> v_omim;
  const auto w0 = 0.01;
  const auto wmax = 100.0;
  const auto w_ratio = 2.0;
  for (auto w = w0; w < wmax * w_ratio; w *= w_ratio) {
    v_omim.push_back(w);
  }
  std::cout << "Im(omega) grid: " << v_omim.front() << " - " << v_omim.back()
            << " in " << v_omim.size() << " steps\n";
  const auto dw_factor = 0.5 * (w_ratio - 1.0 / w_ratio);
  // not sure if endpoint helps?
  const auto dw_endpoint = w0 * 0.25 * (1.0 + 1.0 / w_ratio);

  // {
  //   // Test integration grid:
  //   auto f = [](double x) {
  //     return x * std::exp(-(x - 1.0) * (x - 1.0) / 100.0);
  //   };
  //   double intF = 0.0;
  //   for (const auto omim : v_omim) {
  //     const auto dw =
  //         (omim == w0) ? omim * (dw_factor + dw_endpoint) : omim * dw_factor;
  //     intF += f(omim) * dw;
  //   }
  //   auto exact = NumCalc::num_integrate(f, 0.0, 100.0, 1000);
  //   std::cout << intF << ". error: " << intF - exact << " = "
  //             << 100 * (intF - exact) / exact << "%\n";
  //   std::cin.get();
  // }

  // Greens function (at om_re) remains same, so calculate it once only:
  std::vector<GMatrix> gBetas;
  gBetas.reserve(std::size_t(m_maxkindex + 1));
  for (int ibeta = 0; ibeta <= m_maxkindex; ++ibeta) {
    const auto kbeta = Angular::kappaFromIndex(ibeta);
    gBetas.push_back(Green_hf(kbeta, env + omre));
  }

  for (int ia = 0; ia <= m_maxkindex_core; ++ia) {
    printf("Sigma F: %3i/%3i\r", ia, m_maxkindex_core);
    std::cout << std::flush;
    const auto ka = Angular::kappaFromIndex(ia);
#pragma omp parallel for
    for (int ialpha = 0; ialpha <= m_maxkindex; ++ialpha) {
      const auto kalpha = Angular::kappaFromIndex(ialpha);

      for (const auto omim : v_omim) { // for omega integral
        // nb: do w loop outside beta, since Pi depends only on a,A, and w

        // From middle: (not including x2 from -ve w)
        const auto dw =
            (omim == w0) ? omim * dw_factor + dw_endpoint : omim * dw_factor;

        const auto pi_aalpha = Polarisation(ka, kalpha, omre, omim);

        for (int ibeta = 0; ibeta <= m_maxkindex; ++ibeta) {
          const auto kbeta = Angular::kappaFromIndex(ibeta);

          const auto &g_beta_re = gBetas[std::size_t(ibeta)];
          const auto g_beta = ComplexG(g_beta_re, omim);

          // // sum over k:
          const auto sum_qpq = sumk_cQPQ(kv, ka, kalpha, kbeta, pi_aalpha);
          // // nb: imag part is huge!? why!
          const auto dSdw = mult_elements(g_beta, sum_qpq).get_real();
          // const auto dSdw = (g_beta * sum_qpq).get_real();
#pragma omp critical(sumSig)
          { Sigma_kappa[0] += (2.0 * dw / (2.0 * M_PI)) * dSdw; } // 2.0 for -ve

        } // beta
      }   // w
    }     // alpha
  }       // a
  std::cout << "\n";

  // devide through by dri, drj [these included in q's, but want differential
  // operator for sigma]
  // or.. include one of these in definition of opertion S|v> ?
  for (auto i = 0ul; i < stride_points; ++i) {
    const auto si = onto_fullGrid(i);
    const auto dri = p_gr->drdu[si] * p_gr->du * double(stride);
    for (auto j = 0ul; j < stride_points; ++j) {
      const auto sj = onto_fullGrid(j);
      const auto drj = p_gr->drdu[sj] * p_gr->du * double(stride);
      Sigma_kappa[0].ff[i][j] /= (dri * drj);
      if (include_G) {
        Sigma_kappa[0].fg[i][j] /= (dri * drj);
        Sigma_kappa[0].gf[i][j] /= (dri * drj);
        Sigma_kappa[0].gg[i][j] /= (dri * drj);
      }
    }
  }
}

//******************************************************************************
ComplexGMatrix
CorrelationPotential::sumk_cQPQ(int kv, int ka, int kalpha, int kbeta,
                                const ComplexGMatrix &pi_aalpha) const {
  auto sp = IO::Profile::safeProfiler(__func__);
  // sum over k:
  // sum_k [ck qk * pi(w) * qk], ck angular factor

  auto sum_qpq = ComplexGMatrix(stride_points, include_G);
  // min/max k to include in sum [just to save calculating terms=0]
  // based on two Ck angular factors
  const auto [kmin1, kmax1] = Angular::kminmax_Ck(kv, kbeta);  // Ck_vB
  const auto [kmin2, kmax2] = Angular::kminmax_Ck(ka, kalpha); // Ck_aA
  const auto kmin = std::max(kmin1, kmin2);
  const auto kmax = std::min({kmax1, kmax2, m_maxk});
  const auto tjvp1 = Angular::twoj_k(kv) + 1; //[jv]
  for (int k = kmin; k <= kmax; ++k) {
    const auto c1 =
        Angular::Ck_kk(k, kv, kbeta) * Angular::Ck_kk(k, ka, kalpha);
    if (c1 == 0.0)
      continue;
    const double c_ang = c1 * c1 / double(tjvp1 * (2 * k + 1));
    const auto &q = m_qhat[std::size_t(k)];
    sum_qpq += (c_ang) * (q * (pi_aalpha * q));
  } // k
  // nb: imag part is huge!? why!
  return sum_qpq;
}

} // namespace MBPT
