#include "MBPT/FeynmanSigma.hpp"
#include "Angular/Angular_tables.hpp"
#include "Coulomb/Coulomb.hpp"
#include "Coulomb/YkTable.hpp"
#include "DiracODE/DiracODE.hpp"
#include "HF/HartreeFock.hpp"
#include "IO/FRW_fileReadWrite.hpp"
#include "IO/SafeProfiler.hpp"
#include "MBPT/CorrelationPotential.hpp"
#include "MBPT/GreenMatrix.hpp"
#include "Maths/Grid.hpp"
#include "Maths/LinAlg_MatrixVector.hpp"
#include <algorithm>
#include <cassert>
#include <numeric>

//! Many-body perturbation theory
namespace MBPT {

//******************************************************************************
FeynmanSigma::FeynmanSigma(const HF::HartreeFock *const in_hf,
                           const std::vector<DiracSpinor> &basis,
                           const Sigma_params &sigp,
                           const rgrid_params &subgridp,
                           const std::vector<double> &en_list,
                           const std::string &atom)
    : CorrelationPotential(in_hf, basis, sigp, subgridp),
      m_screen_Coulomb(sigp.screenCoulomb),
      m_omre(sigp.real_omega),
      basis_for_Green(sigp.GreenBasis),
      basis_for_Pol(sigp.PolBasis),
      m_Green_method(sigp.GreenBasis ? GrMethod::basis : GrMethod::Green),
      m_Pol_method(sigp.PolBasis ? GrMethod::basis : GrMethod::Green),
      p_hf(in_hf),
      m_min_core_n(sigp.min_n_core),
      m_max_kappaindex_core(2 * DiracSpinor::max_l(p_hf->get_core())),
      m_max_kappaindex(2 * sigp.max_l_excited) {

  std::cout << "\nCorrelation potential (Sigma): Feynman\n";

  // note: if l_max (for Green's fn) larger than l in basis, need to increase
  // m_kmax, and therefore extend the 6j and Ck tables
  const auto kmax_ex = m_max_kappaindex + 1;
  if (kmax_ex > m_maxk)
    m_maxk = kmax_ex;

  // io file name:
  const std::string ext = ".SigmaF";
  const auto fname =
      atom == "" ? "" : atom + "_" + std::to_string(p_gr->num_points) + ext;
  const bool read_ok = read_write(fname, IO::FRW::read);

  // if (en_list.empty())
  //   return;

  if (!read_ok) {
    // Extand 6j and Ck
    m_6j.fill(m_maxk, m_maxk);
    m_yeh.extend_Ck(m_maxk, m_maxk);
    prep_Feynman();
    form_Sigma(en_list, fname);
  }
}

//******************************************************************************
//******************************************************************************
void FeynmanSigma::fill_Sigma_k(GMatrix *Sigma, const int kappa,
                                const double en) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Find lowest "valence" state from the basis (for energy shift)
  const auto find_kappa = [=](const auto &f) { return f.k == kappa; };
  const auto Fk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);

  // Direct:
  *Sigma = FeynmanDirect(kappa, en);

  // Print out the direct energy shift:
  if (Fk != cend(m_excited)) {
    const auto deD = *Fk * act_G_Fv(*Sigma, *Fk);
    printf("de= %.4f + ", deD);
    std::cout << std::flush;
  }

  const auto exch = m_ex_method == ExchangeMethod::Goldstone
                        ? Exchange_Goldstone(kappa, en)
                        : m_ex_method == ExchangeMethod::w1
                              ? FeynmanEx_1(kappa, en)
                              : FeynmanEx_w1w2(kappa, en);

  // Print out the exchange energy shift:
  if (Fk != cend(m_excited)) {
    const auto deX = *Fk * act_G_Fv(exch, *Fk);
    printf("%.5f = ", deX);
  }

  *Sigma += exch;
}

//******************************************************************************
void FeynmanSigma::prep_Feynman() {

  std::cout << "lmax = " << Angular::lFromIndex(m_max_kappaindex) << "\n";
  std::cout << "Using " << ParseEnum(m_Green_method)
            << " method for Green's functions\n";
  std::cout << "Using " << ParseEnum(m_Pol_method)
            << " method for Polarisation operator\n";
  std::cout << "Including from core n=" << m_min_core_n << "\n";

  if (m_screen_Coulomb)
    std::cout << "Including Coulomb screening\n";

  std::cout << "Using " << ParseEnum(m_ex_method) << " for exchange\n";

  // Print effective screening factors (if used)
  if (m_ex_method == ExchangeMethod::Goldstone && !m_fk.empty()) {
    std::cout << "w/ fk=";
    std::for_each(cbegin(m_fk), cend(m_fk),
                  [](auto fk) { std::cout << fk << ", "; });
    std::cout << "\n";
  }

  if (m_include_G)
    std::cout << "(Including FG/GF and GG)\n";

  form_Q_dr();
  form_Vx();
  form_Pa_core();
  setup_omega_grid();

  const auto max_k = std::min(m_maxk, m_k_cut);
  m_qpq_wk = form_QPQ_wk(max_k, m_Pol_method, m_omre, *m_wgridD);
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Q_dr() {

  std::vector<LinAlg::SqMatrix> qhat(std::size_t(m_maxk + 1), m_subgrid_points);

  LinAlg::SqMatrix t_dri{m_subgrid_points};

  // calculate (real) Qhat and Jacobian
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto sj = ri_subToFull(j);
      const auto drj = dr_subToFull(j);
      const auto rl = p_gr->r[std::min(si, sj)];
      const auto rm = p_gr->r[std::max(si, sj)];
      const auto ratio = rl / rm; // = r_< / r_>
      // q0 = (1.0 / rm) = r_<^k / r_>^k+1 ,  for k=0
      qhat[0][i][j] = (1.0 / rm) * dri * drj;
      t_dri[i][j] = dri;
      for (int k = 1; k <= m_maxk; ++k) {
        qhat[std::size_t(k)][i][j] = qhat[std::size_t(k - 1)][i][j] * ratio;
      }
    } // j
  }   // i

  // fill (complex) Qhat
  m_qhat.resize(std::size_t(m_maxk) + 1, {m_subgrid_points, m_include_G});
  for (auto k = 0; k <= m_maxk; ++k) {
    m_qhat[std::size_t(k)].ff =
        LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, qhat[std::size_t(k)]);
    if (m_include_G) {
      m_qhat[std::size_t(k)].gg = m_qhat[std::size_t(k)].ff;
    }
  }

  // fill (complex) Jacobian
  m_dri = std::make_unique<ComplexGMatrix>(m_subgrid_points, m_include_G);
  m_drj = std::make_unique<ComplexGMatrix>(m_subgrid_points, m_include_G);
  m_dri->ff = LinAlg::ComplexSqMatrix::make_complex({1.0, 0.0}, t_dri);
  m_drj->ff = m_dri->ff.transpose();
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Vx() {
  // Fill Vx:
  const auto max_kappa_index = std::size_t(m_max_kappaindex);
  m_Vxk.resize(max_kappa_index + 1, {m_subgrid_points, m_include_G});

  // If using a local method, exclude exchange here. Note: Only used for
  // testing!
  if (p_hf->excludeExchangeQ())
    return;

  for (auto ik = 0ul; ik <= max_kappa_index; ik++) {
    m_Vxk[ik] = calculate_Vx_kappa(Angular::kappaFromIndex(int(ik)));
  }
}

//------------------------------------------------------------------------------
GMatrix FeynmanSigma::calculate_Vx_kappa(int kappa) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  assert(m_dri != nullptr && "⚠️ form_Q_dr() must be run before form_Vx()");

  GMatrix Vx(m_subgrid_points, m_include_G);
  const auto tj = Angular::twoj_k(kappa);
  const auto kmax = Angular::twojFromIndex(m_max_kappaindex_core);
  const auto &core = p_hf->get_core();

  for (int k = 0; k <= kmax; ++k) {
    GMatrix Vx_k(m_subgrid_points, m_include_G);
    for (const auto &a : core) {
      const auto ck = Angular::Ck_kk(k, kappa, a.k);
      if (ck == 0.0)
        continue;
      const auto c_ang = -1.0 * ck * ck / double(tj + 1);
      addto_G(&Vx_k, a, a, c_ang);
    }
    Vx += Vx_k.mult_elements_by(get_qk(k).get_real());
    // Vx += Vx_k;
  }

  return Vx;
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Pa_core() {
  // Fill core |a><a|
  const auto &core = p_hf->get_core();
  m_Pa.resize(core.size(), {m_subgrid_points, m_include_G});
  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_Pa[ia] = G_single(core[ia], core[ia], ComplexDouble{1.0, 0.0});
  }
}

//------------------------------------------------------------------------------
void FeynmanSigma::setup_omega_grid() {
  // Set up imaginary frequency grid:
  std::cout << "Re(w) = " << m_omre << "\n";

  // used for direct:
  auto wmax_core = 0.0;
  const auto &core = p_hf->get_core();
  for (const auto &Fc : core) {
    if (Fc.n < m_min_core_n)
      continue;
    if (std::abs(Fc.en) > wmax_core)
      wmax_core = std::abs(Fc.en);
  }

  {
    const auto w0 = 0.01;
    // const auto wmax = 2000.0;
    const auto wmax = 150.0; // 1.2 * wmax_core;
    const double wratio = 1.5;
    const std::size_t wsteps = Grid::calc_num_points_from_du(
        w0, wmax, std::log(wratio), GridType::logarithmic);
    m_wgridD = std::make_unique<Grid>(w0, wmax, wsteps, GridType::logarithmic);
    std::cout << "Im(w) " << m_wgridD->gridParameters();
    printf(". r=%.2f\n", m_wgridD->r[1] / m_wgridD->r[0]);
  }
  m_wX_stride = 2; // stride for exchange XXX Make input!
}

//******************************************************************************
const GMatrix &FeynmanSigma::get_Vx_kappa(int kappa) const {
  const auto kappa_index = std::size_t(Angular::indexFromKappa(kappa));
  assert(kappa_index < m_Vxk.size());
  return m_Vxk[kappa_index];
}

const ComplexGMatrix &FeynmanSigma::get_qk(int k) const {
  assert(std::size_t(k) < m_qhat.size());
  return m_qhat[std::size_t(k)];
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
ComplexGMatrix FeynmanSigma::Green(int kappa, ComplexDouble en, States states,
                                   GrMethod method) const {
  if (states == States::core) {
    return Green_core(kappa, en);
  } else if (states == States::excited) {
    return Green_ex(kappa, en, method);
  }
  return (method == GrMethod::Green) ? Green_hf(kappa, en)
                                     : Green_hf_basis(kappa, en);
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::G_single(const DiracSpinor &ket,
                                      const DiracSpinor &bra,
                                      const ComplexDouble f) const {
  ComplexGMatrix Gmat(m_subgrid_points, m_include_G);
  const auto [x, iy] = f.unpack();

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto sj = ri_subToFull(j);
      Gmat.ff[i][j] =
          ComplexDouble(x * ket.f[si] * bra.f[sj], iy * ket.f[si] * bra.f[sj])
              .val;
    } // j
  }   // i

  if (m_include_G) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      const auto si = ri_subToFull(i);
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        const auto sj = ri_subToFull(j);
        Gmat.fg[i][j] =
            ComplexDouble(x * ket.f[si] * bra.g[sj], iy * ket.f[si] * bra.g[sj])
                .val;
        Gmat.gf[i][j] =
            ComplexDouble(x * ket.g[si] * bra.f[sj], iy * ket.g[si] * bra.f[sj])
                .val;
        Gmat.gg[i][j] =
            ComplexDouble(x * ket.g[si] * bra.g[sj], iy * ket.g[si] * bra.g[sj])
                .val;
      } // j
    }   // i
  }

  return Gmat;
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::Green_core(int kappa, ComplexDouble en) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "complex");
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.k=k
  ComplexGMatrix Gcore(m_subgrid_points, m_include_G);

  // loop over HF core, not Sigma core (used in subtraction to get G^excited)
  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.k != kappa)
      continue;
    const auto inv_de = (en - ComplexDouble{a.en}).inverse();
    Gcore += inv_de * m_Pa[ia]; // Pa = |a><a|
  }
  return Gcore;
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_ex(int kappa, ComplexDouble en,
                                      GrMethod method) const {
  ComplexGMatrix Gk(m_subgrid_points, m_include_G);

  if (method == GrMethod::basis) {
    Gk = Green_hf_basis(kappa, en, true);
  } else {
    // Subtract core states, by forceing Gk to be orthogonal to core:
    // Gk -> Gk - \sum_a|a><a|G
    Gk = Green_hf(kappa, en); // - Green_core(kappa, en);
    const auto &core = p_hf->get_core();
    const auto &drj = get_drj();
    for (auto ia = 0ul; ia < core.size(); ++ia) {
      if (core[ia].k != kappa)
        continue;
      if (core[ia].n < m_min_core_n)
        continue;
      auto pa = m_Pa[ia];       // copy
      pa.mult_elements_by(drj); // pa * Gk includes 'r_j' integral
      Gk -= (pa * Gk);
    }
  }

  return Gk;
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_hf_basis(int kappa, ComplexDouble en,
                                            bool ex_only) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  ComplexGMatrix Gc(m_subgrid_points, m_include_G);
  // const auto &core = p_hf->get_core(); // ?????
  const auto &core = m_holes; // p_hf->get_core(); // ?????
  // XXX Should include all states, even below n_min_core?
  // XXX BUT, total greens fn should be orthog!
  // No matter, since we don't actually use this method for green.....
  // (except for excited)
  const auto &ex = m_excited;
  for (const auto orbs : {&core, &ex}) {
    if (ex_only && orbs == &core)
      continue;
    for (const auto &a : *orbs) {
      if (a.k != kappa)
        continue;

      const auto inv_de = (en - ComplexDouble{a.en}).inverse();
      Gc += G_single(a, a, inv_de);
    }
  }
  return Gc;
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_hf(int kappa, ComplexDouble en) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Solve DE (no exchange), regular at 0, infinity ("pinf")
  DiracSpinor x0(0, kappa, p_gr);
  DiracSpinor xI(0, kappa, p_gr);
  const auto alpha = p_hf->get_alpha();
  const auto vl = p_hf->get_vlocal(x0.l());
  const auto &Hmag = p_hf->get_Hrad_mag(x0.l());
  DiracODE::regularAtOrigin(x0, en.re(), vl, Hmag, alpha);
  DiracODE::regularAtInfinity(xI, en.re(), vl, Hmag, alpha);

  // Evaluate Wronskian at ~65% of the way to pinf. Should be inependent of r
  const auto pp = std::size_t(0.65 * double(xI.pinf));
  const auto w = -1.0 * (xI.f[pp] * x0.g[pp] - x0.f[pp] * xI.g[pp]) / alpha;
  // Not sure why -ve sign here... ??? But needed to agree w/ basis version

  // Get G0 (Green's function, without exchange):
  const auto g0 = MakeGreensG0(x0, xI, w);
  // XXX NB: Have to modify Vx also in Breit case. XXX
  const auto &Vx = get_Vx_kappa(kappa);

  static const ComplexDouble one{1.0, 0.0}; // to convert real to complex

  // Include exchange, and imaginary energy part:

  if (en.im() == 0.0) {
    // G = [1 - G0*Vx]^{-1} * G0 = -[G0*Vx-1]^{-1} * G0
    // nb: much faster to invert _before_ make complex!
    // (but, only if imag. part is zero)
    return (-1 * one) * ((g0 * Vx).plusIdent(-1.0).invert() * g0);
  }

  // G0 := G0(re{e}) - no exchange, only real part
  // G(e) = [1 + i*Im{e}*G0 - G0*Vx]^{-1} * G0
  // Note: differential dr is included in Vx (via Q)
  ComplexDouble iw{0.0, en.im()};
  return ((iw * g0).mult_elements_by(*m_drj) - one * (g0 * Vx))
             .plusIdent(1.0)
             .invert() *
         (one * g0);
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::GreenAtComplex(const ComplexGMatrix &Gr,
                                            double om_imag) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Given G(wr) and wi, returns G(wr+i*wi)
  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = Gr * [1 + i*wi*Gr]^-1 = Gr * iX
  const ComplexDouble iw{0.0, om_imag};
  return ((iw * Gr).mult_elements_by(*m_drj).plusIdent(1.0).invert()) * Gr;
}

//------------------------------------------------------------------------------
GMatrix FeynmanSigma::MakeGreensG0(const DiracSpinor &x0, const DiracSpinor &xI,
                                   const double w) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Takes sub-grid into account; ket,bra are on full grid, G on sub-grid
  // G(r1,r2) = x0(rmin)*xI(imax)/w
  GMatrix g0I(m_subgrid_points, m_include_G);

  const auto winv = 1.0 / w;

  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto si = ri_subToFull(i);
    for (auto j = 0ul; j <= i; ++j) { // j <= i
      const auto sj = ri_subToFull(j);
      g0I.ff[i][j] = x0.f[sj] * xI.f[si] * winv;
      // g0I is symmetric
      g0I.ff[j][i] = g0I.ff[i][j];
    } // j
  }   // i

  if (m_include_G) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      const auto si = ri_subToFull(i);
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        const auto sj = ri_subToFull(j);
        const auto irmin = std::min(sj, si);
        const auto irmax = std::max(sj, si);
        g0I.fg[i][j] = x0.f[irmin] * xI.g[irmax] * winv;
        // fg = gf?
        g0I.gf[i][j] = x0.g[irmin] * xI.f[irmax] * winv;
        g0I.gg[i][j] = x0.g[irmin] * xI.g[irmax] * winv;
      } // j
    }   // i
  }

  return g0I;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
ComplexGMatrix FeynmanSigma::Polarisation_k(int k, ComplexDouble omega,
                                            GrMethod method) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  ComplexGMatrix pi_k(m_subgrid_points, m_include_G);

  static const auto Iunit = ComplexDouble{0.0, 1.0};
  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.n < m_min_core_n) // ?? ok?
      continue;
    const auto &pa = m_Pa[ia]; // |a><a|

    for (int ialpha = 0; ialpha <= m_max_kappaindex; ++ialpha) {
      const auto kA = Angular::kappaFromIndex(ialpha);
      const auto ck_aA = Angular::Ck_kk(k, a.k, kA);
      if (ck_aA == 0.0)
        continue;
      const double c_ang = ck_aA * ck_aA / double(2 * k + 1);
      // XXX Inlcude the 2 * k + 1 here ??

      const auto ea_minus_w = ComplexDouble{a.en} - omega;
      const auto ea_plus_w = ComplexDouble{a.en} + omega;

      pi_k += c_ang * (Green_ex(kA, ea_minus_w, method) +
                       Green_ex(kA, ea_plus_w, method))
                          .mult_elements_by(pa);
    }
  }
  return Iunit * pi_k;
}

//******************************************************************************
std::vector<ComplexGMatrix> FeynmanSigma::OneMinusPiQInv(
    const std::vector<std::vector<ComplexGMatrix>> &pi_wk) const {
  // Function of omega!
  // XXX Note sure if right.. depends on pi angular reduction.
  // X_PiQ = [1-PiQ]^{-1}

  std::vector<ComplexGMatrix> X_PiQ(pi_wk.size(),
                                    {m_subgrid_points, m_include_G});

#pragma omp parallel for
  for (auto iw = 0ul; iw < pi_wk.size(); ++iw) {
    auto &X_PiQ_w = X_PiQ[iw];
    for (auto k = 0ul; k < pi_wk[iw].size(); ++k) {
      const auto &pik = pi_wk[iw][k];
      const auto &qk = get_qk(int(k));
      X_PiQ_w -= pik * qk; //-ve, since [1-PiQ]
    }
    // [1-PiQ]^{-1}
    X_PiQ_w.plusIdent(1.0).invert();
  }
  return X_PiQ;
}

//******************************************************************************
std::vector<std::vector<ComplexGMatrix>>
FeynmanSigma::make_pi_wk(int max_k, GrMethod pol_method, double omre,
                         const Grid &wgrid) const {

  const auto num_ks = std::size_t(max_k + 1);
  std::vector<std::vector<ComplexGMatrix>> pi_wk(wgrid.num_points);

#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points; ++iw) {
    const auto omega = ComplexDouble{omre, wgrid.r[iw]};
    pi_wk[iw].reserve(num_ks);
    for (auto k = 0ul; k < num_ks; ++k) {
      pi_wk[iw].push_back(Polarisation_k(int(k), omega, pol_method));
    }
  }
  return pi_wk;
}

//******************************************************************************
std::vector<std::vector<ComplexGMatrix>>
FeynmanSigma::form_QPQ_wk(int max_k, GrMethod pol_method, double omre,
                          const Grid &wgrid) const {
  // Returns QPQ, function of w and k

  std::vector<std::vector<ComplexGMatrix>> qpq(wgrid.num_points);
  std::cout << "Forming QPQ(w,k) matrix.." << std::flush;

  const auto pi_wk = make_pi_wk(max_k, pol_method, omre, wgrid);
  std::cout << "." << std::flush;
  // X_PiQ = [1-PiQ]^{-1}
  // Only calculate X_PiQ if doing screening
  const auto X_PiQ =
      m_screen_Coulomb ? OneMinusPiQInv(pi_wk) : std::vector<ComplexGMatrix>{};
  std::cout << "." << std::flush;

  const auto num_ks = std::size_t(max_k + 1);
#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points; ++iw) {
    qpq[iw].reserve(num_ks);
    for (auto k = 0ul; k < num_ks; ++k) {
      // q*p*q => q*X*p*q, x = [1-Pi*Q]^(-1)
      const auto &qk = get_qk(int(k));
      const auto &pi = pi_wk[iw][k];
      if (m_screen_Coulomb)
        qpq[iw].emplace_back(qk * X_PiQ[iw] * pi * qk);
      else
        qpq[iw].emplace_back(qk * pi * qk);
    }
  }
  std::cout << "..done\n";
  return qpq;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
GMatrix FeynmanSigma::FeynmanDirect(int kv, double env) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  GMatrix Sigma(m_subgrid_points, m_include_G);

  static const ComplexDouble I{0.0, 1.0};

  // // Set up imaginary frequency grid:
  const double omre = m_omre;
  const auto &wgrid = *m_wgridD;

  const auto max_k = std::min(m_maxk, m_k_cut);

  // Store gAs in advance
  // Need to do for each kappa, so fine to do in here..
  // (Unless we use same grid for exchange!)
  const auto num_kappas = std::size_t(m_max_kappaindex + 1);
  std::vector<std::vector<ComplexGMatrix>> gBs(num_kappas);
#pragma omp parallel for
  for (auto ik = 0ul; ik < num_kappas; ++ik) {
    const auto kappa = Angular::kappaFromIndex(int(ik));
    gBs[ik].reserve(wgrid.num_points);
    for (auto iw = 0ul; iw < wgrid.num_points; iw++) {
      ComplexDouble evpw{env + omre, wgrid.r[iw]};
      gBs[ik].push_back(Green(kappa, evpw, States::both, m_Green_method));
    }
  }

#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points; iw++) { // for omega integral
    // I, since dw is on imag. grid; 2 from symmetric +/- w
    const auto dw = I * (2.0 * wgrid.drdu[iw] * wgrid.du / (2 * M_PI));

    for (auto k = 0ul; k <= std::size_t(max_k); k++) {

      const auto qpq_dw = dw * m_qpq_wk[iw][k];

      for (auto iB = 0ul; iB < num_kappas; ++iB) {
        const auto kB = Angular::kappaFromIndex(int(iB));
        const auto ck_vB = Angular::Ck_kk(int(k), kv, kB);
        if (ck_vB == 0.0)
          continue;
        const double c_ang = ck_vB * ck_vB / double(Angular::twoj_k(kv) + 1);

        auto gB = gBs[iB][iw]; // copy

#pragma omp critical(sum_sigma_d)
        { Sigma += c_ang * (gB.mult_elements_by(qpq_dw)).get_real(); }

      } // k
    }   // omega
  }     // beta

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sigma.ff[i][j] /= (dri * drj);
      if (m_include_G) {
        Sigma.fg[i][j] /= (dri * drj);
        Sigma.gf[i][j] /= (dri * drj);
        Sigma.gg[i][j] /= (dri * drj);
      }
    }
  }

  return Sigma;
} // namespace MBPT

//******************************************************************************
GMatrix FeynmanSigma::FeynmanEx_w1w2(int kv, double en_r) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Note: Takes really long time, and doesn't work! XXX

  GMatrix Sx(m_subgrid_points, m_include_G);

  // Set up imaginary frequency grid:
  const double omre = m_omre;
  const auto &wgrid = *m_wgridD;
  static const ComplexDouble I{0.0, 1.0};
  const ComplexDouble en{en_r, 0.0};
  const auto tjvp1 = Angular::twoj_k(kv) + 1;

  const auto max_k = std::min(m_maxk, m_k_cut);

  // Store gs in advance
  // Note: gB depends on w1 and w2!
  const auto num_kappas = std::size_t(m_max_kappaindex + 1);
  std::vector<std::vector<ComplexGMatrix>> Gs(num_kappas);
  std::vector<std::vector<std::vector<ComplexGMatrix>>> GBs(num_kappas);
#pragma omp parallel for
  for (auto ik = 0ul; ik < num_kappas; ++ik) {
    const auto kappa = Angular::kappaFromIndex(int(ik));
    Gs[ik].reserve(wgrid.num_points);
    GBs[ik].resize(wgrid.num_points);
    for (auto iw1 = 0ul; iw1 < wgrid.num_points; iw1++) {
      ComplexDouble evpw{en_r + omre, wgrid.r[iw1]};
      Gs[ik].push_back(Green(kappa, evpw, States::both, m_Green_method));
      GBs[ik][iw1].reserve(wgrid.num_points);
      for (auto iw2 = 0ul; iw2 < wgrid.num_points; iw2++) {
        ComplexDouble evpw1pw2{en_r + 2.0 * omre, wgrid.r[iw1] + wgrid.r[iw2]};
        GBs[ik][iw1].push_back(
            Green(kappa, evpw1pw2, States::both, m_Green_method));
      }
    }
  }

  // Store parts of Sx seperately, for more efficient parallelisation
  std::vector<GMatrix> Sxs(num_kappas, {m_subgrid_points, m_include_G});

#pragma omp parallel for
  for (auto iA = 0ul; iA < num_kappas; ++iA) { // alpha
    auto &ScA = Sxs[iA];
    for (auto iB = 0ul; iB < num_kappas; ++iB) {   // beta
      for (auto iG = 0ul; iG < num_kappas; ++iG) { // gamma
        const auto kA = Angular::kappaFromIndex(int(iA));
        const auto kB = Angular::kappaFromIndex(int(iB));
        const auto kG = Angular::kappaFromIndex(int(iG));

        for (auto iw1 = 0ul; iw1 < wgrid.num_points; iw1 += m_wX_stride) { // w1
          const auto dw1 = 2.0 * double(m_wX_stride) * wgrid.drdu[iw1] *
                           wgrid.du / (2 * M_PI);

          for (auto iw2 = 0ul; iw2 < wgrid.num_points;
               iw2 += m_wX_stride) { // w2
            const auto dw2 = 2.0 * double(m_wX_stride) * wgrid.drdu[iw2] *
                             wgrid.du / (2 * M_PI);

            // get the Green's functions:
            const auto &gA = Gs[iA][iw1];
            const auto &gB = GBs[iB][iw1][iw2];
            const auto &gG = Gs[iG][iw2];

            const auto gqgqg = sumkl_gqgqg(gA, gB, gG, kv, kA, kB, kG, max_k);
            ScA -= (dw1 * dw2 / tjvp1) * gqgqg;

          } // w2
        }   // w1
      }     // gamma
    }       // beta
  }         // alpha

  for (const auto &SxA : Sxs) {
    Sx += SxA;
  }

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sx.ff[i][j] /= (dri * drj);
      // if (m_include_G) {
      //   Sx.fg[i][j] /= (dri * drj);
      //   Sx.gf[i][j] /= (dri * drj);
      //   Sx.gg[i][j] /= (dri * drj);
      // }
    }
  }

  return Sx;
}

//******************************************************************************
GMatrix FeynmanSigma::FeynmanEx_1(int kv, double env) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  GMatrix Sx1(m_subgrid_points, m_include_G);

  // Set up imaginary frequency grid:
  const double omre = m_omre; // does seem to depend on this..
  const auto &wgrid = *m_wgridD;
  const auto &core = p_hf->get_core();
  const auto tjvp1 = Angular::twoj_k(kv) + 1;

  // Store gAs in advance, since these depend only on w and kA, not a, kB or k
  const auto num_kappas = std::size_t(m_max_kappaindex + 1);
//   std::vector<std::vector<ComplexGMatrix>> gAs(num_kappas);
// #pragma omp parallel for
//   for (auto ik = 0ul; ik < num_kappas; ++ik) {
//     const auto kappa = Angular::kappaFromIndex(int(ik));
//     gAs[ik].reserve(wgrid.num_points);
//     for (auto iw = 0ul; iw < wgrid.num_points; iw++) {
//       ComplexDouble evpw{env + omre, wgrid.r[iw]};
//       gAs[ik].push_back(Green(kappa, evpw, States::both, m_Green_method));
//     }
//   }

// #pragma omp parallel for collapse(2)
#pragma omp parallel for
  for (auto iB = 0ul; iB < num_kappas; ++iB) {
    const auto kB = Angular::kappaFromIndex(int(iB));
    for (auto iw = 0ul; iw < wgrid.num_points; iw += m_wX_stride) {

      auto omim = wgrid.r[iw]; // XXX Symmetric?? Or Not??
      const auto omega = ComplexDouble{omre, omim};
      const auto dw1 =
          -2.0 * double(m_wX_stride) * wgrid.drdu[iw] * wgrid.du / (2.0 * M_PI);

      const auto *const qpqw_k = m_screen_Coulomb ? &m_qpq_wk[iw] : nullptr;

      for (auto ia = 0ul; ia < core.size(); ++ia) {
        const auto &Fa = core[ia];
        if (Fa.n < m_min_core_n)
          continue;
        const auto &pa = m_Pa[ia];
        const auto ea = ComplexDouble{Fa.en, 0.0};

        const auto gxBm = Green_ex(kB, ea - omega, m_Green_method);
        const auto gxBp = Green_ex(kB, ea + omega, m_Green_method);

        for (auto iA = 0ul; iA < num_kappas; ++iA) {
          const auto kA = Angular::kappaFromIndex(int(iA));

          // const auto &gA = gAs[iA][iw];
          const auto gA =
              Green(kA, {env + omre, omim}, States::both, m_Green_method);
          const auto gqpg =
              sumkl_GQPGQ(gA, gxBm, gxBp, pa, kv, kA, kB, Fa.k, qpqw_k);

#pragma omp critical(sum_x1)
          { Sx1 += (dw1 / tjvp1) * gqpg; }
        } // alpha
      }   // a
    }
  }

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sx1.ff[i][j] /= (dri * drj);
      // No g
    }
  }

  return Sx1;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
GMatrix FeynmanSigma::sumkl_gqgqg(const ComplexGMatrix &gA,
                                  const ComplexGMatrix &gB,
                                  const ComplexGMatrix &gG, int kv, int kA,
                                  int kB, int kG, int kmax) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // EXCHANGE part, used in w1w2 version

  auto gqgqg = GMatrix(m_subgrid_points, m_include_G);

  for (int k = 0; k <= kmax; k++) { // k (k1)
    const auto &qk = get_qk(k);

    for (int l = 0; l <= kmax; l++) { // l (k2)
      const auto &ql = get_qk(l);

      // tensor_5_product adds the real part of below to result
      // Sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ]
      const auto Lkl = Lkl_abcd(k, l, kv, kB, kA, kG);
      if (Lkl != 0.0) {
        const auto s = Angular::neg1pow(k + l);
        const auto sLkl = ComplexDouble{s * Lkl, 0.0};
        tensor_5_product(&gqgqg, sLkl, qk, gB, gG, gA, ql);
      }
    } // l
  }   // k

  return gqgqg;
}

//******************************************************************************
GMatrix FeynmanSigma::sumkl_GQPGQ(
    const ComplexGMatrix &gA, const ComplexGMatrix &gxBm,
    const ComplexGMatrix &gxBp, const ComplexGMatrix &pa, int kv, int kA,
    int kB, int ka, const std::vector<ComplexGMatrix> *const qpqw_k) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // EXCHANGE part, used in w1 version
  // (w2 was integrated over analytically)

  GMatrix sum_GQPG{m_subgrid_points, m_include_G};

  const auto kmax = std::min(m_maxk, m_k_cut);

  for (auto k = 0; k <= kmax; ++k) {
    // qk -> qk + 2 * QPxQ ??
    const auto &tqk = get_qk(k);
    // Try to screen q^k(w1) - wrong sign?? Anti-screening??
    const auto qk = qpqw_k != nullptr ? tqk + (*qpqw_k)[std::size_t(k)] : tqk;

    for (auto l = 0; l <= kmax; ++l) {
      const auto &ql = get_qk(l);

      const auto s0 = Angular::neg1pow(k + l);

      const auto L1 = Lkl_abcd(k, l, kv, kB, kA, ka);
      const auto L2 = Lkl_abcd(k, l, kv, ka, kA, kB);

      const auto ic1 = ComplexDouble(s0 * L1, 0.0);
      const auto ic2 = ComplexDouble(s0 * L2, 0.0);

      // tensor_5_product adds the real part of below to result
      // Sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ]

      if (L1 != 0.0)
        tensor_5_product(&sum_GQPG, ic1, qk, gxBp, pa, gA, ql);

      if (L2 != 0.0)
        tensor_5_product(&sum_GQPG, ic2, qk, pa, gxBm, gA, ql);

    } // l
  }   // k

  return sum_GQPG;
}

//******************************************************************************
void FeynmanSigma::tensor_5_product(
    GMatrix *result, const ComplexDouble &factor, const ComplexGMatrix &a,
    const ComplexGMatrix &b, const ComplexGMatrix &c, const ComplexGMatrix &d,
    const ComplexGMatrix &e) const {
  // Adds real part of below to result
  // Sum_ij [ factor * a1j * bij * cj2 * (d_1i * e_i2) ]
  const auto size = result->size;
  ComplexDouble sum_j, sum_ij;

  // The a*c part depends only on j (not i).
  // Doing this mult early saves factor of 'size' complex multiplications
  // Leads to a >2x speed-up!
  std::vector<ComplexDouble> ac(size); // see below

  for (auto r1 = 0ul; r1 < size; ++r1) {
    for (auto r2 = 0ul; r2 < size; ++r2) {

      for (auto j = 0ul; j < size; ++j) {
        // early a*c mult
        ac[j] = ComplexDouble(a.ff[r1][j]) * c.ff[j][r2];
      }
      sum_ij = {0.0, 0.0};
      for (auto i = 0ul; i < size; ++i) {
        sum_j = {0.0, 0.0};
        for (auto j = 0ul; j < size; ++j) {
          sum_j += ac[j] * b.ff[i][j];
        }
        sum_ij += sum_j * d.ff[r1][i] * e.ff[i][r2];
      }
      result->ff[r1][r2] += (factor * sum_ij).cre();

    } // r2
  }   // r1
}

//******************************************************************************
double FeynmanSigma::Lkl_abcd(int k, int l, int ka, int kb, int kc,
                              int kd) const {

  const auto &Ck = m_yeh.Ck();

  const auto Ckac = Ck(k, ka, kc);
  const auto Ckbd = Ck(k, kb, kd);
  const auto Clad = Ck(l, ka, kd);
  const auto Clbc = Ck(l, kb, kc);
  if (Ckac == 0.0 || Ckbd == 0.0 || Clad == 0.0 || Clbc == 0.0)
    return 0.0;

  const auto tja = Angular::twoj_k(ka);
  const auto tjb = Angular::twoj_k(kb);
  const auto tjc = Angular::twoj_k(kc);
  const auto tjd = Angular::twoj_k(kd);
  return m_6j(tja, tjc, tjb, tjd, k, l) * (Ckac * Ckbd * Clad * Clbc);
}

//******************************************************************************
GMatrix FeynmanSigma::Exchange_Goldstone(const int kappa,
                                         const double en) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // XXX Update Goldstone to work the same way as this one.
  // And, make me able to call it!!

  // Two second-order exchange diagrams:
  // Diagram (b) (exchange):
  // |Q^k_amn><P^k_amn| / de_amn / [k][j]
  // Diagram (d) (exchange):
  // |Q^k_nba><P^k_nba| / de_nba / [k][j]
  // where:
  // All indeces are summed over,
  // a & b are core states, n & m are virtual excited states,
  // k is multipolarity [Coloulmb expansion]
  // de_xyz = e_v + e_x - e_y - e_z

  GMatrix SxG(m_subgrid_points, m_include_G);

  const auto &Ck = m_yeh.Ck();

  std::vector<GMatrix> Gxs(m_holes.size(), {m_subgrid_points, m_include_G});
#pragma omp parallel for
  for (auto ia = 0ul; ia < m_holes.size(); ia++) {
    const auto &a = m_holes[ia];
    auto &Ga_x = Gxs[ia];
    auto Qkv = DiracSpinor(0, kappa, p_gr); // re-use to reduce alloc'ns
    auto Pkv = DiracSpinor(0, kappa, p_gr); // re-use to reduce alloc'ns
    for (const auto &n : m_excited) {
      const auto [kmin_nb, kmax_nb] = m_yeh.k_minmax(n, a);
      for (int k = kmin_nb; k <= kmax_nb; ++k) {
        if (Ck(k, a.k, n.k) == 0)
          continue;
        const auto f_kkjj = (2 * k + 1) * (Angular::twoj_k(kappa) + 1);
        const auto &yknb = m_yeh(k, n, a);

        // Effective screening parameter:
        const auto fk = get_fk(k);

        // Diagrams (b) [exchange]
        for (const auto &m : m_excited) {
          if (Ck(k, kappa, m.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, a, m, n, k, yknb, Ck);
          // Pkv_bcd_2 allows different screening factor for each 'k2' in exch.
          Coulomb::Pkv_bcd_2(&Pkv, a, m, n, k, m_yeh(m, a), Ck, m_6j, m_fk);
          const auto dele = en + a.en - m.en - n.en;
          const auto factor = fk / (f_kkjj * dele);
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // m

        // Diagrams (d) [exchange]
        for (const auto &b : m_holes) {
          if (Ck(k, kappa, b.k) == 0)
            continue;
          Coulomb::Qkv_bcd(&Qkv, n, b, a, k, yknb, Ck);
          Coulomb::Pkv_bcd_2(&Pkv, n, b, a, k, m_yeh(n, b), Ck, m_6j, m_fk);
          const auto dele = en + n.en - b.en - a.en;
          const auto factor = fk / (f_kkjj * dele); // XXX
          addto_G(&Ga_x, Qkv, Pkv, factor);
        } // b

      } // k
    }   // n
  }     // a

  for (const auto &Gx : Gxs)
    SxG += Gx;

  return SxG;
}

} // namespace MBPT
