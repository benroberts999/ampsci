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
      screen_Coulomb(sigp.screenCoulomb),
      m_omre(sigp.real_omega),
      basis_for_Green(sigp.GreenBasis),
      basis_for_Pol(sigp.PolBasis),
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
  //   return; //?

  if (!read_ok) {
    std::cout << "Form correlation potential: ";
    std::cout << "Feynman method\n";
    // Extand 6j and Ck
    m_6j.fill(m_maxk, m_maxk);
    m_yeh.extend_Ck(m_maxk, m_maxk);
    prep_Feynman();
    if (m_include_G)
      std::cout << "(Including FG/GF and GG)\n";
    form_Sigma(en_list, fname);
  }
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
void FeynmanSigma::prep_Feynman() {

  // Do intial setup:

  std::cout << "lmax = " << Angular::lFromIndex(m_max_kappaindex) << ", ";
  if (basis_for_Green)
    std::cout << "Using basis for Green's fns; ";
  else
    std::cout << "Using Green method for Green's fn; ";
  if (basis_for_Pol)
    std::cout << "basis for Polarisation op.\n";
  else
    std::cout << "Green method for Polarisation op.\n";
  std::cout << "\n";

  form_Q_dr();
  form_Vx();
  form_Pa_core();
  setup_omega_grid();
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

  const auto tj = Angular::twoj_k(kappa);
  GMatrix Vx(m_subgrid_points, m_include_G);
  const auto &core = p_hf->get_core();

  const auto kmax = Angular::twojFromIndex(m_max_kappaindex_core);
  for (int k = 0; k <= kmax; ++k) {
    GMatrix Vx_k(m_subgrid_points, m_include_G);
    for (const auto &a : core) {
      const auto ck = Angular::Ck_kk(k, kappa, a.k); // XXX lookup
      if (ck == 0.0)
        continue;
      const auto c_ang = -1.0 * ck * ck / double(tj + 1);
      addto_G(&Vx_k, a, a, c_ang);
    }
    Vx_k.mult_elements_by(get_qk(k).get_real());
    Vx += Vx_k;
  }

  return Vx;
}

//------------------------------------------------------------------------------
void FeynmanSigma::form_Pa_core() {
  // Fill core |a><a|
  const auto &core = p_hf->get_core();
  m_Pa.resize(core.size(), {m_subgrid_points, m_include_G});
  for (auto ia = 0ul; ia < core.size(); ia++) {
    m_Pa[ia] = G_single(core[ia], core[ia], {1.0, 0.0});
  }
}

//------------------------------------------------------------------------------
void FeynmanSigma::setup_omega_grid() {
  // Set up imaginary frequency grid:
  std::cout << "Re(w) = " << m_omre << "\n";
  {
    const auto w0 = 0.05;
    const auto wmax = 50.0;
    const double wratio = 2.0;
    const std::size_t wsteps = Grid::calc_num_points_from_du(
        w0, wmax, std::log(wratio), GridType::logarithmic);
    // const auto wgrid = Grid(w0, wmax, wsteps, GridType::logarithmic);
    m_wgridD = std::make_unique<Grid>(w0, wmax, wsteps, GridType::logarithmic);
    std::cout << "Dir | Im(w) " << m_wgridD->gridParameters();
    printf(". r=%.2f\n", m_wgridD->r[1] / m_wgridD->r[0]);
  }
  {
    // exchange seems quite insensitive to grid!
    const auto w0 = 0.05;
    const auto wmax = 30.0;
    const double wratio = 3.3;
    const std::size_t wsteps = Grid::calc_num_points_from_du(
        w0, wmax, std::log(wratio), GridType::logarithmic);
    m_wgridX = std::make_unique<Grid>(w0, wmax, wsteps, GridType::logarithmic);
    std::cout << "Exch| Im(w) " << m_wgridX->gridParameters();
    printf(". r=%.2f\n", m_wgridX->r[1] / m_wgridX->r[0]);
  }
}

//******************************************************************************
//******************************************************************************
void FeynmanSigma::fill_Sigma_k(GMatrix *Sigma, const int kappa,
                                const double en) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  auto find_kappa = [=](const auto &f) { return f.k == kappa; };
  const auto Fk = std::find_if(cbegin(m_excited), cend(m_excited), find_kappa);

  // const auto dir
  const auto tmp = FeynmanDirect(kappa, en);
  *Sigma = tmp;

  auto deD = 0.0;
  if (Fk != cend(m_excited)) {
    deD = *Fk * Sigma_G_Fv(*Sigma, *Fk);
    printf("de= %.4f + ", deD);
    std::cout << std::flush;
  }

  const auto exch = FeynmanEx_1(kappa, en);

  if (Fk != cend(m_excited)) {
    auto deX = *Fk * Sigma_G_Fv(exch, *Fk);
    printf("%.5f = ", deX);
  }

  *Sigma += exch;
}

//******************************************************************************
//******************************************************************************

//******************************************************************************
ComplexGMatrix FeynmanSigma::Green(int kappa, double en_re, double en_im,
                                   States states, GrMethod method) const {
  if (states == States::core) {
    return Green_core(kappa, en_re, en_im);
  } else if (states == States::excited) {
    Green_ex(kappa, en_re, en_im, method);
  }
  return (method == GrMethod::Green) ? Green_hf(kappa, en_re, en_im)
                                     : Green_hf_basis(kappa, en_re, en_im);
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
ComplexGMatrix FeynmanSigma::Green_core(int kappa, double en_re,
                                        double en_im) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "complex");
  // G_core = \sum_a |a><a|/(e_r + i*e_i-ea), for all a with a.k=k
  ComplexGMatrix Gcore(m_subgrid_points, m_include_G);

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
ComplexGMatrix FeynmanSigma::Green_ex(int kappa, double en_re, double en_im,
                                      GrMethod method) const {
  if (method == GrMethod::basis) {
    return Green_hf_basis(kappa, en_re, en_im, true);
  }
  return (Green_hf(kappa, en_re, en_im) - Green_core(kappa, en_re, en_im));
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_hf_basis(int kappa, double en_re,
                                            double en_im, bool ex_only) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  ComplexGMatrix Gc(m_subgrid_points, m_include_G);
  // const auto &core = p_hf->get_core(); // ?????
  const auto &core = m_holes; // p_hf->get_core(); // ?????
  // XXX Should include all states, even below n_min_core?
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
ComplexGMatrix FeynmanSigma::Green_hf(int kappa, double en_re,
                                      double en_im) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__, "complexE");
  return en_im == 0.0 ? Green_hf(kappa, en_re)
                      : GreenAtComplex(Green_hf(kappa, en_re), en_im);
}

//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Green_hf(int kappa, double en) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // Solve DE (no exchange), regular at 0, infinity ("pinf")
  DiracSpinor x0(0, kappa, p_gr);
  DiracSpinor xI(0, kappa, p_gr);
  const auto alpha = p_hf->get_alpha();
  const auto vl = p_hf->get_vlocal(x0.l());
  const auto &Hmag = p_hf->get_Hrad_mag(x0.l());
  DiracODE::regularAtOrigin(x0, en, vl, Hmag, alpha);
  DiracODE::regularAtInfinity(xI, en, vl, Hmag, alpha);

  // Evaluate Wronskian at ~65% of the way to pinf. Should be inependent of r
  const auto pp = std::size_t(0.65 * double(xI.pinf));
  const auto w = -1.0 * (xI.f[pp] * x0.g[pp] - x0.f[pp] * xI.g[pp]) / alpha;
  // not sure why -ve sign.. is sign even defined by DiracODE?

  // Get G0 (Green's function, without exchange):
  const auto g0 = MakeGreensG0(x0, xI, w);
  const auto &Vx = get_Vx_kappa(kappa);

  // Include exchange:
  // G = [1 - G0*Vx]^{-1} * G0 = -[G0*Vx-1]^{-1} * G0
  // nb: much faster to invert _before_ make complex!
  return ((g0 * Vx).plusIdent(-1.0).invert() * g0).make_complex({-1.0, 0.0});
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
    } // j
  }   // i
  g0I.ff.fill_symmetric_lower();

  if (m_include_G) {
    for (auto i = 0ul; i < m_subgrid_points; ++i) {
      const auto si = ri_subToFull(i);
      for (auto j = 0ul; j < m_subgrid_points; ++j) {
        const auto sj = ri_subToFull(j);
        const auto irmin = std::min(sj, si);
        const auto irmax = std::max(sj, si);
        g0I.fg[i][j] = x0.f[irmin] * xI.g[irmax] * winv;
        g0I.gf[i][j] = x0.g[irmin] * xI.f[irmax] * winv;
        g0I.gg[i][j] = x0.g[irmin] * xI.g[irmax] * winv;
      } // j
    }   // i
  }

  return g0I;
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
ComplexGMatrix FeynmanSigma::GreenAtComplexShift(const ComplexGMatrix &Gr,
                                                 double de_re,
                                                 double de_im) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // Given G(wr) and wi, returns G(wr+i*wi)
  // G(w) =  G(re(w)+im(w)) ;  Gr = G(re(w)), G = G(w),   im(w) = wi
  // G = Gr * [1 + i*wi*Gr]^-1 = Gr * iX
  // Note: only works when de_re=0.0 ??
  const ComplexDouble dele{de_re, de_im};
  return ((dele * Gr).mult_elements_by(*m_drj).plusIdent(1.0).invert()) * Gr;
}

//******************************************************************************
//******************************************************************************

ComplexGMatrix FeynmanSigma::Polarisation_k(int k, double omre, double omim,
                                            GrMethod method) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  ComplexGMatrix pi_k(m_subgrid_points, m_include_G);
  static const auto Iunit = ComplexDouble{0.0, 1.0};
  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.n < m_min_core_n)
      continue;
    const auto &pa = m_Pa[ia]; // |a><a|

    for (int ialpha = 0; ialpha <= m_max_kappaindex; ++ialpha) {
      const auto kA = Angular::kappaFromIndex(ialpha);
      const auto ck_aA = Angular::Ck_kk(k, a.k, kA);
      if (ck_aA == 0.0)
        continue;
      const double c_ang = ck_aA * ck_aA / double(2 * k + 1);
      // XXX Inlcude the 2 * k + 1 ??

      // pi_k += c_ang * Polarisation_a(pa, a.en, kA, om_re, om_im, method);
      pi_k += c_ang * (Green_ex(kA, a.en - omre, -omim, method) +
                       Green_ex(kA, a.en + omre, omim, method))
                          .mult_elements_by(pa);
    }
  }
  return Iunit * pi_k;
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::Polarisation(int k_a, int k_alpha, double om_re,
                                          double om_im, GrMethod method) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  ComplexGMatrix pi(m_subgrid_points, m_include_G);

  const auto &core = p_hf->get_core();
  for (auto ia = 0ul; ia < core.size(); ++ia) {
    const auto &a = core[ia];
    if (a.n < m_min_core_n || a.k != k_a)
      continue;
    const auto &pa = m_Pa[ia]; // |a><a|
    pi += Polarisation_a(pa, a.en, k_alpha, om_re, om_im, method);
  }
  return pi;
}
//------------------------------------------------------------------------------
ComplexGMatrix FeynmanSigma::Polarisation_a(const ComplexGMatrix &pa,
                                            double ena, int kA, double omre,
                                            double omim,
                                            GrMethod method) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  static const auto Iunit = ComplexDouble{0.0, 1.0};

  return Iunit * ((Green_ex(kA, ena - omre, -omim, method) +
                   Green_ex(kA, ena + omre, omim, method))
                      .mult_elements_by(pa));

  // // mutliply i unit through to demoninator (* and / by (-i))
  // // NO! Probably works for basis? But not otherwise (poles!)
  // return (Green_ex(kA, -omim, -ena + omre, method) +
  //         Green_ex(kA, omim, -ena - omre, method))
  //     .mult_elements_by(pa);
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::screenedCoulomb(const ComplexGMatrix &q,
                                             const ComplexGMatrix &pi) const {
  // not checked!
  // XXX Are the dr_i, dr_j parts correct??
  // Qscr = [1 - Q*Pi]^{-1} * Q = -[Q*Pi - 1]^{-1} * Q
  return -1 * ((q * pi).plusIdent(-1.0).invert() * q);
  // return -1 * q * ((pi * q).plusIdent(-1.0).invert());
}

//******************************************************************************
ComplexGMatrix FeynmanSigma::sum_qpq(int k, double om_re, double om_im) const {

  const auto pol_method = basis_for_Pol ? GrMethod::basis : GrMethod::Green;
  auto pi = Polarisation_k(k, om_re, om_im, pol_method);
  const auto &q = get_qk(k);
  const auto &q_scr = screen_Coulomb ? screenedCoulomb(q, pi) : q; //???
  return q * pi * q_scr;
}

//******************************************************************************
GMatrix FeynmanSigma::FeynmanDirect(int kv, double env) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  GMatrix Sigma(m_subgrid_points, m_include_G);

  // // Set up imaginary frequency grid:
  const double omre = m_omre; // does seem to depend on this..
  const auto &wgrid = *m_wgridD;
  const auto pol_method = basis_for_Pol ? GrMethod::basis : GrMethod::Green;

  // Greens function (at om_re) remains same, so calculate it once only:
  std::vector<ComplexGMatrix> gBetas;
  if (!basis_for_Green) {
    gBetas.reserve(std::size_t(m_max_kappaindex + 1));
    for (int ibeta = 0; ibeta <= m_max_kappaindex; ++ibeta) {
      const auto kB = Angular::kappaFromIndex(ibeta);
      gBetas.push_back(Green_hf(kB, env + omre));
    }
  }

  std::vector<GMatrix> Sigma_Ws(std::size_t(wgrid.num_points),
                                {m_subgrid_points, m_include_G});
#pragma omp parallel for
  for (auto iw = 0ul; iw < wgrid.num_points; iw++) { // for omega integral
    const auto omim = wgrid.r[iw];
    const auto dw = wgrid.drdu[iw] * wgrid.du;
    auto &Sigma_w = Sigma_Ws[iw];

    for (int k = 0; k <= m_maxk; k++) {

      const auto pi = Polarisation_k(k, omre, omim, pol_method);
      const auto &q = get_qk(k);
      const auto &q_scr = screen_Coulomb ? screenedCoulomb(q, pi) : q; //???
      const auto qpq_dw = (dw / M_PI) * (q * pi * q_scr);

      for (int ibeta = 0; ibeta <= m_max_kappaindex; ++ibeta) {
        const auto kB = Angular::kappaFromIndex(ibeta);
        const auto ck_vB = Angular::Ck_kk(k, kv, kB);
        if (ck_vB == 0.0)
          continue;
        const double c_ang = ck_vB * ck_vB / double(Angular::twoj_k(kv) + 1);

        // Calc and store first? Don't depend on k (but much faster than pi)
        auto g_beta =
            basis_for_Green
                ? Green(kB, env + omre, omim, States::both, GrMethod::basis)
                : GreenAtComplex(gBetas[std::size_t(ibeta)], omim);

        // sum over k:
        Sigma_w += c_ang * (g_beta.mult_elements_by(qpq_dw)).get_real();

      } // k
    }   // omega
  }     // beta

  for (const auto &sw : Sigma_Ws) {
    Sigma += sw;
  }

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
GMatrix FeynmanSigma::FeynmanEx_1(int kv, double env) {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  // XXX Sort of works? Out by factor? Unstable?

  GMatrix Sx_e1(m_subgrid_points, m_include_G);
  std::vector<GMatrix> Sx_e1s(std::size_t(m_max_kappaindex + 1),
                              {m_subgrid_points, m_include_G});

  if (exclude_exchange)
    return Sx_e1;

  // // Set up imaginary frequency grid:
  const double omre = m_omre; // does seem to depend on this..
  const auto &wgrid = *m_wgridX;

  // Greens function (at om_re) remains same, so calculate it once only:
  std::vector<ComplexGMatrix> gws;
  for (int ik = 0; ik <= m_max_kappaindex; ++ik) {
    const auto kappa = Angular::kappaFromIndex(ik);
    gws.push_back(Green_hf(kappa, env + omre));
  }

  const auto &core = p_hf->get_core();

  const auto gr_meth =
      basis_for_Green ? GrMethod::basis : GrMethod::Green; // XXX

#pragma omp parallel for
  for (int ibeta = 0; ibeta <= m_max_kappaindex; ++ibeta) {
    const auto kB = Angular::kappaFromIndex(ibeta);
    auto &Sx_e1B = Sx_e1s[std::size_t(ibeta)];

    for (auto ia = 0ul; ia < core.size(); ++ia) {
      const auto &Fa = core[ia];
      if (Fa.n < m_min_core_n)
        continue;
      const auto ka = Fa.k;

      for (auto iw = 0ul; iw < wgrid.num_points; iw++) {
        const auto omim = wgrid.r[iw];
        const auto dw1 = wgrid.drdu[iw] * wgrid.du;

        const auto gxBm = Green_ex(kB, Fa.en - omre, -omim, gr_meth); // B,a,w
        const auto gxBp = Green_ex(kB, Fa.en + omre, omim, gr_meth);  // B,a,w

        for (int ialpha = 0; ialpha <= m_max_kappaindex; ++ialpha) {
          const auto kA = Angular::kappaFromIndex(ialpha);
          const auto &gA_re = gws[std::size_t(ialpha)];

          // const auto gA = GreenAtComplex(gA_re, omim); // A, w
          const auto gA = basis_for_Green ? Green(kA, env + omre, omim,
                                                  States::both, GrMethod::basis)
                                          : GreenAtComplex(gA_re, omim);

          const auto &pa = m_Pa[ia];

          const auto gqpg = sumkl_GQPGQ(gA, gxBm, gxBp, pa, kv, kA, kB, ka);

          const auto fudge_factor = 1.0 / (2 * M_PI); //???
          Sx_e1B += (fudge_factor * dw1 / M_PI) * gqpg;

        } // alpha
      }   // w
    }     // beta
  }       // a

  for (const auto &sB : Sx_e1s) {
    Sx_e1 += sB;
  }

  // devide through by dri, drj [these included in q's, but want
  // differential operator for sigma] or.. include one of these in
  // definition of opertion S|v> ?
  for (auto i = 0ul; i < m_subgrid_points; ++i) {
    const auto dri = dr_subToFull(i);
    for (auto j = 0ul; j < m_subgrid_points; ++j) {
      const auto drj = dr_subToFull(j);
      Sx_e1.ff[i][j] /= (dri * drj);
      // No g
    }
  }

  return Sx_e1;
}

//******************************************************************************
GMatrix FeynmanSigma::sumk_cGQPQ(int kv, int ka, int kA, int kB,
                                 const ComplexGMatrix &g_beta,
                                 const ComplexGMatrix &pi_aalpha) const {
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);
  // sum over k: (for DIRECT part:)
  // sum_k [ck qk * pi(w) * qk], ck angular factor

  auto gqpq = ComplexGMatrix(m_subgrid_points, m_include_G);
  // min/max k to include in sum [just to save calculating terms=0]
  // based on two Ck angular factors
  const auto [kmin1, kmax1] = Angular::kminmax_Ck(kv, kB); // Ck_vB
  const auto [kmin2, kmax2] = Angular::kminmax_Ck(ka, kA); // Ck_aA
  const auto kmin = std::max(kmin1, kmin2);
  const auto kmax = std::min(kmax1, kmax2);
  const auto tjvp1 = Angular::twoj_k(kv) + 1; //[jv]
  for (int k = kmin; k <= kmax; ++k) {
    const auto c1 = Angular::Ck_kk(k, kv, kB);
    if (c1 == 0.0)
      continue;
    const auto c2 = Angular::Ck_kk(k, ka, kA);
    if (c2 == 0.0)
      continue;
    const double c_ang = c1 * c1 * c2 * c2 / double(tjvp1 * (2 * k + 1));
    const auto &q = get_qk(k);
    const auto &q_scr = screen_Coulomb ? screenedCoulomb(q, pi_aalpha) : q;
    gqpq += c_ang * q * pi_aalpha * q_scr;
    // nb: can screen either q, but not both!
  } // k
  return gqpq.mult_elements_by(g_beta).get_real();
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
GMatrix FeynmanSigma::sumkl_GQPGQ(const ComplexGMatrix &gA,
                                  const ComplexGMatrix &gxBm,
                                  const ComplexGMatrix &gxBp,
                                  const ComplexGMatrix &pa, int kv, int kA,
                                  int kB, int ka) const {
  // EXCHANGE part:
  [[maybe_unused]] auto sp = IO::Profile::safeProfiler(__func__);

  auto sum_GQPG = GMatrix(m_subgrid_points, m_include_G);

  const auto &Ck = m_yeh.Ck();

  const auto kmin = 0;
  const int k_cut = 6; // XXX TEMP!
  const auto kmax = std::min(m_maxk, k_cut);

  const auto tjv = Angular::twoj_k(kv);
  const auto tjA = Angular::twoj_k(kA);
  const auto tjB = Angular::twoj_k(kB);
  const auto tja = Angular::twoj_k(ka);
  const auto tjvp1 = tjv + 1; //[jv]

  for (auto k = kmin; k <= kmax; ++k) {
    const auto CkvA = Ck(k, kv, kA);
    if (CkvA == 0.0)
      continue;
    const auto CkaB = Ck(k, ka, kB);
    if (CkaB == 0.0)
      continue;
    const auto ck1 = CkvA * CkaB;
    // nb: differ only by (at most) sign
    const auto ck2 = CkvA * Ck(k, kB, ka);

    const auto &qk = get_qk(k); // XXX this one

    for (auto l = kmin; l <= kmax; ++l) {
      const auto cl1 = Ck(l, kv, kB) * Ck(l, ka, kA);
      const auto cl2 = Ck(l, kv, ka) * Ck(l, kB, kA);
      if (cl1 == 0.0 && cl2 == 0.0)
        continue;

      const auto sj1 = cl1 == 0.0 ? 0.0 : m_6j(tjv, tjA, tja, tjB, k, l);
      const auto sj2 =
          cl2 == 0.0 ? 0.0 : tja == tjB ? sj1 : m_6j(tjv, tjA, tjB, tja, k, l);
      if (sj1 == 0.0 && sj2 == 0.0)
        continue;

      const auto &ql = get_qk(l);

      const auto s0 = Angular::evenQ(k + l) ? 1 : -1;
      const auto cang1 = s0 * ck1 * cl1 * sj1 / tjvp1;
      const auto cang2 = s0 * ck2 * cl2 * sj2 / tjvp1;
      const auto ic1 = ComplexDouble(0.0, cang1);
      const auto ic2 = ComplexDouble(0.0, cang2);

      if (sj1 != 0.0)
        tensor_5_product(&sum_GQPG, ic1, qk, pa, gxBm, gA, ql);
      if (sj2 != 0.0)
        tensor_5_product(&sum_GQPG, ic2, qk, gxBp, pa, gA, ql);

    } // l
  }   // k

  return sum_GQPG;
}

} // namespace MBPT
